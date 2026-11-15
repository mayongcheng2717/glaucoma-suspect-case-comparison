"""Resolve references against Europe PMC / Crossref and format them.

References are never written from memory: each entry is matched to a real
record and the DOI, title, journal, year and authors are taken from it.
Europe PMC is tried first because its title index is exact; Crossref is the
fallback.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.parse
from difflib import SequenceMatcher
from typing import List, Optional

sys.path.insert(0, r"D:\eye100\code")
from dse.netsession import get_json  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CACHE = r"D:\eye100\cache\refs"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF = "https://api.crossref.org/works"
BAD_CONTAINERS = ("yearbook", "scivve", "sci vee", "reference reviews", "faculty opinions")


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def _epmc(title: str) -> Optional[dict]:
    q = f'TITLE:"{title}"'
    d = get_json(EPMC, params={"query": q, "format": "json", "pageSize": 8,
                               "resultType": "core"})
    if not d:
        return None
    best, best_sc = None, 0.0
    for it in d.get("resultList", {}).get("result", []):
        t = it.get("title", "")
        sc = SequenceMatcher(None, _norm(title), _norm(t)).ratio()
        if sc > best_sc:
            best, best_sc = it, sc
    if best is None or best_sc < 0.90:
        return None
    return {"src": "europepmc", "title": best.get("title", ""),
            "journal": (best.get("journalInfo", {}) or {}).get("journal", {}).get("title")
                       or best.get("journalTitle", ""),
            "year": best.get("pubYear"),
            "doi": best.get("doi"),
            "pmid": best.get("pmid"),
            "authorString": best.get("authorString", ""),
            "volume": (best.get("journalInfo", {}) or {}).get("volume"),
            "issue": (best.get("journalInfo", {}) or {}).get("issue"),
            "page": best.get("pageInfo", ""),
            "score": round(best_sc, 3)}


def _crossref(title: str) -> Optional[dict]:
    d = get_json(CROSSREF, params={"query.bibliographic": title, "rows": 10,
                                   "select": "DOI,title,container-title,issued,author,volume,issue,page,type"})
    if not d:
        return None
    best, best_sc = None, 0.0
    for it in d.get("message", {}).get("items", []):
        t = (it.get("title") or [""])[0]
        cont = ((it.get("container-title") or [""])[0] or "").lower()
        if any(b in cont for b in BAD_CONTAINERS):
            continue
        if it.get("type") not in ("journal-article", None):
            continue
        sc = SequenceMatcher(None, _norm(title), _norm(t)).ratio()
        if sc > best_sc:
            best, best_sc = it, sc
    if best is None or best_sc < 0.85:
        return None
    yr = None
    try:
        yr = str(best["issued"]["date-parts"][0][0])
    except Exception:  # noqa: BLE001
        pass
    return {"src": "crossref", "title": (best.get("title") or [""])[0],
            "journal": (best.get("container-title") or [""])[0], "year": yr,
            "doi": best.get("DOI"), "pmid": None,
            "authorString": "; ".join(
                f"{a.get('family','')} {a.get('given','')}" for a in (best.get("author") or [])),
            "volume": best.get("volume"), "issue": best.get("issue"),
            "page": best.get("page"), "score": round(best_sc, 3)}


def _authors_from(s: str, max_n: int = 3) -> str:
    if not s:
        return ""
    parts = [p.strip() for p in re.split(r"[;,]", s) if p.strip()]
    out = []
    for p in parts[:max_n]:
        toks = p.split()
        if len(toks) >= 2 and toks[0].isupper() and len(toks[0]) <= 4:
            out.append(f"{toks[0]} {''.join(t[0] for t in toks[1:])}")   # PubMed style
        else:
            fam = toks[0] if toks else ""
            ini = "".join(t[0] for t in toks[1:]) if len(toks) > 1 else ""
            out.append(f"{fam} {ini}".strip())
    s2 = ", ".join(out)
    if len(parts) > max_n:
        s2 += ", et al."
    return s2


JOURNAL_FIX = {
    "plos genetics": "PLoS Genet",
    "nature genetics": "Nat Genet",
    "nature reviews. genetics": "Nat Rev Genet",
    "science (new york, n.y.)": "Science",
    "proceedings of the national academy of sciences of the united states of america": "Proc Natl Acad Sci USA",
    "experimental eye research": "Exp Eye Res",
    "genetics in medicine : official journal of the american college of medical genetics": "Genet Med",
    "royal society open science": "R Soc Open Sci",
    "journal of the royal statistical society. series b, statistical methodology": "J R Stat Soc Series B",
    "taiwan journal of ophthalmology": "Taiwan J Ophthalmol",
    "bmc medical genomics": "BMC Med Genomics",
    "ophthalmology": "Ophthalmology",
    "nature": "Nature",
}


def _fix_journal(j: str) -> str:
    if not j:
        return ""
    k = j.strip().lower()
    return JOURNAL_FIX.get(k, j.strip())


def vancouver(r: dict) -> str:
    au = _authors_from(r.get("authorString", ""))
    t = (r.get("title") or "").strip().rstrip(".")
    bits = []
    if au:
        bits.append(au if au.endswith(".") else au + ".")
    bits.append(t + ".")
    j = _fix_journal(r.get("journal") or "")
    if j:
        bits.append(j.rstrip(".") + ".")
    yv = r.get("year") or ""
    loc = yv
    if r.get("volume"):
        loc += f";{r['volume']}"
    if r.get("issue"):
        loc += f"({r['issue']})"
    if r.get("page"):
        loc += f":{r['page']}"
    if loc.strip(";"):
        bits.append(loc.strip(";") + ".")
    if r.get("doi"):
        bits.append(f"doi:{r['doi']}")
    elif r.get("pmid"):
        bits.append(f"PMID:{r['pmid']}")
    return " ".join(b for b in bits if b)


def resolve(title: str) -> dict:
    os.makedirs(CACHE, exist_ok=True)
    key = re.sub(r"[^A-Za-z0-9]+", "_", title)[:110]
    fp = os.path.join(CACHE, key + ".json")
    if os.path.exists(fp):
        with open(fp, encoding="utf-8") as f:
            return json.load(f)
    rec = _epmc(title) or _crossref(title)
    out = {"query": title, "found": rec is not None, "record": rec,
           "ref": vancouver(rec) if rec else None}
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False)
    return out


def build(titles: List[str]) -> List[dict]:
    out = []
    for t in titles:
        out.append(resolve(t))
        time.sleep(0.15)
    return out


if __name__ == "__main__":
    titles = [l.strip() for l in open(sys.argv[1], encoding="utf-8") if l.strip()]
    res = build(titles)
    ok = 0
    for i, r in enumerate(res, 1):
        flag = "OK  " if r["found"] else "MISS"
        if r["found"]:
            ok += 1
        print(f"{flag} [{i}] {r['ref'] or r['query']}")
    print(f"\nmatched {ok}/{len(res)}")
    out = sys.argv[2] if len(sys.argv) > 2 else r"D:\eye100\results\refs_auto.json"
    json.dump(res, open(out, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("saved", out)
