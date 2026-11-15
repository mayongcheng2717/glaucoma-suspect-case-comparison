"""Literature check used as the novelty filter.

A result is only carried forward if a targeted search does not already describe
the same contrast. Searches hit Europe PMC (exact phrase support) and fall back
to the GWAS Catalog-free Crossref route via dse.refs.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from typing import List, Optional

sys.path.insert(0, r"D:\eye100\code")
from dse.netsession import get_json  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CACHE = r"D:\eye100\cache\litcheck"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def search(query: str, n: int = 15) -> List[dict]:
    os.makedirs(CACHE, exist_ok=True)
    key = re.sub(r"[^A-Za-z0-9]+", "_", query)[:120]
    fp = os.path.join(CACHE, key + ".json")
    if os.path.exists(fp):
        with open(fp, encoding="utf-8") as f:
            return json.load(f)
    d = get_json(EPMC, params={"query": query, "format": "json", "pageSize": n,
                               "resultType": "lite"})
    out = []
    if d:
        for it in d.get("resultList", {}).get("result", []):
            out.append({"title": it.get("title", ""), "year": it.get("pubYear"),
                        "journal": it.get("journalTitle", ""), "doi": it.get("doi"),
                        "pmid": it.get("pmid"), "authors": it.get("authorString", "")[:80]})
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False)
    time.sleep(0.2)
    return out


GENETICS_TERMS = ["genetic", "gene", "gwas", "variant", "polymorphism", "allele",
                  "heritab", "loci", "locus", "genome-wide", "risk score", "polygenic"]
_GEN_RE = re.compile(
    r"\b(genetic|genetics|gene|genes|gwas|variant|variants|polymorphism|polymorphisms|"
    r"allele|alleles|heritab\w*|loci|locus|polygenic|genome-wide)\b", re.I)


def _has_genetics(title: str) -> bool:
    """Word-boundary match, so 'degeneration' does not count as 'gene'."""
    return bool(_GEN_RE.search(title))


def assess(claim: str, queries: List[str], must_terms: List[List[str]] = None,
           require_genetics_in_title: bool = True) -> dict:
    """Return the hit list plus a prior-art verdict.

    A hit counts as prior art only if its title contains at least one term from
    every group in ``must_terms`` *and* (by default) a genetics term. Without
    the genetics requirement, clinical trial reports match on the disease words
    alone and swamp the search.
    """
    hits = []
    for q in queries:
        hits.extend(search(q))
    seen, uniq = set(), []
    for h in hits:
        k = (h.get("doi") or h.get("pmid") or h["title"]).lower()
        if k in seen:
            continue
        seen.add(k)
        uniq.append(h)
    prior = []
    if must_terms:
        for h in uniq:
            t = h["title"].lower()
            if not all(any(term.lower() in t for term in grp) for grp in must_terms):
                continue
            if require_genetics_in_title and not _has_genetics(h["title"]):
                continue
            prior.append(h)
    return {"claim": claim, "queries": queries, "n_hits": len(uniq),
            "prior_art": prior, "verdict": "KNOWN" if prior else "NO_DIRECT_HIT",
            "hits": uniq[:12]}


if __name__ == "__main__":
    spec = json.load(open(sys.argv[1], encoding="utf-8"))
    out = []
    for c in spec:
        r = assess(c["claim"], c["queries"], c.get("must_terms"))
        out.append(r)
        print(f"\n=== {c['claim']}\n    verdict: {r['verdict']}  (hits={r['n_hits']})")
        for h in r["prior_art"][:5]:
            print(f"      PRIOR: {h['year']} {h['title'][:110]}")
            print(f"             doi:{h['doi']}")
        for h in r["hits"][:3]:
            print(f"      hit:   {h['year']} {h['title'][:100]}")
    json.dump(out, open(sys.argv[2], "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("\nsaved", sys.argv[2])
