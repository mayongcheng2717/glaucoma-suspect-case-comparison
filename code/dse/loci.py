"""Disease locus catalogues from GWAS Catalog (EFO-filtered associations).

``/api/v2/associations?efoTrait=<text>`` is the only filter that actually
restricts results; ``locations`` come back as ``"chr:pos"`` strings (GRCh38).
"""

from __future__ import annotations

import json
import os
import urllib.parse
from typing import Dict, List

import sys
sys.path.insert(0, r"D:\eye100\code")
from dse.netsession import get_json  # noqa: E402

CACHE = r"D:\eye100\cache\loci"
BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2/"

# EFO trait labels as used by the GWAS Catalog v2 search
TRAITS = {
    "AMD": ["age-related macular degeneration"],
    "DR": ["diabetic retinopathy"],
    "GLAUCOMA": ["open-angle glaucoma", "intraocular pressure measurement"],
    "MYOPIA": ["myopia", "refractive error"],
    "RVO": ["retinal vein occlusion"],
}
PMax = 5e-8


def fetch_efo(trait: str, size: int = 500, max_pages: int = 6) -> List[dict]:
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, f"efo_{urllib.parse.quote(trait)}.json")
    if os.path.exists(fp):
        with open(fp, encoding="utf-8") as f:
            return json.load(f)
    out, page = [], 0
    total = None
    while page < max_pages:
        d = get_json(f"{BASE}associations?efoTrait={urllib.parse.quote(trait)}"
                     f"&size={size}&page={page}")
        if not d:
            break
        if total is None:
            total = d.get("page", {}).get("totalElements", 0)
            print(f"  [{trait}] {total} associations")
        items = d.get("_embedded", {}).get("associations", [])
        if not items:
            break
        for a in items:
            try:
                p = float(a.get("p_value"))
            except (TypeError, ValueError):
                continue
            if p > PMax:
                continue
            locs = a.get("locations") or []
            for loc in locs:
                if ":" not in str(loc):
                    continue
                ch, pos = str(loc).split(":")[:2]
                if not ch.isdigit() or not pos.isdigit():
                    continue
                out.append({"chrom": ch, "pos": int(pos), "p": p,
                            "genes": a.get("mapped_genes") or [],
                            "pubmed": a.get("pubmed_id"),
                            "trait": (a.get("reported_trait") or [None])[0]})
        page += 1
        if page * size >= (total or 0):
            break
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(out, f)
    print(f"  [{trait}] kept {len(out)} GWS associations")
    return out


def cluster(assocs: List[dict], flank: int = 500_000) -> List[dict]:
    by: Dict[str, List[dict]] = {}
    for a in assocs:
        by.setdefault(a["chrom"], []).append(a)
    loci = []
    for ch, lst in by.items():
        lst.sort(key=lambda x: x["pos"])
        cur = None
        for a in lst:
            if cur is None or a["pos"] - cur["end"] > flank:
                if cur:
                    loci.append(cur)
                cur = {"chrom": ch, "start": a["pos"], "end": a["pos"], "members": [a]}
            else:
                cur["end"] = max(cur["end"], a["pos"])
                cur["members"].append(a)
        if cur:
            loci.append(cur)
    for l in loci:
        best = min(l["members"], key=lambda m: m["p"])
        l["lead_pos"] = best["pos"]
        l["lead_p"] = best["p"]
        l["genes"] = sorted({g for m in l["members"] for g in (m["genes"] or [])})
        l["pubmed"] = sorted({str(m["pubmed"]) for m in l["members"] if m["pubmed"]})
        l["mid"] = (l["start"] + l["end"]) // 2
        l["id"] = f"{l['chrom']}:{l['mid'] / 1e6:.2f}Mb"
        l["n_assoc"] = len(l["members"])
    loci.sort(key=lambda x: (int(x["chrom"]), x["mid"]))
    return loci


def locus_set(disease: str, flank: int = 500_000) -> List[dict]:
    assocs = []
    for t in TRAITS[disease]:
        assocs.extend(fetch_efo(t))
    seen, uniq = set(), []
    for a in assocs:
        k = (a["chrom"], a["pos"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(a)
    loci = cluster(uniq, flank)
    print(f"[{disease}] {len(uniq)} unique GWS positions -> {len(loci)} loci")
    return loci


if __name__ == "__main__":
    import sys
    want = sys.argv[1:] or list(TRAITS)
    bank = {}
    for d in want:
        bank[d] = locus_set(d)
        for l in bank[d][:8]:
            print(f"   {l['id']:14s} p={l['lead_p']:.1e} n={l['n_assoc']:3d} {','.join(l['genes'][:4])[:56]}")
    os.makedirs(CACHE, exist_ok=True)
    with open(os.path.join(CACHE, "locus_bank.json"), "w", encoding="utf-8") as f:
        json.dump(bank, f, indent=0)
    print("\nsummary:", {k: len(v) for k, v in bank.items()})
