"""Targeted prior-art check on the two-stage glaucoma hypothesis."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, r"D:\eye100\code")
from dse import litcheck as LC  # noqa: E402

OUT = r"D:\eye100\results\gwcas"
CASES = [
    {"claim": "TMCO1 and CDKN2B act at different stages of glaucoma: TMCO1 enriches for "
              "being a glaucoma suspect, CDKN2B for having confirmed glaucoma",
     "queries": ['("TMCO1" AND "CDKN2B") AND glaucoma',
                 '"TMCO1" AND glaucoma AND ("intraocular pressure" OR suspect)',
                 '"CDKN2B" AND glaucoma AND (conversion OR progression OR advanced OR severity)'],
     "must_terms": [["tmco1", "cdkn2b"]]},
    {"claim": "Glaucoma suspects and confirmed glaucoma cases carry opposite alleles at "
              "intraocular-pressure versus neurodegeneration loci",
     "queries": ['("glaucoma suspect" OR "ocular hypertension") AND ("TMCO1" OR "CDKN2B")',
                 '"ocular hypertension" "TMCO1"'],
     "must_terms": [["tmco1", "cdkn2b"]]},
    {"claim": "Stage-specific genetic architecture of glaucoma revealed by case-case "
              "comparison of suspects and cases",
     "queries": ['("stage-specific" OR "two-stage" OR "distinct stages") AND glaucoma AND '
                 '(genetic OR loci)',
                 'glaucoma genetic loci "risk of conversion" versus "risk of disease"'],
     "must_terms": [["stage", "conversion", "converter", "two-stage"]]},
]

res = []
for c in CASES:
    r = LC.assess(c["claim"], c["queries"], c.get("must_terms"))
    res.append(r)
    print(f"\n{'=' * 78}\n{r['claim'][:140]}")
    print(f"  verdict: {r['verdict']}   (hits {r['n_hits']})")
    for h in r["prior_art"][:6]:
        print(f"    PRIOR {h['year']}  {h['title'][:105]}")
        print(f"          {str(h['journal'])[:46]}  doi:{h['doi']}")
    # show top hits regardless, so the closest work is visible
    print("  closest titles:")
    for h in r["hits"][:6]:
        print(f"      {h['year']}  {h['title'][:100]}")

json.dump(res, open(os.path.join(OUT, "prior_art_twostage.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
clean = [r for r in res if r["verdict"] == "NO_DIRECT_HIT"]
print(f"\n\n{len(clean)}/{len(res)} claims have no direct prior-art title")
