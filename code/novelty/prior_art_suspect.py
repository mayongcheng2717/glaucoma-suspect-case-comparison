"""Prior-art check: glaucoma suspects versus confirmed open-angle glaucoma."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, r"D:\eye100\code")
from dse import litcheck as LC  # noqa: E402

OUT = r"D:\eye100\results\gwcas"
CASES = [
    {"claim": "Genome-wide comparison of glaucoma suspects with confirmed primary "
              "open-angle glaucoma",
     "queries": ['("glaucoma suspect" OR "glaucoma suspects" OR "ocular hypertension") AND '
                 '("open-angle glaucoma") AND (genome-wide OR GWAS OR polygenic)',
                 '"ocular hypertension" "primary open angle glaucoma" genetic risk conversion'],
     "must_terms": [["glaucoma suspect", "ocular hypertension", "hypertension"],
                    ["open-angle", "open angle", "poag"]]},

    {"claim": "Genetic determinants of conversion from ocular hypertension to glaucoma",
     "queries": ['("ocular hypertension") AND (conversion OR converter OR progression) AND '
                 '(genetic OR polygenic OR GWAS)',
                 '"glaucoma suspect" conversion "polygenic risk"'],
     "must_terms": [["convert", "progress", "incident", "polygenic"]]},

    {"claim": "Case-case genome-wide contrast between glaucoma suspects and glaucoma cases",
     "queries": ['"case-case" AND glaucoma',
                 '("cases-only" OR "case only") AND glaucoma AND genetic'],
     "must_terms": [["glaucoma"]]},

    {"claim": "Genome-wide case-case analysis of retinal or ophthalmic disease phenotypes",
     "queries": ['"case-case genome-wide" OR "case-case GWAS"',
                 '"CC-GWAS" AND (ophthalm* OR retinal OR glaucoma OR macular)'],
     "must_terms": [["case-case", "cc-gwas", "case-case genome-wide"]]},
]

res = []
for c in CASES:
    r = LC.assess(c["claim"], c["queries"], c.get("must_terms"))
    res.append(r)
    print(f"\n{'=' * 78}\n{r['claim'][:130]}")
    print(f"  verdict: {r['verdict']}   (hits {r['n_hits']})")
    for h in r["prior_art"][:5]:
        print(f"    PRIOR {h['year']}  {h['title'][:100]}")
        print(f"          {str(h['journal'])[:48]}  doi:{h['doi']}")
    if not r["prior_art"]:
        for h in r["hits"][:5]:
            print(f"    near  {h['year']}  {h['title'][:100]}")

json.dump(res, open(os.path.join(OUT, "prior_art_suspect.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
clean = [r for r in res if r["verdict"] == "NO_DIRECT_HIT"]
print(f"\n\n{len(clean)}/{len(res)} claims have no direct prior-art title")
