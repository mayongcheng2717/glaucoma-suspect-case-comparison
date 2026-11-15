"""Prior-art check on the genome-wide one-locus claim and the LOXL1 contamination."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, r"D:\eye100\code")
from dse import litcheck as LC  # noqa: E402

OUT = r"D:\eye100\results\gwcas"
os.makedirs(OUT, exist_ok=True)

CASES = [
    {"claim": "Genome-wide, neovascular and atrophic AMD differ at only one locus (10q26)",
     "queries": ['("neovascular" OR "wet") AND ("atrophic" OR "geographic atrophy") AND '
                 '"macular degeneration" AND ("genome-wide" OR GWAS) AND (differ OR difference '
                 'OR subtype)',
                 '"case-case" "macular degeneration" neovascular atrophic',
                 'CC-GWAS "macular degeneration" subtype'],
     "must_terms": [["macular"], ["genome-wide", "gwas", "case-case", "loci", "locus"]]},

    {"claim": "The 10q26 locus is the only genome-wide significant determinant of AMD subtype",
     "queries": ['"10q26" AND (neovascular OR atrophic OR subtype) AND (genome-wide OR GWAS)',
                 '"ARMS2" "HTRA1" neovascular atrophic "genome-wide" difference'],
     "must_terms": [["10q26", "arms2", "htra1"]]},

    {"claim": "The FinnGen primary open-angle glaucoma endpoint contains exfoliation cases, "
              "detectable because LOXL1 separates NTG from POAG genome-wide",
     "queries": ['"LOXL1" AND ("open-angle glaucoma" OR POAG) AND ("exfoliation" OR '
                 'pseudoexfoliation) AND (misclassif OR contamina OR definition OR coding)',
                 '"exfoliation glaucoma" "primary open angle glaucoma" genetic overlap '
                 'misclassification'],
     "must_terms": [["loxl1"], ["exfoliation", "pseudoexfoliation"]]},

    {"claim": "A genome-wide case-case analysis of retinal disease subtypes",
     "queries": ['"case-case" AND ("retinal" OR "macular" OR "glaucoma") AND '
                 '(genome-wide OR GWAS)',
                 '"cases-only" GWAS retinal disease subtype comparison'],
     "must_terms": [["case-case", "cases-only", "case only"]]},
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
        for h in r["hits"][:4]:
            print(f"    near  {h['year']}  {h['title'][:100]}")

json.dump(res, open(os.path.join(OUT, "prior_art_onelocus.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
clean = [r for r in res if r["verdict"] == "NO_DIRECT_HIT"]
print(f"\n\n{len(clean)}/{len(res)} claims have no direct prior-art title")
