"""Recompute the case-control OR table for the two lead variants."""
import sys

import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(0, r"D:\eye100\code")

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

spec = pd.read_csv(r"D:\eye100\results\gwcas\lead_specificity.tsv", sep="\t")
N_TOTAL = 500186
rows = []
for r in spec.itertuples():
    ncase = float(r.n_cases)
    n = float(N_TOTAL)
    afc, af = float(r.af_cases), float(r.af_all)
    nctrl = n - ncase
    afctl = (af * n - afc * ncase) / nctrl
    lo = np.log((afc / (1 - afc)) / (afctl / (1 - afctl)))
    se = np.sqrt(1 / (2 * ncase * afc * (1 - afc))
                 + 1 / (2 * nctrl * afctl * (1 - afctl)))
    rows.append({"lead": r.lead, "phenotype": r.phenotype, "label": r.label,
                 "or": float(np.exp(lo)), "lo": float(np.exp(lo - 1.96 * se)),
                 "hi": float(np.exp(lo + 1.96 * se)),
                 "p": float(2 * st.norm.sf(abs(lo / se)))})
orr = pd.DataFrame(rows)
orr.to_csv(r"D:\eye100\results\gwcas\lead_case_control.tsv", sep="\t", index=False)
print(orr.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print("saved", len(orr), "rows")
