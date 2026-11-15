"""Specificity of the two stage-specific lead variants across ocular endpoints."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:\eye100\code")
from dse import finngen as fg  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CACHE = r"D:\eye100\cache\casearm"
OUT = r"D:\eye100\results\gwcas"
PHENOS = ["H7_GLAUCSUSP", "H7_GLAUCOMA_POAG", "H7_GLAUCOMA_NTG", "H7_GLAUCOMA_XFG",
          "H7_GLAUCOMA", "WET_AMD", "DRY_AMD"]
LEADS = [("9p21_CDKN2B", "9", 22_041_444), ("1q24_TMCO1", "1", 165_754_533)]

rows = []
for name, ch, pos in LEADS:
    for ph in PHENOS:
        p = os.path.join(CACHE, f"{ph}.npz")
        if not os.path.exists(p):
            continue
        z = np.load(p)
        chrom = (z["key"] >> np.int64(57)).astype(np.int8)
        ppos = ((z["key"] >> np.int64(26)) & np.int64(0x7FFFFFFF)).astype(int)
        sel = (chrom == int(ch)) & (ppos == pos)
        if not sel.any():
            continue
        afc = float(z["afc"][sel][0])
        af = float(z["af"][sel][0])
        try:
            n = fg.n_cases(ph)
            nco = fg.n_controls(ph)
        except Exception:  # noqa: BLE001
            continue
        # case-control effect at this variant (log OR), for context only
        lo = np.log((afc / (1 - afc)) / (af / (1 - af)))
        rows.append({"lead": name, "phenotype": ph, "label": fg.label(ph)[:44],
                     "n_cases": n, "af_cases": afc, "af_all": af,
                     "log_or_vs_all": float(lo), "or_vs_all": float(np.exp(lo))})
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "lead_specificity.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200)
print("=== case-arm allele frequency of each lead across ocular endpoints ===")
for name, sub in df.groupby("lead"):
    print(f"\n{name}")
    print(sub[["phenotype", "n_cases", "af_cases", "af_all", "or_vs_all"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
print("\nsaved")
