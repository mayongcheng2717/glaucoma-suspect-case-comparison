"""Functional anchoring for the two stage-specific glaucoma loci.

For each lead variant: genes in the window (Ensembl), whether the variant is a
cis-eQTL for those genes in GTEx v8 (via the eQTL Catalogue tabix files), and
its phenome-wide profile in FinnGen.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:\eye100\code")
from dse import finngen as fg, genes as GS, stats  # noqa: E402
from dse.netsession import get_json  # noqa: E402
from dse.remoteio import Tabix  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = r"D:\eye100\results\gwcas"
TBI = r"D:\eye100\cache\tbi"
FT = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/"
TISSUES = ["Whole_Blood", "Brain_Frontal_Cortex_BA9", "Artery_Tibial", "Cells_Cultured_fibroblasts"]

LEADS = [
    {"name": "9p21 (CDKN2B)", "chrom": "9", "pos": 22_041_444},
    {"name": "1q24 (TMCO1)", "chrom": "1", "pos": 165_754_533},
]

out = {}
for L in LEADS:
    ch, pos = L["chrom"], L["pos"]
    print(f"\n=== {L['name']}  chr{ch}:{pos:,}")
    g = GS.genes_in(ch, max(1, pos - 300_000), pos + 300_000)
    g = sorted(g, key=lambda x: abs((x["start"] + x["end"]) / 2 - pos))
    print("  genes within 300 kb:", ", ".join(x["symbol"] for x in g[:12]))
    rec = {"lead": L, "genes": g[:12], "eqtl": []}
    # nearest gene and the two canonical candidates
    test_genes = [x["symbol"] for x in g[:4]]
    if ch == "9":
        test_genes = ["CDKN2B", "CDKN2B-AS1", "CDKN2A", "MTAP"]
    if ch == "1":
        test_genes = ["TMCO1", "ALDH9A1", "TMEM132B", "MAEL"]
    for sym in test_genes:
        gi = GS.gene(sym)
        if not gi:
            continue
        ens = gi["id"]
        for tis in TISSUES:
            url = f"{FT}{tis}.tsv.gz"
            try:
                tb = Tabix(url, TBI)
            except Exception as e:  # noqa: BLE001
                print(f"    {tis}: tabix failed {type(e).__name__}")
                continue
            lo, hi = pos - 500_000, pos + 500_000
            rows = []
            for line in tb.query(ch, lo, hi):
                f = line.split("\t")
                if len(f) < 16:
                    continue
                if f[6].split(".")[0] != ens.split(".")[0]:
                    continue
                try:
                    rows.append((int(f[13]), f[14], f[15], float(f[8]), float(f[9]),
                                 float(f[2])))
                except Exception:  # noqa: BLE001
                    continue
            if not rows:
                continue
            df = pd.DataFrame(rows, columns=["pos", "ref", "alt", "beta", "se", "p"])
            hit = df[df["pos"] == pos]
            best = df.loc[df["p"].idxmin()]
            rec["eqtl"].append({
                "gene": sym, "ensembl": ens, "tissue": tis, "n_variants": int(len(df)),
                "lead_is_eqtl": bool(len(hit) and hit["p"].iloc[0] < 1e-4),
                "lead_p": float(hit["p"].iloc[0]) if len(hit) else None,
                "lead_beta": float(hit["beta"].iloc[0]) if len(hit) else None,
                "best_variant_pos": int(best["pos"]), "best_p": float(best["p"]),
                "best_beta": float(best["beta"]),
            })
            if len(hit) and hit["p"].iloc[0] < 1e-4:
                print(f"    {sym:12s} {tis:28s} LEAD IS eQTL  beta={hit['beta'].iloc[0]:+.3f} "
                      f"p={hit['p'].iloc[0]:.1e}")
            else:
                print(f"    {sym:12s} {tis:28s} best p={best['p']:.1e} at {int(best['pos']):,}")
    out[L["name"]] = rec

json.dump(out, open(os.path.join(OUT, "eqtl_anchoring.json"), "w", encoding="utf-8"),
          indent=1, default=str)
print("\nsaved", os.path.join(OUT, "eqtl_anchoring.json"))
