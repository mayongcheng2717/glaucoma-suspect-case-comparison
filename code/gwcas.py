"""Genome-wide case-case scan from extracted case-arm allele frequencies."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:\eye100\code")
from dse import finngen as fg  # noqa: E402

CACHE = r"D:\eye100\cache\casearm"
OUT = r"D:\eye100\results\gwcas"
os.makedirs(OUT, exist_ok=True)


def load(pheno):
    p = os.path.join(CACHE, f"{pheno}.npz")
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    z = np.load(p)
    return {k: z[k] for k in z.files}


def scan(pa, pb, maf_min=0.01, out_tag=None):
    A = load(pa)
    B = load(pb)
    nA, nB = fg.n_cases(pa), fg.n_cases(pb)
    print(f"\n=== {pa} ({nA:,}) vs {pb} ({nB:,})", flush=True)
    print(f"    {pa}: {len(A['key']):,} variants; {pb}: {len(B['key']):,}", flush=True)
    common, ia, ib = np.intersect1d(A["key"], B["key"], assume_unique=True,
                                   return_indices=True)
    print(f"    shared: {len(common):,}", flush=True)
    pa_af = A["afc"][ia].astype(np.float64)
    pb_af = B["afc"][ib].astype(np.float64)
    ok = ((pa_af > maf_min) & (pa_af < 1 - maf_min)
          & (pb_af > maf_min) & (pb_af < 1 - maf_min)
          & (2 * nA * pa_af >= 60) & (2 * nB * pb_af >= 60))
    common, pa_af, pb_af = common[ok], pa_af[ok], pb_af[ok]
    chrom = (common >> np.int64(57)).astype(np.int8)
    pos = ((common >> np.int64(26)) & np.int64(0x7FFFFFFF)).astype(np.int64)
    print(f"    tested: {len(common):,}", flush=True)

    lo = np.log(pa_af / (1 - pa_af)) - np.log(pb_af / (1 - pb_af))
    var = 1.0 / (2 * nA * pa_af * (1 - pa_af)) + 1.0 / (2 * nB * pb_af * (1 - pb_af))
    z = lo / np.sqrt(var)

    maf = np.minimum(np.minimum(pa_af, 1 - pa_af), np.minimum(pb_af, 1 - pb_af))
    lam_raw = float(np.median(z[maf > 0.10] ** 2) / 0.4549364)
    # A deflated statistic is not rescaled upward: report the raw estimate and test on z.
    lam = max(lam_raw, 1.0)
    zc = z / np.sqrt(lam)
    from scipy import stats as st
    p = 2 * st.norm.sf(np.abs(zc))
    print(f"    lambda = {lam:.3f};  min P = {p.min():.3e};  n(P<5e-8) = {(p < 5e-8).sum()}",
          flush=True)

    df = pd.DataFrame({"chrom": chrom, "pos": pos, "af_a": pa_af, "af_b": pb_af,
                       "or": np.exp(lo), "z": zc, "p": p})
    # lead variants of independent loci
    sig = df[df["p"] < 5e-8].sort_values("p")
    leads = []
    for _, r in sig.iterrows():
        if any(r["chrom"] == L["chrom"] and abs(r["pos"] - L["pos"]) < 1_000_000 for L in leads):
            continue
        leads.append({"chrom": int(r["chrom"]), "pos": int(r["pos"]),
                      "or": float(r["or"]), "z": float(r["z"]), "p": float(r["p"]),
                      "af_a": float(r["af_a"]), "af_b": float(r["af_b"])})
    print(f"    independent loci at P<5e-8: {len(leads)}", flush=True)
    tag = out_tag or f"{pa}_vs_{pb}"
    df.to_parquet(os.path.join(OUT, f"{tag}_all.parquet")) if False else None
    df.to_csv(os.path.join(OUT, f"{tag}_top.tsv.gz"), sep="\t", index=False,
              compression="gzip",
              header=True) if False else None
    df.nlargest(200000, "z", keep="all")
    df.reindex(df["p"].sort_values().index[:200000]).to_csv(
        os.path.join(OUT, f"{tag}_top200k.tsv.gz"), sep="\t", index=False,
        compression="gzip")
    json.dump({"pair": [pa, pb], "n_a": nA, "n_b": nB, "n_tested": int(len(df)),
               "lambda": lam, "lambda_raw": lam_raw, "n_gws": int((p < 5e-8).sum()),
               "n_independent": len(leads), "leads": leads[:60]},
              open(os.path.join(OUT, f"{tag}_leads.json"), "w", encoding="utf-8"),
              indent=1)
    print("    saved", os.path.join(OUT, f"{tag}_leads.json"), flush=True)
    return df, leads


if __name__ == "__main__":
    pairs = [("WET_AMD", "DRY_AMD"), ("H7_GLAUCOMA_NTG", "H7_GLAUCOMA_POAG")]
    for pa, pb in pairs:
        try:
            scan(pa, pb)
        except FileNotFoundError as e:
            print("missing input:", e, flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
