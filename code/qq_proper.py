"""Proper QQ data: random subsample of the full genome-wide z distribution."""
from __future__ import annotations

import json
import sys

import numpy as np

sys.path.insert(0, r"D:\eye100\code")
from dse import finngen as fg  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CACHE = r"D:\eye100\cache\casearm"
OUT = r"D:\eye100\results\gwcas"


def load(p):
    z = np.load(f"{CACHE}\\{p}.npz")
    return {k: z[k] for k in z.files}


A, B = load("H7_GLAUCSUSP"), load("H7_GLAUCOMA_POAG")
nA, nB = fg.n_cases("H7_GLAUCSUSP"), fg.n_cases("H7_GLAUCOMA_POAG")
common, ia, ib = np.intersect1d(A["key"], B["key"], assume_unique=True, return_indices=True)
pa = A["afc"][ia].astype(np.float64)
pb = B["afc"][ib].astype(np.float64)
ok = ((pa > 0.01) & (pa < 0.99) & (pb > 0.01) & (pb < 0.99)
      & (2 * nA * pa >= 60) & (2 * nB * pb >= 60))
pa, pb = pa[ok], pb[ok]
lo = np.log(pa / (1 - pa)) - np.log(pb / (1 - pb))
var = 1 / (2 * nA * pa * (1 - pa)) + 1 / (2 * nB * pb * (1 - pb))
z = lo / np.sqrt(var)
maf = np.minimum(np.minimum(pa, 1 - pa), np.minimum(pb, 1 - pb))
lam = float(np.median(z[maf > 0.10] ** 2) / 0.4549364)
lam = max(lam, 1.0)
print(f"n tested {len(z):,}   lambda {lam:.4f}   sd(z) {z.std():.4f}")

rng = np.random.default_rng(20261115)
idx = rng.choice(len(z), size=min(400_000, len(z)), replace=False)
np.save(f"{OUT}\\qq_sample.npy", z[idx].astype(np.float32))
json.dump({"n_tested": int(len(z)), "lambda": lam, "sd_z": float(z.std()),
           "n_qq_sample": int(len(idx)),
           "median_z2_common": float(np.median(z[maf > 0.10] ** 2))},
          open(f"{OUT}\\qq_stats.json", "w"), indent=1)
print("saved qq_sample.npy and qq_stats.json")
