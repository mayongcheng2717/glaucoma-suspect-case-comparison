"""Stream a FinnGen R13 endpoint and keep only suggestive variants.

Memory footprint stays at one chunk (~300 MB) so several endpoints can be
streamed at once on a 16 GB machine.  The output is the complete set of
variants with p < 1e-5, which is all that is needed to define loci and to
compare subtypes at every locus that either subtype nominates.
"""

from __future__ import annotations

import gzip
import io
import os
import sys
import time
import urllib.request

import pandas as pd

BASE = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
OUT = r"D:\eye100\data\tophits"
USE = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "beta", "sebeta",
       "af_alt", "af_alt_cases"]
PMAX = 1e-5


def extract(pheno: str, chunk_rows: int = 300_000):
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, f"{pheno}.tsv.gz")
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        print(f"[{pheno}] already done", flush=True)
        return dst
    url = BASE + f"finngen_R13_{pheno}.gz"
    t0 = time.time()
    kept, total = [], 0
    req = urllib.request.Request(url, headers={"User-Agent": "dse/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        reader = pd.read_csv(
            io.TextIOWrapper(gz, encoding="utf-8", newline=""),
            sep="\t", chunksize=chunk_rows, dtype=str,
            names=USE, header=0, usecols=USE, on_bad_lines="skip",
        )
        for ch in reader:
            total += len(ch)
            p = pd.to_numeric(ch["pval"], errors="coerce")
            sub = ch.loc[p < PMAX]
            if len(sub):
                kept.append(sub)
            if total % 3_000_000 < chunk_rows:
                print(f"[{pheno}] {total:,} rows, kept {sum(len(k) for k in kept):,}, "
                      f"{time.time() - t0:.0f}s", flush=True)
    if not kept:
        print(f"[{pheno}] no variants under {PMAX}", flush=True)
        return None
    df = pd.concat(kept, ignore_index=True)
    for c in ("pos", "pval", "beta", "sebeta", "af_alt", "af_alt_cases"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["pos", "pval"]).sort_values("pos")
    df["pos"] = df["pos"].astype("int64")
    df.to_csv(dst, sep="\t", index=False, compression="gzip")
    print(f"[{pheno}] DONE total={total:,} kept={len(df):,} in {time.time() - t0:.0f}s", flush=True)
    return dst


if __name__ == "__main__":
    for p in sys.argv[1:]:
        try:
            extract(p)
        except Exception as e:  # noqa: BLE001
            print(f"[{p}] FAILED: {type(e).__name__}: {e}", flush=True)
