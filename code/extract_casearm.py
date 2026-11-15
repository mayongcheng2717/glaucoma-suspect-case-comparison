"""Extract case-arm allele frequencies genome-wide for a set of FinnGen endpoints.

Only what the case-case contrast needs is kept: position, chromosome, the
alt-allele frequency among cases, and a hash of the allele pair. No complete
summary-statistics file is retained, and memory stays at one chunk.

Output: D:/eye100/cache/casearm/<PHENO>.npz
    key   (int64)   chrom<<57 | pos<<26 | hash26(ref/alt)
    pos   (int32)
    chrom (int8)
    afc   (float32) alt-allele frequency among cases
    af    (float32) alt-allele frequency overall
"""

from __future__ import annotations

import gzip
import io
import os
import sys
import time
import urllib.request
import zlib

import numpy as np
import pandas as pd

BASE = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
OUT = r"D:\eye100\cache\casearm"
ALL_COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp",
            "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]
USE = ["chrom", "pos", "ref", "alt", "af_alt", "af_alt_cases"]
CHROM_CODE = {str(i): i for i in range(1, 23)}
CHROM_CODE.update({"X": 23, "Y": 24, "MT": 25, "M": 25})
MAF_MIN = 0.005


def extract(pheno: str, chunk_rows: int = 400_000):
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, f"{pheno}.npz")
    if os.path.exists(dst):
        print(f"[{pheno}] already done", flush=True)
        return dst
    url = BASE + f"finngen_R13_{pheno}.gz"
    t0 = time.time()
    keys, poss, chroms, afcs, afs = [], [], [], [], []
    total = 0
    req = urllib.request.Request(url, headers={"User-Agent": "dse/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        reader = pd.read_csv(
            io.TextIOWrapper(gz, encoding="utf-8", newline=""),
            sep="\t", chunksize=chunk_rows, dtype=str,
            names=ALL_COLS, header=0, usecols=USE, on_bad_lines="skip",
        )
        for ch in reader:
            total += len(ch)
            cc = ch["chrom"].map(CHROM_CODE)
            pos = pd.to_numeric(ch["pos"], errors="coerce")
            af = pd.to_numeric(ch["af_alt"], errors="coerce")
            afc = pd.to_numeric(ch["af_alt_cases"], errors="coerce")
            ok = cc.notna() & pos.notna() & af.notna() & afc.notna()
            ok &= (af > MAF_MIN) & (af < 1 - MAF_MIN)
            if not ok.any():
                continue
            cc = cc[ok].astype(np.int64).to_numpy()
            pos = pos[ok].astype(np.int64).to_numpy()
            af = af[ok].to_numpy(np.float32)
            afc = afc[ok].to_numpy(np.float32)
            al = (ch["ref"][ok].astype(str) + "/" + ch["alt"][ok].astype(str)).to_numpy()
            h = np.fromiter((zlib.crc32(a.encode()) for a in al), np.int64, len(al))
            h &= 0x3FFFFFF
            keys.append((cc << 57) | (pos << 26) | h)
            poss.append(pos.astype(np.int32))
            chroms.append(cc.astype(np.int8))
            afcs.append(afc)
            afs.append(af)
            if total % 4_000_000 < chunk_rows:
                print(f"[{pheno}] {total:,} rows, kept {sum(len(k) for k in keys):,}, "
                      f"{time.time() - t0:.0f}s", flush=True)
    key = np.concatenate(keys)
    arrs = {"key": key, "pos": np.concatenate(poss), "chrom": np.concatenate(chroms),
            "afc": np.concatenate(afcs), "af": np.concatenate(afs)}
    order = np.argsort(key, kind="stable")
    for k in arrs:
        arrs[k] = arrs[k][order]
    np.savez(dst, **arrs)
    print(f"[{pheno}] DONE total={total:,} kept={len(key):,} in {time.time() - t0:.0f}s "
          f"({os.path.getsize(dst) / 1e6:.0f} MB)", flush=True)
    return dst


if __name__ == "__main__":
    for p in sys.argv[1:]:
        try:
            extract(p)
        except Exception as e:  # noqa: BLE001
            print(f"[{p}] FAILED {type(e).__name__}: {e}", flush=True)
