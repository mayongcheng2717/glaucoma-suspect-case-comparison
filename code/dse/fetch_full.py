"""Stream a FinnGen R13 endpoint and reduce it to a compact typed array.

No 760 MB file is kept on disk: the HTTP stream is decompressed and parsed in
chunks, and only the columns needed for cross-phenotype contrast are stored.

Output:  D:/eye100/cache/full/<PHENO>.npz
    chrom (int8)  pos (int32)  ref_id/alt_id (int32, vocabulary in .json)
    beta, se, pval, af, afc (float32)
"""

from __future__ import annotations

import gzip
import io
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

BASE = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
OUT = r"D:\eye100\cache\full"
USE = ["chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt", "af_alt_cases"]
CHROM_CODE = {str(i): i for i in range(1, 23)}
CHROM_CODE.update({"X": 23, "Y": 24, "MT": 25, "M": 25})


def fetch(pheno: str, chunk_rows: int = 2_000_000):
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, f"{pheno}.npz")
    if os.path.exists(dst):
        print(f"[{pheno}] already done", flush=True)
        return dst
    url = BASE + f"finngen_R13_{pheno}.gz"
    t0 = time.time()
    parts = {k: [] for k in ("chrom", "pos", "ref", "alt", "beta", "se", "pval", "af", "afc")}
    ref_vocab, alt_vocab = {}, {}
    n_rows = 0
    req = urllib.request.Request(url, headers={"User-Agent": "dse/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        reader = pd.read_csv(
            io.TextIOWrapper(gz, encoding="utf-8", newline=""),
            sep="\t", chunksize=chunk_rows, dtype=str,
            names=USE, header=0, usecols=USE, on_bad_lines="skip",
        )
        for ch in reader:
            ch = ch.rename(columns={"sebeta": "se", "af_alt": "af", "af_alt_cases": "afc"})
            cc = ch["chrom"].map(CHROM_CODE).fillna(0).astype(np.int8)
            pos = pd.to_numeric(ch["pos"], errors="coerce")
            ok = pos.notna() & (cc > 0)
            cc, ch, pos = cc[ok], ch[ok], pos[ok].astype(np.int64)
            if not len(ch):
                continue
            rid = np.array([ref_vocab.setdefault(a, len(ref_vocab)) for a in ch["ref"]], np.int32)
            aid = np.array([alt_vocab.setdefault(a, len(alt_vocab)) for a in ch["alt"]], np.int32)
            parts["chrom"].append(cc.to_numpy())
            parts["pos"].append(pos.to_numpy())
            parts["ref"].append(rid)
            parts["alt"].append(aid)
            for src, dstk in (("beta", "beta"), ("se", "se"), ("pval", "pval"), ("af", "af"), ("afc", "afc")):
                parts[dstk].append(pd.to_numeric(ch[src], errors="coerce").to_numpy(np.float32))
            n_rows += len(ch)
            if n_rows % 5_000_000 < chunk_rows:
                print(f"[{pheno}] {n_rows:,} rows  {time.time() - t0:.0f}s", flush=True)
    arrs = {k: np.concatenate(v) for k, v in parts.items()}
    arrs["chrom"] = arrs["chrom"].astype(np.int8)
    arrs["pos"] = arrs["pos"].astype(np.int32)
    np.savez_compressed(dst, **arrs)
    with open(os.path.join(OUT, f"{pheno}.vocab.json"), "w") as f:
        json.dump({"ref": ref_vocab, "alt": alt_vocab}, f)
    print(f"[{pheno}] DONE {n_rows:,} rows in {time.time() - t0:.0f}s -> {dst} "
          f"({os.path.getsize(dst) / 1e6:.0f} MB)", flush=True)
    return dst


if __name__ == "__main__":
    for p in sys.argv[1:]:
        try:
            fetch(p)
        except Exception as e:  # noqa: BLE001
            print(f"[{p}] FAILED: {e}", flush=True)
