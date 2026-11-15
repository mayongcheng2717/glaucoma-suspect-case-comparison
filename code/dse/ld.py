"""LD reference from 1000 Genomes phase 3, read remotely through tabix.

Only the requested 1-Mb bin is fetched; the EUR subset (GBR, FIN, CEU, TSI, IBS)
is cached as a compact int8 matrix so that LD look-ups are instant afterwards.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from .remoteio import RangedHTTP, Tabix

FT = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
CACHE = r"D:\eye100\cache\ld"
TBIC = r"D:\eye100\cache\tbi"
EUR_POPS = {"GBR", "FIN", "CEU", "TSI", "IBS"}

_panel: Optional[List[str]] = None
_samples: Optional[List[str]] = None
_eur_idx: Optional[np.ndarray] = None
_tabix: Dict[str, Tabix] = {}


def _chrom_file(chrom: str) -> str:
    chrom = str(chrom).replace("chr", "")
    return FT + f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"


def panel():
    """EUR sample names in VCF column order."""
    global _panel, _samples, _eur_idx
    if _panel is not None:
        return _panel, _eur_idx
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, "eur_samples.txt")
    if os.path.exists(fp):
        with open(fp) as f:
            _panel = f.read().split()
        return _panel, np.arange(len(_panel))
    # header
    tb = _get("1")
    data, _ = tb.bg.block(0)
    hdr = None
    for line in data.decode("utf-8", "replace").split("\n"):
        if line.startswith("#CHROM"):
            hdr = line
            break
    assert hdr, "could not read VCF header"
    all_samples = hdr.split("\t")[9:]
    rf = RangedHTTP(FT + "integrated_call_samples_v3.20130502.ALL.panel")
    raw = rf.read(0, rf.size or 200_000).decode("utf-8", "replace")
    pop = {}
    for line in raw.split("\n")[1:]:
        f = line.split("\t")
        if len(f) >= 2:
            pop[f[0]] = f[1]
    keep = [i for i, s in enumerate(all_samples) if pop.get(s) in EUR_POPS]
    _panel = [all_samples[i] for i in keep]
    _eur_idx = np.array(keep)
    with open(fp, "w") as f:
        f.write("\n".join(_panel))
    return _panel, _eur_idx


def _get(chrom: str) -> Tabix:
    if chrom not in _tabix:
        _tabix[chrom] = Tabix(_chrom_file(chrom), TBIC)
    return _tabix[chrom]


def _parse_gt(sub: bytes, n_sample: int) -> np.ndarray:
    need = n_sample * 4
    if len(sub) < need:                      # last field has no trailing tab
        sub = sub + b"\t" * (need - len(sub))
    arr = np.frombuffer(sub[:need], dtype=np.uint8).reshape(n_sample, 4)
    d = (arr[:, 0].astype(np.int16) - 48) + (arr[:, 2].astype(np.int16) - 48)
    d[(arr[:, 0] == ord(".")) | (arr[:, 2] == ord("."))] = -1
    d[(d < 0) | (d > 2)] = -1
    return d.astype(np.int8)


def bin_data(chrom: str, bin_start: int, bin_end: Optional[int] = None) -> Optional[dict]:
    """Genotypes for a 100-kb-aligned window (common SNPs only), cached on disk."""
    chrom = str(chrom).replace("chrom", "").replace("chr", "")
    step = 100_000
    lo = (int(bin_start) // step) * step
    hi = int(bin_end) if bin_end is not None else bin_start + 1_000_000
    hi = min(((int(hi) + step - 1) // step) * step, lo + 5_000_000)
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, f"kg3_{chrom}_{lo}_{hi}.npz")
    if os.path.exists(fp):
        z = np.load(fp, allow_pickle=True)
        return {"pos": z["pos"], "ref": z["ref"], "alt": z["alt"], "G": z["G"]}
    names, eidx = panel()
    tb = _get(chrom)
    pos, ref, alt, rows = [], [], [], []
    n_seen = 0
    for line in tb.query(chrom, lo, hi):
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if len(f) < 10:
            continue
        n_seen += 1
        info = f[7]
        a = f[4]
        if "," in a:                       # multi-allelic
            continue
        if "VT=SNP" not in info:            # SNP only
            continue
        try:
            af = float(info.split("AF=")[1].split(";")[0])
        except Exception:  # noqa: BLE001
            continue
        if af < 0.01 or af > 0.99:          # common variants only
            continue
        sub = f[9].encode()
        if len(sub) < 2504 * 4 - 1:
            continue
        d = _parse_gt(sub, 2504)
        if (d[eidx] < 0).mean() > 0.02:
            continue
        pos.append(int(f[1]))
        ref.append(f[3])
        alt.append(a)
        rows.append(d[eidx])
    print(f"  [1000G] {chrom}:{lo}-{hi} seen={n_seen} kept={len(pos)}", flush=True)
    if not pos:
        return None
    G = np.vstack(rows).astype(np.int8)
    out = {"pos": np.array(pos, np.int64), "ref": np.array(ref, dtype="U48"),
           "alt": np.array(alt, dtype="U48"), "G": G}
    np.savez_compressed(fp, **out)
    return out


def region_genotypes(chrom: str, start: int, end: int) -> Optional[dict]:
    """All 1000G EUR common SNPs in a window."""
    step = 100_000
    lo = (int(start) // step) * step
    hi = ((int(end) + step - 1) // step) * step
    parts = []
    for s in range(lo, hi, 1_000_000):
        d = bin_data(chrom, s, min(s + 1_000_000, hi))
        if d is None:
            continue
        m = (d["pos"] >= start) & (d["pos"] <= end)
        if m.any():
            parts.append({k: d[k][m] for k in d})
    if not parts:
        return None
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def ld_matrix(chrom: str, start: int, end: int,
              targets: List[Tuple[int, str, str]]) -> Tuple[np.ndarray, pd_index]:
    """LD (r) among ``targets`` [(pos, ref, alt)] using 1000G EUR."""
    d = region_genotypes(chrom, start, end)
    if d is None:
        return np.eye(len(targets)), list(range(len(targets)))
    idx = []
    for (p, r, a) in targets:
        hit = np.where((d["pos"] == p) & (d["ref"] == r) & (d["alt"] == a))[0]
        if hit.size == 0:
            hit = np.where((d["pos"] == p) & (d["ref"] == a) & (d["alt"] == r))[0]
        idx.append(int(hit[0]) if hit.size else -1)
    G = d["G"].astype(float)
    gmu = G.mean(1, keepdims=True)
    Gc = G - gmu
    sd = Gc.std(1, keepdims=True)
    sd[sd == 0] = np.nan
    Z = Gc / sd
    out = []
    for i in idx:
        row = []
        for j in idx:
            if i < 0 or j < 0:
                row.append(np.nan)
            else:
                row.append(float(np.nanmean(Z[i] * Z[j])))
        out.append(row)
    return np.array(out), idx


def r2_with(chrom: str, pos_lead: int, ref: str, alt: str,
            positions: List[int], start: int, end: int) -> np.ndarray:
    """r^2 of the lead variant with a list of positions."""
    d = region_genotypes(chrom, start, end)
    if d is None:
        return np.full(len(positions), np.nan)
    h = np.where((d["pos"] == pos_lead) & (d["ref"] == ref) & (d["alt"] == alt))[0]
    if h.size == 0:
        h = np.where(d["pos"] == pos_lead)[0]
    if h.size == 0:
        return np.full(len(positions), np.nan)
    G = d["G"].astype(float)
    Gc = G - G.mean(1, keepdims=True)
    sd = Gc.std(1)
    gl = Gc[h[0]]
    out = []
    for p in positions:
        j = np.where(d["pos"] == p)[0]
        if j.size == 0:
            out.append(np.nan)
            continue
        j = j[0]
        if sd[j] == 0 or sd[h[0]] == 0:
            out.append(np.nan)
            continue
        r = float(np.mean(gl * Gc[j]) / (sd[h[0]] * sd[j]))
        out.append(r * r)
    return np.array(out)


class pd_index(list):
    pass
