"""FinnGen R13 access layer: endpoint registry, metadata, remote region queries.

All access goes through :mod:`dse.remoteio` (HTTP-range tabix), so a locus costs
a few hundred kilobytes instead of a 760 MB download.  Region pulls are cached
as gzipped TSVs so that re-running an analysis is free.
"""

from __future__ import annotations

import io
import json
import os
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

from .remoteio import Tabix

BASE = "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
ROOT = r"D:\eye100"
CACHE = os.path.join(ROOT, "cache")
TBI_CACHE = os.path.join(CACHE, "tbi")
REGION_CACHE = os.path.join(CACHE, "regions")

COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp",
        "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]

# ---------------------------------------------------------------------------
# Endpoint registry — only endpoints actually used by this project
# ---------------------------------------------------------------------------
ENDPOINTS: Dict[str, str] = {
    # ---- age-related macular degeneration ----
    "WET_AMD": "Wet (neovascular) age-related macular degeneration",
    "DRY_AMD": "Dry AMD, includes geographic atrophy",
    "H7_AMD": "AMD, dry or wet",
    "H7_MACULADEGEN": "Degeneration of macula and posterior pole",
    "MACULA_NAS": "Degenerative macular diseases",
    "MACULAR_CYST": "Macular cyst (DME proxy)",
    "MACULAR_HOLE": "Macular hole (non-exudative traction control)",
    "MACULAR_PUCKER": "Epiretinal membrane / macular pucker",
    # ---- diabetic retinopathy ----
    "DM_RETINOPATHY": "Diabetic retinopathy, all",
    "DM_RETINOPATHY_STRICT": "Diabetic retinopathy, strict definition",
    "DM_RETINOPATHY_EXMORE": "Diabetic retinopathy, extended control exclusions",
    "H7_PROLIFRETINOPATHOTH": "Other proliferative retinopathy",
    "H7_BICKRNDRETINOPAT": "Background retinopathy",
    "DM_NEOVASCULAR_GLAUCOMA": "Neovascular glaucoma",
    "T2D": "Type 2 diabetes",
    "DM": "Diabetes mellitus, broad",
    # ---- glaucoma ----
    "H7_GLAUCOMA_NTG": "Normal-tension glaucoma",
    "H7_GLAUCOMA_POAG": "Primary open-angle glaucoma, strict",
    "H7_GLAUCOMA": "Glaucoma, any",
    "H7_GLAUCOMAMEDICATION": "Use of antiglaucoma preparations (IOP-treatment proxy)",
    "H7_GLAUCSUSP": "Glaucoma suspect",
    "H7_GLAUCCLOSEPRIM": "Primary angle-closure glaucoma",
    "H7_GLAUCOMA_XFG": "Exfoliation glaucoma",
    # ---- myopia ----
    "H7_MYOPIA": "Myopia",
    "H7_DEGENMYOPIA": "Degenerative myopia",
    "H7_HYPERMETRO": "Hypermetropia",
    "H7_ASTIGMATISM": "Astigmatism",
    "H7_REFRACCOMMODIS": "Disorders of refraction and accommodation",
    # ---- retinal vascular occlusion ----
    "H7_RETIVASCOCCLUSION": "Retinal vascular occlusion (H34, artery+vein)",
    "H7_CENTRRETARTOCC": "Central retinal artery occlusion",
    "H7_OTHRETARTOCC": "Other retinal artery occlusion",
    "H7_RETINAVASC": "Retinal vascular disorders",
    # ---- comparators / safety ----
    "H7_CATARACTSENILE": "Senile cataract",
    "H7_RETINALDETACH": "Retinal detachment",
    "H7_VISIMPBILATERAL": "Blindness or visual impairment, binocular",
    "H7_OPTIC_ATROPHY": "Optic atrophy",
}

_meta_lock = threading.Lock()
_meta: Optional[dict] = None
_tabix: Dict[str, Tabix] = {}
_tabix_lock = threading.Lock()


def meta() -> dict:
    """Endpoint metadata (case/control counts) from the R13 API, cached on disk."""
    global _meta
    with _meta_lock:
        if _meta is not None:
            return _meta
        fp = os.path.join(CACHE, "finngen_r13_phenos.json")
        if not os.path.exists(fp):
            raise RuntimeError("run probe_finngen_meta.py first to cache endpoint metadata")
        with open(fp, encoding="utf-8") as f:
            _meta = json.load(f)
        return _meta


def n_cases(pheno: str) -> int:
    m = meta().get(pheno)
    if not m or not m.get("num_cases"):
        raise KeyError(f"no case count for {pheno}")
    return int(m["num_cases"])


def n_controls(pheno: str) -> int:
    return int(meta()[pheno]["num_controls"])


def label(pheno: str) -> str:
    return ENDPOINTS.get(pheno) or meta().get(pheno, {}).get("phenostring", pheno)


def tabix(pheno: str) -> Tabix:
    with _tabix_lock:
        if pheno not in _tabix:
            _tabix[pheno] = Tabix(BASE + f"finngen_R13_{pheno}.gz", TBI_CACHE)
        return _tabix[pheno]


def _to_frame(lines: Iterable[str]) -> pd.DataFrame:
    rows = [l.split("\t") for l in lines]
    rows = [r for r in rows if len(r) >= 13 and r[1].isdigit()]
    if not rows:
        return pd.DataFrame(columns=COLS).astype({"pos": "int64"})
    df = pd.DataFrame([r[:13] for r in rows], columns=COLS)
    for c in ("pos",):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["pos"].astype("int64")
    return df


def region(pheno: str, chrom: str, start: int, end: int, use_cache: bool = True) -> pd.DataFrame:
    """Variants for ``pheno`` in the 1-based inclusive window [start, end]."""
    chrom = str(chrom).replace("chr", "")
    key = f"{pheno}_{chrom}_{start}_{end}.tsv.gz"
    fp = os.path.join(REGION_CACHE, key)
    os.makedirs(REGION_CACHE, exist_ok=True)
    if use_cache and os.path.exists(fp):
        df = pd.read_csv(fp, sep="\t", compression="gzip")
        # the cache round-trip loses the string dtype of ``chrom``
        df["chrom"] = df["chrom"].astype(str)
        df["pos"] = df["pos"].astype("int64")
        return df
    lines = tabix(pheno).query(chrom, int(start) - 1, int(end))
    df = _to_frame(lines)
    if use_cache and len(df):
        df.to_csv(fp, sep="\t", index=False, compression="gzip")
    return df


def region_many(phenos: List[str], chrom: str, start: int, end: int,
                threads: int = 6) -> Dict[str, pd.DataFrame]:
    """Parallel region fetch across endpoints."""
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {p: ex.submit(region, p, chrom, start, end) for p in phenos}
        return {p: f.result() for p, f in futs.items()}


def variant(pheno: str, chrom: str, pos: int, pad: int = 1) -> pd.DataFrame:
    """Rows at one position (small padded window), returned exactly."""
    df = region(pheno, chrom, pos - pad, pos + pad)
    if not len(df):
        return df
    return df[df["pos"] == pos]


def variants_many(pheno: str, chrom: str, positions: List[int], pad: int = 1,
                  threads: int = 8) -> pd.DataFrame:
    """Fetch several positions of one endpoint efficiently (merged window)."""
    if not positions:
        return pd.DataFrame(columns=COLS)
    lo, hi = min(positions) - pad, max(positions) + pad
    span = hi - lo
    if span <= 400_000:
        df = region(pheno, chrom, lo, hi)
        return df[df["pos"].isin(positions)]
    frames = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = [ex.submit(region, pheno, chrom, p - pad, p + pad) for p in positions]
        for f in futs:
            d = f.result()
            if len(d):
                frames.append(d[d["pos"].isin(positions)])
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLS)


def all_endpoints() -> List[str]:
    return sorted(meta().keys())
