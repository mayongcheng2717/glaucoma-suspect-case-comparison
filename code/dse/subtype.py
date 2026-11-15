"""Cross-phenotype (subtype) contrast engine.

Design note
-----------
FinnGen gives every endpoint its own control set, so comparing two GWAS betas
(b1 vs b2) mixes a true subtype effect difference with a difference in the
control definition, and the shared samples make the two estimates correlated in
a way that a naive ``se = sqrt(se1^2 + se2^2)`` ignores.  We therefore test
subtype specificity *within cases*: allele frequencies are compared between two
case groups, which is a logistic model of subtype membership on genotype and is
invariant to the control definition.

For nested endpoints (child ⊂ parent) the allele frequency of the residual
group ``parent \\ child`` is obtained exactly from the two reported case
frequencies, so the same contrast applies without any approximation.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import finngen as fg
from . import stats

DATA = r"D:\eye100\data\tophits"


# ---------------------------------------------------------------------------
# locus catalogue
# ---------------------------------------------------------------------------
def load_tophits(pheno: str) -> pd.DataFrame:
    fp = os.path.join(DATA, f"{pheno}.tsv.gz")
    if not os.path.exists(fp):
        raise FileNotFoundError(f"run extract_top.py {pheno}")
    df = pd.read_csv(fp, sep="\t", compression="gzip")
    df["chrom"] = df["chrom"].astype(str)
    return df


def loci_from_tophits(df: pd.DataFrame, flank: int = 500_000,
                      pmax: float = 5e-8) -> List[dict]:
    """Cluster genome-wide significant variants into loci."""
    d = df[df["pval"] <= pmax].copy()
    if not len(d):
        return []
    loci = []
    for ch, sub in d.groupby("chrom"):
        sub = sub.sort_values("pos")
        cur = None
        for _, r in sub.iterrows():
            if cur is None or r["pos"] - cur["end"] > flank:
                if cur:
                    loci.append(cur)
                cur = {"chrom": ch, "start": int(r["pos"]), "end": int(r["pos"]),
                       "rows": [r]}
            else:
                cur["end"] = int(r["pos"])
                cur["rows"].append(r)
        if cur:
            loci.append(cur)
    out = []
    for l in loci:
        rows = pd.DataFrame(l["rows"])
        best = rows.loc[rows["pval"].idxmin()]
        out.append({
            "chrom": l["chrom"], "start": l["start"], "end": l["end"],
            "lead_pos": int(best["pos"]), "lead_ref": best["ref"], "lead_alt": best["alt"],
            "lead_rsid": best["rsids"], "lead_p": float(best["pval"]),
            "lead_beta": float(best["beta"]), "lead_gene": best["nearest_genes"],
            "mid": (l["start"] + l["end"]) // 2,
        })
    out.sort(key=lambda x: (x["chrom"], x["mid"]))
    for o in out:
        o["id"] = f"{o['chrom']}:{o['mid'] / 1e6:.2f}Mb"
    return out


def union_loci(*locus_lists: Sequence[dict], flank: int = 500_000) -> List[dict]:
    """Merge locus catalogues from several endpoints into independent regions."""
    allv = []
    for lst in locus_lists:
        for l in lst:
            allv.append(l)
    if not allv:
        return []
    bychr: Dict[str, List[dict]] = {}
    for l in allv:
        bychr.setdefault(l["chrom"], []).append(l)
    merged = []
    for ch, lst in bychr.items():
        lst.sort(key=lambda x: x["mid"])
        cur = None
        for l in lst:
            if cur is None or l["mid"] - cur["mid"] > flank:
                if cur:
                    merged.append(cur)
                cur = {"chrom": ch, "start": l["start"], "end": l["end"],
                       "mid": l["mid"], "members": [l]}
            else:
                cur["start"] = min(cur["start"], l["start"])
                cur["end"] = max(cur["end"], l["end"])
                cur["members"].append(l)
        if cur:
            merged.append(cur)
    merged.sort(key=lambda x: (x["chrom"], x["mid"]))
    for m in merged:
        m["id"] = f"{m['chrom']}:{m['mid'] / 1e6:.2f}Mb"
        m["genes"] = sorted({str(x.get("lead_gene")) for x in m["members"]
                             if x.get("lead_gene") not in (None, "NA", "nan")})
        m["best_p"] = min(x["lead_p"] for x in m["members"])
    return merged


# ---------------------------------------------------------------------------
# contrast computation
# ---------------------------------------------------------------------------
def _wide(a: pd.DataFrame, b: pd.DataFrame, na: str, nb: str) -> pd.DataFrame:
    m = a.merge(b, on=["chrom", "pos", "ref", "alt"], suffixes=(f"_{na}", f"_{nb}"))
    return m


def _add_contrast(m: pd.DataFrame, af_a: np.ndarray, n_a: int,
                  af_b: np.ndarray, n_b: int,
                  maf_min: float = 0.01, min_case_ac: float = 60.0) -> pd.DataFrame:
    """Attach case-only contrast columns for allele frequencies af_a vs af_b.

    Rare variants are excluded because the case-only variance ``1/(2 n p)`` is
    dominated by a handful of alleles and the log-odds ratio is unbounded as
    ``p -> 0``; such variants otherwise win any scan on p-value alone.
    """
    af_a = np.asarray(af_a, float)
    af_b = np.asarray(af_b, float)
    lo_af = np.minimum(af_a, 1 - af_a)
    lo_bf = np.minimum(af_b, 1 - af_b)
    ok = (np.isfinite(af_a) & np.isfinite(af_b)
          & (lo_af >= maf_min) & (lo_bf >= maf_min)
          & (2 * n_a * af_a >= min_case_ac) & (2 * n_b * af_b >= min_case_ac))
    m = m.loc[ok].reset_index(drop=True)
    af_a, af_b = af_a[ok], af_b[ok]
    lo = stats.logit(af_a) - stats.logit(af_b)
    # Var(logit p) = 1 / (2 n p (1-p)) for an allele frequency estimated from
    # 2n alleles; the (1-p) term is required and makes the test conservative
    # relative to dosage shrinkage in imputed data (empirical null sd 0.85).
    se = np.sqrt(1.0 / (2 * n_a * af_a * (1 - af_a))
                 + 1.0 / (2 * n_b * af_b * (1 - af_b)))
    z = lo / se
    if "af_alt" in m.columns:
        m["af_alt_total"] = m["af_alt"]
    m["af_case_a"] = af_a
    m["af_case_b"] = af_b
    m["maf_case"] = np.minimum(np.minimum(af_a, 1 - af_a), np.minimum(af_b, 1 - af_b))
    m["n_a"] = n_a
    m["n_b"] = n_b
    m["log_or_sub"] = lo
    m["se_sub"] = se
    m["or_sub"] = np.exp(lo)
    m["z_sub"] = z
    m["p_sub"] = stats.z_to_p(z)
    m["q_sub"] = stats.bh_fdr(m["p_sub"].to_numpy())
    return m


def contrast(a: pd.DataFrame, b: pd.DataFrame, n_a: int, n_b: int,
             na: str = "a", nb: str = "b") -> pd.DataFrame:
    """Case-only subtype contrast at every shared variant (disjoint case sets)."""
    if not len(a) or not len(b):
        return pd.DataFrame()
    m = _wide(a, b, na, nb)
    if not len(m):
        return m
    return _add_contrast(m, m[f"af_alt_cases_{na}"].to_numpy(float), n_a,
                         m[f"af_alt_cases_{nb}"].to_numpy(float), n_b)


def nested_residual(parent: pd.DataFrame, child: pd.DataFrame,
                    n_parent: int, n_child: int) -> pd.DataFrame:
    """Allele frequencies of the residual case group ``parent \\ child``."""
    m = _wide(parent, child, "p", "c")
    if not len(m):
        return m
    afp = m["af_alt_cases_p"].to_numpy(float)
    afc = m["af_alt_cases_c"].to_numpy(float)
    afr = (n_parent * afp - n_child * afc) / (n_parent - n_child)
    good = np.isfinite(afr) & (afr > 1e-6) & (afr < 1 - 1e-6)
    m = m.loc[good].reset_index(drop=True)
    m["af_case_res"] = afr[good]
    m["n_res"] = n_parent - n_child
    return m


def contrast_nested(parent: pd.DataFrame, child: pd.DataFrame,
                    n_parent: int, n_child: int) -> pd.DataFrame:
    """Compare the child case group against the residual of the parent group."""
    m = nested_residual(parent, child, n_parent, n_child)
    if not len(m):
        return m
    return _add_contrast(m, m["af_alt_cases_c"].to_numpy(float), n_child,
                         m["af_case_res"].to_numpy(float), n_parent - n_child)


def lead_of(df: pd.DataFrame, by: str = "p_sub") -> Optional[pd.Series]:
    if df is None or not len(df):
        return None
    return df.loc[df[by].idxmin()]


def calibrate(m: pd.DataFrame, p_a: str, p_b: str, null_p: float = 0.10,
              min_null: int = 150) -> dict:
    """Empirical null calibration of the case-only z statistic.

    The binomial variance of a case allele frequency does not include
    imputation error, which is not reported in the FinnGen summary statistics.
    Variants with no association in either subtype form an empirical null; the
    median of their z^2 gives the variance inflation factor lambda, and the
    contrast z is rescaled by sqrt(lambda).  This is the same logic as genomic
    control, applied to the contrast rather than to either GWAS.
    """
    lam = 1.0
    n_null = 0
    if len(m) and p_a in m.columns and p_b in m.columns:
        null = (m[p_a].to_numpy(float) > null_p) & (m[p_b].to_numpy(float) > null_p)
        n_null = int(null.sum())
        if n_null >= min_null:
            lam = float(np.median(m.loc[null, "z_sub"].to_numpy(float) ** 2) / 0.4549364)
    lam = max(lam, 1.0)
    m["z_sub_cal"] = m["z_sub"] / np.sqrt(lam)
    m["p_sub_cal"] = stats.z_to_p(m["z_sub_cal"])
    m["q_sub_cal"] = stats.bh_fdr(m["p_sub_cal"].to_numpy())
    return {"lambda": lam, "n_null": n_null}


def run_contrast(pheno_a: str, pheno_b: str, loci: Sequence[dict],
                 pad: int = 250_000, nested: bool = False,
                 threads: int = 4) -> Dict[str, pd.DataFrame]:
    """Compute the case-only subtype contrast locus by locus.

    ``nested=True`` treats ``pheno_a`` as a *subset* of ``pheno_b`` and compares
    the subset against the residual of the parent case group.
    """
    n_a = fg.n_cases(pheno_a)
    n_b = fg.n_cases(pheno_b)
    out = {}
    for l in loci:
        s, e = max(1, l["start"] - pad), l["end"] + pad
        fa = fg.region(pheno_a, l["chrom"], s, e)
        fb = fg.region(pheno_b, l["chrom"], s, e)
        if not len(fa) or not len(fb):
            continue
        if nested:
            d = contrast_nested(fb, fa, n_b, n_a)
        else:
            d = contrast(fa, fb, n_a, n_b, na=pheno_a, nb=pheno_b)
        if len(d):
            d["locus"] = l["id"]
            out[l["id"]] = d
    return out
