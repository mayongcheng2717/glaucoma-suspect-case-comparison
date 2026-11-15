"""Statistics for the eye-genetics pipeline.

Everything is implemented from first principles on summary statistics so that
the analyses do not depend on an R installation:

* z / p conversion, BH-FDR
* case-only subtype contrast (immune to control-set changes and to the
  sample-overlap error that inflates naive two-GWAS "effect difference" tests)
* residual allele frequency for a nested case definition
* Cochran's Q heterogeneity, IVW / Wald-ratio Mendelian randomisation
* single-causal-variant colocalisation (Wakefield approximate Bayes factors)
* approximate conditional analysis given an LD matrix
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats as st


# ---------------------------------------------------------------------------
# basic conversions
# ---------------------------------------------------------------------------
def p_to_z(p) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1.0)
    return st.norm.isf(p / 2.0)


def z_to_p(z) -> np.ndarray:
    return 2.0 * st.norm.sf(np.abs(np.asarray(z, dtype=float)))


def bh_fdr(p) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def logit(x):
    x = np.clip(np.asarray(x, dtype=float), 1e-12, 1 - 1e-12)
    return np.log(x / (1 - x))


def expit(x):
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=float)))


# ---------------------------------------------------------------------------
# case-only subtype contrast
# ---------------------------------------------------------------------------
def case_only_contrast(af_a, n_a, af_b, n_b) -> dict:
    """Compare allele frequencies between two *disjoint* case groups.

    Uses the 2x2 table of allele counts implied by (n, af) under HWE, i.e. a
    logistic model of subtype membership on genotype.  Because it conditions on
    case status it is unaffected by the control definition — the main reason
    "effect differences" between two FinnGen endpoints usually cannot be trusted.

    Returns OR (A vs B), its SE, z, p and the two allele counts.
    """
    af_a = float(np.clip(af_a, 1e-9, 1 - 1e-9))
    af_b = float(np.clip(af_b, 1e-9, 1 - 1e-9))
    ac_a = 2 * n_a * af_a
    ac_b = 2 * n_b * af_b
    lo = math.log(af_a / (1 - af_a)) - math.log(af_b / (1 - af_b))
    var = 1.0 / ac_a + 1.0 / (2 * n_a) + 1.0 / ac_b + 1.0 / (2 * n_b)
    # Hessian of the logistic likelihood gives 1/ac + 1/(2n) - (1/n) for the
    # logit of a binomial proportion; the -1/n terms cancel in the contrast.
    var = 1.0 / ac_a + 1.0 / ac_b
    se = math.sqrt(var)
    z = lo / se if se > 0 else np.nan
    return {"log_or": lo, "se": se, "or": math.exp(lo), "z": z,
            "p": float(2 * st.norm.sf(abs(z))) if se > 0 else np.nan,
            "ac_a": ac_a, "ac_b": ac_b, "n_a": n_a, "n_b": n_b,
            "af_a": af_a, "af_b": af_b}


def residual_case_af(af_parent, n_parent, af_child, n_child):
    """Allele frequency in ``parent \\ child`` for a nested case definition."""
    if n_parent <= n_child:
        return np.nan
    v = (n_parent * af_parent - n_child * af_child) / (n_parent - n_child)
    if not np.isfinite(v) or v < -0.02 or v > 1.02:
        return np.nan
    return float(np.clip(v, 1e-6, 1 - 1e-6))


# ---------------------------------------------------------------------------
# heterogeneity / effect-difference tests
# ---------------------------------------------------------------------------
def cochran_q(betas, ses) -> Tuple[float, float, int]:
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = b.size
    if k < 2:
        return np.nan, np.nan, k
    w = 1.0 / s ** 2
    mu = np.sum(w * b) / np.sum(w)
    q = float(np.sum(w * (b - mu) ** 2))
    return q, float(st.chi2.sf(q, k - 1)), k


def effect_difference(b1, s1, b2, s2, rho: float = 0.0) -> dict:
    """Naive two-estimate difference with an optional assumed correlation."""
    d = b1 - b2
    v = s1 ** 2 + s2 ** 2 - 2 * rho * s1 * s2
    se = math.sqrt(v) if v > 0 else np.nan
    z = d / se if se and np.isfinite(se) else np.nan
    return {"diff": d, "se": se, "z": z, "p": float(2 * st.norm.sf(abs(z))) if np.isfinite(z) else np.nan}


# ---------------------------------------------------------------------------
# Mendelian randomisation
# ---------------------------------------------------------------------------
def wald_ratio(bx, by, sy, sx: Optional[float] = None) -> dict:
    sx = 0.0 if sx is None else sx
    ratio = by / bx if bx != 0 else np.nan
    se = abs(sy / bx) if bx != 0 else np.nan
    z = ratio / se if se and np.isfinite(se) else np.nan
    return {"ratio": ratio, "se": se, "z": z,
            "p": float(2 * st.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
            "F": (bx / sx) ** 2 if sx and sx > 0 else np.nan}


def ivw(bx, by, sy, sx: Optional[np.ndarray] = None) -> dict:
    bx, by, sy = map(lambda a: np.asarray(a, float), (bx, by, sy))
    ok = np.isfinite(bx) & np.isfinite(by) & np.isfinite(sy) & (sy > 0) & (bx != 0)
    bx, by, sy = bx[ok], by[ok], sy[ok]
    k = bx.size
    if k == 0:
        return {"b": np.nan, "se": np.nan, "p": np.nan, "k": 0,
                "q": np.nan, "qp": np.nan, "F": np.nan, "i2": np.nan}
    ratio = by / bx
    w = (bx / sy) ** 2
    if k == 1:
        b = float(ratio[0])
        se = float(abs(sy[0] / bx[0]))
        q, qp, fs = np.nan, np.nan, (bx[0] / sx[0]) ** 2 if sx is not None and sx[0] > 0 else np.nan
        i2 = np.nan
    else:
        b = float(np.sum(w * ratio) / np.sum(w))
        se = float(math.sqrt(1.0 / np.sum(w)))
        q = float(np.sum(w * (ratio - b) ** 2))
        qp = float(st.chi2.sf(q, k - 1))
        i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
        fs = float(np.mean((bx / sx) ** 2)) if sx is not None else np.nan
    z = b / se if se else np.nan
    return {"b": b, "se": se, "z": z, "p": float(2 * st.norm.sf(abs(z))),
            "k": k, "q": q, "qp": qp, "F": fs, "i2": i2}


def mvmr(bx_exposures: np.ndarray, by, bz, sx=None, sz=None) -> dict:
    """Multivariable IVW for k exposures (weighted least squares)."""
    Bx = np.asarray(bx_exposures, float)          # k instruments x p exposures
    by = np.asarray(by, float)
    bz = np.asarray(bz, float)                    # outcome SE
    if Bx.ndim == 1:
        Bx = Bx[:, None]
    if sz is None:
        sz = np.full_like(by, np.nan)
    ok = np.isfinite(by) & (np.abs(by) > 0) & np.lexsort((np.arange(by.size),))[0] >= 0
    ok = np.all(np.isfinite(Bx), axis=1) & np.isfinite(by) & np.isfinite(bz) & (bz > 0)
    Bx, by, bz = Bx[ok], by[ok], bz[ok]
    k, p = Bx.shape
    if k <= p:
        return {"coef": [np.nan] * p, "se": [np.nan] * p, "p": [np.nan] * p, "k": k}
    W = np.diag(1.0 / bz ** 2)
    A = Bx.T @ W @ Bx
    try:
        Ainv = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return {"coef": [np.nan] * p, "se": [np.nan] * p, "p": [np.nan] * p, "k": k}
    coef = Ainv @ (Bx.T @ W @ by)
    resid = by - Bx @ coef
    dof = k - p
    sigma2 = float(resid @ W @ resid / dof) if dof > 0 else 1.0
    sigma2 = max(sigma2, 1.0)
    se = np.sqrt(np.diag(Ainv) * sigma2)
    z = coef / se
    return {"coef": coef, "se": se, "z": z, "p": 2 * st.norm.sf(np.abs(z)), "k": k,
            "q": float(resid @ W @ resid), "dof": dof}


# ---------------------------------------------------------------------------
# colocalisation (single causal variant)
# ---------------------------------------------------------------------------
def _logsum(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return -np.inf
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


def abf(z, v, w: float = 0.15 ** 2) -> np.ndarray:
    """Wakefield approximate Bayes factor for a single variant."""
    z = np.asarray(z, float)
    v = np.asarray(v, float)
    r = w / (v + w)
    return 0.5 * np.log(1 - r) + 0.5 * z ** 2 * r


def coloc_abf(b1, s1, b2, s2, p1=1e-4, p2=1e-4, p12=5e-6, w1=0.15 ** 2, w2=0.15 ** 2) -> dict:
    b1, s1, b2, s2 = map(lambda a: np.asarray(a, float), (b1, s1, b2, s2))
    ok = np.isfinite(b1) & np.isfinite(s1) & np.isfinite(b2) & np.isfinite(s2) & (s1 > 0) & (s2 > 0)
    if ok.sum() < 2:
        return {"PP.H4": np.nan, "PP.H3": np.nan, "PP.H1": np.nan, "PP.H2": np.nan,
                "PP.H0": np.nan, "nsnps": int(ok.sum())}
    z1 = b1[ok] / s1[ok]
    z2 = b2[ok] / s2[ok]
    l1 = abf(z1, s1[ok] ** 2, w1)
    l2 = abf(z2, s2[ok] ** 2, w2)
    h0 = 0.0
    h1 = _logsum(l1) + math.log(p1)
    h2 = _logsum(l2) + math.log(p2)
    h3 = _logsum(l1) + _logsum(l2) + math.log(p1 * p2)
    h4 = _logsum(l1 + l2) + math.log(p12)
    allh = np.array([h0, h1, h2, h3, h4], float)
    m = np.max(allh)
    pp = np.exp(allh - m)
    pp /= pp.sum()
    return {"PP.H4": pp[4], "PP.H3": pp[3], "PP.H1": pp[1], "PP.H2": pp[2],
            "PP.H0": pp[0], "nsnps": int(ok.sum()),
            "SNP.PP.H4_top": float(np.exp(l1 + l2 - np.max(l1 + l2)).max() /
                                   np.exp(l1 + l2 - np.max(l1 + l2)).sum())}


# ---------------------------------------------------------------------------
# approximate conditional analysis
# ---------------------------------------------------------------------------
def conditional_beta(beta, R, keep):
    """Approximate joint effect of SNP ``keep`` conditioning on others."""
    R = np.asarray(R, float)
    j = list(keep)
    o = [i for i in range(R.shape[0]) if i not in j]
    if not o:
        return np.asarray(beta, float)[j], None
    try:
        Roo_inv = np.linalg.inv(R[np.ix_(o, o)])
    except np.linalg.LinAlgError:
        Roo_inv = np.linalg.pinv(R[np.ix_(o, o)])
    return np.asarray(beta, float)[j] - R[np.ix_(j, o)] @ (Roo_inv @ np.asarray(beta, float)[o]), Roo_inv


def forward_select(beta, se, R, p_thresh=1e-5, max_steps=10):
    """Greedy approximate conditional analysis; returns selected indices."""
    beta = np.asarray(beta, float)
    se = np.asarray(se, float)
    chosen = []
    cur = beta.copy()
    Rc = np.asarray(R, float).copy()
    for _ in range(max_steps):
        z = np.abs(cur / se)
        order = np.argsort(-z)
        best = None
        for i in order:
            if i in chosen:
                continue
            best = i
            break
        if best is None or z[best] < p_to_z(0) * 0 + st.norm.isf(p_thresh / 2):
            break
        chosen.append(int(best))
        keep = [i for i in range(Rc.shape[0]) if i not in chosen]
        if not keep:
            break
        try:
            Roo_inv = np.linalg.inv(Rc[np.ix_(keep, keep)])
        except np.linalg.LinAlgError:
            Roo_inv = np.linalg.pinv(Rc[np.ix_(keep, keep)])
        cur[keep] = beta[keep] - Rc[np.ix_(keep, chosen)] @ (np.linalg.pinv(Rc[np.ix_(chosen, chosen)]) @ beta[chosen])
        cur[chosen] = beta[chosen]
    return chosen


def credible_set(z, r2, coverage=0.95, pp_prior=None):
    """Approximate credible set from marginal z and LD r^2 (per-SNP PP)."""
    z = np.asarray(z, float)
    bf = z ** 2
    k = bf.size
    pp = np.exp(bf / 2 - bf.max() / 2)
    pp /= pp.sum()
    if pp_prior is not None:
        pp = pp * np.asarray(pp_prior, float)
        pp /= pp.sum()
    order = np.argsort(-pp)
    cum = np.cumsum(pp[order])
    n = int(np.searchsorted(cum, coverage) + 1)
    return order[:n], pp, float(cum[min(n, k) - 1])
