"""All seven glaucoma figures, journal style: uppercase panels, no titles, correct QQ."""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st

sys.path.insert(0, r"D:\eye100\code")
from dse import plot as P  # noqa: E402

RES = r"D:\eye100\results\gwcas"
DEL = r"D:\eye100\deliverables\GLAUC_continuum"
FIG = os.path.join(DEL, "figures")
os.makedirs(FIG, exist_ok=True)

d = pd.read_csv(os.path.join(RES, "H7_GLAUCSUSP_vs_H7_GLAUCOMA_POAG_top200k.tsv.gz"),
                sep="\t")
d["mlp"] = -np.log10(np.clip(d["p"].to_numpy(float), 1e-300, None))
LEADS = json.load(open(os.path.join(RES, "H7_GLAUCSUSP_vs_H7_GLAUCOMA_POAG_leads.json")))
QQ = json.load(open(os.path.join(RES, "qq_stats.json")))

# ============================================================ Figure 1
fig = plt.figure(figsize=(180 / 25.4, 62 / 25.4))
ax = fig.add_axes([0.095, 0.18, 0.885, 0.71])
d2 = d.sort_values(["chrom", "pos"])
offsets, tick, pos = {}, [], 0
for c in sorted(d2["chrom"].unique()):
    m = d2["chrom"] == c
    span = d2.loc[m, "pos"].max() / 1e6 + 20
    offsets[c] = pos
    tick.append(pos + span / 2)
    pos += span
xd = np.array([p / 1e6 + offsets[c] for p, c in zip(d2["pos"], d2["chrom"])])
cols = [P.PAL["blue"] if int(c) % 2 else P.PAL["teal"] for c in d2["chrom"]]
ax.scatter(xd, d2["mlp"], s=1.8, c=cols, alpha=0.75, linewidths=0, rasterized=True)
ax.axhline(-np.log10(5e-8), ls="--", lw=0.55, c="k")
for L, nm in zip(LEADS["leads"][:2], ["9p21", "1q24"]):
    xx = L["pos"] / 1e6 + offsets.get(L["chrom"], 0)
    ax.scatter([xx], [-np.log10(max(L["p"], 1e-300))], marker="v", s=26,
               c=P.PAL["orange"], zorder=5, edgecolor="k", linewidth=0.3)
    ax.annotate(nm, xy=(xx, -np.log10(max(L["p"], 1e-300))), xytext=(0, 7),
                textcoords="offset points", fontsize=7, ha="center")
ax.set_ylim(0, d2["mlp"].max() * 1.30)
ax.set_xticks(tick)
ax.set_xticklabels([str(int(c)) for c in sorted(d2["chrom"].unique())], fontsize=5.2)
ax.set_xlim(xd.min(), xd.max())
ax.set_xlabel("Chromosome")
ax.set_ylabel(r"$-\log_{10}$ $P$")
P.save(fig, os.path.join(FIG, "Fig1_manhattan_suspect"))
print("Fig1")

# ============================================================ Figure 2
spec = pd.read_csv(os.path.join(RES, "lead_specificity.tsv"), sep="\t")
ORDER = ["WET_AMD", "DRY_AMD", "H7_GLAUCOMA_XFG", "H7_GLAUCSUSP",
         "H7_GLAUCOMA_NTG", "H7_GLAUCOMA_POAG"]
NICE = {"WET_AMD": "wet AMD", "DRY_AMD": "dry AMD", "H7_GLAUCOMA_XFG": "exfoliation",
        "H7_GLAUCSUSP": "suspect", "H7_GLAUCOMA_NTG": "normal tension",
        "H7_GLAUCOMA_POAG": "open angle"}
fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 76 / 25.4))
for ax, (lead, name, lab) in zip(axes, [("1q24_TMCO1", "TMCO1", "A"),
                                        ("9p21_CDKN2B", "CDKN2B", "B")]):
    sub = spec[spec["lead"] == lead].set_index("phenotype")
    xs = [NICE[p] for p in ORDER if p in sub.index]
    yv = [sub.loc[p, "af_cases"] for p in ORDER if p in sub.index]
    base = float(sub["af_all"].iloc[0])
    ax.axhline(base, ls="--", lw=0.7, c="k")
    ax.plot(range(len(xs)), yv, marker="o", ms=3.6, lw=0.9, color=P.PAL["blue"])
    for i, v in enumerate(yv):
        ax.annotate(f"{v:.3f}", (i, v), textcoords="offset points",
                    xytext=(0, -11 if i == 0 else 6), fontsize=5.4, ha="center")
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels(xs, rotation=32, fontsize=6.0, ha="right")
    ax.set_ylabel("Alt-allele frequency among cases")
    span = max(yv) - min(yv)
    ax.set_ylim(min(yv) - 0.25 * span, max(yv) + 0.30 * span)
    ax.text(0.03, 0.95, name, transform=ax.transAxes, fontsize=7.5, fontweight="bold",
            va="top")
    ax.text(0.03, 0.03, f"population {base:.3f}", transform=ax.transAxes, fontsize=6,
            style="italic")
    P.panel(ax, lab, x=-0.17, y=1.02)
plt.tight_layout()
P.save(fig, os.path.join(FIG, "Fig2_gradient"))
print("Fig2")

# ============================================================ Figure 3
z = np.load(os.path.join(RES, "qq_sample.npy")).astype(float)
fig = plt.figure(figsize=(178 / 25.4, 74 / 25.4))
ax1 = fig.add_axes([0.085, 0.18, 0.35, 0.70])
ax2 = fig.add_axes([0.60, 0.18, 0.35, 0.70])
qq = np.sort(np.abs(z))
theo = st.norm.ppf(0.5 + 0.5 * (np.arange(1, qq.size + 1) - 0.5) / qq.size)
ax1.scatter(theo, qq, s=1.2, c=P.PAL["grey"], alpha=0.45, linewidths=0, rasterized=True)
lim = max(theo.max(), qq.max()) * 1.03
ax1.plot([0, lim], [0, lim], ls="--", lw=0.6, c="k")
ax1.set_xlabel(r"Expected $|z|$")
ax1.set_ylabel(r"Observed $|z|$")
ax1.text(0.05, 0.93,
         f"$\\lambda$ = {QQ['lambda_raw']:.3f}\nSD(z) = {QQ['sd_z']:.3f}",
         transform=ax1.transAxes, fontsize=6.6, va="top")
P.panel(ax1, "A", x=-0.22, y=1.02)
counts, lams, names = [], [], []
for tag, nice in (("H7_GLAUCSUSP_vs_H7_GLAUCOMA_POAG", "suspect\nvs POAG"),
                  ("H7_GLAUCOMA_NTG_vs_H7_GLAUCOMA_POAG", "NTG\nvs POAG"),
                  ("H7_GLAUCOMA_XFG_vs_H7_GLAUCOMA_POAG", "XFG\nvs POAG"),
                  ("WET_AMD_vs_DRY_AMD", "wet vs\ndry AMD")):
    pp = os.path.join(RES, f"{tag}_leads.json")
    if os.path.exists(pp):
        j = json.load(open(pp))
        counts.append(j["n_independent"])
        lams.append(j.get("lambda_raw", j["lambda"]))
        names.append(nice)
ax2.bar(range(len(counts)), counts, color=P.PAL["blue"], width=0.55)
for i, (c, lm) in enumerate(zip(counts, lams)):
    ax2.text(i, c + 0.08, f"{c}", ha="center", fontsize=7)
    ax2.text(i, 0.15, f"$\\lambda$ {lm:.2f}", ha="center", fontsize=5.6, rotation=90,
             color="white")
ax2.set_xticks(range(len(names)))
ax2.set_xticklabels(names, fontsize=6.2)
ax2.set_ylabel("Independent loci at $P < 5\\times10^{-8}$")
ax2.set_ylim(0, max(counts) + 1.2)
P.panel(ax2, "B", x=-0.22, y=1.02)
P.save(fig, os.path.join(FIG, "Fig3_calibration"))
print("Fig3")

# ============================================================ Figure 4
fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 76 / 25.4))
for ax, (ch, pos_, name, lab) in zip(axes, [(9, 22_041_444, "CDKN2B", "A"),
                                            (1, 165_754_533, "TMCO1", "B")]):
    sub = d[(d["chrom"] == ch) & (np.abs(d["pos"] - pos_) < 400_000)].sort_values("pos")
    ax.scatter(sub["pos"] / 1e6, -np.log10(np.clip(sub["p"], 1e-300, None)),
               s=6, c=P.PAL["blue"], alpha=0.8, linewidths=0)
    ax.axvline(pos_ / 1e6, ls=":", lw=0.6, c=P.PAL["red"])
    ax.axhline(-np.log10(5e-8), ls="--", lw=0.5, c="k")
    ax.set_xlabel("Position (Mb, GRCh38)")
    ax.set_ylabel(r"$-\log_{10}$ $P$")
    ax.text(0.03, 0.94, name, transform=ax.transAxes, fontsize=7.5, fontweight="bold",
            va="top")
    P.panel(ax, lab, x=-0.18, y=1.02)
plt.tight_layout()
P.save(fig, os.path.join(FIG, "Fig4_regional"))
print("Fig4")

# ============================================================ Figure 5
eq = json.load(open(os.path.join(RES, "eqtl_anchoring.json")))
rows = []
for locus, rec in eq.items():
    for e in rec.get("eqtl", []):
        if e.get("lead_is_eqtl"):
            rows.append({"gene": e["gene"], "tissue": e["tissue"].replace("_", " "),
                         "beta": e["lead_beta"], "p": e["lead_p"]})
eqd = pd.DataFrame(rows)
fig, ax = P.style(108, 60)
ypos = np.arange(len(eqd))[::-1]
ax.barh(ypos, eqd["beta"], color=P.PAL["teal"], height=0.6)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{r.gene}, {r.tissue}" for r in eqd.itertuples()], fontsize=6.2)
for y, r in zip(ypos, eqd.itertuples()):
    ax.text(r.beta + 0.006, y, f"$P$ = {r.p:.0e}", va="center", fontsize=6)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("cis-eQTL effect size of the suspect-enriched allele")
ax.set_xlim(0, float(eqd["beta"].max()) * 1.45)
P.save(fig, os.path.join(FIG, "Fig5_eqtl"))
print("Fig5")

# ============================================================ Figure 6
fig, ax = P.style(94, 82)
ax.scatter(np.log2(d["or"]), d["mlp"], s=1.6, c=P.PAL["grey"], alpha=0.45,
           linewidths=0, rasterized=True)
sig = d["p"] < 5e-8
ax.scatter(np.log2(d.loc[sig, "or"]), d.loc[sig, "mlp"], s=2.6, c=P.PAL["red"],
           alpha=0.85, linewidths=0, rasterized=True)
ax.axhline(-np.log10(5e-8), ls="--", lw=0.6, c="k")
ax.set_xlabel("log2 odds ratio, suspect vs confirmed POAG")
ax.set_ylabel(r"$-\log_{10}$ $P$")
ax.set_xlim(-1.1, 1.1)
P.save(fig, os.path.join(FIG, "Fig6_volcano"))
print("Fig6")

# ============================================================ Figure 7
orr = pd.read_csv(os.path.join(RES, "lead_case_control.tsv"), sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 72 / 25.4), sharey=True)
NICE2 = {"WET_AMD": "wet AMD", "DRY_AMD": "dry AMD", "H7_GLAUCOMA_XFG": "exfoliation",
         "H7_GLAUCSUSP": "suspect", "H7_GLAUCOMA_NTG": "normal tension",
         "H7_GLAUCOMA_POAG": "open angle"}
for ax, (lead, name, lab) in zip(axes, [("1q24_TMCO1", "TMCO1", "A"),
                                        ("9p21_CDKN2B", "CDKN2B", "B")]):
    sub = orr[orr["lead"] == lead].set_index("phenotype")
    seq = [p for p in ORDER if p in sub.index][::-1]
    for i, p in enumerate(seq):
        r = sub.loc[p]
        col = P.PAL["grey"] if p in ("WET_AMD", "DRY_AMD") else P.PAL["blue"]
        ax.plot([r["lo"], r["hi"]], [i, i], color=col, lw=1.0, solid_capstyle="round")
        ax.scatter([r["or"]], [i], s=18, color=col, zorder=3, edgecolor="white",
                   linewidth=0.3)
    ax.axvline(1.0, color="k", lw=0.6, ls="--")
    ax.set_yticks(range(len(seq)))
    ax.set_yticklabels([NICE2[p] for p in seq], fontsize=6.4)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.text(0.97, 0.95, name, transform=ax.transAxes, fontsize=7.5, fontweight="bold",
            ha="right", va="top")
    if lead.startswith("1q"):
        ax.set_xlim(0.70, 1.10)
        ax.set_xticks([0.75, 0.85, 0.95, 1.05])
    else:
        ax.set_xlim(0.95, 1.72)
        ax.set_xticks([1.0, 1.1, 1.25, 1.5])
    ax.set_xticklabels([f"{v:g}" for v in ax.get_xticks()])
    ax.minorticks_off()
    P.panel(ax, lab, x=-0.20, y=1.02)
plt.tight_layout()
P.save(fig, os.path.join(FIG, "Fig7_forest"))
print("Fig7")

for f in sorted(os.listdir(FIG)):
    if f.endswith(".png"):
        print("  ", f, os.path.getsize(os.path.join(FIG, f)) // 1024, "KB")
