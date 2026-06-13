# -*- coding: utf-8 -*-
"""
make_mr_forest.py
Clean forest plot of two-sample MR causal estimates (hub-gene eQTL exposure ->
IVD / cervical-spondylosis outcomes). Reads the pre-computed estimate table only;
does NOT re-run run_MR.R or touch OpenGWAS.

Fixes vs. old Fig_MR_forest:
  - outcome labels were truncated to 'interverteb' and used ASCII '->'  -> mapped to
    clean human-readable names and a real Unicode arrow.
  - keep only the PRIMARY method per gene x outcome row (IVW where multi-SNP,
    Wald ratio where single-SNP).
  - LOG x-axis (OR is multiplicative) with sensible ticks + OR=1 reference.
  - points coloured by direction of effect (risk OR>1 / protective OR<1), Okabe-Ito.
  - right-hand text column 'OR (95% CI), P=...'.
  - rows grouped/sorted by gene; label 'GENE -> outcome'.

Run:  cd "H:/毕业设计/网药部分/expanded_analysis/MR" && py make_mr_forest.py
Requires: pandas, numpy, matplotlib. (pandas 2.x, matplotlib 3.8+ tested)
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

# --- shared house style (font, palette, panel, save) ---
# walk up from this file until figstyle.py is found (script lives 2 dirs below root)
_here = os.path.dirname(os.path.abspath(__file__))
_d = _here
for _ in range(5):
    if os.path.exists(os.path.join(_d, "figstyle.py")):
        sys.path.insert(0, _d)
        break
    _d = os.path.dirname(_d)
import figstyle as F  # noqa: E402

F.apply()

# ----------------------------------------------------------------------------
# 1. Data preparation
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "MR_results")
CSV = os.path.join(OUTDIR, "MR_estimates.csv")

df = pd.read_csv(CSV)

# keep only the primary causal-estimate method for each row
PRIMARY = {"Inverse variance weighted", "Wald ratio"}
df = df[df["method"].isin(PRIMARY)].copy()

# clean human-readable outcome names (fixes truncated 'interverteb')
OUTCOME_NAME = {
    "cervical_disc": "Cervical disc disorders",
    "spondylosis": "Spondylosis",
    "interverteb": "Other IVD disorders",
}
df["outcome_clean"] = df["outcome_name"].map(OUTCOME_NAME).fillna(df["outcome_name"])

# real arrow:  GENE -> outcome
ARROW = "→"  # U+2192 RIGHTWARDS ARROW
df["row_label"] = df["gene"].astype(str) + f"  {ARROW}  " + df["outcome_clean"]

# ----------------------------------------------------------------------------
# 2. Ordering: group by gene (sorted), and a fixed outcome order within each gene
# ----------------------------------------------------------------------------
OUTCOME_ORDER = ["Cervical disc disorders", "Spondylosis", "Other IVD disorders"]
df["outcome_clean"] = pd.Categorical(
    df["outcome_clean"], categories=OUTCOME_ORDER, ordered=True
)
# genes alphabetically; within gene, the fixed outcome order
df = df.sort_values(["gene", "outcome_clean"], ascending=[True, True]).reset_index(drop=True)

# plot bottom-to-top so the first gene (A...) sits at the TOP of the figure
df = df.iloc[::-1].reset_index(drop=True)
df["y"] = np.arange(len(df))

# ----------------------------------------------------------------------------
# 3. Colour by direction of effect (Okabe-Ito, matches F.PAL)
# ----------------------------------------------------------------------------
RISK = "#D55E00"        # OR > 1  -> risk-increasing (Okabe-Ito vermillion)
PROTECT = "#0072B2"     # OR < 1  -> protective (Okabe-Ito blue)
df["color"] = np.where(df["or"] >= 1.0, RISK, PROTECT)
# significance flag (primary-method p) for optional emphasis
df["sig"] = df["pval"] < 0.05

# right-hand annotation text: OR (95% CI), P=...
def fmt_p(p):
    if p < 0.001:
        return "P<0.001"
    return f"P={p:.3f}"

df["anno"] = df.apply(
    lambda r: f"{r['or']:.2f} ({r['or_lci95']:.2f}–{r['or_uci95']:.2f}), {fmt_p(r['pval'])}",
    axis=1,
)

# ----------------------------------------------------------------------------
# 4. Figure
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 6.5))
# leave room on the right for the OR/CI/P text column
fig.subplots_adjust(left=0.30, right=0.70, top=0.93, bottom=0.10)

XTICKS = [0.5, 0.7, 1.0, 1.4, 2.0]
XLIM = (0.45, 2.15)

# faint gene-band shading (alternate genes) for readability
genes_in_order = list(dict.fromkeys(df["gene"]))  # bottom-to-top order as plotted
for i, g in enumerate(genes_in_order):
    if i % 2 == 0:
        ys = df.loc[df["gene"] == g, "y"]
        ax.axhspan(ys.min() - 0.5, ys.max() + 0.5, color="#F2F2F2", zorder=0)

# OR = 1 reference line
ax.axvline(1.0, color="0.35", lw=0.9, ls="--", zorder=1)

# CI whiskers + points
for _, r in df.iterrows():
    ax.plot(
        [r["or_lci95"], r["or_uci95"]],
        [r["y"], r["y"]],
        color=r["color"],
        lw=1.3,
        solid_capstyle="round",
        zorder=2,
    )
# significant points filled & slightly larger, non-sig open
sig = df[df["sig"]]
nsig = df[~df["sig"]]
ax.scatter(sig["or"], sig["y"], s=34, c=sig["color"], edgecolors="white",
           linewidths=0.6, zorder=3)
ax.scatter(nsig["or"], nsig["y"], s=26, facecolors="white",
           edgecolors=nsig["color"], linewidths=1.2, zorder=3)

# ---- axes cosmetics ----
ax.set_xscale("log")
ax.set_xlim(*XLIM)
ax.set_xticks(XTICKS)
ax.set_xticklabels([f"{t:g}" for t in XTICKS])
# log scale auto-adds minor ticks/labels (e.g. 6x10^-1) -> suppress them
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_yticks(df["y"])
ax.set_yticklabels(df["row_label"])
ax.set_ylim(-0.6, len(df) - 0.4)
ax.set_xlabel("Odds ratio per SD increase in gene expression (95% CI)")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", length=3)

# ---- right-hand text column: OR (95% CI), P ----
x_text = 1.015  # axes-fraction x just right of the plot rectangle
ax.text(x_text, 1.035, "OR (95% CI), P value", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=7.5, fontweight="bold")
for _, r in df.iterrows():
    yfrac = (r["y"] - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.text(x_text, yfrac, r["anno"], transform=ax.transAxes,
            ha="left", va="center", fontsize=7,
            color="0.15", family="DejaVu Sans Mono")

# ---- direction legend (custom handles) ----
from matplotlib.lines import Line2D  # noqa: E402

legend_handles = [
    Line2D([0], [0], marker="o", color=RISK, markerfacecolor=RISK,
           markeredgecolor="white", markersize=6, lw=0, label="Risk (OR > 1)"),
    Line2D([0], [0], marker="o", color=PROTECT, markerfacecolor=PROTECT,
           markeredgecolor="white", markersize=6, lw=0, label="Protective (OR < 1)"),
    Line2D([0], [0], marker="o", color="0.4", markerfacecolor="white",
           markeredgecolor="0.4", markersize=6, lw=0, label="P ≥ 0.05 (open)"),
]
# legend below the x-axis label -> keeps the whole top strip free for the
# title (left) and the OR/CI/P column header (right), no collisions
ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.075),
          ncol=3, frameon=False, fontsize=7, handletextpad=0.4,
          columnspacing=1.4, borderaxespad=0.0)

# short title (no full sentence burned in), left-aligned over the plot rectangle
fig.suptitle("Mendelian randomization: hub-gene expression → IVD risk",
             x=0.30, y=1.0, ha="left", fontsize=9.5, fontweight="bold")

F.save(fig, os.path.join(OUTDIR, "Fig_MR_forest"))
print(f"Rows plotted: {len(df)}  |  genes: {genes_in_order}")
print("Saved:", os.path.join(OUTDIR, "Fig_MR_forest.png/.pdf"))
