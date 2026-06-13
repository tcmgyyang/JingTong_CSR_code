import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------- shared style ----------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

# Okabe-Ito
VERMILLION = "#D55E00"
GREY = "#999999"
BLACK = "#000000"

CSV = r"H:/毕业设计/网药部分/JingTong_CSR_paper/tables/MR/MR_estimates.csv"
OUT_BASE = r"H:/毕业设计/网药部分/JingTong_CSR_paper/figures/composite/Fig5_MR"

df = pd.read_csv(CSV)

# Keep ONLY primary-method row per gene-outcome:
#   nsnp > 1  -> Inverse variance weighted
#   nsnp == 1 -> Wald ratio
def is_primary(row):
    if row["nsnp"] > 1:
        return row["method"] == "Inverse variance weighted"
    else:
        return row["method"] == "Wald ratio"

prim = df[df.apply(is_primary, axis=1)].copy()

outcomes = [
    ("cervical_disc", "Cervical disc disorders"),
    ("spondylosis",   "Spondylosis"),
    ("interverteb",   "Other IVD"),
]

genes = sorted(prim["gene"].unique())  # alphabetical
n_genes = len(genes)
# y positions: alphabetical TOP -> down (so first gene at top)
ypos = {g: n_genes - 1 - i for i, g in enumerate(genes)}

# Determine shared x-limits across all panels (log scale)
all_lo = prim["or_lci95"].min()
all_hi = prim["or_uci95"].max()
xmin = 10 ** np.floor(np.log10(all_lo) * 2) / 1.0  # placeholder
xmin = max(all_lo * 0.92, 1e-3)
xmax = all_hi * 1.08

fig = plt.figure(figsize=(9, 8), facecolor="white")
gs = fig.add_gridspec(
    nrows=1, ncols=3, wspace=0.12,
    left=0.13, right=0.985, top=0.90, bottom=0.10,
)

axes = []
for j, (okey, otitle) in enumerate(outcomes):
    ax = fig.add_subplot(gs[0, j])
    axes.append(ax)
    sub = prim[prim["outcome_name"] == okey]

    for _, r in sub.iterrows():
        g = r["gene"]
        y = ypos[g]
        sig = r["pval"] < 0.05
        col = VERMILLION if sig else GREY
        # CI line
        ax.plot([r["or_lci95"], r["or_uci95"]], [y, y],
                color=col, lw=1.2, solid_capstyle="round", zorder=2)
        # point: filled vermillion if sig, else filled grey
        ax.plot(r["or"], y, marker="o", ms=5.0,
                mfc=col, mec=BLACK if sig else col,
                mew=0.6 if sig else 0.0, zorder=3)

        # Annotate nominally significant ones with OR (95% CI)
        if sig:
            txt = f"{r['or']:.2f} ({r['or_lci95']:.2f}-{r['or_uci95']:.2f})"
            ax.annotate(
                txt, xy=(r["or_uci95"], y),
                xytext=(4, 0), textcoords="offset points",
                ha="left", va="center", fontsize=6.0,
                color=VERMILLION, fontweight="bold", zorder=4,
                clip_on=False,
            )

    # reference line OR = 1
    ax.axvline(1.0, color=BLACK, ls="--", lw=0.8, zorder=1)

    ax.set_xscale("log")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.6, n_genes - 0.4)

    # x ticks
    xticks = [0.5, 0.75, 1.0, 1.5, 2.0]
    xticks = [t for t in xticks if xmin <= t <= xmax]
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:g}" for t in xticks])
    ax.minorticks_off()
    ax.tick_params(axis="x", labelsize=7)

    ax.set_title(otitle, fontsize=9, pad=8)
    ax.set_xlabel("Odds ratio (95% CI)", fontsize=8)

    if j == 0:
        ax.set_yticks([ypos[g] for g in genes])
        ax.set_yticklabels(genes, fontsize=7.5)
        ax.set_ylabel("Gene", fontsize=8)
    else:
        ax.set_yticks([ypos[g] for g in genes])
        ax.set_yticklabels([])
        # keep left spine for uniformity but no labels
    ax.tick_params(axis="y", length=2)

    # faint horizontal guide rows
    for g in genes:
        ax.axhline(ypos[g], color="0.92", lw=0.5, zorder=0)

# panel letter
axes[0].text(-0.42, 1.045, "A", transform=axes[0].transAxes,
             fontsize=15, fontweight="bold", ha="left", va="bottom")

# Legend (shared, below)
legend_handles = [
    Line2D([0], [0], marker="o", color="none", mfc=VERMILLION, mec=BLACK,
           mew=0.6, ms=6, label="p < 0.05 (nominal)"),
    Line2D([0], [0], marker="o", color="none", mfc=GREY, mec=GREY,
           ms=6, label="p ≥ 0.05"),
    Line2D([0], [0], color=BLACK, ls="--", lw=0.8, label="OR = 1 (null)"),
]
fig.legend(handles=legend_handles, loc="lower center",
           ncol=3, fontsize=7.5, bbox_to_anchor=(0.55, 0.005),
           handletextpad=0.5, columnspacing=1.6)

fig.suptitle(
    "Mendelian randomization: causal effect of target-gene expression on IVD-related outcomes",
    fontsize=10, y=0.965,
)

# count significant
nsig = int((prim["pval"] < 0.05).sum())
print("primary rows:", len(prim))
print("nominally significant (p<0.05):", nsig)
print("genes:", n_genes, genes)

for ext, dpi in [("png", 300), ("pdf", None), ("svg", None)]:
    kw = dict(bbox_inches="tight", facecolor="white")
    if dpi:
        kw["dpi"] = dpi
    fig.savefig(f"{OUT_BASE}.{ext}", **kw)

plt.close(fig)

# report pixel dims of PNG
from PIL import Image
im = Image.open(f"{OUT_BASE}.png")
print("PNG pixel dimensions:", im.size)
