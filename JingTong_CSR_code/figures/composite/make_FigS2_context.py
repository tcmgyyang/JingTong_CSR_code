# -*- coding: utf-8 -*-
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm

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

# Okabe-Ito palette
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERM = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
YELLOW = "#F0E442"
GREY = "#999999"
BLACK = "#000000"

BASE = r"H:\\毕业设计\\网药部分\\JingTong_CSR_paper"
TBL = os.path.join(BASE, "tables", "ML")
OUT_DIR = os.path.join(BASE, "figures", "composite")
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, "FigS2_context")


def panel_letter(ax, letter):
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes,
            fontsize=15, fontweight="bold", ha="left", va="bottom")


# ===================== load data =====================
# (A) tissue corroboration
A = pd.read_csv(os.path.join(TBL, "tissue_hubgene_GSE153761.csv"))
A.columns = [c.strip() for c in A.columns]
A = A.sort_values("log2FC_deg_vs_nor").reset_index(drop=True)

# (B) sample meta + age correlation
META = pd.read_csv(os.path.join(TBL, "sample_meta.csv"))
META.columns = [c.strip() for c in META.columns]
AGE = pd.read_csv(os.path.join(TBL, "hubgene_age_correlation.csv"))
AGE.columns = [c.strip() for c in AGE.columns]

# (C) hub x immune spearman r
C = pd.read_csv(os.path.join(TBL, "hub_immune_spearman_r.csv"), index_col=0)

# ===================== figure layout =====================
fig = plt.figure(figsize=(11, 9))
# 2 rows: top row A | B (equal heights), bottom row C spanning full width
gs = gridspec.GridSpec(
    2, 2, figure=fig,
    height_ratios=[1.0, 1.15],
    width_ratios=[1.0, 1.0],
    hspace=0.55, wspace=0.42,
    left=0.085, right=0.965, top=0.93, bottom=0.085,
)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, :])

# ---------------------------------------------------------------
# Panel A: diverging horizontal bar of log2FC (degenerate vs normal)
# ---------------------------------------------------------------
genes = A["hub_gene"].tolist()
lfc = A["log2FC_deg_vs_nor"].to_numpy(dtype=float)
pvals = A["p"].to_numpy(dtype=float)
ypos = np.arange(len(genes))

colors = [VERM if v > 0 else BLUE for v in lfc]
axA.barh(ypos, lfc, color=colors, edgecolor="black", linewidth=0.5, height=0.68)
axA.axvline(0, color=BLACK, linewidth=0.8)
axA.set_yticks(ypos)
axA.set_yticklabels(genes, fontsize=8)
axA.set_xlabel("log$_2$FC (degenerate vs normal)", fontsize=8)
axA.set_title("Tissue corroboration (GSE153761)", fontsize=9, pad=6)

# mark p<0.05 with "*"
xmax = np.nanmax(np.abs(lfc))
axA.set_xlim(-xmax * 1.30, xmax * 1.30)
for yi, v, p in zip(ypos, lfc, pvals):
    if p < 0.05:
        off = 0.06 * xmax
        if v >= 0:
            axA.text(v + off, yi, "*", ha="left", va="center",
                     fontsize=13, fontweight="bold", color=BLACK)
        else:
            axA.text(v - off, yi, "*", ha="right", va="center",
                     fontsize=13, fontweight="bold", color=BLACK)

# legend proxies
from matplotlib.patches import Patch
legA = [Patch(facecolor=VERM, edgecolor="black", label="Up in degenerate"),
        Patch(facecolor=BLUE, edgecolor="black", label="Down in degenerate")]
axA.legend(handles=legA, loc="lower right", fontsize=7, handlelength=1.2,
           borderpad=0.3, labelspacing=0.3)
axA.text(0.985, 0.04, "* p < 0.05 (TP53)", transform=axA.transAxes,
         ha="right", va="bottom", fontsize=7, color=GREY)
panel_letter(axA, "A")

# ---------------------------------------------------------------
# Panel B: age by group boxplot + strip, with inset lollipop of hub~age r
# ---------------------------------------------------------------
group_order = ["CON", "CSR", "DCM"]
group_colors = {"CON": GREEN, "CSR": ORANGE, "DCM": BLUE}
data_by_group = [META.loc[META["group"] == g, "age"].to_numpy(dtype=float)
                 for g in group_order]

bp = axB.boxplot(data_by_group, positions=range(len(group_order)),
                 widths=0.55, patch_artist=True, showfliers=False,
                 medianprops=dict(color=BLACK, linewidth=1.2),
                 whiskerprops=dict(color=BLACK, linewidth=0.8),
                 capprops=dict(color=BLACK, linewidth=0.8),
                 boxprops=dict(linewidth=0.8))
for patch, g in zip(bp["boxes"], group_order):
    patch.set_facecolor(group_colors[g])
    patch.set_alpha(0.35)
    patch.set_edgecolor(BLACK)

rng = np.random.default_rng(7)
for i, (g, vals) in enumerate(zip(group_order, data_by_group)):
    jitter = rng.uniform(-0.16, 0.16, size=len(vals))
    axB.scatter(i + jitter, vals, s=14, color=group_colors[g],
                edgecolor=BLACK, linewidth=0.35, alpha=0.9, zorder=3)

axB.set_xticks(range(len(group_order)))
axB.set_xticklabels([f"{g}\n(n={len(v)})" for g, v in zip(group_order, data_by_group)],
                    fontsize=8)
axB.set_ylabel("Age (years)", fontsize=8)
axB.set_title("Age confound across groups", fontsize=9, pad=6)
axB.set_ylim(18, 75)

# annotate non-overlapping control ages
con_max = np.nanmax(data_by_group[0])
axB.axhline(con_max, color=GREEN, linestyle="--", linewidth=0.8, alpha=0.7)
axB.text(0.02, con_max + 1.0, "Controls strictly younger (non-overlapping)",
         transform=axB.get_yaxis_transform(), fontsize=6.6, color=GREEN,
         ha="left", va="bottom")
panel_letter(axB, "B")

# inset: hub-gene ~ age Spearman r lollipop
axin = axB.inset_axes([0.56, 0.10, 0.40, 0.42])
ag = AGE.copy()
rcol = "Spearman r vs age"
ag = ag.sort_values(rcol).reset_index(drop=True)
ygenes = ag["hub_gene"].tolist()
rvals = ag[rcol].to_numpy(dtype=float)
yy = np.arange(len(ygenes))
axin.hlines(yy, 0, rvals, color=GREY, linewidth=1.0, zorder=1)
axin.scatter(rvals, yy, s=22, color=PURPLE, edgecolor=BLACK,
             linewidth=0.4, zorder=2)
axin.axvline(0, color=BLACK, linewidth=0.7)
axin.set_yticks(yy)
axin.set_yticklabels(ygenes, fontsize=6.2)
axin.set_xlim(-0.5, 0.5)
axin.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
axin.tick_params(axis="x", labelsize=6)
axin.tick_params(axis="y", labelsize=6.2)
axin.set_title("hub-gene ~ age Spearman r  (|r| ≤ 0.20)", fontsize=6.6, pad=3)
axin.axvspan(-0.20, 0.20, color=GREEN, alpha=0.10, zorder=0)
for s in ["top", "right"]:
    axin.spines[s].set_visible(False)

# ---------------------------------------------------------------
# Panel C: hub-gene x immune-cell Spearman r heatmap (RdBu_r, centred 0)
# ---------------------------------------------------------------
Cmat = C.copy()
# keep cell columns in original order; mark empty (all-NaN) columns
cell_cols = list(Cmat.columns)
M = Cmat.to_numpy(dtype=float)

vmax = np.nanmax(np.abs(M))
vmax = max(vmax, 0.6)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad(color="#E6E6E6")  # NaN / not-estimable cells in light grey

im = axC.pcolormesh(M, cmap=cmap, norm=norm, edgecolors="white", linewidth=0.6)

axC.set_xticks(np.arange(len(cell_cols)) + 0.5)
axC.set_xticklabels(cell_cols, rotation=45, ha="right", fontsize=6.8)
axC.set_yticks(np.arange(len(Cmat.index)) + 0.5)
axC.set_yticklabels(list(Cmat.index), fontsize=8)
axC.set_ylim(len(Cmat.index), 0)  # genes top-to-bottom
axC.set_title("Hub-gene vs immune-cell Spearman correlation", fontsize=9, pad=6)
axC.tick_params(length=0)
for s in ["top", "right", "left", "bottom"]:
    axC.spines[s].set_visible(False)

# annotate notable cells
def mark_cell(gene, cell, label):
    if gene in list(Cmat.index) and cell in cell_cols:
        gi = list(Cmat.index).index(gene)
        ci = cell_cols.index(cell)
        val = M[gi, ci]
        if not np.isnan(val):
            axC.text(ci + 0.5, gi + 0.5, "●", ha="center", va="center",
                     fontsize=6.5, color=BLACK)
            txtcol = "white" if abs(val) > 0.40 else BLACK
            axC.text(ci + 0.5, gi + 0.5, f"{val:+.2f}", ha="center", va="center",
                     fontsize=6.0, color=txtcol, fontweight="bold")

mark_cell("IL1B", "Neutrophils", "IL1B~Neutrophils")
mark_cell("MMP9", "Macrophages M2", "MMP9~Macrophages M2")

# textual callouts under the heatmap
axC.text(0.0, -0.52,
         "Notable: IL1B – Neutrophils r = +0.60   |   MMP9 – Macrophages M2 r = −0.52",
         transform=axC.transAxes, fontsize=7.2, color=BLACK, ha="left", va="top")
axC.text(1.0, -0.52,
         "grey = not estimable (constant fraction)",
         transform=axC.transAxes, fontsize=6.5, color=GREY, ha="right", va="top")

# colorbar
cbar = fig.colorbar(im, ax=axC, orientation="vertical", fraction=0.018,
                    pad=0.012, aspect=18)
cbar.set_label("Spearman r", fontsize=7.5)
cbar.ax.tick_params(labelsize=6.5)
cbar.outline.set_linewidth(0.6)

panel_letter(axC, "C")

# ===================== save =====================
import gc
# Compute the tight bbox once (with a small pad) so each save does a single
# render -- avoids holding two large dpi=300 RGBA buffers at the same time,
# which triggers MemoryError under the constrained BLAS/memory environment.
pad = 0.15  # inches
tb = fig.get_tightbbox(fig.canvas.get_renderer())
tb = tb.padded(pad)

# vector formats first (cheap on memory)
fig.savefig(OUT + ".pdf", bbox_inches=tb, facecolor="white")
gc.collect()
fig.savefig(OUT + ".svg", bbox_inches=tb, facecolor="white")
gc.collect()
# raster last
fig.savefig(OUT + ".png", dpi=300, bbox_inches=tb, facecolor="white")
gc.collect()

# report pixel size of png
from PIL import Image
with Image.open(OUT + ".png") as im_:
    print("PNG_SIZE", im_.size)
print("SAVED", OUT + ".png")
