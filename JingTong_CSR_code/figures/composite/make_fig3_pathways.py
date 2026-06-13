import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm

# ----------------------------------------------------------------------
# Shared style
# ----------------------------------------------------------------------
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
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
YELLOW = "#F0E442"
GREY = "#999999"
BLACK = "#000000"

BASE = "H:/毕业设计/网药部分/JingTong_CSR_paper"
GSVA_CSV = BASE + "/tables/enrichment/GSVA_key_pathways.csv"
GSEA_CSV = BASE + "/tables/enrichment/GSEA_key_pathways_summary.csv"
OUT = BASE + "/figures/composite/Fig3_pathways"

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
gsva = pd.read_csv(GSVA_CSV)
gsea = pd.read_csv(GSEA_CSV)

# ----------------------------------------------------------------------
# Figure layout: 1 row x 2 cols
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(12, 5))
gs = GridSpec(1, 2, figure=fig, wspace=0.55, width_ratios=[1.0, 1.05],
              left=0.07, right=0.94, top=0.90, bottom=0.12)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# ----------------------------------------------------------------------
# Panel A : GSVA pathway activity (horizontal lollipop)
# ----------------------------------------------------------------------
gA = gsva.sort_values("delta_CERVvsCON", ascending=True).reset_index(drop=True)
# ascending=True so that largest delta ends up at the TOP of horizontal axis
y = np.arange(len(gA))
delta = gA["delta_CERVvsCON"].values
pvals = gA["p"].values
labels = gA["pathway"].values

colors = [VERMILLION if d >= 0 else BLUE for d in delta]

# stems
for yi, di, ci in zip(y, delta, colors):
    axA.plot([0, di], [yi, yi], color=ci, lw=2.2, zorder=2, solid_capstyle="round")
# markers
axA.scatter(delta, y, color=colors, s=55, zorder=3, edgecolor="white", linewidth=0.6)

axA.axvline(0, color=BLACK, lw=0.8, zorder=1)

axA.set_yticks(y)
axA.set_yticklabels(labels, fontsize=7.5)
axA.set_xlabel("GSVA Δ (degeneration − control)", fontsize=8)

# annotate p-value at each point, placed on the side away from the bar
xpad = (delta.max() - delta.min()) * 0.04 if len(delta) else 0.01
xpad = max(xpad, 0.006)
for yi, di, pv in zip(y, delta, pvals):
    txt = f"p={pv:.3g}"
    if di >= 0:
        axA.text(di + xpad, yi, txt, va="center", ha="left", fontsize=6.5, color=GREY)
    else:
        axA.text(di - xpad, yi, txt, va="center", ha="right", fontsize=6.5, color=GREY)

# expand x-limits to give room for annotations
xmin = min(delta.min(), 0)
xmax = max(delta.max(), 0)
span = xmax - xmin
axA.set_xlim(xmin - span * 0.42, xmax + span * 0.42)
axA.set_ylim(-0.7, len(gA) - 0.3)
axA.tick_params(axis="both", length=3, width=0.8)
axA.set_title("GSVA pathway activity", fontsize=10, fontweight="bold", pad=8)
axA.text(-0.12, 1.06, "A", transform=axA.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")

# ----------------------------------------------------------------------
# Panel B : Single-gene GSEA NES heatmap
# ----------------------------------------------------------------------
nes_mat = gsea.pivot(index="gene", columns="Term", values="NES")
fdr_mat = gsea.pivot(index="gene", columns="Term", values="FDR")

# order genes by mean NES (descending -> strongest at top)
gene_order = nes_mat.mean(axis=1).sort_values(ascending=False).index
nes_mat = nes_mat.loc[gene_order]
fdr_mat = fdr_mat.loc[gene_order]

# order pathways by mean NES (descending columns left->right)
term_order = nes_mat.mean(axis=0).sort_values(ascending=False).index
nes_mat = nes_mat[term_order]
fdr_mat = fdr_mat[term_order]

data = nes_mat.values.astype(float)
fdr = fdr_mat.values.astype(float)

# diverging norm centred at 1.0 (GSEA null = NES 1)
vmin = np.nanmin(data)
vmax = np.nanmax(data)
# symmetric range around 1.0 for balanced diverging colours
half = max(1.0 - vmin, vmax - 1.0)
vmin_n = 1.0 - half
vmax_n = 1.0 + half
norm = TwoSlopeNorm(vmin=vmin_n, vcenter=1.0, vmax=vmax_n)

nrows, ncols = data.shape
# Use pcolormesh (flat shading) instead of imshow to avoid imshow's large
# internal resampling buffer, which overruns this constrained-memory env.
xedges = np.arange(ncols + 1)
yedges = np.arange(nrows + 1)
im = axB.pcolormesh(xedges, yedges, data, cmap="RdBu_r", norm=norm,
                    edgecolors="white", linewidth=1.0, shading="flat")
axB.set_aspect("auto")
axB.set_xlim(0, ncols)
axB.set_ylim(0, nrows)
axB.invert_yaxis()  # row 0 (first gene) at TOP, imshow-like orientation

# cell centres for ticks/annotations
xcent = np.arange(ncols) + 0.5
ycent = np.arange(nrows) + 0.5

axB.set_xticks(xcent)
axB.set_yticks(ycent)
axB.set_xticklabels(nes_mat.columns, fontsize=7, rotation=40, ha="right",
                    rotation_mode="anchor")
axB.set_yticklabels(nes_mat.index, fontsize=7.5)
axB.tick_params(axis="both", length=0)
for spine in axB.spines.values():
    spine.set_visible(False)

# annotate FDR<0.05 with *  (row i -> centre ycent[i])
for i in range(nrows):
    for j in range(ncols):
        if not np.isnan(fdr[i, j]) and fdr[i, j] < 0.05:
            axB.text(xcent[j], ycent[i], "*", ha="center", va="center",
                     fontsize=11, color=BLACK, fontweight="bold")

axB.set_title("Single-gene GSEA (NES)", fontsize=10, fontweight="bold", pad=8)
axB.text(-0.12, 1.06, "B", transform=axB.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")

# colorbar matched to panel height
cbar = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.02)
cbar.set_label("NES", fontsize=8)
cbar.ax.tick_params(labelsize=7, length=2)
cbar.outline.set_linewidth(0.6)
# mark NES=1 reference on colorbar
cbar.ax.axhline(norm(1.0), color=BLACK, lw=0.8)

# ----------------------------------------------------------------------
# Make both panels the SAME plotting height
# ----------------------------------------------------------------------
fig.canvas.draw()
posA = axA.get_position()
posB = axB.get_position()
# align top and bottom of both axes (use the common min top / max bottom)
top = min(posA.y1, posB.y1)
bottom = max(posA.y0, posB.y0)
axA.set_position([posA.x0, bottom, posA.width, top - bottom])
axB.set_position([posB.x0, bottom, posB.width, top - bottom])

# ----------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------
for ext in ("png", "pdf", "svg"):
    kw = {"dpi": 300} if ext == "png" else {}
    fig.savefig(f"{OUT}.{ext}", bbox_inches="tight", facecolor="white", **kw)

# report PNG pixel size
from PIL import Image
with Image.open(f"{OUT}.png") as im_:
    print("PNG_SIZE", im_.size)
print("DONE")
