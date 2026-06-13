# Figure 3 — native, side-by-side: A GSVA pathway activity (left), B single-gene GSEA NES heatmap (right).
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["OMP_NUM_THREADS"] = "1"
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
E = os.path.join(ROOT, "tables", "enrichment")
OUT = os.path.join(ROOT, "figures", "composite", "Fig3_pathways.png")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42})
VERM, BLUE = "#D55E00", "#0072B2"

SHORT = {"NF-kappa B signaling pathway": "NF-κB", "MAPK signaling pathway": "MAPK",
         "AGE-RAGE signaling pathway in diabetic complications": "AGE-RAGE",
         "TNF signaling pathway": "TNF", "IL-17 signaling pathway": "IL-17",
         "PI3K-Akt signaling pathway": "PI3K-Akt", "Lipid and atherosclerosis": "Lipid/athero."}

fig = plt.figure(figsize=(13, 5.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.32)

# ---- A: GSVA lollipop ----
rows = list(csv.DictReader(open(os.path.join(E, "GSVA_key_pathways.csv"), encoding="utf-8-sig")))
rows.sort(key=lambda r: float(r["delta_CERVvsCON"]))
path = [r["pathway"] for r in rows]
delta = [float(r["delta_CERVvsCON"]) for r in rows]
pval = [float(r["p"]) for r in rows]
y = np.arange(len(path))
axA = fig.add_subplot(gs[0])
cols = [VERM if d >= 0 else BLUE for d in delta]
axA.hlines(y, 0, delta, color=cols, lw=2.2, zorder=2)
axA.scatter(delta, y, color=cols, s=42, zorder=3)
for i, (d, p) in enumerate(zip(delta, pval)):
    axA.text(d + (0.004 if d >= 0 else -0.004), i, f"p={p:.2f}",
             va="center", ha="left" if d >= 0 else "right", fontsize=6.3, color="#444")
axA.axvline(0, color="0.6", lw=0.8)
axA.set_yticks(y); axA.set_yticklabels([SHORT.get(p, p) for p in path], fontsize=8)
axA.set_xlabel("GSVA Δ (cervical degeneration − control)", fontsize=8.5)
axA.set_title("GSVA pathway activity (none p<0.05)", fontsize=9.5)
axA.set_xlim(min(delta) * 1.5 - 0.02, max(delta) * 1.5 + 0.02)
axA.spines[["top", "right"]].set_visible(False)
axA.text(-0.30, 1.05, "A", transform=axA.transAxes, fontsize=16, fontweight="bold", va="bottom")

# ---- B: single-gene GSEA NES heatmap ----
recs = list(csv.DictReader(open(os.path.join(E, "GSEA_key_pathways_summary.csv"), encoding="utf-8-sig")))
genes, terms = [], []
for r in recs:
    if r["gene"] not in genes: genes.append(r["gene"])
    if r["Term"] not in terms: terms.append(r["Term"])
NES = np.full((len(genes), len(terms)), np.nan); FDR = np.full_like(NES, np.nan)
for r in recs:
    i, j = genes.index(r["gene"]), terms.index(r["Term"])
    NES[i, j] = float(r["NES"]); FDR[i, j] = float(r["FDR"])
axB = fig.add_subplot(gs[1])
vmin, vmax = np.nanmin(NES), np.nanmax(NES)
norm = TwoSlopeNorm(vmin=min(vmin, 0.99), vcenter=1.0, vmax=max(vmax, 1.01))
im = axB.imshow(NES, cmap="RdBu_r", norm=norm, aspect="auto")
axB.set_xticks(range(len(terms))); axB.set_xticklabels([SHORT.get(t, t) for t in terms],
                                                       rotation=40, ha="right", fontsize=8)
axB.set_yticks(range(len(genes))); axB.set_yticklabels(genes, fontsize=8)
for i in range(len(genes)):
    for j in range(len(terms)):
        if not np.isnan(NES[i, j]):
            star = "*" if FDR[i, j] < 0.05 else ""
            axB.text(j, i, f"{NES[i,j]:.2f}{star}", ha="center", va="center", fontsize=5.6,
                     color="white" if abs(NES[i, j] - 1) > 0.55 else "black")
cb = fig.colorbar(im, ax=axB, fraction=0.045, pad=0.02); cb.set_label("NES (null = 1.0)", fontsize=8)
cb.ax.tick_params(labelsize=7)
axB.set_title("Single-gene GSEA (NES; * FDR<0.05)", fontsize=9.5)
axB.text(-0.16, 1.05, "B", transform=axB.transAxes, fontsize=16, fontweight="bold", va="bottom")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT.replace(".png", f".{ext}"), dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig3_pathways (native, side-by-side)")
