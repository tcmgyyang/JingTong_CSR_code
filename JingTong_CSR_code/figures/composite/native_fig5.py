# Figure 5 — Mendelian randomization forest, redrawn from raw MR estimates.
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["OMP_NUM_THREADS"] = "1"
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
OUT = os.path.join(ROOT, "figures", "composite", "Fig5_MR.png")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42})
RISK, PROT, NS = "#D55E00", "#0072B2", "#c4c4c4"

# primary-method estimate per (gene, outcome)
prim = {}
for r in csv.DictReader(open(os.path.join(ROOT, "tables", "MR", "MR_estimates.csv"), encoding="utf-8-sig")):
    nsnp = int(float(r["nsnp"])); meth = r["method"]
    want = "Inverse variance weighted" if nsnp > 1 else "Wald ratio"
    if meth != want:
        continue
    try:
        prim[(r["gene"], r["outcome_name"])] = (float(r["or"]), float(r["or_lci95"]),
                                                float(r["or_uci95"]), float(r["pval"]))
    except ValueError:
        pass

OUTC = [("cervical_disc", "Cervical disc disorders"), ("spondylosis", "Spondylosis"),
        ("interverteb", "Other IVD")]
genes = sorted({g for (g, _) in prim})
ny = len(genes)
fig, axes = plt.subplots(1, 3, figsize=(11.5, 7.4), sharey=True)

for ax, (ocode, otitle) in zip(axes, OUTC):
    for i in range(ny):                                 # alternating row bands
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color="#f6f6f6", zorder=0)
    ax.axvline(1.0, color="#888888", ls="--", lw=0.9, zorder=1)
    for i, g in enumerate(genes):
        v = prim.get((g, ocode))
        if v is None:
            continue
        orr, lo, hi, p = v
        sig = p < 0.05
        col = (RISK if orr >= 1 else PROT) if sig else NS
        ax.plot([lo, hi], [i, i], color=col, lw=2.2 if sig else 1.0, zorder=3,
                solid_capstyle="round")
        ax.scatter([orr], [i], s=58 if sig else 20, color=col, edgecolors="black" if sig else "none",
                   linewidths=0.6, zorder=4)
        if sig:
            ax.text(2.08, i, f"{orr:.2f} ({lo:.2f}–{hi:.2f})", va="center", ha="left",
                    fontsize=6.2, color=col, fontweight="bold")
    ax.set_xscale("log"); ax.set_xlim(0.5, 2.0)
    ax.set_xticks([0.5, 0.7, 1.0, 1.4, 2.0]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"], fontsize=7.5)
    ax.set_xlabel("Odds ratio (95% CI)", fontsize=8.5)
    ax.set_title(otitle, fontsize=10, fontweight="bold", pad=6)
    ax.set_ylim(-0.6, ny - 0.4); ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_yticks(range(ny)); axes[0].set_yticklabels(genes, fontsize=8)
axes[0].set_ylabel("Hub-gene exposure", fontsize=9)
fig.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=RISK, markeredgecolor="black", markersize=8, label="Risk-increasing (p<0.05)"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=PROT, markeredgecolor="black", markersize=8, label="Protective (p<0.05)"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=NS, markersize=6, label="Not significant")],
           loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Two-sample Mendelian randomization (eQTLGen → FinnGen): 17 genes × 3 outcomes; 7 nominally significant",
             fontsize=10.5, y=0.99)
fig.tight_layout(rect=[0, 0.04, 1, 0.96])
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT.replace(".png", f".{ext}"), dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig5_MR (native, polished forest)")
