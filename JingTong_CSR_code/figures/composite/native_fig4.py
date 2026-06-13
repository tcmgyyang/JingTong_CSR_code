# Figure 4 — native, recomputed from the expression matrix.
# A AUC lollipop (cervical degeneration), B AUC lollipop (CSR), C hub-gene expression boxplots, D ROC curves.
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["OMP_NUM_THREADS"] = "1"
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
ML = os.path.join(ROOT, "tables", "ML")
OUT = os.path.join(ROOT, "figures", "composite", "Fig4_diagnostics.png")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42})
OK = {"CON": "#56B4E9", "CSR": "#E69F00", "DCM": "#D55E00"}
BLUE, VERM, GREEN, PURPLE = "#0072B2", "#D55E00", "#009E73", "#CC79A7"

expr = pd.read_csv(os.path.join(ROOT, "data", "processed", "expr_symbol.csv.gz"), index_col=0, compression="gzip")
lab = pd.read_csv(os.path.join(ML, "labels.csv"), index_col=0).iloc[:, 0]
expr = expr[lab.index]
cerv = lab.isin(["CSR", "DCM"]).values        # positive class for cervical degeneration
csr = (lab == "CSR").values
con = (lab == "CON").values

# ---- per-gene AUC/CI tables ----
def read_blood_auc():
    rows = list(csv.reader(open(os.path.join(ML, "blood_hub_AUC.csv"), encoding="utf-8-sig")))
    out = {}
    for r in rows[3:]:
        if r and r[0]:
            out[r[0]] = {"cerv": (float(r[1]), float(r[2]), float(r[3])),
                         "csr": (float(r[4]), float(r[5]), float(r[6]))}
    return out

def read_roc(fn):
    return {r["feature"]: (float(r["AUC"]), float(r["CI_low"]), float(r["CI_high"]))
            for r in csv.DictReader(open(os.path.join(ML, fn), encoding="utf-8-sig"))}

bh = read_blood_auc()
roc_c = read_roc("ROC_CERVvsCON.csv"); roc_s = read_roc("ROC_CSRvsCON.csv")
auc_cerv = {g: bh[g]["cerv"] for g in bh}; auc_cerv["TP53"] = roc_c["TP53"]
auc_csr = {g: bh[g]["csr"] for g in bh};  auc_csr["TP53"] = roc_s["TP53"]

def lollipop(ax, aucs, title, L):
    items = sorted(aucs.items(), key=lambda kv: kv[1][0])
    g = [k for k, _ in items]; a = [v[0] for _, v in items]
    lo = [v[1] for _, v in items]; hi = [v[2] for _, v in items]
    y = np.arange(len(g))
    ax.hlines(y, lo, hi, color="#bbbbbb", lw=2, zorder=1)
    cols = [VERM if av >= 0.7 else BLUE for av in a]
    ax.scatter(a, y, color=cols, s=44, zorder=3)
    ax.axvline(0.5, color="0.6", ls="--", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(g, fontsize=8)
    ax.set_xlim(0.2, 1.02); ax.set_xlabel("AUC (95% CI)", fontsize=8.5)
    ax.set_title(title, fontsize=9.5); ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.22, 1.04, L, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")

def manual_roc(scores, pos):
    if np.corrcoef(scores, pos.astype(float))[0, 1] < 0:
        scores = -scores
    order = np.argsort(-scores); s = pos[order].astype(float)
    tps = np.cumsum(s); fps = np.cumsum(1 - s); P = s.sum(); N = len(s) - P
    tpr = np.concatenate([[0], tps / P]); fpr = np.concatenate([[0], fps / N])
    return fpr, tpr

fig = plt.figure(figsize=(12, 9.6))
gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.30)
lollipop(fig.add_subplot(gs[0, 0]), auc_cerv, "Cervical degeneration vs control", "A")
lollipop(fig.add_subplot(gs[0, 1]), auc_csr, "CSR vs control", "B")

# ---- C: expression boxplots ----
axC = fig.add_subplot(gs[1, 0])
genes = ["MMP9", "AKT1", "IL1B", "TP53"]; grp_order = ["CON", "CSR", "DCM"]
width = 0.24
for gi, g in enumerate(genes):
    vals = np.log2(expr.loc[g].values.astype(float) + 1)
    for ki, k in enumerate(grp_order):
        d = vals[(lab == k).values]
        bp = axC.boxplot([d], positions=[gi + (ki - 1) * width], widths=width * 0.85,
                         patch_artist=True, showfliers=False, manage_ticks=False)
        for b in bp["boxes"]:
            b.set(facecolor=OK[k], alpha=0.85, edgecolor="black", linewidth=0.5)
        for med in bp["medians"]:
            med.set(color="black", linewidth=0.8)
axC.set_xticks(range(len(genes))); axC.set_xticklabels(genes, fontsize=8.5)
axC.set_ylabel("log2(FPKM + 1)", fontsize=8.5); axC.set_title("Hub-gene expression by group", fontsize=9.5)
axC.spines[["top", "right"]].set_visible(False)
from matplotlib.patches import Patch
axC.legend(handles=[Patch(facecolor=OK[k], label=k) for k in grp_order], fontsize=7.5, loc="upper right")
axC.text(-0.18, 1.04, "C", transform=axC.transAxes, fontsize=16, fontweight="bold", va="bottom")

# ---- D: ROC curves (cervical degeneration) ----
axD = fig.add_subplot(gs[1, 1])
rocgenes = [("IL1B", BLUE), ("MMP9", VERM), ("AKT1", GREEN), ("TP53", PURPLE)]
for g, c in rocgenes:
    fpr, tpr = manual_roc(expr.loc[g].values.astype(float), cerv)
    a = roc_c.get(g, (np.nan,))[0]
    axD.plot(fpr, tpr, color=c, lw=1.8, label=f"{g}: AUC={a:.2f}")
axD.plot([0, 1], [0, 1], color="0.6", ls="--", lw=0.8)
axD.set_xlabel("1 − specificity", fontsize=8.5); axD.set_ylabel("Sensitivity", fontsize=8.5)
axD.set_title("ROC: cervical degeneration vs control", fontsize=9.5)
axD.legend(fontsize=7.5, loc="lower right"); axD.spines[["top", "right"]].set_visible(False)
axD.set_xlim(-0.02, 1.02); axD.set_ylim(-0.02, 1.02)
axD.text(-0.18, 1.04, "D", transform=axD.transAxes, fontsize=16, fontweight="bold", va="bottom")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT.replace(".png", f".{ext}"), dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig4_diagnostics (native, recomputed)")
