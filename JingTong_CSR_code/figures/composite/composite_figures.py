# Rebuild consolidated figures:
#  - uniform-width PIL montages (consistent panel sizes, tidy spacing, bold labels)
#  - Figure 6 re-rendered natively from data with a clean palette (no rainbow/jet)
# Run: py composite_figures.py
import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

FIG = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ROOT = os.path.abspath(os.path.join(FIG, ".."))
OUT = os.path.join(FIG, "composite")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.linewidth": 0.8,
                     "svg.fonttype": "none", "pdf.fonttype": 42})

OKABE = {"AKT1": "#0072B2", "IL6": "#E69F00", "IL1B": "#009E73",
         "TP53": "#CC79A7", "MMP9": "#D55E00"}
try:
    FONT_PATH = "C:/Windows/Fonts/arialbd.ttf"
    ImageFont.truetype(FONT_PATH, 40)
except Exception:
    FONT_PATH = None

# ---------------- PIL montage (uniform panel widths) ----------------
def _rw(im, w):
    return im.resize((w, max(1, round(im.height * w / im.width))), Image.LANCZOS)

def montage(rows, out, panel_w=980, gap=50, pad=50, lab_h=74, font_sz=60, bg=(255, 255, 255)):
    font = ImageFont.truetype(FONT_PATH, font_sz) if FONT_PATH else ImageFont.load_default()
    grid = [[(_rw(Image.open(os.path.join(FIG, p)).convert("RGB"), panel_w), L) for (p, L) in row]
            for row in rows]
    ncols = max(len(r) for r in grid)
    rowh = [max(im.height for im, _ in r) for r in grid]
    W = pad * 2 + ncols * panel_w + (ncols - 1) * gap
    H = pad * 2 + sum(rowh) + len(grid) * lab_h + (len(grid) - 1) * gap
    canvas = Image.new("RGB", (W, H), bg)
    d = ImageDraw.Draw(canvas)
    y = pad
    for ri, row in enumerate(grid):
        nr = len(row)
        roww = nr * panel_w + (nr - 1) * gap
        x = (W - roww) // 2
        for im, L in row:
            if L:
                d.text((x, y), L, font=font, fill=(0, 0, 0))
            yy = y + lab_h + (rowh[ri] - im.height) // 2
            canvas.paste(im, (x, yy))
            x += panel_w + gap
        y += lab_h + rowh[ri] + gap
    sz = canvas.size
    canvas.save(out, dpi=(170, 170))
    canvas.close()
    print("montage", os.path.basename(out), sz)

# ---------------- Figure 6 (native: docking matrix + MM-GBSA + per-residue) ----------------
def fig6():
    # docking matrix
    with open(os.path.join(ROOT, "tables", "docking_matrix", "docking_affinity_matrix.csv"),
              encoding="utf-8-sig") as f:
        rd = list(csv.reader(f))
    cols = rd[0][1:]                                   # AKT1,IL1B,IL6,MMP9,TP53
    names, vals = [], []
    for r in rd[1:]:
        row = []
        ok = False
        for v in r[1:]:
            if v.strip() == "":
                row.append(np.nan)
            else:
                row.append(float(v)); ok = True
        if ok:
            names.append(r[0]); vals.append(row)
    M = np.array(vals)
    order_c = [cols.index(c) for c in ["MMP9", "AKT1", "IL6", "IL1B", "TP53"]]
    M = M[:, order_c]; cols = ["MMP9", "AKT1", "IL6", "IL1B", "TP53"]
    strongest = np.nanmin(M, axis=1)
    ridx = np.argsort(strongest)                        # strongest binders first
    M = M[ridx]; names = [names[i] for i in ridx]
    def pretty(n):
        n = n.replace("_", " ").strip()
        return (n[:22] + "…") if len(n) > 23 else n
    names = [pretty(n) for n in names]

    # MM-GBSA
    with open(os.path.join(ROOT, "tables", "docking_MD", "MD_MMGBSA_summary.csv"),
              encoding="utf-8-sig") as f:
        md = list(csv.DictReader(f))
    md.sort(key=lambda r: float(r["mean_dG_GB_kcal_mol"]))   # most negative first
    cplx = [r["complex"] for r in md]
    dG = [float(r["mean_dG_GB_kcal_mol"]) for r in md]
    sd = [float(r["SD"]) for r in md]
    bcol = [OKABE[c.split("-")[0]] for c in cplx]

    # per-residue AKT1-quercetin
    with open(os.path.join(ROOT, "tables", "docking_MD", "MD_binding_residues.csv"),
              encoding="utf-8-sig") as f:
        pr = [r for r in csv.DictReader(f) if r["complex"] == "AKT1-quercetin"]
    pr = sorted(pr, key=lambda r: int(r["rank"]))[:8]
    res = [r["residue"].replace(":", "") for r in pr]
    rE = [float(r["mean_TOTAL_kcal_mol"]) for r in pr]

    fig = plt.figure(figsize=(11, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.18], height_ratios=[1, 1])
    # A: matrix
    from matplotlib.colors import TwoSlopeNorm
    axA = fig.add_subplot(gs[:, 0])
    finite = M[~np.isnan(M)]; vcen = -7.0
    norm = TwoSlopeNorm(vmin=min(finite.min(), vcen - 0.1), vcenter=vcen, vmax=max(finite.max(), vcen + 0.1))
    im = axA.imshow(M, aspect="auto", cmap="RdBu", norm=norm, interpolation="nearest")
    axA.set_xticks(range(len(cols))); axA.set_xticklabels(cols, fontsize=9, fontweight="bold")
    axA.set_yticks(range(len(names))); axA.set_yticklabels(names, fontsize=5.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isnan(v):
                axA.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=4.2,
                         color="white" if abs(v - vcen) > 2.0 else "black")
    axA.set_title("Docking affinity matrix (kcal/mol; red = stronger)", fontsize=10.5, pad=6)
    cb = fig.colorbar(im, ax=axA, fraction=0.045, pad=0.02); cb.ax.tick_params(labelsize=7)
    cb.set_label("Affinity (kcal/mol)", fontsize=8.5)
    axA.text(-0.18, 1.01, "A", transform=axA.transAxes, fontsize=22, fontweight="bold", va="bottom")
    # B: MM-GBSA
    axB = fig.add_subplot(gs[0, 1])
    yb = np.arange(len(cplx))
    axB.barh(yb, dG, xerr=sd, color=bcol, edgecolor="black", linewidth=0.6,
             error_kw=dict(ecolor="0.3", capsize=3, lw=1))
    axB.set_yticks(yb); axB.set_yticklabels(cplx, fontsize=9); axB.invert_yaxis()
    axB.set_xlabel("MM-GBSA ΔG$_{bind}$ (kcal/mol)", fontsize=10)
    axB.set_title("Binding free energy (mean ± SD, n=3)", fontsize=11, pad=6)
    axB.axvline(0, color="0.6", lw=0.8); axB.grid(axis="x", ls=":", alpha=0.5)
    from matplotlib.patches import Patch
    axB.legend(handles=[Patch(facecolor=OKABE[t], edgecolor="k", label=t) for t in ["AKT1", "IL6", "IL1B"]],
               fontsize=8, loc="lower left", frameon=False)
    axB.text(-0.16, 1.04, "B", transform=axB.transAxes, fontsize=22, fontweight="bold", va="bottom")
    # C: per-residue AKT1-quercetin
    axC = fig.add_subplot(gs[1, 1])
    yc = np.arange(len(res))
    axC.barh(yc, rE, color=OKABE["AKT1"], edgecolor="black", linewidth=0.6)
    axC.set_yticks(yc); axC.set_yticklabels(res, fontsize=9); axC.invert_yaxis()
    axC.set_xlabel("Per-residue ΔG (kcal/mol)", fontsize=10)
    axC.set_title("AKT1–quercetin ATP-pocket hot-spot residues", fontsize=11, pad=6)
    axC.grid(axis="x", ls=":", alpha=0.5)
    axC.text(-0.16, 1.04, "C", transform=axC.transAxes, fontsize=22, fontweight="bold", va="bottom")
    fig.savefig(os.path.join(OUT, "Fig6_docking_md.png"), dpi=150, facecolor="white")
    plt.close(fig); print("native Fig6_docking_md.png")

# ---------------- build all ----------------
fig6()

# NOTE: Fig 2/3/4/7 and S1/S2/S4 are now drawn NATIVELY by native_fig2/3/4/5.py,
# native_supp.py and rebuild_fig7.py. Their old montage calls were removed so this
# script only (re)builds Fig 6 (native) and the MD-trajectory supplementary (S5).

montage([[("MD/fig_rmsd_backbone.png", "A"), ("MD/fig_rmsd_ligand.png", "B")],
         [("MD/fig_rmsf.png", "C"), ("MD/fig_rgyr.png", "D")],
         [("MD/fig_sasa.png", "E"), ("MD/fig_hbond.png", "F")]],
        os.path.join(OUT, "FigS5_md.png"))

print("ALL DONE ->", OUT)
