# Disease-target UpSet plot (panel D of Figure 2): intersection of cervical-spondylosis
# disease-target databases. Pure-csv + matplotlib (low memory; no pandas/upsetplot).
import os, csv
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(os.path.dirname(__file__))
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.linewidth": 0.8})

# ---- load per-database disease genes ----
src = os.path.join(ROOT, "data", "source", "Diease-gene.csv")
rows = list(csv.DictReader(open(src, encoding="utf-8-sig")))
gene_src = {}
for r in rows:
    g = (r["Gene"] or "").strip()
    db = (r["Database"] or "").strip()
    if g and db:
        gene_src.setdefault(g, set()).add(db)

# canonical order + display labels (largest sets first)
SRC = ["GEO", "GeneCards", "Disgenet", "DrugBank", "TTD", "OMIM"]
LABEL = {"GEO": "GEO DEG", "GeneCards": "GeneCards", "Disgenet": "DisGeNET",
         "DrugBank": "DrugBank", "TTD": "TTD", "OMIM": "OMIM"}
setsize = {s: sum(1 for v in gene_src.values() if s in v) for s in SRC}

combo = Counter(frozenset(v) for v in gene_src.values())
items = sorted(combo.items(), key=lambda kv: (-kv[1]))[:13]   # top intersections
K = len(items)
nS = len(SRC)

BLUE = "#0072B2"; GREY_F = "#3b3b3b"; GREY_E = "#d9d9d9"; SETBAR = "#9ecae1"

fig = plt.figure(figsize=(8.4, 4.8))
gs = GridSpec(2, 2, width_ratios=[1.0, 4.4], height_ratios=[2.7, 2.0],
              hspace=0.06, wspace=0.34, figure=fig)
ax_bar = fig.add_subplot(gs[0, 1])
ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
ax_set = fig.add_subplot(gs[1, 0], sharey=ax_dot)

# --- top: intersection-size bars ---
sizes = [c for _, c in items]
ax_bar.bar(range(K), sizes, color=BLUE, width=0.62, edgecolor="none")
for i, v in enumerate(sizes):
    ax_bar.text(i, v + max(sizes) * 0.015, str(v), ha="center", va="bottom", fontsize=6.5)
ax_bar.set_ylabel("Intersection size", fontsize=9)
ax_bar.set_ylim(0, max(sizes) * 1.16)
ax_bar.spines[["top", "right"]].set_visible(False)
ax_bar.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_bar.tick_params(axis="y", labelsize=7.5)

# --- bottom-right: membership dot matrix ---
for c, (members, _) in enumerate(items):
    filled = [r for r, s in enumerate(SRC) if s in members]
    for r, s in enumerate(SRC):
        ax_dot.scatter(c, r, s=46, color=(GREY_F if s in members else GREY_E),
                       zorder=3, edgecolors="none")
    if len(filled) > 1:
        ax_dot.plot([c, c], [min(filled), max(filled)], color=GREY_F, lw=1.6, zorder=2)
# row shading for readability
for r in range(nS):
    if r % 2 == 0:
        ax_dot.axhspan(r - 0.5, r + 0.5, color="#f4f4f4", zorder=0)
ax_dot.set_xlim(-0.6, K - 0.4)
ax_dot.set_ylim(-0.6, nS - 0.4)
ax_dot.invert_yaxis()
ax_dot.set_yticks(range(nS))
ax_dot.set_yticklabels([LABEL[s] for s in SRC], fontsize=8.5, ha="right")
ax_dot.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
ax_dot.tick_params(axis="y", length=0, pad=3)
for sp in ax_dot.spines.values():
    sp.set_visible(False)

# --- bottom-left: per-database set sizes (bars grow leftward) ---
ax_set.barh(range(nS), [setsize[s] for s in SRC], color=SETBAR,
            edgecolor="#6baed6", height=0.6, zorder=3)
mx = max(setsize.values())
for r, s in enumerate(SRC):
    ax_set.text(setsize[s] + mx * 0.03, r, str(setsize[s]),
                va="center", ha="right", fontsize=6.5, color="#333333")
ax_set.set_xlim(mx * 1.15, -mx * 0.02)
ax_set.set_xlabel("Set size", fontsize=9)
ax_set.tick_params(axis="x", labelsize=7)
ax_set.tick_params(axis="y", which="both", left=False, labelleft=False)
ax_set.spines[["top", "left", "right"]].set_visible(False)

fig.suptitle("Cervical-spondylosis disease targets across databases (UpSet)",
             fontsize=10, y=0.99)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUT, f"Fig2_disease_upset.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig2_disease_upset  set sizes:", {LABEL[s]: setsize[s] for s in SRC},
      "| top combos:", [(sorted(m), c) for m, c in items[:6]])
