import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx
from itertools import combinations

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

ROOT = "H:/毕业设计/网药部分/JingTong_CSR_paper"

# ============================================================
# Generic UpSet drawing routine on a given subplot spec
# ============================================================
def draw_upset(fig, spec, membership, set_order, set_sizes, title, top_n=10,
               bar_color=BLUE, set_bar_color=GREY):
    """
    membership: dict combo(frozenset of set names) -> intersection size
    set_order : list of set display names (top->bottom in matrix)
    set_sizes : dict set name -> total size
    """
    n_sets = len(set_order)
    # nested gridspec: matrix occupies bottom, top bars on top, left set bars
    # layout columns: [left set-size bars | matrix/intersection columns]
    # layout rows:    [top intersection bars / matrix]
    inner = GridSpecFromSubplotSpec(
        2, 2, subplot_spec=spec,
        width_ratios=[0.9, 3.0], height_ratios=[2.0, 1.25],
        wspace=0.05, hspace=0.06)

    ax_top = fig.add_subplot(inner[0, 1])      # intersection bars
    ax_mat = fig.add_subplot(inner[1, 1], sharex=ax_top)  # dot matrix
    ax_set = fig.add_subplot(inner[1, 0], sharey=ax_mat)  # set-size bars

    # ---- pick top intersections ----
    combos = sorted(membership.items(), key=lambda kv: kv[1], reverse=True)
    combos = [(c, v) for c, v in combos if v > 0][:top_n]
    n_comb = len(combos)
    x = np.arange(n_comb)

    # row y positions: set_order[0] at top -> highest y
    y_pos = {name: (n_sets - 1 - i) for i, name in enumerate(set_order)}

    # ---- top intersection bar chart ----
    sizes = [v for _, v in combos]
    ax_top.bar(x, sizes, color=bar_color, width=0.6, zorder=3)
    for xi, s in zip(x, sizes):
        ax_top.text(xi, s + max(sizes) * 0.015, str(s), ha="center", va="bottom",
                    fontsize=6.3, color=BLACK)
    ax_top.set_ylim(0, max(sizes) * 1.18)
    ax_top.set_ylabel("Intersection\nsize", fontsize=7.5)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.tick_params(axis="y", labelsize=6.5)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.set_title(title, fontsize=9.5, pad=6)

    # ---- dot membership matrix ----
    ax_mat.set_xlim(-0.5, n_comb - 0.5)
    ax_mat.set_ylim(-0.5, n_sets - 0.5)
    # alternating row shading
    for name in set_order:
        yy = y_pos[name]
        if (n_sets - 1 - yy) % 2 == 0:
            ax_mat.axhspan(yy - 0.5, yy + 0.5, color="#f2f2f2", zorder=0)
    # light grey dots everywhere
    for xi, (combo, _) in enumerate(combos):
        members = set(combo)
        ys_member = [y_pos[s] for s in set_order if s in members]
        for name in set_order:
            yy = y_pos[name]
            is_mem = name in members
            ax_mat.scatter(xi, yy, s=58,
                           color=(VERM if is_mem else "#d9d9d9"),
                           edgecolors="none", zorder=3)
        # vertical link line between members
        if len(ys_member) > 1:
            ax_mat.plot([xi, xi], [min(ys_member), max(ys_member)],
                        color=VERM, lw=1.6, zorder=2)
    ax_mat.set_yticks([y_pos[n] for n in set_order])
    ax_mat.set_yticklabels(set_order, fontsize=7)
    ax_mat.set_xticks([])
    for sp in ["top", "right", "bottom", "left"]:
        ax_mat.spines[sp].set_visible(False)
    ax_mat.tick_params(axis="both", length=0)

    # ---- left set-size horizontal bars (grow leftward) ----
    set_vals = [set_sizes[n] for n in set_order]
    ys = [y_pos[n] for n in set_order]
    ax_set.barh(ys, set_vals, color=set_bar_color, height=0.55, zorder=3)
    ax_set.set_xlim(max(set_vals) * 1.05, 0)  # invert -> bars grow left
    ax_set.set_ylim(-0.5, n_sets - 0.5)
    ax_set.set_xlabel("Set size", fontsize=7.5)
    ax_set.tick_params(axis="x", labelsize=6.3)
    ax_set.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax_set.spines[sp].set_visible(False)
    # annotate set sizes at bar tips
    for yy, v in zip(ys, set_vals):
        ax_set.text(v + max(set_vals) * 0.04, yy, str(v), ha="left", va="center",
                    fontsize=6.2, color=BLACK)

    return ax_top  # for panel-letter placement


def build_membership_from_bool(df, set_cols):
    """df rows = genes, set_cols = list of boolean columns. Return membership dict."""
    membership = {}
    for _, row in df.iterrows():
        combo = frozenset(c for c in set_cols if bool(row[c]))
        if combo:
            membership[combo] = membership.get(combo, 0) + 1
    set_sizes = {c: int(df[c].sum()) for c in set_cols}
    return membership, set_sizes


def build_membership_from_long(df, db_col, item_col, db_order):
    """Long df with database/item rows. Return membership keyed by db-combos per item."""
    item_to_dbs = {}
    for _, row in df.iterrows():
        item_to_dbs.setdefault(row[item_col], set()).add(row[db_col])
    membership = {}
    for item, dbs in item_to_dbs.items():
        combo = frozenset(d for d in db_order if d in dbs)
        if combo:
            membership[combo] = membership.get(combo, 0) + 1
    set_sizes = {d: int((df[db_col] == d).sum()) for d in db_order}
    # use unique items per db for set size (dedupe)
    set_sizes = {d: df.loc[df[db_col] == d, item_col].nunique() for d in db_order}
    return membership, set_sizes


# ============================================================
# Load data
# ============================================================
# (A) drug targets 4 source
drug = pd.read_csv(f"{ROOT}/tables/targets/drug_targets_4source.csv")
drug_sets = ["TCMSP", "STITCH", "SwissTarget", "HERB"]
drug_mem, drug_sizes = build_membership_from_bool(drug, drug_sets)

# (B) disease targets long
dis = pd.read_csv(f"{ROOT}/data/source/Diease-gene.csv")
# relabel
relabel = {"Disgenet": "DisGeNET", "GEO": "GEO DEG"}
dis["DBlabel"] = dis["Database"].map(lambda d: relabel.get(d, d))
db_order = ["GEO DEG", "GeneCards", "DisGeNET", "DrugBank", "TTD", "OMIM"]
dis_mem, dis_sizes = build_membership_from_long(dis, "DBlabel", "Gene", db_order)

# (C) PPI
ppi_nodes = pd.read_csv(f"{ROOT}/tables/targets/PPI_nodes.csv")
ppi_edges = pd.read_csv(f"{ROOT}/tables/targets/PPI_edges.csv")
HUBS = ["AKT1", "IL6", "TNF", "IL1B", "ESR1", "MYC", "EGFR", "JUN", "HIF1A",
        "MMP9", "CCL2", "CXCL8", "BDNF", "IL10", "IL2", "IL4", "BCL2L1", "NFKBIA"]

# (D) herb-compound-target
hct = pd.read_csv(f"{ROOT}/tables/targets/herb_compound_target_edges.csv")
HERB_SET = {"baishao", "chuanxiong", "gegen", "qianghuo", "sanqi",
            "weilingxian", "yanhusuo"}
HERB_DISPLAY = {
    "baishao": "Baishao", "chuanxiong": "Chuanxiong", "gegen": "Gegen",
    "qianghuo": "Qianghuo", "sanqi": "Sanqi", "weilingxian": "Weilingxian",
    "yanhusuo": "Yanhusuo"}

# ============================================================
# Figure
# ============================================================
fig = plt.figure(figsize=(11, 10))
gs = GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.22,
              left=0.06, right=0.97, top=0.94, bottom=0.05)

# ---------- Panel A ----------
axA = draw_upset(fig, gs[0, 0], drug_mem, drug_sets, drug_sizes,
                 "Drug targets across 4 databases", top_n=10,
                 bar_color=BLUE, set_bar_color=GREY)
axA.text(-0.20, 1.10, "A", transform=axA.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")

# ---------- Panel B ----------
axB = draw_upset(fig, gs[0, 1], dis_mem, db_order, dis_sizes,
                 "Disease targets across databases", top_n=10,
                 bar_color=GREEN, set_bar_color=GREY)
axB.text(-0.20, 1.10, "B", transform=axB.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")

# ---------- Panel C : PPI ----------
axC = fig.add_subplot(gs[1, 0])
G = nx.Graph()
for _, r in ppi_nodes.iterrows():
    G.add_node(r["gene"], degree=int(r["degree"]), is_hub=bool(r["is_hub"]))
for _, r in ppi_edges.iterrows():
    if r["node1"] in G and r["node2"] in G:
        G.add_edge(r["node1"], r["node2"], w=float(r["combined_score"]))
deg = dict(G.degree())
pos = nx.spring_layout(G, seed=42, k=0.55, iterations=200)
hubset = set(HUBS)
node_colors = [VERM if (n in hubset) else "#cfcfcf" for n in G.nodes()]
# node size proportional to degree (use stored degree attr)
node_sizes = [40 + G.nodes[n]["degree"] * 13 for n in G.nodes()]
nx.draw_networkx_edges(G, pos, ax=axC, edge_color=GREY, width=0.4, alpha=0.5)
nx.draw_networkx_nodes(G, pos, ax=axC, node_color=node_colors,
                       node_size=node_sizes, linewidths=0.4,
                       edgecolors="white")
hub_labels = {n: n for n in G.nodes() if n in hubset}
nx.draw_networkx_labels(G, pos, labels=hub_labels, ax=axC, font_size=6.0,
                        font_color=BLACK)
axC.set_title("PPI network (18 hubs)", fontsize=9.5, pad=6)
axC.axis("off")
axC.text(-0.06, 1.06, "C", transform=axC.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")
legC = [Line2D([0], [0], marker="o", color="none", markerfacecolor=VERM,
               markersize=8, label="Hub gene"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#cfcfcf",
               markersize=7, label="Other target")]
axC.legend(handles=legC, loc="lower left", fontsize=6.8,
           bbox_to_anchor=(-0.02, -0.02), handletextpad=0.3)

# ---------- Panel D : herb-compound-target tripartite ----------
axD = fig.add_subplot(gs[1, 1])
T = nx.Graph()
node_class = {}
for _, r in hct.iterrows():
    s, t, et = r["source"], r["target"], r["edge_type"]
    if et in ("H-C", "C-H"):
        herb = s if s in HERB_SET else (t if t in HERB_SET else s)
        comp = t if herb == s else s
        node_class[herb] = "herb"
        node_class.setdefault(comp, "compound")
        T.add_edge(herb, comp)
    elif et == "C-T":
        node_class.setdefault(s, "compound")
        node_class[t] = "target"  # targets are definitive
        T.add_edge(s, t)
# ensure herbs override compound if mislabeled
for n in HERB_SET:
    if n in T:
        node_class[n] = "herb"

herbs = [n for n in T if node_class.get(n) == "herb"]
comps = [n for n in T if node_class.get(n) == "compound"]
targs = [n for n in T if node_class.get(n) == "target"]

# three-column layout: herbs left, compounds middle, targets right
def column_positions(nodes, xval, spread=2.0):
    nodes = sorted(nodes)
    n = len(nodes)
    if n == 1:
        ys = [0.0]
    else:
        ys = np.linspace(spread, -spread, n)
    return {node: (xval, y) for node, y in zip(nodes, ys)}

posT = {}
posT.update(column_positions(herbs, -1.0, spread=1.7))
posT.update(column_positions(comps, 0.0, spread=2.05))
posT.update(column_positions(targs, 1.0, spread=2.05))

class_color = {"herb": ORANGE, "compound": GREEN, "target": SKY}
class_marker = {"herb": "h", "compound": "o", "target": "s"}
class_size = {"herb": 230, "compound": 60, "target": 55}

nx.draw_networkx_edges(T, posT, ax=axD, edge_color=GREY, width=0.3, alpha=0.45)
for cls in ["target", "compound", "herb"]:
    nodes = [n for n in T if node_class.get(n) == cls]
    nx.draw_networkx_nodes(
        T, posT, nodelist=nodes, ax=axD,
        node_color=class_color[cls], node_shape=class_marker[cls],
        node_size=class_size[cls], linewidths=0.4, edgecolors="white")
# label herbs only
herb_labels = {n: HERB_DISPLAY.get(n, n) for n in herbs}
for n, (xx, yy) in posT.items():
    if n in herb_labels:
        axD.text(xx - 0.08, yy, herb_labels[n], ha="right", va="center",
                 fontsize=6.6, fontweight="bold", color=BLACK)
axD.set_title("Herb-compound-target network", fontsize=9.5, pad=6)
axD.axis("off")
axD.margins(0.12)
axD.text(-0.06, 1.06, "D", transform=axD.transAxes, fontsize=15,
         fontweight="bold", ha="left", va="bottom")
legD = [Patch(facecolor=ORANGE, label=f"Herb (n={len(herbs)})"),
        Patch(facecolor=GREEN, label=f"Compound (n={len(comps)})"),
        Patch(facecolor=SKY, label=f"Target (n={len(targs)})")]
axD.legend(handles=legD, loc="lower center", ncol=3, fontsize=6.8,
           bbox_to_anchor=(0.5, -0.04), handletextpad=0.3, columnspacing=1.0)

# ============================================================
# Save
# ============================================================
outdir = f"{ROOT}/figures/composite"
os.makedirs(outdir, exist_ok=True)
png_path = f"{outdir}/Fig2_network.png"
for ext in ["png", "pdf", "svg"]:
    fig.savefig(f"{outdir}/Fig2_network.{ext}", dpi=300,
                bbox_inches="tight", facecolor="white")
plt.close(fig)

# report pixel dims
from PIL import Image
im = Image.open(png_path)
print("SAVED", png_path)
print("PIXELS", im.size)
print("counts A combos:", len(drug_mem), "B combos:", len(dis_mem))
print("PPI nodes:", G.number_of_nodes(), "edges:", G.number_of_edges())
print("herbs/comps/targs:", len(herbs), len(comps), len(targs))
