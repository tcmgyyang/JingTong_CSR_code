# Figure 2 — native (no stitching): A drug UpSet, B disease UpSet, C PPI network, D herb-compound network.
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["OMP_NUM_THREADS"] = "1"
import csv
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
import networkx as nx

ROOT = r"H:\毕业设计\网药部分\JingTong_CSR_paper"
T = os.path.join(ROOT, "tables", "targets")
OUT = os.path.join(ROOT, "figures", "composite", "Fig2_network.png")
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8, "svg.fonttype": "none", "pdf.fonttype": 42})
BLUE, VERM, GREEN, SKY, ORANGE, GREY = "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#E69F00", "#bdbdbd"

def letter(ax, L, x=-0.08, y=1.04):
    ax.text(x, y, L, transform=ax.transAxes, fontsize=16, fontweight="bold", ha="left", va="bottom")

# ---------- compact UpSet panel ----------
def upset(fig, cell, sources, labels, gene_src, title, L, topk):
    setsize = {s: sum(1 for v in gene_src.values() if s in v) for s in sources}
    combo = Counter(frozenset(v & set(sources)) for v in gene_src.values() if v & set(sources))
    items = [kv for kv in sorted(combo.items(), key=lambda kv: -kv[1]) if kv[0]][:topk]
    K = len(items); nS = len(sources)
    gs = GridSpecFromSubplotSpec(2, 1, cell, height_ratios=[2.0, 1.7], hspace=0.04)
    axb = fig.add_subplot(gs[0]); axd = fig.add_subplot(gs[1], sharex=axb)
    sizes = [c for _, c in items]
    axb.bar(range(K), sizes, color=BLUE, width=0.62)
    for i, v in enumerate(sizes):
        axb.text(i, v + max(sizes) * 0.02, str(v), ha="center", va="bottom", fontsize=6)
    axb.set_ylim(0, max(sizes) * 1.18); axb.set_ylabel("Intersection", fontsize=7.5)
    axb.tick_params(labelbottom=False, labelsize=7); axb.spines[["top", "right"]].set_visible(False)
    axb.set_title(title, fontsize=9, pad=4)
    for c, (members, _) in enumerate(items):
        filled = [r for r, s in enumerate(sources) if s in members]
        for r in range(nS):
            axd.scatter(c, r, s=34, color=("#3b3b3b" if sources[r] in members else "#dcdcdc"), zorder=3)
        if len(filled) > 1:
            axd.plot([c, c], [min(filled), max(filled)], color="#3b3b3b", lw=1.4, zorder=2)
    for r in range(nS):
        if r % 2 == 0:
            axd.axhspan(r - 0.5, r + 0.5, color="#f5f5f5", zorder=0)
    axd.set_xlim(-0.6, K - 0.4); axd.set_ylim(-0.6, nS - 0.4); axd.invert_yaxis()
    axd.set_yticks(range(nS))
    axd.set_yticklabels([f"{labels[s]} ({setsize[s]})" for s in sources], fontsize=7.5)
    axd.tick_params(axis="x", bottom=False, labelbottom=False); axd.tick_params(axis="y", length=0)
    for sp in axd.spines.values():
        sp.set_visible(False)
    letter(axb, L, x=-0.13)

# ---------- concentric ring layout helper ----------
def ring(nodes, r, start=90.0):
    p = {}; n = max(len(nodes), 1)
    for i, nd in enumerate(nodes):
        a = np.deg2rad(start - 360.0 * i / n)
        p[nd] = (r * np.cos(a), r * np.sin(a))
    return p

def rlabel(ax, name, x, y, off, fs, color):
    ang = np.degrees(np.arctan2(y, x)); r = np.hypot(x, y) + off
    rot = ang if -90 <= ang <= 90 else ang + 180
    ha = "left" if -90 <= ang <= 90 else "right"
    ax.text(r * np.cos(np.radians(ang)), r * np.sin(np.radians(ang)), name, fontsize=fs,
            rotation=rot, rotation_mode="anchor", ha=ha, va="center", color=color)

def cat_grad(nodes, deg, cmap_name, lo=0.3):
    """colour each node within a category by its degree (within-class gradient)."""
    cmap = matplotlib.colormaps[cmap_name]
    vals = [deg[n] for n in nodes]
    lo_v, hi_v = (min(vals), max(vals)) if vals else (0, 1)
    span = (hi_v - lo_v) or 1
    return [cmap(lo + (1 - lo) * ((deg[n] - lo_v) / span)) for n in nodes]

# ---------- C: PPI — hubs on OUTER ring coloured by degree gradient, others grey inside ----------
def ppi_net(ax):
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize
    G = nx.Graph(); deg, hub = {}, {}
    for r in csv.DictReader(open(os.path.join(T, "PPI_nodes.csv"), encoding="utf-8-sig")):
        deg[r["gene"]] = int(r["degree"]); hub[r["gene"]] = (r["is_hub"].strip().lower() == "true")
        G.add_node(r["gene"])
    for r in csv.DictReader(open(os.path.join(T, "PPI_edges.csv"), encoding="utf-8-sig")):
        G.add_edge(r["node1"], r["node2"])
    hubs = sorted([n for n in G if hub.get(n)], key=lambda n: -deg[n])   # high degree first
    oth = sorted([n for n in G if not hub.get(n)], key=lambda n: -deg[n])
    pos = {}; pos.update(ring(hubs, 1.0)); pos.update(ring(oth, 0.46))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#e6e6e6", width=0.45)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=oth, node_color=cat_grad(oth, deg, "Greys", 0.30),
                           node_size=[deg[n] * 5 + 12 for n in oth], edgecolors="white", linewidths=0.3)
    dv = [deg[n] for n in hubs]; norm = Normalize(min(dv), max(dv)); cmap = matplotlib.colormaps["YlOrRd"]
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=hubs, node_color=[cmap(norm(deg[n])) for n in hubs],
                           node_size=[deg[n] * 17 for n in hubs], edgecolors="#333333", linewidths=0.5)
    import matplotlib.patheffects as pe
    for n in hubs:
        x, y = pos[n]
        ax.text(x, y, n, fontsize=5.5, ha="center", va="center", fontweight="bold",
                color="white", path_effects=[pe.withStroke(linewidth=1.5, foreground="black")])
    sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = ax.figure.colorbar(sm, ax=ax, fraction=0.042, pad=0.04); cb.set_label("Hub degree", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_xlim(-1.18, 1.18); ax.set_ylim(-1.15, 1.15); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("PPI network (hubs on outer ring, by degree)", fontsize=9, pad=2)
    letter(ax, "C", x=-0.02)

# ---------- D: herb (inner) – compound (middle) – target (outer), concentric ----------
def herb_net(ax):
    herbs, targets, comps = set(), set(), set(); edges = []
    for r in csv.DictReader(open(os.path.join(T, "herb_compound_target_edges.csv"), encoding="utf-8-sig")):
        s, t, et = r["source"], r["target"], r["edge_type"].upper(); edges.append((s, t))
        if et == "H-C":
            herbs.add(s); comps.add(t)
        elif et == "C-T":
            comps.add(s); targets.add(t)
        elif et == "C-H":                 # artifact: compound -> herb (target IS a herb)
            comps.add(s); herbs.add(t)
    targets -= herbs                       # keep herbs strictly in the inner ring
    comps -= herbs | targets
    G = nx.Graph(); G.add_edges_from(edges)
    HL = {"sanqi": "Sanqi", "chuanxiong": "Chuanxiong", "yanhusuo": "Yanhusuo", "qianghuo": "Qianghuo",
          "baishao": "Baishao", "weilingxian": "Weilingxian", "gegen": "Gegen", "notoginseng": "Sanqi"}
    pos = {}
    pos.update(ring(sorted(herbs), 0.20)); pos.update(ring(sorted(comps), 0.55)); pos.update(ring(sorted(targets), 1.0))
    deg = dict(G.degree())
    tl, cl, hl = sorted(targets), sorted(comps), sorted(herbs)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#eeeeee", width=0.22)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=tl, node_color=cat_grad(tl, deg, "Blues", 0.35),
                           node_size=[16 + deg[n] * 9 for n in tl], edgecolors="white", linewidths=0.2)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=cl, node_color=cat_grad(cl, deg, "Greens", 0.35),
                           node_size=[20 + deg[n] * 4 for n in cl], edgecolors="white", linewidths=0.2)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=hl, node_color=cat_grad(hl, deg, "Oranges", 0.45),
                           node_size=430, node_shape="h", edgecolors="black", linewidths=0.7)
    def trunc(s, n=14):
        return s if len(s) <= n else s[:n - 1] + "…"
    for t in targets:
        rlabel(ax, t, pos[t][0], pos[t][1], 0.02, 4.2, "#222222")
    for c in comps:
        rlabel(ax, trunc(c), pos[c][0], pos[c][1], 0.015, 3.4, "#0a5c3e")
    for h in herbs:
        ax.text(pos[h][0], pos[h][1], HL.get(h, h), fontsize=6.0, ha="center", va="center", fontweight="bold")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="h", color="w", markerfacecolor=ORANGE, markersize=10, label="Herb (7, inner)"),
                       Line2D([0], [0], marker="o", color="w", markerfacecolor=GREEN, markersize=8, label="Compound (44, middle)"),
                       Line2D([0], [0], marker="o", color="w", markerfacecolor=SKY, markersize=7, label="Target (43, outer)")],
              fontsize=6.8, loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.07))
    ax.set_xlim(-1.26, 1.26); ax.set_ylim(-1.24, 1.24); ax.set_aspect("equal"); ax.axis("off")
    ax.margins(0)
    ax.set_title("Herb–compound–target network", fontsize=9, pad=2)
    letter(ax, "D", x=-0.02)

# ---------- data for UpSets ----------
def drug_gene_src():
    gs = {}
    for r in csv.DictReader(open(os.path.join(T, "drug_targets_4source.csv"), encoding="utf-8-sig")):
        s = {k for k in ["TCMSP", "STITCH", "SwissTarget", "HERB"] if r[k].strip().lower() == "true"}
        if s:
            gs[r["gene"]] = s
    return gs

def disease_gene_src():
    gs = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "data", "source", "Diease-gene.csv"), encoding="utf-8-sig")):
        g, db = (r["Gene"] or "").strip(), (r["Database"] or "").strip()
        if g and db:
            gs.setdefault(g, set()).add(db)
    return gs

fig = plt.figure(figsize=(15, 14.5))
top = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.18], hspace=0.14, wspace=0.10)
upset(fig, top[0, 0], ["TCMSP", "STITCH", "SwissTarget", "HERB"],
      {"TCMSP": "TCMSP", "STITCH": "STITCH", "SwissTarget": "SwissTarget", "HERB": "HERB"},
      drug_gene_src(), "Drug targets across 4 databases", "A", 8)
upset(fig, top[0, 1], ["GEO", "GeneCards", "Disgenet", "DrugBank", "TTD", "OMIM"],
      {"GEO": "GEO DEG", "GeneCards": "GeneCards", "Disgenet": "DisGeNET",
       "DrugBank": "DrugBank", "TTD": "TTD", "OMIM": "OMIM"},
      disease_gene_src(), "Disease targets across 6 databases", "B", 10)
ppi_net(fig.add_subplot(top[1, 0]))
herb_net(fig.add_subplot(top[1, 1]))
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT.replace(".png", f".{ext}"), dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote Fig2_network (native 4-panel)")
