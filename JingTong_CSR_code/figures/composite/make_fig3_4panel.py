# -*- coding: utf-8 -*-
"""Fig 3 (4-panel): a = GO-BP ORA, b = KEGG ORA, c = GSVA, d = single-gene GSEA.
GO/KEGG come from Enrichr ORA on the 43 intersection targets (tables/enrichment/*_ORA.csv);
GSVA/GSEA reuse the existing key-pathway summaries. Reproduces Results 3.3."""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
    "font.size":8,"axes.linewidth":0.8,"axes.spines.top":False,"axes.spines.right":False,
    "legend.frameon":False,"svg.fonttype":"none","pdf.fonttype":42})
BLUE="#0072B2"; ORANGE="#E69F00"; GREEN="#009E73"; VERM="#D55E00"; GREY="#999999"; BLACK="#000000"

BASE = r"O:\JingTong_CSR_paper"
ENR = os.path.join(BASE,"tables","enrichment")
OUT = os.path.join(BASE,"figures","composite","Fig3_pathways_4panel")

keg = pd.read_csv(os.path.join(ENR,"KEGG_ORA.csv"))
gsva= pd.read_csv(os.path.join(ENR,"GSVA_key_pathways.csv"))
gsea= pd.read_csv(os.path.join(ENR,"GSEA_key_pathways_summary.csv"))
from matplotlib.patches import Patch
SKY="#56B4E9"
GO_FILES={"BP":"GO_BP_ORA.csv","CC":"GO_CC_ORA.csv","MF":"GO_MF_ORA.csv"}
GO_COLORS={"BP":GREEN,"CC":SKY,"MF":ORANGE}
GO_NAMES={"BP":"Biological process","CC":"Cellular component","MF":"Molecular function"}

# KEGG: drop promiscuous disease/infection/cancer gene-sets -> keep signaling/mechanism pathways
BAD = ("infection","disease","cancer","sarcoma","leukemia","carcinoma","hepatitis",
       "influenza","tuberculosis","measles","Epstein","virus","Yersinia","Chagas",
       "Salmonella","Shigellosis","Amoebiasis","Toxoplasmosis","Malaria","Pertussis",
       "Legionellosis","Leishmaniasis","platelet")
keg_sig = keg[~keg["term"].str.lower().apply(lambda t: any(b.lower() in t for b in BAD))].copy()

def clean_go(t): return t.rsplit(" (GO:",1)[0]

fig = plt.figure(figsize=(12.6, 10.4))
gs = GridSpec(2, 2, figure=fig, wspace=0.62, hspace=0.42,
              left=0.085, right=0.93, top=0.93, bottom=0.08,
              height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.05])
axA = fig.add_subplot(gs[0,0]); axB = fig.add_subplot(gs[0,1])
axC = fig.add_subplot(gs[1,0]); axD = fig.add_subplot(gs[1,1])

def wrap_term(t, width=34):
    """Soft-wrap a long term name onto at most two lines for tick labels."""
    if len(t) <= width:
        return t
    words = t.split(" ")
    line1, line2, cur = [], [], 0
    for w in words:
        if cur + len(w) + (1 if line1 else 0) <= width and not line2:
            line1.append(w); cur += len(w) + 1
        else:
            line2.append(w)
    out = " ".join(line1)
    if line2:
        out += "\n" + " ".join(line2)
    return out

def size_from_count(n, nmin, nmax, smin=24, smax=190):
    """Map gene count -> marker area (pts^2), linear over the observed range."""
    if nmax == nmin:
        return np.full_like(np.asarray(n, float), (smin + smax) / 2.0)
    return smin + (np.asarray(n, float) - nmin) / (nmax - nmin) * (smax - smin)

# ---------------------------------------------------------------------------
# Panel A: GO enrichment as a faceted bubble (lollipop-dot) plot.
# Three ontology blocks (BP/CC/MF), top 4 each, kept visually separated and
# colour-coded by ontology; x = -log10(adj P); dot SIZE encodes gene count.
# A thin stem links each dot to the axis so the ranking still reads as a bar
# chart would, but the figure now carries an extra (count) dimension.
# ---------------------------------------------------------------------------
rows = []  # (ontology, label, nlp, n_overlap)
for onto in ("BP", "CC", "MF"):
    d = pd.read_csv(os.path.join(ENR, GO_FILES[onto])).sort_values("adj_p").head(4)
    for _, r in d.iterrows():
        rows.append((onto, clean_go(r["term"]),
                     -np.log10(max(r["adj_p"], 1e-300)), int(r["n_overlap"])))
# y positions top->bottom (BP block on top) with a gap between ontology blocks
cur, yv, last = 0.0, [], None
for onto, *_ in rows:
    if last is not None and onto != last:
        cur += 0.9
    yv.append(cur); cur += 1.0; last = onto
ymax = max(yv); yv = [ymax - v for v in yv]
vals = np.array([r[2] for r in rows])
cols = [GO_COLORS[r[0]] for r in rows]
ncs  = np.array([r[3] for r in rows])
labels = [wrap_term(r[1]) for r in rows]
# shared gene-count scale across panels A and B so dot sizes are comparable
allcnt = np.concatenate([ncs, keg_sig.sort_values("adj_p").head(10)["n_overlap"].values])
CMIN, CMAX = int(allcnt.min()), int(allcnt.max())
szA = size_from_count(ncs, CMIN, CMAX)
xmaxA = vals.max() * 1.20
# subtle ontology background bands + stems
for onto in ("BP", "CC", "MF"):
    idx = [i for i, r in enumerate(rows) if r[0] == onto]
    ylo, yhi = min(yv[i] for i in idx) - 0.5, max(yv[i] for i in idx) + 0.5
    axA.axhspan(ylo, yhi, color=GO_COLORS[onto], alpha=0.06, zorder=0)
    # ontology tag on the right margin of each band
    axA.text(xmaxA*0.992, yhi-0.18, GO_NAMES[onto], ha="right", va="top",
             fontsize=6.4, color=GO_COLORS[onto], fontweight="bold", style="italic", zorder=4)
for yi, v, c in zip(yv, vals, cols):
    axA.plot([0, v], [yi, yi], color=c, lw=1.1, alpha=0.55, zorder=2, solid_capstyle="round")
axA.scatter(vals, yv, s=szA, c=cols, edgecolor="white", linewidth=0.8, zorder=3)
# gene-count number printed inside / beside each dot
for yi, v, nc, s in zip(yv, vals, ncs, szA):
    axA.text(v + xmaxA*0.018, yi, str(int(nc)), va="center", ha="left",
             fontsize=6.2, color=GREY, zorder=4)
axA.set_yticks(yv); axA.set_yticklabels(labels, fontsize=6.8)
axA.set_xlim(0, xmaxA); axA.set_ylim(min(yv)-0.7, max(yv)+0.7)
axA.set_xlabel("$-\\log_{10}$ adjusted $P$", fontsize=8)
axA.tick_params(axis="both", length=3, width=0.8)
axA.set_title("GO enrichment", fontsize=10, fontweight="bold", pad=8)
axA.text(-0.02, 1.05, "A", transform=axA.transAxes, fontsize=15, fontweight="bold", ha="right", va="bottom")
# (ontology identity is already carried by the per-band italic tags + colour;
#  a separate colour legend would be redundant, so it is intentionally omitted)

# ---------------------------------------------------------------------------
# Panel B: KEGG signalling pathways as a continuous-gradient bubble plot.
# x = -log10(adj P); dot SIZE = gene count (n_overlap); dot COLOUR = a
# sequential significance gradient (more significant = deeper). Companion
# colourbar + a size legend make this the classic Cell/Nature enrichment dot
# plot, distinct in visual language from the categorical Panel A.
# ---------------------------------------------------------------------------
dB = keg_sig.sort_values("adj_p").head(10).copy()
dB["nlp"] = -np.log10(dB["adj_p"].clip(lower=1e-300))
dB = dB.iloc[::-1].reset_index(drop=True)   # most significant at top after plotting
yB = np.arange(len(dB))
nlpB = dB["nlp"].values
ncB = dB["n_overlap"].values.astype(int)
szB = size_from_count(ncB, CMIN, CMAX)
xmaxB = nlpB.max() * 1.20
normB = Normalize(vmin=nlpB.min(), vmax=nlpB.max())
cmapB = plt.get_cmap("YlOrRd")
# clamp away from the very pale tail so small-significance dots stay visible
colB = cmapB(0.30 + 0.65 * normB(nlpB))
# faint guide stems
for yi, v in zip(yB, nlpB):
    axB.plot([0, v], [yi, yi], color=GREY, lw=0.7, alpha=0.35, zorder=1, solid_capstyle="round")
scB = axB.scatter(nlpB, yB, s=szB, c=colB, edgecolor="white", linewidth=0.8, zorder=3)
for yi, v, nc in zip(yB, nlpB, ncB):
    axB.text(v + xmaxB*0.016, yi, str(int(nc)), va="center", ha="left",
             fontsize=6.2, color=GREY, zorder=4)
axB.set_yticks(yB); axB.set_yticklabels([wrap_term(t) for t in dB["term"]], fontsize=6.8)
axB.set_xlim(0, xmaxB); axB.set_ylim(-0.7, len(dB)-0.3)
axB.set_xlabel("$-\\log_{10}$ adjusted $P$", fontsize=8)
axB.tick_params(axis="both", length=3, width=0.8)
axB.set_title("KEGG pathway", fontsize=10, fontweight="bold", pad=8)
axB.text(-0.02, 1.05, "B", transform=axB.transAxes, fontsize=15, fontweight="bold", ha="right", va="bottom")
# significance colourbar (slim, to the right of panel B)
smB = ScalarMappable(norm=normB, cmap=cmapB)
cbB = fig.colorbar(smB, ax=axB, fraction=0.040, pad=0.14, aspect=22)
cbB.set_label("$-\\log_{10}$ adj $P$", fontsize=6.8)
cbB.ax.tick_params(labelsize=6, length=2); cbB.outline.set_linewidth(0.6)
# gene-count size legend (lower-right, inside plot)
cnt_ticks = [CMIN, int(round((CMIN+CMAX)/2)), CMAX]
size_handles = [Line2D([0], [0], marker="o", linestyle="none",
                       markerfacecolor=GREY, markeredgecolor="white", markeredgewidth=0.6,
                       markersize=np.sqrt(size_from_count(c, CMIN, CMAX)),
                       label=str(c)) for c in cnt_ticks]
leg_size = axB.legend(handles=size_handles, loc="lower right",
                      bbox_to_anchor=(1.0, 0.02), fontsize=6.2,
                      frameon=False, title="Gene count", title_fontsize=6.4,
                      labelspacing=1.0, handletextpad=0.6, borderpad=0.3)
leg_size._legend_box.align = "left"
axB.add_artist(leg_size)

# Panel C: GSVA lollipop
gC = gsva.sort_values("delta_CERVvsCON", ascending=True).reset_index(drop=True)
yc = np.arange(len(gC)); delta = gC["delta_CERVvsCON"].values; pv = gC["p"].values
cols = [VERM if v>=0 else BLUE for v in delta]
for yi,di,ci in zip(yc,delta,cols): axC.plot([0,di],[yi,yi],color=ci,lw=2.2,zorder=2,solid_capstyle="round")
axC.scatter(delta,yc,color=cols,s=52,zorder=3,edgecolor="white",linewidth=0.6)
axC.axvline(0,color=BLACK,lw=0.8)
axC.set_yticks(yc); axC.set_yticklabels(gC["pathway"].values,fontsize=7.3)
axC.set_xlabel("GSVA Δ (degeneration − control)",fontsize=8)
xpad=max((delta.max()-delta.min())*0.04,0.006)
for yi,di,p in zip(yc,delta,pv):
    t=f"p={p:.3g}"
    axC.text(di+(xpad if di>=0 else -xpad),yi,t,va="center",ha=("left" if di>=0 else "right"),fontsize=6.3,color=GREY)
sp=delta.max()-delta.min()
axC.set_xlim(min(delta.min(),0)-sp*0.42,max(delta.max(),0)+sp*0.42); axC.set_ylim(-0.7,len(gC)-0.3)
axC.tick_params(axis="both",length=3,width=0.8)
axC.set_title("GSVA pathway activity",fontsize=10,fontweight="bold",pad=8)
axC.text(-0.02,1.05,"C",transform=axC.transAxes,fontsize=15,fontweight="bold",ha="right",va="bottom")

# Panel D: single-gene GSEA NES heatmap
nes=gsea.pivot(index="gene",columns="Term",values="NES"); fdr=gsea.pivot(index="gene",columns="Term",values="FDR")
go_order=nes.mean(1).sort_values(ascending=False).index; nes=nes.loc[go_order]; fdr=fdr.loc[go_order]
to=nes.mean(0).sort_values(ascending=False).index; nes=nes[to]; fdr=fdr[to]
data=nes.values.astype(float); fdrv=fdr.values.astype(float)
half=max(1.0-np.nanmin(data),np.nanmax(data)-1.0); norm=TwoSlopeNorm(vmin=1.0-half,vcenter=1.0,vmax=1.0+half)
nr,nc=data.shape
im=axD.pcolormesh(np.arange(nc+1),np.arange(nr+1),data,cmap="RdBu_r",norm=norm,edgecolors="white",linewidth=1.0,shading="flat")
axD.set_aspect("auto"); axD.set_xlim(0,nc); axD.set_ylim(0,nr); axD.invert_yaxis()
xc=np.arange(nc)+0.5; yc2=np.arange(nr)+0.5
axD.set_xticks(xc); axD.set_yticks(yc2)
axD.set_xticklabels([t.replace(" signaling pathway","").replace(" in diabetic complications","") for t in nes.columns],fontsize=6.8,rotation=40,ha="right",rotation_mode="anchor")
axD.set_yticklabels(nes.index,fontsize=7.3); axD.tick_params(length=0)
for s in axD.spines.values(): s.set_visible(False)
for i in range(nr):
    for j in range(nc):
        if not np.isnan(fdrv[i,j]) and fdrv[i,j]<0.05:
            axD.text(xc[j],yc2[i],"*",ha="center",va="center",fontsize=11,color=BLACK,fontweight="bold")
axD.set_title("Single-gene GSEA (NES)",fontsize=10,fontweight="bold",pad=8)
axD.text(-0.02,1.05,"D",transform=axD.transAxes,fontsize=15,fontweight="bold",ha="right",va="bottom")
cb=fig.colorbar(im,ax=axD,fraction=0.046,pad=0.02); cb.set_label("NES",fontsize=8)
cb.ax.tick_params(labelsize=7,length=2); cb.outline.set_linewidth(0.6); cb.ax.axhline(norm(1.0),color=BLACK,lw=0.8)

for ext in ("png","pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
from PIL import Image
with Image.open(f"{OUT}.png") as im_: print("PNG_SIZE",im_.size)
print("KEGG panel terms:", keg_sig.sort_values('adj_p').head(10)['term'].tolist())
print("DONE ->", OUT+".png")
