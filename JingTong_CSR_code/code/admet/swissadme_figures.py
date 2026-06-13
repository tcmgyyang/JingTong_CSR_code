# -*- coding: utf-8 -*-
"""Publication ADMET figures from SwissADME output.

Two figures:
  Fig_ADMET_SwissADME : two-panel drug-likeness profile
      (a) CONTINUOUS physchem block  -> each column min-max normalised for COLOUR
          (sequential, colourblind-safe) but annotated with the RAW value that the
          same normalisation maps to (colour and text encode the SAME mapping).
      (b) CATEGORICAL drug-likeness flags -> discrete ListedColormap with Yes/No / counts.
  Fig_BOILEDegg       : BOILED-Egg WLOGP-vs-TPSA scatter, GI by colour + BBB by SHAPE,
                        de-collided labels, egg regions in legend.

Reads ../swissadme.csv  (12 active compounds, Jingtong-Granules).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import figstyle as F
F.apply()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
os.makedirs(FIG, exist_ok=True)

d = pd.read_csv(os.path.join(ROOT, 'swissadme.csv'))
d['Molecule'] = d['Molecule'].astype(str)
# tidy compound labels: split camelCase joined names so ticks read cleanly, no underscores
LABEL = {'ginsenosideRh2': 'Ginsenoside Rh2', 'tetrahydropalmatine': 'Tetrahydropalmatine',
         'tetrahydroberberine': 'Tetrahydroberberine'}
d['Label'] = d['Molecule'].map(lambda s: LABEL.get(s, s.capitalize()))

# =====================================================================================
# Figure 1 : two-panel ADMET / drug-likeness profile (subplot_mosaic)
# =====================================================================================
# ---- (a) CONTINUOUS block --------------------------------------------------------
cont_cols = ['MW', 'Consensus Log P', 'TPSA', '#Rotatable bonds', '#H-bond donors',
             '#H-bond acceptors', 'Synthetic Accessibility']
cont_hdr  = ['MW\n(Da)', 'cLogP', 'TPSA\n(Å$^2$)', 'Rot.\nbonds', 'HBD', 'HBA',
             'Synthetic\naccessibility']
INT_COLS  = {'MW', '#Rotatable bonds', '#H-bond donors', '#H-bond acceptors'}
C = d[cont_cols].copy(); C.index = d['Label']
# per-column min-max -> this single normalisation drives BOTH colour AND the raw annotation,
# so a deep-colour cell is genuinely the column max (no more "deep cell shows 0.0").
Cn = (C - C.min()) / (C.max() - C.min())          # 0..1, used for colour only
Craw = C.copy()                                    # raw, used for text only

# ---- (b) CATEGORICAL block -------------------------------------------------------
cat = pd.DataFrame(index=d['Label'])
cat['GI\nabsorption']              = (d['GI absorption'] == 'High').astype(int).values
cat['BBB\npermeant']               = (d['BBB permeant'] == 'Yes').astype(int).values
cat['Lipinski\nviolations']        = d['Lipinski #violations'].values
cat['Veber\nviolations']           = d['Veber #violations'].values
cat['PAINS\nalerts']               = d['PAINS #alerts'].values
cat['Brenk\nalerts']               = d['Brenk #alerts'].values
cat['Bioavailability\nscore']      = d['Bioavailability Score'].values

# text shown in each categorical cell
GI_TXT  = {1: 'High', 0: 'Low'}
BBB_TXT = {1: 'Yes', 0: 'No'}

fig, axd = plt.subplot_mosaic(
    [['a', 'b']], figsize=(7.5, 4.4),
    gridspec_kw={'width_ratios': [7, 7], 'wspace': 0.30, 'left': 0.15,
                 'right': 0.99, 'top': 0.84, 'bottom': 0.30})

n = len(C)

# ---------- panel a : continuous, sequential colourblind-safe (cividis) ----------
axa = axd['a']
im = axa.imshow(Cn.values, cmap=F.SEQUENTIAL, aspect='auto', vmin=0, vmax=1)
axa.set_xticks(range(len(cont_cols))); axa.set_xticklabels(cont_hdr, fontsize=7)
axa.set_yticks(range(n)); axa.set_yticklabels(C.index, fontsize=7)
axa.set_xticks(np.arange(-.5, len(cont_cols), 1), minor=True)
axa.set_yticks(np.arange(-.5, n, 1), minor=True)
axa.grid(which='minor', color='white', lw=0.8); axa.tick_params(which='minor', length=0)
for sp in axa.spines.values(): sp.set_visible(False)
# annotate RAW value; text colour chosen from the SAME normalised value (contrast on cividis)
for i in range(n):
    for j in range(len(cont_cols)):
        val = Craw.values[i, j]
        txt = f'{val:.0f}' if cont_cols[j] in INT_COLS else f'{val:.1f}'
        tc = 'white' if Cn.values[i, j] > 0.55 else 'black'   # Blues: dark=high
        axa.text(j, i, txt, ha='center', va='center', fontsize=6.3, color=tc)
axa.set_title('Physicochemical properties', fontsize=8, pad=6)
# slim colourbar for the continuous block
cb = fig.colorbar(im, ax=axa, fraction=0.040, pad=0.025)
cb.set_label('Per-column min–max scaled', fontsize=6.5)
cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)
F.panel(axa, 'a', x=-0.30, y=1.10)

# ---------- panel b : categorical / drug-likeness flags --------------------------
axb = axd['b']
cat_cols = list(cat.columns)
# Build a colour matrix manually: GI/BBB are binary (favourable vs not),
# violation/alert counts use a 3-step "0 / 1 / >=2" discrete scale,
# bioavailability score is the SwissADME 0.55 (good) vs 0.17 (lower) flag.
# Discrete colourblind-safe set: favourable = teal-blue, neutral = light, unfavourable = vermillion.
GOOD   = '#0072B2'   # favourable (Okabe-Ito blue)
MILD   = '#E69F00'   # one violation / mild caution (Okabe-Ito amber - harmonises better than yellow)
BAD    = '#D55E00'   # >=2 violations / unfavourable (Okabe-Ito vermillion)
LIGHT  = '#ECECEC'   # absent / zero alerts (neutral light grey)

cmat = np.zeros((n, len(cat_cols), 3))  # RGB per cell

def hex2rgb(h):
    h = h.lstrip('#'); return tuple(int(h[k:k+2], 16) / 255 for k in (0, 2, 4))

for i in range(n):
    for j, c in enumerate(cat_cols):
        v = cat.values[i, j]
        if c.startswith('GI'):
            col = GOOD if v == 1 else BAD            # High GI = good, Low = unfavourable
        elif c.startswith('BBB'):
            col = GOOD if v == 1 else LIGHT          # BBB+ highlighted, BBB- neutral (not "bad")
        elif c.startswith('Bioavailability'):
            col = GOOD if v >= 0.5 else BAD          # 0.55 good vs 0.17 low
        else:  # violation / alert counts
            col = LIGHT if v == 0 else (MILD if v == 1 else BAD)
        cmat[i, j] = hex2rgb(col)

axb.imshow(cmat, aspect='auto')
axb.set_xticks(range(len(cat_cols))); axb.set_xticklabels(cat_cols, fontsize=7)
axb.set_yticks(range(n)); axb.set_yticklabels([])           # share rows with panel a
axb.set_xticks(np.arange(-.5, len(cat_cols), 1), minor=True)
axb.set_yticks(np.arange(-.5, n, 1), minor=True)
axb.grid(which='minor', color='white', lw=0.8); axb.tick_params(which='minor', length=0)
for sp in axb.spines.values(): sp.set_visible(False)
# cell text: Yes/No / High/Low / counts / score
for i in range(n):
    for j, c in enumerate(cat_cols):
        v = cat.values[i, j]
        if c.startswith('GI'):
            txt = GI_TXT[v]
        elif c.startswith('BBB'):
            txt = BBB_TXT[v]
        elif c.startswith('Bioavailability'):
            txt = f'{v:.2f}'
        else:
            txt = f'{int(v)}'
        # dark text on light cells, white on saturated cells
        lum = cmat[i, j] @ np.array([0.299, 0.587, 0.114])
        tc = 'black' if lum > 0.6 else 'white'
        axb.text(j, i, txt, ha='center', va='center', fontsize=6.3, color=tc)
axb.set_title('Drug-likeness & alerts', fontsize=8, pad=6)
F.panel(axb, 'b', x=-0.08, y=1.10)

# discrete legend below the categorical panel
handles = [Patch(facecolor=GOOD,  edgecolor='none', label='Favourable (High GI / BBB+ / good score)'),
           Patch(facecolor=MILD,  edgecolor='none', label='1 violation or alert'),
           Patch(facecolor=BAD,   edgecolor='none', label='Unfavourable (Low GI / ≥2 violations / low score)'),
           Patch(facecolor=LIGHT, edgecolor='grey', lw=0.4, label='Absent / 0 alerts')]
axb.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, -0.30),
           ncol=2, fontsize=6.3, frameon=False, handlelength=1.1,
           columnspacing=1.2, handletextpad=0.5)

fig.suptitle('ADMET / drug-likeness profile (SwissADME)', fontsize=9.5, y=0.97)
F.save(fig, os.path.join(FIG, 'Fig_ADMET_SwissADME'))

# =====================================================================================
# Figure 2 : BOILED-Egg  (WLOGP vs TPSA), GI = colour, BBB = SHAPE, de-collided labels
# =====================================================================================
from matplotlib.patches import Polygon
import egg_coords as EGG

fig, ax = plt.subplots(figsize=(6.8, 5.0))

# Egg regions drawn from the AUTHENTIC Daina & Zoete (2016) boundary coordinates
# (digitised via pyBOILEDegg) -> smooth, correct egg, not a hand-estimated ellipse.
# Position alone encodes the prediction: inside white = high GI absorption (HIA),
# inside yolk = BBB permeant. So points need NO redundant colour/shape coding.
ax.add_patch(Polygon(EGG.GIA_COORDS, closed=True, facecolor='#F5F3EC',
                     edgecolor='#C7C2B6', lw=1.0, zorder=0))                 # white (HIA)
ax.add_patch(Polygon(EGG.BBB_COORDS, closed=True, facecolor='#FBE08A',
                     edgecolor='#E0B84A', lw=1.0, zorder=1, alpha=0.92))     # yolk (BBB)

# colour points by chemical class (richer + adds a structural dimension); position still
# carries GI/BBB, so this is informative, not redundant.
CLASS = {'quercetin': 'Flavonoid', 'kaempferol': 'Flavonoid', 'catechin': 'Flavonoid',
         'tetrahydropalmatine': 'Alkaloid', 'tetrahydroberberine': 'Alkaloid',
         'cavidine': 'Alkaloid', 'cryptopine': 'Alkaloid',
         'ginsenosideRh2': 'Triterpenoid saponin', 'stigmasterol': 'Phytosterol',
         'ligustilide': 'Phthalide', 'puerarin': 'Glycoside', 'paeoniflorin': 'Glycoside'}
CLASS_PAL = {'Flavonoid': '#D55E00', 'Alkaloid': '#0072B2', 'Glycoside': '#009E73',
             'Triterpenoid saponin': '#CC79A7', 'Phytosterol': '#E69F00', 'Phthalide': '#56B4E9'}
xs = d['TPSA'].values.astype(float); ys = d['WLOGP'].values.astype(float)
pcols = [CLASS_PAL[CLASS.get(m, 'Flavonoid')] for m in d['Molecule']]
ax.scatter(xs, ys, s=62, c=pcols, edgecolors='black', linewidths=0.6, zorder=3)

ax.set_xlabel('TPSA (Å$^2$)'); ax.set_ylabel('WLOGP')
ax.set_title('BOILED-Egg: gastrointestinal absorption & BBB permeation', fontsize=9)
ax.set_xlim(-15, 178); ax.set_ylim(-3, 8)

# clean automatic label placement with thin leader lines
texts = [ax.text(x, y, lab, fontsize=6.6, color='#222222', zorder=4)
         for x, y, lab in zip(xs, ys, d['Label'])]
try:
    from adjustText import adjust_text
    adjust_text(texts, x=list(xs), y=list(ys), ax=ax,
                arrowprops=dict(arrowstyle='-', color='#9AA0A6', lw=0.5),
                expand=(1.3, 1.6))
except Exception as e:
    print('adjustText fallback:', e)
    for t, x, y in zip(texts, xs, ys):
        t.set_position((x + 3, y + 0.18))

# two-block legend OUTSIDE the plot (keeps the egg uncluttered):
# (i) egg regions, (ii) chemical class colours
reg = [Patch(facecolor='#F5F3EC', edgecolor='#C7C2B6', label='HIA: high GI absorption'),
       Patch(facecolor='#FBE08A', edgecolor='#E0B84A', label='BBB permeant')]
cls = [Line2D([0], [0], marker='o', color='none', markerfacecolor=CLASS_PAL[k],
              markeredgecolor='black', markeredgewidth=0.5, markersize=7, label=k)
       for k in CLASS_PAL]
leg1 = ax.legend(handles=reg, fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                 frameon=False, title='Egg region', title_fontsize=7.5, borderpad=0.4)
ax.add_artist(leg1)
ax.legend(handles=cls, fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 0.74),
          frameon=False, title='Chemical class', title_fontsize=7.5,
          borderpad=0.4, labelspacing=0.5)

F.save(fig, os.path.join(FIG, 'Fig_BOILEDegg'))

print('Saved Fig_ADMET_SwissADME (2-panel) + Fig_BOILEDegg')
print('NOTE: CYP inhibitor columns were "n/d" in this CSV export (read from web table if needed).')
