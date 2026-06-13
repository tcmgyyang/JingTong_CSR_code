# -*- coding: utf-8 -*-
"""
Shared publication figure style for the Jingtong-Granules / cervical-spondylosis manuscript
(target: OMICS / SCI). Single source of truth for font, palette, sizing, panel labels, export.
Import from any script in a subfolder of the project root:
    import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import figstyle as F
    F.apply()                       # sets rcParams
    ... ax uses F.PAL['CON'] ...
    F.panel(ax,'A');  F.save(fig, os.path.join(figdir,'Fig_x'))
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# ---- Okabe-Ito colourblind-safe palette (single source of truth) ----
PAL   = {'CON': '#0072B2', 'CSR': '#D55E00', 'DCM': '#009E73'}          # group triad
PAL2  = {'Normal': '#0072B2', 'Degenerate': '#D55E00'}                   # two-group
NEUTRAL = '#BBBBBB'                                                      # use when colour is redundant with x labels
OKABE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
         '#56B4E9', '#F0E442', '#999999', '#000000']
DIVERGING = 'RdBu_r'      # correlation heatmaps (always center=0, symmetric)
SEQUENTIAL = 'Blues'      # sequential magnitude heatmaps (harmonises with Okabe-Ito blue)

# print sizes (inches): single column ~85 mm, double column ~178 mm
SINGLE = 3.35
DOUBLE = 7.0

def apply():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
        'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.8, 'figure.dpi': 120, 'savefig.dpi': 300,
        'legend.frameon': False,
    })

def panel(ax, letter, x=-0.16, y=1.06):
    """Bold uppercase panel label in axes coords."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='bottom', ha='right')

def star(p):
    """Significance: * <0.05, ** <0.01, *** <0.001; '' (nothing shown) if not significant."""
    try:
        p = float(p)
    except (TypeError, ValueError):
        return ''
    return '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))

def stroke(txt_artist, lw=1.4, fg='white'):
    """White stroke so asterisks/labels read on any cell colour."""
    txt_artist.set_path_effects([pe.withStroke(linewidth=lw, foreground=fg)])

def boxstrip(ax, data, x, y, order, palette, width=0.62, pt=2.9):
    """Polished journal-style boxplot: group-coloured boxes (semi-transparent) + clean
    jittered black points. Replaces dull uniform-grey boxes."""
    import seaborn as sns
    sns.boxplot(data=data, x=x, y=y, order=order, hue=x, hue_order=order, palette=palette,
                legend=False, ax=ax, fliersize=0, width=width, linewidth=0.8,
                medianprops={'color': 'black', 'linewidth': 1.3},
                whiskerprops={'color': 'black', 'linewidth': 0.8},
                capprops={'color': 'black', 'linewidth': 0.8})
    for patch in ax.patches:
        patch.set_edgecolor('black'); patch.set_alpha(0.80)
    sns.stripplot(data=data, x=x, y=y, order=order, ax=ax, color='black',
                  size=pt, alpha=0.45, jitter=0.18, linewidth=0)
    ax.set_xlabel('')

def save(fig, basepath_noext, tiff=False):
    """Export PNG (300 dpi) + PDF (vector, fonts embedded) [+ optional 600-dpi TIFF]."""
    fig.savefig(basepath_noext + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(basepath_noext + '.pdf', bbox_inches='tight')
    if tiff:
        fig.savefig(basepath_noext + '.tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
