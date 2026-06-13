# -*- coding: utf-8 -*-
"""
STEP 1 (run on a memory-rich machine/server) -- scRNA QC + clustering + cell-type
annotation for the cervical/disc endplate dataset (GSE242040 cartilage endplate, or the
GSE165722 disc atlas). Output: adata_annotated.h5ad (raw counts kept in a layer for
CellOracle) + UMAP / hub-gene figures.

Usage:  py 01_annotate.py <10x_mtx_dir | path.h5ad> [out_dir]
Env:    pip install scanpy leidenalg igraph
"""
import os, sys
import scanpy as sc, numpy as np, pandas as pd
sc.settings.verbosity = 1

INP = sys.argv[1] if len(sys.argv) > 1 else 'DATA'        # 10x dir or .h5ad
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))
sc.settings.figdir = OUT

HUBS = ['AKT1','IL6','TNF','IL1B','ESR1','MYC','EGFR','JUN','HIF1A','MMP9','CCL2',
        'CXCL8','BDNF','IL10','IL2','IL4','BCL2L1','NFKBIA']
# marker panels for spinal cartilage-endplate / disc cell types
MARKERS = {
    'Chondrocyte':        ['COL2A1','ACAN','SOX9','COL9A1','COL11A1'],
    'Hypertrophic/Osteo': ['COL10A1','RUNX2','SPP1','IBSP','ALPL','MMP13'],
    'Fibrochondrocyte':   ['COL1A1','COL1A2','COL3A1','DCN','LUM'],
    'Endothelial':        ['PECAM1','VWF','CLDN5','CDH5'],
    'Immune/Macrophage':  ['PTPRC','CD68','LYZ','CD14','AIF1','CD163'],
    'Proliferating':      ['MKI67','TOP2A','PCNA','CCNB1'],
    'Nucleus pulposus':   ['KRT19','KRT8','CD24','T','CHST3'],
}

import glob, re
print('loading', INP)
def _fix_genes(a):
    for gc in ('gene_names', 'Gene', 'GeneName', 'gene', 'var_names', 'features'):  # GSE160756 uses 'gene_names'
        if gc in a.var.columns:
            a.var_names = a.var[gc].astype(str).values; break
    for cc in ('cell_names', 'CellID', 'obs_names'):                                 # GSE160756 uses 'cell_names'
        if cc in a.obs.columns:
            a.obs_names = a.obs[cc].astype(str).values; break
    a.var_names_make_unique(); a.obs_names_make_unique(); return a
looms = sorted(glob.glob(os.path.join(INP, '*.loom'))) if os.path.isdir(INP) else []
if looms:                                                  # GSE160756 atlas: one .loom per sample
    ads = []
    for f in looms:
        a = _fix_genes(sc.read_loom(f, sparse=True)); base = os.path.basename(f).split('.')[0]
        a.obs['sample'] = base
        a.obs['tissue'] = ('NP' if re.search('NP', base, re.I) else
                           'CEP' if re.search(r'C?EP', base, re.I) else
                           'AF' if re.search('AF', base, re.I) else 'NA')
        a.obs_names = [base + '_' + bc for bc in a.obs_names]
        ads.append(a); print('  ', base, a.shape, a.obs['tissue'][0])
    adata = ads[0].concatenate(*ads[1:], join='outer', index_unique=None) if len(ads) > 1 else ads[0]
elif INP.endswith('.h5ad'):
    adata = sc.read_h5ad(INP)
else:
    adata = sc.read_10x_mtx(INP, var_names='gene_symbols', cache=True)
adata.var_names_make_unique()
print('combined:', adata.shape, '| tissues:', adata.obs['tissue'].value_counts().to_dict() if 'tissue' in adata.obs else 'n/a')

# ---- QC ----
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)
adata = adata[(adata.obs.n_genes_by_counts > 200) & (adata.obs.n_genes_by_counts < 7000) &
              (adata.obs.pct_counts_mt < 20)].copy()
sc.pp.filter_genes(adata, min_cells=3)
print('after QC:', adata.shape)

adata.layers['raw_count'] = adata.X.copy()                # keep raw counts for CellOracle
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2500)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata); sc.tl.leiden(adata, resolution=0.6)
print('clusters:', adata.obs.leiden.nunique())

# ---- annotate clusters by marker score ----
adata_full = adata.raw.to_adata()
for ct, gs in MARKERS.items():
    g = [x for x in gs if x in adata_full.var_names]
    if g: sc.tl.score_genes(adata_full, g, score_name='sc_' + ct)
score_cols = ['sc_' + c for c in MARKERS if 'sc_' + c in adata_full.obs]
adata.obs = adata.obs.join(adata_full.obs[score_cols])
mean_by_cl = adata.obs.groupby('leiden')[score_cols].mean()
cl2type = {cl: row.idxmax().replace('sc_', '') for cl, row in mean_by_cl.iterrows()}
adata.obs['cell_type'] = adata.obs.leiden.map(cl2type).astype('category')
print('cluster -> cell type:', cl2type)

# ---- figures ----
_cols = ['leiden', 'cell_type'] + (['tissue'] if 'tissue' in adata.obs else [])
sc.pl.umap(adata, color=_cols, save='_clusters.png', show=False)
hub_present = [h for h in HUBS if h in adata.raw.var_names]
sc.pl.umap(adata, color=hub_present[:12], use_raw=True, save='_hub_expr.png', show=False, ncols=4)
sc.pl.dotplot(adata, hub_present, groupby='cell_type', use_raw=True,
              standard_scale='var', save='_hub_dotplot.png', show=False)

adata.write(os.path.join(OUT, 'adata_annotated.h5ad'))
adata.obs[['leiden', 'cell_type']].to_csv(os.path.join(OUT, 'cell_annotation.csv'))
print('SAVED adata_annotated.h5ad + UMAP/dotplot. Next: py 02_celloracle_dualvirtual.py adata_annotated.h5ad')
