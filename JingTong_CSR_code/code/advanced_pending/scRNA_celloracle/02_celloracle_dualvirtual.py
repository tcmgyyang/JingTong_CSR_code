# -*- coding: utf-8 -*-
"""
STEP 2 (server) -- CellOracle GRN + DUAL VIRTUAL perturbation (virtual knockout +
virtual over-expression) of the perturbable hub genes, on the annotated scRNA.

NOTE on biology: CellOracle simulates perturbation of *regulatory* genes (TFs in the base
GRN). Among the 18 multi-source hubs the TF/regulator ones (e.g. JUN, MYC, ESR1, HIF1A) can
be virtually KO'd / over-expressed; effector hubs (MMP9, IL1B, IL6, TNF, CCL2, CXCL8 ...)
are not GRN regulators -> they are reported by *cell-type localisation* (step 1), not
virtually perturbed. The script auto-selects 18hub INTERSECT oracle-regulators.

For each perturbable hub it runs BOTH:
  - virtual knockout  : expression -> 0
  - virtual over-expr : expression -> ~2x its observed max
and saves the perturbation-vector (quiver) on the UMAP + a perturbation-score summary.

Usage:  py 02_celloracle_dualvirtual.py adata_annotated.h5ad [out_dir]
Env:    pip install celloracle   (also needs scanpy; see README for the conda recipe)
"""
import os, sys
import numpy as np, pandas as pd
import scanpy as sc, celloracle as co
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

ADATA = sys.argv[1] if len(sys.argv) > 1 else 'adata_annotated.h5ad'
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(OUT, 'perturb_figs'); os.makedirs(FIG, exist_ok=True)
HUBS = ['AKT1','IL6','TNF','IL1B','ESR1','MYC','EGFR','JUN','HIF1A','MMP9','CCL2',
        'CXCL8','BDNF','IL10','IL2','IL4','BCL2L1','NFKBIA']

adata = sc.read_h5ad(ADATA)
if 'raw_count' in adata.layers:                      # CellOracle wants raw counts in .X
    adata.X = adata.layers['raw_count'].copy()
adata = adata.raw.to_adata() if adata.raw is not None and adata.n_vars < 4000 else adata

# ---- build Oracle + GRN ----
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name='cell_type', embedding_name='X_umap')
base_GRN = co.data.load_human_promoter_base_GRN()     # built-in human base GRN (promoter/TSS)
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()
n_comps = min(50, np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])) > 0.002)[0][0] if False else 50)
oracle.knn_imputation(n_pca_dims=50, k=max(25, int(0.025 * oracle.adata.shape[0])),
                      balanced=True, b_sight=int(0.1 * oracle.adata.shape[0]),
                      b_maxl=int(0.05 * oracle.adata.shape[0]), n_jobs=4)

links = oracle.get_links(cluster_name_for_GRN_unit='cell_type', alpha=10, n_jobs=4, verbose_level=0)
links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

regulators = set(oracle.active_regulatory_genes) if hasattr(oracle, 'active_regulatory_genes') \
             else set(oracle.adata.var_names)
perturbable = [g for g in HUBS if g in oracle.adata.var_names and g in regulators]
print('perturbable hub TFs:', perturbable)
pd.Series(perturbable, name='perturbable_hub').to_csv(os.path.join(OUT, 'perturbable_hubs.csv'), index=False)

def run_perturb(gene, value, tag):
    oracle.simulate_shift(perturb_condition={gene: value}, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    fig, ax = plt.subplots(figsize=(5, 5))
    oracle.plot_simulation_flow_on_grid(scale=25, ax=ax)
    ax.set_title(f'{tag}: {gene}', fontsize=11); ax.axis('off')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f'{gene}_{tag}.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    return float(np.abs(oracle.delta_embedding).sum())

summary = []
for g in perturbable:
    mx = float(np.asarray(oracle.adata[:, g].layers['imputed_count']).max()) if 'imputed_count' in oracle.adata.layers else 1.0
    ko = run_perturb(g, 0.0, 'KO')                      # virtual knockout
    oe = run_perturb(g, mx * 2.0, 'OE')                # virtual over-expression
    summary.append([g, round(ko, 3), round(oe, 3)])
    print(f'{g}: KO shift={ko:.2f}  OE shift={oe:.2f}')

pd.DataFrame(summary, columns=['gene', 'KO_total_shift', 'OE_total_shift']).to_csv(
    os.path.join(OUT, 'dualvirtual_summary.csv'), index=False, encoding='utf-8-sig')
# GRN cluster degree (which hubs are central regulators)
try:
    links.plot_degree_distributions(plot_model=True, save=os.path.join(FIG, 'GRN_degree'))
except Exception:
    pass
print('DONE -> perturb_figs/<gene>_{KO,OE}.png + dualvirtual_summary.csv')
