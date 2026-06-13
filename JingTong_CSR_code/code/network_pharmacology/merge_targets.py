# -*- coding: utf-8 -*-
"""
MERGE the 4 drug-target sources for Jingtong Granules (7 herbs) -> drug-disease
intersection -> STRING PPI -> hub genes. Runs with whatever sources are present:
  1. TCMSP            : ../allTargets.symbol.xls   (Drug,MolId,MolName,Symbol)   [always]
  2. STITCH           : STITCH_targets.csv         (herb,gene,MolName)           [auto, done]
  3. SwissTargetPred. : swisstarget_results/*.csv  (per-compound SwissADME export)[user, tomorrow]
  4. HERB             : HERB_targets.csv           (herb,gene)                   [user, tomorrow]
Outputs source-annotated targets, union & >=2-source sets, disease intersections,
STRING-derived hub genes. Re-run as sources arrive.
"""
import os, glob, urllib.request, urllib.parse
import pandas as pd, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.dirname(os.path.abspath(__file__))
SEVEN = ['baishao', 'chuanxiong', 'gegen', 'qianghuo', 'sanqi', 'weilingxian', 'yanhusuo']
U = lambda s: str(s).strip().upper()

src = {}    # source -> set(gene)

# 1. TCMSP -------------------------------------------------------------------
t = pd.read_csv(os.path.join(ROOT, 'allTargets.symbol.xls'), sep='\t')
src['TCMSP'] = set(t[t.Drug.isin(SEVEN)]['Symbol'].dropna().map(U))

# 2. STITCH ------------------------------------------------------------------
p = os.path.join(OUT, 'STITCH_targets.csv')
if os.path.exists(p):
    s = pd.read_csv(p)
    if len(s): src['STITCH'] = set(s['gene'].dropna().map(U))

# 3. SwissTargetPrediction (pre-built swiss_targets.csv via build_herb_swiss.py) ---
p = os.path.join(OUT, 'swiss_targets.csv')
if os.path.exists(p):
    s = pd.read_csv(p)
    if len(s): src['SwissTarget'] = set(s['gene'].dropna().map(U))

# 4. HERB --------------------------------------------------------------------
for name in ('HERB_targets.csv', 'HERB_targets.txt'):
    p = os.path.join(OUT, name)
    if os.path.exists(p):
        h = pd.read_csv(p)
        gcol = next((c for c in h.columns if c.strip().lower() in ('gene', 'target', 'gene symbol', 'symbol')), h.columns[-1])
        src['HERB'] = set(h[gcol].dropna().map(U)); break

print('sources present:', {k: len(v) for k, v in src.items()})

# ---- source-annotated union ----
all_genes = set().union(*src.values())
rows = [{'gene': g, **{k: (g in v) for k, v in src.items()},
         'n_sources': sum(g in v for v in src.values())} for g in sorted(all_genes)]
ann = pd.DataFrame(rows).sort_values(['n_sources', 'gene'], ascending=[False, True])
ann.to_csv(os.path.join(OUT, 'drug_targets_4source.csv'), index=False, encoding='utf-8-sig')
union = set(ann.gene); ge2 = set(ann[ann.n_sources >= 2].gene)
print(f'drug targets: union={len(union)}  >=2 sources={len(ge2)}')

# ---- disease targets (reconstructed: DB + GEO DEGs) ----
dis = set()
try: dis |= set(pd.read_csv(os.path.join(ROOT, 'Diease-gene.csv'), encoding='utf-8-sig')['Gene'].dropna().map(U))
except Exception: pass
try: dis |= set(pd.read_csv(os.path.join(ROOT, 'GeneCards_1.csv'))['Gene Symbol'].dropna().map(U))
except Exception: pass
try: dis |= set(pd.read_csv(os.path.join(ROOT, '87geoPharm(1)', '07.diffPvalue', 'diff.xls'), sep='\t')['id'].dropna().map(U))
except Exception: pass
print('disease targets:', len(dis))

inter_union = sorted(union & dis); inter_ge2 = sorted(ge2 & dis)
pd.Series(inter_union).to_csv(os.path.join(OUT, 'intersection_union.csv'), index=False, header=['gene'])
pd.Series(inter_ge2).to_csv(os.path.join(OUT, 'intersection_2source.csv'), index=False, header=['gene'])
print(f'drug-disease intersection: union={len(inter_union)}  >=2 sources={len(inter_ge2)}')

# ---- STRING PPI on the (union) intersection -> hub genes ----
def string_ppi(genes):
    genes = [g for g in genes if g]
    url = 'https://string-db.org/api/tsv/network'
    data = urllib.parse.urlencode({'identifiers': '%0d'.join(genes), 'species': '9606',
                                   'caller_identity': 'jingtong_netpharm'}).encode()
    req = urllib.request.Request(url, data=data, headers={'User-Agent': 'jingtong-netpharm/1.0'})
    txt = urllib.request.urlopen(req, timeout=90).read().decode('utf-8', 'replace')
    edges = []
    for ln in txt.splitlines()[1:]:
        f = ln.split('\t')
        if len(f) >= 6: edges.append((f[2], f[3], float(f[5])))   # preferredName_A,B, combined_score
    return edges

def string_hub(genes, tag):
    if len(genes) < 3: return
    try:
        import networkx as nx
        G = nx.Graph()
        for a, b, sc in string_ppi(genes):
            if sc >= 0.4: G.add_edge(a, b, weight=sc)
        if not G.number_of_nodes(): return
        deg = dict(G.degree()); bet = nx.betweenness_centrality(G)
        md, mb = np.median(list(deg.values())), np.median(list(bet.values()))
        hub = sorted([n for n in G if deg[n] >= md and bet[n] >= mb], key=lambda n: -deg[n])
        pd.DataFrame([(n, deg[n], round(bet[n], 4)) for n in hub],
                     columns=['gene', 'degree', 'betweenness']).to_csv(
            os.path.join(OUT, f'hub_genes_{tag}.csv'), index=False, encoding='utf-8-sig')
        print(f'[{tag}] STRING PPI: {G.number_of_nodes()} nodes / {G.number_of_edges()} edges; '
              f'HUB (deg & betw >= median): {len(hub)}')
        print('   ', ', '.join(hub[:20]))
    except Exception as e:
        print(f'[{tag}] STRING PPI failed:', e)

print('\n== PRIMARY: high-confidence >=2-source intersection ==')
string_hub(inter_ge2, '2source')
print('== SENSITIVITY: full union intersection ==')
string_hub(inter_union, 'union')
print('\nWrote drug_targets_4source.csv, intersection_union.csv, intersection_2source.csv, hub_genes_4source.csv')
