# -*- coding: utf-8 -*-
"""
ceRNA step 2: build a lncRNA-miRNA-mRNA network anchored on the disease hub genes,
using experimentally-supported interactions from ENCORI/starBase. Cross-references the
network lncRNAs/circRNAs with the GSE153761 DE ncRNA lists (step 1).
Outputs: ceRNA_edges.csv, ceRNA_nodes.csv (Cytoscape-ready) + Fig_ceRNA_network.png
"""
import os, time, urllib.request, io, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import networkx as nx
OUT=os.path.dirname(os.path.abspath(__file__))

# anchors = (18 multi-source hub genes) INTERSECT (DE in GSE153761 endplate tissue)
ANCHORS=['IL1B','IL6','JUN','BCL2L1','CCL2','CXCL8','TNF']
UA={'User-Agent':'Mozilla/5.0 jingtong-netpharm/1.0'}
def encori(geneType, miRNA='all', target='all', prog=2):
    url=(f'https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType={geneType}'
         f'&miRNA={miRNA}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum={prog}'
         f'&program=None&target={target}&cellType=all')
    for k in range(4):
        try:
            req=urllib.request.Request(url,headers=UA)
            r=urllib.request.urlopen(req,timeout=50).read().decode('utf-8',errors='replace')
            lines=[l for l in r.splitlines() if l and not l.startswith('#')]
            if len(lines)<2: return pd.DataFrame()
            hdr=lines[0].split('\t')
            return pd.DataFrame([x.split('\t') for x in lines[1:]],columns=hdr)
        except Exception:
            time.sleep(2*(k+1))
    return pd.DataFrame()

# 1. miRNAs targeting each anchor mRNA (high confidence: >=3 programs)
edges=[]; mirnas=set()
for g in ANCHORS:
    df=encori('mRNA', target=g, prog=3)
    if df.empty: print('no miRNA for',g); continue
    df['clipExpNum']=pd.to_numeric(df['clipExpNum'],errors='coerce').fillna(0)
    top=df.sort_values('clipExpNum',ascending=False)['miRNAname'].drop_duplicates().head(8).tolist()
    for m in top: edges.append(('miRNA',m,'mRNA',g)); mirnas.add(m)
    print(f'{g}: {len(top)} miRNAs')
    time.sleep(0.4)

# 2. lncRNAs sponging those miRNAs
for m in sorted(mirnas):
    df=encori('lncRNA', miRNA=m, target='all', prog=2)
    if df.empty: continue
    df['clipExpNum']=pd.to_numeric(df['clipExpNum'],errors='coerce').fillna(0)
    top=df.sort_values('clipExpNum',ascending=False)['geneName'].drop_duplicates().head(4).tolist()
    for ln in top: edges.append(('lncRNA',ln,'miRNA',m))
    time.sleep(0.4)

E=pd.DataFrame(edges,columns=['srcType','source','tgtType','target']).drop_duplicates()
E.to_csv(os.path.join(OUT,'ceRNA_edges.csv'),index=False,encoding='utf-8-sig')

# 3. cross-reference network lncRNAs with array DE lncRNAs (symbol-level, best effort)
try:
    de_lnc=pd.read_csv(os.path.join(OUT,'ceRNA_DE_lncRNA.csv'))
    de_names=set(de_lnc['ACC'].astype(str).str.upper())|set(de_lnc.iloc[:,0].astype(str).str.upper())
except Exception: de_names=set()
net_lnc=set(E[E.srcType=='lncRNA'].source)
matched=sorted({l for l in net_lnc if l.upper() in de_names or any(l.upper() in d for d in de_names)})
print(f'\nNetwork: {len(net_lnc)} lncRNAs, {len(mirnas)} miRNAs, {len(ANCHORS)} mRNA hubs')
print('lncRNAs also DE on the GSE153761 array (best-effort match):', matched if matched else 'none cleanly matched (ID-namespace differs; see note)')

# 4. node table + figure
nodes={}
for _,r in E.iterrows():
    nodes[r['source']]=r['srcType']; nodes[r['target']]=r['tgtType']
pd.DataFrame([(n,t) for n,t in nodes.items()],columns=['node','type']).to_csv(os.path.join(OUT,'ceRNA_nodes.csv'),index=False,encoding='utf-8-sig')

G=nx.Graph()
for _,r in E.iterrows(): G.add_edge(r['source'],r['target'])
col={'mRNA':'#C44E52','miRNA':'#55A868','lncRNA':'#4C72B0'}
ncol=[col[nodes[n]] for n in G.nodes()]
nsz=[900 if nodes[n]=='mRNA' else (350 if nodes[n]=='miRNA' else 180) for n in G.nodes()]
plt.figure(figsize=(15,11))
pos=nx.spring_layout(G,k=0.55,seed=42,iterations=80)
nx.draw_networkx_edges(G,pos,alpha=0.25,width=0.6)
nx.draw_networkx_nodes(G,pos,node_color=ncol,node_size=nsz,linewidths=0.4,edgecolors='k')
lab={n:n for n in G.nodes() if nodes[n] in ('mRNA','miRNA')}
nx.draw_networkx_labels(G,pos,labels=lab,font_size=7)
import matplotlib.patches as mp
plt.legend(handles=[mp.Patch(color=c,label=t) for t,c in col.items()],loc='upper right',frameon=False)
plt.title('ceRNA network of Jingtong-Granules hub genes (lncRNA-miRNA-mRNA, ENCORI; GSE153761 DE-anchored)')
plt.axis('off'); plt.tight_layout()
for ext in ('png','pdf'): plt.savefig(os.path.join(OUT,'figures','Fig_ceRNA_network.'+ext),dpi=300,bbox_inches='tight')
print('\nSaved ceRNA_edges.csv, ceRNA_nodes.csv (Cytoscape-ready), Fig_ceRNA_network.png')
print('network edges:',len(E))
