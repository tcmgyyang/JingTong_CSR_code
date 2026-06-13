# -*- coding: utf-8 -*-
"""Build HERB_targets.csv and swiss_targets.csv from the user's downloads, with the
correct file->herb / file->compound mapping, ready for merge_targets.py."""
import os, glob, re
import pandas as pd
OUT = os.path.dirname(os.path.abspath(__file__))
U = lambda s: str(s).strip().upper()

# ---------- HERB: 7 files -> herbs (table order; 三七 = 6_9 file) ----------
HERB_MAP = {
    'herb_target2026_6_9.xlsx':       'sanqi',        # 三七 (06-09)
    'herb_target2026_6_10.xlsx':      'chuanxiong',   # 川芎 (no suffix)
    'herb_target2026_6_10 (1).xlsx':  'yanhusuo',     # 延胡索
    'herb_target2026_6_10 (2).xlsx':  'qianghuo',     # 羌活
    'herb_target2026_6_10 (3).xlsx':  'baishao',      # 白芍
    'herb_target2026_6_10 (4).xlsx':  'weilingxian',  # 威灵仙
    'herb_target2026_6_10 (5).xlsx':  'gegen',        # 葛根
}
herb_rows = []
for fn, herb in HERB_MAP.items():
    p = os.path.join(OUT, 'herb', fn)
    if not os.path.exists(p):
        print('MISSING HERB file:', fn); continue
    df = pd.read_excel(p)
    gcol = 'Target name' if 'Target name' in df.columns else df.columns[1]
    genes = df[gcol].dropna().map(U).unique()
    for g in genes:
        herb_rows.append([herb, g])
    print(f'  HERB {herb:12s} <- {fn:34s}  {len(genes)} targets')
herb_df = pd.DataFrame(herb_rows, columns=['herb', 'gene']).drop_duplicates()
herb_df.to_csv(os.path.join(OUT, 'HERB_targets.csv'), index=False, encoding='utf-8-sig')
print('HERB total unique genes:', herb_df.gene.nunique(), '\n')

# ---------- SwissTarget: 43 files in order -> compound -> herb(s) ----------
inp = pd.read_csv(os.path.join(OUT, 'swisstarget_input.csv'))   # ordered rows = file order
inp = inp.reset_index(drop=True)
def filenum(path):
    m = re.search(r'\((\d+)\)', os.path.basename(path))
    return int(m.group(1)) if m else 0          # no-suffix -> 0 (first)
files = sorted(glob.glob(os.path.join(OUT, 'swisstarget', 'SwissTargetPrediction*.csv')), key=filenum)
print(f'SwissTarget files: {len(files)}  | input compounds: {len(inp)}')
sw_rows = []
PROB = 0.0
for i, f in enumerate(files):
    if i >= len(inp): print('extra file (no matching compound):', os.path.basename(f)); continue
    comp = inp.iloc[i]
    herbs = str(comp['herbs']).split('|')
    df = pd.read_csv(f)
    gcol = 'Common name' if 'Common name' in df.columns else df.columns[1]
    pcol = next((c for c in df.columns if 'probab' in c.lower()), None)
    sub = df if pcol is None else df[pd.to_numeric(df[pcol], errors='coerce').fillna(0) > PROB]
    for g in sub[gcol].dropna().map(U).unique():
        for h in herbs:
            sw_rows.append([h, g, comp['MolName']])
sw_df = pd.DataFrame(sw_rows, columns=['herb', 'gene', 'compound']).drop_duplicates()
sw_df.to_csv(os.path.join(OUT, 'swiss_targets.csv'), index=False, encoding='utf-8-sig')
print('SwissTarget total unique genes (prob>0):', sw_df.gene.nunique())
print('per-source gene counts -> HERB:', herb_df.gene.nunique(), '| SwissTarget:', sw_df.gene.nunique())
