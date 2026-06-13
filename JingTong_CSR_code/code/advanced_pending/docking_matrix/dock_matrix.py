# -*- coding: utf-8 -*-
"""
Drug VIRTUAL SCREENING (focused, on-formula) -- AutoDock Vina docking MATRIX:
44 Jingtong-Granules active compounds  x  5 hub receptors (AKT1/IL1B/IL6/MMP9/TP53),
reusing the corrected structures + grid boxes in ../../docking_redo. Produces an
affinity matrix (kcal/mol) + heatmap to systematically rank compound-target pairs.

Run on a server (Vina + Open Babel). Receptor *.pdbqt must be prepared first (see
docking_redo/run_docking.sh, which builds <HUB>_receptor.pdbqt).

Usage:  py dock_matrix.py
Env:    AutoDock Vina (vina or vina.exe on PATH), Open Babel (obabel), rdkit, pandas
"""
import os, re, subprocess, shutil, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import figstyle as F; F.apply()

HERE = os.path.dirname(os.path.abspath(__file__))
EXP  = os.path.dirname(HERE)
REDO = os.path.join(os.path.dirname(EXP), 'docking_redo')      # corrected receptors + configs
LIGDIR = os.path.join(HERE, 'ligands_pdbqt'); os.makedirs(LIGDIR, exist_ok=True)
VINA = shutil.which('vina') or shutil.which('vina.exe') or os.path.join(os.path.dirname(EXP), 'docking', 'vina.exe')
OBABEL = shutil.which('obabel') or 'obabel'
HUBS = ['AKT1', 'IL1B', 'IL6', 'MMP9', 'TP53']                # receptors with structures+configs

def parse_cfg(hub):
    cfg = {}
    with open(os.path.join(REDO, f'config_{hub}.txt'), encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            m = re.match(r'\s*(center_[xyz]|size_[xyz])\s*=\s*([-\d.]+)', ln)
            if m: cfg[m.group(1)] = float(m.group(2))
    return cfg

def smi_to_pdbqt(name, smi):
    out = os.path.join(LIGDIR, f'{name}.pdbqt')
    if os.path.exists(out): return out
    # 3D + Gasteiger + pdbqt via Open Babel (RDKit alternative noted in README)
    cmd = [OBABEL, f'-:{smi}', '-opdbqt', '-O', out, '--gen3d', '--partialcharge', 'gasteiger']
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=180)
        return out if os.path.exists(out) and os.path.getsize(out) > 0 else None
    except Exception as e:
        print('  ligand prep failed', name, e); return None

def dock(receptor_pdbqt, ligand_pdbqt, cfg):
    out = ligand_pdbqt.replace('.pdbqt', '_out.pdbqt')
    cmd = [VINA, '--receptor', receptor_pdbqt, '--ligand', ligand_pdbqt, '--out', out,
           '--center_x', cfg['center_x'], '--center_y', cfg['center_y'], '--center_z', cfg['center_z'],
           '--size_x', cfg['size_x'], '--size_y', cfg['size_y'], '--size_z', cfg['size_z'],
           '--exhaustiveness', 16, '--num_modes', 5]
    cmd = [str(x) for x in cmd]
    try:
        r = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1200)
        for ln in (open(out).read().splitlines() if os.path.exists(out) else r.stdout.splitlines()):
            m = re.search(r'REMARK VINA RESULT:\s*([-\d.]+)', ln)
            if m: return float(m.group(1))
        for ln in r.stdout.splitlines():               # fallback: parse stdout table
            m = re.match(r'\s*1\s+([-\d.]+)', ln)
            if m: return float(m.group(1))
    except Exception as e:
        print('  dock failed', os.path.basename(ligand_pdbqt), e)
    return np.nan

# ---- ligands ----
lig = pd.read_csv(os.path.join(EXP, 'swisstarget_input.csv'))
lig = lig[['MolName', 'SMILES']].dropna()
lig['lid'] = lig['MolName'].str.replace(r'[^A-Za-z0-9]', '_', regex=True).str[:24]
print(f'{len(lig)} ligands x {len(HUBS)} receptors')
cfgs = {h: parse_cfg(h) for h in HUBS}

mat = pd.DataFrame(index=lig['lid'], columns=HUBS, dtype=float)
for _, r in lig.iterrows():
    lp = smi_to_pdbqt(r['lid'], r['SMILES'])
    if not lp: continue
    for h in HUBS:
        rec = os.path.join(REDO, f'{h}_receptor.pdbqt')
        if not os.path.exists(rec): print('  missing receptor', rec); continue
        mat.loc[r['lid'], h] = dock(rec, lp, cfgs[h])
    print(f"  {r['lid']:24s} " + ' '.join(f'{h}={mat.loc[r["lid"],h]}' for h in HUBS))

mat.to_csv(os.path.join(HERE, 'docking_affinity_matrix.csv'), encoding='utf-8-sig')

# ---- heatmap (more negative = stronger; darker = stronger) ----
m = mat.astype(float).dropna(how='all')
m = m.loc[m.mean(axis=1).sort_values().index]                 # strongest binders on top
fig, ax = plt.subplots(figsize=(4.6, max(4, 0.22 * len(m))))
sns.heatmap(m, cmap='RdYlBu', center=-6, annot=True, fmt='.1f', annot_kws={'size': 5.5},
            linewidths=0.4, linecolor='white', cbar_kws={'label': 'Vina affinity (kcal/mol)', 'shrink': 0.4}, ax=ax)
ax.set_xlabel('Hub target'); ax.set_ylabel('Active compound')
plt.setp(ax.get_yticklabels(), fontsize=5.5); plt.setp(ax.get_xticklabels(), fontsize=8)
ax.set_title('Compound–target docking matrix (AutoDock Vina)', fontsize=9)
fig.tight_layout(); F.save(fig, os.path.join(EXP, 'figures', 'Fig_docking_matrix'))
print('SAVED docking_affinity_matrix.csv + Fig_docking_matrix')
print('strongest pairs (<= -8 kcal/mol):')
s = m.stack(); print(s[s <= -8].sort_values().head(20).to_string())
