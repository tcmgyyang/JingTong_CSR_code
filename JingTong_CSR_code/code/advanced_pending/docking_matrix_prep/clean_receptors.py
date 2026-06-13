# -*- coding: utf-8 -*-
"""Tool-free receptor cleaning: for each hub PDB keep the single protein chain whose atoms
sit in/near the docking box (from config_<HUB>.txt), drop waters/heteroatoms/altlocs.
Output <HUB>_receptor_clean.pdb -> ready for a 1-command obabel/prepare_receptor step on
the server (see prepare_receptors.sh). No external tools needed here."""
import os, re, math
REDO = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'docking_redo')
OUT  = os.path.dirname(os.path.abspath(__file__))
PDB  = {'AKT1': '4EKL', 'IL1B': '8RYS', 'IL6': '1ALU', 'MMP9': '1GKC', 'TP53': '2OCJ'}
STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE'.split())

def box_center(hub):
    c = {}
    for ln in open(os.path.join(REDO, f'config_{hub}.txt'), encoding='utf-8', errors='replace'):
        m = re.match(r'\s*(center_[xyz])\s*=\s*([-\d.]+)', ln)
        if m: c[m.group(1)] = float(m.group(2))
    return (c['center_x'], c['center_y'], c['center_z'])

for hub, pdbid in PDB.items():
    src = os.path.join(REDO, f'{pdbid}.pdb')
    if not os.path.exists(src):
        print('MISSING', src); continue
    cx, cy, cz = box_center(hub)
    lines = open(src, encoding='utf-8', errors='replace').read().splitlines()
    # score each chain by #protein atoms within 16 A of the box centre
    chain_hits, chain_atoms = {}, {}
    for ln in lines:
        if ln.startswith('ATOM') and ln[17:20].strip() in STD:
            ch = ln[21]
            try: x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            except ValueError: continue
            chain_atoms.setdefault(ch, []).append(ln)
            if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= 16**2:
                chain_hits[ch] = chain_hits.get(ch, 0) + 1
    if not chain_atoms:
        print(hub, 'no protein atoms?!'); continue
    keep = max(chain_hits or {c: len(a) for c, a in chain_atoms.items()},
               key=lambda c: chain_hits.get(c, 0))
    # write kept chain, primary altloc only
    out = os.path.join(OUT, f'{hub}_receptor_clean.pdb')
    n = 0
    with open(out, 'w') as fh:
        for ln in chain_atoms[keep]:
            alt = ln[16]
            if alt not in (' ', 'A'): continue
            fh.write(ln + '\n'); n += 1
        fh.write('TER\nEND\n')
    print(f'{hub} ({pdbid}) chain {keep}: {n} atoms (box hits {chain_hits.get(keep,0)}) -> {os.path.basename(out)}')
print('\nDone. Next (on server): bash prepare_receptors.sh  ->  <HUB>_receptor.pdbqt')
