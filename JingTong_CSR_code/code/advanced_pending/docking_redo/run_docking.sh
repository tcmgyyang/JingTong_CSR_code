#!/bin/bash
# =============================================================================
# run_docking.sh — AutoDock Vina docking of quercetin & kaempferol to the
# corrected/added targets (AKT1=4EKL, IL6=1ALU, MMP9=1GKC, TP53=2OCJ).
# Prereq: *_receptor.pdbqt + quercetin.pdbqt + kaempferol.pdbqt prepared (see README §3).
# Configs (grid boxes) are pre-computed in config_<TARGET>.txt.
# Usage: bash run_docking.sh
# =============================================================================
set -uo pipefail
VINA=${VINA:-vina}     # path to AutoDock Vina (your docking/vina.exe on Windows, or 'vina')
LIGS=(quercetin kaempferol)
TARGETS=(AKT1 IL6 MMP9 TP53)
mkdir -p results

for T in "${TARGETS[@]}"; do
  for L in "${LIGS[@]}"; do
    echo ">>> docking $L -> $T"
    $VINA --config "config_${T}.txt" \
          --receptor "${T}_receptor.pdbqt" \
          --ligand   "${L}.pdbqt" \
          --out      "results/${T}_${L}_out.pdbqt" \
          --log      "results/${T}_${L}.log"
    # best affinity (mode 1) from log
    grep -A3 "^-----" "results/${T}_${L}.log" | awk 'NR==2{print "    best affinity:",$2,"kcal/mol"}'
  done
done

# ---- re-docking positive control (RMSD<2A) for targets with a co-crystal ligand ----
# Convert *_ref_ligand.pdb -> pdbqt, dock back into same box, compare to crystal pose.
for T in AKT1 MMP9; do
  if [ -f "${T}_ref_ligand.pdb" ]; then
    echo ">>> redocking control for $T"
    obabel "${T}_ref_ligand.pdb" -O "${T}_ref.pdbqt" -p 7.4 --partialcharge gasteiger 2>/dev/null || true
    $VINA --config "config_${T}.txt" --receptor "${T}_receptor.pdbqt" \
          --ligand "${T}_ref.pdbqt" --out "results/${T}_redock.pdbqt" --log "results/${T}_redock.log"
    # RMSD vs crystal pose (needs obrms from Open Babel)
    obabel "results/${T}_redock.pdbqt" -O "results/${T}_redock_best.pdb" -f 1 -l 1 2>/dev/null || true
    obrms "${T}_ref_ligand.pdb" "results/${T}_redock_best.pdb" 2>/dev/null \
      && echo "    ^ redocking RMSD (want < 2.0 A)" || echo "    (install obrms to report RMSD)"
  fi
done
echo "DONE. Update heat.csv with best affinities; run binding_residues.py on the best poses."
