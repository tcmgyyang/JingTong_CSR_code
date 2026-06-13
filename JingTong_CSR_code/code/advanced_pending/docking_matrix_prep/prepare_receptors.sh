#!/bin/bash
# =============================================================================
# prepare_receptors.sh  (run on the SERVER, where Open Babel / ADFR are available)
# Converts the tool-free-cleaned receptor PDBs -> AutoDock Vina rigid receptor .pdbqt
# Cleaned PDBs (chain A, protein only, box-validated) were produced locally by
# clean_receptors.py. Output <HUB>_receptor.pdbqt is what dock_matrix.py / run_docking.sh need.
#
# Usage:  bash prepare_receptors.sh        # then copy *_receptor.pdbqt into docking_redo/
# =============================================================================
set -uo pipefail
HUBS=(AKT1 IL1B IL6 MMP9 TP53)

for H in "${HUBS[@]}"; do
  IN="${H}_receptor_clean.pdb"
  OUT="${H}_receptor.pdbqt"
  [ -f "$IN" ] || { echo "skip $H (no $IN)"; continue; }
  echo ">>> $H : $IN -> $OUT"
  # --- Option A: Open Babel (rigid receptor, protonate at pH 7.4, Gasteiger) ---
  obabel "$IN" -O "$OUT" -xr -p 7.4 --partialcharge gasteiger 2>/dev/null \
    && echo "    OK (obabel)" && continue
  # --- Option B (fallback): AutoDockTools / ADFR ---
  #   prepare_receptor -r "$IN" -o "$OUT" -A hydrogens -U nphs_lps_waters
  echo "    obabel failed -> try: prepare_receptor -r $IN -o $OUT -A hydrogens -U nphs_lps_waters"
done

echo ""
echo "Done. Sanity check (each file should start with REMARK / ROOT and contain ATOM lines):"
for H in "${HUBS[@]}"; do
  [ -f "${H}_receptor.pdbqt" ] && echo "  ${H}_receptor.pdbqt : $(grep -c '^ATOM' "${H}_receptor.pdbqt") atoms"
done
echo "Then: cp *_receptor.pdbqt ../docking_redo/   (so dock_matrix.py finds them)"
