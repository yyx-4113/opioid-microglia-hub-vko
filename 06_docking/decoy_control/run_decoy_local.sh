#!/bin/bash
# Local runner: downloaded Vina 1.2.7 windows binary; all paths local (D:)
# No set -e: a single failed docking must not abort the 630-job sweep.
ROOT="D:/2026.9/极速交付9月会员日优惠套路/05_多组学+虚拟敲除药物发现/方案二_阿片耐受与痛觉过敏小胶质枢纽基因虚拟敲除"
# vina called by absolute path inside each command (vina_cmds_abs.txt); PATH left for safety
export REC_PQ="$ROOT/06_docking/hpc_bundle_autodl/06_docking/rec_pdbqt"
export LIG_PQ="$ROOT/06_docking/hpc_bundle_autodl/06_docking/lig_pdbqt"
export DEC_PQ="$ROOT/06_docking/decoy_control/lig_decoy_pdbqt"
export OUTDIR="$ROOT/06_docking/decoy_control"
mkdir -p "$OUTDIR/results_realbox" "$OUTDIR/results_randbox" "$OUTDIR/results_decoy"
NP=${1:-8}
echo "Launching 630 vina jobs with -P $NP on local Vina 1.2.7 ($(date))"
cat "$OUTDIR/vina_cmds_abs.txt" | xargs -P "$NP" -I{} bash -c '{}' && true
echo "DONE_ALL_DECOY_RUNS $(date)"
