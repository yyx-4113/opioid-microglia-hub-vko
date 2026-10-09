#!/bin/bash
#SBATCH --job-name=lincs_vina
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=vina_%A.out
module load vina 2>/dev/null || true
bash $SLURM_SUBMIT_DIR/run_vina.sh
