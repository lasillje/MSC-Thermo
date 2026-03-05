#!/bin/bash
#SBATCH --job-name=TFS_DRG_FLUXES
#SBATCH --cpus-per-task=16
#SBATCH --time=48:05:00
#SBATCH --partition=regular
#SBATCH --mem-per-cpu=1G


IN_DIR="$1"
IN_POINTS="$2"
OUT_DIR="$3"
NUM_SAMPLES_IN="$4"
NUM_SAMPLES_OUT="$5"

module purge

CONDA_BASE=$(conda info --base)

echo "Conda base: $CONDA_BASE"

if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

conda activate /home5/s3997200/.conda/envs/pta-clone/

python tfs_drg_and_fluxes.py "$IN_DIR" "$IN_POINTS" "$OUT_DIR" "$NUM_SAMPLES_IN" "$NUM_SAMPLES_OUT"