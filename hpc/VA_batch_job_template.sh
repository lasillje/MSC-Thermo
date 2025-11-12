#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --cpus-per-task=CPUS
#SBATCH --time=HOURS:02:00
#SBATCH --partition=regular
#SBATCH --mem-per-cpu=MEMORY

## Load Python and activate environment
module purge
module load Python/3.10.8-GCCcore-12.2.0
module load Gurobi
source ~/venvs/tf_env/bin/activate

##lines to be appended
