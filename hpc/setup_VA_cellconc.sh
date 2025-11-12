#!/bin/bash
#
#Script to setup simulations on the Habrok cluster
##launch lots of small va 
#
#NF template José Losa
#2025.10


# Store the input arguments in variables with understandable names:
# NAME=$1
CONDITION=$1
CPU=$2
HOURS=$3
MEMORY=$4
NRXNS=$5
OUTPUT_FOLDER=$6
BNAME=$7 
##batchname like __modxx_
# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f todolist.sh
touch todolist.sh

SLURM_MINS=$((MINS+1)) #add 1 min to job to allow gurobi to get to time limit and save .sol file 


# Setup batch scripts for each optimization:
for i in $(seq 0 $NRXNS)
do
#   FNAME=${NAME}_rxn${i}
    FNAME="202509_13C_${BNAME}_${CONDITION}_rxn${i}_fluxVA_"
    FILE="${HOME}/20250913C_files/2513C_VA_SCRIPTS/MPS_FOLDER/cellconc/${CONDITION}/${FNAME}.mps"
    echo "$FILE"
    # Change the settings for the job submission:
    cat ~/20250913C_files/2513C_VA_SCRIPTS/VA_batch_job_template.sh | sed "s/JOBNAME/$FNAME/g" | sed "s/CPUS/$CPU/g"  | sed "s/HOURS/$HOURS/g" | sed "s/MEMORY/$MEMORY/g" > batch_script_${FNAME}.slurm

    # Append line to call Gurobi:
    echo "gurobi.sh ~/20250913C_files/2513C_VA_SCRIPTS/gurobi_optimization_VA.py $FILE $FNAME $HOURS $CPU $OUTPUT_FOLDER" >> batch_script_${FNAME}.slurm

    # Append job to the 'todolist':
    echo "sbatch  batch_script_${FNAME}.slurm" >> todolist.sh


done

# Turn the 'todolist' into an executable:
chmod +x todolist.sh
