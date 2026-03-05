#!/bin/bash
#
#Script to setup simulations on the cluster
##launch lots of small va 
#


# Store the input arguments in variables with understandable names:
CONDITION=$1
CPU=$2
HOURS=$3
MEMORY="${4}G"
INPUT_DIR=$5
OUTPUT_FOLDER=$6
SCRIPT_FOLDER=$7
NAME=$8
MIPGAP=$9
CHECK_FOLDER=${10}

# Setup an empty 'todolist' to store all the jobs to be submitted:
rm -f "${SCRIPT_FOLDER}/todolist.sh"
touch "${SCRIPT_FOLDER}/todolist.sh"

SLURM_MINS=$((MINS+1)) #add 1 min to job to allow gurobi to get to time limit and save .sol file 

# Setup batch scripts for each optimization:
for FILE in "${INPUT_DIR}"/*.mps.gz;
do
	[[ -f "$FILE" ]] || { echo "No .mps files found in $INPUT_DIR"; exit 1; }
	
    FULL_BASENAME=$(basename "$FILE")
    FNAME="${FULL_BASENAME%.mps}"
    REACTION_INDEX=$(echo "$FULL_BASENAME" | cut -d'_' -f1)
    
    if ls "${CHECK_FOLDER}/${REACTION_INDEX}"_* >/dev/null 2>&1; then
        echo "Skipping $REACTION_INDEX: Output already exists in $CHECK_FOLDER"
        continue
    fi
	
	  echo "Processing: $FILE -> FNAME=$FNAME"
	
    echo "$FILE"
    # Change the settings for the job submission:
    cat ~/hpc/VA_batch_job_template.sh | sed "s/JOBNAME/${NAME}_${FNAME}/g" | sed "s/CPUS/$CPU/g"  | sed "s/HOURS/$HOURS/g" | sed "s/MEMORY/$MEMORY/g" > "${SCRIPT_FOLDER}/batch_script_${FNAME}.slurm"

    # Append line to call Gurobi:
    echo "python ~/hpc/gurobi_optimization_VA.py $FILE $FNAME $HOURS $CPU $OUTPUT_FOLDER $MIPGAP" >> "${SCRIPT_FOLDER}/batch_script_${FNAME}.slurm"

    # Append job to the 'todolist':
    echo "sbatch  batch_script_${FNAME}.slurm" >> "${SCRIPT_FOLDER}/todolist.sh"
	
done

# Turn the 'todolist' into an executable:
# chmod +x todolist.sh
