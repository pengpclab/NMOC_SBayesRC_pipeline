#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm array launcher for step1_SBayesRC_v2.r. Each array task reads one line
# from step1_AllAnnotation_input_file.txt and passes its four whitespace-
# separated fields to the R script. The --array range must equal the number of
# input lines. Analysis parameters are unchanged; script paths reflect the
# cleaned repository layout.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p defq
#SBATCH -J S1_SBayesRC
#SBATCH -t 24:00:00
#SBATCH --mem=100000MB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=your_email_address
#SBATCH -o ./PRS_STEP1_LOG/STEP1_%A_%a.out
#SBATCH -e ./PRS_STEP1_LOG/STEP1_%A_%a.err
#SBATCH --array=1-10                         # Remember to modify it based on the number of lines in input file


# Create log directory if it doesn't exist
mkdir -p ./PRS_STEP1_LOG

echo "===== JOB STARTED ====="
echo "Date         : $(date)"
echo "Array Job ID : ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Array Task   : $SLURM_ARRAY_TASK_ID"
echo "Node         : $(hostname)"
echo "Working Dir  : $(pwd)"
echo "========================"


# Load conda environment
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env


# Read the input file and split the line into arguments
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/SBayesRC/data/step1_AllAnnotation_input_file.txt)  #Remember to modify the input file
set -- $LINE
echo "Parsed input line: $LINE"


# Run the R script with those arguments
Rscript /path/to/SBayesRC/scripts/step1_sbayesrc/step1_SBayesRC_v2.r "$1" "$2" "$3" "$4"
conda deactivate
echo "===== JOB FINISHED $(date) ====="







# Note: DO NOT USE ONLY 12 CPU, CAUSE CORE DUMP
#Rscript /path/to/SBayesRC/scripts/step1_sbayesrc/step1_SBayesRC_v2.r 0 ukbEUR_HM3 TRUE annotation_file.txt   #example row in "step1_AllAnnotation_input_file.txt"
