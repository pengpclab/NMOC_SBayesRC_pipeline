#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm launcher for explicitly selected step-1 SBayesRC runs. In the deposited
# study version, the active commands run annot_ALL.txt for both HapMap3 and
# Imputed LD panels with group 0; the commented commands record earlier targeted
# annotation runs. Analysis parameters are unchanged; script paths reflect the
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
#SBATCH -o ./PRS_STEP1_LOG/STEP1_%j.out
#SBATCH -e ./PRS_STEP1_LOG/STEP1_%j.err


# Create log directory if it doesn't exist
mkdir -p ./PRS_STEP1_LOG


# Load conda environment
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env


# Run the R script with those arguments
#Rscript /path/to/SBayesRC/scripts/step1_sbayesrc/step1_SBayesRC_v2.r 0 ukbEUR_Imputed TRUE annotation_file.txt
Rscript /path/to/SBayesRC/scripts/step1_sbayesrc/step1_SBayesRC_v2.r 0 ukbEUR_HM3 TRUE annotation_file.txt

conda deactivate


# Note: DO NOT USE ONLY 12 CPU, CAUSE CORE DUMP
