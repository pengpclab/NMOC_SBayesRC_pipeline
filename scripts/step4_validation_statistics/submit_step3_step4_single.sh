#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm launcher that runs held-out data preparation (Python step 3) followed by
# PRS performance statistics (R step 4) for one output folder/group pair. The
# sequential environment activation, failure handling, and arguments are
# unchanged; script paths reflect the cleaned repository layout.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p defq
#SBATCH -J S34_SBayesRC
#SBATCH -t 24:00:00
#SBATCH --mem=100000MB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=your_email_address
#SBATCH -o ./PRS_STEP3_STEP4_LOG/STEP34_%j.out
#SBATCH -e ./PRS_STEP3_STEP4_LOG/STEP34_%j.err

# Example usage:
# sbatch /path/to/SBayesRC/scripts/step4_validation_statistics/submit_step3_step4_single.sh prs_result_NMOC_GROUP_0_HapMap3_sbrc_annotation 0


# Create log directory if it doesn't exist
mkdir -p ./PRS_STEP3_STEP4_LOG


# Check if arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Must provide output_folder and group number"
    echo "Usage: sbatch $0 <output_folder> <group_num>"
    exit 1
fi

OUTPUT_FOLDER=$1
GROUP_NUM=$2


# Activate conda environment
source /path/to/anaconda3/etc/profile.d/conda.sh


# Run Step 3 (Python script)
echo "Running Step 3: Generating PRS case and control files"
conda activate python_env
python /path/to/SBayesRC/scripts/step3_validation_data/step3_file_modification.py "$OUTPUT_FOLDER" "$GROUP_NUM"
if [ $? -ne 0 ]; then
    echo "Error: Step 3 (Python script) failed"
    conda deactivate
    exit 1
fi
conda deactivate


# Run Step 4 (R script)
echo "Running Step 4: Statistics calculation"
conda activate r_env
Rscript /path/to/SBayesRC/scripts/step4_validation_statistics/step4_stat_calculation_v3.r "$OUTPUT_FOLDER" "$GROUP_NUM"
if [ $? -ne 0 ]; then
    echo "Error: Step 4 (R script) failed"
    conda deactivate
    exit 1
fi
conda deactivate

echo "Both steps completed successfully"
