#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm array launcher for validation across multiple result-folder/group rows.
# It runs steps 3 and 4, then uses flock to serialize appends to the aggregate
# statistics CSV. Version 2 was used to prevent the race condition described at
# the end of this file. Analysis and aggregation commands are unchanged.
#SBATCH --job-name=S34_STAT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p defq
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email_address
#SBATCH --output=/path/to/SBayesRC/script/PRS_STEP3_STEP4_LOG/STEP34_%A_%a.out
#SBATCH --error=/path/to/SBayesRC/script/PRS_STEP3_STEP4_LOG/STEP34_%A_%a.err
#SBATCH --array=1-10

# Example usage:
# sbatch /path/to/SBayesRC/scripts/step4_validation_statistics/submit_step3_step4_v2_array.sh /path/to/SBayesRC/data/step3_input_file_HM3.txt
# sbatch /path/to/SBayesRC/scripts/step4_validation_statistics/submit_step3_step4_v2_array.sh /path/to/SBayesRC/data/step3_AllAnnotation_input_file.txt


# Create log directory if it doesn't exist
mkdir -p ./PRS_STEP3_STEP4_LOG

echo "===== JOB STARTED ====="
echo "Date         : $(date)"
echo "Array Job ID : ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Array Task   : $SLURM_ARRAY_TASK_ID"
echo "Node         : $(hostname)"
echo "Working Dir  : $(pwd)"
echo "========================"


echo "PRS result list is: " $1
RESULT_LIST=$1
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${RESULT_LIST})
set -- $LINE
echo "Parsed input line: $LINE"


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



# Export statistic results according to the input result list
OUTPUT_DIR="/path/to/SBayesRC/result"
RESULT_LIST_BASENAME=$(basename "$RESULT_LIST" .txt)
STAT_RESULT="${OUTPUT_DIR}/${RESULT_LIST_BASENAME}_stat_result.csv"
STAT_RESULT_SORTED="${OUTPUT_DIR}/${RESULT_LIST_BASENAME}_stat_result_sorted.csv"
LOCKFILE="${STAT_RESULT}.lock"


SRC_FILE="${OUTPUT_DIR}/${OUTPUT_FOLDER}/nmoc_grp${GROUP_NUM}_prs_results.csv"

if [[ ! -s "${SRC_FILE}" ]]; then
	echo "WARNING: source result missing or empty: ${SRC_FILE}" >&2
else
	{
		# Open (or create) the lock file on FD 200 and take an exclusive lock
		exec 200>"${LOCKFILE}"
		# Wait indefinitely
		flock -x 200

		if [[ ! -f "${STAT_RESULT}" ]]; then
			# Keep header + all rows
			cat "${SRC_FILE}" > "${STAT_RESULT}"
		else
			# Append only the last row
			tail -n 1 "${SRC_FILE}" >> "${STAT_RESULT}"
		fi

		# Release the lock
		flock -u 200
	} 2>> "${OUTPUT_DIR}/aggregation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
fi

cat ${STAT_RESULT} | sort > ${STAT_RESULT_SORTED}

echo "All the Statistics results are now stored in this single csv file: " ${STAT_RESULT_SORTED}
echo "===== JOB FINISHED $(date) ====="


# Note: In version one, due to race conditions, multiple tasks may reach the final step at the same time. Lines can be dropped, duplicated, or truncated. Version two had fixed this issue.
