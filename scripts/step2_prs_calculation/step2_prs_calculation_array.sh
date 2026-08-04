#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm array job for PRS calculation. Each task selects one SBayesRC weight
# file, runs PLINK2 scoring for every chromosome, joins chromosome-level scores
# by IID, and sums them into score_file/prs_result_sum.txt. The scoring columns,
# centering option, file naming, paths, and merge/sum operations are unchanged.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p defq
#SBATCH -J S2_PLINK
#SBATCH -t 48:00:00
#SBATCH --mem=100000MB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=your_email_address
#SBATCH -o ./PRS_STEP2_LOG/STEP2_%A_%a.out
#SBATCH -e ./PRS_STEP2_LOG/STEP2_%A_%a.err
#SBATCH --array=1-10                         # Remember to modify it based on the number of lines in input file



# Usage:
#sbatch /path/to/SBayesRC/scripts/step2_prs_calculation/step2_prs_calculation_array.sh <WEIGHT_LIST_FILE.txt>
#sbatch /path/to/SBayesRC/scripts/step2_prs_calculation/step2_prs_calculation_array.sh /path/to/SBayesRC/data/step2_AllAnnotation_input_file.txt


#
# Step 0 Create log directory if it doesn't exist
#
mkdir -p ./PRS_STEP2_LOG

echo "===== JOB STARTED ====="
echo "Date         : $(date)"
echo "Array Job ID : ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Array Task   : $SLURM_ARRAY_TASK_ID"
echo "Node         : $(hostname)"
echo "Working Dir  : $(pwd)"
echo "========================"


#
# Step 1 Setup
#
if [ $# -lt 1 ]; then
    echo "Usage: $0 WEIGHT_LIST_FILE"
    exit 1
fi

WEIGHT_LIST=$1
WEIGHT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${WEIGHT_LIST})
WEIGHT_BASENAME=$(basename "$WEIGHT_FILE" _ID.txt)

echo "Task ${SLURM_ARRAY_TASK_ID} using weight file: $WEIGHT_FILE"

source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate linux_env


#
# Step 2 Directories and resources
#
PGEN_DIR="/path/to/EOC_Data/Genotype_PGEN"
WORKING_DIR="/path/to/SBayesRC/data"
OUT_DIR="/path/to/SBayesRC/result"
WEIGHT_DIR="/path/to/SBayesRC/data/SBayesRC_PGM"

OUT_FOLDER="prs_result_${WEIGHT_BASENAME}"
mkdir -p ${OUT_DIR}/${OUT_FOLDER}
echo "Output folder: $OUT_FOLDER"
SCORE_FOLDER="${OUT_DIR}/${OUT_FOLDER}/score_file"
mkdir -p ${SCORE_FOLDER}

THREADS=12
MEMORY=800000


#
# Step 3 Calculate PRS per chromosome
#
for DATA in ${PGEN_DIR}/ALL_CHR*.pvar;
do
	echo "--------------------START ONE CHROMOSOME--------------------"
	PLINK_FILE=$(echo ${DATA} | sed 's/.pvar$//')
	CHR=$(basename ${PLINK_FILE} | cut -d'_' -f2)
	RESULT_FILE=$(echo ${OUT_DIR}/${OUT_FOLDER}/${OUT_FOLDER}_${CHR})

	echo "Running: ${PLINK_FILE}, Chr: ${CHR}"

	# RUN PLINK SCORE FUNCTION
	plink2 \
	--pfile ${PLINK_FILE} \
	--score ${WEIGHT_DIR}/${WEIGHT_FILE} 1 2 3 header center list-variants \
	--out ${RESULT_FILE} \
	--threads ${THREADS} \
	--memory ${MEMORY}

	echo "--------------------FINISH ONE CHROMOSOME--------------------"
done

mv ${OUT_DIR}/${OUT_FOLDER}/*.sscore ${SCORE_FOLDER}


#
# Step 4 Merge all the prs data from each chromosomes into one file
#
FILE_NUM=0

for PROFILE in ${SCORE_FOLDER}/*.sscore;
do
	FILE_NUM=$(( FILE_NUM+1 ))
	TITLE=$(basename ${PROFILE} | cut -d '.' -f 1 | rev | cut -d '_' -f -1 | rev)
	if [ ${FILE_NUM} == 1 ]
	then
		cut -f 2,5 ${PROFILE} | tail -n +2 > ${OUT_DIR}/${OUT_FOLDER}/final.txt
		sed -i "1i IID $TITLE" ${OUT_DIR}/${OUT_FOLDER}/final.txt

	else
		cut -f 2,5 ${PROFILE} | tail -n +2 > ${OUT_DIR}/${OUT_FOLDER}/add_column.txt
		sed -i "1i IID $TITLE" ${OUT_DIR}/${OUT_FOLDER}/add_column.txt
		join ${OUT_DIR}/${OUT_FOLDER}/final.txt ${OUT_DIR}/${OUT_FOLDER}/add_column.txt > ${OUT_DIR}/${OUT_FOLDER}/tmp && mv ${OUT_DIR}/${OUT_FOLDER}/tmp ${OUT_DIR}/${OUT_FOLDER}/final.txt
	fi
done

mv ${OUT_DIR}/${OUT_FOLDER}/final.txt ${SCORE_FOLDER}/prs_result_ALL.txt
rm ${OUT_DIR}/${OUT_FOLDER}/add_column.txt


#
# Step 5 Sum up the PRS from each chromosomes and extract the result
#

old_filename=${SCORE_FOLDER}/prs_result_ALL.txt
new_filename="${old_filename/_ALL.txt/_sum.txt}"

awk '
{
	if (NR == 1) {
		print $1, "prs_score"
	} else {
		sum = 0
		for (i = 2; i <= NF; i++) {
			sum += $i
		}
		print $1, sum
	}
}' "$old_filename" > "$new_filename"


conda deactivate
echo "===== JOB FINISHED $(date) ====="
