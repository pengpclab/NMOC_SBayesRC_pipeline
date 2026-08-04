#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm launcher for the PRS percentile-bin analysis. It passes one PRS result
# folder and group number to prs_percentile_bin_analysis.R. Scheduler resources,
# environment activation, and analysis arguments are unchanged; the script path
# reflects the cleaned repository layout and renamed analysis file.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -p defq
#SBATCH -J PRS_BINS
#SBATCH -t 12:00:00
#SBATCH --mem=20000MB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=your_email_address
#SBATCH -o ./PRS_PERCENTILE_BIN_LOG/PRS_BINS_%j.out
#SBATCH -e ./PRS_PERCENTILE_BIN_LOG/PRS_BINS_%j.err

#sbatch submit_prs_percentile_bin_analysis.sh /path/to/SBayesRC/result/prs_result_NMOC_GROUP_0_HapMap3_sbrc_annotation 0

mkdir -p ./PRS_PERCENTILE_BIN_LOG

echo "===== JOB STARTED ====="
echo "Date       : $(date)"
echo "Job ID     : $SLURM_JOB_ID"
echo "Node       : $(hostname)"
echo "Working Dir: $(pwd)"
echo "========================"

  if [ $# -ne 2 ]; then
      echo "Error: Must provide output_folder and group number"
      echo "Usage: sbatch $0 <output_folder> <group_num>"
      exit 1
  fi

  OUTPUT_FOLDER=$1
  GROUP_NUM=$2

  source /path/to/miniconda/etc/profile.d/conda.sh
  conda activate r_env

  echo "Running PRS percentile-bin analysis"
  Rscript /path/to/SBayesRC/scripts/additional_analyses/prs_percentile_bin_analysis.R "$OUTPUT_FOLDER" "$GROUP_NUM"
  if [ $? -ne 0 ]; then
      echo "Error: PRS percentile-bin analysis (R script) failed"
      conda deactivate
      exit 1
  fi
  conda deactivate

  echo "===== JOB FINISHED $(date) ====="
