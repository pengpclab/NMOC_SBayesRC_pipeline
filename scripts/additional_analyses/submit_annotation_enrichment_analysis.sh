#!/bin/bash
# PUBLIC RELEASE DOCUMENTATION
# Slurm launcher for the annotation enrichment summary and figure analysis. It
# supplies the study SBayesRC data directory and NMOC_enrichment_summary output
# directory to annotation_enrichment_analysis.R. Resources, environment, and
# command-line arguments are unchanged; the path reflects the renamed file.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p defq
#SBATCH -J SBayesRC_Enrich
#SBATCH -t 24:00:00
#SBATCH --mem=100000MB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=your_email_address
#SBATCH -o ./SBayesRC_Enrich_Log/Enrich_%j.out
#SBATCH -e ./SBayesRC_Enrich_Log/Enrich_%j.err


# Create log directory if it doesn't exist
mkdir -p ./SBayesRC_Enrich_Log


# Activate conda environment
source /path/to/anaconda3/etc/profile.d/conda.sh


# Run Enrichment Analysis (R script)
echo "Running Enrichment Analysis"
conda activate r_env
Rscript /path/to/SBayesRC/scripts/additional_analyses/annotation_enrichment_analysis.R /path/to/SBayesRC/data /path/to/SBayesRC/result/NMOC_enrichment_summary
if [ $? -ne 0 ]; then
    echo "Error: Enrichment Analysis failed"
    conda deactivate
    exit 1
fi
conda deactivate

echo "Analysis completed successfully"
