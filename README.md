# NMOC SBayesRC analysis scripts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21796787.svg)](https://doi.org/10.5281/zenodo.21796787)

This repository contains the scripts used for SBayesRC-based polygenic risk
score analyses in the non-mucinous ovarian cancer (NMOC) study.

The main workflow includes SBayesRC model fitting, PLINK2-based PRS
calculation, preparation of validation datasets, and evaluation of PRS
performance. Additional scripts are included for PRS percentile-bin and
functional annotation enrichment analyses.

Individual-level genotype and phenotype data, intermediate files, and
study-specific annotation source files are not included.

## Repository structure

```text
scripts/
├── step1_sbayesrc/
├── step2_prs_calculation/
├── step3_validation_data/
├── step4_validation_statistics/
└── additional_analyses/
```
## Study design and input preparation

### Five-fold GWAS and cross-validation design

The primary cohort was divided into five mutually exclusive subsets, denoted
A, B, C, D, and E. In each cross-validation fold, four subsets were combined
to form the training dataset, representing approximately 80% of the cohort.
The remaining subset was retained as an independent testing dataset,
representing approximately 20% of the cohort.

The five training and testing configurations were:

| Fold | Training subsets | Testing subset |
|:----:|------------------|:--------------:|
| 1 | A, B, C, D | E |
| 2 | A, B, C, E | D |
| 3 | A, B, D, E | C |
| 4 | A, C, D, E | B |
| 5 | B, C, D, E | A |

For each fold, GWAS summary statistics were generated de novo using only the
individuals in the four training subsets. These fold-specific GWAS summary
statistics were used as input for SBayesRC model fitting. The resulting SNP
weights were then used to calculate and evaluate PRS in the held-out testing
subset.

Each participant therefore contributed to the testing dataset once and to the
training dataset in the remaining four folds. Identical training and testing
assignments were used across all annotation models and SNP-density panels to
ensure direct comparisons between models.

The individual-level fold assignment file is not included because it contains
study participant identifiers.

### Functional annotation preparation

Study-specific functional annotations were prepared using the SNP identifiers
and genomic positions provided in the SBayesRC `snp.info` files.

For interval-based annotations, including epigenomic regions, active
susceptibility-gene regions, and transcription factor binding sites, the SNP
positions were represented as single-base genomic intervals. The SNP intervals
and annotation BED files were compared using `bedtools intersect`.

For each interval-based annotation:

- a SNP was assigned a value of `1` when its genomic position overlapped an
  annotation interval;
- a SNP was assigned a value of `0` when no overlap was observed; and
- the resulting binary values formed one annotation column in the
  SNP-by-annotation matrix.

Precomputed variant-level quantitative annotations, such as annotation scores,
were mapped directly to matching SNP identifiers or genomic positions rather
than converted to binary values.

All SNP and annotation files were required to use a consistent genome build,
chromosome naming convention, and coordinate convention.

After annotation assignment, the study-specific annotation matrix was joined
to the SBayesRC BaselineLD v2.2 annotation matrix by SNP identifier using
`xan join`. The duplicate SNP identifier column introduced by the join was
removed, and the final matrix was written in the tab-delimited format required
by SBayesRC.

The annotation source datasets are described in the manuscript and
supplementary materials. Source BED files, intermediate overlap files, and
final annotation matrices are not distributed in this repository.

## Main workflow: Steps 1–4

### Step 1: SBayesRC analysis

Directory: [`scripts/step1_sbayesrc/`](scripts/step1_sbayesrc/)

- `step1_SBayesRC_v2.r` runs the SBayesRC `tidy`, `impute`, and `sbayesrc`
  functions for one cross-validation group and LD reference panel. It then adds
  chromosome, position, and allele information to the SBayesRC weight file for
  subsequent PLINK2 scoring.
- `submit_step1_SBayesRC_v2_single.sh` submits selected analyses as individual Slurm jobs.
- `submit_step1_SBayesRC_v2_array.sh` submits multiple runs simultaneously.

The R script accepts:

```text
Rscript step1_SBayesRC_v2.r GROUP_NUM LD_FOLDER WITH_ANNOTATION [ANNOTATION_FILE]
```

The functional annotation matrix was prepared separately as a study-specific input
in the format required by SBayesRC. Annotation-preparation scripts and source
files are not included in this repository.

### Step 2: PRS calculation

Directory:
[`scripts/step2_prs_calculation/`](scripts/step2_prs_calculation/)

- `step2_prs_calculation_single.sh` calculates PRS for SBayesRC weight file as a single job
- `step2_prs_calculation_array.sh` calculates PRS for one SBayesRC weight file per Slurm
  array task.

For each chromosome, PLINK2 calculates centered scores using columns
1–3 of the SBayesRC weight file:

```text
variant ID, effect allele, SNP weight
```

The chromosome-level scores are joined by `IID` and summed to produce
`score_file/prs_result_sum.txt`.

### Step 3: Prepare validation data

Directory:
[`scripts/step3_validation_data/`](scripts/step3_validation_data/)

- `step3_file_modification.py` joins the phenotype table, summed PRS, and
  `groups_info.txt`. It selects NMOC cases (`all_non_mucinous == 1`) and
  controls (`all_non_mucinous == 0`) from the requested validation group and
  writes separate case and control files.

The script accepts:

```text
python step3_file_modification.py OUTPUT_FOLDER GROUP_NUM
```

### Step 4: Validation statistics

Directory:
[`scripts/step4_validation_statistics/`](scripts/step4_validation_statistics/)

- `step4_stat_calculation_v3.r` calculates AUC with a DeLong 95% interval, the
  likelihood-ratio chi-square comparing PRS and intercept-only logistic models,
  and the log odds ratio per standard deviation of PRS with its 95% interval.
- `submit_step3_step4_single.sh` runs Steps 3 and 4 for one result
  folder/group pair.
- `submit_step3_step4_v2_array.sh` runs Steps 3 and 4 as a Slurm array and
  uses `flock` when aggregating result rows.

The R script accepts:

```text
Rscript step4_stat_calculation_v3.r OUTPUT_FOLDER GROUP_NUM
```

## Additional analyses

Directory:
[`scripts/additional_analyses/`](scripts/additional_analyses/)

### PRS percentile-bin analysis

- `prs_percentile_bin_analysis.R` defines PRS percentile cutoffs from controls,
  uses the 40–60% group as the reference, and estimates odds ratios, 95%
  intervals, and p-values for the remaining percentile groups.
- `submit_prs_percentile_bin_analysis.sh` is the corresponding Slurm launcher.

### Annotation enrichment analysis

- `annotation_enrichment_analysis.R` summarizes SBayesRC per-SNP heritability
  enrichment across five cross-validation groups and the HapMap3 and Imputed
  reference panels. It writes annotation-level summary tables and the corresponding
  supplementary figures.
- `submit_annotation_enrichment_analysis.sh` is the corresponding Slurm
  launcher.

For the annotation-specific figures, error bars represent the minimum and
maximum fold-specific enrichment point estimates. They are not confidence
intervals or posterior credible intervals.

## Study-specific paths and dependencies

The scripts contain the HPC paths, Conda environment names, resource
requests, and file naming conventions used for the study. These paths must be
updated if the scripts are run in another environment.

Required software and packages include:

- R and the SBayesRC package
- R packages: `tidyverse`, `pROC`, `lmtest`, `patchwork`, and `scales`
- Python with `pandas`
- PLINK2
- BEDTools
- `xan`
- Slurm

## Data availability

Individual-level genotype and phenotype data are not included because they are
controlled-access data. Participant-level fold assignments are also not
distributed because they contain study identifiers.

Functional annotation sources are listed in the manuscript and supplementary
materials. The corresponding source files, intermediate annotation files, and
final SBayesRC annotation matrices are not included.

The SBayesRC software is available from its original repository:

https://github.com/zhilizheng/SBayesRC

## Archived release

Version `v1.0.0` of this repository is archived on Zenodo:

**DOI:** [10.5281/zenodo.21796787](https://doi.org/10.5281/zenodo.21796787)
