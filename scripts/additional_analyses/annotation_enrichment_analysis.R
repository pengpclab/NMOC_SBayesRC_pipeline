#!/usr/bin/env Rscript

# PUBLIC RELEASE DOCUMENTATION
# This is the annotation enrichment and figure analysis used in the study,
# given a clear version-free name for the simplified repository. Analysis logic
# and command-line arguments are unchanged.

# SBayesRC enrichment analysis
#
# Summarize SBayesRC per-SNP heritability enrichment from five CV folds
# and compare HapMap3 versus Imputed reference panels.
#
#
# Default input structure:
# /path/to/SBayesRC/data/
#   NMOC_GROUP_0_HapMap3_sbrc_ALL/*.AnnoPerSnpHsqEnrichment
#   ...
#   NMOC_GROUP_4_Imputed_sbrc_ALL/*.AnnoPerSnpHsqEnrichment
#
# Usage:
#   Rscript annotation_enrichment_analysis.R
# or
#   Rscript annotation_enrichment_analysis.R /path/to/data /path/to/output

options(stringsAsFactors = FALSE)

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "forcats",
  "ggplot2", "patchwork", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, for example with conda:\n",
    "conda install -c conda-forge r-tidyverse r-patchwork r-scales"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "/path/to/SBayesRC/data"
out_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "NMOC_enrichment_summary")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure settings --------------------------------------------------------------
top_n_all <- 20
top_n_per_category <- 8

# Keep this empty to rank all annotations. Add a source name only when exclusion
# is scientifically pre-specified, not merely because its estimates are large.
exclude_sources_from_top20 <- character(0)

panel_labels <- c(
  "HapMap3" = "HapMap3 variants",
  "Imputed" = "Imputed variants"
)
panel_colors <- c(
  "HapMap3" = "#3B8EC2",
  "Imputed" = "#E69F00"
)

# Baseline-LD v2.2 annotation names used in the SBayesRC release.
# Explicit names are safer than assigning the first N rows to Baseline-LD.
baseline_ld_annotations <- c(
  "Coding_UCSC", "Coding_UCSC.flanking.500",
  "Conserved_LindbladToh", "Conserved_LindbladToh.flanking.500",
  "CTCF_Hoffman", "CTCF_Hoffman.flanking.500",
  "DGF_ENCODE", "DGF_ENCODE.flanking.500",
  "DHS_peaks_Trynka", "DHS_Trynka", "DHS_Trynka.flanking.500",
  "Enhancer_Andersson", "Enhancer_Andersson.flanking.500",
  "Enhancer_Hoffman", "Enhancer_Hoffman.flanking.500",
  "FetalDHS_Trynka", "FetalDHS_Trynka.flanking.500",
  "H3K27ac_Hnisz", "H3K27ac_Hnisz.flanking.500",
  "H3K27ac_PGC2", "H3K27ac_PGC2.flanking.500",
  "H3K4me1_peaks_Trynka", "H3K4me1_Trynka",
  "H3K4me1_Trynka.flanking.500",
  "H3K4me3_peaks_Trynka", "H3K4me3_Trynka",
  "H3K4me3_Trynka.flanking.500",
  "H3K9ac_peaks_Trynka", "H3K9ac_Trynka",
  "H3K9ac_Trynka.flanking.500",
  "Intron_UCSC", "Intron_UCSC.flanking.500",
  "PromoterFlanking_Hoffman", "PromoterFlanking_Hoffman.flanking.500",
  "Promoter_UCSC", "Promoter_UCSC.flanking.500",
  "Repressed_Hoffman", "Repressed_Hoffman.flanking.500",
  "SuperEnhancer_Hnisz", "SuperEnhancer_Hnisz.flanking.500",
  "TFBS_ENCODE", "TFBS_ENCODE.flanking.500",
  "Transcr_Hoffman", "Transcr_Hoffman.flanking.500",
  "TSS_Hoffman", "TSS_Hoffman.flanking.500",
  "UTR_3_UCSC", "UTR_3_UCSC.flanking.500",
  "UTR_5_UCSC", "UTR_5_UCSC.flanking.500",
  "WeakEnhancer_Hoffman", "WeakEnhancer_Hoffman.flanking.500",
  "GERP.NS", "GERP.RSsup4",
  paste0("MAFbin", 1:10),
  "MAF_Adj_Predicted_Allele_Age", "MAF_Adj_LLD_AFR",
  "Recomb_Rate_10kb", "Nucleotide_Diversity_10kb",
  "Backgrd_Selection_Stat", "CpG_Content_50kb", "MAF_Adj_ASMC",
  "GTEx_eQTL_MaxCPP", "BLUEPRINT_H3K27acQTL_MaxCPP",
  "BLUEPRINT_H3K4me1QTL_MaxCPP",
  "BLUEPRINT_DNA_methylation_MaxCPP",
  "synonymous", "non_synonymous",
  "Conserved_Vertebrate_phastCons46way",
  "Conserved_Vertebrate_phastCons46way.flanking.500",
  "Conserved_Mammal_phastCons46way",
  "Conserved_Mammal_phastCons46way.flanking.500",
  "Conserved_Primate_phastCons46way",
  "Conserved_Primate_phastCons46way.flanking.500",
  "BivFlnk", "BivFlnk.flanking.500",
  "Human_Promoter_Villar", "Human_Promoter_Villar.flanking.500",
  "Human_Enhancer_Villar", "Human_Enhancer_Villar.flanking.500",
  "Ancient_Sequence_Age_Human_Promoter",
  "Ancient_Sequence_Age_Human_Promoter.flanking.500",
  "Ancient_Sequence_Age_Human_Enhancer",
  "Ancient_Sequence_Age_Human_Enhancer.flanking.500",
  "Human_Enhancer_Villar_Species_Enhancer_Count",
  "Human_Promoter_Villar_ExAC",
  "Human_Promoter_Villar_ExAC.flanking.500"
)

# FAVOR annotations appended to BaselineLD v2.2 in annotation_baseline_FAVOR.
# These are grouped with Baseline in Supplementary Table 1.
favor_annotations <- c(
  "apc_conservation_v2",
  "apc_epigenetics_active",
  "apc_epigenetics_repressed",
  "apc_epigenetics_transcription",
  "apc_local_nucleotide_diversity_v3",
  "apc_mappability",
  "apc_mutation_density",
  "apc_protein_function_v3",
  "apc_transcription_factor",
  "linsight",
  "fathmm_xf",
  "cadd_rawscore",
  "cadd_phred"
)

# Broad categories follow Supplementary Table 1 exactly.
classify_supp_category <- function(annotation) {
  dplyr::case_when(
    annotation == "Intercept" ~ "Intercept",
    annotation %in% c(baseline_ld_annotations, favor_annotations) ~ "Baseline",
    str_detect(
      annotation,
      "^(C04_IOSE|C05_FTSEC|E01_HGSOC|E02_high|F04_endo|G02_clear)"
    ) ~ "Epigenomics",
    annotation == "raw_score" ~ "Epigenomics",
    str_detect(
      annotation,
      "^(CCOC|EEC|FT|HGSOC|IOSE|LGSOC|MOC):E[1-7]$"
    ) ~ "Epigenomics",
    annotation %in% c("MAGMA", "MultiXcan", "chromMAGMA") ~
      "Germline susceptible genes",
    str_detect(annotation, "^H3K27Ac_") ~ "Transcription factor",
    TRUE ~ "Other / review"
  )
}

# Detailed annotation-set labels preserve the organization of Supplementary
# Table 1. Chromatin states are labeled only by ovarian subtype; E1-E7 remain
# identifiers in the original Annotation column and are not translated here.
classify_annotation_set <- function(annotation) {
  dplyr::case_when(
    annotation == "Intercept" ~ "Intercept",
    annotation %in% baseline_ld_annotations ~ "BaselineLD v2.2",
    annotation %in% favor_annotations ~ "FAVOR",

    str_detect(annotation, "^C05_FTSEC") ~ "H3K27Ac FTSEC (precursor)",
    str_detect(annotation, "^C04_IOSE") ~ "H3K27Ac IOSE (precursor)",
    str_detect(annotation, "^E01_HGSOC") ~
      "H3K27Ac HGSOC (cancer cell line)",
    str_detect(annotation, "^G02_clear") ~ "H3K27Ac CCOC (tumor)",
    str_detect(annotation, "^F04_endo") ~ "H3K27Ac ENOC (tumor)",
    str_detect(annotation, "^E02_high") ~ "H3K27Ac HGSOC (tumor)",
    annotation == "raw_score" ~ "ATAC-seq FTSEC (precursor cell line)",

    str_detect(annotation, "^CCOC:E[1-7]$") ~ "Chromatin state CCOC",
    str_detect(annotation, "^EEC:E[1-7]$") ~ "Chromatin state EEC",
    str_detect(annotation, "^FT:E[1-7]$") ~ "Chromatin state FTSEC",
    str_detect(annotation, "^HGSOC:E[1-7]$") ~ "Chromatin state HGSOC",
    str_detect(annotation, "^IOSE:E[1-7]$") ~ "Chromatin state IOSE",
    str_detect(annotation, "^LGSOC:E[1-7]$") ~ "Chromatin state LGSOC",
    str_detect(annotation, "^MOC:E[1-7]$") ~ "Chromatin state MOC",

    annotation == "MAGMA" ~ "MAGMA gene list",
    annotation == "chromMAGMA" ~ "chromMAGMA gene list",
    annotation == "MultiXcan" ~ "MultiXcan gene list",

    # These are active TFBS annotations.
    str_detect(annotation, "^H3K27Ac_") ~ str_replace(
      annotation,
      "^H3K27Ac_",
      "H3K27Ac "
    ),
    TRUE ~ "Other / review"
  )
}

extract_chromatin_subtype <- function(annotation) {
  dplyr::if_else(
    str_detect(annotation, "^(CCOC|EEC|FT|HGSOC|IOSE|LGSOC|MOC):E[1-7]$"),
    str_extract(annotation, "^[^:]+"),
    NA_character_
  )
}

classify_annotation <- function(annotation) {
  dplyr::case_when(
    annotation == "Intercept" ~ "Intercept",
    annotation %in% baseline_ld_annotations ~ "BaselineLD v2.2",
    annotation %in% favor_annotations ~ "FAVOR",
    str_detect(
      annotation,
      "^(C04_IOSE|C05_FTSEC|E01_HGSOC|E02_high|F04_endo|G02_clear)"
    ) ~ "H3K27Ac epigenomics",
    annotation == "raw_score" ~ "ATAC-seq FTSEC",
    str_detect(
      annotation,
      "^(CCOC|EEC|FT|HGSOC|IOSE|LGSOC|MOC):E[1-7]$"
    ) ~ "Chromatin-state annotations",
    annotation %in% c("MAGMA", "MultiXcan", "chromMAGMA") ~
      "Germline susceptible genes",
    str_detect(annotation, "^H3K27Ac_") ~ "Transcription factor",
    TRUE ~ "Other / review"
  )
}

chromatin_state_names <- c(
  "E1" = "Weak promoter",
  "E2" = "Active promoter",
  "E3" = "Active region",
  "E4" = "Active enhancer",
  "E5" = "Weak enhancer",
  "E6" = "Insulator",
  "E7" = "Transcribed"
)

make_pretty_label <- function(annotation) {
  labels <- annotation |>
    str_replace("\\.flanking\\.500$", " (flanking 500 bp)") |>
    str_replace("\\.extend\\.500$", " (extended 500 bp)") |>
    str_replace_all("_", " ")

  chromatin_pattern <- "^(CCOC|EEC|FT|HGSOC|IOSE|LGSOC|MOC):E[1-7]$"
  is_chromatin_state <- str_detect(annotation, chromatin_pattern)

  if (any(is_chromatin_state)) {
    subtype <- str_extract(annotation[is_chromatin_state], "^[^:]+")
    subtype <- dplyr::recode(subtype, "FT" = "FTSEC")
    state_code <- str_extract(annotation[is_chromatin_state], "E[1-7]$")
    state_name <- unname(chromatin_state_names[state_code])

    # Retain the ovarian subtype while replacing the internal E1-E7 code with
    # a reader-facing biological state label. The original Annotation column
    # remains unchanged in all output tables.
    labels[is_chromatin_state] <- paste0(
      subtype,
      " - ",
      state_name
    )
  }

  labels
}

is_ovarian_source <- function(source) {
  source %in% c(
    "H3K27Ac epigenomics",
    "ATAC-seq FTSEC",
    "Chromatin-state annotations",
    "Germline susceptible genes",
    "Transcription factor"
  )
}

# Locate and validate input files ----------------------------------------------
input_files <- Sys.glob(
  file.path(
    base_dir,
    "NMOC_GROUP_*_*_sbrc_ALL",
    "*.AnnoPerSnpHsqEnrichment"
  )
)

if (length(input_files) == 0) {
  stop("No .AnnoPerSnpHsqEnrichment files found under: ", base_dir)
}

file_index <- tibble(
  File = normalizePath(input_files, mustWork = TRUE),
  Run = basename(dirname(input_files))
) |>
  extract(
    Run,
    into = c("Fold", "Panel"),
    regex = "^NMOC_GROUP_([0-9]+)_(HapMap3|Imputed)_sbrc_ALL$",
    remove = FALSE,
    convert = TRUE
  )

if (any(is.na(file_index$Fold)) || any(is.na(file_index$Panel))) {
  bad <- file_index |> filter(is.na(Fold) | is.na(Panel))
  stop(
    "Could not parse fold/panel from these directories:\n",
    paste(bad$Run, collapse = "\n")
  )
}

expected <- tidyr::crossing(Fold = 0:4, Panel = c("HapMap3", "Imputed"))
found <- file_index |> distinct(Fold, Panel)
missing_runs <- anti_join(expected, found, by = c("Fold", "Panel"))
if (nrow(missing_runs) > 0) {
  stop(
    "Missing expected runs:\n",
    paste0("Fold ", missing_runs$Fold, ", ", missing_runs$Panel,
           collapse = "\n")
  )
}

if (nrow(file_index) != 10) {
  warning(
    "Expected 10 files (5 folds x 2 panels), but found ", nrow(file_index),
    ". The script will continue after checking duplicates."
  )
}

if (any(duplicated(file_index[c("Fold", "Panel")]))) {
  stop("More than one enrichment file was found for at least one fold/panel.")
}

read_one_enrichment <- function(path, fold, panel, run) {
  x <- read_tsv(path, show_col_types = FALSE, progress = FALSE, trim_ws = TRUE)

  required_columns <- c("Annotation", "Enrich", "SD")
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0) {
    stop(
      "File is missing columns [",
      paste(missing_columns, collapse = ", "),
      "]: ", path
    )
  }

  x |>
    transmute(
      Fold = as.integer(fold),
      Panel = panel,
      Run = run,
      File = path,
      Annotation = as.character(Annotation),
      Enrich = as.numeric(Enrich),
      Posterior_SD = as.numeric(SD)
    )
}

raw_list <- lapply(seq_len(nrow(file_index)), function(i) {
  read_one_enrichment(
    file_index$File[[i]],
    file_index$Fold[[i]],
    file_index$Panel[[i]],
    file_index$Run[[i]]
  )
})
raw <- bind_rows(raw_list)

if (anyDuplicated(raw[c("Fold", "Panel", "Annotation")]) > 0) {
  stop("Duplicate annotation rows were found within at least one fold/panel.")
}

if (any(!is.finite(raw$Enrich)) || any(!is.finite(raw$Posterior_SD))) {
  stop("Non-finite Enrich or SD values were found.")
}
if (any(raw$Posterior_SD < 0)) {
  stop("Negative posterior SD values were found.")
}

annotation_coverage <- raw |>
  distinct(Fold, Panel, Annotation) |>
  count(Annotation, name = "N_runs")
if (any(annotation_coverage$N_runs != 10)) {
  warning(
    "Some annotations are not present in all 10 files. See: ",
    file.path(out_dir, "annotation_coverage.tsv")
  )
}
write_tsv(annotation_coverage, file.path(out_dir, "annotation_coverage.tsv"))

# Annotation metadata ----------------------------------------------------------
# Optional manual metadata file:
#   /path/to/SBayesRC/data/annotation_metadata.tsv
# with columns:
#   Annotation    Source    Ovarian_relevant
# Ovarian_relevant should be TRUE/FALSE.
metadata_path <- file.path(base_dir, "annotation_metadata.tsv")

auto_metadata <- tibble(Annotation = sort(unique(raw$Annotation))) |>
  mutate(
    Source = classify_annotation(Annotation),
    Ovarian_relevant = is_ovarian_source(Source),
    Plot_label = make_pretty_label(Annotation)
  )

write_tsv(
  auto_metadata,
  file.path(out_dir, "annotation_metadata_auto_review.tsv")
)

if (file.exists(metadata_path)) {
  manual_metadata <- read_tsv(
    metadata_path,
    show_col_types = FALSE,
    progress = FALSE
  )
  required_meta <- c("Annotation", "Source", "Ovarian_relevant")
  missing_meta <- setdiff(required_meta, names(manual_metadata))
  if (length(missing_meta) > 0) {
    stop(
      "annotation_metadata.tsv is missing columns: ",
      paste(missing_meta, collapse = ", ")
    )
  }

  metadata <- auto_metadata |>
    select(Annotation, Auto_Source = Source,
           Auto_Ovarian_relevant = Ovarian_relevant, Plot_label) |>
    left_join(
      manual_metadata |>
        select(Annotation, Source, Ovarian_relevant),
      by = "Annotation"
    ) |>
    mutate(
      Source = coalesce(Source, Auto_Source),
      Ovarian_relevant = coalesce(
        as.logical(Ovarian_relevant),
        Auto_Ovarian_relevant
      )
    ) |>
    select(Annotation, Source, Ovarian_relevant, Plot_label)
} else {
  metadata <- auto_metadata
}

# Add Supplementary Table 1 grouping independently of the Figure 2 source
# grouping. This keeps Figure 1B faithful to the manuscript table while
# retaining the biologically focused Figure 2 facets.
metadata <- metadata |>
  mutate(
    Category = classify_supp_category(Annotation),
    Annotation_set = classify_annotation_set(Annotation),
    Chromatin_subtype = extract_chromatin_subtype(Annotation)
  )

write_tsv(
  metadata,
  file.path(out_dir, "annotation_metadata_supp_table1.tsv")
)

unclassified <- metadata |>
  filter(Source == "Other / review" | Category == "Other / review")
if (nrow(unclassified) > 0) {
  write_tsv(
    unclassified,
    file.path(out_dir, "annotations_needing_source_review.tsv")
  )
  message(
    nrow(unclassified),
    " annotations were assigned to 'Other / review'. Review ",
    file.path(out_dir, "annotation_metadata_auto_review.tsv"),
    " and optionally create ", metadata_path
  )
}

# Fold-level approximate posterior intervals ----------------------------------
# Each file gives an SBayesRC posterior mean and posterior SD for one CV run.
# We approximate that fold-specific posterior as Normal(Enrich, Posterior_SD^2).
# These are descriptive Bayesian/model-based intervals, not sampling CIs,
# S-LDSC standard errors, or frequentist significance tests.
raw <- raw |>
  filter(Annotation != "Intercept") |>
  left_join(metadata, by = "Annotation") |>
  mutate(
    Panel = factor(Panel, levels = c("HapMap3", "Imputed")),
    Approx_posterior_95_lower = pmax(0, Enrich - 1.96 * Posterior_SD),
    Approx_posterior_95_upper = Enrich + 1.96 * Posterior_SD,
    Approx_posterior_prob_gt1 = if_else(
      Posterior_SD > 0,
      pnorm((Enrich - 1) / Posterior_SD),
      if_else(Enrich > 1, 1, 0)
    )
  )

write_tsv(raw, file.path(out_dir, "enrichment_fold_level.tsv"))

# Aggregate five folds ---------------------------------------------------------
# The primary display summaries are the five-fold mean plus the observed
# minimum and maximum fold-specific point estimates. These ranges directly
# describe fold-to-fold stability and are used in Figures 1A and 2.
#
# Posterior-SD-based pooled model columns are still calculated below for
# transparency and backward compatibility with earlier output tables. They are
# not used as error bars in the main figures. The pooled model SD follows the
# variance of an equal-weight mixture of the five approximate fold-specific
# normal posterior distributions:
# Var = mean(SD_i^2 + mean_i^2) - pooled_mean^2.
summary_by_panel <- raw |>
  group_by(Annotation, Plot_label, Source, Category, Annotation_set, Chromatin_subtype, Ovarian_relevant, Panel) |>
  summarise(
    N_folds = n(),
    Mean_enrichment = mean(Enrich),
    Median_enrichment = median(Enrich),
    Min_enrichment = min(Enrich),
    Max_enrichment = max(Enrich),
    Between_fold_SD = if (n() > 1) sd(Enrich) else NA_real_,
    RMS_posterior_SD = sqrt(mean(Posterior_SD^2)),
    Pooled_model_variance = max(
      0,
      mean(Posterior_SD^2 + Enrich^2) - mean(Enrich)^2
    ),
    Mean_posterior_prob_gt1 = mean(Approx_posterior_prob_gt1),
    Folds_point_estimate_gt1 = sum(Enrich > 1),
    Folds_approx_posterior_95_excludes_1 = sum(Approx_posterior_95_lower > 1),
    .groups = "drop"
  ) |>
  mutate(
    Pooled_model_SD = sqrt(Pooled_model_variance),
    Pooled_interval_lower = pmax(
      0,
      Mean_enrichment - 1.96 * Pooled_model_SD
    ),
    Pooled_interval_upper = Mean_enrichment + 1.96 * Pooled_model_SD,
    # These are the endpoints used for the annotation-specific error bars in
    # Figures 1A and 2. They are actual fold-level point estimates, not
    # posterior or confidence limits.
    CV_range_lower = Min_enrichment,
    CV_range_upper = Max_enrichment,
    Panel = factor(Panel, levels = c("HapMap3", "Imputed"))
  )

write_tsv(
  summary_by_panel,
  file.path(out_dir, "enrichment_summary_by_panel.tsv")
)

# Paired panel comparison, including within-fold Imputed-HapMap3 differences.
paired_folds <- raw |>
  select(Annotation, Plot_label, Source, Category, Annotation_set, Chromatin_subtype, Ovarian_relevant, Fold, Panel, Enrich) |>
  pivot_wider(names_from = Panel, values_from = Enrich) |>
  filter(!is.na(HapMap3), !is.na(Imputed)) |>
  mutate(
    Difference_Imputed_minus_HapMap3 = Imputed - HapMap3,
    Ratio_Imputed_over_HapMap3 = if_else(
      HapMap3 != 0,
      Imputed / HapMap3,
      NA_real_
    )
  )

paired_difference_summary <- paired_folds |>
  group_by(Annotation, Plot_label, Source, Category, Annotation_set, Chromatin_subtype, Ovarian_relevant) |>
  summarise(
    N_paired_folds = n(),
    Mean_difference = mean(Difference_Imputed_minus_HapMap3),
    SD_difference = if (n() > 1) sd(Difference_Imputed_minus_HapMap3) else NA_real_,
    Min_difference = min(Difference_Imputed_minus_HapMap3),
    Max_difference = max(Difference_Imputed_minus_HapMap3),
    Folds_Imputed_greater = sum(Difference_Imputed_minus_HapMap3 > 0),
    Mean_ratio = mean(Ratio_Imputed_over_HapMap3, na.rm = TRUE),
    .groups = "drop"
  )

panel_means_wide <- summary_by_panel |>
  select(
    Annotation, Plot_label, Source, Category, Annotation_set,
    Chromatin_subtype, Ovarian_relevant, Panel,
    Mean_enrichment, Min_enrichment, Max_enrichment,
    CV_range_lower, CV_range_upper,
    Pooled_interval_lower, Pooled_interval_upper
  ) |>
  pivot_wider(
    names_from = Panel,
    values_from = c(
      Mean_enrichment,
      Min_enrichment,
      Max_enrichment,
      CV_range_lower,
      CV_range_upper,
      Pooled_interval_lower,
      Pooled_interval_upper
    )
  )

panel_comparison <- panel_means_wide |>
  left_join(
    paired_difference_summary,
    by = c(
      "Annotation", "Plot_label", "Source", "Category", "Annotation_set",
      "Chromatin_subtype", "Ovarian_relevant"
    )
  ) |>
  arrange(desc(Mean_enrichment_Imputed))

write_tsv(paired_folds, file.path(out_dir, "paired_fold_differences.tsv"))
write_tsv(
  panel_comparison,
  file.path(out_dir, "enrichment_panel_comparison.tsv")
)

# Plot theme ------------------------------------------------------------------
# Figures 1 and 2 use Arial for all text, with every text element forced to black.
plot_theme <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(family = "Arial", color = "black"),
    plot.subtitle = element_text(family = "Arial", color = "black"),
    plot.caption = element_text(family = "Arial", color = "black"),
    axis.title = element_text(family = "Arial", color = "black"),
    axis.text = element_text(family = "Arial", color = "black"),
    legend.title = element_text(family = "Arial", color = "black"),
    legend.text = element_text(family = "Arial", color = "black"),
    strip.text = element_text(family = "Arial", color = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_blank()
  )

dodge <- position_dodge(width = 0.55)

# Figure 1A: top annotations across all annotation sources ---------------------
top20_candidates <- summary_by_panel |>
  filter(Panel == "Imputed")
if (length(exclude_sources_from_top20) > 0) {
  top20_candidates <- top20_candidates |>
    filter(!Source %in% exclude_sources_from_top20)
}

top20_names <- top20_candidates |>
  arrange(desc(Mean_enrichment)) |>
  slice_head(n = top_n_all) |>
  pull(Annotation)

top20_order <- summary_by_panel |>
  filter(Panel == "Imputed", Annotation %in% top20_names) |>
  arrange(Mean_enrichment) |>
  pull(Plot_label)

plot_top20_data <- summary_by_panel |>
  filter(Annotation %in% top20_names) |>
  mutate(Plot_label = factor(Plot_label, levels = top20_order))

p_top20 <- ggplot(
  plot_top20_data,
  aes(x = Plot_label, y = Mean_enrichment, color = Panel)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(
    aes(
      ymin = CV_range_lower,
      ymax = CV_range_upper
    ),
    width = 0.18,
    position = dodge,
    linewidth = 0.45
  ) +
  geom_point(position = dodge, size = 2.1) +
  coord_flip() +
  scale_color_manual(values = panel_colors, labels = panel_labels) +
  labs(
    title = "A  Top annotations across all sources",
    x = NULL,
    y = "Per-SNP heritability enrichment"
  ) +
  plot_theme

# Figure 1B-C: Supplementary Table 1 category summaries -----------------------
# The germline susceptibility gene-list annotations have much larger enrichment
# values than the other categories. They are therefore displayed in a dedicated
# panel (B), while Baseline, Epigenomics, and Transcription factor annotations
# are summarized together in panel C. This is a presentation choice based on the
# pre-specified annotation class, not on statistical significance.
category_levels <- c(
  "Baseline",
  "Epigenomics",
  "Germline susceptible genes",
  "Transcription factor"
)

unexpected_categories <- summary_by_panel |>
  filter(!Category %in% c(category_levels, "Intercept")) |>
  distinct(Category)
if (nrow(unexpected_categories) > 0) {
  warning(
    "Annotations outside Supplementary Table 1 categories were found: ",
    paste(unexpected_categories$Category, collapse = ", ")
  )
}

plot_distribution_data <- summary_by_panel |>
  filter(Category %in% category_levels) |>
  mutate(Category = factor(Category, levels = category_levels))

# Panel B: summarize the three germline susceptibility gene-list annotations
# using the same boxplot-plus-points design as panel C. The individual gene
# lists are already displayed in panel A; panel B emphasizes their category-
# level distribution on a dedicated scale.
plot_germline_data <- summary_by_panel |>
  filter(Category == "Germline susceptible genes") |>
  mutate(
    Category = factor(
      Category,
      levels = "Germline susceptible genes"
    )
  )

p_germline <- ggplot(
  plot_germline_data,
  aes(x = Category, y = Mean_enrichment, color = Panel)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_boxplot(
    aes(group = interaction(Category, Panel)),
    position = position_dodge(width = 0.75),
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    size = 1.25,
    alpha = 0.75
  ) +
  coord_flip() +
  scale_color_manual(values = panel_colors, labels = panel_labels) +
  labs(
    title = "B  Germline susceptibility gene-list enrichment",
    x = NULL,
    y = "Fold-mean per-SNP heritability enrichment"
  ) +
  plot_theme

# Panel C: order the remaining three categories by their median enrichment in
# the Imputed panel, from largest at the top to smallest at the bottom.
remaining_categories <- c(
  "Baseline",
  "Epigenomics",
  "Transcription factor"
)

remaining_category_order <- summary_by_panel |>
  filter(
    Panel == "Imputed",
    Category %in% remaining_categories
  ) |>
  group_by(Category) |>
  summarise(
    Imputed_median_enrichment = median(Mean_enrichment),
    .groups = "drop"
  ) |>
  arrange(Imputed_median_enrichment) |>
  pull(Category)

plot_remaining_data <- summary_by_panel |>
  filter(Category %in% remaining_categories) |>
  mutate(Category = factor(Category, levels = remaining_category_order))

p_remaining <- ggplot(
  plot_remaining_data,
  aes(x = Category, y = Mean_enrichment, color = Panel)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_boxplot(
    aes(group = interaction(Category, Panel)),
    position = position_dodge(width = 0.75),
    width = 0.58,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.13,
      dodge.width = 0.75
    ),
    size = 0.85,
    alpha = 0.55
  ) +
  coord_flip() +
  scale_color_manual(values = panel_colors, labels = panel_labels) +
  labs(
    title = "C  Enrichment across remaining annotation categories",
    x = NULL,
    y = "Fold-mean per-SNP heritability enrichment"
  ) +
  plot_theme

figure1_caption <- paste(
  "Annotations are grouped according to the four categories defined in",
  "Supplementary Table 1. Panel A shows mean SBayesRC-derived per-SNP",
  "heritability enrichment across five cross-validation models; error bars",
  "span the minimum and maximum fold-specific point estimates. Panels B and C summarize",
  "category-level distributions: boxes show the median and interquartile",
  "range, and points represent individual annotations. Germline susceptibility",
  "gene lists are displayed separately in panel B because their enrichment",
  "values are substantially larger. The dashed line denotes enrichment = 1."
)

right_column <- p_germline / p_remaining +
  plot_layout(heights = c(0.62, 1.38), guides = "collect")

figure1 <- p_top20 | right_column
figure1 <- figure1 +
  plot_layout(widths = c(1.25, 1), guides = "collect") +
  plot_annotation(caption = figure1_caption)
figure1 <- figure1 &
  theme(
    text = element_text(family = "Arial", color = "black"),
    legend.position = "bottom",
    plot.caption = element_text(
      family = "Arial",
      color = "black",
      size = 9,
      hjust = 0,
      margin = margin(t = 8)
    )
  )

ggsave(
  file.path(out_dir, "Figure1_all_annotations.pdf"),
  figure1,
  width = 13.5,
  height = 9.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(out_dir, "Figure1_all_annotations.png"),
  figure1,
  width = 13.5,
  height = 9.0,
  units = "in",
  dpi = 300
)

# Figure 2: focused annotation-specific enrichment -----------------------------
# Up to eight annotations are selected from Baseline, Epigenomics, and
# Transcription factor categories according to mean enrichment in the Imputed
# panel. Germline susceptibility gene lists are omitted here because they are
# already shown individually in Figure 1A and summarized in Figure 1B.
figure2_category_levels <- c(
  "Baseline",
  "Epigenomics",
  "Transcription factor"
)

figure2_top_names <- summary_by_panel |>
  filter(
    Panel == "Imputed",
    Category %in% figure2_category_levels
  ) |>
  mutate(Category = factor(Category, levels = figure2_category_levels)) |>
  group_by(Category) |>
  arrange(desc(Mean_enrichment), .by_group = TRUE) |>
  slice_head(n = top_n_per_category) |>
  ungroup() |>
  pull(Annotation)

figure2_order_table <- summary_by_panel |>
  filter(
    Panel == "Imputed",
    Annotation %in% figure2_top_names
  ) |>
  mutate(Category = factor(Category, levels = figure2_category_levels)) |>
  arrange(Category, Mean_enrichment) |>
  mutate(Facet_label = paste(Plot_label, Category, sep = "___"))

figure2_levels <- figure2_order_table$Facet_label

plot_figure2_data <- summary_by_panel |>
  filter(Annotation %in% figure2_top_names) |>
  mutate(
    Category = factor(Category, levels = figure2_category_levels),
    Facet_label = paste(Plot_label, Category, sep = "___"),
    Facet_label = factor(Facet_label, levels = figure2_levels)
  )

p_category_annotations <- ggplot(
  plot_figure2_data,
  aes(x = Facet_label, y = Mean_enrichment, color = Panel)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(
    aes(
      ymin = CV_range_lower,
      ymax = CV_range_upper
    ),
    width = 0.18,
    position = dodge,
    linewidth = 0.45
  ) +
  geom_point(position = dodge, size = 2.0) +
  coord_flip() +
  facet_wrap(~Category, scales = "free", ncol = 3) +
  scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_color_manual(values = panel_colors, labels = panel_labels) +
  labs(
    title = "Focused annotation-specific per-SNP heritability enrichment",
    subtitle = paste0(
      "Top ", top_n_per_category,
      " annotations within Baseline, Epigenomics, and Transcription factor categories"
    ),
    caption = paste(
      "This focused view excludes germline susceptibility gene lists, which are",
      "shown in Figure 1. Within each category, annotations are ordered from",
      "highest to lowest mean enrichment in the Imputed panel. Chromatin-state",
      "labels show ovarian subtype and biological state; original E1-E7 codes are",
      "retained in the accompanying data table. Error bars span the minimum and",
      "maximum fold-specific point estimates. Numerical axes vary by facet to improve visibility."
    ),
    x = NULL,
    y = "Per-SNP heritability enrichment"
  ) +
  plot_theme +
  theme(
    plot.caption = element_text(
      family = "Arial",
      color = "black",
      size = 9,
      hjust = 0,
      margin = margin(t = 8)
    )
  )

ggsave(
  file.path(out_dir, "Figure2_focused_regulatory_annotations.pdf"),
  p_category_annotations,
  width = 15,
  height = 7.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(out_dir, "Figure2_focused_regulatory_annotations.png"),
  p_category_annotations,
  width = 15,
  height = 7.2,
  units = "in",
  dpi = 300
)

# Supplementary figure: panel-to-panel enrichment comparison ------------------
scatter_data <- panel_comparison |>
  filter(
    is.finite(Mean_enrichment_HapMap3),
    is.finite(Mean_enrichment_Imputed)
  )

scatter_limit <- max(
  c(
    scatter_data$Mean_enrichment_HapMap3,
    scatter_data$Mean_enrichment_Imputed
  ),
  na.rm = TRUE
)

p_scatter <- ggplot(
  scatter_data,
  aes(
    x = Mean_enrichment_HapMap3,
    y = Mean_enrichment_Imputed,
    color = Category
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 1.8, alpha = 0.78) +
  coord_equal(xlim = c(0, scatter_limit), ylim = c(0, scatter_limit)) +
  labs(
    title = "Comparison of enrichment between reference panels",
    subtitle = "Each point is one annotation averaged over five CV folds",
    x = "Per-SNP enrichment, HapMap3",
    y = "Per-SNP enrichment, Imputed",
    color = "Supplementary Table 1 category"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "FigureS1_panel_comparison_scatter.pdf"),
  p_scatter,
  width = 9.5,
  height = 7.5,
  units = "in"
)
ggsave(
  file.path(out_dir, "FigureS1_panel_comparison_scatter.png"),
  p_scatter,
  width = 9.5,
  height = 7.5,
  units = "in",
  dpi = 300
)

# Supplementary figure: fold-level stability for the overall top annotations ---
fold_top_data <- raw |>
  filter(Annotation %in% top20_names) |>
  mutate(Plot_label = factor(Plot_label, levels = top20_order))

p_fold_stability <- ggplot(
  fold_top_data,
  aes(
    x = factor(Fold),
    y = Enrich,
    group = Panel,
    color = Panel
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey65") +
  geom_line(linewidth = 0.45, alpha = 0.75) +
  geom_point(size = 1.35) +
  facet_wrap(~Plot_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = panel_colors, labels = panel_labels) +
  labs(
    title = "Fold-level stability of top annotation enrichment",
    x = "Cross-validation fold",
    y = "Per-SNP heritability enrichment"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7.5),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, "FigureS2_fold_stability.pdf"),
  p_fold_stability,
  width = 12,
  height = 10,
  units = "in"
)
ggsave(
  file.path(out_dir, "FigureS2_fold_stability.png"),
  p_fold_stability,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300
)

# Tables used directly by figures ---------------------------------------------
write_tsv(
  plot_top20_data |> arrange(Panel, desc(Mean_enrichment)),
  file.path(out_dir, "Figure1A_top20_data.tsv")
)
write_tsv(
  plot_germline_data |> arrange(Panel, desc(Mean_enrichment)),
  file.path(out_dir, "Figure1B_germline_gene_list_data.tsv")
)
write_tsv(
  plot_remaining_data |> arrange(Category, Panel, desc(Mean_enrichment)),
  file.path(out_dir, "Figure1C_remaining_category_distribution_data.tsv")
)
write_tsv(
  plot_distribution_data |> arrange(Category, Panel, desc(Mean_enrichment)),
  file.path(out_dir, "Figure1BC_all_category_data.tsv")
)
write_tsv(
  plot_figure2_data |> arrange(Category, Panel, desc(Mean_enrichment)),
  file.path(out_dir, "Figure2_focused_regulatory_annotations_data.tsv")
)

write_tsv(
  plot_distribution_data |>
    distinct(Annotation, Category) |>
    count(Category, name = "N_annotations"),
  file.path(out_dir, "supp_table1_category_counts.tsv")
)

# Compact text summary ---------------------------------------------------------
comparison_summary <- paired_folds |>
  summarise(
    N_annotations = n_distinct(Annotation),
    N_fold_annotation_pairs = n(),
    Percent_pairs_Imputed_greater = 100 * mean(
      Difference_Imputed_minus_HapMap3 > 0
    )
  )

sink(file.path(out_dir, "analysis_notes.txt"))
cat("SBayesRC enrichment summary\n")
cat("===========================\n\n")
cat("Input directory: ", base_dir, "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Files analyzed: ", nrow(file_index), "\n", sep = "")
cat("Annotations analyzed (excluding Intercept): ",
    n_distinct(raw$Annotation), "\n", sep = "")
cat("Fold-annotation pairs with Imputed > HapMap3: ",
    round(comparison_summary$Percent_pairs_Imputed_greater, 1), "%\n", sep = "")
cat("\nImportant interpretation note:\n")
cat(paste(
  "In Figures 1A and 2, each point is the mean enrichment across the five CV",
  "models and each error bar spans the minimum to maximum of the five",
  "fold-specific enrichment point estimates. These ranges describe fold-level",
  "stability; they are not confidence intervals, posterior credible intervals,",
  "S-LDSC standard errors, or formal significance tests. Posterior-SD-based",
  "approximate intervals are retained only in the output tables for reference.\n"
))
sink()

message("Completed. Results written to: ", out_dir)
message("Main figures:")
message("  ", file.path(out_dir, "Figure1_all_annotations.pdf"))
message("  ", file.path(out_dir, "Figure2_focused_regulatory_annotations.pdf"))
message("Supplementary figures:")
message("  ", file.path(out_dir, "FigureS1_panel_comparison_scatter.pdf"))
message("  ", file.path(out_dir, "FigureS2_fold_stability.pdf"))
