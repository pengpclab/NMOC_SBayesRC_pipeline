# PUBLIC RELEASE DOCUMENTATION
# Study role: estimate odds ratios across control-defined PRS percentile bins.
# Cut points are the 10th, 20th, 30th, 40th, 60th, 70th, 80th, and 90th
# percentiles among controls. The 40-60% bin is the logistic-regression
# reference group. The output reports case/control counts, ORs, 95% intervals,
# and p-values for the remaining bins.
# Positional arguments are the PRS result folder and group number. All analysis
# code, paths, bin boundaries, regression calls, and outputs below are unchanged
# from the study script.

# Load packages
library(tidyverse)

  # Set working directory
setwd("/path/to/SBayesRC/result")
  
  # Check arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Error: Must provide the folder path and group number (0-4)")
}
  
FOLDER    <- args[1]
GROUP_NUM <- as.integer(args[2])

  # Read files
  PRS_CASE    <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_case.txt")
  PRS_CONTROL <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_control.txt")

  pheno_case <- read.csv(PRS_CASE, sep = '\t') %>%
    select('IID', 'prs_score') %>%
    mutate(PHENO = 2)

  pheno_control <- read.csv(PRS_CONTROL, sep = '\t') %>%
    select('IID', 'prs_score') %>%
    mutate(PHENO = 1)

  pheno <- rbind(pheno_control, pheno_case) %>%
    mutate(PHENO_binary = if_else(PHENO == 2, 1, 0))

  # Compute cutoffs from controls only
  p10 <- quantile(pheno_control$prs_score, 0.10)
  p20 <- quantile(pheno_control$prs_score, 0.20)
  p30 <- quantile(pheno_control$prs_score, 0.30)
  p40 <- quantile(pheno_control$prs_score, 0.40)
  p60 <- quantile(pheno_control$prs_score, 0.60)
  p70 <- quantile(pheno_control$prs_score, 0.70)
  p80 <- quantile(pheno_control$prs_score, 0.80)
  p90 <- quantile(pheno_control$prs_score, 0.90)

  # Assign everyone to a bin
  pheno <- pheno %>%
    mutate(prs_bin = case_when(
      prs_score <  p10                   ~ "bin1_0_10",
      prs_score >= p10 & prs_score < p20 ~ "bin2_10_20",
      prs_score >= p20 & prs_score < p30 ~ "bin3_20_30",
      prs_score >= p30 & prs_score < p40 ~ "bin4_30_40",
      prs_score >= p40 & prs_score < p60 ~ "ref_40_60",
      prs_score >= p60 & prs_score < p70 ~ "bin5_60_70",
      prs_score >= p70 & prs_score < p80 ~ "bin6_70_80",
      prs_score >= p80 & prs_score < p90 ~ "bin7_80_90",
      prs_score >= p90                   ~ "bin8_90_100"
    )) %>%
    mutate(prs_bin = factor(prs_bin, levels = c(
      "ref_40_60", "bin1_0_10", "bin2_10_20", "bin3_20_30", "bin4_30_40",
      "bin5_60_70", "bin6_70_80", "bin7_80_90", "bin8_90_100"
    )))

  # Print case/control counts per bin
  print("Case/Control counts per bin:")
  print(table(pheno$prs_bin, pheno$PHENO_binary))

  # Logistic regression — reference group is 40-60%
  model        <- glm(PHENO_binary ~ prs_bin, data = pheno, family = "binomial")
  coef_summary <- summary(model)$coefficients

  # Extract OR and CI for each non-reference bin
  bins       <- c("bin1_0_10", "bin2_10_20", "bin3_20_30", "bin4_30_40",
                  "bin5_60_70", "bin6_70_80", "bin7_80_90", "bin8_90_100")
  bin_labels <- c("0-10%", "10-20%", "20-30%", "30-40%",
                  "60-70%", "70-80%", "80-90%", "90-100%")

  results_list <- list()

  # Reference row
  results_list[["ref"]] <- data.frame(
    Folder     = FOLDER,
    Group      = GROUP_NUM,
    Bin        = "40-60% (ref)",
    N_Cases    = sum(pheno$prs_bin == "ref_40_60" & pheno$PHENO_binary == 1),
    N_Controls = sum(pheno$prs_bin == "ref_40_60" & pheno$PHENO_binary == 0),
    OR         = 1.0,
    OR_LCI     = NA,
    OR_UCI     = NA,
    Pvalue     = NA
  )

  # All other bins
  for (i in seq_along(bins)) {
    bin       <- bins[i]
    label     <- bin_labels[i]
    coef_name <- paste0("prs_bin", bin)

    beta  <- coef_summary[coef_name, "Estimate"]
    se    <- coef_summary[coef_name, "Std. Error"]
    pval  <- coef_summary[coef_name, "Pr(>|z|)"]
    OR       <- exp(beta)
    OR_lower <- exp(beta - 1.96 * se)
    OR_upper <- exp(beta + 1.96 * se)

    print(paste(label, "OR:", round(OR, 3), "95% CI:", round(OR_lower, 3), "-", round(OR_upper, 3), "P:", pval))

    results_list[[bin]] <- data.frame(
      Folder     = FOLDER,
      Group      = GROUP_NUM,
      Bin        = label,
      N_Cases    = sum(pheno$prs_bin == bin & pheno$PHENO_binary == 1),
      N_Controls = sum(pheno$prs_bin == bin & pheno$PHENO_binary == 0),
      OR         = OR,
      OR_LCI     = OR_lower,
      OR_UCI     = OR_upper,
      Pvalue     = pval
    )
  }

  results <- do.call(rbind, results_list)

  output_file <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_results_deciles.csv")
  write.csv(results, output_file, row.names = FALSE, quote = FALSE)
  print(paste("Results saved to:", output_file))
