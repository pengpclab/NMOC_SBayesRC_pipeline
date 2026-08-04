# PUBLIC RELEASE DOCUMENTATION
# Study role: evaluate a held-out fold using AUC with a DeLong 95% interval, a
# likelihood-ratio test comparing PRS and intercept-only logistic models, and
# the log odds ratio per standard deviation of PRS with its 95% interval.
# Positional arguments are the PRS result folder and group number. The script
# consumes the case/control files written by step3_file_modification.py and
# writes nmoc_grp<GROUP_NUM>_prs_results.csv in the same result folder.
# Statistical calculations, paths, variable coding, and output columns below
# are unchanged from the script used for the reported analysis.

# Load the packages
library(tidyverse)
library(pROC)
library(lmtest)

# Set working directory
setwd("/path/to/SBayesRC/result")


# Check if argument is provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Error: Must provide the folder path and group number (0-4)")
}

FOLDER <- args[1]
GROUP_NUM <- as.integer(args[2])


# Set up the variables
PRS_CASE <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_case.txt")
PRS_CONTROL <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_control.txt")


# Read files
pheno_case <- read.csv(PRS_CASE, sep = '\t') %>%
select('IID', 'prs_score')

pheno_control <- read.csv(PRS_CONTROL, sep = '\t') %>%
select('IID', 'prs_score')

pheno_case['PHENO'] <- 2
pheno_control['PHENO'] <- 1


# AUC calculation
pheno = rbind(pheno_control, pheno_case)

roc_pheno <- roc(
  response = pheno$PHENO,
  predictor = pheno$prs_score,
  levels = c(1, 2),   # 1 = control, 2 = case
  direction = "<",    # higher PRS is expected in cases
  quiet = TRUE
)

auc_pheno <- as.numeric(auc(roc_pheno))


# 95% confidence interval for AUC using DeLong's method
auc_ci <- ci.auc(
  roc_pheno,
  conf.level = 0.95,
  method = "delong"
)

auc_lower <- as.numeric(auc_ci[1])
auc_upper <- as.numeric(auc_ci[3])

print(paste("AUC:", auc_pheno))
print(paste(
  "95% CI for AUC:",
  auc_lower,
  "-",
  auc_upper
))


# log-likelihood chi-square calculation
norm_pheno <- pheno %>%
mutate(prs_standardized = scale(prs_score, scale = TRUE, center = TRUE)[,1]) %>% # mean = 0; sd = 1
mutate(PHENO = if_else(PHENO == 2, 1, 0))

model_with_PRS <- glm(PHENO ~ prs_standardized, data = norm_pheno, family = "binomial")
model_without_PRS <- glm(PHENO ~ 1, data = norm_pheno, family = "binomial")

LLR <- lrtest(model_without_PRS, model_with_PRS)
LLR_Chisq <- LLR$Chisq[2]
print(paste("Log-likelihood chi-square:", LLR_Chisq))


# log(Odds Ratio per SD)
logoddsratio <- model_with_PRS$coefficients[-1]
sd_prs <- sd(norm_pheno$prs_standardized)
log_oddsratio_per_sd <- logoddsratio*sd_prs
print(paste("Log Odds ratio per standard deviation:", log_oddsratio_per_sd ))


# 95% CI (Confidence Interval) (added in the version 2 script)
coef_summary <- summary(model_with_PRS)$coefficients
beta <- coef_summary["prs_standardized", "Estimate"] # same as logoddsratio
se <- coef_summary["prs_standardized", "Std. Error"]
lower <- beta - 1.96 * se
upper <- beta + 1.96 * se


# Export results to a CSV file
results <- data.frame(
	Folder = FOLDER,
	Group = GROUP_NUM,
	AUC = auc_pheno,
	AUC_LCI = auc_lower,
	AUC_UCI = auc_upper,
	Log_Likelihood_Chisq = LLR_Chisq,
	Log_Odds_Ratio_per_SD = log_oddsratio_per_sd,
	Standard_Error = se,
	Log_OR_LCI = lower,
	Log_OR_UCI = upper
)
output_file <- paste0(FOLDER, "/nmoc_grp", GROUP_NUM, "_prs_results.csv")
write.csv(results, output_file, row.names = FALSE, quote = FALSE)
print(paste("Results saved to:", output_file))
