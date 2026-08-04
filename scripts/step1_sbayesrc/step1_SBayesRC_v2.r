# PUBLIC RELEASE DOCUMENTATION
# Study role: fit one fold-specific SBayesRC model and create the corresponding
# PLINK2 weight file with chromosome, position, reference, and alternate alleles.
# Positional arguments are: GROUP_NUM, LD_FOLDER, WITH_ANNOTATION, and (when
# WITH_ANNOTATION is TRUE) ANNOTATION_FILE. GROUP_NUM identifies the held-out
# group. LD_FOLDER selects the HapMap3 or Imputed reference directory.
# The random seed, working directory, file naming, SBayesRC function calls, and
# output transformations below are unchanged from the script used in the study.
# The absolute paths are therefore historical, site-specific requirements.

#
# Part 1: SBayesRC Main Functions
#

# Library loading
library(SBayesRC)


# For reproducibility
set.seed(16888)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript step1_SBayesRC_v2.r <GROUP_NUM> <LD_FOLDER> <WITH_ANNOTATION> [ANNOTATION_FILE]")
}


# Parse arguments
GROUP_NUM <- as.integer(args[1])
LD_FOLDER <- args[2]
WITH_ANNOTATION <- as.logical(args[3])
ANNOT_FILE <- ifelse(WITH_ANNOTATION, args[4], NA)


# Working directory
WK_DIR <- "/path/to/SBayesRC/data"
setwd(WK_DIR)


# Define GWAS input file
ma_file <- sprintf("all_non_mucinous_group%d.rsID.cojoformat.txt", GROUP_NUM)


# Define output prefix
out_prefix <- ifelse(LD_FOLDER == "ukbEUR_HM3",
	sprintf("NMOC_GROUP_%d_HapMap3", GROUP_NUM),
	sprintf("NMOC_GROUP_%d_Imputed", GROUP_NUM))


# Run tidy
OUTPUT_TIDY <- paste0(out_prefix, "_tidy.ma")
tidy(mafile = ma_file, LDdir = LD_FOLDER, output = OUTPUT_TIDY, log2file = TRUE)


# Run impute
OUTPUT_IMPUTE <- paste0(out_prefix, "_imp.ma")
impute(mafile = OUTPUT_TIDY, LDdir = LD_FOLDER, output = OUTPUT_IMPUTE, log2file = TRUE)


# Run SBayesRC
if (WITH_ANNOTATION) {
	ANNOT_LABEL <- sub("^annot_", "", ANNOT_FILE)
	ANNOT_LABEL <- sub("\\.txt$", "", ANNOT_LABEL)
	OUTPUT_SBRC <- paste0(out_prefix, "_sbrc_", ANNOT_LABEL)
	sbayesrc(mafile = OUTPUT_IMPUTE, LDdir = LD_FOLDER, outPrefix = OUTPUT_SBRC, annot = ANNOT_FILE, log2file = TRUE)
} else {
	OUTPUT_SBRC <- paste0(out_prefix, "_sbrc_noAnnot")
	sbayesrc(mafile = OUTPUT_IMPUTE, LDdir = LD_FOLDER, outPrefix = OUTPUT_SBRC, log2file = TRUE)
}

# Create output folder
system(paste("mkdir -p", OUTPUT_SBRC), wait = TRUE)


#
# Part 2: Add SNP information for PRS Calculation
#

# Library loading
library(tidyverse)


# Define file paths
INPUT_FILE <- paste0(OUTPUT_SBRC, ".txt")
SNP_INFO_FILE <- "ukbEUR_Imputed/snp.info"
OUTPUT_DIR <- file.path(WK_DIR, "SBayesRC_PGM")
OUTPUT_FILE <- file.path(OUTPUT_DIR, paste0(OUTPUT_SBRC, "_ID.txt"))


# Create output directory if it doesn't exist
system(paste("mkdir -p", OUTPUT_DIR), wait = TRUE)


# Load weight file
weight_rsid <- read.csv(INPUT_FILE, sep = "\t", header = TRUE)


# Load SNP information file from LD folder
snp_info <- read.csv(SNP_INFO_FILE, sep = "\t", header = TRUE) %>%
setNames(c("CHR", "SNP", "Index", "GenPos", "POSITION", "ALT", "REF", "AltFreq", "N", "Block")) %>%
select(c(CHR, SNP, POSITION, ALT, REF))


# Merge and reformat
new_weight_file <- weight_rsid %>%
left_join(snp_info, by = 'SNP') %>%
unite('tmp', c('POSITION', 'REF', 'ALT'), sep = "_", remove = FALSE) %>%
unite('ID', c('CHR', 'tmp'), sep = ":", remove = FALSE) %>%
select(!tmp) %>%
relocate(ID, .before = SNP) %>%
relocate(SNP, .after = ALT) %>%
relocate(REF, .before = ALT)


# Save outputs
write.table(new_weight_file, file = OUTPUT_FILE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
message("Wrote merged file: ", OUTPUT_FILE)

# Move SBayesRC outputs into their own folder
system(paste("mv", paste0(OUTPUT_SBRC, ".*"), OUTPUT_SBRC), wait = TRUE)
