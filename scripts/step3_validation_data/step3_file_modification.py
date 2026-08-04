#!/bin/env python3

# PUBLIC RELEASE DOCUMENTATION
# Study role: combine phenotype records, summed PLINK2 PRS values, and the
# pre-existing five-group assignment table; then extract held-out NMOC cases
# and controls for the group supplied on the command line.
# Usage: python step3_file_modification.py <output_folder> <group_num>
# Case/control definitions are all_non_mucinous == 1 and == 0, respectively.
# All paths, merge behavior, column names, and outputs below are unchanged from
# the study script. The individual-level inputs are controlled-access data and
# are not distributed with the public code.

import pandas as pd
import os
import sys


# Check if arguments are provided
if len(sys.argv) != 3:
	print("Error: Must provide output_folder and group number")
	print("Usage: python step3_file_modification.py <output_folder> <group_num>")
	sys.exit(1)


# Get arguments
OUTPUT_FOLDER = sys.argv[1]
GROUP_NUM = int(sys.argv[2])


# Define paths
PHENO_FILE = "/path/to/data/Genotype_EOC/joined_phenotypes.txt"
PRS_RESULT = f"/path/to/SBayesRC/result/{OUTPUT_FOLDER}/score_file/prs_result_sum.txt"
GROUPS_FILE = "/path/to/EOC_Data/groups_info.txt"
OUTPUT_DIR = f"/path/to/SBayesRC/result/{OUTPUT_FOLDER}"

# Check if files exist
for file_path in [PHENO_FILE, PRS_RESULT, GROUPS_FILE]:
	if not os.path.isfile(file_path):
		print(f"Error: File not found: {file_path}")
		sys.exit(1)

print("START")

# Read files
pheno = pd.read_csv(PHENO_FILE, sep = '\t')
prs = pd.read_csv(PRS_RESULT, sep = ' ', header = 0)
group = pd.read_csv(GROUPS_FILE, sep = '\t', names = ['GROUP', 'IID'])


# Merg phenotype info and PRS data
pheno_prs = pd.merge(pheno, prs, left_on = 'OCACID', right_on = 'IID')


# Process NMOC histotype
prs_case = pheno_prs[pheno_prs['all_non_mucinous'] == 1]
prs_control = pheno_prs[pheno_prs['all_non_mucinous'] == 0]


# Merge with group info
prs_case_group_info = pd.merge(prs_case, group, on = 'IID')
prs_control_group_info = pd.merge(prs_control, group, on = 'IID')


# Filter for specified group
prs_case_group = prs_case_group_info[prs_case_group_info['GROUP'] == GROUP_NUM]
prs_control_group = prs_control_group_info[prs_control_group_info['GROUP'] == GROUP_NUM]


# Save output files
case_file = os.path.join(OUTPUT_DIR, f"nmoc_grp{GROUP_NUM}_prs_case.txt")
control_file = os.path.join(OUTPUT_DIR, f"nmoc_grp{GROUP_NUM}_prs_control.txt")
prs_case_group.to_csv(case_file, sep='\t', index = False, na_rep = 'NA')
prs_control_group.to_csv(control_file, sep='\t', index = False, na_rep = 'NA')

print(f"Saved {case_file} and {control_file}")
print("DONE")
