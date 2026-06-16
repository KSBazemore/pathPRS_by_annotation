#!/bin/bash

set -euo pipefail

####################################
#Setup

home_dir="/path/to/home_directory"

#Pre-QC'd IGAP summary statistics. Full summary statistics available at https://dss.niagads.org/datasets/ng00075/
input_sumstats="IGAP_summary_statistics_QCd.txt"

#Path to magma executable
magma="/path/to/magma"

#Gene boundaries file, same build as summary statistics
gene_loc="/path/to/gene_bounds"

#Created files
sans_apoe="${home_dir}/IGAP_sansAPOE.txt"
snp_locs="${home_dir}/IGAP_snplocs.txt"


####################################
#Filter summary statistics and perform annotation

#Assumed column order
#1 = chromosome
#2 = base-pair position
#3 = SNP ID

#Remove APOE region, for hg19
awk 'NR==1 || (($1 != 19 && $1 != "chr19") || $2 < 44377739 || $2 > 46125148)' "${home_dir}/${input_sumstats}" > "${sans_apoe}"

#Prepare snp location file
awk 'NR > 1 {print $3, $1, $2}' "${sans_apoe}" > "${snp_locs}"

#S-1 SNP to gene annotation
"$magma" --annotate window=35,10 --snp-loc "${snp_locs}" --gene-loc "${gene_loc}" --out "${home_dir}/S_1_annotation"

