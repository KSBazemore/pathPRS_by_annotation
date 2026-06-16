#!/bin/bash

set -euo pipefail

##############################
#Set up
inc_set="${1:-}"

PRSiceR="path/to/PRSice.R"
PRSice_linux="path/to/PRSice_linux"

home_dir="path/to/home_directory"

UKBdat="path/to/UKB_data_prefix"
base_file="path/to/IGAP_summary_statistics_QCd.txt"
qcd_snplist="path/to/qcd_snplist.txt"
snp_sets="${home_dir}/cluster_snp_sets.txt"
 
LDref="path/to/LD_reference_prefix"


###############################
#Set inputs for training

    case ${inc_set} in
    	A) 
    	  cov_data="${home_dir}/SetA_PRSice_ready.txt"
    		covs="age,sex,e4_count,@f.22009.0.[1-5]"
    		outcomes="AD_A_Tr_1,AD_A_Tr_2,AD_A_Tr_3,AD_A_Tr_4,AD_A_Tr_5,AD_A_Tr_6,AD_A_Tr_7,AD_A_Tr_8,AD_A_Tr_9,AD_A_Tr_10"
    		;;
    	B) 
    	  cov_data="${home_dir}/SetB_PRSice_ready.txt"
    		covs="age,sex,e4_count,@f.22009.0.[1-5]"
    		outcomes="AD_B_Tr_1,AD_B_Tr_2,AD_B_Tr_3,AD_B_Tr_4,AD_B_Tr_5,AD_B_Tr_6,AD_B_Tr_7,AD_B_Tr_8,AD_B_Tr_9,AD_B_Tr_10"
    		;;
    	C) 
    	  cov_data="${home_dir}/SetC_PRSice_ready.txt"
    		covs="sex,e4_count,f.34.0.0,f.22000.0.0,f_age,m_age,@f.22009.0.[1-5]"
    		outcomes="AD_C_Tr_1,AD_C_Tr_2,AD_C_Tr_3,AD_C_Tr_4,AD_C_Tr_5,AD_C_Tr_6,AD_C_Tr_7,AD_C_Tr_8,AD_C_Tr_9,AD_C_Tr_10"
    		;;
  		*)
    		echo "Error: inc_set must be A, B, or C"
    		exit 1
    		;;
    esac


###############################
#Set output directory for training results

out_dir="${home_dir}/set_${inc_set}_multi"
mkdir -p "${out_dir}"


##############################
#Training

Rscript "${PRSiceR}" \
	--prsice "${PRSice_linux}" \
	--base "${base_file}" \
	--beta \
	--extract "${qcd_snplist}" \
	--snp-set "${snp_sets}" \
	--target "${UKBdat}" \
	--type bed \
	--binary-target T,T,T,T,T,T,T,T,T,T \
	--proxy 0.8 \
	--ld "${LDref}" \
	--pheno "${cov_data}" \
	--pheno-col "${outcomes}" \
	--ignore-fid \
	--cov "${cov_data}" \
	--cov-col "${covs}" \
	--cov-factor sex,e4_count \
	--bar-levels 5e-08,5e-06,5e-04,0.05,0.1,0.5 \
	--print-snp \
	--fastscore \
	--score std \
	--out "${out_dir}/training"
  	



