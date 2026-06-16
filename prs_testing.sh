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


###############################
#Set input parameters for testing    

    case "${inc_set}" in
            A) 
              cov_data="${home_dir}/SetA_PRSice_ready.txt"
          		covs="age,sex,e4_count,@f.22009.0.[1-5]"
          		;;
            B) 
              cov_data="${home_dir}/SetB_PRSice_ready.txt"
          		covs="age,sex,e4_count,@f.22009.0.[1-5]"
          		;;
            C) 
              cov_data="${home_dir}/SetC_PRSice_ready.txt"
          		covs="sex,e4_count,f.34.0.0,f.22000.0.0,f_age,m_age,@f.22009.0.[1-5]"
          		;;
        		*)
          		echo "Error: inc_set must be A, B, or C"
          		exit 1
          		;;
    esac
    

    
###############################
#Set output directory for testing results

out_dir="${home_dir}/set_${inc_set}_multi"
mkdir -p "${out_dir}"


##############################
#Testing

#For testing, loop over each holdout because the snp-set will vary based on training results
for i in {1..10}; do 

  phenotype="AD_${inc_set}_Ts_${i}"
  snp_sets="${home_dir}/set${inc_set}_trained_snp_sets_${i}.txt"
  
    Rscript "${PRSiceR}" \
  	--prsice "${PRSice_linux}" \
  	--base "${base_file}" \
  	--snp-set "${snp_sets}" \
  	--beta \
  	--target "${UKBdat}" \
  	--type bed \
  	--binary-target T \
  	--no-clump \
  	--pheno "${cov_data}" \
  	--pheno-col "${phenotype}" \
  	--ignore-fid \
  	--cov "${cov_data}" \
  	--cov-col "${covs}" \
  	--cov-factor sex,e4_count \
  	--bar-levels 1 \
  	--print-snp \
  	--fastscore \
  	--score std \
  	--out "${out_dir}/testing_${i}"
  	
done
  	

