####################################
#Setup

library(dplyr)
library(data.table)

home_dir <- "/path/to/home_directory"

source("path/to/functions.R")

####################################
#Inputs

inc_set <- "A"

trained_snps_file <- file.path(home_dir, paste0("set_", inc_set, "_multi/training.snp"))
trained_model_file <- file.path(home_dir, paste0("set_", inc_set, "_multi/training.summary"))


###################################
#Filter SNP sets to only those in best-fit PRS for each holdout

trained_snps <- fread(trained_snps_file)
trained_model <- fread(trained_model_file)

phenotypes <- unique(trained_model$Phenotype)

for (k in 1:length(phenotypes)) {
  
  pheno <- phenotypes[k]
  
  this_trained_model <- trained_model[which(trained_model$Phenotype == pheno)]
  
  out_file <- file.path(home_dir, paste0("set", inc_set, "_trained_snp_sets_", k, ".txt"))
  
  var_sets <- vector("list", 61)
  annots <- c("S1", "S2", "S3")
  counter <- 1
  
  for (i in 0:20) {
    
    for (j in annots) {
      
      if (i == 0) {
        cluster_name <- "Base"
        if (counter > 1) {
          next
        }  
      } else {
        cluster_name <- paste0(j, "_cluster_", i)
      }
      
      cols <- c("SNP", cluster_name, "P")
      trained_snp_set <- trained_snps[ , ..cols]
      
      this_threshold <- this_trained_model[Set == cluster_name, Threshold]
      
      these_vars <- trained_snp_set[which(trained_snp_set[,2]==1 & trained_snp_set$P <= this_threshold)]$SNP
      
      if (cluster_name == "Base") {
        cluster_name <- "trained_base"
      }
      
      var_sets[[counter]] <- list(cluster_name, these_vars)
      counter <- counter + 1
      
    }
  }
  
  write_varset(var_sets, out_file)
}