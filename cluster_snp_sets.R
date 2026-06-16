####################################
#Setup

library(ActivePathways)
library(purrr)
library(dplyr)
library(data.table)

home_dir <- "/path/to/home_directory"

source("path/to/functions.R")

####################################
#Inputs

#Results of EnrichmentMap / AutoAnnotate pathway clustering
cluster_file <- file.path(home_dir, "pathway_clusters.csv")

#Gene Ontology pathway annotation file containing all pathways of interest
#sourced from https://download.baderlab.org/EM_Genesets/ and filtered to exclude pathways with <10 or >1000 genes
gmt_file <- file.path(home_dir, "GO_filtered.gmt")

#S-1 annotation file
s1_file <- file.path(home_dir, "S_1_annotation.txt")

#S-2 annotation file
s2_file <- file.path(home_dir, "S_2_annotation.txt")

#S-3 annotation file
#annotation file created as in PMID: 36289406 with the addition of eQTL SNPs
s3_file <- file.path(home_dir, "S_3_annotation.txt")

#List of QC'd SNPs
qcd_snplist <- "path/to/qcd_snplist.txt"

#Specify output filename
out_file <- file.path(home_dir, "cluster_snp_sets.txt")


####################################
#Read inputs

clusters <- read.csv(cluster_file, header = T) %>%
  dplyr::rename(names = Pathway_Name, DESCR = Pathway_Description, Cluster_Num = X__mclCluster)

gmt <- read.GMT(gmt_file)


####################################
#Merge pathway clustering results and pathway annotations

#Format gmt file as a data frame
names <- names(gmt)
gmt_frame <- as.data.frame(do.call(rbind, gmt)) %>% cbind(., names)
rownames(gmt_frame) <- NULL

#Merging clusters with gmt to get list of all genes in pathways of interest
clusters <- left_join(clusters, gmt_frame, by = "names")
names(clusters$genes) <- NULL


####################################
#Filter the V2G annotation to genes of interest

#Non-overlapping list of all genes across all clusters
all_genes <- unique(unlist(clusters$genes)) %>% gsub(" ", "", .)


######################################
#For each annotation strategy, create sets of SNPs annotated each cluster

annotations <- list(
  S1 = s1_file,
  S2 = s2_file,
  S3 = s3_file
)

for (annotation_name in names(annotations)) {

  this_v2g <- read.GMT(annotations[[annotation_name]])
  
  #Filtering v2g annotation while removing any spaces in x$id
  filtered_v2g <- purrr::map(this_v2g, function(x) {
    x$id <- gsub(" ", "", x$id)
    x$name <- gsub(" ", "", x$name)
    x$genes <- gsub(" ", "", x$genes)
    x
  }) %>%
    purrr::keep(~ .x$id %in% all_genes)
  
  #Output SNPs annotated to clusters by cluster
  for (i in sort(unique(clusters$Cluster_Num))) {
    
    selected_cluster <- filter(clusters, Cluster_Num == i)
    
    genes_in_cluster <- unique(unlist(selected_cluster$genes))
    genes_in_cluster <- gsub(" ", "", genes_in_cluster)
    
    cluster_genes <- purrr::keep(filtered_v2g, ~.x$id %in% genes_in_cluster)
    snps <- unique(unlist(lapply(cluster_genes, function(x) x$genes)))
    
    output_file <- file.path(home_dir, paste0("annot_", annotation_name, "_cluster", i, "_SNPs.txt"))
    write.table(data.frame(snps), output_file, quote = F, row.names = F, col.names = F)

  }
}


######################################
#Format annotated SNPs for PRSice input

#Read in list of QC'd SNPs
qcsnps <- read.table(qcd_snplist)$V1

#Initialize a list to read in variant sets
var_sets <- vector("list", 60)
annots <- c("S1", "S2", "S3")

counter <- 1

for (i in 1:20) {
  
  for (j in annots) {
    
    cluster_name <- paste0(j, "_cluster_", i)
    these_vars <- read.table(paste0(home_dir, "/annot_", j, "_cluster", i, "_SNPs.txt"))$V1
    these_vars <- these_vars[these_vars %in% qcsnps]
    
    var_sets[[counter]] <- list(cluster_name, these_vars)
    counter <- counter + 1
  }
}

write_varset(var_sets, out_file)



