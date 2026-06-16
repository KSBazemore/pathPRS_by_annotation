####################################
#Setup

library(ActivePathways)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

home_dir <- "/path/to/home_directory"

source(file.path(home_dir, "functions.R")) 

####################################
#Inputs

#Positional annotation file from MAGMA
annot_file <- file.path(home_dir, "S_1_annotation.genes.annot")

#Chromatin interaction mapping file, from FUMA on QC'd summary statistics
ci_file <- "/path/to/ci.txt"

#eQTL mapping file, from FUMA on QC'd summary statistics 
eqtl_file <- "/path/to/eqtl.txt"

#Created file
annotation_output <- file.path(home_dir, "S_2_annotation.txt")

#####################################
#Read inputs

annot <- read.GMT(annot_file)

ci <- fread(ci_file, header = TRUE)
eqtl <- fread(eqtl_file, header = TRUE)


#####################################
#Process the chromatin interaction file to yield a list of genes, each with a list of corresponding snps

#Individual rows in the source file contain multiple genes, splitting these up 
ci <- separate_rows(ci, genes, sep = ":")

#Reformat SNP list in each row so it can be read as a vector
ci$rsids <- str_split(ci$SNPs, ';')

#Get a list of unique genes in the ci file
ci <- filter(ci, !is.na(genes))
ci_genes <- unique(as.vector(ci$genes))

#Make the final nested list of genes, each with a list of their mapped snps
ci_v2g <- list()

for (i in seq_along(ci_genes)) {
  ci_v2g[[i]] <- list()
  ci_v2g[[i]]$gene <- ci_genes[i]
  filter_ci <- filter(ci, genes == ci_genes[i])
  ci_v2g[[i]]$snps <- unique(unlist(filter_ci$rsids))
}

names(ci_v2g) <- sapply(ci_v2g, function(x) x$gene)


#####################################
#Add SNPs from CI to the appropriate genes in the annotation file
#Skip over the first two elements as they specify windows

for (i in seq_along(annot)[-c(1, 2)]) {
  gene <- annot[[i]]$id
  
  if (length(gene) ==1 && (gene %in% names(ci_v2g)) == T) {
    annot[[i]]$genes <- unique(c(annot[[i]]$genes, ci_v2g[[gene]]$snps))
  } 
}


#####################################
#Process the eQTL file to yield a list of genes, each with a list of corresponding snps

#Identify list of unique genes with eqtl snps
eqtl_genes <- unique(as.vector(eqtl$gene))

#Make the final nested list of genes, each with a list of their mapped snps
eqtl_v2g <- list()

for (i in seq_along(eqtl_genes)) {
  eqtl_v2g[[i]] <- list()
  eqtl_v2g[[i]]$gene <- eqtl_genes[i]
  filter_eqtl <- filter(eqtl, gene == eqtl_genes[i])
  eqtl_v2g[[i]]$snps <- unique(unlist(filter_eqtl$uniqID))
}

names(eqtl_v2g) <- sapply(eqtl_v2g, function(x) x$gene)


#####################################
#Add snps from eqtl to the appropriate genes in gmt and above vector
#Skip over the first two elements as they specify windows

for (i in seq_along(annot)[-c(1, 2)]) {
  gene <- annot[[i]]$id
  
  if (length(gene) ==1 && (gene %in% names(eqtl_v2g)) == T) {
    annot[[i]]$genes <- unique(c(annot[[i]]$genes, eqtl_v2g[[gene]]$snps)) 
  } 
}

write_v2g_file(annot, annotation_output)
