#############################################################
#
write_v2g_file <- function(gmt_data, file_path) {
  if (length(gmt_data) == 0) {
    cat("Error: No data provided.\n")
    return()
  }
  
  # Open the file connection for writing
  file_conn <- file(file_path, "w")
  
  # Write each gene set to the file
  for (i in seq_along(gmt_data)) {
    # Extracting gene set information
    id <- gsub(" ", "_", gmt_data[[i]][["id"]])  
    name <- gsub(" ", "_", gmt_data[[i]][["name"]]) 
    genes <- gmt_data[[i]][["genes"]]
    
    # Writing gene set information to file
    cat(id, "\t", name, "\t", paste(genes, collapse = "\t"), "\n", file = file_conn)
  }
  
  # Close the file connection
  close(file_conn)
  
  cat("V2G file written successfully.\n")
}






############################################################
#Function to write variant sets in format needed for PRSice
write_varset <- function(gmt_data, file_path) {
  if (length(gmt_data) == 0) {
    cat("Error: No data provided.\n")
    return()
  }
  
  # Open the file connection for writing
  file_conn <- file(file_path, "w")
  
  # Write each gene set to the file
  for (i in seq_along(gmt_data)) {
    
    snp_field <- length(gmt_data[[1]])
    # Extracting gene set information
    id <- gmt_data[[i]][[1]]  
    variants <- gmt_data[[i]][[2]]  
    
    # Writing gene set information to file
    cat(id, "\t", paste(variants, collapse = "\t"), "\n", file = file_conn)
  }
  
  # Close the file connection
  close(file_conn)
  
  cat("Variant-set file written successfully.\n")
}
