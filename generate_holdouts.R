###############################
#Setup

library(data.table)
library(dplyr)
library(splitTools)
library(tidysmd)


home_dir <- "/path/to/home_directory"

#################################
#Inputs

#UK Biobank included participants & covariates for the proxy, true case, and Wu/Marioni analyses
proxy_data <- "/path/to/proxy_data.txt"
true_data <- "/path/to/true_case_data.txt"
wumarioni_data <- "/path/to/wu_marioni_data.txt"

#UK Biobank bim file
#This should be the bim file to be used in PRSice jobs
bim_file <- "/path/to/ukb.bim"


#################################
#Specify inclusion set
inc_set <- "A"


#################################
#Read in set data & filter to included individuals

set_file <- case_when(inc_set == "A" ~ proxy_data,
                      inc_set == "B" ~ true_data,
                      inc_set == "C" ~ wumarioni_data)

set_col <- case_when(inc_set == "A" ~ "proxy_status",
                     inc_set == "B" ~"true_status",
                     inc_set == "C" ~"wu_marioni_status")

set_data <- fread(set_file) %>%
  filter(!is.na(!!sym(set_col))) %>%
  rename(status = !!sym(set_col))

#############################
#Vector of seeds
set.seed(2025)
seeds <- sample.int(.Machine$integer.max, 10, replace = FALSE)

############################
#Make train and test samples

#Set data stratification parameters
d <- set_data %>% select(status, sex)
k_num <- length(unique(paste0(d$status, d$sex)))
y <- multi_strata(d, k = k_num)

#Initialize data.frame to store SMDs for each split
all_smds <- data.frame(variable = character(),
                       method = character(),
                       group = character(),
                       smd = numeric()
                       )

outcome_names <- colnames(set_data)
start_cols <- length(outcome_names)


#Generate 10 training and testing sets
for (i in 1:10) {
  
  #Get seed for this run
  this_seed <- seeds[i]
  
  #Split data into training and testing using this seed
  inds <- partition(y, p = c(train = 0.8, test = 0.2), type = "stratified", split_into_list = T, seed = this_seed)
  set.train <- set_data[inds$train, ]
  set.test <- set_data[inds$test, ]
  
  #Add column names for these training and testing sets
  outcome_names <- c(outcome_names, paste0("AD_", inc_set, "_Tr", "_", i), paste0("AD_", inc_set, "_Ts", "_", i))
  
  
    outcome_train_vector <- ifelse(set_data$ID_1 %in% set.train$ID_1,
                                 set_data$AD,
                                 NA)
    
    outcome_test_vector <- ifelse(set_data$ID_1 %in% set.test$ID_1,
                                       set_data$AD,
                                       NA)

  set_data <- cbind(set_data, outcome_train_vector, outcome_test_vector)
  
  #Get SMDs of covariates across training and testing sets
  smd_data <- set_data[, 1:start_cols]
  
  smd_data <- smd_data %>%
    mutate(group = case_when(ID_1 %in% set.train$ID_1 ~ "Training",
                             ID_1 %in% set.test$ID_1 ~ "Testing"))
  
  #If inc_set is C then we have a different set of covariates to compare
  if (inc_set == "C") {
    this_tidy_smd <- as.data.frame(tidy_smd(smd_data, c(status, f.34.0.0, sex, e4_count, parents_affected), .group = group, na.rm = T))
  } else {
    this_tidy_smd <- as.data.frame(tidy_smd(smd_data, c(status, age, sex, e4_count, parents_affected), .group = group, na.rm = T))
  }
  
  rownames(this_tidy_smd) <- paste0(rownames(this_tidy_smd), ".", i)
  
  all_smds <- rbind(all_smds, this_tidy_smd)

}

names(set_data) <- outcome_names

split_cols <- paste0(
  "AD_", inc_set, "_",
  rep(c("Tr", "Ts"), times = 10),
  "_",
  rep(1:10, each = 2)
)

if (inc_set == "C") {
  set_data <- set_data %>%
    select(ID_1, ID_2, AD, sex, e4_count, f.22000.0.0, f.34.0.0, f_age, m_age, f.22009.0.1:f.22009.0.5, all_of(split_cols))
} else {
  set_data <- set_data %>%
    select(ID_1, ID_2, AD, age, sex, e4_count, f.22009.0.1:f.22009.0.5, all_of(split_cols))
} 

write.table(all_smds, file.path(home_dir, paste0("Set", inc_set, "_All_SMDs.txt")), quote = F, col.names = F, row.names = F)


######################
#Join with UKB bim file to make sure IDs are in the right order for PRSice

ukb_bim <- fread(bim_file, header = T)

bim_pheno <- left_join(ukb_bim, set_data, by = c("ID_1", "ID_2"))

write.table(bim_pheno, file.path(home_dir, paste0("Set", inc_set, "_PRSice_ready.txt")), quote = F, col.names = T, row.names = F, sep = "\t")




