## ---------------------------
## Purpose of script: 
##  Calculate statistics on important predictors chosen in elastic net.
## Author: Henri Chung
## ---------------------------

library(tidyverse)
library(glmnet)
library(parallel)
rm(list = ls())

# Set up data folder
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

sample_metadata <- mydata[["metadata"]]
sample_genotypes <- mydata[["genotypes"]] 
sample_phenotypes <- mydata[["phenotypes"]]

# gene metadata
gene_metadata <- reference$gene_metadata %>%
  mutate(reference$gene_metadata, combo = paste(gene, resistance_drug)) 

# tidy gene names to formal gene names
gene_key <- gene_metadata %>%
  select(gene_identifier) %>%
  unique() %>%
  mutate(gene = janitor::make_clean_names(gene_identifier))

# table for antibiotic class
drug_classes <- read_csv("data/reference/drug_classes_new.csv") %>%
  rename(mic_id = antibiotic) %>%
  mutate(class = ifelse(class == "sulfonamide", "sulphonamide", class),) %>%
  mutate(class = ifelse(class == "folate pathway antagonist", "FPA", class)) %>%
  mutate(class = gsub("beta lactam combo", "beta-lactam", class)) 

# stability selection
ss_data <- readRDS("outputs/ss_models.RDS")

# vector of host animals
host_animals <- c("cat", "cattle", "dog", "swine", "horse", "chicken", "turkey")

# Predictor Statistics
#####################

# the possible features fed into the model.
possible_predictors <- ss_data %>%
  dplyr::select(mic_id, breakpoint, training) %>%
  mutate(possible = map(training, function(.x){
    res <- .x %>%
      colnames() %>%
      as_tibble() 
    })) %>%
  dplyr::select(-training) %>%
  unnest(possible) %>%
  unique() %>%
  filter(value != "custom_phenotype") %>%
  rename(gene = "value") %>%
  mutate(possible = TRUE)


# identify important predictors (selected for in 2/3rds of 1000 replicates for the models.)
important_predictors <- ss_data %>%
  mutate(imp = purrr::map(coefs, function(.x){
    res <- .x %>%
      as.matrix() %>% 
      as.data.frame() %>%
      dplyr::filter(values >= 0.66) 
  })) %>%
  dplyr::select(c("mic_id","breakpoint", "imp")) %>%
  unnest(imp) %>%
  dplyr::select(-values) %>%
  unique() %>%
  filter(ind != "(Intercept)") %>%
  rename(gene = "ind") %>%
  mutate(important = TRUE)

# coeficient for important predictors
model_list <- readRDS("outputs/model_list.RDS")
bc_model <- model_list[["Binary-Count"]]
coef_predictors  <- bc_model %>%
  mutate(coef = purrr::map(en, function(.x){
    res <- .x %>% 
      coef() %>%
      as.matrix() %>% 
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      filter(gene != "(Intercept)") %>%
      mutate(s0 = abs(s0))
  })) %>%
  dplyr::select(c("mic_id","breakpoint", "coef")) %>%
  unnest(coef) 

# sum of predictor importance
coef_predictors_sum <- coef_predictors %>%
  group_by(mic_id, breakpoint) %>%
  summarize(coef_sum = sum(s0))

# sum of predictor importance by type of predictor
coef_predictors_by_group <- coef_predictors %>%
  left_join(important_predictors) %>%
  filter(important == TRUE) %>%
  select(-important) %>% 
  mutate(group = case_when(
    grepl(":", gene) ~ 'interactions',
    grepl(paste(host_animals, collapse = "|"), gene) ~ 'animal_main',
    gene == "binary" ~ 'binary',
    gene == "count" ~ 'count')) %>%
  mutate(group = ifelse(is.na(group), "gene_main", group)) %>%
  group_by(mic_id, breakpoint, group) %>% 
  summarize(coef = sum(s0)) %>% 
  pivot_wider(names_from = "group", values_from = "coef", values_fill = 0) %>%
  rowwise() %>%
  mutate(animal_all = sum(interactions, animal_main), gene_all = sum(gene_main, binary, count)) %>% 
  select(breakpoint, mic_id, gene_main, animal_main, binary, count, interactions, animal_all, gene_all) %>%
  arrange(breakpoint, mic_id) %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim_sulfamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>% 
  mutate(mic_id = ifelse(mic_id == "ticarcillin_clavulanic_acid", "co-ticarclav", mic_id)) %>%
  mutate(breakpoint = ifelse(breakpoint == "clsi", "CLSI", "ECOFF")) %>%
  mutate_if(is.numeric, signif, 2) %>%
  mutate_all(as.character)
coef_predictors_by_group[coef_predictors_by_group == 0] <- ""
colnames(coef_predictors_by_group) <- c("BP", "Antibiotic", "Genes", "Animal Main", "Binary", "Count", "Interactions", "Animal All", "Gene All")
write_csv(coef_predictors_by_group, "outputs/predictor_importance_by_group.csv")

# calculate the animal related and proportion of important predictors
animal_predictors <- important_predictors %>% 
  mutate(animal =  str_extract(pattern = paste(host_animals, collapse = "|"), string = gene)) %>%
  mutate(animal = ifelse(is.na(animal), "main", animal)) %>%
  group_by(mic_id, breakpoint) %>% 
  nest() %>%
  mutate(prop = purrr::map(data, function(.x){
    temp <- .x %>% mutate(animal2 = ifelse(animal == "main", "main", "animal")) %>%
    pull(animal2) %>% table()
    return(as_tibble(temp/sum(temp)))
    })) %>%
  unnest(prop) %>%
  pivot_wider(names_from = ".", values_from = "n") %>% 
  arrange(desc(animal))
animal_predictors

# calculate the animal related and proportion of possible predictors
animal_possible <- possible_predictors %>% 
  mutate(animal =  str_extract(pattern = paste(host_animals, collapse = "|"), string = gene)) %>%
  mutate(animal = ifelse(is.na(animal), "main", animal)) %>%
  group_by(mic_id, breakpoint) %>% 
  nest() %>%
  mutate(prop = purrr::map(data, function(.x){
    temp <- .x %>% mutate(animal2 = ifelse(animal == "main", "main", "animal")) %>%
    pull(animal2) %>% table()
    return(as_tibble(temp/sum(temp)))
    })) %>%
  unnest(prop) %>%
  pivot_wider(names_from = ".", values_from = "n") %>% 
  rename(animal_possible = "animal", main_possible = "main", data_possible = "data")

# combine possible and important animal proportions and compare
animal_pred_data <- left_join(animal_possible, animal_predictors, by = c("mic_id", "breakpoint")) %>%
  mutate(n_possible = purrr::map(data_possible, nrow), n = purrr::map(data, nrow)) %>%
  unnest(c(n_possible, n)) %>%
  select(-c("data", "data_possible"))

animal_pred_data %>% filter(is.na(animal)) %>% arrange(animal_possible) 
animal_pred_data %>% filter(!is.na(animal)) %>% arrange(animal)

# replace the "tidy" version of genes text with regular version.
important_predictors_table <- important_predictors %>% 
  left_join(gene_key) %>% 
  filter(!grepl(":", gene)) %>% 
  mutate(gene_identifier = ifelse(is.na(gene_identifier), gene, gene_identifier)) %>%
  select(-c("important", "gene")) %>% 
  group_by(mic_id, breakpoint) %>% 
  arrange(gene_identifier) %>%
  summarize(gene = paste(gene_identifier, collapse = ", ")) %>%
  arrange(breakpoint, mic_id) %>%
  select(breakpoint, mic_id, gene)
write_csv(important_predictors_table, "outputs/important_predictors.csv")

# identify predictors unique to an antibiotic relative to other models in its class.
unique_predictors <- important_predictors %>%
  left_join(drug_classes, by = "mic_id") %>%
  group_by(class, breakpoint, gene) %>%
  filter(n() == 1) %>%
  mutate(unique = TRUE) %>%
  ungroup() %>%
  dplyr::select(-c("class", "important")) 

# identify shared predictors between models of the same antibiotic class.
possible_shared_predictors <-  possible_predictors %>%
  left_join(drug_classes, by = "mic_id") %>%
  group_by(class, breakpoint) %>%
  mutate(l = length(unique(mic_id))) %>%
  group_by(class, breakpoint, gene) %>%
  mutate(n = n()) %>%
  filter(l == n & l != 1) %>%
  ungroup() %>%
  mutate(shared = TRUE) %>%
  group_by(class, breakpoint) %>% 
  summarize(possible = length(unique(gene)))

# identify shared predictors between models of the same antibiotic class.
shared_predictors <-  important_predictors %>%
  left_join(drug_classes, by = "mic_id") %>%
  group_by(class, breakpoint) %>%
  mutate(l = length(unique(mic_id))) %>%
  group_by(class, breakpoint, gene) %>%
  mutate(n = n()) %>%
  filter(l == n & l != 1) %>%
  ungroup() %>%
  mutate(shared = TRUE) %>%
  group_by(class, breakpoint) %>% 
  summarize(shared = length(unique(gene))) %>%
  left_join(possible_shared_predictors, by = c("class", "breakpoint")) %>%
  mutate(shared = paste(shared, possible, sep = "/"))

# aggregate the power values from power analysis into tables
#####################

# read the power data in.
power_sim <- read_csv("outputs/power_vals.csv")

# extract power values from power sim.
power_vals <- power_sim %>%
  dplyr::select(mic_id, breakpoint,gene, power, ef) %>%
  filter(gene != "(Intercept)") %>%
  filter(ef == 10) %>% dplyr::select(-ef)

power_predictors <- power_vals %>%
  mutate(power = ifelse(power >= 0.80, TRUE, FALSE)) 

# join predictor information into single df.
predictor_data <- possible_predictors %>%
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  left_join(unique_predictors, by = c("mic_id", "breakpoint", "gene")) %>% 
  left_join(power_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  left_join(coef_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  group_by(mic_id, breakpoint) %>%
  summarize(
    n_possible = sum(possible),
    n_important = sum(important, na.rm = TRUE),
    unique = sum(unique, na.rm = TRUE),
    n_power = sum(power, na.rm = TRUE),
    importance = signif(sum(s0, na.rm = TRUE),3)
  ) %>%
  mutate(percent_unique = signif(unique/n_possible, 2)) %>%
  left_join(drug_classes, by = "mic_id") %>%
  left_join(shared_predictors, by = c("class", "breakpoint")) 

# function to remove duplicates beyond first for latex conversion
replace_duplicates <- function(x) {
  seen <- c()
  output <- c()
  for (i in seq_along(x)) {
    if (x[i] %in% seen) {
      output[i] <- ""
    } else {
      seen <- c(seen, x[i])
      output[i] <- x[i]
    }
  }
  return(output)
}


# format table for import to latex.
predictor_table <- predictor_data %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim_sulfamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>% 
  mutate(mic_id = ifelse(mic_id == "ticarcillin_clavulanic_acid", "co-ticarclav", mic_id)) %>%
  mutate(class = ifelse(class == "folate pathway antagonist", "FPA", class)) %>% 
  mutate(breakpoint = ifelse(breakpoint == "clsi", "CLSI", "ECOFF")) %>%
  dplyr::select(breakpoint, class, mic_id, shared, unique, n_important, importance, n_power, n_possible) %>%
  unique() %>%
  arrange(breakpoint, class, mic_id, shared, unique, n_possible, n_power, n_important) %>%
  replace(is.na(.), "0") %>%
  mutate_all(as.character) 
predictor_table$breakpoint <- replace_duplicates(predictor_table$breakpoint)
predictor_table$class <- replace_duplicates(predictor_table$class)
predictor_table$shared <- replace_duplicates(predictor_table$shared)
predictor_table[predictor_table == 0] <- ""
colnames(predictor_table) <- c("BP", "Class", "Antibiotic", "Shared",  "Unique", "X", "Xpow", "Xhat", "importance")
write_csv(predictor_table,  "outputs/predictor_table.csv", na = "")

# Calculate the proportion of power values greater than or equal to .80 and were important.
power_dist <- power_vals %>% 
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  mutate(important = ifelse(is.na(important), FALSE, TRUE), power = ifelse(power >= 0.80, TRUE, FALSE)) %>%
  group_by(mic_id, breakpoint,power, important) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(ind = paste(power, important)) %>%
  dplyr::select(-c("power", "important")) %>%
  pivot_wider(names_from = "ind", values_from = n) %>%
  replace(is.na(.), NA) %>%
  rowwise() %>%
  mutate(n = sum(`FALSE FALSE`, `FALSE TRUE`, `TRUE FALSE`, `TRUE TRUE`, na.rm = TRUE)) 

# reformat table for import to latex.
power_table <- power_dist %>%
  left_join(drug_classes, by = "mic_id") %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim_sulfamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>%
  mutate(class = ifelse(mic_id == "folate pathway antagonist", "FPA", class)) %>%
  mutate(breakpoint = ifelse(breakpoint == "clsi", "CLSI", "ECOFF")) %>%
  dplyr::select(breakpoint, class, mic_id, `FALSE FALSE`, `FALSE TRUE`, `TRUE FALSE`, `TRUE TRUE`, everything()) %>%
  arrange(breakpoint, class, mic_id) %>%
  replace(is.na(.), 0) %>%
  mutate_all(as.character)
colnames(power_table) <- c("Breakpoint", "Class", "Antibiotic", "-/-", "-/+", "+/-", "+/+", "n")
write_csv(power_table, "outputs/power_table.csv", na = "") 

# calculate chi-square to test for effect of power on predictor importance.
#####################
cq_data_grad <- power_sim %>%
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  mutate(important = ifelse(is.na(important), FALSE, TRUE), power = ifelse(power >= 0.80, TRUE, FALSE))  %>%
  group_by(power, important, ef) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "important", values_from = "n") %>%
  group_by(ef) %>% group_split() %>%
  lapply(function(.x){chisq.test(select(.x, -ef))})

cq_data <-  power_vals %>% 
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  mutate(important = ifelse(is.na(important), FALSE, TRUE), power = ifelse(power >= 0.80, TRUE, FALSE))  %>%
  group_by(power, important) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "important", values_from = "n") 
chisq.test(cq_data)

# Calculate differences in predictor importance from WT and CLSI breakpoints.
#####################
breakpoint_differences <- important_predictors %>%
  group_by(mic_id) %>%
  mutate(l = length(unique(breakpoint))) %>%
  filter(l == 2) %>%
  group_by(mic_id, gene) %>%
  mutate(n = n()) %>%
  mutate(status = ifelse(n == 1, "unique", "shared")) %>%
  select(-c("n", "important")) %>% 
  group_by(mic_id, breakpoint, status) %>%
  summarize(n2 = n()) %>%
  rowwise() %>%
  mutate(status = ifelse(status == "unique", paste(breakpoint, status, sep = "_"), status)) %>%
  ungroup() %>% select(-breakpoint) %>%
  unique() %>%
  pivot_wider(names_from = "status", values_from = "n2")

# reformat breakpoint differences into table for import to latex.
breakpoint_table <- breakpoint_differences %>%
  left_join(drug_classes, by = "mic_id") %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim_sulfamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(class = ifelse(mic_id == "folate pathway antagonist", "FPA", class)) %>%
  select(class, mic_id, clsi_unique, ecoff_unique, shared) %>%
  arrange(class, mic_id) 
colnames(breakpoint_table) <- c("Class", "Antibiotic", "CLSI", "ECOFF", "Shared")
write_csv(breakpoint_table, "outputs/breakpoint_table.csv", na = "") 

# write important predictors (not-including interaction terms) in wide format
important_predictors_table <- important_predictors %>% 
  left_join(gene_key) %>% 
  filter(!grepl(":", gene)) %>% 
  mutate(gene_identifier = ifelse(is.na(gene_identifier), gene, gene_identifier)) %>%
  select(-c("important", "gene")) %>% 
  group_by(mic_id, breakpoint) %>% 
  arrange(gene_identifier) %>%
  summarize(gene = paste(gene_identifier, collapse = ", ")) %>%
  arrange(breakpoint, mic_id) %>%
  select(breakpoint, mic_id, gene)
write_csv(important_predictors_table, "outputs/important_predictors_table.csv")

# predictor table combining importance in model, 
power_wide <- power_sim %>% 
  mutate(ef = paste0("ef_", round(log(ef), 2))) %>%
  pivot_wider(names_from = "ef", values_from = power) %>%
  rename(beta = "s0") %>%
  mutate(abs_beta = abs(beta))

importance_power <- power_wide %>%
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  mutate(important = ifelse(is.na(important), FALSE, important))
write_csv(importance_power, "outputs/importance_power.csv")

model_res <- readRDS("outputs/model_res.RDS") %>% select(model, mic_id, breakpoint, .metric, .estimate)
write_csv(model_res, "outputs/model_res.csv")
quit()

model_res %>% 
  group_by(model, breakpoint, mic_id, .metric) %>% 
  summarize(.estimate = mean(.estimate, na.rm = TRUE)) %>% 
  filter(model == "Count" | model == "Only-Count" | model == "Base")  %>% 
  pivot_wider(names_from = .metric, values_from = .estimate) %>% 
  select(model, breakpoint, f_meas) %>% 
  arrange(mic_id) %>%
  pivot_wider(names_from = "model", values_from = "f_meas") %>% 
  mutate(diff = abs(`Only-Count`-Count)) %>%
  arrange(desc(diff))


model_res %>% 
  filter(mic_id %in% c("amoxicillin_clavulanic_acid", "ampicillin", "cephalexin")) %>% 
  filter(.metric == "f_meas" & breakpoint == "ecoff") %>% 
  select(model, mic_id, breakpoint, .estimate) %>%
  pivot_wider(names_from = "model", values_from = ".estimate")

model_res %>% 
  filter(model == "Binary-Count") %>%
  filter(.metric == "f_meas") %>% 
  group_by(mic_id) %>% filter(n() == 1) %>%
  filter(breakpoint == "ecoff") %>%
  select(model, mic_id, breakpoint, .estimate) %>%
  pivot_wider(names_from = "model", values_from = ".estimate") 