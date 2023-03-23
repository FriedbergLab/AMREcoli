# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(tidymodels)
library(janitor)
library(keras)
rm(list = ls())

#source("code/helper_functions.R")
message(Sys.time(), "| Setting up folder structure")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
message(Sys.time(), "| Reading in data")
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

# remove ducks
sample_metadata <- mydata[["metadata"]] %>%
  filter(host_animal_common != "ducks") %>% mutate(host_animal_common = as.character(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) 

sample_genotypes <- mydata[["genotypes"]] %>%
  filter(host_animal_common != "duck") %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog")))

sample_phenotypes <- mydata[["phenotypes"]] %>%
  filter(sample_id != "FL34741PPY20064")

identifiers <- read_csv(paste(dataFolder, "identifier_edited.csv", sep = ""))

gene_metadata <- reference$gene_metadata

pval = 0.05


# reformat genotype data
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_identifier) %>%
  unique() %>%
  mutate(values = 1)


# Check which method uncovers the most "expected" resistances given prior information
custom_gene_subset <- function(.x, gene_group){
  res <- .x$resistance %>%
    select(-card_id) %>%
    left_join(gene_metadata, by = "gene") %>%
    rename(gene_group = gene_group) %>%
    select(gene_group, gene_resistance) %>%
    rename(cardab_id = "gene_resistance") %>%
    left_join(reference$classes) %>%
    select(gene_group, mic_id) %>%
    mutate(combo = paste(gene_group, mic_id, sep = "_")) %>%
    unique() %>%
    filter(!is.na(mic_id))
}

gene_resistance <- custom_gene_subset(reference, "gene_identifier")
card_lookups_post <- read_csv("data/card_lookups_post.csv") 


# 1. Univariate modeling of genotype to binarized phenotypic output of clsi data.
# N = 7364 
# Subset dataset such that all I -> R
one_pheno <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# Same samples are significant in I > R and I > S
# what about for 50th quantile?
two_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# what about for 90th quantile?
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic90") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno, three_pheno)
names(data_phenotypes) <- c("ItoR", "ic90")

# Join and nest phenotype/genotype data for analysis
custom_nest <- function(pheno, genos, gene_group, animal = FALSE){
  
  res1 <- genos %>%
    rename(gene_group = gene_group) %>%
    mutate(gene_group = paste("g_", gene_group, sep = "")) %>%
    pivot_wider(names_from = "gene_group", values_from = "values") %>%
    replace(is.na(.), 0 ) %>% 
    pivot_longer(cols = contains("g_"), values_to = "values") %>%
    mutate(name = gsub("g_", "", name)) %>%
    rename(gene_group = "name") %>%
    left_join(pheno) %>%
    select(sample_id, host_animal_common, mic_id, phenotype, gene_group, values) %>%
    group_by(mic_id, gene_group) %>%
    nest() 
  
}
message(Sys.time(), "| Nesting model data")
model_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = FALSE)})

temp <- model_data[[1]] %>%
  unnest(data) %>% 
  unique() %>%
  filter(!is.na(mic_id))

model_data0 <- temp %>%
  group_by(sample_id) %>%
  rowwise() %>%
  pivot_wider(names_from = "gene_group", values_from = "values", values_fill = 0) %>%
  pivot_wider(names_from = "mic_id", values_from = "phenotype", values_fill = 0) %>%
  mutate(animal_val = 1) %>%
  pivot_wider(names_from = "host_animal_common", values_from = "animal_val", values_fill = 0) %>%
  clean_names() %>%
  mutate_if(is.numeric, as.factor)

ubi_genes <- names(model_data0[-1,])[apply(model_data0[-1,], 2, function(.x){length(unique(.x))}) == 1]

model_data1 <- select(model_data0, -all_of(ubi_genes))

# Fit model
message(Sys.time(), "| Starting ML workflow.")
set.seed(12345)


message(Sys.time(), "| Splitting data.")
data_split <- initial_split(model_data1, prop = 0.80)
train_data <- training(data_split)
test_data <- testing(data_split)

animals <- as.character(unique(temp$host_animal_common)) %>% make_clean_names()
genes <- colnames(model_data1)[colnames(model_data1) %in% (unique(temp$gene_group) %>% make_clean_names() )]
antibiotics <- unique(temp$mic_id) %>% make_clean_names()


message(Sys.time(), "| Crossfolding data.")
train_vfold <- vfold_cv(train_data, v = 5) %>% 
  mutate(df_ana = map(splits, analysis),
         df_ass = map(splits, assessment))

responses <- paste0(c(antibiotics, antibiotics), collapse = "+")
predictors <- paste0(genes, collapse = "+")
f <- paste0(responses, " ~ ", predictors, sep = "")  

message(Sys.time(), "| Writing recipe.")
data_rec <- recipes::recipe(x = train_data) %>%
  update_role(sample_id, new_role = "ID") %>%
  update_role(all_of(genes), new_role = "predictor") %>%
  update_role(all_of(c(animals, antibiotics)), new_role = "outcome")

data_rec <- prep(data_rec)

message(Sys.time(), "| Preparing model.")
data_model <-
  mlp(epochs = 100, hidden_units = 5, dropout = 0.1) %>%
  set_mode("classification") %>% 
  set_engine("keras", verbose = 0) 

message(Sys.time(), "| Combining workflow.")
data_workflow <- 
  workflow() %>% 
  add_model(data_model) %>% 
  add_recipe(data_rec)

message(Sys.time(), "| Fitting workflow to folds.")
model_fit <- train_vfold %>% 
    mutate(recipe = map (df_ana, ~prep(data_rec, training = .x)),
      df_ana = map (recipe,  juice),
        df_ass = map2(recipe, df_ass, ~bake(.x, new_data = .y))) %>% 
    mutate(model_fit  = map(df_ana, function(.x){data_workflow %>% fit(data = .x)})) %>% 
    mutate(model_pred = map2(model_fit, df_ass, ~predict(.x, new_data = .y)))