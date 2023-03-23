# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
library(lme4)
rm(list = ls())

#source("code/helper_functions.R")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
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

# Analysis Tests
# How many samples are there with an intermediate designation?
sample_intermediate <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  unique() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "phenotype", values_from = "value") %>%
  replace(is.na(.), 0 ) %>%
  left_join(select(sample_metadata, c("sample_id", "host_animal_common"))) %>%
  group_by(mic_id, host_animal_common) %>%
  summarize(S = sum(S), I = sum(I), R = sum(R)) %>%
  filter(I != 0) %>%
  arrange(desc(I)) 
#write_csv(sample_intermediate, "sample_intermediates.csv")


# reformat genotype data
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_identifier) %>%
  unique() %>%
  mutate(values = 1)
head(logreg_genos)

logreg_genos_family <- logreg_genos %>%
  left_join(select(gene_metadata, c("gene_identifier", "gene_name_family", "gene_type"))) %>%
  filter(gene_type != "plasmid") %>%
  unique() %>%
  group_by(host_animal_common, sample_id, gene_name_family) %>%
  summarize(values = sum(values)) %>%
  ungroup() %>%
  mutate(values = as.factor(values)) 
head(logreg_genos_family)

gene_metadata %>% select(-gene_resistance) %>% filter(gene_type != "plasmid") %>%unique() %>% View()
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
gene_resistance_family <- custom_gene_subset(reference, "gene_name_family")
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

# what about for 50th quantile?
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno,two_pheno, three_pheno)
names(data_phenotypes) <- c("ItoR","ic50", "ic50")


# Univariate Logistic Regression w/ Pooled animal groups

# fit a univariate logistic regression to data


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

model_data <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = FALSE)})
model_data_family <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos_family, "gene_name_family", animal = FALSE)})

# filter to only relationships found in CARD
custom_filter_nest <- function(.x, gene_resistance){
  res <- .x %>%
    ungroup() %>%
    mutate(combo = paste(gene_group, mic_id, sep = "_")) %>%
    filter(combo %in% gene_resistance$combo)
  return(res)
}



card_model_data <- lapply(model_data, function(.x){custom_filter_nest(.x, gene_resistance = gene_resistance)})
card_model_data_family <- lapply(model_data_family, function(.x){custom_filter_nest(.x, gene_resistance = gene_resistance_family)})
# Custom filter data to fit model requirements  
custom_filter <- function(x, filter = TRUE){
  res <- x %>%
    filter(!is.na(mic_id)) %>%
    mutate(p = map(data, function(.x){length(unique(.x$phenotype))})) %>% # data must have at least 2 phenotypes
    mutate(l = map(data, function(.x){length(unique(.x$values))})) %>% # data must have at least 2 different gene values
    mutate(m  = map(data, function(.x){mean(.x$values)})) %>% # data must have gene values that are not too close to 0 or 1
    mutate(a = map(data, function(.x){length(unique(.x$host_animal_common))})) %>% # there must be at least more than one animal in a group for animal effect
    unnest(p) %>% unnest(l) %>% unnest(m) %>% unnest(a)
  
  if(filter == TRUE){
    res <- res %>% filter(p > 1 & l > 1 & m < 0.98 & m > 0.02 & a > 1)
  }else{
    res <- res %>% filter(p <= 1 | l <= 1 | m >= 0.98 | m <= 0.02 | a <= 1)
  }
}

filtered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = TRUE)})
unfiltered_model_data <- lapply(card_model_data, function(.x){custom_filter(.x, filter = FALSE)})

filtered_model_data_family <- lapply(card_model_data_family, function(.x){custom_filter(.x, filter = TRUE)})
unfiltered_model_data_family <- lapply(card_model_data_family, function(.x){custom_filter(.x, filter = FALSE)})

# fit model to data
custom_test <- function(x){
  res <- x %>%
    mutate(fit = map(data, function(.x){return(lme4::glmer(phenotype~as.factor(values) + (1|host_animal_common), data = .x, family = binomial(link = "logit")))})) %>%
    mutate(summ = map(fit, function(.x){summary(.x)})) %>%
    mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
    mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
    mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
    unnest(pval.factor) %>%
    mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

filtered_results <- lapply(filtered_model_data, function(.x){custom_test(.x)})
filtered_results_family <- lapply(filtered_model_data_family, function(.x){custom_test(.x)})
# fit a separate model to data that doesn't have an animal random effect

custom_test_noanimal <- function(x){
  res <- x %>%
    mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
    mutate(summ = map(fit, function(.x){summary(.x)})) %>%
    mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
    mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
    mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
    unnest(pval.factor) %>%
    mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

unfiltered_model_data2 <- lapply(unfiltered_model_data, function(.x){filter(.x, p > 1 & l > 1 & m < 0.98 & m > 0.02 & a == 1)})
unfiltered_results <- lapply(unfiltered_model_data2, function(.x){custom_test_noanimal(.x)})

unfiltered_model_data2_family <- lapply(unfiltered_model_data_family, function(.x){filter(.x, p > 1 & l > 1 & m < 0.98 & m > 0.02 & a == 1)})
unfiltered_results_family <- lapply(unfiltered_model_data2_family, function(.x){custom_test_noanimal(.x)})

filtered_results
unfiltered_results

lapply(filtered_results, dim)
lapply(unfiltered_results, dim)

filtered_sig_results <- lapply(filtered_results, function(.x){filter(.x, pval.adj < 0.05)})
unfiltered_sig_results <- lapply(unfiltered_results, function(.x){filter(.x, pval.adj < 0.05)})

lapply(filtered_sig_results, dim)  
lapply(unfiltered_sig_results, dim)  


filtered_results_family
unfiltered_results_family

lapply(filtered_results_family, dim)
lapply(unfiltered_results_family, dim)

filtered_sig_results_family <- lapply(filtered_results_family, function(.x){filter(.x, pval.adj < 0.05)})
unfiltered_sig_results_family <- lapply(unfiltered_results_family, function(.x){filter(.x, pval.adj < 0.05)})

lapply(filtered_sig_results_family, dim)  
lapply(unfiltered_sig_results_family, dim) 


family_sort <- function(x){
  res <- x %>%
    mutate(pvals = map(pval, function(.x){as.data.frame(.x)})) %>%
    unnest(pvals) %>%
    pivot_longer(contains("("), names_to = "factor_level", values_to = "pvals", names_repair = "unique") %>%
    filter(!is.na(pvals)) %>%
    mutate(pval.adj = p.adjust(pval.factor, method = "holm"))
  return(res)
}

test <- lapply(filtered_results_family, family_sort)
filtered_sig_results_family <- lapply(test, function(.x){filter(.x, pval.adj < 0.05)})

filtered_sig_results_family[[1]] %>%
  select(mic_id, gene_group, factor_level, pval.adj, pvals) %>% View()
