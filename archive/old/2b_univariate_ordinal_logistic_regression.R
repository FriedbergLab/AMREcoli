# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
rm(list = ls())

source("code/helper_functions.R")
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

# reformat genotype data
logreg_genos <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_identifier) %>%
  unique() %>%
  mutate(values = 1)

# reformat genotype data
logreg_genos_family <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene) %>%
  left_join(gene_metadata) %>%
  select(host_animal_common, sample_id, gene_name_family, gene_identifier) %>%
  unique() %>%
  mutate(values = 1) %>%
  group_by(host_animal_common, sample_id, gene_name_family) %>%
  summarize(values = sum(values)) %>%
  filter(gene_name_family != "")
head(logreg_genos_family)

# Check which method uncovers the most "expected" resistances given prior information
gene_resistance <- custom_gene_subset(reference, "gene_identifier")
gene_resistance_family <- custom_gene_subset(reference, "gene_name_family")

head(gene_resistance)
head(gene_resistance_family)

# Ordinal Logistic Regression w/ Pooling animal groups
ord_phenotypes <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique()

pval = 0.05

data_gene_unpooled <- custom_nest(phenotypes, logreg_genos, "gene_identifier", animal = TRUE) 
test_gene_unpooled <- custom_test(data_gene_unpooled, ordinal = TRUE)
sig_gene_unpooled <- filter(test_gene_unpooled, pval.adj < !!pval)

data_gene_unpooled_filtered <- custom_filter_nest(data_gene_unpooled, gene_resistance) 
test_gene_unpooled_filtered <- custom_test(data_gene_unpooled_filtered, ordinal = TRUE)
sig_gene_unpooled_filtered <- filter(test_gene_unpooled_filtered, pval.adj < !!pval)

# Univariate Logistic Regression w/ Pooled animal groups
data_gene_pooled <- custom_nest(phenotypes, logreg_genos, "gene_identifier", animal = FALSE) 
test_gene_pooled <- custom_test(data_gene_pooled, ordinal = TRUE)
sig_gene_pooled <- filter(test_gene_pooled, pval.adj < !!pval)

data_gene_pooled_filtered <- custom_filter_nest(data_gene_pooled, gene_resistance) 
test_gene_pooled_filtered <- custom_test(data_gene_pooled_filtered, ordinal = TRUE)
sig_gene_pooled_filtered <- filter(test_gene_pooled_filtered, pval.adj < !!pval)


# Univariate Logistic Regression + Gene Families
data_family_pooled <- custom_nest(phenotypes, logreg_genos, "gene_identifier", animal = FALSE) 
test_family_pooled <- custom_test(data_family_pooled, ordinal = TRUE)
sig_family_pooled <- filter(test_family_pooled, pval.adj < !!pval)

data_family_pooled_filtered <- custom_filter_nest(data_family_pooled, gene_resistance_family) 
test_family_pooled_filtered <- custom_test(data_family_pooled_filtered, ordinal = TRUE)
sig_family_pooled_filtered <- filter(test_family_pooled_filtered, pval.adj < !!pval)

# Univariate Logistic Regression w/ Pooled animal groups + Gene Families
data_family_pooled <- custom_nest(phenotypes, logreg_genos, "gene_identifier", animal = FALSE) 
test_family_pooled <- custom_test(data_family_pooled, ordinal = TRUE)
sig_family_pooled <- filter(test_family_pooled, pval.adj < !!pval)

data_family_pooled_filtered <- custom_filter_nest(data_family_pooled, gene_resistance_family) 
test_family_pooled_filtered <- custom_test(data_family_pooled_filtered, ordinal = TRUE)
sig_family_pooled_filtered <- filter(test_family_pooled_filtered, pval.adj < !!pval)


my_data <- logreg_genos %>%
  mutate(gene_identifier = paste("g_", gene_identifier, sep = "")) %>%
  pivot_wider(names_from = "gene_identifier", values_from = "values") %>%
  replace(is.na(.), 0 ) %>% 
  pivot_longer(cols = contains("g_"), values_to = "values") %>%
  mutate(name = gsub("g_", "", name)) %>%
  rename(gene_identifier = "name") %>%
  left_join(five_pheno) %>%
  select(sample_id, host_animal_common, mic_id, phenotype, gene_identifier, values) %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  group_by(mic_id, gene_identifier) %>%
  nest() %>%
  filter(!is.na(mic_id)) %>%
  mutate(l  = map(data, function(.x){length(unique(.x$values))})) %>%
  unnest(l) %>%
  filter(l > 1)

pos_polr = possibly(.f = MASS::polr, otherwise = NA)

five_test <- five_data %>%
  mutate(fit = map(data, function(.x){pos_polr(formula = phenotype ~ as.factor(values), data = .x, Hess = TRUE)})) %>%
  filter(!is.na(fit)) %>%
  mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
  mutate(summ = map(fit, function(.x){summary(.x)})) %>%
  mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
  mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
  mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
  unnest(pval.factor) %>%
  mutate(pval.adj = p.adjust(pval.factor, method = "holm"))

# How many tests failed for ordinal regression?
nrow(five_data) - nrow(five_test)

# How many significant tests?
five_sig <- five_test %>% filter(pval.adj < pval)
nrow(five_sig)

# How many significant relationships recovered?
five_combo <- five_sig %>%
  ungroup() %>%
  select(gene_identifier, mic_id) %>%
  mutate(combo = paste(gene_identifier, mic_id, sep = "_")) %>%
  pull(combo) 
e <- five_combo %in% gene_resistance$combo
five_combo[e] %>% sort()
length(five_combo[e])

# filter and test only relevant relationships?
five_test_filtered <- five_data %>%
  ungroup() %>%
  mutate(combo = paste(gene_identifier, mic_id, sep = "_")) %>%
  filter(combo %in% gene_resistance$combo)  %>%
  mutate(fit = map(data, function(.x){pos_polr(formula = phenotype ~ as.factor(values), data = .x, Hess = TRUE)})) %>%
  filter(!is.na(fit)) %>%
  mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
  mutate(summ = map(fit, function(.x){summary(.x)})) %>%
  mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
  mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
  mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
  unnest(pval.factor) %>%
  mutate(pval.adj = p.adjust(pval.factor, method = "holm"))

five_sig_filtered <- five_test_filtered %>%
  filter(pval.adj < pval)  

nrow(five_test_filtered); nrow(five_sig_filtered)

# Ordinal Logistic Regression w/ Pooling animal groups

five_pheno <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique()

five_data_family <- logreg_genos_family %>%
  mutate(gene_name_family = paste("g_", gene_name_family, sep = "")) %>%
  pivot_wider(names_from = "gene_name_family", values_from = "values") %>%
  replace(is.na(.), 0 ) %>% 
  pivot_longer(cols = contains("g_"), values_to = "values") %>%
  mutate(name = gsub("g_", "", name)) %>%
  rename(gene_name_family = "name") %>%
  left_join(five_pheno) %>%
  select(sample_id, host_animal_common, mic_id, phenotype, gene_name_family, values) %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  group_by(mic_id, gene_name_family) %>%
  nest() %>%
  filter(!is.na(mic_id)) %>%
  mutate(l  = map(data, function(.x){length(unique(.x$values))})) %>%
  unnest(l) %>%
  filter(l > 1)

pos_polr = possibly(.f = MASS::polr, otherwise = NA)

five_test_family <- five_data_family %>%
  mutate(fit = map(data, function(.x){pos_polr(formula = phenotype ~ as.factor(values), data = .x, Hess = TRUE)})) %>%
  filter(!is.na(fit)) %>%
  mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
  mutate(summ = map(fit, function(.x){summary(.x)})) %>%
  mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
  mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
  mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
  unnest(pval.factor) %>%
  mutate(pval.adj = p.adjust(pval.factor, method = "holm"))

# How many tests failed for ordinal regression?
nrow(five_data_family) - nrow(five_test_family)

# How many significant tests?
five_sig_family <- five_test_family %>% filter(pval.adj < pval)
nrow(five_sig_family)

# How many significant relationships recovered?
five_combo_family <- five_sig_family %>%
  ungroup() %>%
  select(gene_name_family, mic_id) %>%
  mutate(combo = paste(gene_name_family, mic_id, sep = "_")) %>%
  pull(combo) 
e_family <- five_combo_family %in% gene_resistance_family$combo
five_combo_family[e_family] %>% sort()
length(five_combo_family[e_family])

# filter and test only relevant relationships?
five_test_filtered_family <- five_data_family %>%
  ungroup() %>%
  mutate(combo = paste(gene_name_family, mic_id, sep = "_")) %>%
  filter(combo %in% gene_resistance_family$combo)  %>%
  mutate(fit = map(data, function(.x){pos_polr(formula = phenotype ~ as.factor(values), data = .x, Hess = TRUE)})) %>%
  filter(!is.na(fit)) %>%
  mutate(fit = map(data, function(.x){return(glm(phenotype~as.factor(values), data = .x, family = binomial(link = "logit")))})) %>%
  mutate(summ = map(fit, function(.x){summary(.x)})) %>%
  mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
  mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% select('Pr(>|z|)') %>% t())})) %>%
  mutate(pval.factor = map(pval, function(.x){.x[,2]})) %>%
  unnest(pval.factor) %>%
  mutate(pval.adj = p.adjust(pval.factor, method = "holm"))

five_sig_filtered_family <- five_test_filtered_family %>%
  filter(pval.adj < pval)  

nrow(five_test_filtered_family); nrow(five_sig_filtered_family)