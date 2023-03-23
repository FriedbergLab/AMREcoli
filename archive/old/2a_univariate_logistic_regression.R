# Load required packages and set up folder directory
setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(ontologyIndex)
library(lme4)
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
write_csv(sample_intermediate, "sample_intermediates.csv")


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

# I -> S
two_pheno <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = ifelse(phenotype == "I", "S", phenotype)) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

# Same samples are significant in I > R and I > S
# what about for 50th quantile?
three_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic50") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))


# what about for 90th quantile?
four_pheno <- sample_phenotypes %>%
  filter(breakpoint == "ic90") %>%
  select(sample_id, mic_id, phenotype) %>%
  unique() %>%
  mutate(phenotype = ifelse(phenotype == "R", 1, 0))

data_phenotypes <- list(one_pheno, two_pheno, three_pheno, four_pheno)
names(data_phenotypes) <- c("ItoR", "ItoS", "ic50", "ic90")


# Univariate Logistic Regression 
data_gene_unpooled <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = TRUE)})
test_gene_unpooled <- lapply(data_gene_unpooled, function(.x){custom_test(.x, ordinal = FALSE)})
sig_gene_unpooled <- lapply(test_gene_unpooled, function(.x){filter(.x, pval.adj < !!pval)})

data_gene_unpooled_filtered <- lapply(data_gene_unpooled, function(.x){custom_filter_nest(.x, gene_resistance)})
test_gene_unpooled_filtered <- lapply(data_gene_unpooled_filtered, function(.x){custom_test(.x, ordinal = FALSE)})
sig_gene_unpooled_filtered <- lapply(test_gene_unpooled_filtered, function(.x){filter(.x, pval.adj < !!pval)})

# Univariate Logistic Regression w/ Pooled animal groups
data_gene_pooled <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = FALSE)})
test_gene_pooled <- lapply(data_gene_pooled, function(.x){custom_test(.x, ordinal = FALSE)})
sig_gene_pooled <- lapply(test_gene_pooled, function(.x){filter(.x, pval.adj < !!pval)})

data_gene_pooled_filtered <- lapply(data_gene_pooled, function(.x){custom_filter_nest(.x, gene_resistance)})
test_gene_pooled_filtered <- lapply(data_gene_pooled_filtered, function(.x){custom_test(.x, ordinal = FALSE)})
sig_gene_pooled_filtered <- lapply(test_gene_pooled_filtered, function(.x){filter(.x, pval.adj < !!pval)})

# Univariate Logistic Regression + Gene Families
data_family_unpooled <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos_family, "gene_name_family", animal = TRUE)})
test_family_unpooled <- lapply(data_family_unpooled, function(.x){custom_test(.x, ordinal = FALSE)})
sig_family_unpooled <- lapply(test_family_unpooled, function(.x){filter(.x, pval.adj < !!pval)})

data_family_unpooled_filtered <- lapply(data_family_unpooled, function(.x){custom_filter_nest(.x, gene_resistance_family)})
test_family_unpooled_filtered <- lapply(data_family_unpooled_filtered, function(.x){custom_test(.x, ordinal = FALSE)})
sig_family_unpooled_filtered <- lapply(test_family_unpooled_filtered, function(.x){filter(.x, pval.adj < !!pval)})

# Univariate Logistic Regression w/ Pooled animal groups + Gene Families
data_family_pooled <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos_family, "gene_name_family", animal = FALSE)})
test_family_pooled <- lapply(data_family_pooled, function(.x){custom_test(.x, ordinal = FALSE)})
sig_family_pooled <- lapply(test_family_pooled, function(.x){filter(.x, pval.adj < !!pval)})

data_family_pooled_filtered <- lapply(data_family_pooled, function(.x){custom_filter_nest(.x, gene_resistance_family)})
test_family_pooled_filtered <- lapply(data_family_pooled_filtered, function(.x){custom_test(.x, ordinal = FALSE)})
sig_family_pooled_filtered <- lapply(test_family_pooled_filtered, function(.x){filter(.x, pval.adj < !!pval)})

# Ordinal 
ord_phenotypes <- sample_phenotypes %>%
  filter(breakpoint == "clsi") %>%
  mutate(phenotype = factor(phenotype, levels = c("S", "I", "R"))) %>%
  select(sample_id, mic_id, phenotype) %>%
  unique()

# Univariate Ordinal Logistic Regression 
ord_data_gene_unpooled <- custom_nest(ord_phenotypes, logreg_genos, "gene_identifier", animal = TRUE) 
ord_test_gene_unpooled <- custom_test(ord_data_gene_unpooled, ordinal = TRUE)
ord_sig_gene_unpooled <- filter(ord_test_gene_unpooled, pval.adj < !!pval)

ord_data_gene_unpooled_filtered <- custom_filter_nest(ord_data_gene_unpooled, gene_resistance) 
ord_test_gene_unpooled_filtered <- custom_test(ord_data_gene_unpooled_filtered, ordinal = TRUE)
ord_sig_gene_unpooled_filtered <- filter(ord_test_gene_unpooled_filtered, pval.adj < !!pval)

# Univariate Ordinal Logistic Regression w/ Pooled animal groups
ord_data_gene_pooled <- custom_nest(ord_phenotypes, logreg_genos, "gene_identifier", animal = FALSE) 
ord_test_gene_pooled <- custom_test(ord_data_gene_pooled, ordinal = TRUE)
ord_sig_gene_pooled <- filter(ord_test_gene_pooled, pval.adj < !!pval)

ord_data_gene_pooled_filtered <- custom_filter_nest(ord_data_gene_pooled, gene_resistance) 
ord_test_gene_pooled_filtered <- custom_test(ord_data_gene_pooled_filtered, ordinal = TRUE)
ord_sig_gene_pooled_filtered <- filter(ord_test_gene_pooled_filtered, pval.adj < !!pval)

# Univariate Ordinal Logistic Regression + Gene Families
ord_data_family_unpooled <- custom_nest(ord_phenotypes, logreg_genos_family, "gene_name_family", animal = TRUE) 
ord_test_family_unpooled <- custom_test(ord_data_family_unpooled, ordinal = TRUE)
ord_sig_family_unpooled <- filter(ord_test_family_unpooled, pval.adj < !!pval)

ord_data_family_unpooled_filtered <- custom_filter_nest(ord_data_family_unpooled, gene_resistance_family) 
ord_test_family_unpooled_filtered <- custom_test(ord_data_family_unpooled_filtered, ordinal = TRUE)
ord_sig_family_unpooled_filtered <- filter(ord_test_family_unpooled_filtered, pval.adj < !!pval)

# Univariate Ordinal Logistic Regression w/ Pooled animal groups + Gene Families
ord_data_family_pooled <- custom_nest(ord_phenotypes, logreg_genos_family, "gene_name_family", animal = FALSE) 
ord_test_family_pooled <- custom_test(ord_data_family_pooled, ordinal = TRUE)
ord_sig_family_pooled <- filter(ord_test_family_pooled, pval.adj < !!pval)

ord_data_family_pooled_filtered <- custom_filter_nest(ord_data_family_pooled, gene_resistance_family) 
ord_test_family_pooled_filtered <- custom_test(ord_data_family_pooled_filtered, ordinal = TRUE)
ord_sig_family_pooled_filtered <- filter(ord_test_family_pooled_filtered, pval.adj < !!pval)

# compile results into table

results <- as.data.frame(matrix(nrow = 24, ncol = 5))
colnames(results) <- c("ItoR", "ItoS", "ic50", "ic90", "ord")

results[1,] <- c(lapply(data_gene_unpooled, nrow) %>% unlist(), nrow(ord_data_gene_unpooled))
results[2,] <- c(lapply(test_gene_unpooled, nrow) %>% unlist(), nrow(ord_test_gene_unpooled))
results[3,] <- c(lapply(sig_gene_unpooled, nrow) %>% unlist(), nrow(ord_sig_gene_unpooled))

results[4,] <- c(lapply(data_gene_unpooled_filtered, nrow) %>% unlist(), nrow(ord_data_gene_unpooled_filtered))
results[5,] <- c(lapply(test_gene_unpooled_filtered, nrow) %>% unlist(), nrow(ord_test_gene_unpooled_filtered))
results[6,] <- c(lapply(sig_gene_unpooled_filtered, nrow) %>% unlist(), nrow(ord_sig_gene_unpooled_filtered))

results[7,] <- c(lapply(data_gene_pooled, nrow) %>% unlist(), nrow(ord_data_gene_pooled))
results[8,] <- c(lapply(test_gene_pooled, nrow) %>% unlist(), nrow(ord_test_gene_pooled))
results[9,] <- c(lapply(sig_gene_pooled, nrow) %>% unlist(), nrow(ord_sig_gene_pooled))

results[10,] <- c(lapply(data_gene_pooled_filtered, nrow) %>% unlist(), nrow(ord_data_gene_pooled_filtered))
results[11,] <- c(lapply(test_gene_pooled_filtered, nrow) %>% unlist(), nrow(ord_test_gene_pooled_filtered))
results[12,] <- c(lapply(sig_gene_pooled_filtered, nrow) %>% unlist(), nrow(ord_sig_gene_pooled_filtered))

results[13,] <- c(lapply(data_family_unpooled, nrow) %>% unlist(), nrow(ord_data_family_unpooled))
results[14,] <- c(lapply(test_family_unpooled, nrow) %>% unlist(), nrow(ord_test_family_unpooled))
results[15,] <- c(lapply(sig_family_unpooled, nrow) %>% unlist(), nrow(ord_sig_family_unpooled))

results[16,] <- c(lapply(data_family_unpooled_filtered, nrow) %>% unlist(), nrow(ord_data_family_unpooled_filtered))
results[17,] <- c(lapply(test_family_unpooled_filtered, nrow) %>% unlist(), nrow(ord_test_family_unpooled_filtered))
results[18,] <- c(lapply(sig_family_unpooled_filtered, nrow) %>% unlist(), nrow(ord_sig_family_unpooled_filtered))

results[19,] <- c(lapply(data_family_pooled, nrow) %>% unlist(), nrow(ord_data_family_pooled))
results[20,] <- c(lapply(test_family_pooled, nrow) %>% unlist(), nrow(ord_test_family_pooled))
results[21,] <- c(lapply(sig_family_pooled, nrow) %>% unlist(), nrow(ord_sig_family_pooled))

results[22,] <- c(lapply(data_family_pooled_filtered, nrow) %>% unlist(), nrow(ord_data_family_pooled_filtered))
results[23,] <- c(lapply(test_family_pooled_filtered, nrow) %>% unlist(), nrow(ord_test_family_pooled_filtered))
results[24,] <- c(lapply(sig_family_pooled_filtered, nrow) %>% unlist(), nrow(ord_sig_family_pooled_filtered))

results$class <- rep(c("data", "test", "sig"), 8)
results$gene_type <- c(rep("gene", 12), rep("family", 12))
results$pool <- rep(c(rep("unpooled", 6), rep("pooled", 6)), 2)
results$filtered <- rep(c(rep("unfiltered",3), rep("filtered", 3)), 4)

custom_compare <- function(x, y){
  names <- c("ItoR", "ItoS", "ic50", "ic90", "ord")
  a <- c(lapply(x, nrow) %>% unlist(), nrow(y))
  b <- c(lapply(x, function(.x){custom_combo(.x, gene_resistance) %>% length()}) %>% unlist(),
         custom_combo(y, gene_resistance) %>% length())
  res <- bind_rows(a, b)
  colnames(res) <- names
  res$filter <- c(FALSE, TRUE)
  return(res)
}

custom_compare_family <- function(x, y){
  names <- c("ItoR", "ItoS", "ic50", "ic90", "ord")
  a <- c(lapply(x, nrow) %>% unlist(), nrow(y))
  b <- c(lapply(x, function(.x){custom_combo(.x, gene_resistance_family) %>% length()}) %>% unlist(),
         custom_combo(y, gene_resistance_family) %>% length())
  res <- bind_rows(a, b)
  colnames(res) <- names
  res$filter <- c(FALSE, TRUE)
  return(res)
}

combo_sig_gene_unpooled <- custom_compare(sig_gene_unpooled, ord_sig_gene_unpooled)
combo_sig_gene_pooled <- custom_compare(sig_gene_pooled, ord_sig_gene_pooled)

combo_sig_family_unpooled <- custom_compare_family(sig_family_unpooled, ord_sig_family_unpooled)
combo_sig_family_pooled <- custom_compare_family(sig_family_pooled, ord_sig_family_pooled)


# Univariate Mixed Effects Logistic Regression model with Animal Random Effects
data_gene_pooled_mixed <- lapply(data_gene_pooled, custom_filter_animal)
test_gene_pooled_mixed <- lapply(data_gene_pooled_mixed, function(.x){custom_test_mixed(.x)})
sig_gene_pooled_mixed <- lapply(test_gene_pooled_mixed, function(.x){filter(.x, pval.adj < !!pval)})

data_gene_pooled_mixed_filtered <- lapply(data_gene_pooled, function(.x){custom_filter_nest_animal(.x, gene_resistance)})
test_gene_pooled_mixed_filtered <- lapply(data_gene_pooled_mixed_filtered, function(.x){custom_test_mixed(.x)})
sig_gene_pooled_mixed_filtered <- lapply(test_gene_pooled_mixed_filtered, function(.x){filter(.x, pval.adj < !!pval)})


# Univariate Mixed Effects Logistic Regression model with Animal Random Effects w/ Gene Families
data_family_pooled_mixed <- lapply(data_family_pooled, custom_filter_animal)
test_family_pooled_mixed <- lapply(data_family_pooled_mixed, function(.x){custom_test_mixed(.x)})
sig_family_pooled_mixed <- lapply(test_family_pooled_mixed, function(.x){filter(.x, pval.adj < !!pval)})

data_family_pooled_mixed_filtered <- lapply(data_family_pooled, function(.x){custom_filter_nest_animal(.x, gene_resistance_family)})
test_family_pooled_mixed_filtered <- lapply(data_family_pooled_mixed_filtered, function(.x){custom_test_mixed(.x)})
sig_family_pooled_mixed_filtered <- lapply(test_family_pooled_mixed_filtered, function(.x){filter(.x, pval.adj < !!pval)})


#compare models by summing AIC and BIC values from significant results
rm(list = ls(pattern = "data_.*"))
rm(all_sigs)
all_sigs <- mget(ls(pattern = "^sig.*"), envir = globalenv())
names(all_sigs)

custom_aic <- function(x){
  res <- lapply(x[["fit"]], function(.x){AIC(logLik(.x))})
  return(res)
}

all_aic <- lapply(all_sigs, function(.y){lapply(.y, custom_aic)}) %>%
  lapply(function(.z){lapply(.z, unlist) %>% stack()}) %>%
  bind_rows(.id = "model")

all_aic %>%
  group_by(model, ind) %>%
  summarize(avg = mean(values), n = n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = "ind", values_from = "avg") %>%
  arrange(avg)

all_ords <- mget(ls(pattern = "^ord_sig.*"), envir = globalenv())

all_ord_aic <- lapply(all_ords, custom_aic) %>%
  stack() %>% 
  group_by(ind) %>%
  summarize(avg = mean(values)) %>%
  arrange(avg)
head(all_ord_aic)

#

all_tests <-  mget(ls(pattern = "test_.*"), envir = globalenv())[-1]
names(all_tests) %>% sort()

all_tests <- list(test_gene_unpooled, 
                  test_gene_pooled, 
                  test_family_unpooled, 
                  test_family_pooled, 
                  test_gene_pooled_mixed, 
                  test_family_pooled_mixed)

names(all_tests) <- c("test_gene_unpooled", "test_gene_pooled", "test_family_unpooled", 
                               "test_family_pooled", "test_gene_pooled_mixed", "test_family_pooled_mixed")

all_tests_filtered <- list(test_gene_unpooled_filtered, 
                  test_gene_pooled_filtered, 
                  test_family_unpooled_filtered, 
                  test_family_pooled_filtered, 
                  test_gene_pooled_mixed_filtered, 
                  test_family_pooled_mixed_filtered)

names(all_tests_filtered) <- c("test_gene_unpooled_filtered", "test_gene_pooled_filtered", "test_family_unpooled_filtered", 
                               "test_family_pooled_filtered", "test_gene_pooled_mixed_filtered", "test_family_pooled_mixed_filtered")

all_tests_ord <- list(ord_test_gene_unpooled, 
                      ord_test_family_unpooled, 
                      ord_test_gene_pooled, 
                      ord_test_family_unpooled)

names(all_tests_ord) <- c("ord_tests_gene_unpooled", "ord_test_family_unpooled", "ord_test_gene_pooled", "ord_test_family_unpooled")

all_tests_ord_filtered <- list(ord_test_gene_unpooled_filtered, 
                      ord_test_family_unpooled_filtered, 
                      ord_test_gene_pooled_filtered, 
                      ord_test_family_pooled_filtered)

names(all_tests_ord_filtered) <- c("ord_tests_gene_unpooled_filtered", "ord_test_family_unpooled_filtered", 
                          "ord_test_gene_pooled_filtered", "ord_test_family_pooled_filtered")

# number of tests
lapply(all_tests, function(.x){
  lapply(.x, nrow) 
}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_filtered, function(.x){
  lapply(.x, nrow) 
}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_ord, nrow) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

lapply(all_tests_ord_filtered, nrow) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

# AIC values (mean)
lapply(all_tests, function(.x){
  lapply(.x, function(.y){
    .y %>%
      mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
      unnest(AIC) %>%
      pull(AIC) %>%
      mean()
  })
}) %>%
  unlist() %>%
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_filtered, function(.x){
  lapply(.x, function(.y){
    .y %>%
      mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
      unnest(AIC) %>%
      pull(AIC) %>%
      mean()
  })
}) %>%
  unlist() %>%
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_ord,function(.y){
  .y %>%
    mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
    unnest(AIC) %>%
    pull(AIC) %>%
    mean()}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

lapply(all_tests_ord_filtered,function(.y){
  .y %>%
    mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
    unnest(AIC) %>%
    pull(AIC) %>%
    mean()}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

# AIC values (median)
lapply(all_tests, function(.x){
  lapply(.x, function(.y){
    .y %>%
      mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
      unnest(AIC) %>%
      pull(AIC) %>%
      median()
  })
}) %>%
  unlist() %>%
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_filtered, function(.x){
  lapply(.x, function(.y){
    .y %>%
      mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
      unnest(AIC) %>%
      pull(AIC) %>%
      median()
  })
}) %>%
  unlist() %>%
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_tests_ord,function(.y){
  .y %>%
    mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
    unnest(AIC) %>%
    pull(AIC) %>%
    median()}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

lapply(all_tests_ord_filtered,function(.y){
  .y %>%
    mutate(AIC = map(fit, function(.z){AIC(logLik(.z))})) %>%
    unnest(AIC) %>%
    pull(AIC) %>%
    median()}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

# 
all_tests <-  mget(ls(pattern = "test_.*"), envir = globalenv())[-1]
names(all_tests) %>% sort()

all_sigs <- list(sig_gene_unpooled, 
                  sig_gene_pooled, 
                  sig_family_unpooled, 
                  sig_family_pooled, 
                  sig_gene_pooled_mixed, 
                  sig_family_pooled_mixed)

names(all_sigs) <- c("sig_gene_unpooled", "sig_gene_pooled", "sig_family_unpooled", 
                      "sig_family_pooled", "sig_gene_pooled_mixed", "sig_family_pooled_mixed")

all_sigs_filtered <- list(sig_gene_unpooled_filtered, 
                           sig_gene_pooled_filtered, 
                           sig_family_unpooled_filtered, 
                           sig_family_pooled_filtered, 
                           sig_gene_pooled_mixed_filtered, 
                           sig_family_pooled_mixed_filtered)

names(all_sigs_filtered) <- c("sig_gene_unpooled_filtered", "sig_gene_pooled_filtered", "sig_family_unpooled_filtered", 
                               "sig_family_pooled_filtered", "sig_gene_pooled_mixed_filtered", "sig_family_pooled_mixed_filtered")

all_sigs_ord <- list(ord_sig_gene_unpooled, 
                      ord_sig_family_unpooled, 
                      ord_sig_gene_pooled, 
                      ord_sig_family_unpooled)

names(all_sigs_ord) <- c("ord_sigs_gene_unpooled", "ord_sig_family_unpooled", "ord_sig_gene_pooled", "ord_sig_family_unpooled")

all_sigs_ord_filtered <- list(ord_sig_gene_unpooled_filtered, 
                               ord_sig_family_unpooled_filtered, 
                               ord_sig_gene_pooled_filtered, 
                               ord_sig_family_unpooled_filtered)

names(all_sigs_ord_filtered) <- c("ord_sigs_gene_unpooled_filtered", "ord_sig_family_unpooled_filtered", 
                                   "ord_sig_gene_pooled_filtered", "ord_sig_family_unpooled_filtered")

# significant tests
lapply(all_sigs, function(.x){
  lapply(.x, nrow) 
}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_sigs_filtered, function(.x){
  lapply(.x, nrow) 
}) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() %>%
  separate(ind, into = c("model", "data"), sep = "\\.") %>%
  pivot_wider(names_from = "data", values_from = "values")

lapply(all_sigs_ord, nrow) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

lapply(all_sigs_ord_filtered, nrow) %>% 
  unlist() %>% 
  stack() %>% 
  as.data.frame() 

all_sigs_ord[[1]]

names(all_sigs_ord)
custom_combo(all_sigs_ord[[1]], gene_resistance) %>% length()
custom_combo(all_sigs_ord[[2]], gene_resistance_family) %>% length()
custom_combo(all_sigs_ord[[3]], gene_resistance) %>% length()
custom_combo(all_sigs_ord[[4]], gene_resistance_family) %>% length()

custom_combo(all_sigs_ord_filtered[[1]], gene_resistance) %>% length()
custom_combo(all_sigs_ord_filtered[[2]], gene_resistance_family) %>% length()
custom_combo(all_sigs_ord_filtered[[3]], gene_resistance) %>% length()
custom_combo(all_sigs_ord_filtered[[4]], gene_resistance_family) %>% length()

all_sigs %>% names()

lapply(all_sigs[[1]], function(.y){custom_combo(.y, gene_resistance) %>% length()}) %>% stack()
lapply(all_sigs[[2]], function(.y){custom_combo(.y, gene_resistance) %>% length()}) %>% stack()
lapply(all_sigs[[3]], function(.y){custom_combo(.y, gene_resistance_family) %>% length()}) %>% stack()
lapply(all_sigs[[4]], function(.y){custom_combo(.y, gene_resistance_family) %>% length()}) %>% stack()
lapply(all_sigs[[5]], function(.y){custom_combo(.y, gene_resistance) %>% length()}) %>% stack()
lapply(all_sigs[[6]], function(.y){custom_combo(.y, gene_resistance_family) %>% length()}) %>% stack()
