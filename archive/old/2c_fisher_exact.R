# Fisher exact is not appropriate when Expected Value > 5
# Bower, Keith M. 2003. "When to Use Fisher's Exact Test." In American Society for Quality, Six Sigma Forum Magazine, 2:35-37. 4.
# McCrum-Gardner, Evie. 2008. "Which Is the Correct Statistical Test to Use?" British Journal of Oral and Maxillofacial Surgery 46 (1): 38-41.

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

data_phenotypes

# Univariate Logistic Regression 
data_gene_unpooled <- lapply(data_phenotypes, function(.x){custom_nest(.x, logreg_genos, "gene_identifier", animal = TRUE)})

custom_contigency <- function(x){
  res <- x %>%
    mutate(phenotype = ifelse(phenotype == 0, "S", "R")) %>%
    rename("gene" = values) %>%
    mutate(gene = ifelse(gene == 0, "0", "1")) 
  return(res)
}

chi_test_gene_unpooled <- lapply(data_gene_unpooled, function(.x){
  res <- mutate(.x, data = map(data, function(.y){custom_contigency(.y)})) %>%
    mutate(chi_square = map(data, function(.y){return(chisq.test(table(.y$phenotype, .y$gene)))})) %>%
    mutate(pval = map(chi_square, function(.y){return(.y$p.value)})) %>%
    mutate(pval.adj = p.adjust(pval, method = "holm"))
  })

chi_sig_gene_unpooled <- lapply(chi_test_gene_unpooled, function(.x){filter(.x, pval.adj < 0.05)})

