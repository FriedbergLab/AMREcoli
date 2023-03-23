## ---------------------------
## Purpose of script: 
##  Generate plots and tables describing experimental data.
##  Each plot aims to answer a basic question about the dataset.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("~/amr_analysis")
library(tidyverse)
library(glmnet)
library(rstatix)
library(yardstick)
rm(list = ls())
set.seed(123)

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


# What are the most prevalent genes among samples? (n = 846)
t1d <- sample_genotypes %>%
  select(host_animal_common, sample_id, gene, gene_type) %>%
  unique() 

t1 <- t1d %>%
  select(sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum)) %>%
  mutate(prop = sum/nrow(sample_metadata))
t1
write_csv(t1, "outputs/total_gene_counts.csv")

# Most prevalent genes per animal.
t1b <- t1d %>%
  dplyr::select(host_animal_common, sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(host_animal_common, gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum)) %>%
  pivot_wider(names_from = host_animal_common, values_from = sum)
t1b
write_csv(t1b, "outputs/total_gene_counts_animal.csv")

# Most prevalent genes among samples by animal as a percentage.
sample_counts <- sample_genotypes %>%
  select(host_animal_common, sample_id) %>%
  unique() %>%
  pull(host_animal_common) %>%
  table()

t2 <- t1d %>%
  select(host_animal_common, sample_id, gene, gene_type) %>%
  unique() %>%
  group_by(host_animal_common, gene, gene_type) %>%
  summarize(sum = n()) %>%
  ungroup() %>%
  arrange(desc(sum)) %>% 
  filter(!is.na(gene)) %>%
  pivot_wider(names_from = host_animal_common, values_from = sum) %>%
  replace(is.na(.), 0 ) %>%
  mutate(cattle = cattle/sample_counts[["cattle"]],
         swine = swine/sample_counts[["swine"]], 
         chicken = chicken/sample_counts[["chicken"]], 
         turkey = turkey/sample_counts[["turkey"]],
         horse = horse/sample_counts[["horse"]],
         dog = dog/sample_counts[["dog"]],
         cat = cat/sample_counts[["cat"]])
t2  
write_csv(t2, "outputs/gene_presence_per_animal.csv")

# Number of samples of each animal species.
p1d <- sample_metadata %>% 
  select(sample_id, host_animal_common) %>%
  unique() %>%
  group_by(host_animal_common) %>%
  summarize(n = n()) %>%
  ungroup()

p1 <- p1d %>%
  ggplot(aes(x = reorder(host_animal_common, desc(n)), y = n)) + 
  geom_bar(stat = "identity")  + 
  xlab("") + ylab("count") + 
  ylim(0, 250) + 
  labs(title = "Samples per Host Animal", subtitle = "N = 980", fill = "Host Animal") +
  geom_text(aes(label=n), vjust = -2) +
  theme_classic()
p1

pdf("outputs/sample_histogram_animal.pdf")
p1
dev.off()

# Number of samples of each animal species. (species names)
p2d <- sample_metadata %>% 
  select(sample_id, host_animal_species) %>%
  unique() %>%
  group_by(host_animal_species) %>%
  summarize(n = n()) %>%
  ungroup()

p2 <- p2d %>%
  ggplot(aes(x = reorder(host_animal_species, desc(n)), y = n, fill = host_animal_species)) + 
  geom_bar(stat = "identity")  + 
  xlab("") + ylab("count") + 
  ylim(0, 250) + 
  labs(title = "Samples per Host Animal", subtitle = "N = 980", fill = "Host Animal") +
  geom_text(aes(label=n), vjust = -2) +
  theme_classic()
p2

pdf("outputs/sample_histogram_animal_species.pdf", height = 6, width = 10)
p2
dev.off()

# How many genes were found per sample?
p3d <- sample_genotypes %>%
  select(sample_id, host_animal_common, gene)  %>%
  group_by(sample_id, host_animal_common) %>%
  summarize(n = n())
p3d

p3 <- p3d %>%
  ggplot(aes(x = host_animal_common, y = n)) +
  geom_boxplot()+
  labs(fill = "Host Animal" , title = "# of AMR Genes") +
  xlab("Host Animal") + ylab("Genes") + labs(title = "") + 
  theme_classic()
p3

pdf("outputs/unique_genes_animal.pdf")
p3
dev.off()

quit()