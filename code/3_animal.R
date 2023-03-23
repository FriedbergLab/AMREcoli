## ---------------------------
## Purpose of script: 
##  Plot and test differences in resistant prevalence.
## Author: Henri Chung
## Efficienc:
## ---------------------------

# Load required packages and set up folder directory
#setwd("~/amr_analysis")
library(tidyverse)
library(glmnet)
library(rstatix)
library(yardstick)
library(grid)
library(gridExtra)

rm(list = ls())
set.seed(123)

message(Sys.time(), " Loading model data.")

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

# vector of animals
host_animals <- c("cat", "cattle", "dog", "swine", "horse", "chicken", "turkey")

# table for antibiotic class
drug_classes <- read_csv("data/reference/drug_classes_new.csv") %>%
  rename(mic_id = antibiotic) %>%
  mutate(class = ifelse(class == "sulfonamide", "sulphonamide", class),) %>%
  mutate(class = ifelse(class == "folate pathway antagonist", "FPA", class)) %>%
  mutate(class = gsub("beta lactam combo", "beta-lactam", class)) %>% 
  mutate(class = gsub("aminoglycoside", "  aminoglycoside  ", class)) %>%
  mutate(class = gsub("fluoroquinolone", "  fluoroquinolone  ", class))

# Is there a difference in phenotype proportions between different animals?

# Filter data 
# test antimicrobial, breakpoint, and phenotype combinations with more than 5 samples
# test antimicrobial and breakpoint combinations that show both phenotypes
# label intermediate phenotypes as resistant

p1_data <- sample_phenotypes %>%
  left_join(unique(dplyr::select(sample_metadata, c("sample_id", "host_animal_common"))), by = "sample_id") %>%
  filter(phenotype != "NI") %>% 
  mutate(custom_phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  dplyr::select(-phenotype) %>%
  group_by(mic_id, breakpoint, custom_phenotype) %>%
  filter(n() > 5) %>%
  ungroup() %>%
  unique() %>%
  mutate(combo = paste(mic_id, breakpoint))
  
keep <- p1_data %>% group_by(mic_id, breakpoint, combo) %>% 
  summarize(k = length(unique(custom_phenotype))) %>%
  filter(k > 1)

p1_data <- p1_data %>%
  filter(combo %in% keep$combo) %>% select(-combo) %>%
  nest(data = c("sample_id", "host_animal_common", "custom_phenotype"))

mic_values <- read_csv("outputs/mic_values.csv")
mic_levels <- mic_values %>% arrange(num, equal) %>% pull(value) %>% unique() 
mic_values <- mutate(mic_values, value = factor(value, levels = mic_levels))


# Fisher test for animal resistance differences

# calculate the proportion of animals which exhibit a resistance phenotype.
p1_data <- p1_data %>%
  mutate(prop_data = map(data, function(.x){
      .y <- .x %>% 
        group_by(host_animal_common, custom_phenotype) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(host_animal_common = as.character(host_animal_common), custom_phenotype = factor(custom_phenotype, levels = c("R", "S"))) %>%
        ungroup() %>%
        tidyr::complete(host_animal_common, custom_phenotype,fill =  list(n = 0)) %>%
        arrange(host_animal_common, custom_phenotype, n)}))

# calculate a fisher's exact test for difference in proportions. Adjust P-value using BH
p1_data <- p1_data %>%
  mutate(prop_test = map(prop_data, function(.z){
      if(nrow(.z) < 3){return(NULL)}else{
        r = filter(.z, custom_phenotype == "R") %>% pull(n)
        s = filter(.z, custom_phenotype == "S") %>% pull(n)
        m <- matrix(c(r, s), ncol = 2)
        rownames(m) = unique(.z$host_animal_common)
        res <- fisher.test(m, simulate.p.value=TRUE)
        return(res)
      } 
    })) %>%
  mutate(pval = map(prop_test, function(.a){.a$p.value})) %>% 
  unnest(pval) %>%
  mutate(pval.adj = p.adjust(pval, method = "holm")) 

# calculate pairwise comparisons between interactions
p1_data <- p1_data %>%  
  mutate(pwc = map(prop_data, function(.b){
    res <- .b %>% 
      pivot_wider(names_from = "custom_phenotype", values_from = "n") %>% column_to_rownames("host_animal_common") %>% 
      as.matrix() %>%
      pairwise_fisher_test(p.adjust.method = "holm")
    return(res)
    })) 


# Plot significant differences

# filter to significant difference between animal groups
p1_sig <- p1_data %>% filter(pval.adj <= 0.05) %>%
  mutate(sig_labels = "") %>%
  mutate(sig_labels = ifelse(pval.adj <= 0.05, "*", sig_labels)) %>%
  mutate(sig_labels = ifelse(pval.adj <= 0.01, "**", sig_labels)) %>%
  mutate(sig_labels = ifelse(pval.adj <= 0.001, "***", sig_labels)) 

# Color palette for plotting.
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Plot significant differences in animal resistance.
p1_plot <- p1_sig %>%
  mutate(prop_data = map(prop_data, function(.c){
    .c %>% 
      pivot_wider(names_from = "custom_phenotype", values_from = "n") %>%
      mutate(prop = R / (R + S))
    })) %>%
  unnest(prop_data) %>%
  mutate(breakpoint = ifelse(breakpoint == "ecoff", "ECOFF", "CLSI")) %>%
  left_join(unique(select(reference$classes, test_type_antibiotic_class, mic_id)), by = "mic_id") %>%
  ungroup() %>%
  mutate(host_animal_common = case_when( # convert species names to common name
  grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "Ducks",
  grepl("cat\\>", .$host_animal_common, ignore.case = TRUE) ~ "Cats",
  grepl("cattle", .$host_animal_common, ignore.case = TRUE) ~ "Cattle",
  grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Dogs",
  grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "Horses",
  grepl("chicken", .$host_animal_common, ignore.case = TRUE) ~ "Chickens",
  grepl("turkey", .$host_animal_common, ignore.case = TRUE) ~ "Turkeys",
  grepl("swine", .$host_animal_common, ignore.case = TRUE) ~ "Swine")) %>%
  mutate(mic_id = gsub("_sulfamethoxazole", "-sulphamethoxazole", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim-sulphamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>%
  mutate(test_type_antibiotic_class = gsub("folate pathway antagonist", "FPA", test_type_antibiotic_class)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("Swine", "Cattle", "Chickens", "Turkeys", "Horses", "Dogs", "Cats"))) %>%
  mutate(label = ifelse(breakpoint == "ECOFF", paste(mic_id, sig_labels, " (Non-WT)", sep = ""), paste(mic_id, sig_labels, " (R)", sep = ""))) %>%
  ggplot(aes(x = label, y = prop*100, shape = host_animal_common, size = 5)) + 
  geom_point() + 
  coord_flip() + 
  facet_grid(test_type_antibiotic_class~., scales = "free_y", space = "free_y",  labeller = labeller(group = label_wrap_gen(width = 10, multi_line = TRUE))) +
  labs(shape = "Host Animal ", title = "B") + 
  xlab("") + ylab("\nPercent Resistant") +
  theme_bw() +
  theme(legend.position="right", 
    strip.text.y.right = element_text(angle = 0, margin = margin(r = 20, l = 20)), 
    text = element_text(size=20), 
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 0)))  +
  guides(size = "none", shape = guide_legend(override.aes = list(size = 5))) +
  scale_shape_manual(values = c(0, 5, 1, 10, 7, 3, 4)) 

p1b_plot <- p1_plot + labs(title = "")
ggsave(p1b_plot, file = "outputs/Fig1.tiff", height = 6, width = 8, units = "cm", scale = 4, dpi = 300, compression = "lzw")
ggsave(p1b_plot, file = "outputs/animal_comparison.png", height = 6, width = 8, units = "cm", scale = 4, dpi = 300)


# Plot the proportion of samples that were resistant to a given antibiotic.
p2 <- sample_phenotypes %>%
  filter(phenotype != "NI") %>% 
  mutate(phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  group_by(mic_id, breakpoint) %>%
  filter(n() > 5) %>%
  group_by(mic_id, breakpoint, phenotype) %>%
  summarize(n = n()) %>%
  group_by(mic_id, breakpoint) %>%
  mutate(prop = n / sum(n)) 

p2_plot <- p2 %>%
  select(-n) %>%
  pivot_wider(names_from = "phenotype", values_from = "prop") %>% 
  replace(is.na(.), 0) %>% 
  rename(prop = "R") %>% 
  select(-S) %>%
  left_join(drug_classes, by = "mic_id") %>%
  mutate(mic_id = gsub("_", "-", mic_id)) %>%
  mutate(mic_id = gsub("-acid", " acid", mic_id)) %>%
  mutate(mic_id = gsub("-sulfamethoxazole", "-sulphamethoxazole", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim-sulphamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin-clavulanic acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin-tazobactam", "TZP", mic_id)) %>%
  mutate(breakpoint = ifelse(breakpoint == "ecoff", "Non-WT", "Resistant")) %>% 
  ggplot(aes(x = fct_reorder(mic_id, prop), y = prop*100, shape = breakpoint, size = 5)) + 
  geom_point() +
  coord_flip() +
  facet_grid(class~., scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 100, multi_line = TRUE)) +
  labs(title = "A", shape = "Phenotype") + xlab("Antimicrobial") + ylab("\nPercent") +
  scale_shape_manual(values = c(0, 16)) +
  guides(size = "none", shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() + 
  theme(
    legend.position="right", 
    strip.text.y.right = element_text(angle = 0, margin = margin(r = 20, l = 20)), 
    text = element_text(size=20),
    plot.title = element_text(face = "bold"),
   axis.title.y = element_text(angle = 0, vjust = 0.5)) 

p2b_plot <- p2_plot + labs(title = "")   
ggsave(p2b_plot, file = "outputs/prevalence.tiff", height = 6, width = 8, units = "cm", scale = 4, dpi = 300, compression = "lzw")
ggsave(p2b_plot, file = "outputs/prevalence.png", height = 6, width = 8, units = "cm", scale = 4, dpi = 300)

p3 <- grid.arrange(p2_plot, p1_plot, nrow = 1, vp=viewport(width=1, height=1))
ggsave(p3, file = "outputs/prevalence_animal.tiff", height = 6, width = 16, units="cm", scale=4, compression = "lzw")
ggsave(p3, file = "outputs/prevalence_animal.png", height = 6, width = 16, units="cm", scale=4)