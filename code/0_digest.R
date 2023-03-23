## ---------------------------
## Purpose of script: 
##  Reshape raw experiment data into tidy data and reference 
##  tables for later analysis.
## Author: Henri Chung
## ---------------------------
#module load gsl/2.7.1-uuykddp
# Load required packages
#setwd("~/amr_analysis")
library(tidyverse)
rm(list = ls())

# Set up data folders
dataFolder <- "data/"
dataFiles <- list.files(dataFolder)

require(readxl)

# function to read all of the sheets from an excel document.
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

# Parse through pp# files for sample data.
# Read in the number of individual gene identifier 
# hits per database for each sample.asda
# ===================================
data_files_list <- list()
data_files <- list.files(dataFolder, pattern = ".*ppy.*")
for(j in 1:length(data_files)){
  temp <- read_excel_allsheets(paste(dataFolder, data_files[j], sep = "")) %>%
    .[which(sapply(., nrow) > 0)]
  data_files_list[[j]] <- bind_rows(temp, .id = "sampleID")
}
names(data_files_list) <- data_files

# Separate ppy files based data suffix
#  ast - phenotype data
#  abricate/amrfinder - genotype data 
# ===================================

# classes object is a reference table for the name of the same antibiotic agent across different databases.
classes <- read_csv(paste(dataFolder, "/reference/antibiotic classes.csv", sep = "")) %>%
  mutate(test_type_antibiotic_class = ifelse(test_type_antibiotic_class == "sulfonamide", "sulphonamide", test_type_antibiotic_class))

head(classes)

# Read in information about the S,I,R MIC breakpoints 3
# for different animal species/AB resistances.

# Palanti Urinary breakpoints 
urine_breakpoints <- read_csv(paste(dataFolder, "palantir_mic_sir.csv", sep = "")) %>%
  mutate(host_animal_species = trimws(gsub("\\(organism\\)", "", animal_species_scientific_name))) %>%
  mutate(breakpoint_id = gsub(" MIC", "", test_type_desc)) %>%
  dplyr::select(-c("animal_species_scientific_name", "test_type_desc"))

# Read in Palantir breakpoints
breakpoints <- read_excel(paste(dataFolder, "reference/CLSIbreakpoint_refTablePalantir.xlsx", sep = ""))

# Reshape Breakpoint information
breakpoints_clean <- breakpoints %>%
  filter(bacteria_species_scientific_name_general == "Escherichia coli") %>% # filter to just E.Coli information
  select(c("animal_species_scientific_name", "test_type_desc", "animal_infection_site_uti", contains("test_result"))) %>% # select relevant columns
  mutate(animal_species_scientific_name = trimws(gsub("\\(organism\\)", "", animal_species_scientific_name))) %>% # remove "(organism) from species name column
  unique() %>% # remove duplicate observations
  # change animal species names to alternative spelling
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name == "Genus Felis", "Felis catus", animal_species_scientific_name)) %>% 
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Genus Canis", "Canis familiaris", animal_species_scientific_name)) %>%
  mutate(animal_species_scientific_name = ifelse( animal_species_scientific_name =="Sus scrofa", "Sus domesticus", animal_species_scientific_name)) %>%
  mutate(test_type_desc = gsub(" MIC", "", test_type_desc)) %>% # remove MIC from antibiotic type
  rename(host_animal_species = "animal_species_scientific_name", breakpoint_id = "test_type_desc") %>% # rename
  rbind(urine_breakpoints) %>% #combine CLSI breakpoints with urine specific ones.
  unique()
head(breakpoints_clean)



# Phenotype Data
# reshape and clean phenotype data from "AST" sheets
# ===================================
ast_phenotype_raw <- bind_rows(data_files_list[grepl("ast", names(data_files_list))])

# Check whether samples had a UTI infection diagnosis
urine_diagnosis <- ast_phenotype_raw %>%
  janitor::clean_names() %>%
  select(c("final_diagnosis")) %>%
  mutate(final_diagnosis = tolower(final_diagnosis), animal_infection_site_uti = ifelse(grepl("uti|pyuria|urinary.*infection", final_diagnosis), "UTI/urine/urinary tract", "Not UTI/urine/urinary tract")) %>%
  unique()

write_csv(urine_diagnosis, "data/reference/urine_lookup_post.csv")

# Reshape and clean AST phenotype data into tidy format
ast_phenotype <- ast_phenotype_raw %>%
  janitor::clean_names() %>%  # clean column names
  select(-c("laboratory_name", "unique_specimen_id")) %>%
  mutate(animal_species = tolower(animal_species), # reduce column values to lowercase
         specimen_source_tissue = tolower(specimen_source_tissue), 
         final_diagnosis = tolower(final_diagnosis)) %>%
  mutate("animal_species" = gsub("poultry-domestic ", "", animal_species)) %>% # remove "poultry-domestic" from animal species columns
  mutate(animal_species = factor(animal_species, levels = c("cattle", "swine", "chicken", "turkey", "equine", "dog", "cat", "ducks"))) %>% #convert animals to factor to order
  rename(host_animal_common = "animal_species") %>% # change animal_species to host_animal_comon
  mutate(host_animal_species = case_when( # convert species names to common name
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "Anas platyrhynchos",
    grepl("cat\\>", .$host_animal_common, ignore.case = TRUE) ~ "Felis catus",
    grepl("cattle", .$host_animal_common, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chicken", .$host_animal_common, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turkey", .$host_animal_common, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("swine", .$host_animal_common, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  mutate(host_animal_species = ifelse(host_animal_common == "equine", "Equus caballus", host_animal_species)) %>%
  mutate(animal_infection_site_uti = ifelse(grepl("uti|pyuria|urinary.*infection", final_diagnosis), "UTI/urine/urinary tract", "Not UTI/urine/urinary tract"))

# separate metadata from all data
sample_metadata <- ast_phenotype %>% 
  mutate(sample_id = gsub("\\.1", "", sample_id)) %>%
  select(-contains("mic")) %>% # remove MIC values
  select(sample_id, host_animal_species, host_animal_common, everything())  %>%
  mutate(host_animal_common = as.character(host_animal_common)) %>%
  filter(host_animal_common != "ducks") %>% mutate(host_animal_common = as.character(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) %>%
  unique() 

# For samples that do not have a reference breakpoint value,
# determine a value by checking what quantile breakpoint best classifies known samples
# ===================================asd


# read in ecoff values
ecoff_values <- read_csv("data/ecoli_eucast_ecoff_2023-02-02.csv") %>%
  janitor::clean_names() %>%
  dplyr::select(mic_id, t_ecoff) %>%
  unique() %>%
  mutate(t_ecoff = as.numeric(t_ecoff)) %>%
  filter(!is.na(t_ecoff))


# reshape ast phenotype matrix into long format
phenotypes_clean <- ast_phenotype %>% 
  select(sample_id, host_animal_species, animal_infection_site_uti,  contains("mic")) %>%
  reshape2::melt(id.vars = c("sample_id", "host_animal_species", "animal_infection_site_uti")) %>% # reshape data from wide to long
  mutate(variable = gsub("_mic", "", variable)) %>% # remove _mic suffix 
  rename(mic_id = "variable") %>% # change variable to mic_id
  left_join(ecoff_values, by = "mic_id") %>%
  left_join(select(classes, c("mic_id", "breakpoint_id")), by = "mic_id") %>%
  left_join(breakpoints_clean, by = c("host_animal_species", "breakpoint_id", "animal_infection_site_uti")) %>% # join with breakpoints data
  mutate(num = as.character(str_extract(value, "[0-9/.]+"))) %>%
  separate(num, into = c("a", "b"), sep = "/") %>%
  pivot_longer(cols = c("a", "b"), names_to = "blank", values_to = "num") %>%
  filter(blank != "b") %>%
  select(-c("blank")) %>%
  filter(!is.na(num)) %>%
  mutate(num = as.numeric(num)) %>%
  mutate(equal = str_extract(value, "[<=>]+")) %>%
  mutate(equal = ifelse(is.na(equal), "=", equal)) %>%
  mutate(phenotype = NA)  %>%
  mutate(mic_id = gsub("enrofloxacin_l", "enrofloxacin", mic_id)) 

# separate mic
temp_phenotypes <- phenotypes_clean %>%
  rowwise() %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num == test_result_threshold_mic_resistant_min) & (equal == ">" | equal == ">=" | equal == "=")), "R", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num > test_result_threshold_mic_resistant_min) & (equal != "<=")), "R", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & ((num == test_result_threshold_mic_intermediate_min) | (num == test_result_threshold_mic_intermediate_max)) & (equal == "=")), "I", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num <= test_result_threshold_mic_sensitive_max) & (equal == "<" | equal == "<=" | equal == "=")), "S", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num <= test_result_threshold_mic_intermediate_min | num <= test_result_threshold_mic_resistant_min) & (num > test_result_threshold_mic_sensitive_max) & (equal == "<=")), "NI", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num > test_result_threshold_mic_resistant_min) & (equal == "<=")), "NI", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (!is.na(test_result_threshold_mic_sensitive_max)) & (is.na(test_result_threshold_mic_intermediate_max))), "NI", phenotype))

# label phenotypes assigned with clsi breakpoints
clsi_phenotypes <- temp_phenotypes %>% filter(!is.na(phenotype)) %>%
  mutate(breakpoint = "clsi")

# label phenotypes assigned with ecoff breakpoints
ecoff_phenotypes <- temp_phenotypes %>% filter(is.na(phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num == t_ecoff) & (equal == ">" | equal == ">=")), "R", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num == t_ecoff) & (equal == "<" | equal == "<=" | equal == "=")), "S", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num > t_ecoff) & (equal != "<=" & equal != "<")), "R", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype) & (num < t_ecoff & (equal != ">=" & equal != ">"))), "S", phenotype)) %>%
  mutate(phenotype = ifelse((is.na(phenotype)), "NI", phenotype)) %>%
  mutate(breakpoint = "ecoff") %>%
  left_join(unique(select(sample_metadata, host_animal_common, host_animal_species)), by = "host_animal_species")

# Abs and animal combinations with no ECOFF 
e0_data <- ecoff_phenotypes %>% 
  filter(is.na(t_ecoff)) %>% 
  group_by(mic_id, host_animal_common) %>% 
  summarize(n = n()) %>% 
  group_by(mic_id) %>% 
  summarize(n = sum(n), host_animal_common = paste(sort(host_animal_common), collapse = ", ")) %>% 
  select(mic_id, host_animal_common, n) %>%
  arrange(mic_id) 
write_csv(e0_data, "outputs/ecoff_excluded.csv")

# short-hand table to review assigned phenotypes.
review <- temp_phenotypes %>%
  filter(!is.na(phenotype)) %>%
  group_by(num, equal, phenotype, test_result_threshold_mic_sensitive_max, test_result_threshold_mic_intermediate_min, test_result_threshold_mic_intermediate_max, test_result_threshold_mic_resistant_min, t_ecoff) %>%
  summarize(n = n())
write.csv(review, "outputs/review_phenotypes.csv")

# record distribution of MIC values
mic_values <- temp_phenotypes %>%
  mutate(value = paste0(equal, num)) %>% 
  mutate(equal = factor(equal, levels = c("<=", "=", ">"))) %>%
  select(sample_id, host_animal_species, animal_infection_site_uti, mic_id, num, equal, value)
mic_levels <- mic_values %>% arrange(num, equal) %>% pull(value) %>% unique() 
mic_values <- mutate(mic_values, value = factor(value, levels = mic_levels))
write_csv(mic_values, "outputs/mic_values.csv")

# label phenotypes
phenotypes <- rbind(clsi_phenotypes, ecoff_phenotypes) %>%
  as_tibble() %>%
 filter(!is.na(value)) 


# Separate abricate and amrfinder files based 
# genotype data 
# ===================================

# reshape and clean abricate data
abricate_raw <- bind_rows(data_files_list[grepl("abricate", names(data_files_list))])

abricate_clean <- abricate_raw %>%
  janitor::clean_names() %>% # clean column names
  select(c("sample_id", "strand", "gene", "database", "accession", "product", "resistance", "x_coverage", "x_identity")) %>% # select relevant columns
  rename(coverage = "x_coverage", identity = "x_identity") %>% # rename columns
  mutate(tool = "abricate") %>% # label all observations with abricate
  mutate(sample_id = ifelse(sample_id == "EC-CowMN55108PPY30052", "EC-Cow-MN55108PPY30052", sample_id)) %>% # adjust two errored sample_id
  separate(sample_id, into = c("EC", "host_animal_common", "sample_id"), sep = "-") %>%  # separate sample_id to include animal
  select(-c("EC")) %>% # remove EC
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate(host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "equine",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine",
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "duck"))
head(abricate_clean)


# reshape and clean amrfinder data
amrfinder_raw <-  bind_rows(data_files_list[grepl("amrfinder", names(data_files_list))])

amrfinder_clean <- bind_rows(data_files_list[grepl("amrfinder", names(data_files_list))]) %>% 
  janitor::clean_names() %>%  # clean column names
  select(c("sample_id", "strand", "gene_symbol", "sequence_name", "class", "subclass", "accession_of_closest_sequence", 
           "x_coverage_of_reference_sequence", "x_identity_to_reference_sequence")) %>% # select relevant columns
  rename(coverage = "x_coverage_of_reference_sequence", identity = "x_identity_to_reference_sequence", accession = "accession_of_closest_sequence", gene = "gene_symbol") %>% # rename columns
  mutate(tool = "amrfinder") %>% # label observations with amrfinder
  separate(sample_id, into = c("EC", "host_animal_common", "sample_id"), sep = "-") %>%
  select(-c("EC")) %>%
  mutate(host_animal_common = tolower(host_animal_common)) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "cowmn55108ppy30052", "cow", host_animal_common)) %>%
  mutate(host_animal_common = case_when( #change common animal names to standard format.
    grepl("cat", .$host_animal_common, ignore.case = TRUE) ~ "cat",
    grepl("cow", .$host_animal_common, ignore.case = TRUE) ~ "cattle",
    grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "dog",
    grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "horse",
    grepl("chick", .$host_animal_common, ignore.case = TRUE) ~ "chicken",
    grepl("turk", .$host_animal_common, ignore.case = TRUE) ~ "turkey",
    grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "duck", 
    grepl("pig", .$host_animal_common, ignore.case = TRUE) ~ "swine")) %>%
  mutate(database = "amrfinder") %>%
  rename(product = "sequence_name", resistance = "class")
head(amrfinder_clean)


# combine abricate and amrfinder data
genotypes <- bind_rows(abricate_clean, amrfinder_clean) %>%
  mutate(gene = gsub("_[0-9]", "", gene)) %>%
  select(sample_id, everything()) %>%
  mutate(sample_id = ifelse(sample_id == "N47907PPY10019", "IN47907PPY10019", sample_id)) 

# save gene identifiers to lookup in CARD database
identifiers_lookup <- genotypes %>%
  select(gene, database, accession) %>%
  unique() %>%
  arrange(gene)
head(identifiers_lookup)
write_csv(identifiers_lookup, "data/reference/identifiers_lookup.csv")  

drug_classes <- read.csv("data/reference/drug_classes_new.csv", encoding = "UTF-8") %>%
  rename("resistance_class" = "class", "resistance_drug2" = "antibiotic") %>%
  mutate(resistance_class = ifelse(resistance_class == "beta lactam combo", "beta-lactam", resistance_class)) %>%
  mutate(resistance_class = ifelse(resistance_class == "rifamicin", "rifamycin", resistance_class)) %>%
  mutate(resistance_class = ifelse(resistance_class == "lincomycin", "lincosamide", resistance_class))

# Reference table which labels which genes / gene identifiers are supposed to confer resistance to drug classes / specific drugs.
gene_metadata <- read.csv("data/reference/identifiers_lookup_post_new.csv", encoding = "UTF-8") %>%
  select(-c("database", "accession")) %>%
  mutate(gene_type = ifelse(gene_type == "chromosomal", "chromosome", gene_type)) %>%
  separate_rows(resistance_drug, resistance_class, sep = "[,/]") %>%
  left_join(drug_classes, by = "resistance_class") %>%
  mutate(resistance_drug = ifelse(resistance_drug == "", resistance_drug2, resistance_drug)) %>%
  select(-c("resistance_drug2")) %>%
  mutate(resistance_drug = trimws(resistance_drug), resistance_class = trimws(resistance_class)) %>%
  mutate(gene = ifelse(gene == "tet(M)2", "tet(M)", gene)) %>%
  mutate(combo = paste(resistance_drug, resistance_class)) %>%
  filter(combo != "gentamicin fluoroquinolone") %>%
  select(-combo) 

# Combine relevant data together
sample_genotypes <- genotypes %>%
  mutate(gene_type = ifelse(database == "plasmidfinder", "plasmid", "chromosome")) %>%
  mutate(sample_id = gsub("_2", "", sample_id)) %>%
  select(host_animal_common, sample_id, gene, gene_type) %>%
  mutate(gene = ifelse(gene == "tet(M)2", "tet(M)", gene)) %>%
  left_join(select(gene_metadata, gene, gene_identifier,  gene_type, gene_name_family, resistance_class), by = c("gene", "gene_type")) %>%
  mutate(host_animal_common = ifelse(host_animal_common == "equine", "horse", host_animal_common)) %>%
  mutate(sample_id = ifelse(sample_id == "IA50014PP30306", "IA50014PPY30306", sample_id)) %>%
  unique() %>%
  filter(sample_id %in% sample_metadata$sample_id) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("cattle", "swine", "chicken", "turkey", "horse", "cat", "dog"))) 

# reformat phenotypes into long format, while removing threshold values.
sample_phenotypes <- phenotypes %>%
  mutate(sample_id = gsub("\\.1", "", sample_id)) %>%
  select(sample_id, mic_id, phenotype, breakpoint) %>%
  filter(!is.na(phenotype)) %>%
  filter(sample_id %in% sample_metadata$sample_id) %>%
  unique() 

# Preview data before exporting
# ===================================

# sample genotypes.
sample_genotypes %>% head()
# sample metadata.
sample_metadata %>% head()
# sample phenotypes.
sample_phenotypes %>% head()

# check equal sample information in each category.
a <- unique(sample_genotypes$sample_id) 
b <- unique(sample_phenotypes$sample_id) 
c <- unique(sample_metadata$sample_id) 
d <- c(a,b,c)
table(d) %>% as.numeric() %>% summary()
# breakpoint values
breakpoints_clean %>% head()
# reference table of classes.
classes %>% head()
# table of genes, gene family, and conferred resistance.
gene_metadata %>% head()

# Export data.
sample_data <- list(sample_metadata, sample_genotypes, sample_phenotypes)
names(sample_data) <- c("metadata", "genotypes", "phenotypes")
saveRDS(sample_data, paste(dataFolder, "tidy/samples.RDS", sep = ""))

reference_data <- list(breakpoints_clean, classes, gene_metadata, genotypes)
names(reference_data) <- c("breakpoints", "classes", "gene_metadata", "genotypes")
saveRDS(reference_data, paste(dataFolder, "tidy/reference.RDS", sep = ""))
quit()

