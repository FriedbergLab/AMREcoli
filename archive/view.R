library(tidyverse)
library(glmnet)
library(rstatix)
library(yardstick)

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
model_res <- readRDS("outputs/model_res.RDS")
# rank models by highest average f_meas and accuracy
model_ranks <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = mean(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))

# rank models by highest median f_meas and accuracy
model_median <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = median(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))

# rank models by highest sd for f_meas and accuracy
model_sd <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = sd(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))

# is predicting resistance harder than non-wt?
model_res %>% group_by(model, breakpoint, .metric) %>% summarize(mean = mean(.estimate), sd = sd(.estimate)) %>%
    mutate(value = paste(round(mean,3), " (", round(sd,2), ")", sep = "")) %>%
    select(-c("mean", "sd")) %>%
    pivot_wider(names_from = ".metric", values_from = "value") %>%
    arrange()

# is predicting resistance harder than non-wt?
model_res %>% group_by(model, breakpoint, mic_id, .metric) %>% summarize(mean = mean(.estimate), sd = sd(.estimate)) %>%
    mutate(value = paste(round(mean,3), " (", round(sd,2), ")", sep = "")) %>%
    select(-c("mean", "sd")) %>%
    pivot_wider(names_from = ".metric", values_from = "value") %>%
    arrange(mic_id)

# comparing addition of the count.
model_res %>% filter(.metric == "f_meas") %>%
    unite(value, mic_id:breakpoint) %>%
    select(-c("pr", "roc", ".metric", ".estimator")) %>%
    pivot_wider(names_from = "model", values_from = ".estimate") %>%
    mutate(diff = abs(Count - Base)) %>%
    arrange(desc(diff))

model_res %>% filter(.metric == "f_meas") %>%
    unite(value, mic_id:breakpoint) %>%
    select(-c("pr", "roc", ".metric", ".estimator")) %>%
    pivot_wider(names_from = "model", values_from = ".estimate") %>%
    mutate(diff = abs(`Only-Count` - Count)) %>%
    arrange(desc(diff))

#
model_res %>% group_by(model, breakpoint, .metric) %>% summarize(mean = mean(.estimate), sd = sd(.estimate)) %>%
    mutate(value = paste(round(mean,3), " (", round(sd,2), ")", sep = "")) %>%
    select(-c("mean", "sd")) %>%
    pivot_wider(names_from = ".metric", values_from = "value") %>%
    select(sens)


power_vals <- read_csv("outputs/power_vals.csv")
ss_data <- readRDS("outputs/ss_models.RDS")

animal_metrics <- read_csv("outputs/animal_metrics.csv")

animal_metrics %>% select(-.estimator) %>%  pivot_wider(names_from = ".metric", values_from = ".estimate") %>% 
  select(mic_id, breakpoint, model, host_animal_common, f_meas) %>% 
  pivot_wider(names_from = "host_animal_common", values_from = "f_meas") %>%
  replace(is.na(.), 0) %>% filter(model == "Binary-Count") %>% tail(15)