## ---------------------------
## Purpose of script: 
##  Tune and model genotypes/phenotype data using Random Forest Models.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(tidymodels)
library(janitor)
library(ranger)
library(themis)
library(gridExtra)
rm(list = ls())


#source("code/helper_functions.R")
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

wt_tune <- readRDS("outputs/iq75_tune_list.RDS")
clsi_tune <- readRDS("outputs/ItoR_tune_list.RDS")

wt_res <- readRDS("outputs/iq75_res_list.RDS")
clsi_res <- readRDS("outputs/ItoR_res_list.RDS")

    
# ROC/PR Curves
collect_rocpr <- function(res_list){
  preds <- lapply(res_list, collect_predictions)
  names <- names(res_list)
  new_names <- list() ; roc <- list(); prc <- list()
  for(i in 1:length(res_list)){
    tryCatch({
    roc[[i]] <- roc_curve(preds[[i]], .pred_0, truth = !!names[[i]])
    prc[[i]] <- pr_curve(preds[[i]], .pred_0, truth = !!names[[i]])
    new_names[[i]] <- names[[i]]
    }, error=function(e){cat("ERROR :\n", i)})
}
  names(prc) <- names(roc) <- new_names
  res <- list(prc, roc)
  names(res) <- c("prc", "roc")
  return(res)

}
# ROC / PR

clsi_rocpr <- collect_rocpr(clsi_res)
wt_rocpr <- collect_rocpr(wt_res)

# PR
prc_df <- bind_rows(clsi_rocpr[[1]], .id = "mic_id") %>%
  mutate(mic_id = gsub("p_", "", mic_id)) %>% 
  left_join(select(reference$classes, test_type_antibiotic_class, mic_id))

prc_groups <- split(prc_df, f = prc_df$test_type_antibiotic_class)
for(i in 1:length(prc_groups)){
  temp <- ggplot(prc_groups[[i]], aes(x = recall, y = precision, color = mic_id)) + 
    geom_path() + 
    facet_wrap(~test_type_antibiotic_class) + ylim(0, 1) +
    labs(title = "Precision-Recall Curve for Antimicrobial Agent", legend = "antimicrobial agent")
  assign(value = temp,x = paste("r", i , sep = ""))
}
length(prc_groups)
prc_plot <- grid.arrange(r1, r2, r3, r4, nrow = 2)

pdf("clsi_pr_curves.pdf", width = 16)
plot(prc_plot)
dev.off()

roc_df <- bind_rows(clsi_rocpr[[2]], .id = "mic_id") %>%
  mutate(mic_id = gsub("p_", "", mic_id)) %>%
  left_join(select(reference$classes, test_type_antibiotic_class, mic_id))

roc_groups <- split(roc_df, f = roc_df$test_type_antibiotic_class)
for(i in 1:length(roc_groups)){
  temp <- ggplot(roc_groups[[i]], aes(x = 1-specificity, y = sensitivity, color = mic_id)) + 
    geom_path() + 
    facet_wrap(~test_type_antibiotic_class) +
    ylim(0, 1) + 
    labs(title = "ROC Curve for Antimicrobial Agent", legend = "antimicrobial agent")
  assign(value = temp,x = paste("u", i , sep = ""))
}
length(roc_groups)
roc_plot <- grid.arrange(u1, u2, u3, u4, nrow = 2)

pdf("clsi_roc_curves.pdf", width = 16)
plot(roc_plot )
dev.off()

# PR
prc_df <- bind_rows(wt_rocpr[[1]], .id = "mic_id") %>%
  mutate(mic_id = gsub("p_", "", mic_id)) %>% 
  left_join(select(reference$classes, test_type_antibiotic_class, mic_id))

prc_groups <- split(prc_df, f = prc_df$test_type_antibiotic_class)
for(i in 1:length(prc_groups)){
  temp <- ggplot(prc_groups[[i]], aes(x = recall, y = precision, color = mic_id)) + 
    geom_path() + 
    facet_wrap(~test_type_antibiotic_class) + ylim(0, 1) +
    labs(title = "Precision-Recall Curve for Antimicrobial Agent", legend = "antimicrobial agent")
  assign(value = temp,x = paste("r", i , sep = ""))
}
length(prc_groups)
prc_plot  <- grid.arrange(r1, r2, r3, r4, r5, r6, r7,
                          r8, r9, r10, r11, r12, r13, nrow = 4)


pdf("wt_pr_curves.pdf", width = 16)
plot(prc_plot )
dev.off()

roc_df <- bind_rows(wt_rocpr[[2]], .id = "mic_id") %>%
  mutate(mic_id = gsub("p_", "", mic_id)) %>%
  left_join(select(reference$classes, test_type_antibiotic_class, mic_id))

roc_groups <- split(roc_df, f = roc_df$test_type_antibiotic_class)
for(i in 1:length(roc_groups)){
  temp <- ggplot(roc_groups[[i]], aes(x = 1-specificity, y = sensitivity, color = mic_id)) + 
    geom_path() + 
    facet_wrap(~test_type_antibiotic_class) +
    ylim(0, 1) + 
    labs(title = "ROC Curve for Antimicrobial Agent", legend = "antimicrobial agent")
  assign(value = temp,x = paste("u", i , sep = ""))
}
length(roc_groups)
roc_plot <- grid.arrange(u1, u2, u3, u4, u5, u6, u7,
                          u8, u9, u10, u11, u12, u13, nrow = 4)

pdf("wt_roc_curves.pdf", width = 16)
plot(roc_plot)
dev.off()

# Metrics
drug_classes <- read_csv("data/reference/drug_classes_new.csv") %>%
  rename(mic_id = antibiotic)

custom_metrics <- function(res_list, drug_classes = drug_classes){
  multi_metric <- metric_set(spec, sens, kap, ppv, npv, accuracy, mcc, f_meas)
  metrics <- list()
  preds <- lapply(res_list, collect_predictions)
  names <- names(res_list)
  new_names <- list()
  for(i in 1:length(res_list)){
    tryCatch({
      metrics[[i]] <- multi_metric(preds[[i]], estimate = .pred_class, truth = !!names[[i]])
      new_names[[i]] <- names[[i]]
    }, error=function(e){cat("ERROR :\n", i)})
  }
  names(metrics) <- new_names

  metrics_df <- bind_rows(metrics, .id = "mic_id") %>% 
    mutate(mic_id = gsub("p_", "", mic_id)) %>%
    left_join(drug_classes) %>%
    unique() %>%
    filter(!is.na(mic_id)) %>%
    filter(mic_id != "NULL")

  return(metrics_df)
}

clsi_metrics <- custom_metrics(clsi_res, drug_classes)
wt_metrics <- custom_metrics(wt_res, drug_classes)

metrics_list <- list(clsi_metrics, wt_metrics)
names(metrics_list) <- c("clsi", "wt")
metrics_df <- bind_rows(metrics_list, .id = "model")

avg <- metrics_df %>%
	group_by(model, .metric) %>%
  summarise(mean = mean(.estimate), sd = sd(.estimate)) %>%
	pivot_wider(names_from = "model", values_from = "mean")

avg_class <- metrics_df %>%
	group_by(model, .metric, class) %>%
	summarise(mean = mean(.estimate), sd = sd(.estimate)) %>%
	pivot_wider(names_from = "model", values_from = "mean")

p1 <- metrics_df %>%
  filter(!(.metric %in% c("kap", "mcc"))) %>%
  mutate(.metric2 = case_when( # convert species names to common name
    grepl("accuracy", .$.metric, ignore.case = TRUE) ~ "Accuracy",
    grepl("f_meas", .$.metric, ignore.case = TRUE) ~ "F1-score",
    grepl("npv", .$.metric, ignore.case = TRUE) ~ "NPV",
    grepl("ppv", .$.metric, ignore.case = TRUE) ~ "PPV",
    grepl("sens", .$.metric, ignore.case = TRUE) ~ "Sensitivity",
    grepl("spec", .$.metric, ignore.case = TRUE) ~ "Specificity")) %>%
  mutate(mic_id = factor(mic_id, levels = unique(drug_classes$mic_id))) %>% 
  filter(!is.nan(.estimate)) %>%
  ggplot(aes(y = mic_id, x = .metric2, fill = .estimate)) + 
  geom_tile(aes(width = .9, height = .9)) + 
  scale_fill_gradient(low = "tomato1", high = "lightgreen") +
  facet_grid(class~model, scales = 'free', space = 'free') +
  geom_text(aes(label=round(.estimate, 2))) +
  theme_bw() + xlab("") + ylab("") + 
  labs(fill = "") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.25, vjust=0.25),
          strip.text.y = element_text(angle = 0))
pdf("outputs/rf_metrics.pdf", height = 10, width = 12)
p1
dev.off()


##



# variable importance
variable_importance <- function(x){
  res <- x$.workflow[[1]]$fit$fit$fit$variable.importance %>%
    stack() %>%
    arrange(desc(values)) 
  return(res)
}

genes <- reference$genotypes$gene 
clean_genes <- janitor::make_clean_names(genes)

genes_df <- cbind(genes, clean_genes) %>%
  as_tibble() %>%
  mutate(clean_genes = paste("g_", clean_genes, sep = ""))

res_list <- list(clsi_res, wt_res)
names(res_list) <- c("clsi", "wt")

res_varimp <- lapply(res_list, function(x){
  res <- lapply(x, variable_importance) %>%
  bind_rows(.id = "antibiotic") %>%
  group_by(antibiotic) %>%
  arrange(desc(values), .by_group = TRUE) %>% 
  slice_head(n = 5) %>%
  filter(values != 0)
  return(res)
  })

varimp_df <- bind_rows(res_varimp, .id = "model") %>%
  rename(clean_genes = "ind") %>%
  left_join(genes_df) %>%
  mutate(genes = trimws(genes)) %>%
  mutate(desc = paste(genes, " [", round(values,2), "]", sep = "")) %>%
  mutate(antibiotic = gsub("p_", "", antibiotic )) %>%
  select(-c(genes, values)) %>%
  group_by(model, antibiotic) %>%
  mutate(desc = paste(desc, collapse = ", ")) %>%
  select(model, antibiotic, desc) %>%
  unique() %>%
  rename(mic_id = "antibiotic") %>%
  left_join(drug_classes) %>%
  select(model, class, mic_id, desc) %>%
  arrange(model, class, mic_id)
colnames(varimp_df) <- c("Model", "Class", "Antibiotic", "Important Variables (n = 5)")
write_csv(varimp_df, "outputs/rf_variable_importance.csv")


