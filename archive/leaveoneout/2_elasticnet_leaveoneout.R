## ---------------------------
## Purpose of script: 
##  Model elastic net over AMR data
##  Testing multiple models.
## Author: Henri Chung
## Efficienc: 16.5 hours 2 GB
## ---------------------------

# Load required packages and set up folder directory
#setwd("~/amr_analysis")
library(tidyverse)
library(glmnet)
library(rstatix)
library(yardstick)

rm(list = ls())
set.seed(123)
message(Sys.time(), " Testing leave one out models.")
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
host_animals <- c("cat", "dog", "cattle", "swine", "horse", "chicken", "turkey")
host_animal_interactions <- c(paste(host_animals, "binary", sep = ":"), paste(host_animals, "count", sep = ":"))

drug_classes <- read_csv("data/reference/drug_classes_new.csv") %>%
  rename(mic_id = antibiotic) %>%
  mutate(class = ifelse(class == "sulfonamide", "sulphonamide", class),) %>%
  mutate(class = ifelse(class == "folate pathway antagonist", "FPA", class)) %>%
  mutate(class = gsub("beta lactam combo", "beta-lactam", class)) 


K = 500 # number of data rows to use. Low number for testing. 
ssn = 1000 # number of iterations for stability selection.

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

phenotypes_wide <- sample_phenotypes %>%
  left_join(unique(dplyr::select(sample_metadata, c("sample_id", "host_animal_common"))), by = "sample_id") %>%
  group_by(mic_id, breakpoint, phenotype) %>% 
  summarize(n = n()) %>%
  mutate(mic_id = ifelse(mic_id == "trimethoprim_sulfamethoxazole", "TMP/SMX", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "amoxicillin_clavulanic_acid", "co-amoxiclav", mic_id)) %>%
  mutate(mic_id = ifelse(mic_id == "piperacillin_tazobactam", "TZP", mic_id)) %>%
  mutate(breakpoint = ifelse(breakpoint == "clsi", "CLSI", "ECOFF")) %>%
  pivot_wider(names_from = "phenotype", values_from = "n") %>%
  arrange(breakpoint, mic_id) %>%
  replace(is.na(.), 0) %>%
  select(breakpoint, everything())
write_csv(phenotypes_wide, "outputs/phenotypes_wide_leaveoneout.csv")
# Elastic Net 
#################

# function to prepare phenotype data for elastic net equation
prep_data <- function(.x, .y, filter = TRUE){
  # join data with genotypes
  temp <- .x %>% left_join(sample_genotypes, by = c("sample_id", "host_animal_common")) %>%  
      dplyr::select(sample_id, host_animal_common, custom_phenotype, gene_identifier) %>%
      unique() %>%
      mutate(value = 1) %>%
      pivot_wider(names_from = "gene_identifier", values_from = "value") %>%
      replace(is.na(.), 0) %>%
      janitor::clean_names()

  # option to filter variables to only a priori genes
  if(filter == TRUE){
  vars <- gene_metadata %>%
    filter(resistance_drug == .y) %>%
    pull(gene_identifier) %>% unique() %>%
    janitor::make_clean_names()
  temp <- temp %>% select(c("sample_id", "host_animal_common", "custom_phenotype"), any_of(vars))
  }
    
  # reshape data to long format
  res <- temp %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = "host_animal_common", values_from = "value") %>%
    replace(is.na(.), 0) %>%
    unique() %>%
    dplyr::select(-c("sample_id")) %>%
    mutate(custom_phenotype = ifelse(custom_phenotype == "R", 1, 0)) 

  # remove most common animal as baseline
  present_animals <- select(res, any_of(host_animals)) %>% colSums() %>% sort()
  if(length(present_animals) > 1){res <- select(res, -names(present_animals)[1])}

  # add binary and count columns columns
  temp2 <- dplyr::select(as.data.frame(res), -any_of(c("duck", "cat", "cattle", "dog", "horse", "chicken", "turkey", "swine")))
  binary = as.numeric(rowSums(temp2) > 1) == 1
  count = rowSums(temp2)/max(rowSums(temp2))
  res <- as.data.frame(res == 1)
  if(filter == TRUE){res$binary <- binary; res$count <- count}

  # calculate interaction terms
  gene_preds <- colnames(select(res, -any_of(host_animals), -custom_phenotype))
  animal_preds <- colnames(select(res, any_of(host_animals)))
  int_terms <- expand.grid(animal_preds, gene_preds)
  int_mat <- bind_cols(apply(int_terms, 1 , function(.x){int <- list(res[[.x[1]]] * res[[.x[1]]]); names(int) <- paste(.x[1], .x[2], sep = ":"); return(int)}))
  res2 <- cbind(res, int_mat)
  return(res2)
}

# function to model the elasticnet
en <- function(.z){

    # format x and y
    x <- as.matrix(dplyr::select(.z, -custom_phenotype))
    y <- .z$custom_phenotype 

    # calculate phenotype weights
    Y = as_tibble(y)
    fraction_0 <- rep(1 - sum(Y == 0) / nrow(Y), sum(Y == 0))
    fraction_1 <- rep(1 - sum(Y == 1) / nrow(Y), sum(Y == 1))
    # assign that value to a "weights" vector
    weights <- numeric(nrow(Y))
    weights[Y == 0] <- fraction_0
    weights[Y == 1] <- fraction_1

    # use 10 fold cross validation to tune the alpha and lambda parameters
    alpha_seq <- seq(0.1, 0.9, 0.1)
    alpha_tune <- map(alpha_seq, function(.a){
       cv <- cv.glmnet(x, y, family = "binomial", type.measure = "class", alpha = .a, standardize = FALSE, weights = weights)
       data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], alpha = .a)
    }) %>% 
    bind_rows() %>%
    arrange(cvm) %>% 
    slice_head(n = 1)

    # cross validation
    cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class", alpha = alpha_tune$alpha, weight = weights, nfolds = 10)
    
    # select best model
    best_model <- glmnet(x, y, lambda = cvfit$lambda.min, alpha = alpha_tune$alpha, family = "binomial",weight = weights)
    return(best_model)
}

# Evaluate model performance
#################

training_prop = 0.80

# define evaluation metrics
multi_metric <- metric_set(spec, sens, ppv, npv, accuracy, f_meas)

model_evaluate <- function(.d){
  model <- .d %>% 
    # split data into training and testing
    mutate(data_split = purrr::map(data, rsample::initial_split, prop = training_prop, strata = "custom_phenotype")) %>%
    mutate(training = purrr::map(data_split, rsample::training), testing = purrr::map(data_split, rsample::testing)) %>%
    # fit model to training data
    mutate(en = purrr::map(training, en)) %>%
    # calculate predictions
    mutate(predictions = purrr::map2(.x = en, .y = testing, function(.x, .y){
      preds <- predict(.x, type = "response",  newx = as.matrix(.y[,-1]), s = .x$lambda)
      res <- tibble(truth = factor(.y[,1], levels = c(FALSE, TRUE)), pred_class1 = as.numeric(preds), pred_class2 = 1-as.numeric(preds), predicted = factor(ifelse(pred_class1 > 0.5, TRUE, FALSE), levels = c(FALSE, TRUE)))
      })) %>%
    # evaluate predictions
    mutate(pr = purrr::map(predictions, pr_auc, truth = truth, estimate = pred_class1)) %>% # PRC
    mutate(roc = purrr::map(predictions, roc_auc, truth = truth, estimate = pred_class1)) %>% # ROC
    mutate(metrics = purrr::map(predictions, multi_metric, truth = truth, estimate =  predicted)) # metrics 
  return(model)
}


# Fit multiple models
p2_data <- p1_data %>% dplyr::select(mic_id, breakpoint, data)
p2_data <- head(p2_data, min(nrow(p1_data), K))

# Base: apriori genes + animal
message(Sys.time(), " Fitting base model.")
base_model <- p2_data %>% 
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("count", "binary"))})) %>%
  model_evaluate()

# Binary: apriori genes + animal + binary
message(Sys.time(), " Fitting binary model.")
binary_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("count"))})) %>%
  model_evaluate()

# Count: apriori genes + animal + count
message(Sys.time(), " Fitting count model.")
count_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("binary"))})) %>%
  model_evaluate()

# Only-Count: count + animal
set.seed(123)
message(Sys.time(), " Fitting count only model.")
purecount_model <- p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, c("custom_phenotype", "count"))})) %>%
  mutate(data_split = purrr::map(data, rsample::initial_split, prop = training_prop)) %>%
  mutate(training = purrr::map(data_split, rsample::training), testing = purrr::map(data_split, rsample::testing)) %>%
  mutate(en = purrr::map(training, function(.d){glm(custom_phenotype ~ count, data = .d, family = "binomial")})) %>%
  mutate(predictions = purrr::map2(.x = en, .y = testing, function(.x, .y){
      preds <- predict(.x, newdata = select(.y, count), type = "response")
      res <- tibble(truth = factor(.y[,1], levels = c(FALSE, TRUE)), pred_class1 = as.numeric(preds), pred_class2 = 1-as.numeric(preds), predicted = factor(ifelse(pred_class1 > 0.5, TRUE, FALSE), levels = c(FALSE, TRUE)))
      })) %>%
  mutate(pr = purrr::map(predictions, pr_auc, truth, pred_class1)) %>%
  mutate(roc = purrr::map(predictions, roc_auc, truth, pred_class1)) %>%
  mutate(metrics = purrr::map(predictions, multi_metric, truth = truth, estimate =  predicted)) 

# Only-Binary: binary only
message(Sys.time(), " Fitting binary only model.")
purebinary_model <-  p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x,c("custom_phenotype", "binary"), any_of(host_animals), any_of(host_animal_interactions))})) %>%
  mutate(data_split = purrr::map(data, rsample::initial_split, prop = training_prop)) %>%
  mutate(training = purrr::map(data_split, rsample::training), testing = purrr::map(data_split, rsample::testing)) %>%
  mutate(en = purrr::map(training, function(.d){glm(custom_phenotype ~ binary, data = .d, family = "binomial")})) %>%
  mutate(predictions = purrr::map2(.x = en, .y = testing, function(.x, .y){
      preds <- predict(.x, newdata = select(.y, binary), type = "response")
      res <- tibble(truth = factor(.y[,1], levels = c(FALSE, TRUE)), pred_class1 = as.numeric(preds), pred_class2 = 1-as.numeric(preds), predicted = factor(ifelse(pred_class1 > 0.5, TRUE, FALSE), levels = c(FALSE, TRUE)))
      })) %>%
  mutate(pr = purrr::map(predictions, pr_auc, truth, pred_class1)) %>%
  mutate(roc = purrr::map(predictions, roc_auc, truth, pred_class1)) %>%
  mutate(metrics = purrr::map(predictions, multi_metric, truth = truth, estimate =  predicted)) 

# Binary-Count: apriori genes + animal + binary + count
message(Sys.time(), " Fitting binary-count model.")
bc_model <- p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  model_evaluate()

# Full : all genes + animal
message(Sys.time(), " Fitting all model.")
full_model <- p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data, filter = FALSE)) %>%
  model_evaluate()


# Summarize Results
#################

# bind models tibbles together list
model_list <-list(base_model, binary_model, count_model, bc_model, full_model, purecount_model, purebinary_model)
names(model_list) <- c("Base", "Binary", "Count", "Binary-Count", "Full", "Only-Count", "Only-Binary")

saveRDS(model_list, "outputs/model_list_leaveoneout.RDS")
model_list <- readRDS("outputs/model_list_leaveoneout.RDS")

# unnest evaluation metrics, convert to df
model_res <- lapply(model_list, function(.l){
  .l %>%
    select(mic_id, breakpoint, pr, roc, metrics) %>%
    unnest(metrics)
  }) %>% bind_rows(.id = "model")

saveRDS(model_res, "outputs/model_res_leaveoneout.RDS")
model_res <- readRDS("outputs/model_res_leaveoneout.RDS")

# plot model evaluation results
p2_plot <- model_res %>%
  mutate(model = factor(model, levels = c("Base", "Binary", "Count", "Binary-Count", "Full", "Only-Count" ,"Only-Binary"))) %>%
  mutate(.metric2 = case_when( # convert species names to common name
    grepl("accuracy", .$.metric, ignore.case = TRUE) ~ "Accuracy",
    grepl("f_meas", .$.metric, ignore.case = TRUE) ~ "F1-score",
    grepl("npv", .$.metric, ignore.case = TRUE) ~ "NPV",
    grepl("ppv", .$.metric, ignore.case = TRUE) ~ "PPV",
    grepl("sens", .$.metric, ignore.case = TRUE) ~ "Sensitivity",
    grepl("spec", .$.metric, ignore.case = TRUE) ~ "Specificity")) %>%
  mutate(breakpoint = ifelse(breakpoint == "ecoff", "ECOFF", "CLSI")) %>%
  ggplot(aes(x = .metric2, y = .estimate, fill = model)) +  
  geom_boxplot() + 
  ylab("") + xlab("Metric") + labs(fill = "Model") + 
  facet_wrap(~breakpoint, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Set1") + 
  theme_bw()

ggsave(p2_plot, file = "outputs/en_models_leaveoneout.tiff", unit = "cm", width = 16, height =  6, dpi = 300, scale = 2)
ggsave(p2_plot, file = "outputs/en_models_leaveoneout.png", unit = "cm", width = 16, height =  6, dpi = 300, scale = 2)

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


# identify best model
best_model <- model_ranks$model[1]

# table of model results by animal
animal_metrics <- lapply(model_list, function(.x){
   .x %>%
    mutate(animal_preds = purrr::map2(testing, predictions, function(.a, .b){
      temp <- cbind(.a, .b) %>%
        pivot_longer(
          cols =any_of(c("cat", "cattle", "dog", "swine", "horse", "chicken", "turkey")), 
          names_to = "host_animal_common", values_to = "value") %>% 
        filter(value == TRUE) %>% 
        select(-value) %>%
        group_by(host_animal_common) %>%
        group_split()
      names <- lapply(temp, function(.x){unique(.x$host_animal_common)}) %>% unlist()
      metrics <- lapply(temp, multi_metric, truth = truth, estimate = predicted)
      names(metrics) <- names
      return(bind_rows(metrics, .id = "host_animal_common"))
    }))
  }) %>%
  bind_rows(.id = "model") %>%
  select(mic_id, breakpoint, model, animal_preds) %>%
  unnest(animal_preds)

write_csv(animal_metrics, "outputs/animal_metrics_leaveoneout.csv")

# determine model accuracy for individual models
p3_data <- animal_metrics %>% filter(.metric == "accuracy") %>%
  mutate(host_animal_common = case_when( # convert species names to common name
  grepl("duck", .$host_animal_common, ignore.case = TRUE) ~ "Ducks",
  grepl("cat\\>", .$host_animal_common, ignore.case = TRUE) ~ "Cats",
  grepl("cattle", .$host_animal_common, ignore.case = TRUE) ~ "Cattle",
  grepl("dog", .$host_animal_common, ignore.case = TRUE) ~ "Dogs",
  grepl("horse", .$host_animal_common, ignore.case = TRUE) ~ "Horses",
  grepl("chicken", .$host_animal_common, ignore.case = TRUE) ~ "Chickens",
  grepl("turkey", .$host_animal_common, ignore.case = TRUE) ~ "Turkeys",
  grepl("swine", .$host_animal_common, ignore.case = TRUE) ~ "Swine")) 

# Count the best model for each animal and breakpoint
p3a_data <- p3_data  %>%
  group_by(mic_id, breakpoint, host_animal_common) %>%
  filter(.estimate == max(.estimate)) %>%
  mutate(breakpoint = ifelse(breakpoint == "clsi", "CLSI", "ECOFF")) %>%
  select(-c(".estimator", ".estimate", ".metric")) %>%
  group_by(host_animal_common, breakpoint, model) %>% 
  summarize(n = n()) %>%
  pivot_wider(names_from = "model", values_from = "n") 

# list the best model for each animal and breakpoint in wide format, show the range between highest accuracy models.
p3b_data <- p3_data %>% ungroup() %>%
  group_by(mic_id, breakpoint, host_animal_common) %>%
  arrange(mic_id, breakpoint, host_animal_common) %>% 
  mutate(range = abs(max(.estimate)-min(.estimate))) %>%
  filter(.estimate == max(.estimate)) %>%
  group_by(mic_id, breakpoint, host_animal_common, range) %>%
  summarize(model = paste(sort(model), collapse = ", "))

# write to table.
p3a_table <- p3a_data %>%
  select(breakpoint, host_animal_common, everything()) %>%
  arrange(breakpoint, host_animal_common) %>%
  rename("Host Animal" = host_animal_common, "BP" = breakpoint) 
write_csv(p3a_table, "outputs/best_model_by_animal_leaveoneout.csv")

p3b_table <- p3b_data %>%
  select(breakpoint, host_animal_common, everything()) %>%
  arrange(breakpoint, host_animal_common) %>%
  rename("Host Animal" = host_animal_common, "BP" = breakpoint) 
write_csv(p3b_table, "outputs/range_model_by_animal_leaveoneout.csv")

# View models sorted by f_meas
model_res %>% filter(.metric == "f_meas") %>% 
select(model, mic_id, breakpoint, .estimate, .estimate) %>% 
pivot_wider(names_from = "model", values_from = ".estimate") %>%
rename(only_binary = "Only-Binary", only_count = "Only-Count", binary_count = "Binary-Count") %>% 
arrange(Base) 

# Stability selection of best elastic net model
#################

en_selection <- function(.z){

    # format x and y
    x <- as.matrix(dplyr::select(.z, -custom_phenotype))
    y <- .z$custom_phenotype 

    # randomly select 80% of columns from x
    n_col = round(ncol(x) * 0.80)
    ind <- sample(1:ncol(x), n_col)
    new_x <- x[,ind]

    # calculate phenotype weights
    Y = as_tibble(y)
    fraction_0 <- rep(1 - sum(Y == 0) / nrow(Y), sum(Y == 0))
    fraction_1 <- rep(1 - sum(Y == 1) / nrow(Y), sum(Y == 1))
    # assign that value to a "weights" vector
    weights <- numeric(nrow(Y))
    weights[Y == 0] <- fraction_0
    weights[Y == 1] <- fraction_1

    # use 10 fold cross validation to tune the alpha and lambda parameters
    alpha_seq <- seq(0.1, 0.9, 0.1)
    alpha_tune <- map(alpha_seq, function(.a){
       cv <- cv.glmnet(new_x, y, family = "binomial", type.measure = "class", alpha = .a, standardize = FALSE, weights = weights)
       data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], alpha = .a)
    }) %>% 
    bind_rows() %>%
    arrange(cvm) %>% 
    slice_head(n = 1)

    # cross validation
    cvfit <- cv.glmnet(new_x, y, family = "binomial", type.measure = "class", alpha = alpha_tune$alpha, weight = weights, nfolds = 10)
    
    # select best model
    best_model <- glmnet(x, y, lambda = cvfit$lambda.min, alpha = alpha_tune$alpha, family = "binomial",weight = weights)
    coefs <- coef(best_model) %>% as.matrix() %>% as.data.frame() %>% filter(abs(s0) != 0) %>% rownames()
    return(coefs)
}

# number of replicates for stability selection
library(furrr)
library(future)
plan(multicore)
tictoc::tic()
ss_data <- model_list[[best_model]] %>%
  mutate(coefs = furrr::future_map(training, function(.z){stack(table(unlist(replicate(en_selection(.z), n = ssn)))/ssn)}, .options = furrr_options(seed = TRUE )))
tictoc::toc()
saveRDS(ss_data, "outputs/ss_models_leaveoneout.RDS")

#quit()


ss_data_leaveoneout <- readRDS("outputs/ss_models_leaveoneout.RDS")
ss_data <- readRDS("outputs/ss_models.RDS")
# the possible features fed into the model.
possible_predictors_loo <- ss_data_leaveoneout  %>%
  dplyr::select(mic_id, breakpoint, training) %>%
  mutate(possible = map(training, function(.x){
    res <- .x %>%
      colnames() %>%
      as_tibble() 
    })) %>%
  dplyr::select(-training) %>%
  unnest(possible) %>%
  unique() %>%
  filter(value != "custom_phenotype") %>%
  rename(gene = "value") %>%
  mutate(possible_loo = TRUE)

important_predictors_loo <- ss_data_leaveoneout %>%
  mutate(imp = purrr::map(coefs, function(.x){
    res <- .x %>%
      as.matrix() %>% 
      as.data.frame() %>%
      dplyr::filter(values >= 0.66) 
  })) %>%
  dplyr::select(c("mic_id","breakpoint", "imp")) %>%
  unnest(imp) %>%
  dplyr::select(-values) %>%
  unique() %>%
  filter(ind != "(Intercept)") %>%
  rename(gene = "ind") %>%
  mutate(important_loo = TRUE)

possible_predictors <- ss_data  %>%
  dplyr::select(mic_id, breakpoint, training) %>%
  mutate(possible = map(training, function(.x){
    res <- .x %>%
      colnames() %>%
      as_tibble() 
    })) %>%
  dplyr::select(-training) %>%
  unnest(possible) %>%
  unique() %>%
  filter(value != "custom_phenotype") %>%
  rename(gene = "value") %>%
  mutate(possible = TRUE)

important_predictors <- ss_data %>%
  mutate(imp = purrr::map(coefs, function(.x){
    res <- .x %>%
      as.matrix() %>% 
      as.data.frame() %>%
      dplyr::filter(values >= 0.66) 
  })) %>%
  dplyr::select(c("mic_id","breakpoint", "imp")) %>%
  unnest(imp) %>%
  dplyr::select(-values) %>%
  unique() %>%
  filter(ind != "(Intercept)") %>%
  rename(gene = "ind") %>%
  mutate(important = TRUE)

# Leave one out comparison.
power_vals <- read_csv("outputs/power_vals.csv") %>%
  mutate(ef = paste0("ef_", ef)) %>%
  pivot_wider(names_from = "ef", values_from = "power")
power_vals_loo <- read_csv("outputs/power_vals_leaveoneout.csv") %>%
  rename(power_loo = "power", s0_loo = "s0") %>%
  mutate(ef = paste0("ef_", ef)) %>%
  pivot_wider(names_from = "ef", values_from = "power_loo")

loo_comparison <- possible_predictors %>% 
  left_join(important_predictors, by = c("mic_id", "breakpoint", "gene")) %>%
  left_join(possible_predictors_loo, by = c("mic_id", "breakpoint", "gene")) %>%
  left_join(important_predictors_loo, by = c("mic_id", "breakpoint", "gene")) %>%
  replace(is.na(.), FALSE) %>%
  left_join(power_vals, by = c("mic_id", "breakpoint", "gene"))
write_csv(loo_comparison, "outputs/comparison_leaveoneout.csv")

loo_comparison_n <- loo_comparison %>% 
  group_by(mic_id, breakpoint) %>% 
  summarize(possible_n = sum(possible), important_n = sum(important), possible_loo_n = sum(possible_loo), important_loo_n = sum(important_loo))  %>%
  left_join(drug_classes, by = "mic_id") %>%
  arrange(breakpoint, class) %>%
  mutate(change = important_loo_n/important_n) %>% 
  mutate(change_type = case_when(
    change < 0.90 ~ "negative",
    change >= 1.10 ~ "positive",
    change >= 0.90 & change < 1.10 ~ "neutral")
  )
write_csv(loo_comparison_n, "outputs/comparison_leaveoneout_n.csv")

# predictors that were true in first model but no longer true in important, but were possible.
loo_change_neg <- loo_comparison %>% filter(important == TRUE & important_loo == FALSE & possible_loo == TRUE) 
write_csv(loo_change_neg, "outputs/negative_leaveoneout.csv")

loo_change_pos <- loo_comparison %>% filter(important == FALSE & important_loo == TRUE & possible_loo == TRUE) 
write_csv(loo_change_pos, "outputs/positive_leaveoneout.csv")


summary(loo_change_neg$ef_2)
summary(loo_change_pos$ef_2)

# calculate the animal related and proportion of important predictors
animal_predictors_loo <- important_predictors_loo %>% 
  mutate(animal =  str_extract(pattern = paste(host_animals, collapse = "|"), string = gene)) %>%
  mutate(animal = ifelse(is.na(animal), "main", animal)) %>%
  group_by(mic_id, breakpoint) %>% 
  nest() %>%
  mutate(prop = purrr::map(data, function(.x){
    temp <- .x %>% mutate(animal2 = ifelse(animal == "main", "main", "animal")) %>%
    pull(animal2) %>% table()
    return(as_tibble(temp/sum(temp)))
    })) %>%
  unnest(prop) %>%
  pivot_wider(names_from = ".", values_from = "n") %>% 
  arrange(desc(animal))
animal_predictors_loo

# calculate the animal related and proportion of possible predictors
animal_possible_loo <- possible_predictors_loo %>% 
  mutate(animal =  str_extract(pattern = paste(host_animals, collapse = "|"), string = gene)) %>%
  mutate(animal = ifelse(is.na(animal), "main", animal)) %>%
  group_by(mic_id, breakpoint) %>% 
  nest() %>%
  mutate(prop = purrr::map(data, function(.x){
    temp <- .x %>% mutate(animal2 = ifelse(animal == "main", "main", "animal")) %>%
    pull(animal2) %>% table()
    return(as_tibble(temp/sum(temp)))
    })) %>%
  unnest(prop) %>%
  pivot_wider(names_from = ".", values_from = "n") %>% 
  rename(animal_possible = "animal", main_possible = "main", data_possible = "data")

# combine possible and important animal proportions and compare
animal_pred_data_loo <- left_join(animal_possible_loo, animal_predictors_loo, by = c("mic_id", "breakpoint")) %>%
  mutate(n_possible = purrr::map(data_possible, nrow), n = purrr::map(data, nrow)) %>%
  unnest(c(n_possible, n)) %>%
  select(-c("data", "data_possible"))

animal_pred_data_loo %>% filter(is.na(animal)) %>% arrange(animal_possible) 
animal_pred_data_loo %>% filter(!is.na(animal)) %>% arrange(animal)

model_res <- readRDS("outputs/model_res")
model_res_loo <- readRDS("outputs/model_res_loo")

model_diffs <- left_join(model_res, model_res_loo, by = c("model", "mic_id", "breakpoint", ".metric", ".estimator")) %>% 
  filter(.estimate.x != .estimate.y) %>%
  select(.estimate.x, .estimate.y, model, .metric, model, mic_id, breakpoint)

filter(model_diffs, .metric == "f_meas" & model == "Binary-Count")