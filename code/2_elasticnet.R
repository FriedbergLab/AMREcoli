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
seed = 123

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
host_animal_interactions <- c(paste(host_animals, "binary", sep = ":"), paste(host_animals, "count", sep = ":"))

K = 500 # number of data rows to use. Low number for testing, K > n defaults to n. 
ssn = 1000 # number of iterations for stability selection.

# Is there a difference in phenotype proportions between different animals?

# Filter data 
# test antimicrobial, breakpoint, and phenotype combinations with more than 5 samples
# test antimicrobial and breakpoint combinations that show both phenotypes
# label intermediate phenotypes as resistant

p0_data <- sample_phenotypes %>%
  left_join(unique(dplyr::select(sample_metadata, c("sample_id", "host_animal_common"))), by = "sample_id") %>%
  filter(phenotype != "NI") %>% 
  mutate(custom_phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  dplyr::select(-phenotype) %>%
  group_by(mic_id, breakpoint, custom_phenotype) %>%
  filter(n() > 5) %>%
  ungroup() %>%
  unique() %>%
  mutate(combo = paste(mic_id, breakpoint))

nrow(sample_phenotypes) # 18437
# == NI, 7200
# n() <= 5, 17
# k <= 1, 1048

keep <- p0_data %>% group_by(mic_id, breakpoint, combo) %>% 
  summarize(k = length(unique(custom_phenotype))) %>%
  filter(k > 1)

p1_data <- p0_data %>%
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
write_csv(phenotypes_wide, "outputs/phenotypes_wide.csv")
# Elastic Net 
#################

# function to prepare phenotype data for elastic net equation
prep_data <- function(.x, .y, filter = TRUE){
  # join data with genotypes
  temp <- .x %>% left_join(sample_genotypes, by = c("sample_id", "host_animal_common"), multiple = "all") %>%  
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

    # manually decide folds for reproducibility
    set.seed(seed)
    fold_ids <- sample(rep(seq(10), length.out = nrow(x)))

    # use 10 fold cross validation to tune the alpha and lambda parameters
    alpha_seq <- seq(0.1, 0.9, 0.1)
    alpha_tune <- map(alpha_seq, function(.a){
       cv <- cv.glmnet(x, y, family = "binomial", type.measure = "class", alpha = .a, standardize = FALSE, weights = weights, foldid = fold_ids)
       data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], alpha = .a)
    }) %>% 
      bind_rows() %>%
      arrange(cvm) %>% 
      slice_head(n = 1)
      
    # cross validation
    cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class", alpha = alpha_tune$alpha, weight = weights, foldid = fold_ids)
    
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
  set.seed(seed)
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

message(Sys.time(), " Fitting binary no animal model.")
binary_na_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("count"))})) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -contains(host_animals))})) %>%
  model_evaluate()


# Count: apriori genes + animal + count
message(Sys.time(), " Fitting count model.")
count_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("binary"))})) %>%
  model_evaluate()


# Only-Count: count + animal
message(Sys.time(), " Fitting count only model.")
purecount_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, c("custom_phenotype"), contains("count"), any_of(host_animals))})) %>%
  model_evaluate()

message(Sys.time(), " Fitting count only model.")
purecount_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, c("custom_phenotype", "count"))})) %>%
  model_evaluate()

# Only-Binary: binary only
message(Sys.time(), " Fitting binary only model.")
purebinary_model <-  p2_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, c("custom_phenotype"), contains("binary"), any_of(host_animals))})) %>%
  model_evaluate()

# Binary-Count: apriori genes + animal + binary + count
message(Sys.time(), " Fitting binary-count model.")
bc_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  model_evaluate()

# Binary-Count (no animal): apriori genes + binary + count
message(Sys.time(), " Fitting binary-count no animal model.")
bc_na_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -contains(host_animals))})) %>%
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
model_list <-list(base_model, binary_model, count_model,  bc_model, bc_na_model, full_model, purecount_model, purebinary_model)
names(model_list) <- c("Gene", "Gene+Binary", "Gene+Count",  "Gene+Binary+Count", "BC-noanimal",  "Full", "Count", "Binary")

saveRDS(model_list, "outputs/model_list.RDS")
model_list <- readRDS("outputs/model_list.RDS")

# unnest evaluation metrics, convert to df
model_res <- lapply(model_list, function(.l){
  .l %>%
    select(mic_id, breakpoint, pr, roc, metrics) %>%
    unnest(metrics)
  }) %>% bind_rows(.id = "model")

saveRDS(model_res, "outputs/model_res.RDS")
model_res <- readRDS("outputs/model_res.RDS")

# plot model evaluation results
p2_plot <- model_res %>%
  filter(model != "BC-noanimal" & model != "Binary-na" & model != "Count-na") %>%
  mutate(model = ifelse(model == "Base", "Gene", model)) %>%
  mutate(model = factor(model, levels = c("Gene", "Gene+Binary", "Gene+Count",  "Gene+Binary+Count", "Full", "Count", "Binary"))) %>%
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

ggsave(p2_plot, file = "outputs/Fig2.tiff", unit = "cm", width = 16, height =  6, dpi = 300, scale = 2)
ggsave(p2_plot, file = "outputs/en_models.png", unit = "cm", width = 16, height =  6, dpi = 300, scale = 2)

# rank models by highest average f_meas and accuracy
model_ranks <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = mean(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))
model_ranks
# rank models by highest median f_meas and accuracy
model_median <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = median(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))
model_median
# rank models by highest sd for f_meas and accuracy
model_sd <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = sd(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(accuracy), desc(f_meas))
model_sd

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

write_csv(animal_metrics, "outputs/animal_metrics.csv")

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
write_csv(p3a_table, "outputs/best_model_by_animal.csv")

p3b_table <- p3b_data %>%
  select(breakpoint, host_animal_common, everything()) %>%
  arrange(breakpoint, host_animal_common) %>%
  rename("Host Animal" = host_animal_common, "BP" = breakpoint) 
write_csv(p3b_table, "outputs/range_model_by_animal.csv")

# View models sorted by f_meas
model_res %>% filter(.metric == "f_meas") %>% 
select(model, mic_id, breakpoint, .estimate, .estimate) %>% 
pivot_wider(names_from = "model", values_from = ".estimate") %>%
arrange(Gene) 

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
plan(multicore)
library(furrr)
tictoc::tic()
ss_data <- model_list[[best_model]] %>%
  mutate(coefs = furrr::future_map(training, function(.z){stack(table(unlist(replicate(en_selection(.z), n = ssn)))/ssn)}, .options = furrr_options(seed = TRUE )))
tictoc::toc()
saveRDS(ss_data, "outputs/ss_models.RDS")

# write samples excluded from analysis

e1_data <- sample_phenotypes %>%
  left_join(unique(dplyr::select(sample_metadata, c("sample_id", "host_animal_common"))), by = "sample_id") %>%
  filter(phenotype == "NI") 

e2_data <- sample_phenotypes %>%
  left_join(unique(dplyr::select(sample_metadata, c("sample_id", "host_animal_common"))), by = "sample_id") %>%
  filter(phenotype != "NI") %>% 
  mutate(custom_phenotype = ifelse(phenotype == "I", "R", phenotype)) %>%
  dplyr::select(-phenotype) %>%
  group_by(mic_id, breakpoint, custom_phenotype) %>%
  filter(n() <= 5) %>%
  ungroup() %>%
  unique() %>%
  mutate(combo = paste(mic_id, breakpoint))
e2_data %>% group_by(mic_id, breakpoint) %>% summarize(n = n())
nrow(sample_phenotypes) # 18437
# == NI, 7200
# n() <= 5, 17
# k <= 1, 1048

excluded <- p0_data %>% group_by(mic_id, breakpoint, combo) %>% 
  summarize(k = length(unique(custom_phenotype))) %>%
  filter(k <= 1)

e3_data <- p0_data %>%
  filter(combo %in% excluded$combo) %>% select(-combo) 
e3_data %>% group_by(mic_id, breakpoint) %>% summarize(n = n())
excluded_ecoffs <- read_csv("outputs/ecoff_excluded.csv")

#quit()


