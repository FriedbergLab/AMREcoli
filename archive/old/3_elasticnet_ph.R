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

gene_metadata <- reference$gene_metadata %>%
  mutate(reference$gene_metadata, combo = paste(gene, resistance_drug)) 

gene_key <- gene_metadata %>%
  select(gene_identifier) %>%
  unique() %>%
  mutate(gene = janitor::make_clean_names(gene_identifier))

drug_classes <- read_csv("data/reference/drug_classes_new.csv") %>%
  rename(mic_id = antibiotic) %>%
  mutate(class = ifelse(class == "sulfonamide", "sulphonamide", class))


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


# calculate the proportion of animals which exhibit a resistance phenotype.
message(Sys.time(), " TESTING ANIMAL INTERACTIONS")
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

# The palette with black:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1_plot <- p1_sig %>%
  mutate(prop_data = map(prop_data, function(.c){
    .c %>% 
      pivot_wider(names_from = "custom_phenotype", values_from = "n") %>%
      mutate(prop = R / (R + S))
    })) %>%
  unnest(prop_data) %>%
  mutate(breakpoint = ifelse(breakpoint == "iq", "WT", "CLSI")) %>%
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
  mutate(test_type_antibiotic_class = gsub("folate pathway antagonist", "FPA", test_type_antibiotic_class)) %>%
  mutate(host_animal_common = factor(host_animal_common, levels = c("Swine", "Cattle", "Chickens", "Turkeys", "Horses", "Dogs", "Cats"))) %>%
  mutate(label = paste(mic_id, sig_labels, " (", breakpoint, ")", sep = "")) %>%
  ggplot(aes(x = label, y = prop*100, shape = host_animal_common, size = 5)) + 
  geom_point() + 
  coord_flip() + 
  facet_grid(test_type_antibiotic_class~., scales = "free_y", space = "free_y") + 
  labs(shape = "Host Animal") + 
  xlab("Antimicrobial") + ylab("Percent Resistant") +
  theme(legend.position="right", 
    strip.text.y.right = element_text(angle = 0), 
    text = element_text(size=20),
    axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  guides(size = "none", shape = guide_legend(override.aes = list(size = 5))) +
  scale_shape_manual(values = c(0, 5, 1, 10, 7, 3, 4)) 
pdf("outputs/animal_comparison.pdf", height = 12, width = 13)
p1_plot
dev.off()


# Elastic Net 

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

  # add binary and count columns columns
  temp2 <- dplyr::select(as.data.frame(res), -any_of(c("duck", "cat", "cattle", "dog", "horse", "chicken", "turkey", "swine")))
  binary = as.numeric(rowSums(temp2) > 1) == 1
  count = rowSums(temp2)/max(rowSums(temp2))
  res <- as.data.frame(res == 1)
  if(filter == TRUE){res$binary <- binary; res$count <- count}
  return(res)
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
    mutate(pr = purrr::map(predictions, pr_auc, truth, pred_class1)) %>% # PRC
    mutate(roc = purrr::map(predictions, roc_auc, truth, pred_class1)) %>% # ROC
    mutate(metrics = purrr::map(predictions, multi_metric, truth = truth, estimate =  predicted)) # metrics 
  return(model)
}


# Test multiple models
p2_data <- p1_data %>% dplyr::select(mic_id, breakpoint, data)

# base: apriori genes + animal
message(Sys.time(), " TESTING BASE MODEL")
base_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("count", "binary"))})) %>%
  model_evaluate()

# binary: apriori genes + animal + binary
message(Sys.time(), " TESTING BINARY MODEL")
binary_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("count"))})) %>%
  model_evaluate()

# count: apriori genes + animal + count
message(Sys.time(), " TESTING COUNT MODEL")
count_model <- p2_data %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  mutate(data = purrr::map(data, function(.x){select(.x, -c("binary"))})) %>%
  model_evaluate()

# purecount : count only
message(Sys.time(), " TESTING PURECOUNT MODEL")
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

# bc: apriori genes + animal + binary + count
message(Sys.time(), " TESTING BC MODEL")
bc_model <- p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data)) %>%
  model_evaluate()

# full : all genes + animal
message(Sys.time(), " TESTING FULL MODEL")
full_model <- p1_data %>%
  dplyr::select(mic_id, breakpoint, data) %>%
  mutate(data = purrr::map2(.x = data, .y = mic_id, prep_data, filter = FALSE)) %>%
  model_evaluate()


# Summarize Results

# bind models
model_list <-list(base_model, binary_model, count_model, bc_model, full_model, purecount_model)
names(model_list) <- c("Base", "Binary", "Count", "Binary-Count", "Full", "Only-Count")
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

pair_rank<- model_res %>% 
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>% 
  group_by(model) %>%
  arrange(desc(f_meas), .by_group=TRUE) %>%
  group_split() %>%
  lapply(function(.x){rownames_to_column(.x, "rank")}) %>%
  bind_rows() %>%
  mutate(name = paste(mic_id, breakpoint)) %>%
  select(model, name, rank) %>%
  mutate(rank = as.numeric(rank)) %>%
  arrange(model, name) %>%
  pivot_wider(names_from = "model", values_from = "rank") %>%
  select(-name)

 model_res %>% 
  filter(model == "Count" | model == "Binary-Count") %>% 
  select(model, mic_id, breakpoint, .metric, .estimate) %>% 
  group_by(model, mic_id, .metric) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  pivot_wider(names_from = "model", values_from = ".estimate_mean") %>% 
  mutate(diff = `Binary-Count`-Count) %>% 
  arrange(desc(diff)) %>% filter(.metric == "spec")

model_res %>%
  group_by(model, .metric, mic_id) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  filter(model == "Base" & .metric == "f_meas") %>%
  arrange(.estimate_mean)

model_res %>%
  group_by(model, .metric, mic_id) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  filter(model == "Only-Count" & .metric == "f_meas") %>%
  arrange(.estimate_mean)

pair_rank_cor <- cor(pair_rank, method = "spearman") 
pair_rank_cor[upper.tri(pair_rank_cor, diag = TRUE)] <- NA
pair_rank_cor %>%
  reshape2::melt() %>%
  filter(!is.na(value)) %>%
  arrange(value)
rstatix::cor_test(pair_rank, rank,  method = "spearman")

# plot model evaluation results
p2_plot <- model_res %>%
  mutate(model = factor(model, levels = c("Base", "Binary", "Count", "Binary-Count", "Full", "Only-Count"))) %>%
  mutate(.metric2 = case_when( # convert species names to common name
    grepl("accuracy", .$.metric, ignore.case = TRUE) ~ "Accuracy",
    grepl("f_meas", .$.metric, ignore.case = TRUE) ~ "F1-score",
    grepl("npv", .$.metric, ignore.case = TRUE) ~ "NPV",
    grepl("ppv", .$.metric, ignore.case = TRUE) ~ "PPV",
    grepl("sens", .$.metric, ignore.case = TRUE) ~ "Sensitivity",
    grepl("spec", .$.metric, ignore.case = TRUE) ~ "Specificity")) %>%
  mutate(breakpoint = ifelse(breakpoint == "iq", "WT", "CLSI")) %>%
  ggplot(aes(x = .metric2, y = .estimate, fill = model)) +  
  geom_boxplot() + 
  ylab("") + xlab("Metric") + labs(fill = "Model") + 
  facet_wrap(~breakpoint, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Set1") + theme_bw()

model_res %>% group_by(model, .metric) %>% summarize(mean = mean(.estimate, na.rm = TRUE)) %>% pivot_wider(names_from = ".metric", values_from = "mean")


pdf("outputs/en_models.pdf", width = 11, height = 5)
p2_plot
dev.off()
png(file="outputs/en_models.png", width=1100, height=500)
p2_plot
dev.off()

# what was the best model by F1 score and then accuracy?
model_ranks <- model_res %>%
  group_by(model, .metric) %>%
  summarize(.estimate = median(.estimate, na.rm = TRUE)) %>%
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>%
  arrange(desc(f_meas), desc(accuracy))

best_model <- model_ranks$model[1]
message(Sys.time(), " ", best_model, " is the best model")
# Stability selection of best elastic net model
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
message(Sys.time(), " TESTING STABILITY SELECTION MODEL")
ssn = 1000
ss_data <- model_list[[best_model]] %>%
  mutate(coefs = purrr::map(training, function(.z){stack(table(unlist(replicate(en_selection(.z), n = ssn)))/ssn)}))
saveRDS(ss_data, "outputs/ss_models.RDS")

message(Sys.time(), " COMPLETE")
quit()






power_vals2 <- power_vals %>% 
  filter(power >= 0.80) %>%
  left_join(drug_classes, by = "mic_id") %>%
  group_by(class) %>%
  mutate(l = length(unique(mic_id))) %>%
  group_by(class, gene) %>%
  mutate(n = n()) %>%
  filter(l == n) %>% 
  ungroup() %>%
  select(mic_id, breakpoint, gene, power, s0) %>%
  unique()

# Examine model coefficients 
message(Sys.time(), " TESTING EXTRACTING SIGNIFICANT MODEL FEATURES")
model_coefs <- ss_data %>%
  mutate(imp = purrr::map(coefs, function(.x){
    res <- .x %>%
      as.matrix() %>% 
      as.data.frame() %>%
      dplyr::filter(values > 0.66) 
  })) %>%
  select(c("mic_id","breakpoint", "imp")) %>%
  unnest(imp) %>%
  left_join(drug_classes, by = "mic_id") %>%
  select(-values) %>%
  unique() %>%
  filter(ind != "(Intercept)") %>%
  rename(gene = "ind") %>%
  left_join(power_vals, by = c("mic_id", "breakpoint", "gene")) %>%
  filter(power >= 0.80) %>%
  select(-power)



# determine number of common features shared between models of the same class
common_features <- model_coefs %>%  
  left_join(gene_key, by = "gene") %>%
  mutate(gene_identifier = ifelse(is.na(gene_identifier), gene, gene_identifier)) %>%
  select(-gene) %>%
  group_by(mic_id) %>% 
  select(-breakpoint) %>%
  unique() 

# determine number of unique features for each model
individual_features <- common_features %>%
  mutate(s0 = ifelse(s0 == 0, FALSE, TRUE)) %>%
  group_by(mic_id, class, s0) %>%
  nest(data = gene_identifier)
write_csv(individual_features,  "outputs/individual_features.csv")

# determine number of features in total for each class.
class_features <- common_features %>%
  unnest(features) %>%
  group_by(class) %>% 
  nest(features = c(mic_id, values, ind)) %>%
  mutate(features = purrr::map(features, function(.x){
    l <- length(unique(.x$mic_id))

    res <- .x %>% group_by(ind) %>%
      summarize(n = n()) %>%
      mutate(prop = n / l) %>%
      filter(prop == 1) %>%
      pull(ind) %>%
      unique() %>%
      as.character() %>%
      sort() %>%
      paste(collapse = ", ")
    return(res)
    })) %>%
  unnest(features) %>%
  select(class, features) %>%
  arrange(class) 
write_csv(class_features,  "outputs/class_features.csv")

# calculate number of total features for each antibiotic
total_features <- common_features %>%
  unnest(features) %>% 
  unique() %>%
  group_by(mic_id, class) %>%
  summarize(total_features = n())

# compare power of features
average_power <- model_coefs %>%
  right_join(power_vals2, by = c("mic_id", "breakpoint", "gene")) %>%
  group_by(mic_id, breakpoint) %>% 
  summarize(average_power = mean(power, na.rm = TRUE)) %>%
  group_by(mic_id) %>% 
  summarize(average_power = mean(average_power, na.rm = TRUE)) 

possible_features_n <- model_coefs %>% group_by(mic_id) %>% summarize(possible_features = length(unique(gene)))

unique_features <- common_features %>%
  unnest(features) %>%
  group_by(class) %>% 
  nest(features = c(mic_id, values, ind)) %>%
  mutate(features = purrr::map(features, function(.x){
    l <- length(unique(.x$mic_id))

    res <- .x %>% group_by(ind) %>%
      mutate(n = n()) %>%
      mutate(prop = n / l) %>%
      filter(prop != 1) %>%
      select(-c("values", "n", "prop"))
    return(res)
    })) %>%
  unnest(features) %>% 
  unique() %>%
  group_by(class, mic_id) %>%
  summarize(unique_features = n()) %>%
  left_join(total_features, by = c("mic_id", "class")) %>%
  left_join(average_power, by = "mic_id") %>%
  left_join(possible_features_n, by = c("mic_id")) %>% 
  select(class, mic_id, everything()) %>%
  arrange(class, mic_id,) %>%
  mutate(percent_unique = signif(unique_features/total_features, 2)) %>%
  mutate(average_power = signif(average_power, 2)) %>%
  dplyr::select(class, mic_id, unique_features, total_features, possible_features, percent_unique) %>%
  replace(is.na(.), 0)
write_csv(unique_features,  "outputs/unique_features.csv", na = "")


 model_res %>% 
  filter(model == "Count" | model == "Binary-Count") %>% 
  select(model, mic_id, breakpoint, .metric, .estimate) %>% 
  group_by(model, mic_id, .metric) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  pivot_wider(names_from = "model", values_from = ".estimate_mean") %>% 
  mutate(diff = `Binary-Count`-Count) %>% 
  arrange(desc(diff)) %>% filter(.metric == "spec")

model_res %>%
  group_by(model, .metric, mic_id) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  filter(model == "Base" & .metric == "f_meas") %>%
  arrange(.estimate_mean)

model_res %>%
  group_by(model, .metric, mic_id) %>%
  summarize(.estimate_mean = mean(.estimate)) %>%
  filter(model == "Only-Count" & .metric == "f_meas") %>%
  arrange(.estimate_mean)

pair_rank_cor <- cor(pair_rank, method = "spearman") 
pair_rank_cor[upper.tri(pair_rank_cor, diag = TRUE)] <- NA
pair_rank_cor %>%
  reshape2::melt() %>%
  filter(!is.na(value)) %>%
  arrange(value)
rstatix::cor_test(pair_rank, rank,  method = "spearman")

pair_rank<- model_res %>% 
  pivot_wider(names_from = ".metric", values_from = ".estimate") %>% 
  group_by(model) %>%
  arrange(desc(f_meas), .by_group=TRUE) %>%
  group_split() %>%
  lapply(function(.x){rownames_to_column(.x, "rank")}) %>%
  bind_rows() %>%
  mutate(name = paste(mic_id, breakpoint)) %>%
  select(model, name, rank) %>%
  mutate(rank = as.numeric(rank)) %>%
  arrange(model, name) %>%
  pivot_wider(names_from = "model", values_from = "rank") %>%
  select(-name)

# how many predictors were selected just because they have power?

# models with no predictor power >= 0.80
power_vals %>% group_by(mic_id, breakpoint) %>% filter(sum(power >= 0.80) == 0) %>%
  summarize(power = mean(power)) %>% arrange(mic_id) 

# models with sum predictor's power >= 0.80
power_vals %>% group_by(mic_id, breakpoint) %>% filter(sum(power >= 0.80) != 0) %>%
  summarize(power = mean(power)) %>% arrange(mic_id) %>% arrange(power) %>% tail()
