library(tidyverse)
library(glmnet)
library(parallel)
rm(list = ls())

nsim = 1000
ncores = 36


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

# elastic net model function
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


# take training data from model output
model_data <- readRDS("outputs/ss_models.RDS") %>%  select(mic_id, breakpoint, training, en)

.x <- model_data$training[[1]]

function(.x){
  resp <- .x$custom_phenotype
  preds <- select(.x, -custom_phenotype)
  apply(preds, 2, function(.y){
      no_y <- resp[!.y]
      base <- sum(no_y)/length(no_y)
      logit <- glm(resp ~ .y, family = "binomial")
    })
}
# calculate power with simulation
custom_powersim <- function(.x, n = 10, ef = 0.7){

  # get predictors
  preds <- .x[,-1]

  # assign effect size to positive predictors
  effect_size <-  rep(ef, times = ncol(preds))
  sim_preds <- preds * effect_size

  # probability of resistance 
  p <- 1 / (1 + exp(-rowSums(sim_preds)))

  # determine the significant predictors using the model for 
  # the number of predictors and the given effect size
  # return NULL if their is less than 6 samples with resistance (indeterminate model)
  sig_list <- list()
  seq <- 1:n
  res <- mclapply(seq, function(.x){
    phenotypes <- replicate(expr = rbinom(n = nrow(preds), size = 1, prob = p), n = 1000)
    phenotypes_dist <- apply(phenotypes, 2,  function(.y){res <- (min(table(.y)) > 5) & (length(unique(.y)) > 1)})
    if(sum(phenotypes_dist) < 1){return(NULL)}else{preds$custom_phenotype <- phenotypes[,which(phenotypes_dist == TRUE)[1]]}
    temp_model <- en(preds)
    coefs <- t(as.matrix(coef(temp_model)))
    sig <- colnames(coefs)[abs(coefs) > 0]
    return(sig)
  }, mc.cores = ncores) %>%
  unlist() %>% table() 

  # return the percent of replicates for which a predictor was significant.
  return(res / n)
}

message(Sys.time(), " nsim = ", nsim, ",  ncores = ", ncores)
# calculate power for 2x effect.
message(Sys.time(), " Start power simulation effect size 2x")
power_sim <- model_data %>% 
  mutate(power2 = purrr::map(.x = training, custom_powersim, ef = 0.693, n = nsim)) %>%
  mutate(power1.5 = purrr::map(.x = training, custom_powersim, ef = 0.405, n = nsim)) %>%
  mutate(power1.4 = purrr::map(.x = training, custom_powersim, ef = 0.336, n = nsim)) %>%
  mutate(power1.3 = purrr::map(.x = training, custom_powersim, ef = 0.262, n = nsim)) %>%
  mutate(power1.2 = purrr::map(.x = training, custom_powersim, ef = 0.182, n = nsim)) %>%
  mutate(power1.1 = purrr::map(.x = training, custom_powersim, ef = 0.095, n = nsim)) 

message(Sys.time(), " Saving output")
saveRDS(power_sim, "outputs/power_sim.RDS")

