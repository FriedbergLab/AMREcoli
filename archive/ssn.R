#setwd("~/amr_analysis")
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

K = 500 # number of data rows to use. Low number for testing. 
ssn = 1000 # number of iterations for stability selection.

model_list <- readRDS("outputs/model_list.RDS")
best_model = "Binary-Count"
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
ssn = 1000
library(furrr)
library(future)
plan(multicore, workers = 36)
tictoc::tic()
ss_data <- model_list[[best_model]] %>%
  mutate(coefs = furrr::future_map(training, function(.z){stack(table(unlist(replicate(en_selection(.z), n = ssn)))/ssn)}, .options = furrr_options(seed = TRUE )))
tictoc::toc()
saveRDS(ss_data, "ss_models.RDS")

quit()