### ---------------------------
## Purpose of script: 
##  Post-hoc power analysis of single gene logistic regression 
## Author: Henri Chung
## ---------------------------

#setwd("E:/Projects/amr_analysis")
suppressPackageStartupMessages(library(tidyverse, quiet = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(lme4, quiet = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(simr, quiet = TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(furrr, quiet = TRUE, warn.conflicts = FALSE))
rm(list = ls())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
N = 1000
message("POWER ANALYSIS WITH ", N, " SIMULATIONS", Sys.time())
message("READING DATA ", Sys.time())
results <- readRDS("logreg_results.RDS")
data_split_name <- case_when(
  args[1] == 1 ~ "usda_singlegene_results",
  args[1] == 2 ~ "usda_geneclass_results",
  args[1] == 3 ~ "usda_binclass_results")

if(args[1] == 2){
    test_name = simr::fixed("values.L")
    factor_name = "values.L"
}else if(args[1] == 3){
    test_name <- simr::fixed("values")
    factor_name = "values"
}else{
    test_name <- simr::fixed("values1")
    factor_name = "values1"
}
attr(factor_name, "text") <- function(...) NULL

usda_results <- results[[data_split_name]]

custom_curve <- function(.x, n_sim = 10, test_name = "",factor_name = "",  effect_size = 0.7, n_extend = 100){
    res <- .x %>%
      head() %>%
      rowwise() %>%
      rename(xfit = "fit") %>%
      mutate(temp = ifelse(is.character(unlist(xfit)), TRUE, FALSE)) %>%
      ungroup() %>% filter(temp == FALSE) %>%
      mutate(xfit2 = map(xfit, function(.z){.z@beta[names(fixef(.z)) %in% factor_name] <- effect_size; return(.z)})) %>%
      mutate(power = future_map(xfit2, function(.y){powerSim(.y, nsim =  n_sim, test = test_name)},.options = furrr_options(seed = 123))) %>%
      mutate(power_val = map(power, function(.v){.v$x / 10})) %>%
      mutate(sim_model = future_map2(data,xfit2, create_model)) %>%
      mutate(curve = future_map(sim_model, function(.x){powerCurve(.x, test=simr::fixed("values"), within="values", breaks=c(5,10,20,30,40), nsim = N)})) # 5,10,20,30,40
  return(res)
}

create_model <- function(.x, .y){
    # data
    n = length(unique(.x$values))
    sample_id <- paste("s", seq(1,n*2,1), sep = "_")
    host_animal_common <- factor(rep(c("1", "2"),  n))
    values <- factor(rep(seq(0, n-1,1), each = 2))
    phenotype <- rep(c(0,1), n)

    covars <- data.frame(
        sample_id = sample_id, 
        host_animal_common = host_animal_common, 
        phenotype = phenotype, 
        values = values)
    # parameters
    fixed <- fixef(.y)
    names(fixed)[1:2] <- c("(Intercept)", "values1")
    rand <- as.numeric(VarCorr(.y))
    #res <- attr(VarCorr(.y), "sc")
    if(args[1] == 2){
        fixed <- fixef(.y)
        names(fixed) <- c("(Intercept)", paste("values", 1:(n-1), sep = ""))
    }
     
    model <- makeGlmer(phenotype~values + (1|host_animal_common), family = binomial(link = "logit"), fixef=fixed, VarCorr=rand, data=covars)
    model_ext <- extend(model, within="values", n=40)
    return(model_ext)
}

options( warn = -1 ) # suppress warning messages. Simulations will  repeat "boundary (singular) fit: see ?isSingular"
plan(multisession)

effects <- c(0.7, 1.6, 2.3, 3.9, 4.6)
usda_curve <- list()
for(i in 1:length(effects)){
    message("STARTING POWER ANALYSIS FOR EFFECT SIZE ", effects[i], " ", Sys.time())
    usda_curve[[i]] <- suppressMessages(lapply(usda_results, function(.x){
        custom_curve(.x, n_sim = N, test_name = test_name, factor_name = factor_name, effect_size = effects[i])
    }))
}
names(usda_curve) <- effects

message("USDA SINGLEGENE COMPLETE", Sys.time())
saveRDS(usda_curve, paste(gsub("_results", "_curve.RDS", data_split_name)))
quit()

