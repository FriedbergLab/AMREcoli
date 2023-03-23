## ---------------------------
## Purpose of script: 
##  Perform a logistic regression calculating significance of gene 
##  presence / absence on resistance phenotype.
## Author: Henri Chung
## ---------------------------

# Load required packages and set up folder directory
#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(lme4)
rm(list = ls())

# Set up data folder
dataFolder <- "data/"
dataFile <- "tidy/samples.RDS"
referenceFile <- "tidy/reference.RDS"

#
options(pillar.sigfigs=2)

# Read in tidy data.
mydata <- readRDS(paste(dataFolder, dataFile, sep = ""))
reference <- readRDS(paste(dataFolder, referenceFile, sep = ""))

gene_classes <- reference$classes %>%
  select(test_type_antibiotic_class, mic_id) %>%
  rename("resistance_class" = test_type_antibiotic_class) %>%
  unique() 

raw_results <- readRDS("logreg_results.RDS")

fit_errors <- function(.x, errors = TRUE){
  fits <- unlist(.x$fit)
  ind <- unlist(lapply(fits, is.character))
  if(errors == TRUE){return(.x[!ind,])}else{return(.x[ind,])}
}

results <- lapply(raw_results, function(x){lapply(x, fit_errors)})
errors <- lapply(raw_results, function(x){lapply(x, fit_errors, errors = FALSE)})

usda_singlegene_results <- results[["usda_singlegene_results"]]
usda_geneclass_results <- results[["usda_geneclass_results"]]
usda_binclass_results <- results[["usda_binclass_results"]]

# calculate p-values for fits.
custom_pval <- function(x, y = 2){
  res <- x %>%
      filter(!is.na(fit)) %>%
      mutate(summ = map(fit, function(.x){summary(.x)})) %>%
      mutate(coefs = map(summ, function(.x){.x["coefficients"]})) %>%
      mutate(pval = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% dplyr::select('Pr(>|z|)') %>% t())})) %>%
      mutate(pval.factor = map(pval, function(.x){.x[,y]})) %>%
      unnest(pval.factor) %>%
      mutate(pval.adj = p.adjust(pval.factor, method = "holm")) %>%
      filter(pval.adj < 0.05)
  return(res)
}

# filter to only significant results
usda_singlegene_sig <- lapply(usda_singlegene_results, function(.x){custom_pval(.x)})
usda_geneclass_sig <- lapply(usda_geneclass_results, function(.x){custom_pval(.x)})
usda_binclass_sig <- lapply(usda_binclass_results, function(.x){custom_pval(.x)})


sigs <- list(usda_singlegene_sig, usda_geneclass_sig, usda_binclass_sig)
names(sigs) <- c("usda_singlegene_sig", "usda_geneclass_sig", "usda_binclass_sig")
lapply(results, function(x){lapply(x, dim)})
lapply(sigs, function(x){lapply(x, dim)})

# CLSI Animal Prior / USDA Animal Results 
usda_singlegene_results[[1]] %>% nrow() # total number of tests
usda_singlegene_sig[[1]] %>% nrow() # significant results
usda_singlegene_sig[[1]]$fit[[1]] %>% summary() # example summary results.


#Function to extract model coefficients, extract random effect, calculate odds ratio, calculate confidence intervals.
custom_unpack <- function(.x, gene_class = gene_classes){
  res <- .x %>%
    left_join(gene_class, by = "mic_id") %>% # add gene class labels
    mutate(combo = paste(mic_id, gene_group)) %>%
    mutate(estimate = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% dplyr::select("Estimate") %>% t())}))  %>% # extract coefficient
    mutate(estimate.factor = map(estimate, function(.x){.x[,2]})) %>% #
    mutate(std = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% dplyr::select(contains("Std")) %>% t())})) %>% # extract std of coefficient
    mutate(std.factor = map(std, function(.x){.x[,2]})) %>%
    mutate(ranefs = map(fit, function(.x){temp <- .x %>% ranef() %>% as_tibble() %>% pull(grp); return(paste(temp, collapse = ", "))}))  %>% # extract random effects.
    mutate(varCorr = map(fit, function(.x){temp <- VarCorr(.x) %>% as.data.frame() %>% select(vcov, sdcor); return(temp[2]) })) %>% # extract variance and std of random effect
    unnest(cols = c(ranefs, varCorr, estimate.factor, std.factor)) %>% 
    mutate(exp_coef = exp(estimate.factor), exp_coef_lower = exp(estimate.factor-std.factor), exp_coef_upper = exp(estimate.factor+std.factor)) %>% # exponentiate coefficient, calculate upper and lower bound
    mutate(odds_ratio = paste(round(exp_coef, 2), " (", round(exp_coef_lower,2),  ",", round(exp_coef_upper,2), ")", sep = "")) %>% # reformat values
    unique() %>%
    dplyr::select(resistance_class, mic_id, gene_group, odds_ratio, pval.adj) %>% # select desired columns
    arrange(resistance_class, mic_id, pval.adj) 
    colnames(res) <- c("Class", "Antibiotic", "Gene", "Odds Ratio", "Pval") 
    return(res)
}

# same function as above but reformatted for geneclass models

custom_unpack2 <- function(.x, gene_class = gene_classes){

  # if no rows in input, return empty dataframe
    if(nrow(.x) < 1){
	temp <- as.data.frame(matrix(ncol = 5, nrow = 1)) 
	colnames(temp) <- c("Class", "Antibiotic", "Gene", "Odds Ratio", "Pval") 
	return(temp)
	}
  res <- .x %>%
    left_join(gene_class, by = "mic_id") %>% # add gene class labels
    mutate(combo = paste(mic_id, gene_group)) %>%
    mutate(estimate = map(coefs, function(.x){return(.x[[1]] %>% as.data.frame() %>% dplyr::select("Estimate") %>% t())}))  %>% # extract coefficient
    mutate(estimate.factor = map(estimate, function(.x){.x[,2]})) %>%
    mutate(std = map(fit, function(.x){return(exp(confint(.x, devtol = 1e-7)) %>% t() %>% as_tibble() %>% dplyr::select(matches("^values[0-9]*|\\.L")))})) %>% # extract std of coefficient
    mutate(exp_coef_lower = map(std, function(.x){colnames(.x) <- "exp_coef_lower"; return(.x[1,])}), exp_coef_upper = map(std, function(.x){colnames(.x) <- "exp_coef_upper"; return(.x[2,])})) %>%   # calculate upper and lower bound
    unnest(c(exp_coef_lower, exp_coef_upper)) %>%    
    mutate(ranefs = map(fit, function(.x){temp <- .x %>% ranef() %>% as_tibble() %>% pull(grp); return(paste(temp, collapse = ", "))}))  %>%
    mutate(varCorr = map(fit, function(.x){temp <- VarCorr(.x) %>% as.data.frame() %>% dplyr::select(vcov, sdcor); return(temp[2]) })) %>% # extract variance and std of random effect
    unnest(cols = c(ranefs, varCorr, estimate.factor)) %>%
    mutate(exp_coef = exp(estimate.factor)) %>%
    mutate(odds_ratio = paste(round(exp_coef, 2), " (", round(exp_coef_lower,2),  ",", round(exp_coef_upper,2), ")", sep = "")) %>% # reformat values
    unique() %>%
    dplyr::select(resistance_class, mic_id, gene_group, odds_ratio, pval.adj) %>%
    arrange(resistance_class, mic_id, pval.adj) 
    colnames(res) <- c("Class", "Antibiotic", "Gene", "Odds Ratio", "Pval") 
    return(res)
}

# Reduce Beta, Std, Odds Ratio into One column
#TESTING ZONES###

m1_data <- lapply(usda_singlegene_sig, custom_unpack2)
m2_data <- lapply(usda_geneclass_sig, custom_unpack)
m3_data <- lapply(usda_binclass_sig, custom_unpack2)

names(m1_data) <- names(m2_data) <- names(m3_data) <- c("clsi", "wt")


write_csv(bind_rows(m1_data, .id = "Group"), "outputs/usda_singlegene_sig.csv")
write_csv(bind_rows(m2_data, .id = "Group"), "outputs/usda_geneclass_sig.csv")
write_csv(bind_rows(m3_data, .id = "Group"), "outputs/usda_binclass_sig.csv")

#######

# number of tested results for each model and dataset.
res_df <- stack(unlist(lapply(results, function(x){lapply(x, nrow)}))) %>%
  mutate(ind = gsub("usda_|_results", "" , as.character(ind)))
res_df

# number of signficant results for each model and dataset.
sig_df <- stack(unlist(lapply(sigs, function(x){lapply(x, nrow)}))) %>%
  mutate(ind = gsub("usda_|_sig", "" , as.character(ind)))
sig_df

# combine total, and significant results
m4_data <- res_df %>%
  as_tibble() %>%
  rename(n_test = "values") %>%
  merge(sig_df) %>%
  rename(n_sig = "values") %>%
  mutate(n_insig = n_test - n_sig) %>%
  separate(ind, into = c("model", "group"), sep = "\\.") 

# figure for insignificant and significant results.
p4 <- m4_data %>%
  mutate(model = ifelse(model == "binclass", "Binary Group", model)) %>%
  mutate(model = ifelse(model == "singlegene", "Single Gene", model)) %>%
  mutate(model = ifelse(model == "geneclass", "Group", model)) %>% 
  mutate(model = factor(model, levels = c("Single Gene", "Group", "Binary Group"))) %>%
  mutate(group = ifelse(group == "iq75", "wt", "clsi")) %>%
  #mutate(n_insig = ifelse(n_insig > 0, n_insig, NA), n_sig = ifelse(n_sig > 0, n_sig, NA)) %>%
  dplyr::select(model,group, n_sig, n_insig) %>%
  pivot_longer(cols = c("n_sig", "n_insig"), names_to = "class", values_to =  "value") %>%
  ggplot(aes(x = model, y = value, fill = class, label = value)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(position = position_dodge2(width = 1, preserve = "single"), angle = 0, vjust=0, hjust=0) + 
  labs(fill = "P value") + 
  scale_fill_manual(labels = c(">0.05", "<0.05"), values = c("lightgray", "darkgray")) + 
  ylim(0, 60) + 
  xlab("") + ylab("# Models") +
  facet_wrap(~group) +
  theme_bw()

# plot
pdf("outputs/p4.pdf", height = 6, width = 8)
p4
dev.off()
#######

# function to extract estimate of andom random effect from tested models
custom_animal <- function(.x){
  if(nrow(.x)== 0){return(as_tibble(0))}
  res <- .x %>%
    select(mic_id, gene_group, fit) %>%
    mutate(ranefs = map(fit, function(.x){return(.x %>% ranef() %>% as.data.frame() %>% select(grp, condval, condsd))}))  %>%
    mutate(varCorr = map(fit, function(.x){VarCorr(.x) %>% as.data.frame() %>% select(vcov, sdcor)})) %>% 
    unnest(cols = c(ranefs, varCorr)) %>%
    filter(vcov != 0) %>%
    group_by(mic_id, gene_group, grp) %>%
    arrange(desc(condval))
  return(res)
}

a1_data <- lapply(usda_singlegene_sig, custom_animal) %>% bind_rows(.id = "group")
a2_data <- lapply(usda_geneclass_sig, custom_animal) %>% bind_rows(.id = "group")
a3_data <- lapply(usda_binclass_sig, custom_animal) %>% bind_rows(.id = "group")


# combine animal random effect data for each model into single dataframe.
a_list <- list(a1_data, a2_data, a3_data)
names(a_list) <- c("Single Gene", "Group", "Binary Group")
a_data <- bind_rows(a_list, .id = "model") %>%
  ungroup() %>%
  filter(!is.na(mic_id)) %>%
  mutate(group = ifelse(group == "ItoR", "clsi", "wt")) %>%
  mutate(grp = as.character(grp)) %>%
  mutate(grp = as.factor(grp)) %>%
  mutate(grp2 = case_when( # convert species names to common name
    grepl("duck", .$grp, ignore.case = TRUE) ~ "Anas platyrhynchos",
    grepl("cat\\>", .$grp, ignore.case = TRUE) ~ "Felis catus",
    grepl("cattle", .$grp, ignore.case = TRUE) ~ "Bos taurus",
    grepl("dog", .$grp, ignore.case = TRUE) ~ "Canis familiaris",
    grepl("horse", .$grp, ignore.case = TRUE) ~ "Equus caballus",
    grepl("chicken", .$grp, ignore.case = TRUE) ~ "Gallus gallus",
    grepl("turkey", .$grp, ignore.case = TRUE) ~ "Meleagris gallopavo",
    grepl("swine", .$grp, ignore.case = TRUE) ~ "Sus domesticus")) %>%
  mutate(grp3 = case_when( # convert species names to common name
    grepl("duck", .$grp, ignore.case = TRUE) ~ "Ducks",
    grepl("cat\\>", .$grp, ignore.case = TRUE) ~ "Cats",
    grepl("cattle", .$grp, ignore.case = TRUE) ~ "Cattle",
    grepl("dog", .$grp, ignore.case = TRUE) ~ "Dogs",
    grepl("horse", .$grp, ignore.case = TRUE) ~ "Horses",
    grepl("chicken", .$grp, ignore.case = TRUE) ~ "Chickens",
    grepl("turkey", .$grp, ignore.case = TRUE) ~ "Turkeys",
    grepl("swine", .$grp, ignore.case = TRUE) ~ "Swine")) %>%
  mutate(grp3 = factor(grp3, levels = c("Swine", "Cattle", "Chickens", "Turkeys", "Horses", "Dogs", "Cats"))) %>%
  mutate(grp2 = as.factor(grp)) %>%
  mutate(grp2 = factor(grp, levels = c("Sus domesticus", "Bos taurus", "Gallus gallus", "Meleagris gallopavo", "Equus caballus", "Canis familiaris", "Felis catus"))) %>%
  group_by(model, mic_id, gene_group) %>%
  filter(n() > 5) # 

# check mean random effect for different models
a_data %>% group_by(model, mic_id, grp3) %>% summarize(mean_condval = mean(condval)) %>% arrange(desc(n))


# plot mean animal effect.
pdf("outputs/p5.pdf", height = 6, width = 8)
p5 <- ggplot(a_data, aes(x = grp3, y = condval, fill = grp3)) + 
  geom_boxplot() +
  facet_wrap(~model) +
  labs(fill = "Animal Source") + ylab("Conditional Mean") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p5
dev.off()

# wfunction to compare the fit of models with and without animal random effect
compare_aic <- function(.x, .y){

  # filter model results which returned errors
  new_x <- filter(.x, !map_lgl(fit, is.character)) 
  new_y <- filter(.y, !map_lgl(fit, is.character)) 

  # find only models which were fit to both datasets
  ind <- intersect(new_x$combo, new_y$combo)
  fit_x <- filter(.x, combo %in% ind) %>% pull(fit)
  fit_y <- filter(.y, combo %in% ind) %>% pull(fit)

  # compare two models with chi-square
  aic_list <- mapply(anova, fit_x, fit_y, SIMPLIFY = FALSE) %>%
    lapply(as_tibble)
  names(aic_list) <- unlist(.x$combo)
  aic_df <- bind_rows(aic_list, .id = "combo")

  # calculate significant model differences
  sig_aic <- aic_df %>%
    filter(`Pr(>Chisq)` <= 0.05) %>%
    pull(combo)
    if(length(sig_aic) == 0){message("no significant results)")}

  # return df of results
  res <- aic_df %>%
    #filter(aic_df, combo %in% sig_aic ) %>%
    select(npar, AIC, combo) %>%
    group_by(combo) %>%
    mutate(npar = ifelse(npar != max(npar), "random", "normal")) %>%
    pivot_wider(names_from = npar, values_from = AIC) %>%
    mutate(compare = ifelse(random < normal, TRUE, FALSE))
  return(res)
}

# extract random effect model results
usda_singlegene_results <- results[["usda_singlegene_results"]]
usda_geneclass_results <- results[["usda_geneclass_results"]]
usda_binclass_results <- results[["usda_binclass_results"]]

# extract results of models w/o random effect
results_nore <- readRDS("logreg_results_nore.RDS")
usda_singlegene_results_nore <- results_nore[["usda_singlegene_results_nore"]]
usda_geneclass_results_nore <- results_nore[["usda_geneclass_results_nore"]]
usda_binclass_results_nore <- results_nore[["usda_binclass_results_nore"]]

# compare AIC for each model type.
compare_aic(usda_singlegene_results[[1]], usda_singlegene_results_nore[[1]]) %>%
  filter(compare == FALSE)
compare_aic(usda_singlegene_results[[2]], usda_singlegene_results_nore[[2]]) %>%
  filter(compare == FALSE)

compare_aic(usda_geneclass_results[[1]], usda_geneclass_results_nore[[1]])
compare_aic(usda_geneclass_results[[2]], usda_geneclass_results_nore[[2]])

compare_aic(usda_binclass_results[[1]], usda_binclass_results_nore[[1]])
compare_aic(usda_binclass_results[[2]], usda_binclass_results_nore[[2]])

#
usda_singlegene_sig_nore <- lapply(usda_singlegene_results_nore, function(.x){custom_pval(.x)})
usda_geneclass_sig_nore <- lapply(usda_geneclass_results_nore, function(.x){custom_pval(.x)})
usda_binclass_sig_nore <- lapply(usda_binclass_results_nore, function(.x){custom_pval(.x)})

# compare fits between models
binclass_itor_combos <- usda_binclass_results_nore[[1]]$combo

geneclass_compare <- usda_geneclass_results_nore[[1]] %>% 
  filter(combo %in% binclass_itor_combos) %>%
  filter(!map_lgl(fit, is.character)) 

singlegene_compare <- usda_singlegene_results_nore[[1]] %>%
  left_join(select(rename(reference$gene_metadata, gene_group = "gene"), c(gene_group, resistance_class))) %>%
  mutate(combo = paste(resistance_class, mic_id)) %>%
  filter(combo %in% binclass_itor_combos) %>%
  group_by(combo) %>%
  unique() %>%
  mutate(AIC = map(fit, AIC)) %>%
  unnest(AIC) %>%
  filter(AIC == min(AIC))

binclass_compare <- filter(usda_binclass_results_nore[[1]], combo %in% singlegene_compare$combo)

anova(binclass_compare$fit[[1]], geneclass_compare$fit[[1]], singlegene_compare$fit[[1]])



binclass_itor_combos <- usda_binclass_results_nore[[2]]$combo

singlegene_compare <- usda_singlegene_results_nore[[2]] %>%
  left_join(select(rename(reference$gene_metadata, gene_group = "gene"), c(gene_group, resistance_class))) %>%
  mutate(combo = paste(resistance_class, mic_id)) %>%
  filter(combo %in% binclass_itor_combos) %>%
  unique() %>%
  mutate(AIC = map(fit, AIC)) %>%
  unnest(AIC) %>%
  group_by(combo) %>%
  filter(AIC == min(AIC)) %>% 
  arrange(desc(combo))

binclass_compare <- filter(usda_binclass_results_nore[[2]], combo %in% singlegene_compare$combo) %>% 
  arrange(desc(combo))

geneclass_compare <- usda_geneclass_results_nore[[2]] %>% 
  filter(combo %in% binclass_itor_combos & combo %in% singlegene_compare$combo) %>%
  filter(!map_lgl(fit, is.character)) %>% 
  arrange(desc(combo))


fit_comparison <- list()
for(i in 1:nrow(singlegene_compare)){
  #fit_comparison[[i]] <- anova(binclass_compare$fit[[i]], geneclass_compare$fit[[i]], singlegene_compare$fit[[i]], test = "LRT")
  fit_comparison[[i]] <- tibble(binclass = AIC(binclass_compare$fit[[i]]), geneclass =  AIC(geneclass_compare$fit[[i]]), singlegene =  AIC(singlegene_compare$fit[[i]]))
}
names(fit_comparison) <- singlegene_compare$combo
aic_df <- bind_rows(fit_comparison, .id = "combo") %>%
  pivot_longer(cols = c("binclass", "geneclass", "singlegene"), names_to = "model", values_to = "AIC") %>%
  group_by(combo) %>%
  filter(AIC == min(AIC)) %>%
  arrange(combo)

filter(aic_df, model == "singlegene") %>% left_join(singlegene_compare) %>% custom_pval() %>% select(pval.adj, everything())
filter(aic_df, model == "geneclass") %>% left_join(geneclass_compare) %>% custom_pval() %>% select(pval.adj, everything())
filter(aic_df, model == "binclass") %>% left_join(binclass_compare) %>% custom_pval() %>% select(pval.adj, everything())
write_csv(aic_df, "aic_df.csv")
anova(a, b, c)

