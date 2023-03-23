### ---------------------------
## Purpose of script: 
##  Post-hoc power analysis of single gene logistic regression 
## Author: Henri Chung
## ---------------------------

#setwd("E:/Projects/amr_analysis")
library(tidyverse)
library(lme4)
library(simr)
rm(list = ls())

usda_singlegene_curve <- readRDS("usda_singlegene_curve.RDS")
usda_geneclass_curve <- readRDS("usda_geneclass_curve.RDS")
usda_binclass_curve <- readRDS("usda_binclass_curve.RDS")


unlist2 = function(x, sep = "_") {
  #save top names
  top_names = names(x)
  x = unname(x)

  #flatten
  x2 = unlist(x)

  #add prefix
  #determine how many prefixes to add of each
  lengths_top = sapply(x, length)
  prefixes = rep(top_names, times = lengths_top)
  names(x2) = paste0(prefixes, sep, names(x2))
  x2
}
curve_list <- list(usda_singlegene_curve, usda_geneclass_curve, usda_binclass_curve)
names(curve_list) <- c("singlegene", "geneclass", "binclass")

mean_power <- lapply(curve_list, function(.x){
  lapply(.x, function(.y){
    res2 <- lapply(.y, function(.z){
      temp <- .z %>%
        unnest(power_val) %>%
        mutate(power_val = power_val) %>% 
        pull(power_val)
      res <- lapply(temp, mean)
      return(unlist(res))})
    return(unlist2(res2))
})}) %>%
  lapply(unlist2) %>%
  unlist2() %>%
  stack() %>%
  as_tibble() %>%
  separate(ind, into = c("model", "effect_size","dataset"), sep = "_") 

test <- mean_power %>% group_by(model, effect_size) %>%
  summarize(mean = mean(values))

p1 <- mean_power %>% 
  group_by(model, effect_size, dataset) %>%
  summarize(mean = mean(values)) %>%
  ggplot(aes(x = effect_size, color = dataset, y = mean)) +
  geom_point() +
  facet_wrap(~model)
pdf("test2.pdf")
p1
dev.off()

p1 <- mean_power %>% 
  ggplot(aes(x = effect_size, color = dataset, y = values)) +
  geom_boxplot() +
  facet_wrap(~model)
pdf("test.pdf")
p1
dev.off()

mean_curve <- lapply(curve_list, function(.x){
  df <- lapply(.x, function(.y){
    res2 <- lapply(.y, function(.z){
      temp <- .z %>%
        pull(curve) %>%
        lapply(function(.q){
          res3 <- lapply(.q$ps, function(.t){
            return(.t$x/10)
            }) %>% unlist()
          names(res3) <- .q$nlevels
          return(stack(res3))
          })
        names(temp) <- unlist(.z$mic_id)
        return(bind_rows(temp, .id = "combo"))})
    return(bind_rows(res2, .id = "dataset"))})
  return(bind_rows(df, .id = "effect"))}) %>% 
bind_rows(.id = "model")


test <- mean_curve %>% 
  group_by(model, effect, ind) %>% 
  summarize(mean = mean(values)) %>% 
  filter(mean > 80) %>% 
  arrange(effect)
