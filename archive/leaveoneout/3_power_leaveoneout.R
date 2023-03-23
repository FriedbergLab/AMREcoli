## ---------------------------
## Purpose of script: 
##  Approximate the power to detect effect size for individual predictors.
## Author: Henri Chung
## ---------------------------

library(tidyverse)
library(WebPower)
library(parallel)
library(grid)
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


# Power Calculation
#####################

# specific effect size to detect
B = 0.7

# take training data from stability selection output
model_data <- readRDS("outputs/ss_models_leaveoneout.RDS") %>%  dplyr::select(mic_id, breakpoint, training, en)

# function to iterate over reach column in the predictor matrix.
power_calc <- function(.x, beta = 0.07){
  resp <- .x$custom_phenotype
  preds <- dplyr::select(.x, -custom_phenotype)
  l <- length(resp)
  res <- list()
  for(i in 1:ncol(preds)){res[[i]] <- log_power(preds[,i], resp = resp, l = l, beta = beta)}
  names(res) <- colnames(preds)
  return(stack(res))
}

# function to calculate power for a single predictor and response.
log_power <- function(.y, resp, l, beta){
      if(sum(.y > 0 & .y < 1) > 0){
        mean_y <- mean(.y[.y > 0])

        p0_i <- (sum(resp[.y >= mean_y])/l) + 1e-06
        a0 = log( (p0_i / (1-p0_i) )) 
        p1_i = exp(a0+beta)/ ( 1 + exp(a0+beta) ) + 1e-06

        pow <- wp.logistic(n = l, p0 = p0_i, p1 = p1_i, alpha = 0.05,power = NULL, family = "normal")
      }else{
        p0_i <- (sum(resp[.y == 0])/l) + 1e-06
        a0 = log( (p0_i / (1-p0_i) )) 
        p1_i = (exp(a0+beta)/ ( 1 + exp(a0+beta) ) )+ 1e-06

        pow <- wp.logistic(n = l, p0 = p0_i, p1 = p1_i, alpha = 0.05,power = NULL, family = "Bernoulli", parameter = 0.5) 
    }
    return(pow$power)
  }

# join power results together with predictors
power_join <- function(.x, .y){
  sig <- as.matrix(.x$beta) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(ind = "rowname") %>%
    left_join(.y, by = "ind")
}

# save results of individual effect size power sim
power_sim <- model_data %>% 
  mutate(power_vals = purrr::map(training, power_calc, beta = B)) %>%
  mutate(power_table = purrr::map2(en, power_vals, power_join))

message(Sys.time(), " Saving output")
saveRDS(power_sim, "outputs/power_sim_leaveoneout.RDS")

# measure power over gradient of effect sizes
power_grad <- model_data %>% 
  mutate(power_vals2 = purrr::map(training, power_calc, beta = 0.7)) %>%
  mutate(power_vals3 = purrr::map(training, power_calc, beta = 1.1)) %>%
  mutate(power_vals5 = purrr::map(training, power_calc, beta = 1.6)) %>%
  mutate(power_vals10 = purrr::map(training, power_calc, beta = 2.3)) %>%
  pivot_longer(cols = contains("power_vals"), names_to = "ef", values_to = "power_vals") %>%
  mutate(power_table = purrr::map2(en, power_vals, power_join))

# unnest power values
power_vals_grad <- power_grad %>%
  dplyr::select(mic_id, breakpoint, power_table, ef) %>%
  unnest(power_table) %>%
  rename(power = "values", gene = "ind") %>% 
  filter(gene != "(Intercept)") %>%
  mutate(ef = gsub("power_vals","", ef)) %>%
  mutate(ef = factor(ef, levels = c(2, 3, 5, 10)))
write_csv(power_vals_grad, "outputs/power_vals_leaveoneout.csv")

# plot distribution of power values over different effect sizes.
power_plot <- power_vals_grad %>%
  mutate(sig = ifelse(s0 == 0, "Not Important", "Important")) %>%
  mutate(combo = paste(mic_id, breakpoint)) %>%
  mutate(ef = factor(ef, levels = c(2, 3, 5, 10))) %>%
  ggplot(aes(x = ef, y = power, fill = sig)) + geom_boxplot() +
  xlab("\nEffect Size") + ylab("Power     ") + labs(fill = "") + 
  scale_fill_grey(start=0.8, end=0.4) +  
  #theme_classic() + 
  labs(title = "", fill = "Predictors") + 
  scale_shape_manual(values = c(0, 16)) +
  theme(
    legend.position="right", 
    strip.text.y.right = element_text(angle = 0), 
    text = element_text(size=20),
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(angle = 0, vjust = 0.5))
ggsave(power_plot, file = "outputs/power.tiff", height = 7, width = 7, dpi = 300, compression = "lzw")

# plot distribution of power alues for individual antibiotics.
power_plot_antibiotic <- power_plot$data %>%
  ggplot(aes(x = sig, y = power, fill = ef)) + geom_boxplot() + 
  facet_wrap(~combo) + 
  theme_classic()
ggsave("outputs/power_antibiotic_leaveoneout.tiff", power_plot_antibiotic, height = 7, width = 12, dpi = 300, compression = "lzw")

power_vals_grad %>% group_by(mic_id) %>% summarize(mean = mean(power))
power_vals_grad %>% group_by(mic_id) %>% summarize(mean = mean(power)) %>% arrange(mean)
quit()