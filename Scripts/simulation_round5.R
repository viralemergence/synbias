#### Set up ####
library(poems)
library(tidyverse)
library(furrr)
library(Rcpp)
library(qs)
source("Scripts/simulation_functions.R")
options(future.globals.maxSize= 9967869120)
data_dir <- "Data/simulation_round5/"
#### Simulation ####
synbias_lhs <- LatinHypercubeSampler$new()
synbias_lhs$set_uniform_parameter("interaction_strength", 0, 0.75)
synbias_lhs$set_class_parameter("priority_effects", c(TRUE, FALSE))
synbias_lhs$set_uniform_parameter("cf_ratio")
synbias_lhs$set_uniform_parameter("strains", 10, 200, decimals = 0)
sample_data3 <- synbias_lhs$generate_samples(10000, random_seed = 216) %>%
  rowid_to_column("id")
plan(sequential)
inputs <- sample_data3 %>% rowwise() %>% group_split() %>%
  map(df_to_inputs)
future_walk(1:length(inputs), \(x) {
  filename <- paste0(data_dir, "sim/", sample_data3$id[x], ".qs")
  if (!file.exists(filename))
    sourceCpp("Scripts/ferrari.cpp")
    x <- inputs[[x]]
    mat <- ferrari(x$initial_pop, x$interactions, 2e-04, 100)
    qsave(mat, filename)
  }, .progress = T)
#### 50 percent prevalence points ####
# Try it with one example first
sim1 <- data_dir |> paste0("sim1.qs") |> qread()
sim1 |> apply(c(2, 3), sum) |> array_branch(1) |> map(\(x) x >= 500) |>
  map_int(match, x = TRUE)
# Scale it up
prevalence50 <- function(path) {
  path |> qread() |>
    apply(c(2, 3), sum) |>
    array_branch(1) |>
    map(\(x) x >= 500) |>
    map_int(match, x = TRUE)
}
timepoints <- 1:10000 |> map(\(x) paste0(data_dir, "sim/", x, ".qs")) |>
  map(prevalence50)
