#### Set up ####
library(poems)
library(tidyverse)
library(furrr)
library(Rcpp)
library(qs)
source("Scripts/simulation_functions.R")
sourceCpp("Scripts/ferrari.cpp")
options(future.globals.maxSize= 9967869120)
data_dir <- "/Data/simulation_round5/"
#### Simulation ####
synbias_lhs <- LatinHypercubeSampler$new()
synbias_lhs$set_uniform_parameter("interaction_strength", 0, 0.75)
synbias_lhs$set_class_parameter("priority_effects", c(TRUE, FALSE))
synbias_lhs$set_uniform_parameter("cf_ratio")
synbias_lhs$set_uniform_parameter("strains", 10, 200, decimals = 0)
synbias_lhs$set_uniform_parameter("sample_prop", 0.01, 1)
sample_data3 <- synbias_lhs$generate_samples(10000, random_seed = 216) %>%
  rowid_to_column("id")
plan(sequential)
inputs <- sample_data3 %>% rowwise() %>% group_split() %>%
  map(df_to_inputs)
future_walk(1:length(inputs), \(x) {
  filename <- paste0(getwd(),
                     "/Data/simulation_round5/sim",
                     sample_data3$id[x],
                     ".qs")
  if (!file.exists(filename)) {
    sourceCpp("Scripts/ferrari.cpp")
    x <- inputs[[x]]
    mat <- ferrari(x$initial_pop, x$interactions, 2e-04, 100)
    qsave(mat, filename)
  }
}, .progress = T)
#### 50 percent prevalence points ####
# Try it with one example first
sim1 <- paste0(getwd(), data_dir, "sim1.qs") |> qread()
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
plan(multisession, workers = 8)
timepoints <- 1:10000 |>
  map(\(x) paste0(getwd(), data_dir, "sim", x, ".qs")) |>
  future_map(prevalence50, .progress = T)
timepoints |> flatten_int() |>
  hist(main = "Histogram of 50% prevalence points",
       xlab = "Timestep when prevalence is 50%")
#### Sample diversity at median 50% prevalence point ####
diversity_sample <- function(path, timepoint, proportion) {
  path |> qread() %>%
    .[sample(1:1000, proportion*1000), , timepoint] |>
    colSums() |> map_lgl(\(x) x > 0) |> sum()
}
sample_data3$richness_count <-  1:10000 |>
  map(\(x) paste0(getwd(), data_dir, "sim", x, ".qs")) |>
  future_map2_int(sample_data3$sample_prop,
                  diversity_sample,
                  timepoint = 35,
                  .progress = T)
