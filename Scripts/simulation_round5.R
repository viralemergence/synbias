#### Set up ####
library(poems)
library(tidyverse)
library(furrr)
library(Rcpp)
library(qs)
library(vroom)
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
plan(multisession, workers = 2)
set.seed(306)
inputs <- sample_data3 %>% rowwise() %>% group_split() %>%
  map(df_to_inputs)
strains <- inputs |> 
  map("interactions") |> 
  map(function(x) {
    x <- x - 1
    t(sapply(1:nrow(x), function(i) {
      c(sum(abs(c(x[i, ][x[i, ] < 0], x[, i][x[, i] < 0]))), 
        sum(c(x[i, ][x[i, ] > 0], x[, i][x[, i] > 0])))
    }))
  }) |> 
  map(~ matrix(.x, ncol = 2, byrow = TRUE)) |>
  map(as.data.frame) |>
  map(\(x) cbind(x, Strain = 1:nrow(x))) |>
  list_rbind(names_to = "Simulation") |>
  set_names(c("Simulation", "Competition", "Facilitation", "Strain"))
vroom_write(strains, "Data/simulation_round5/per_strain.csv")
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
sim1000 <- paste0(getwd(), data_dir, "sim1000.qs") |> qread()
sim1000 |> apply(c(2, 3), sum) |> array_branch(1) |> map(\(x) x >= 500) |>
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
write_csv(sample_data3, "Data/simulation_round5/simulation_round5.csv")
strains |> 
  left_join(sample_data3, by = c("Simulation" = "id")) |>
  ggplot(aes(x = Competition/strains, y = Facilitation/strains)) +
  geom_bin_2d() +
  scale_fill_viridis_c() +
  coord_fixed() +
  cowplot::theme_cowplot()
#### Figures ####
all_results <- "Data/simulation_round5/simulation_round5.csv" |> read_csv() |>
  mutate(percent_error = ((strains-richness_count)/strains)*100)
binned_data <- all_results %>%
  mutate(strains_mid = (floor(strains / 5) * 5) + 2.5,
         cf_ratio_mid = (floor(cf_ratio / 0.05) * 0.05) + 0.025) %>%
  group_by(strains_mid, cf_ratio_mid) %>%
  summarize(avg_percent_error = mean(percent_error, na.rm = TRUE))
# Define breaks and labels for the y-axis
y_breaks <- c(min(binned_data$cf_ratio_mid), max(binned_data$cf_ratio_mid))
y_labels <- c("100%\nfacilitative", "100%\ncompetitive")

# Creating the plot
ggplot(binned_data, aes(x = strains_mid, y = cf_ratio_mid, fill = avg_percent_error)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(-0.5, 10)) +
  labs(title = "Heatmap: Strains vs C:F Ratio",
       x = "Strains", fill = "% Error") +
  theme_minimal_grid() +
  scale_x_continuous(labels = scales::label_number()) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(axis.title.y = element_blank())

binned_data <- all_results |>
  mutate(interaction_strength_mid = (floor(interaction_strength / 0.03) * 0.03) + 0.015,
         cf_ratio_mid = (floor(cf_ratio / 0.05) * 0.05) + 0.025) %>%
  group_by(interaction_strength_mid, cf_ratio_mid) %>%
  summarize(avg_percent_error = mean(percent_error, na.rm = TRUE))
y_breaks <- c(min(binned_data$cf_ratio_mid), max(binned_data$cf_ratio_mid))
y_labels <- c("100%\nfacilitative", "100%\ncompetitive")

# Creating the plot
ggplot(binned_data, aes(x = interaction_strength_mid, y = cf_ratio_mid, fill = avg_percent_error)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(-0.5, 10)) +
  labs(title = "Heatmap: Interaction Strength vs CF Ratio",
       x = "Interaction Strength", fill = "% Error") +
  theme_minimal_grid() +
  scale_x_continuous(labels = scales::label_number()) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(axis.title.y = element_blank())

all_results |>
  ggplot(aes(x = sample_prop*100, y = percent_error)) +
  geom_point() +
  theme_minimal_grid() + xlab("% of Population Sampled") +
  ylab("% Error")

binned_data <- all_results |>
  mutate(sample_prop_mid = (floor(sample_prop / 0.04) * 0.04) + 0.02,
         cf_ratio_mid = (floor(cf_ratio / 0.05) * 0.05) + 0.025) %>%
  group_by(sample_prop_mid, cf_ratio_mid) %>%
  summarize(avg_percent_error = mean(percent_error, na.rm = TRUE))
y_breaks <- c(min(binned_data$cf_ratio_mid), max(binned_data$cf_ratio_mid))
y_labels <- c("100%\nfacilitative", "100%\ncompetitive")

# Creating the plot
ggplot(binned_data, aes(x = sample_prop_mid*100, y = cf_ratio_mid, fill = avg_percent_error)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(-0.5, 50)) +
  labs(title = "Heatmap: Sample Proportion vs CF Ratio",
       x = "% of Population Sampled", fill = "% Error") +
  theme_minimal_grid() +
  scale_x_continuous(labels = scales::label_number()) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  theme(axis.title.y = element_blank())

#### Sensitivity Analysis ####
library(tidymodels)
tidymodels_prefer()
time35_split <- initial_split(all_results)
time35_spec <- rand_forest() |> set_mode("regression")
