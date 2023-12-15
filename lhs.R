library(poems)
library(tidyverse)
synbias_lhs <- LatinHypercubeSampler$new()
synbias_lhs$set_uniform_parameter("interaction_strength", 0, 0.75)
synbias_lhs$set_class_parameter("priority_effects", c(TRUE, FALSE))
synbias_lhs$set_uniform_parameter("cf_ratio")
synbias_lhs$set_uniform_parameter("strains", 10, 100, decimals = 0)
sample_data <- synbias_lhs$generate_samples(10000, random_seed = 324) %>%
  rowid_to_column("id")
# add code here to generate interaction matrices based on sample data
df_to_matrix <- function(df) {
  int_matrix <- matrix(nrow = df$strains, ncol = df$strains)
  if (df$priority_effects) {
    comp_int <- round(df$cf_ratio*df$strains^2)
    fac_int <- df$strains^2 - comp_int
    comp_indices <- sample(length(int_matrix), comp_int)
    int_matrix[comp_indices] <- runif(comp_int, 1 - df$interaction_strength, 1)
    int_matrix[-comp_indices] <- runif(fac_int, 1, 1 + df$interaction_strength)
  } else {
    comp_int <- round(df$cf_ratio*sum(lower.tri(int_matrix)))
    fac_int <- sum(lower.tri(int_matrix)) - comp_int
    comp_indices <- sample(sum(lower.tri(int_matrix)), comp_int)
    int_matrix[lower.tri(int_matrix)][comp_indices] <- runif(comp_int,
                                               1 - df$interaction_strength, 1)
    int_matrix[upper.tri(int_matrix)][-comp_indices] <- int_matrix[lower.tri(int_matrix)]
  }
  diag(int_matrix) <- 1
  return(int_matrix)
}
