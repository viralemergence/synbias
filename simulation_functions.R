disease_sample <- function(matrix, prop,
                           imperfect_assays = c(TRUE, FALSE)) {
  timesteps <- dim(matrix)[3]
  random_timestep <- sample(timesteps, 1)
  individuals <- round(prop*nrow(matrix))
  perfect_sample <- matrix |> _[individuals, , random_timestep]
  if (imperfect_assays) {
    false_positive <- runif(ncol(perfect_sample), 0, 0.01)
    false_negative <- runif(ncol(perfect_sample), 0, 0.07)
    strain_positives <- perfect_sample |> colSums()
    strain_negatives <- nrow(perfect_sample) - strain_positives
    for (c in 1:ncol(perfect_sample)) {
      perfect_sample[,c][perfect_sample[,c]==1] <- rbinom(strain_positives, 1,
                                                          1 - false_negative)
      perfect_sample[,c][perfect_sample[,c]==0] <- rbinom(strain_negatives, 1,
                                                          false_positive)
    }
  }
  return(perfect_sample)
}
# add code here to generate interaction matrices based on sample data
df_to_inputs <- function(df) {
  int_matrix <- matrix(1, nrow = df$strains, ncol = df$strains)
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
    int_matrix[lower.tri(int_matrix)][-comp_indices] <- runif(fac_int,
                                                              1,
                                                              1 + df$interaction_strength)
    int_matrix <- int_matrix*t(int_matrix)
  }
  diag(int_matrix) <- 1
  initial_pop <- c(rep(c(1, rep(0, df$strains)), df$strains),
                   rep(0, df$strains*1000 - df$strains*(df$strains+1))) |>
    matrix(nrow = 1000, byrow = T)
  return(list(initial_pop = initial_pop, interactions = int_matrix))
}
