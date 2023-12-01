library(poems)
synbias_lhs <- LatinHypercubeSampler$new()
synbias_lhs$set_uniform_parameter("interaction_strength", 0, 0.75)
synbias_lhs$set_class_parameter("priority_effects", c(TRUE, FALSE))
synbias_lhs$set_uniform_parameter("cf_ratio")
synbias_lhs$set_uniform_parameter("strains", 10, 100, decimals = 0)
sample_data <- synbias_lhs$generate_samples(10000, random_seed = 324) %>%
  rowid_to_column("id")
# add code here to generate interaction matrices based on sample data
