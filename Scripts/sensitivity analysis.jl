using CSV, DataFrames, MLJ, Gadfly, Cairo, Fontconfig
#### Sensitivity analysis: random timestep data ####
# Import data
rt_dat = vcat(CSV.read("Data/simulation_round1.csv", DataFrame) |> 
                x -> insertcols!(x, :sample_prop => 0.1),
              CSV.read("Data/simulation_round2.csv", DataFrame) |> 
                x -> insertcols!(x, :sample_prop => 0.1),
              CSV.read("Data/simulation_round3.csv", DataFrame))
first(rt_dat)
# Prep data
y, X = unpack(rt_dat, ==(:richness_count), colname -> colname != :id, rng=18)
train, test = partition(eachindex(y), 0.7, shuffle=true, rng=18)
X_train, y_train = X[train, :], y[train]
X_test, y_test = X[test, :], y[test]
# Tune the neural net
Nnet = @load NeuralNetworkRegressor pkg=BetaML
nnet = Nnet()
r1 = range(nnet, :batch_size, lower=8, upper=256)
r2 = range(nnet, :epochs, lower=10, upper=200)
tuned_nnet = TunedModel(
  model = nnet,
  tuning = RandomSearch(),
  resampling = CV(nfolds = 6),
  range = [r1, r2],
  measure = rms,
  n = 25
)
mach = machine(tuned_nnet, X_train, y_train) # epochs = 41, batch_size = 22
# Fit and predict on test data
fit!(mach)
yhat = predict(mach, X[test,:])
mape(yhat, y[test])
# Plot feature importance
function permutation_importance(mach, X, y, measure)
  yhat = predict(mach, X)
  baseline = mean(measure(yhat, y))
  importances = zeros(size(X, 2))
  for i in 1:size(X, 2)
      X_permuted = copy(X)
      X_permuted[:, i] = X[shuffle(1:end), i]
      yhat_permuted = predict(mach, X_permuted)
      importances[i] = mean(measure(yhat_permuted, y)) - baseline
  end
  return importances
end
importances = permutation_importance(mach, X_test, y_test, mape)
x_label = ["Interaction\nStrength", "Priority\nEffects", "C:F Ratio", "# Strains", "Proportion\nSampled"]
p = Gadfly.plot(x = x_label, y = importances, color = x_label, Geom.bar, 
                Guide.xlabel("Simulation Settings"),
                Guide.ylabel("Permutation Importance"), 
                Scale.color_discrete_manual("#B25D91", "#CB87B4", "#EFC7E6", "#1BB6AF", "#088BBE"),
                Coord.cartesian(ymin = -0.05, ymax = 1),
                Theme(key_position = :none))
Gadfly.draw(Gadfly.PNG("Figures/sensitivity_analysis_rt.png"), p)
#### Sensitivity analysis: 35th timestep data ####
dat_35t = CSV.read("Data/simulation_round5/simulation_round5.csv", DataFrame)
first(dat_35t)