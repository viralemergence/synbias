using CoinfectionSimulator
using LatinHypercubeSampling
using DataFrames
using Random
using Distributions
using CSV
using Tables
using Tidier
using Gadfly

interaction_scores = CSV.read("Data/interaction_scores.csv", DataFrame)
input_df = @chain CSV.read("Data/simulation_round7_input.csv", DataFrame) begin
    @mutate(sim = 1:15000)
    @left_join(interaction_scores, sim)
    @mutate(
        competition_score = total_competition/(strains^2),
        facilitation_score = total_facilitation/(strains^2),
    )
    @mutate(diff = facilitation_score - competition_score)
    @mutate(
        cf_cat = case_when(
            diff <= -0.15 => "high competition",
            diff > 0.15 => "high facilitation",
            true => missing,
        )
    )
    dropmissing()
end

high_competition = @filter(input_df, cf_cat .== "high competition")
high_facilitation = @filter(input_df, cf_cat .== "high facilitation")
competition_matrices = create_interaction_matrix(high_competition)
facilitation_matrices = create_interaction_matrix(high_facilitation)

# High Competition SEIR Simulations

n_individuals = 1000

# Collect results
n_sims = nrow(high_competition)
results_high_comp_seir = Vector(undef, n_sims)

# Simulate
for sim in 1:n_sims
    # Create SEIR strains
    n_strains = Int(high_competition.strains[sim])
    strains = Vector{SEIRModel}(undef, n_strains)
    for s in 1:n_strains
        strains[s] = SEIRModel(
            rand(Truncated(Normal(high_competition.transmission[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_competition.disease_mortality[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_competition.recovery[sim], 1), 0, 1)),
            rand(Truncated(Poisson(high_competition.latency[sim]), 1, 10))
        )
    end
    # Initialize population - everyone starts susceptible
    initial_pop = Population(Individual[])
    for p in 1:n_individuals
        individual = Individual(n_strains, 1)  # age 1
        push!(initial_pop.individuals, individual)
    end
    # Set simulation parameters
    params = SimulationParameters(
        strains,
        competition_matrices[sim],
        0.0, # base mortality
        0.0, # fecundity
        1, # Age of maturity
        :simultaneous, # introduction of strains
        100, # time steps
        :density # transmission type
    )

    true_pop = simulate(initial_pop, params)

    results_high_comp_seir[sim] = true_pop
end


agg_data = Dict{Tuple{Int,Int}, Int}()

for (i, sim_result) in enumerate(results_high_comp_seir)
    for (j, pop) in enumerate(sim_result)
        case_count = 0
        for v in pop.individuals
            case_count += sum(v[:, 3])
        end
        key = (i, j)
        agg_data[key] = get(agg_data, key, 0) + case_count
    end
end

# Convert to DataFrame
sim_timestep_df = DataFrame(
    sim = [k[1] for k in keys(agg_data)],
    timestep = [k[2] for k in keys(agg_data)],
    cases = [v for v in values(agg_data)]
)

high_comp_seir_df = @chain sim_timestep_df begin
    @group_by(timestep)
    @summarize(mean_cases = mean(cases), sd_cases = std(cases), median_cases = median(cases))
    @mutate(cases_upper = mean_cases + sd_cases, cases_lower = max(mean_cases - sd_cases, 0), type = "High competition")
end

Gadfly.plot(high_comp_seir_df,
    x = :timestep,
    y = :mean_cases,
    ymin = :cases_lower,
    ymax = :cases_upper,
    Geom.line(),
    Geom.ribbon(),
    Guide.xlabel("Timestep"),
    Guide.ylabel("Total Cases"),
    Guide.title("Average High-Competition SEIR Simulation")
)

# High Competition SI Simulations

# Collect results
n_sims = nrow(high_competition)
results_high_comp_si = Vector(undef, n_sims)

# Simulate
for sim in 1:n_sims
    # Create SI strains
    n_strains = Int(high_competition.strains[sim])
    strains = Vector{SIModel}(undef, n_strains)
    for s in 1:n_strains
        strains[s] = SIModel(
            rand(Truncated(Normal(high_competition.transmission[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_competition.disease_mortality[sim], 1), 0, 1))
        )
    end
    # Initialize population - everyone starts susceptible
    initial_pop = Population(Individual[])
    for p in 1:n_individuals
        individual = Individual(n_strains, 1)  # age 1
        push!(initial_pop.individuals, individual)
    end
    # Set simulation parameters
    params = SimulationParameters(
        strains,
        competition_matrices[sim],
        0.0, # base mortality
        0.0, # fecundity
        1, # Age of maturity
        :simultaneous, # introduction of strains
        100, # time steps
        :frequency # transmission type
    )

    true_pop = simulate(initial_pop, params)

    results_high_comp_si[sim] = true_pop
end

agg_data_si = Dict{Tuple{Int,Int}, Int}()

for (i, sim_result) in enumerate(results_high_comp_si)
    for (j, pop) in enumerate(sim_result)
        case_count = 0
        for v in pop.individuals
            case_count += sum(v[:, 3])
        end
        key = (i, j)
        agg_data_si[key] = get(agg_data_si, key, 0) + case_count
    end
end

# Convert to DataFrame
sim_timestep_df_si = DataFrame(
    sim = [k[1] for k in keys(agg_data_si)],
    timestep = [k[2] for k in keys(agg_data_si)],
    cases = [v for v in values(agg_data_si)]
)

high_comp_si_df = @chain sim_timestep_df_si begin
    @group_by(timestep)
    @summarize(mean_cases = mean(cases), sd_cases = std(cases), median_cases = median(cases))
    @mutate(cases_upper = mean_cases + sd_cases, cases_lower = max(mean_cases - sd_cases, 0), type = "High competition")
end

Gadfly.plot(high_comp_si_df,
    x = :timestep,
    y = :mean_cases,
    ymin = :cases_lower,
    ymax = :cases_upper,
    Geom.line(),
    Geom.ribbon(),
    Guide.xlabel("Timestep"),
    Guide.ylabel("Total Cases"),
    Guide.title("Average High-Competition SI Simulation")
)

# High Facilitation SEIR Simulations

# Collect results
n_sims = nrow(high_facilitation)
results_high_fac_seir = Vector(undef, n_sims)

# Simulate
for sim in 1:n_sims
    # Create SEIR strains
    n_strains = Int(high_facilitation.strains[sim])
    strains = Vector{SEIRModel}(undef, n_strains)
    for s in 1:n_strains
        strains[s] = SEIRModel(
            rand(Truncated(Normal(high_facilitation.transmission[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_facilitation.disease_mortality[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_facilitation.recovery[sim], 1), 0, 1)),
            rand(Truncated(Poisson(high_facilitation.latency[sim]), 1, 10))
        )
    end
    # Initialize population - everyone starts susceptible
    initial_pop = Population(Individual[])
    for p in 1:n_individuals
        individual = Individual(n_strains, 1)  # age 1
        push!(initial_pop.individuals, individual)
    end
    # Set simulation parameters
    params = SimulationParameters(
        strains,
        facilitation_matrices[sim],
        0.0, # base mortality
        0.0, # fecundity
        1, # Age of maturity
        :simultaneous, # introduction of strains
        100, # time steps
        :density # transmission type
    )
    true_pop = simulate(initial_pop, params)
    results_high_fac_seir[sim] = true_pop
end

agg_data_fac = Dict{Tuple{Int,Int}, Int}()

for (i, sim_result) in enumerate(results_high_fac_seir)
    for (j, pop) in enumerate(sim_result)
        case_count = 0
        for v in pop.individuals
            case_count += sum(v[:, 3])
        end
        key = (i, j)
        agg_data_fac[key] = get(agg_data_fac, key, 0) + case_count
    end
end

# Convert to DataFrame
sim_timestep_df_fac = DataFrame(
    sim = [k[1] for k in keys(agg_data_fac)],
    timestep = [k[2] for k in keys(agg_data_fac)],
    cases = [v for v in values(agg_data_fac)]
)

high_fac_seir_df = @chain sim_timestep_df_fac begin
    @group_by(timestep)
    @summarize(mean_cases = mean(cases), sd_cases = std(cases), median_cases = median(cases))
    @mutate(cases_upper = mean_cases + sd_cases, cases_lower = max(mean_cases - sd_cases, 0), type = "High facilitation")
    @bind_rows(high_comp_seir_df)
end

p1 = Gadfly.plot(high_fac_seir_df,
    x = :timestep,
    y = :mean_cases,
    ymin = :cases_lower,
    ymax = :cases_upper,
    color = :type,
    Geom.line(),
    Geom.ribbon(),
    Scale.color_discrete_manual("#004488", "#BB5566"),
    Theme(alphas=[0.5], line_width = 2mm),
    Guide.xlabel("Timestep"),
    Guide.ylabel("Total Cases"),
    Guide.title("SEIR Simulations: The most extreme interaction scenarios")
)

# High Facilitation SI Simulations

# Collect results
n_sims = nrow(high_facilitation)
results_high_fac_si = Vector(undef, n_sims)
# Simulate
for sim in 1:n_sims
    # Create SI strains
    n_strains = Int(high_facilitation.strains[sim])
    strains = Vector{SIModel}(undef, n_strains)
    for s in 1:n_strains
        strains[s] = SIModel(
            rand(Truncated(Normal(high_facilitation.transmission[sim], 1), 0, 1)),
            rand(Truncated(Normal(high_facilitation.disease_mortality[sim], 1), 0, 1))
        )
    end
    # Initialize population - everyone starts susceptible
    initial_pop = Population(Individual[])
    for p in 1:n_individuals
        individual = Individual(n_strains, 1)  # age 1
        push!(initial_pop.individuals, individual)
    end
    # Set simulation parameters
    params = SimulationParameters(
        strains,
        facilitation_matrices[sim],
        0.0, # base mortality
        0.0, # fecundity
        1, # Age of maturity
        :simultaneous, # introduction of strains
        100, # time steps
        :frequency # transmission type
    )
    true_pop = simulate(initial_pop, params)
    results_high_fac_si[sim] = true_pop
end

agg_data_fac_si = Dict{Tuple{Int,Int}, Int}()

for (i, sim_result) in enumerate(results_high_fac_si)
    for (j, pop) in enumerate(sim_result)
        case_count = 0
        for v in pop.individuals
            case_count += sum(v[:, 3])
        end
        key = (i, j)
        agg_data_fac_si[key] = get(agg_data_fac_si, key, 0) + case_count
    end
end

# Convert to DataFrame
sim_timestep_df_fac_si = DataFrame(
    sim = [k[1] for k in keys(agg_data_fac_si)],
    timestep = [k[2] for k in keys(agg_data_fac_si)],
    cases = [v for v in values(agg_data_fac_si)]
)

high_fac_si_df = @chain sim_timestep_df_fac_si begin
    @group_by(timestep)
    @summarize(mean_cases = mean(cases), sd_cases = std(cases), median_cases = median(cases))
    @mutate(cases_upper = mean_cases + sd_cases, cases_lower = max(mean_cases - sd_cases, 0), type = "High facilitation")
    @bind_rows(high_comp_si_df)
end

p2 = Gadfly.plot(high_fac_si_df,
    x = :timestep,
    y = :mean_cases,
    ymin = :cases_lower,
    ymax = :cases_upper,
    color = :type,
    Coord.cartesian(xmin=0, xmax=5),
    Geom.line(),
    Geom.ribbon(),
    Scale.color_discrete_manual("#004488", "#BB5566"),
    Theme(alphas=[0.5], line_width = 2mm),
    Guide.xlabel("Timestep"),
    Guide.ylabel("Total Cases"),
    Guide.title("SI Simulations: The most extreme interaction scenarios")
)

p1 |> PNG("Figures/SuppFig2_SEIR_example_sims.png", 6inch, 4.5inch)
p2 |> PNG("Figures/SuppFig2_SI_example_sims.png", 6inch, 4.5inch)