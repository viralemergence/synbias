"""
  coinfection_simulator(;
	initial_pop::Vector{Matrix{Bool}},
	ages::Vector{Int},
	interactions::Matrix{Float64},
	disease_type::Vector{String},
	base_mortality::Float64,
	disease_mortality::Vector{Float64},
	fecundity::Float64,
	transmission::Vector{Float64},
	time_steps::Int,
	age_maturity::Int,
	introduction::String = "simultaneous",
	latency::Union{Vector{Int}, Nothing} = nothing,
	recovery::Union{Vector{Float64}, Nothing} = nothing,
	immunity_loss::Union{Vector{Float64}, Nothing} = nothing
  ) -> Tuple{Vector{Vector{Matrix{Bool}}}, Vector{Vector{Int}}}

Simulates multiple pathogen strains spreading through a host population with possible coinfection dynamics.

# Arguments
- `initial_pop::Vector{Matrix{Bool}}`: Vector of matrices representing the initial population. 
  Each matrix represents an individual, with rows corresponding to strains and columns 
  representing disease states [S,E,I,R] (Susceptible, Exposed, Infected, Recovered).
- `ages::Vector{Int}`: Ages of each individual in the initial population.
- `interactions::Matrix{Float64}}`: Matrix of interaction factors between strains. Values > 1 indicate 
  synergistic interactions, values < 1 indicate antagonistic interactions.
- `disease_type::Vector{String}`: Vector specifying model type for each strain. 
  Must be one of: "si", "sir", "seir", or "seirs".
- `base_mortality::Float64`: Background mortality rate per time step.
- `disease_mortality::Vector{Float64}`: Additional mortality rate for each strain when infected.
- `fecundity::Float64`: Mean number of offspring per mature individual per time step.
- `transmission::Vector{Float64}`: Transmission probability for each strain.
- `time_steps::Int`: Number of time steps to simulate.
- `age_maturity::Int`: Age at which individuals can reproduce.
- `introduction::String`: How strains are introduced: "simultaneous" (all at once) or 
  "random" (randomly throughout simulation) or "none" (no infections are introduced by the function). 
  Default is "simultaneous".
- `latency::Union{Vector{Int}, Nothing}`: Required for SEIR/SEIRS models. Vector of latency periods 
  for each strain (number of time steps in exposed state).
- `recovery::Union{Vector{Float64}, Nothing}`: Required for SIR/SEIR/SEIRS models. 
  Recovery probability for each strain.
- `immunity_loss::Union{Vector{Float64}, Nothing}`: Required for SEIRS models.
  Probability of losing immunity for each strain.

# Returns
A tuple containing:
1. Vector of population states at each time step, where each state is a vector of individual matrices
2. Vector of individual ages at each time step (in the same order as the population states)

# Disease State Format
Each individual is represented by an n×4 matrix where n is the number of strains:
- Column 1: Susceptible state (true if susceptible)
- Column 2: Exposed state (true if exposed, for SEIR/SEIRS models)
- Column 3: Infected state (true if infected)
- Column 4: Recovered state (true if recovered, for SIR/SEIR/SEIRS models)

"""
# cSpell:ignore susceptibles infecteds indiv seir seirs juvs coinfection
using Random, Distributions, StatsBase

function coinfection_simulator(;
	initial_pop::Vector{Matrix{Bool}},
	ages::Vector{Int},
	interactions::Matrix{Float64},
	disease_type::Vector{String},
	base_mortality::Float64,
	disease_mortality::Vector{Float64},
	fecundity::Float64,
	transmission::Vector{Float64},
	time_steps::Int,
	age_maturity::Int,
	introduction::String = "simultaneous",
	latency::Union{Vector{Int}, Nothing} = nothing,
	immunity_loss::Union{Vector{Float64}, Nothing} = nothing,
	recovery::Union{Vector{Float64}, Nothing} = nothing,
)
	# Input validation
	@assert all(m -> size(m, 2) == 4, initial_pop) "All population matrices must have 4 columns"
	@assert length(unique([size(m, 1) for m in initial_pop])) == 1 "All population matrices must have the same number of rows"
	@assert all(a -> a >= 0, ages) "All ages must be 0 or greater"

	n_strains = size(initial_pop[1], 1)
	n_individuals = length(initial_pop)

	@assert length(ages) == n_individuals "There must be an age given for every individual"
	@assert size(interactions, 1) == n_strains && size(interactions, 2) == n_strains "The interactions matrix must have rows and columns equal to the number of strains"
	@assert all(dt -> dt in ["si", "sir", "seir", "seirs"], disease_type) "All disease types must be one of: si, sir, seir, seirs"
	@assert length(disease_type) == n_strains "Disease type vector length must match number of strains"
	@assert base_mortality >= 0 && base_mortality <= 1 "Base mortality must fall between 0 and 1"
	@assert length(disease_mortality) == n_strains "Disease mortality vector length must match number of strains"
	@assert all(m -> m >= 0 && m <= 1, disease_mortality) "Disease mortality must fall between 0 and 1"
	@assert all(base_mortality .+ disease_mortality .<= 1) "Base mortality and disease mortality combined must not exceed 1"
	@assert fecundity >= 0 "Fecundity cannot be negative"
	@assert length(transmission) == n_strains "Transmission vector length must equal number of strains"
	@assert all(t -> t >= 0 && t <= 1, transmission) "Transmission must fall between 0 and 1"

	if any(dt -> dt in ["sir", "seir", "seirs"], disease_type)
		@assert !isnothing(recovery) "Recovery rates must be provided"
		@assert length(recovery) == n_strains "Recovery vector length must equal the number of strains"
		@assert all(r -> r >= 0 && r <= 1, recovery) "All recovery values must fall between 0 and 1"

		if any(dt -> dt in ["seir", "seirs"], disease_type)
			@assert !isnothing(latency) "Latency periods must be provided"
			@assert length(latency) == n_strains "Latency vector length must equal the number of strains"
			@assert all(l -> l > 0, latency) "Latency in days must be greater than 0"

			if any(dt -> dt == "seirs", disease_type)
				@assert !isnothing(immunity_loss) "Immunity loss rates must be provided"
				@assert length(immunity_loss) == n_strains "Immunity loss vector length must equal number of strains"
				@assert all(il -> il >= 0 && il <= 1, immunity_loss) "All immunity loss values must be between 0 and 1"
			end
		end
	end

	@assert time_steps >= 1 "Time steps must be 1 or greater"
	@assert introduction in ["simultaneous", "random", "none"] "Introduction of viruses must be none, simultaneous or random"
	@assert age_maturity > 0 "Age of maturity must be greater than 0"

	# Strain introduction timing
	if introduction == "simultaneous"
		intro_step = ones(Int, n_strains)
	elseif introduction == "random"
		intro_step = rand(1:time_steps, n_strains)
	elseif introduction == "none"
		intro_step = zeros(Int, n_strains)
	else
		error("Invalid introduction type: $introduction")
	end

	# Initialize results
	result_pop = Vector{Vector{Matrix{Bool}}}(undef, time_steps)
	result_pop[1] = copy.(initial_pop)
	result_ages = Vector{Vector{Int}}(undef, time_steps)
	result_ages[1] = copy(ages)

	# Time loop
	for t in 1:(time_steps-1)
		current_pop = deepcopy(result_pop[t])
		current_ages = deepcopy(result_ages[t])

		# Introduce infections
		if any(intro_step .== t)
			infected_indices = sample(1:length(current_pop), sum(intro_step .== t), replace = false)
			for (i, idx) in enumerate(infected_indices)
				strain = findall(intro_step .== t)[i]
				current_pop[idx][strain, 1] = false
				current_pop[idx][strain, 3] = true
			end
		end

		# Breeding
		breeding_age = current_ages .>= age_maturity
		number_births = rand(Poisson(sum(breeding_age) * fecundity))
		if number_births > 0
			new_juvs = [falses(n_strains, 4) for _ in 1:number_births]
			for m in new_juvs
				m[:, 1] .= true
			end
			append!(current_pop, new_juvs)
			append!(current_ages, zeros(Int, number_births))
		end

		# For each strain, detect whether there are any infected/exposed
		# This lets us save time computationally when no infection is active
		# For each strain, detect whether there are any infected/exposed
		active_infections = if isempty(current_pop)
			falses(n_strains)
		else
			infection_status = map(m -> m[:, [2, 3]], current_pop)
			sum_matrix = reduce(hcat, infection_status)
			sum_vec = vec(sum(sum_matrix, dims = 2))
			sum_vec .> 0
		end

		# Death tracker
		dead_indices = Set{Int}()

		# Process each strain
		for strain in 1:n_strains
			# Death of uninfected individuals
			uninfected = findall(m -> any(view(m, strain, [1, 2, 4])), current_pop)
			uninfected_death = isempty(uninfected) ? 0 : rand(Binomial(length(uninfected), base_mortality))
			if uninfected_death > 0
				dead = sample(uninfected, uninfected_death; replace = false)
				union!(dead_indices, dead)
			end

			# Skip if no active infections for this strain
			!active_infections[strain] && continue

			# Get alive individuals
			alive = [i ∉ dead_indices for i in 1:length(current_pop)]

			# Process based on disease type
			if disease_type[strain] == "si"
				handle_si_disease(
					current_pop, alive, strain, transmission[strain], interactions[strain, :],
					base_mortality, disease_mortality[strain], dead_indices,
				)
			elseif disease_type[strain] == "sir"
				handle_sir_disease(
					current_pop, alive, strain, transmission[strain], interactions[strain, :],
					base_mortality, disease_mortality[strain], recovery[strain], dead_indices,
				)
			elseif disease_type[strain] == "seir"
				handle_seir_disease(
					current_pop, alive, strain, transmission[strain], interactions[strain, :],
					base_mortality, disease_mortality[strain], recovery[strain], latency[strain], dead_indices,
				)
			elseif disease_type[strain] == "seirs"
				handle_seirs_disease(
					current_pop, alive, strain, transmission[strain], interactions[strain, :],
					base_mortality, disease_mortality[strain], recovery[strain], latency[strain],
					immunity_loss[strain], dead_indices,
				)
			end
		end

		# Final cleanup after all mortality:
		alive_mask = [i ∉ dead_indices for i in 1:length(current_pop)]
		current_pop = current_pop[alive_mask]
		current_ages = current_ages[alive_mask]

		# Age surviving individuals
		current_ages .+= 1
		result_pop[t+1] = current_pop
		result_ages[t+1] = current_ages
	end

	return (result_pop, result_ages)
end

# Helper functions
function infect(
	susceptibles::Vector{Matrix{Bool}},
	infecteds::Vector{Matrix{Bool}},
	interactions::Vector{Float64},
	beta::Float64,
	strain::Int,
)
	# Convert susceptibles to individual/strain matrix
	strain_matrix = map(m -> sum.(eachrow(m[:, [2, 3]])), susceptibles)
	strain_matrix = vcat(strain_matrix'...)

	indiv = [strain_matrix[i, :] for i in 1:size(strain_matrix, 1)]
	total_strains = sum.(indiv)

	# Infect individuals
	indiv_new = Vector{Int}(undef, length(indiv))
	for i in 1:length(indiv)
		indiv_prev = indiv[i]
		int_indiv = interactions .* indiv_prev

		if total_strains[i] > 1
			prob = prod(int_indiv) * beta
			indiv_prev[strain] = max(indiv_prev[strain], rand(Binomial(length(infecteds), prob)))
		elseif total_strains[i] == 1
			prob = sum(int_indiv) * beta
			indiv_prev[strain] = max(indiv_prev[strain], rand(Binomial(length(infecteds), prob)))
		else
			indiv_prev[strain] = rand(Binomial(length(infecteds), beta))
		end

		indiv_new[i] = indiv_prev[strain]
	end

	return indiv_new .> 0
end

function handle_si_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)
end

function handle_sir_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	dead_indices::Set{Int},
)
	# Infection
	handle_infection(current_pop, alive, strain, transmission, interactions)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

function handle_seir_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	dead_indices::Set{Int},
)
	# Exposure
	handle_exposure(current_pop, alive, strain, transmission, interactions)

	# Infection from exposed
	handle_exposed_infection(current_pop, alive, strain, latency)

	# Death of infected
	handle_infected_death(current_pop, strain, base_mortality, disease_mortality, dead_indices)

	# Recovery
	handle_recovery(current_pop, alive, strain, recovery)
end

function handle_seirs_disease(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
	base_mortality::Float64,
	disease_mortality::Float64,
	recovery::Float64,
	latency::Int,
	immunity_loss::Float64,
	dead_indices::Set{Int},
)
	# SEIR processes
	handle_seir_disease(current_pop, alive, strain, transmission, interactions,
		base_mortality, disease_mortality, recovery, latency, dead_indices)

	# Loss of immunity
	handle_immunity_loss(current_pop, alive, strain, immunity_loss)
end

# Utility functions
function handle_infection(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible = findall(m -> m[strain, 1], current_pop[alive])
	infected = findall(m -> m[strain, 3], current_pop[alive])

	if !isempty(susceptible)
		infection = infect(current_pop[alive][susceptible],
			current_pop[alive][infected],
			interactions,
			transmission, strain)
		alive_idx = findall(identity, alive)
		for idx in alive_idx[susceptible[infection]]
			current_pop[idx][strain, 1] = false
			current_pop[idx][strain, 3] = true
		end
	end
end

function handle_exposure(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	transmission::Float64,
	interactions::Vector{Float64},
)
	susceptible = findall(m -> m[strain, 1], current_pop[alive])
	infected = findall(m -> m[strain, 3], current_pop[alive])

	if !isempty(susceptible)
		infection = infect(current_pop[alive][susceptible],
			current_pop[alive][infected],
			interactions,
			transmission, strain)
		alive_idx = findall(identity, alive)
		for idx in alive_idx[susceptible[infection]]
			current_pop[idx][strain, 1] = false
			current_pop[idx][strain, 2] = true
		end
	end
end

function handle_exposed_infection(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	latency::Int,
)
	exposed = findall(m -> m[strain, 2], current_pop[alive])
	if !isempty(exposed)
		infection = sample(exposed, rand(Binomial(length(exposed), 1 / latency)); replace = false)
		alive_idx = findall(identity, alive)
		for idx in alive_idx[infection]
			current_pop[idx][strain, 2] = false
			current_pop[idx][strain, 3] = true
		end
	end
end

function handle_recovery(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	recovery::Float64,
)
	infected = findall(m -> m[strain, 3], current_pop[alive])
	if !isempty(infected)
		recovered = sample(infected, rand(Binomial(length(infected), recovery)); replace = false)
		alive_idx = findall(identity, alive)
		for idx in alive_idx[recovered]
			current_pop[idx][strain, 3] = false
			current_pop[idx][strain, 4] = true
		end
	end
end

function handle_immunity_loss(
	current_pop::Vector{Matrix{Bool}},
	alive::Vector{Bool},
	strain::Int,
	immunity_loss::Float64,
)
	recovered = findall(m -> m[strain, 4], current_pop[alive])
	if !isempty(recovered)
		susceptible = sample(recovered, rand(Binomial(length(recovered), immunity_loss)); replace = false)
		alive_idx = findall(identity, alive)
		for idx in alive_idx[susceptible]
			current_pop[idx][strain, 4] = false
			current_pop[idx][strain, 1] = true
		end
	end
end

function handle_infected_death(
	current_pop::Vector{Matrix{Bool}},
	strain::Int,
	base_mortality::Float64,
	disease_mortality::Float64,
	dead_indices::Set{Int},
)
	infected = findall(m -> m[strain, 3], current_pop)
	infected_death = isempty(infected) ? 0 : rand(Binomial(length(infected),
		base_mortality + disease_mortality))
	if infected_death > 0
		dead = sample(infected, infected_death; replace = false)
		union!(dead_indices, dead)
	end
end
