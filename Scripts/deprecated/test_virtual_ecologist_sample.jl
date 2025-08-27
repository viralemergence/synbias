using Test
include("virtual_ecologist_sample.jl")
include("coinfection_simulator.jl")

@testset "virtual_ecologist_sample tests" begin
	@testset "Basic functionality" begin
		# Create a simple virtual population
		# 2 timesteps, 2 individuals per timestep, 3 viruses, 4 columns per matrix
		timesteps = 2
		individuals = 2
		n_viruses = 3

		# Create a virtual population where first individual has virus 1, second has virus 2
		virtual_pop = Vector{Vector{Matrix{Bool}}}(undef, timesteps)
		for t in 1:timesteps
			virtual_pop[t] = Vector{Matrix{Bool}}(undef, individuals)
			for i in 1:individuals
				virtual_pop[t][i] = falses(n_viruses, 4)
				virtual_pop[t][i][i, 3] = true
				# Set all rows except i to true in column 1
				rows_except_i = setdiff(1:n_viruses, i)
				virtual_pop[t][i][rows_except_i, 1] .= true
			end
		end

		# Test with no errors in detection
		result = virtual_ecologist_sample(
			virtual_population = virtual_pop,
			proportion_sampled = 1.0,
			false_positive_rate = 0.0,
			false_negative_rate = 0.0,
		)

		# Since we sample all individuals with no errors, we should detect all viruses
		@test all(result[1, 1:2] .== true)
		@test all(result[2, 1:2] .== true)
		@test all(result[:, 3] .== false) # Virus 3 is never present

		# Test with no sampling
		result_no_sample = virtual_ecologist_sample(
			virtual_population = virtual_pop,
			proportion_sampled = 0.0,
			false_positive_rate = 0.0,
			false_negative_rate = 0.0,
		)

		# With 50% sampling, we should detect at most one virus per timestep
		@test sum(result_no_sample[1, :]) == 0
		@test sum(result_no_sample[2, :]) == 0
	end

	@testset "False positive and negative rates" begin
		# Create a test population with no infected individuals
		virtual_pop = [[Matrix{Bool}([true  false false false;]) for _ in 1:100] for _ in 1:5]

		# Test false positives
		result_fp = virtual_ecologist_sample(
			virtual_population = virtual_pop,
			proportion_sampled = 1.0,
			false_positive_rate = 1.0,  # Always give false positives
			false_negative_rate = 0.0,
		)

		# Should detect all viruses due to false positives
		@test all(result_fp .== true)

		# Create a test population with all individuals infected
		virtual_pop_all = [[Matrix{Bool}([false false true false;]) for _ in 1:100] for _ in 1:5]

		# Test false negatives
		result_fn = virtual_ecologist_sample(
			virtual_population = virtual_pop_all,
			proportion_sampled = 1.0,
			false_positive_rate = 0.0,
			false_negative_rate = 1.0, # Always miss viruses
		)

		# Should miss all viruses due to false negatives
		@test all(result_fn .== false)
	end

	@testset "Input validation" begin
		valid_pop = [[Matrix{Bool}([false false true false;]) for _ in 1:100] for _ in 1:5]

		# Test empty population
		@test_throws AssertionError virtual_ecologist_sample(
			virtual_population = Vector{Vector{Matrix{Bool}}}(),
			proportion_sampled = 0.5,
			false_positive_rate = 0.0,
			false_negative_rate = 0.0,
		)

		# Test invalid proportion
		@test_throws AssertionError virtual_ecologist_sample(
			virtual_population = valid_pop,
			proportion_sampled = -0.1,
			false_positive_rate = 0.0,
			false_negative_rate = 0.0,
		)

		@test_throws AssertionError virtual_ecologist_sample(
			virtual_population = valid_pop,
			proportion_sampled = 1.1,
			false_positive_rate = 0.0,
			false_negative_rate = 0.0,
		)

		# Test invalid error rates
		@test_throws AssertionError virtual_ecologist_sample(
			virtual_population = valid_pop,
			proportion_sampled = 0.5,
			false_positive_rate = -0.1,
			false_negative_rate = 0.0,
		)

		@test_throws AssertionError virtual_ecologist_sample(
			virtual_population = valid_pop,
			proportion_sampled = 0.5,
			false_positive_rate = 0.0,
			false_negative_rate = 1.1,
		)
	end
end
