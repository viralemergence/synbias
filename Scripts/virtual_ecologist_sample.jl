virtual_ecologist_sample = function (;
	virtual_population::Vector{Vector{Matrix{Bool}}},
	proportion_sampled::Float64,
	false_positive_rate::Float64,
	false_negative_rate::Float64,
)
	# Assertions
	@assert length(virtual_population) > 0 "The virtual population must not be empty."
	@assert all(t -> all(m -> size(m, 2) == 4, t), virtual_population) "Each matrix in the virtual population must have 4 columns."
	@assert all(t -> all(m -> length(unique([size(m, 1) for m in t])) == 1, t), virtual_population) "All matrices in the virtual population must have the same number of rows."
	@assert proportion_sampled >= 0 && proportion_sampled <= 1 "Proportion sampled must be between 0 and 1."
	@assert false_positive_rate >= 0 && false_positive_rate <= 1 "False positive rate must be between 0 and 1."
	@assert false_negative_rate >= 0 && false_negative_rate <= 1 "False negative rate must be between 0 and 1."

	n_viruses = size(virtual_population[1][1], 1)
	timesteps = length(virtual_population)
	individuals_sampled = round(Int, proportion_sampled * size(virtual_population, 1))
	detect_matrix = falses(timesteps, n_viruses)

	for i in 1:timesteps
		# Sample individuals
		sampled_individuals = randperm(size(virtual_population[i], 1))[1:individuals_sampled]
		perfect_sample = [virtual_population[i][j][:, 3] for j in sampled_individuals]
		imperfect_sample = falses(individuals_sampled, n_viruses)

		# Apply false positive and false negative rates
		for i in 1:individuals_sampled
			for j in 1:n_viruses
				if perfect_sample[i][j] == true && false_negative_rate > 0
					imperfect_sample[i, j] = rand() < (1 - false_negative_rate) ? true : false
				else
					if false_positive_rate > 0
						imperfect_sample[i, j] = rand() < false_positive_rate ? true : false
					end
				end
			end
		end

		# Store result in detect_matrix
		detect_matrix[i, :] = sum(imperfect_sample, dims = 1) .> 0

	end

	# Return the detection matrix
	return detect_matrix
end
