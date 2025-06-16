"""
    prep_interaction_matrix(; df::DataFrame)

Generates an interaction matrix for each row in the provided DataFrame `df`, based on the specified parameters.

# Arguments
- `df::DataFrame`: A DataFrame containing the following required columns:
    - `:interaction_strength` (`Float64`): Standard deviation for the Gaussian distribution used to generate interaction values.
    - `:cf_ratio` (`Float64`): Mean for the Gaussian distribution used to generate interaction values.
    - `:priority_effects` (`Bool`): If `true`, the interaction matrix is asymmetric; otherwise, it is symmetric.
    - `:strains` (`Int`): The number of strains, which determines the size of the interaction matrix.

# Returns
- `int_matrix_list::Vector{Matrix{Float64}}`: A list of interaction matrices, where each matrix corresponds to a row in the DataFrame. Each matrix is of size `strains × strains`.

# Errors
- Throws an error if any required column is missing from the DataFrame.
- Throws an error if the number of strains is not positive.

# Notes
- For each row in the DataFrame, an interaction matrix of size `strains × strains` is created.
- If `priority_effects` is `true`, off-diagonal elements are filled with random values from a Gaussian distribution with mean `cf_ratio` and standard deviation `interaction_strength`.
- If `priority_effects` is `false`, the matrix is filled symmetrically.
"""
using DataFrames
using Distributions
using Random

function prep_interaction_matrix(;
    df::DataFrame
)
    # Ensure the DataFrame has the required columns
    required_columns = [:interaction_strength, :cf_ratio, :priority_effects, :strains]
    for col in required_columns
        if !haskey(df, col)
            error("DataFrame must contain column: $col")
        end
    end

    # Ensure columns are of the correct type
    df.interaction_strength = convert(Vector{Float64}, df.interaction_strength)
    df.cf_ratio = convert(Vector{Float64}, df.cf_ratio)
    df.priority_effects = convert(Vector{Bool}, df.priority_effects)
    df.strains = convert(Vector{Int}, df.strains)

    # Collect matrices in a vector
    int_matrix_list = Vector{Matrix{Float64}}(undef, nrow(df))

    for row in eachrow(df)
        # Create matrix of the correct size
        if row.strains <= 0
            error("Number of strains must be positive, got: $(row.strains)")
        end
        int_matrix = ones(Float64, row.strains, row.strains)
        # Fill the interaction matrix
        for i in 1:row.strains
            for j in 1:row.strains
                if row.priority_effects
                    if i != j
                        int_matrix[i, j] = rand(Gaussian(row.cf_ratio, row.interaction_strength))
                    end
                else
                    if i < j
                        int_matrix[i, j] = rand(Gaussian(row.cf_ratio, row.interaction_strength))
                    end
                    if i > j
                        int_matrix[i, j] = int_matrix[j, i]
                    end
                end
            end
        end

        # Store the interaction matrix in the list
        int_matrix_list[row.index] = int_matrix
    end

    # Return the list of interaction matrices
    return int_matrix_list

end