using Test
using DataFrames
using Distributions
using Random

include("coinfection_data_prep.jl")

@testset "prep_interaction_matrix tests" begin
    
    @testset "Missing columns error" begin
        # Test missing interaction_strength column
        df_missing = DataFrame(
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [3]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_missing)
        
        # Test missing all required columns
        df_empty = DataFrame()
        @test_throws ErrorException prep_interaction_matrix(df=df_empty)
    end
    
    @testset "Invalid strains count" begin
        df_invalid = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [0]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_invalid)
        
        df_negative = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [-1]
        )
        @test_throws ErrorException prep_interaction_matrix(df=df_negative)
    end
    
    @testset "Valid input - basic functionality" begin
        Random.seed!(123)
        df = DataFrame(
            interaction_strength = [0.1, 0.2],
            cf_ratio = [0.5, 0.8],
            priority_effects = [true, false],
            strains = [2, 3]
        )
        
        result = prep_interaction_matrix(df=df)
        
        @test length(result) == 2
        @test size(result[1]) == (2, 2)
        @test size(result[2]) == (3, 3)
        @test all(isa.(result, Matrix{Float64}))
    end
    
    @testset "Priority effects - asymmetric matrix" begin
        Random.seed!(456)
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [3]
        )
        
        result = prep_interaction_matrix(df=df)
        matrix = result[1]
        
        # Diagonal should be 1.0
        @test all(diag(matrix) .== 1.0)
        
        # Off-diagonal elements should not necessarily be symmetric
        @test matrix[1,2] != matrix[2,1] || matrix[1,3] != matrix[3,1]
    end
    
    @testset "No priority effects - symmetric matrix" begin
        Random.seed!(789)
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [false],
            strains = [3]
        )
        
        result = prep_interaction_matrix(df=df)
        matrix = result[1]
        
        # Diagonal should be 1.0
        @test all(diag(matrix) .== 1.0)
        
        # Matrix should be symmetric
        @test matrix == matrix'
    end
    
    @testset "Single strain matrix" begin
        df = DataFrame(
            interaction_strength = [0.1],
            cf_ratio = [0.5],
            priority_effects = [true],
            strains = [1]
        )
        
        result = prep_interaction_matrix(df=df)
        
        @test size(result[1]) == (1, 1)
        @test result[1][1,1] == 1.0
    end
end