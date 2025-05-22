using Test, Random

# Set random seed for reproducible tests
Random.seed!(123)

@testset "Coinfection Simulator Tests" begin

    @testset "infect() Function" begin

        # Test 1: Full susceptibility with 100% transmission
        susceptibles = [Matrix{Bool}(falses(2,4)) for _ in 1:10]
        interactions = [0.5, 0.5]
        infecteds = [Matrix{Bool}(trues(2,4))]
        beta = 1.0
        @test sum(infect(susceptibles, infecteds, interactions, beta, 1)) == length(susceptibles)

        # Test 2: Transmission modulated by interactions
        susceptibles = [Matrix{Bool}([true  false false false;
                                      false false true false]) for _ in 1:10]  # Strain 1 susceptible
        infecteds = [Matrix{Bool}([false false true false;
                                   false false true false])]  # Strain 2 infected
        interactions = [1.0, 0.5]  # Strain 1's interactions
        beta = 1.0

        # Expected probability: 1.0 * 0.5 = 0.5
        # Expected infections: ~50% of 10 = 5
        result = infect(susceptibles, infecteds, interactions, beta, 1)
        @test 2 <= sum(result) <= 8  # Allow probabilistic variation
    end

    @testset "Disease Handlers" begin
        # Initialize test population
        pop = [Matrix{Bool}([true  false false false;
                             true false false false]) for _ in 1:100]
        pop[3] = Matrix{Bool}([false false true false;
                               false false true false])

        @testset "SI Model" begin
            # Test infection progression
            current_pop = deepcopy(pop)
            handle_infection(current_pop, Vector{Bool}(trues(100)), 1, 0.5, [1.0, 0.5])
            @test sum(m[1,3] for m in current_pop) > 0
        end

        @testset "SIR Model" begin
            # Test recovery mechanism
            current_pop = deepcopy(pop)
            current_pop[1][1,3] = true
            current_pop[1][1,1] = false
            handle_recovery(current_pop, Vector{Bool}(trues(100)), 1, 1.0)
            @test current_pop[1][1,4]  # Should be recovered
        end

        @testset "SEIR Model" begin
            # Test exposure mechanism
            current_pop = deepcopy(pop)
            current_pop[1][1,2] = true
            current_pop[1][1,1] = false
            handle_exposed_infection(current_pop, Vector{Bool}(trues(100)), 1, 1)
            @test current_pop[1][1,3] # should be infected
        end

        @testset "SEIRS Model" begin
            # Test immunity loss mechanism
            current_pop = deepcopy(pop)
            current_pop[1][1,4] = true
            current_pop[1][1,1] = false
            handle_immunity_loss(current_pop, Vector{Bool}(trues(100)), 1, 1.0)
            @test current_pop[1][1,1] # should be susceptible            
        end
    end

    @testset "Full Simulation" begin
        # Basic SI model test
        pop = [Matrix{Bool}([true  false false false;
                             true false false false]) for _ in 1:100]
        pop[3] = Matrix{Bool}([false false true false;
                               false false true false])

        results = coinfection_simulator(
            pop,
            ones(Int, 100),
            [1.0 0.9; 0.9 1.0],  # Interaction matrix (1x1)
            ["si", "si"],
            0.0,       # No mortality
            [0.0, 0.0],
            0.0,       # No births
            [1.0, 1.0],     # Full transmission
            2,         # 2 time steps
            1
        )

        # Test infection spread
        final_pop = results[1][end]
        @test sum(m[1,3] for m in final_pop) > 1
    end

    @testset "Edge Cases" begin

        pop = [Matrix{Bool}([true  false false false;]) for _ in 1:100]
        # Test empty population
        @test_throws AssertionError coinfection_simulator(
            Vector{Matrix{Bool}}(),
            Vector{Int64}(),
            zeros(0,0),
            String[],
            0.0,
            Float64[],
            0.0,
            Float64[],
            1,
            1
        )

        # Test full immunity
        immune_pop = deepcopy(pop)
        for m in immune_pop
            m[1,4] = true  # All recovered
        end
        results = coinfection_simulator(
            immune_pop,
            zeros(Int, 100),
            [1.0;;],
            ["sir"],
            0.0,
            [0.0],
            0.0,
            [1.0],
            2,
            1,
            recovery=[1.0]
        )
        @test all(m[1,4] for m in results[1][end])
    end

    @testset "Validation Checks" begin
        pop = [Matrix{Bool}([true  false false false;]) for _ in 1:100]

        # Test invalid disease type
        @test_throws AssertionError coinfection_simulator(
            pop,
            ones(Int, 100),
            [1.0;;],
            ["invalid"],
            0.0,
            [0.0],
            0.0,
            [1.0],
            1,
            1
        )

        # Test mortality bounds
        @test_throws AssertionError coinfection_simulator(
            pop,
            ones(Int, 100),
            [1.0;;],
            ["si"],
            1.5,  # Invalid mortality
            [0.0],
            0.0,
            [1.0],
            1,
            1
        )
    end
end