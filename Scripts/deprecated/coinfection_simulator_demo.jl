include("coinfection_simulator.jl")
using BenchmarkTools, Tidier, CategoricalArrays, DataFrames, Gadfly
pop = [hcat(Matrix{Bool}(trues(4, 1)), fill(false, 4, 3)) for _ in 1:100]
pop[1][1, 1] = false
pop[1][1, 3] = true
pop[2][2, 1] = false
pop[2][2, 3] = true
pop[3][3, 1] = false
pop[3][3, 3] = true
pop[4][4, 1] = false
pop[4][4, 3] = true
interaction_matrix = [i == j ? 1.0 : rand(0.7:0.01:1.3) for i in 1:4, j in 1:4]
ages = collect(rand(1:199) for _ in 1:100)
disease_mortality = collect(rand(0:0.02) for _ in 1:4)
transmission = collect(rand(0.001:0.1) for _ in 1:4)

@benchmark coinfection_simulator(
  initial_pop = $pop,
  ages = $ages,
  interactions = $interaction_matrix,
  disease_type = ["si", "sir", "seir", "seirs"],
  base_mortality = 0.001,
  disease_mortality = $disease_mortality,
  fecundity = 0.01,
  transmission = $transmission,
  time_steps = 100,
  age_maturity = 30,
  introduction = "simultaneous",
  latency = [1, 1, 2, 3],
  immunity_loss = [0, 0, 0, 0.02],
  recovery = [0, 0.1, 0.1, 0.1]
)

sim1 = coinfection_simulator(
  initial_pop = pop,
  ages = ages,
  interactions = interaction_matrix,
  disease_type = ["si", "sir", "seir", "seirs"],
  base_mortality = 0.001,
  disease_mortality = disease_mortality,
  fecundity = 0.01,
  transmission = transmission,
  time_steps = 100,
  age_maturity = 30,
  introduction = "simultaneous",
  latency = [1, 1, 2, 3],
  immunity_loss = [0, 0, 0, 0.02],
  recovery = [0, 0.1, 0.1, 0.1]
)

pop_size = [length(vec) for vec in sim1[2]]
infected1 = [sum(m[1,3] for m in t) for t in sim1[1]]
infected2 = [sum(m[2,3] for m in t) for t in sim1[1]]
infected3 = [sum(m[3,3] for m in t) for t in sim1[1]]
infected4 = [sum(m[4,3] for m in t) for t in sim1[1]]

result_df = @chain DataFrame(Time=1:100, total = pop_size, infected1 = infected1, infected2 = infected2, infected3 = infected3, infected4 = infected4) begin
  @pivot_longer(total:infected4, names_to = "State", values_to = "Count")
end
result_df.State = categorical(result_df.State; levels=["infected1", "infected2", "infected3", "infected4", "total"], ordered=true)

p = Gadfly.plot(
  result_df,
  x=:Time,
  y=:Count,
  color=:State,
  Geom.line,
  Guide.colorkey(labels=["Total", "SI", "SIR", "SEIR", "SEIRS"]),
  Scale.y_continuous(minvalue=0)
)

# Save the plot as an image file in the same directory as the .jl file
output_path = "/Users/caryinstitute/Library/CloudStorage/GoogleDrive-pilowskyj@caryinstitute.org/My Drive/Grants/F32/simulation_example.png"
draw(PNG(output_path, 16cm, 12cm), p)