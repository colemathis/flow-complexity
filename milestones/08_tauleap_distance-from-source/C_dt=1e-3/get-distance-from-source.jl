using Graphs
using JLD2
using DataFrames
using CSV

nchem = 25
sim_numbers = [30, 60, 100]

sim1 = JLD2.load("data/sims/000030.jld2")
sim2 = JLD2.load("data/sims/000060.jld2")
sim3 = JLD2.load("data/sims/000100.jld2")

sim_array = [sim1["sim"], sim2["sim"], sim3["sim"]]

df = DataFrame(sim_number=Int[], chemostat_id=Int[], distance_from_source=Int[])

for i in 1:3

    sim_number = sim_numbers[i]
    
    g = sim_array[i].ensemble.ensemble_graph
    
    dijkstra_result = dijkstra_shortest_paths(g, 1)
    dist = []
    for j in 1:nchem
        d = dijkstra_result.dists[j]
        push!(df, (sim_number=sim_numbers[i], chemostat_id=j, distance_from_source=d))
    end

end

CSV.write("distance_from_source.csv", df)