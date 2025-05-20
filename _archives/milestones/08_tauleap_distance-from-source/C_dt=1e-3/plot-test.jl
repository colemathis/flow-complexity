using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using DataFrames, Statistics, Plots, LaTeXStrings

using JLD2
using Graphs

# sim_array = FC.load_simulation_array()

# nchem = sim_array[1].params[:N_reactors]
nchem = 25

# pop = FC.create_populations_dataframe(sim_array)

pop = CSV.read("data/timeseries.csv", DataFrame)
# @save "data/pop.jld2" pop
# @load "data/pop.jld2" pop

sim_numbers = [30,60,100]
markers = [:circle, :diamond, :cross]
colors = [:blue, :red, :green]

s = scatter(title = "Integer average vs distance from source", legend=nothing)

distances = CSV.read("distance_from_source.csv", DataFrame)

for i in 1:3

    sim_number = sim_numbers[i]
    
    # g = sim_array[sim_number].ensemble.ensemble_graph
    
    # dijkstra_result = dijkstra_shortest_paths(g, 1)
    # dist = []
    # for i in 1:nchem
    #     d = dijkstra_result.dists[i]
    #     push!(dist, d)
    # end

    dist = filter(row -> row.sim_number == sim_number, distances)
    sort!(dist, :chemostat_id)
    dist = dist.distance_from_source
    
    f = filter(row -> row.sim_number == sim_number, pop)
    filter!(row -> row.time == 100, f)
    
    int_avg = []
    for i in 1:nchem
        ff = filter(row -> row.chemostat_id == i, f)
        ai = FC.calculate_assembly_index.(ff.integer)
        a = sum(ai .* ff.frequency) / sum(ff.frequency)
        push!(int_avg, a)
    end

    scatter!(s, dist, int_avg, color=colors[i])

end

mkpath("figs")
savefig("figs/test-julia.pdf")    

