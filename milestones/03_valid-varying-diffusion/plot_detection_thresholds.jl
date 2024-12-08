using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using DataFrames, Statistics, Plots, LaTeXStrings

using JLD2
# using Graphs

sim_array = FC.load_simulation_array()
nchem = sim_array[1].params[:N_reactors]

# pop = FC.create_populations_dataframe(sim_array)
# @save "data/pop.jld2" pop
@load "data/pop.jld2" pop

nsim = length(sim_array)

k_ds = []
for i in 1:nsim
    k_d = sim_array[i].params[:outflow_rate]
    push!(k_ds, k_d)
end

# for each simulation, calculate the total # of molecules
total_populations = zeros(Int, nsim)
f = filter(row -> row.time == 100, pop)
for i in 1:nsim
    ff = filter(row -> row.sim_number == i, f)
    gff = groupby(ff, [:integer])
    gff = combine(gff, :frequency => sum => :frequency)
    total_populations[i] = sum(gff.frequency)
end

# for each simulation, calculate the max detectable integers
thresholds = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
nthresh = length(thresholds)
highest_integer = zeros(Int, nsim, nthresh)
f = filter(row -> row.time == 100, pop)
for i in 1:nsim
    ff = filter(row -> row.sim_number == i, f)
    gff = groupby(ff, [:integer])
    gff = combine(gff, :frequency => sum => :frequency)
    transform!(gff, :frequency => ByRow(x -> x / total_populations[i]) => :fraction)
    for j in 1:nthresh
        fff = filter(row -> row.fraction > thresholds[j], gff)
        sorted_gff = sort(fff, :integer, rev=true)
        highest_integer[i, j] = sorted_gff.integer[1]
    end
end

# convert these integers to AIs
highest_AI = zeros(Int, nsim, nthresh)
for i in 1:nsim
    for j in 1:nthresh
        highest_AI[i, j] = FC.calculate_assembly_index(highest_integer[i, j])
    end
end

# plot integers
s = scatter(xscale = :log10)
for i in 1:nthresh
    thresh = thresholds[i]
    scatter!(s, k_ds, highest_integer[:, i], label="$thresh")
end

mkpath("figs")
# display(s)
savefig("figs/detection_thresholds_int.pdf")

# plot AIs
s = scatter(xscale = :log10)
for i in 1:nthresh
    thresh = thresholds[i]
    scatter!(s, k_ds, highest_AI[:, i], label="$thresh")
end

mkpath("figs")
# display(s)
savefig("figs/detection_thresholds_AI.pdf")

