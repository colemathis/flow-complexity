using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using DataFrames, Statistics, Plots, LaTeXStrings

using JLD2
using Graphs

sim_array = FC.load_simulation_array()
nchem = sim_array[1].params[:N_reactors]

pop = FC.create_populations_dataframe(sim_array)
@save "data/pop.jld2" pop
# @load "data/pop.jld2" pop

nsim = length(sim_array)

k_ds = []
for i in 1:nsim
    k_d = sim_array[i].params[:outflow_rate]
    push!(k_ds, k_d)
end

nspecies = 10
avg_populations = zeros(Float64, nsim, nspecies)
f = filter(row -> row.time == 100, pop)
for i in 1:nsim
    ff = filter(row -> row.sim_number == i, f)
    gff = groupby(ff, [:integer])
    gff = combine(gff, :frequency => mean => :frequency)
    for j in 2:nspecies
        fff = filter(row -> row.integer == j, gff)
        avg_populations[i, j] = fff.frequency[1]
    end
end

s = scatter(xscale = :log10)

for i in 2:nspecies

    scatter!(s, k_ds, avg_populations[:,i], label = "$i")

end

mkpath("figs")
# display(s)
savefig("figs/populations-vs-kd.pdf")
