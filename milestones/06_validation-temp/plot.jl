using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using Statistics, Plots, LaTeXStrings

sim_array = FC.load_simulation_array()

mkpath("figs")

nsim = length(sim_array)
timestamps = sim_array[1].output[:timestamps]
nt = length(timestamps)
nspecies = 100
nchem = 1

pop = zeros(Int, nsim, nt, nchem, nspecies)

for i in 1:nsim
    this_sim = sim_array[i]
    molecules = this_sim.output[:populations]
    for j in 1:nt, k in 1:nchem, l in 1:nspecies
        these_molecules = molecules[j][k]
        these_populations = count(x -> x == l, these_molecules)
        pop[i,j,k,l] = these_populations
    end
end

# average over sims
pop = mean(pop, dims=1)
pop = dropdims(pop, dims=1)

# average over chemostats
pop = mean(pop, dims=2)
pop = dropdims(pop, dims=2)

# Create plot with basic settings
p = plot(title="Time series", xlabel="Time", ylabel="Count",
         legendtitle="Integers", legend=:outerright)

for i in 1:10
    plot!(p, timestamps, pop[:, i], label="$i")
end
savefig("figs/time-series.pdf")

# keep the last (temporal) iteration
pop = pop[end,:]

ns = 25
b = bar(title="Integer Frequency", xlabel="Integers", ylabel="Frequency", legend=false)
bar!(1:ns, pop[1:ns])
savefig("figs/hist.pdf")
