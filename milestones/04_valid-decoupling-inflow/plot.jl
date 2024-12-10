using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using DataFrames, Statistics, Plots, LaTeXStrings

sim_array = FC.load_simulation_array()
pop = FC.create_populations_dataframe(sim_array)

# create a time series

# only keep the first chemostat
ts_allsims = copy(pop)
filter!(row -> row.chemostat_id == 2, ts_allsims)

# take average over simulations
gts_allsims = groupby(ts_allsims, [:time, :integer])
gts_allsims = combine(gts_allsims, :frequency => mean => :frequency)

p = plot(title="Time series for chemostat 1, all sims", xlabel="Time", ylabel="Count",
         legendtitle="Integers", legend=:outerright)

for i in 1:10
    f = filter(row -> row.integer == i, gts_allsims)
    plot!(p, f.time, f.frequency, label="$i")
end

mkpath("figs")
savefig("figs/time-series.pdf")

# create a histogram

filter!(row -> row.time == 100.0, gts_allsims)
filter!(row -> row.integer < 25, gts_allsims)
sort!(gts_allsims, :integer)

b = bar(xlabel="Integers", ylabel="Frequency", legend=false)
bar!(gts_allsims.integer, gts_allsims.frequency)
savefig("figs/hist.pdf")

# create a multiplot for all chemostats of sim 1

ts_multi = copy(pop)
sim_no = 1
filter!(row -> row.sim_number == sim_no, ts_multi)

# p = plot(layout = grid(5, 5), size = (900, 900), legend=nothing, link=:both, ylim = (0, 1000))
p = plot(layout = grid(5, 5), size = (900, 900), legend=nothing, link=:both, ylim = (0, 10000))

nspecies = 10
nchem = sim_array[sim_no].params[:N_reactors]
for i in 1:nchem
    f = filter(row -> row.chemostat_id == i, ts_multi)
    for j in 1:nspecies
        g = filter(row -> row.integer == j, f)
        x = g.time
        y = g.frequency
        plot!(p, x, y, subplot = i)
    end
end

savefig("figs/multi.pdf")

