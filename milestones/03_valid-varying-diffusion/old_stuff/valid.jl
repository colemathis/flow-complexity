# Import necessary packages
using DrWatson
using JLD2

# Activate the project environment
# @quickactivate :FlowComplexity

# Define the file path
file_path = joinpath("./data", "data.jld2")

# Load the data from the specified file
@load file_path sim_array

using StatsBase

function convert_timeseries_to_tidy_df(timeseries)

    # get the recorded variables of the timeseries dict
    recorded_vars = [k for k in keys(timeseries[1])]

    # get the times
    times = [t for t in keys(timeseries[1][recorded_vars[1]])]

    # get the reactors
    reactors = collect(keys(timeseries))

    # create a data array and fill, for every reactor and time frame, the variable and its value
    data = []
    for r in reactors
        for t in times
            for var in recorded_vars
                if var == :complete_timeseries #p3 not too sure what’s the difference here
                    if timeseries[r][:complete_timeseries][t] != []
                        time_counts = countmap(timeseries[r][:complete_timeseries][t])
                        for (v,c) in time_counts
                            push!(data, Dict("reactor"=> r,"time"=> t, "variable"=>string(v), "value"=> c))
                        end
                    end
                else
                    push!(data, Dict("reactor"=> r,"time" => t, "variable" => String(var), "value"=> get(timeseries[r][var],t,0) ))
                end
            end
        end
    end

    # wrap this up in a dataframe
    tidy_df = DataFrame(time = map(x -> x["time"], data),
                        reactor = map(x -> x["reactor"], data),
                        variable = map(x -> x["variable"], data),
                        value = map(x -> x["value"], data))

    # and return it
    return tidy_df
    
end

# pop_array = zeros(Int, nchem, max_t, nspecies)
nsim=100
nspecies=10

avg_pop = zeros(Float64, nsim, nspecies)

using DataFrames

n = 100
kds = exp10.(LinRange(-6,2,n))

for i in 1:nsim
    print(".")
    sim = sim_array[i]
    p = convert_timeseries_to_tidy_df(sim.time_evolution)
    filter!(row -> row.time == 100, p)
    for j in 2:nspecies
        pp = filter(row -> row.variable == string(j), p)
        # println(pp)
        # ppp = mean(pp)
        gpp = groupby(pp, [:variable])
        gpp = combine(gpp, :value => mean => :value)
        avg_pop[i, j] = gpp.value[1]
        # exit()
    end
end
println("")


using Plots

p = scatter(xscale=:log10)
for i in 2:10
    a = avg_pop[:,i]    
    scatter!(p, kds, a)
end

display(p)
savefig("test.pdf")
