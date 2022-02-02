# July 29 2021 
# Exploring parameters of well mixed model 

# Key parameters are Total Mass, Iterations, outflow rate, and epsilion (constructive process - destructive process = 2 epsilion)
include("../src/Chemostats.jl")
include("../src/TimeEvolve.jl")
include("process_bson_to_csv.jl")

using DrWatson
using Random 
using JLD2
using FileIO
using Dates


function test_chemostat()
    tau_max = 100.0

    mass = 1000
    
    outflow = 1.0
    reaction_rates = [10.0*(1.0/(mass)), 1.0, outflow] # Constructive, destructive, outflow

    molecules = repeat([1], mass)

    well_mixed_chemostat = Chemostat(0, reaction_rates, molecules)
    record = [:molecule_count, :average_length, :complete_timeseries]
    evolution_out = evolve_well_mixed(well_mixed_chemostat, tau_max, 1.0, record)
    save("data/raw/test_run.bson", evolution_out)
    tidy_df = bson_to_tidy_df("data/raw/test_run.bson")
    writedlm("data/raw/test_run.csv",Iterators.flatten(([names(tidy_df)], eachrow(tidy_df))), ',')

end

function output_test_run()

end

function complete_well_mixed_parameters(mass, outflow_rates, forward_rates, repetitions)
    # Complete exploration of input parameters 
    for i in 1:repetitions
        for outflow in outflow_rates
            for f in forward_rates
                well_mixed_rates = [f*(1.0/(mass)), 1.0, outflow] # Constructive, destructive, outflow

                molecules = repeat([1], mass)
                well_mixed_chemostat = Chemostat(0, [], [], molecules, well_mixed_rates, mass, mass)

                molecules = repeat([1], mass)

                well_mixed_chemostat = Chemostat(0, [], [], molecules, well_mixed_rates, mass, mass)
                record = [:molecule_count, :average_length, :complete_timeseries]
                evolution_out = evolve_well_mixed(well_mixed_chemostat, 100., 1.0, record)

                params = [i, mass, outflow, f]
                save_data(evolution_out, "data/raw", params)
            end
        end
    end

end

function logunif(min, max)
    scale = log10(max) - log10(min)
    r = 10^(scale*rand() + log10(min))
    return r 
end

function save_data(data, parameters)
    fname = ""
    for p in parameters
        fname = fname*string(p)*"_"
    end
    fname = datadir("sims", savename(parameters, "bson", connector = "^"))
    save(fname, data)
end

function complete_line_reactors_parameters(mass, outflow_rates, forward_rates, repetitions)
    # Complete exploration of input parameters 
    for i in 1:repetitions
        for outflow in outflow_rates
            for f in forward_rates
                line_reactor_rates = [f*(1.0/(mass)), 1.0, outflow] # Constructive, destructive, outflow
                line_reactors = make_line_reactors(3, line_reactor_rates, mass, mass)

                record = [:molecule_count, :average_length, :var_length, :complete_timeseries]
                evolution_out = evolve_distributed(line_reactors, 100., 1.0, record)

                params = [i, mass, outflow, f]
                save_data(evolution_out, "data/raw_line_reactor", params)
            end
        end
    end

end


function complete_line_reactors_n_reactors(mass, outflow_rate, forward_rate, reactors)
    # Complete exploration of input parameters 
    for n in reactors
        seed = parse(Int64, Dates.format(now(), "SSMMHHddmm"))
        line_reactor_rates = [forward_rate*(1.0/(mass)), 1.0, outflow_rate] # Constructive, destructive, outflow
        line_reactors = make_line_reactors(n, line_reactor_rates, mass, mass)

        record = [:molecule_count, :average_length, :var_length, :complete_timeseries]
        evolution_out = evolve_distributed(line_reactors, 100., 1.0, record, seed)

        mode = "line-graph"
        d = @strdict n mass outflow_rate forward_rate mode seed
        save_data(evolution_out, d)
    end

end


