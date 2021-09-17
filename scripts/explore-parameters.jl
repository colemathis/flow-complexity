# July 29 2021 
# Exploring parameters of well mixed model 

# Key parameters are Total Mass, Iterations, outflow rate, and epsilion (constructive process - destructive process = 2 epsilion)
include("../src/Chemostats.jl")
include("../src/TimeEvolve.jl")
using Random 
using JLD2
using FileIO


function test_chemostat()
    iterations = 100000

    mass = 1000
    
    outflow = 0.00
    ϵ = 0.005
    reaction_rates = [((1.0 - outflow)/2.0) + ϵ, ((1.0 - outflow)/2.0) - ϵ , outflow] # Constructive, destructive, outflow

    molecules = repeat([1], mass)

    well_mixed_chemostat = Chemostat(0, [], [], molecules, reaction_rates, mass, mass)
    record = [:molecule_count, :average_length]
    evolution_out = evolve_well_mixed(well_mixed_chemostat, iterations, record)

    params = [iterations, mass, outflow]
    save_data(evolution_out, "data/raw/", params)

end


function complete_well_mixed_parameters(mass, iteration_choices, outflow_rates, epsilions, repetitions)
    # Complete exploration of input parameters 
    for i in 1:repetitions
        for max_iterations in iteration_choices
            for outflow in outflow_rates
                for ϵ in epsilions
                    well_mixed_rates = [((1.0 - outflow)/2.0) + ϵ, ((1.0 - outflow)/2.0) - ϵ , outflow] # Constructive, destructive, outflow

                    molecules = repeat([1], mass)
                    well_mixed_chemostat = Chemostat(0, [], [], molecules, well_mixed_rates, mass, mass)

                    record = [:molecule_count, :average_length]
                    evolution_out = evolve_well_mixed(well_mixed_chemostat, max_iterations, record)

                    params = [i, max_iterations, outflow, ϵ]
                    save_data(evolution_out, "data/raw", params)
                end
            end
        end
    end

end

function logunif(min, max)
    scale = log10(max) - log10(min)
    r = 10^(scale*rand() + log10(min))
    return r 
end

function save_data(data, directory, parameters)
    fname = ""
    for p in parameters
        fname = fname*string(p)*"_"
    end
    fname = directory*"/"*fname*".bson"
    save(fname, data)
end
