"""

##########################################################################

TEST_TIMING



##########################################################################

"""

# include("../src/Chemostats.jl")
# include("../src/TimeEvolve.jl")
#FIXME: to remove

using Random 
using JLD2
using FileIO
using DelimitedFiles
using DataFrames

# using DrWatson
# using Distributed
@quickactivate

@everywhere include("../src/Simulation.jl")

"""

##########################################################################

##########################################################################

Input:


Output:


##########################################################################

"""

function time_well_mixed_parameters()
    # Complete exploration of input parameters 

    #p1: try to reduce the range for parameters here since the code crashes at some point
    # ϵ_choices = [0.0, 0.0001, 0.001, 0.01, 0.1]
    ϵ_choices = [0.0, 0.0001]
    # outflow_choices = [0.0, 0.001]
    outflow_choices = [0.0, 0.001]
    iteration_choices = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]
    mass_choices = [100, 200, 500, 1000, 2000, 5000, 10000, 20000]

    timing_df = DataFrame(ϵ = Float64[], outflow = Float64[], iterations = Int64[], mass = Int64[], time = Float64[])
    ran_first = false
    for mass in mass_choices
        println("mass=", mass)
        for max_iterations in iteration_choices
            println("max_iterations=",max_iterations)
            for outflow in outflow_choices
                println("outflow=",outflow)
                for ϵ in ϵ_choices
                    println("ϵ=",ϵ)
                    well_mixed_rates = [((1.0 - outflow)/2.0) + ϵ, ((1.0 - outflow)/2.0) - ϵ , outflow] # Constructive, destructive, outflow

                    molecules = repeat([1], mass)

                    # this was the original call:
                    # well_mixed_chemostat = Chemostat(0, [], [], molecules, well_mixed_rates, mass, mass)

                    # here’s the constructor and the call being made in Ensemble.jl
                    #=
                    function Chemostat(ID                       ::Int64,
                                       reaction_rate_constants  ::Vector{Float64}; 
                                       molecules                = Vector{Int64}[],
                                       mass_fixed               = false, 
                                       neighbors                = Vector{Int64}[], 
                                       neighbor_flows           = Vector{Float64}[],
                                       stabilized_integers      = Vector{Int64}[]
                                       )

                    this_chemostat = Chemostat(r, reaction_rate_constants,
                                        molecules = molecules,
                                        mass_fixed = mass_fixed,
                                        neighbors = these_neigbors,
                                        stabilized_integers = stabilized_integers)
                    =#

                    #p1: if I’m using the constructor I get the following error:
                    # "type Chemostat has no field mass_fixed"
                    
                    # whereas if I try to instantiate directly the struct, the evolve_well_mixed()
                    # function below fails at random points while iterating over parameter permutations
                    
                    well_mixed_chemostat = Chemostat(0,
                                                    well_mixed_rates,
                                                    molecules=molecules, 
                                                    mass_fixed=true, 
                                                    neighbors=[0],
                                                    neighbor_flows=[1.0], 
                                                    stabilized_integers=[0])

                    record = [:molecule_count, :average_length]
                    if !ran_first

                        #p1 this is the original call:
                        # @elapsed evolution_out = evolve_well_mixed(well_mixed_chemostat, max_iterations, record);

                        # function definition:
                        # function evolve_well_mixed(chemostat::Chemostat, tau_max::Float64, output_freq::Float64, outputs::Array{Symbol,1})

                        # the modified call I’ve tried:
                        @elapsed evolution_out = evolve_well_mixed(well_mixed_chemostat, 100.0, 0.1, record);

                    end

                    #p1 and this fails too
                    # runtime = @elapsed evolution_out = evolve_well_mixed(well_mixed_chemostat, max_iterations, record)
                    runtime = @elapsed evolution_out = evolve_well_mixed(well_mixed_chemostat, 100.0, 0.1, record);
                    
                    results = [ϵ, outflow, max_iterations, mass, runtime]
                    push!(timing_df, results)
                end
            end
        end
    end
    writedlm("data/timing_results.csv", Iterators.flatten(([names(timing_df)], eachrow(timing_df))), ',')
end

time_well_mixed_parameters()