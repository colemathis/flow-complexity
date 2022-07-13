include("Chemostats.jl")
include("Ensemble.jl")
include("SaveData.jl")
include("TimeEvolve.jl")

using Graphs
using Dates
using GraphIO


## Simulation Struct
mutable struct Simulation
    ## This type contains all the information needed to produce and run a simulation
    ### Dynamic information
    total_time::Float64
    output_time::Float64
    output_count::Int64
    random_seed::Int64
    ### Topological information
    graph_type::String
    N_reactors::Int64
    N_inflow::Int64
    ### Rate Constants and mass
    mass::Int64
    all_constants::Array{Float64}
    stablized_integers::Dict{Int64, Tuple{Int64}} # This might be a type issue... Fucking Julia
    ### Timeseries
    recorded_variables::Vector{Symbol}
    time_evolution::Dict{Any, Any}
    ### Simulation Number
    sim_number::Int64
    save_directory::String

    ### The Ensemble
    ensemble::Ensemble
end

# Generator for Simulation Struct
function Simulation(
    mass::Int64,
    graph_type::String,
    N_reactors::Int64,
    forward_rate::Float64,
    outflow_rate::Float64;
    backword_rate = 1.0,
    total_time = 100.0,
    output_time = 1.0,
    output_count = -1,
    N_inflow = 1,
    random_seed = parse(Int64, Dates.format(now(), "SSMMHHddmm")),
    recorded_variables = [:complete_timeseries],
    stabilized_integers = Dict(-1 => Tuple(-1,))
)
    ## Figure out what you need to do here! 

    # Get the simulation number of save name 
    sim_number = get_sim_number()
    save_name = datadir("sims", string(sim_number))

    # Form rates
    single_rate = [forward_rate*(1.0/mass), 1.0, outflow_rate]
    all_constants = zeros(N_reactors, 3)
    for i in 1:N_reactors
        all_constants[i,:] = single_rate
    end
    
    # Check output count
    if output_count > 1
        output_time = total_time / float(output_count)
    end
    time_series = Dict{Any, Any}()
    this_ensemble = generate_ensemble_from_specs(graph_type,
                                                 N_reactors,
                                                 N_inflow,
                                                 forward_rate,
                                                 backword_rate,
                                                 outflow_rate,
                                                 mass)

    return Simulation(total_time, 
                      output_time,
                      output_count,
                      random_seed,
                      graph_type,
                      N_reactors,
                      N_inflow,
                      mass,
                      all_constants,
                      stabilized_integers,
                      recorded_variables,
                      time_series,
                      sim_number,
                      save_name,
                      this_ensemble)
end

function generate_ensemble_from_specs(graph_type, N_reactors, N_inflow, forward_rate, backward_rate, outflow_rate, mass)
    ## Generate the chemostat specifications
    chemostat_list = gen_chemostat_spec_list(N_reactors, forward_rate, backward_rate, outflow_rate)
    ## Generate the ensemble_graph
    this_ensemble = Ensemble(N_reactors, graph_type, N_inflow, mass, chemostat_list=chemostat_list)

    return this_ensemble
end

function get_parameters(sim)
    parameter_fields = [:total_time, :output_time, :output_count,
                        :random_seed, :graph_type, :N_reactors,
                        :N_inflow,:mass, :recorded_variables,
                        :sim_number, :save_directory]
    parameters = Dict(key=>getfield(sim, key) for key ∈ parameter_fields )

    all_rate_constants = sim.all_constants
    forward_rates = all_rate_constants[:,1]
    outflow_rates = all_rate_constants[:,3]
    parameters[:ave_forward_rate] = mean(forward_rates)
    parameters[:ave_outflow_rate] = mean(outflow_rates)
 
    return parameters
end

function RunSimulation(sim)

    time_series_data = evolve_distributed(sim.ensemble,
                                          sim.total_time,
                                          sim.output_time,
                                          sim.recorded_variables,
                                          sim.random_seed)
    sim.time_evolution = time_series_data
    println("Sim Completed")
    parameters = get_parameters(sim)

    save_data(time_series_data, parameters, sim.ensemble, sim.sim_number)
   
end