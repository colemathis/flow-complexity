include("Chemostats.jl")
include("Ensemble.jl")
include("SaveData.jl")

using Graphs
using Dates

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
    ### Rate Constants
    all_constants::Array{Float64}(undef,N_reactors,3)
    stablized_integers::Dict{Int64, Tuple{Int64}} # This might be a type issue... Fucking Julia
    ### Timeseries
    recorded_variables::Vector{Symbol}
    time_evolution::Dict{Any, Any}
    ### Simulation Number
    sim_number::Int64
    save_directory::String
end

# Generator for Simulation Struct
function Simulation(
    mass::Int64,
    graph_type::String,
    N_reactors::Int64,
    forward_rate::Float64,
    outflow_rate::Float64;
    total_time = 100.0,
    ouput_time = 1.0,
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
    single_rate = [forward_rate, 1.0, outflow_rate]
    all_rates = zeros(3, N_reactors)
    for i in 1:N_reactors
        all_rates[i,:] = single_rate
    end
    
    # Check output count
    if output_count > 1
        output_time = total_time / float(output_count)
    end

    return Simulation()
end