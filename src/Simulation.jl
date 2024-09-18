"""

##########################################################################

SIMULATION



##########################################################################

"""

include("Chemostats.jl")
include("Ensemble.jl")
include("SaveData.jl")
include("TimeEvolve.jl")

using Graphs
using Dates
using GraphIO

export Simulation, RunSimulation

"""

##########################################################################
SIMULATION STRUCT

Contains all that is needed to run a simulation.
##########################################################################

dynamic information:
    total_time          (float)         = total physical time of the simulation (e.g., 100.0)
    output_time         (float)         = interval for saving data (e.g., 1.0)
    output_count        (int)           = number of "snapshots" of data to save (e.g., 50)
    random_seed         (int)           = seed for the random number generator

topological information:
    graph_type          (string)        = topology of the reactors ("ER", "BA", "regular", "line", "lattice")
    N_reactors          (int)           = number of reactors (e.g., 1...N)
    N_inflow            (int)           = number of reactors that have an in-flow (e.g., 1...N)

rate constants and mass:
    mass                (int)           = total physical mass of the reactors (e.g., 100, 200, ...)
    all_constants       (float array)   = forward (constructive), backward (destructive) 
                                          and outflow rates (e.g., 1.0, ...)
    stablized_integers  (int,int dict)  = integers stabilized in a reactor #p4: to be implemented

timeseries:
    recorded_variables  (#p1 not sure)  = recorded variables (e.g., "complete_timeseries")
    time_evolution      (dict)          = time series data (e.g., time_series_data)

simulation metadata:
    sim_number          (int)           = simulation number (e.g. 1...N)
    save_time_series    (bool)          = 
    save_parameters     (bool)          = 
    save_graph          (bool)          = 
    save_simulation     (bool)          = 
    save_directory      (string)        = save directory (e.g., "data")
    sim_notes           (string)        = notes 
    
ensemble:
    ensemble            (Ensemble)      = ensemble (e.g., the_ensemble)

##########################################################################

"""

mutable struct Simulation
    
    ### Dynamic information
    total_time          ::Float64
    output_time         ::Float64
    output_count        ::Int64
    random_seed         ::Int64

    ### Topological information
    graph_type          ::String
    N_reactors          ::Int64
    N_inflow            ::Int64

    ### Rate Constants and mass
    mass                ::Int64
    all_constants       ::Array{Float64}
    stablized_integers  ::Dict{Int64, Tuple{Int64}} #p3: (cole) This might be a type issue... Fucking Julia 

    ### Timeseries
    recorded_variables  ::Vector{Symbol}
    time_evolution      ::Dict{Any, Any}

    ### Simulation Number
    sim_number          ::Int64
    save_time_series    ::Bool
    save_parameters     ::Bool
    save_graph          ::Bool
    save_simulation     ::Bool
    save_name           ::String
    sim_notes           ::String

    ### The Ensemble
    ensemble            ::Ensemble

end

"""

##########################################################################
GENERATOR FOR SIMULATION STRUCT

Parametrizes the Simulation struct.
##########################################################################

Input:

    mass                (int)           = (see struct)
    graph_type          (string)        = (see struct)
    N_reactors          (int)           = (see struct)
    forward_rate        (float)         = rate constant for A + B -> C (e.g., 1.0)
    outflow_rate        (float)         = rate constant for G -> ∅ (e.g., 1.0)
    backword_rate       (float)         = rate constant for D -> E + F (e.g., 1.0)
    total_time          (float)         = (see struct)
    output_time         (float)         = (see struct)
    output_count        (int)           = (see struct)
    N_inflow            (int)           = (see struct)
    notes               (string)        = (see struct)
    sim_number          (int)           = (see struct)
    random_seed         (int)           = (see struct)
    recorded_variables  (#p1 not sure)  = (see struct)
    stabilized_integers (int,int dict)  = (see struct)

Output:

    simulation          (Simulation)    = simulation

##########################################################################

"""

function Simulation(
    mass                ::Int64,
    graph_type          ::String,
    N_reactors          ::Int64,
    forward_rate        ::Float64,
    outflow_rate        ::Float64;
    backword_rate       = 1.0,
    total_time          = 100.0,
    output_time         = 1.0,
    output_count        = -1,
    N_inflow            = 1,
    notes               = "None",
    sim_number          = -1,
    save_time_series    = true,
    save_parameters     = true,
    save_graph          = true,
    save_simulation     = true,
    save_directory      = "",
    random_seed         = Random.rand(1:10000000),
    recorded_variables  = [:complete_timeseries],
    stabilized_integers = Dict(-1 => Tuple(-1,))
    )

    # Get the next simulation number available
    if sim_number < 0
        sim_number = get_sim_number()
    end
    sim_number_string = lpad(sim_number,6,"0")
    save_name = datadir("sims", save_directory, sim_number_string)

    # Create vectors containing the forward, outflow and backward rates
    # and assign these to every reactor
    single_rate = [forward_rate*(1.0/mass), 1.0, outflow_rate]
    all_constants = zeros(N_reactors, 3)
    for i in 1:N_reactors
        all_constants[i,:] = single_rate
    end
    
    # If output count specified, define output time
    if output_count > 1
        output_time = total_time / float(output_count)
    end

    # Placeholder for time series
    time_series = Dict{Any, Any}()

    # Generate an ensemble of N reactors with these parameters
    this_ensemble = generate_ensemble_from_specs(graph_type,
                                                 N_reactors,
                                                 N_inflow,
                                                 forward_rate,
                                                 backword_rate,
                                                 outflow_rate,
                                                 mass)

    # return the Simulation struct
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
                      save_time_series,
                      save_parameters,
                      save_graph,
                      save_simulation,
                      save_name,
                      notes,
                      this_ensemble)

end

"""

##########################################################################
GENERATE AN ENSEMBLE FROM GIVEN SPECIFICATIONS

Create an ensemble of chemostats from the given parameters.
##########################################################################

Input:

    graph_type      (string)    = ensemble graph type
    N_reactors      (int)       = number of reactors in the ensemble
    N_inflow        (int)       = #p1 to determine
    forward_rate    (float)     = forward reaction rate
    backward_rate   (float)     = backward reaction rate
    outflow_rate    (float)     = outflow reaction rate #p1 clarify
    mass            (int)       = mass of the ensemble #p1 clarify

Output:

    this_ensemble   (Ensemble)  = ensemble

##########################################################################

"""

function generate_ensemble_from_specs(graph_type, 
                                      N_reactors, 
                                      N_inflow, 
                                      forward_rate, 
                                      backward_rate, 
                                      outflow_rate, 
                                      mass
                                      )

    # generate a list of chemostat specifications which can be passed to the Ensemble constructor
    chemostat_list = gen_chemostat_spec_list(N_reactors, forward_rate, backward_rate, outflow_rate)

    # create an ensemble of chemostats with the list of specifications
    this_ensemble = Ensemble(N_reactors, graph_type, N_inflow, mass, chemostat_list=chemostat_list)

    return this_ensemble

end

"""

##########################################################################
GET ALL SIMULATION PARAMETERS

For a given simulation, get the parameters in a dictionary.
##########################################################################

Input:

    sim         (Simulation)    = simulation

Output:

    parameters  (dict)          = parameters

##########################################################################

"""

function get_parameters(sim)

    parameter_fields = [:total_time, :output_time, :output_count,
                        :random_seed, :graph_type, :N_reactors,
                        :N_inflow,:mass, :recorded_variables,
                        :sim_number, :sim_notes, :save_name]

    # create the dictionary that will hold all parameters
    parameters = Dict(key=>getfield(sim, key) for key ∈ parameter_fields )

    # fetch the individual rate constants from the rate dictionary
    all_rate_constants = sim.all_constants
    forward_rates = all_rate_constants[:,1]
    outflow_rates = all_rate_constants[:,3]

    # and calculate the average rates
    parameters[:ave_forward_rate] = mean(forward_rates)
    parameters[:ave_outflow_rate] = mean(outflow_rates)
 
    return parameters

end

"""

##########################################################################
RUN SIMULATION

Run the simulation, then save the time series and parameters to file.
##########################################################################

Input:

    sim     (Simulation)    = simulation

Output:

    (none)

##########################################################################

"""

function RunSimulation(sim)

    # run the simulation
    # p1: what is "evolved" ?
    time_series_data = evolve_distributed(sim.ensemble,
                                          sim.total_time,
                                          sim.output_time,
                                          sim.recorded_variables,
                                          sim.random_seed)

    # sim completed, save the time series data in the struct
    sim.time_evolution = time_series_data

    println("Sim Completed")

    # save both sim time series data and parameters
    parameters = get_parameters(sim)
    save_data(sim, parameters)
   
end