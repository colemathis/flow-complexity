#==============================================================================#
# IMPORTS
#==============================================================================#

import DataFrames
import Random
import Graphs
import Dates
using SparseArrays

#==============================================================================#
# DATA TYPES
#==============================================================================#

# Set the backward rate to a constant
const backward_rate = 1.0

#==============================================================================#

mutable struct Simulation
    """
    A mutable struct representing a simulation.
    """

    # define dictionaries for parameters, output
    params   ::Dict{Symbol, Any}
    output   ::Dict{Symbol, Any}
    ensemble ::Ensemble

    run             ::Function
    from_sim_array  ::Function
    save_A_matrix   ::Function
    load_A_matrix   ::Function

end

#==============================================================================#

function Simulation(; user_params = Dict{Symbol, Any}())
    """
    Create a Simulation object with default and user-defined parameters.
    """

    # apply default parameter values
    params = apply_params_default_values()

    # apply user-defined parameter values
    params = apply_params_user_values(params, user_params)
    
    # create output dictionary
    output = Dict{Symbol, Any}()

    # generate ensemble of chemostats
    ensemble = generate_ensemble_from_specs(params)

    # variables defined — let’s define the methods
    # in this case though we need to create the object beforehand as we’ll be using 
    # it as arguments to the methods
    sim = Simulation(params, output, ensemble,
                    () -> nothing,
                    () -> nothing,
                    () -> nothing,
                    () -> nothing
                )

    # assign wrapper functions
    sim.run =               ()           -> run_simulation(sim)
    sim.from_sim_array =    (sim_number) -> launch_simulation_from_sim_array(sim_number)
    sim.save_A_matrix =     (fn)         -> save_A_matrix(sim, fn)
    sim.load_A_matrix =     (fn)         -> load_A_matrix!(sim, fn)

    # return the objects we previously created
    return sim
    
end

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function launch_simulation_from_sim_array(sim; dry_run=false)
    """
    Launch a simulation from the simulation array.
    """

    # Load the parameter array from the CSV file
    array_fn = "./data/params.csv"
    
    # Check if sim is a string and convert to Int if necessary
    if typeof(sim) == String
        sim = parse(Int, sim)
    end
    
    println("Loading parameters for simulation $sim from file $array_fn")

    # Read the CSV file into a DataFrame
    df = DataFrames.DataFrame(CSV.File(array_fn))

    # Convert each row to a dictionary, adjust data types if needed and add to array
    params_dict_array = convert_param_df_to_dict(df)

    # Get the parameter dictionary for the specified simulation
    params_dict = params_dict_array[sim]
    println("Parameters loaded.\n")

    # Create a Simulation object with the loaded parameters
    simulation = FlowComplexity.Simulation(; user_params = params_dict)

    # Print the parameters (skipping matrices)
    print_parameters(simulation)

    # Run the simulation
    run_simulation(simulation, dry_run=dry_run)

    return simulation
    
end

#==============================================================================#

function apply_params_default_values()
    """
    Create a dictionary of default parameter values.
    """

    default_values = Dict(
        # hidden parameters - either unused or should not be customizable
        # fix eventually
        :N_inflow           => 1,
        :sim_number         => 1,

        # simulation parameters
        :save_interval    	=> 1,
        :method         	=> "tau-leaping",
        :dt             	=> 1e-3,
        :random_seed    	=> "random",
    
        # physical parameters
        :total_time     	=> 100,
        :initial_mass       => 0,
    
        # topological parameters
        :graph_type         => "lattice-2way",
        :randomize_edges    => false,
        :N_reactors     	=> 25,
    
        # reaction rates
        :inflow_mols        => 0,
        :forward_rate   	=> 1e-3,
        :diffusion_rate   	=> 1e-3,
        :outflow_rate       => 1e-3
    )
    
    return default_values

end

#==============================================================================#

function apply_params_user_values(previous_params, user_params)
    """
    Apply user-defined parameter values.
    """

    # create dict for final parameters
    final_params = copy(previous_params)
    
    # overwrite as needed with the ones provided by the user
    for (user_key, user_value) in user_params
        final_params[user_key] = user_value
    end
    
    return final_params

end

#==============================================================================#

function print_parameters(simulation)
    """
    Print the parameters of a simulation.
    """

    for (k, v) in simulation.params
        if !(v isa SparseMatrixCSC)
            println("   $k = $v")
        else
            println("   $k = <sparse matrix>")
        end
    end
    println("")

end

#==============================================================================#

function convert_param_df_to_dict(df)
    """
    Convert a DataFrame of parameters to an array of dictionaries.
    """

    params_dict_array = []
    for row in eachrow(df)
        d = Dict()
        for (name, value) in pairs(row)
            if (typeof(value) == CSV.String7) | (typeof(value) == CSV.String15)
                value = convert(String, value)
            end
            if name == :A
                value = maybe_restore_A(value)
            end
            merge!(d, Dict(name => value))
        end
        push!(params_dict_array, d)
    end

    return params_dict_array

end

#==============================================================================#

function maybe_restore_R(value)
    """
    Restore the R matrix from a string representation.
    """

    if (typeof(value) == CSV.String7) | (typeof(value) == CSV.String15)
        value = String(value)
    end
    if value isa AbstractString && startswith(value, "sparse(")
        m = match(r"^sparse\(\[(.*?)\],\s*\[(.*?)\],\s*\[(.*?)\](?:,\s*(\d+)\s*,\s*(\d+))?\)$", strip(value))
        if m !== nothing
            I = parse.(Int,     split(m.captures[1], r"\s*,\s*"))
            J = parse.(Int,     split(m.captures[2], r"\s*,\s*"))
            V = parse.(Float64, split(m.captures[3], r"\s*,\s*"))
            value = m.captures[4] === nothing ?
                sparse(I, J, V) :
                sparse(I, J, V, parse(Int, m.captures[4]), parse(Int, m.captures[5]))
        end
    end
    return value

end

#==============================================================================#

function generate_ensemble_from_specs(params)
    """
    Generate an Ensemble of chemostats from the given parameters.
    """

    # generate a list of chemostat specifications which can be passed to the Ensemble constructor
    chemostat_list = gen_chemostat_spec_list(params[:N_reactors],
                                             params[:forward_rate],
                                             backward_rate,
                                             params[:diffusion_rate],
                                             params[:outflow_rate])

    # create an ensemble of chemostats with the list of specifications
    ensemble = Ensemble(params[:N_reactors],
                             params[:graph_type],
                             params[:N_inflow],
                             params[:randomize_edges],
                             params[:initial_mass],
                             chemostat_list=chemostat_list)

    return ensemble

end

#==============================================================================#

function run_simulation(simulation; dry_run=false)
    """
    Run a simulation with the given parameters.
    """

    # set random seed to either a random value or user-defined value
    set_random_seed(simulation)

    # initialize output variables
    simulation.output[:timestamps] = []
    simulation.output[:populations] = []

    # start timing, then run the model
    elapsed_time = @elapsed begin

        # run the simulation using the SSA algorithm
        if simulation.params[:method] == "exact"
            println("Launching simulation using exact algorithm...")
            evolve_distributed_exact(simulation)

        # run the simulation using the tau-leaping algorithm
        elseif simulation.params[:method] == "tau-leaping"
            println("Launching simulation using tau-leaping algorithm...")
            evolve_distributed_tau_leaping(simulation; dry_run=dry_run)

        # error if method is unknown
        else
            println("error: method unknown")
            exit()
        end

    # end timing block
    end

    # calculate time series and delete molecules
    ts = calculate_time_series(simulation)
    simulation.output[:timeseries] = ts
    delete!(simulation.output, :populations)

    # print the elapsed time
    elapsed_time = round(elapsed_time, digits = 2)
    println("Sim Completed. Time taken: $elapsed_time seconds.")
    simulation.output[:total_time] = elapsed_time

    # save the simulation data
    save_data(simulation)
    println("")

    return simulation
    
end

#==============================================================================#

function set_random_seed(simulation)
    """
    Set the random seed for the simulation to either a random or user-defined value.
    """

    if simulation.params[:random_seed] != "random"
        Random.seed!(simulation.params[:random_seed])
        println("Random seed set to custom value: ", simulation.params[:random_seed])
    else
        random_seed = rand(UInt)   # generate a random seed
        Random.seed!(random_seed)  # set the RNG to this new seed
        println("Generated random seed: ", random_seed)
    end

end

#==============================================================================#

function save_checkpoint(sim, tau::Float64)
    """
    Save the current time and populations to the output.
    """
    
    # cap the precision on the time value
    i = round(tau, digits=3)
    println(Dates.now(), " --  saving at t=$i")
    
    # record the current time
    push!(sim.output[:timestamps], i)
    
    # create an empty array for the populations at this time
    push!(sim.output[:populations], [])

    # loop over each reactor
    for i in 1:sim.params[:N_reactors]

        # create a copy of the molecules for reactor i
        copy_of_molecules = copy(sim.ensemble.reactors[i].molecules)

        # push this copy to the populations array
        push!(sim.output[:populations][end], copy_of_molecules)
        # the resulting structure is thus sim.output[:populations][time_index][reactor_index]

    end
    
end

#==============================================================================#

function calculate_time_series(sim)
    """
    Calculate the time series DataFrame from the simulation output.
    """

    # create empty DataFrame for time series
    ts = DataFrames.DataFrame()

    # extract data from simulation parameters & output
    sim_number = sim.params[:sim_number]
    nt = length(sim.output[:timestamps])
    nchem = sim.params[:N_reactors]

    # loop over each time point & chemostat to build time series
    for j in 1:nt
        time = sim.output[:timestamps][j]
        for k in 1:nchem

            # get molecules at this time & chemostat
            mols = sim.output[:populations][j][k]

            # get unique molecules
            unique_mols = unique(mols)

            # loop over unique molecules to get frequencies
            for m in unique_mols

                # count frequency
                f = count(x -> x == m, mols)

                # append row to time series with dimensions
                r = (sim_number = sim_number,
                     time = time,
                     chemostat_id = k,
                     integer = m,
                     frequency = f)
                push!(ts, r)

            end
        end
    end

    return ts
    
end

#==============================================================================#
# END OF FILE
#==============================================================================#

