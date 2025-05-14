#==============================================================================#
# IMPORTS
#==============================================================================#

import DataFrames
import Random
import Graphs
import Dates

#==============================================================================#
# DATA TYPES
#==============================================================================#

# Set the backward rate to a constant
const backward_rate = 1.0

#==============================================================================#

mutable struct Simulation

    params   ::Dict{Symbol, Any}
    output   ::Dict{Symbol, Any}
    ensemble ::Ensemble

end

#==============================================================================#

function Simulation(; user_params = Dict{Symbol, Any}())

    params = apply_params_default_values()
    params = apply_params_user_values(params, user_params)
    
    output = Dict{Symbol, Any}()
    ensemble = generate_ensemble_from_specs(params)
    
    return Simulation(params, output, ensemble)
    
end

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function launch_simulation(sim; dry_run=false)

    array_fn = "./data/params.csv"
    
    if typeof(sim) == String
        sim = parse(Int, sim)
    end
    
    println("Loading parameters for simulation $sim from file $array_fn")

    df = DataFrames.DataFrame(CSV.File(array_fn))

    params_dict_array = []
    for row in eachrow(df)
        d = Dict()
        for (name, value) in pairs(row)
            if (typeof(value) == CSV.String7) | (typeof(value) == CSV.String15)
                value = convert(String, value)
            end
            merge!(d,Dict(name => value))
        end
        push!(params_dict_array,d)
    end

    params_dict = params_dict_array[sim]
    println("Parameters loaded.\n")

    for (k, v) in params_dict
    	println("   $k = $v")
    end
    println("")
    
    simulation = FlowComplexity.Simulation(; user_params = params_dict)
    RunSimulation(simulation, dry_run=dry_run)

    return simulation
    
end

#==============================================================================#

function apply_params_default_values()

    default_values = Dict(
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
        :outflow_rate       => 1e-3,

        # hidden parameters - either unused or should not be customizable
        # fix eventually
        :N_inflow           => 1,
        :sim_number         => 1,
        :save_name          => get_relative_path()
    )
    
    
    return default_values

end

#==============================================================================#

function apply_params_user_values(previous_params, user_params)

    final_params = copy(previous_params)
    
    # overwrite with the ones provided by the user
    for (user_key, user_value) in user_params
        final_params[user_key] = user_value
    end
    
    return final_params

end

#==============================================================================#

function generate_ensemble_from_specs(params)

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

function RunSimulation(sim; dry_run=false)

    if sim.params[:random_seed] != "random"
        Random.seed!(sim.params[:random_seed])
        println("Random seed set to custom value: ", sim.params[:random_seed])
    else
        random_seed = rand(UInt)  # generate a random seed
        Random.seed!(random_seed)  # set the RNG to this new seed
        println("Generated random seed: ", random_seed)
    end

    sim_number = sim.params[:sim_number]

    sim.params[:save_name] = joinpath(dirname(Pkg.project().path), "milestones", sim.params[:save_name], "data", "sims")
    sim.output[:timestamps] = []
    sim.output[:populations] = []

    elapsed_time = @elapsed begin
        if sim.params[:method] == "exact"
            if dry_run == true
                println("Dry run for exact algorithm unsupported.")
                exit()
            end
            println("Launching simulation using exact algorithm...")
            evolve_distributed_exact(sim)
        elseif sim.params[:method] == "tau-leaping"
            println("Launching simulation using tau-leaping algorithm...")
            evolve_distributed_tau_leaping(sim; dry_run=dry_run)
        else
            println("error: method unknown")
            exit()
        end
    end

    # calculate time series and delete molecules
    ts = calculate_time_series(sim)
    sim.output[:timeseries] = ts
    delete!(sim.output, :populations)
    
    elapsed_time = round(elapsed_time, digits = 2)
    println("Sim Completed. Time taken: $elapsed_time seconds.")
    sim.output[:total_time] = elapsed_time

    save_data(sim)
    println("")
    
end

#==============================================================================#

function save_data(sim)

    sim_number = sim.params[:sim_number]
    sim_number_str = sim_number_string(sim_number)
    sim_dir = joinpath(sim.params[:save_name], "$sim_number_str")
    mkpath(sim_dir)
    fn = joinpath(sim_dir, "timeseries.csv")
    CSV.write(fn, sim.output[:timeseries])

    io_nodes = DataFrames.DataFrame(sim_number=Int[], chemostat_in=Int[], chemostat_out=Int[])
    g = sim.ensemble.ensemble_graph
    for e in Graphs.edges(g)
        push!(io_nodes, (
            sim_number = sim_number,
            chemostat_in = Graphs.src(e),
            chemostat_out = Graphs.dst(e)
        ))
    end
    sim_dir = joinpath(sim.params[:save_name], "$sim_number_str")
    mkpath(sim_dir)
    fn = joinpath(sim_dir, "graph.csv")
    CSV.write(fn, io_nodes)

    fn = joinpath(sim_dir, "meta.csv")
    meta = DataFrames.DataFrame(
        sim_number = sim.params[:sim_number],
        total_time = sim.output[:total_time],
        total_constructive_rxn = get(sim.output, :total_constructive_rxn, missing),
        skipped_constructive_rxn = get(sim.output, :skipped_constructive_rxn, missing),
        total_destructive_rxn = get(sim.output, :total_destructive_rxn, missing),
        skipped_destructive_rxn = get(sim.output, :skipped_destructive_rxn, missing),
        total_diffusion_rxn = get(sim.output, :total_diffusion_rxn, missing),
        skipped_diffusion_rxn = get(sim.output, :skipped_diffusion_rxn, missing),
        total_outflow_rxn = get(sim.output, :total_outflow_rxn, missing),
        skipped_outflow_rxn = get(sim.output, :skipped_outflow_rxn, missing),
    )
    CSV.write(fn, meta)

    rel = relpath(fn, pwd())
    println("Data saved in $rel")

end

#==============================================================================#

function save_checkpoint(sim, tau::Float64)
    
    i = round(tau, digits=3)
    println(Dates.now(), " --  saving at t=$i")
    # print(".")
    
    # record the current time
    push!(sim.output[:timestamps], i)
    
    # loop over reactors then and add vector of molecules
    push!(sim.output[:populations], [])
    for i in 1:sim.params[:N_reactors]
        copy_of_molecules = copy(sim.ensemble.reactors[i].molecules)
        push!(sim.output[:populations][end], copy_of_molecules)
    end

    # the resulting ndim array will have dimensions:
    # 1 = time id (NOT time stamp)
    # 2 = reactor id
    # 3 = molecule vector
    #
    # so if time id is "t", reactor id is "i" and we want to access molecule "m",
    # sim.populations[t, i, m]
    
end

#==============================================================================#

function calculate_time_series(sim)

    ts = DataFrames.DataFrame()

    sim_number = sim.params[:sim_number]
    nt = length(sim.output[:timestamps])
    nchem = sim.params[:N_reactors]
    for j in 1:nt
        time = sim.output[:timestamps][j]
        for k in 1:nchem
            mols = sim.output[:populations][j][k]
            unique_mols = unique(mols)
            for m in unique_mols
                f = count(x -> x == m, mols)
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

