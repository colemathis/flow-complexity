include("Pipeline.jl")
include("Chemostats.jl")
include("Ensemble.jl")
include("EvolveStochastic.jl")
include("EvolveTauleaping.jl")
include("Analysis.jl")

using CSV, DataFrames
using Pkg
using Random

# Set the backward rate to a constant
const backward_rate = 1.0

#==============================================================================#
# DATA TYPES
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
                                             params[:diffusion_rate])

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

    # sim_number_string = lpad(sim_number,6,"0")
    # sim.params[:save_name] = projectdir("milestones", sim.params[:save_name], "data", "sims", sim_number_string)
    # sim.params[:save_name] = projectdir("milestones", sim.params[:save_name], "data", "sims")
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

    save_data(sim)
    println("")
    
end

#==============================================================================#

function sim_number_string(sim_number)

    sim_number_string = lpad(sim_number, 6, "0")

    return sim_number_string

end

#==============================================================================#

function save_data(sim)

    # sim_number = string(sim.params[:sim_number])
    sim_number = sim.params[:sim_number]
    sim_number_str = sim_number_string(sim_number)
    # fn = joinpath(sim.params[:save_name], "simulation.jld2")
    # fn = joinpath(sim.params[:save_name], "$sim_number_str.jld2")
    # save(fn, Dict("sim" => sim))

    # save sim.output[:timeseries] as a CSV file
    # fn = joinpath(sim.params[:save_name], "timeseries.csv")
    sim_dir = joinpath(sim.params[:save_name], "$sim_number_str")
    mkpath(sim_dir)
    fn = joinpath(sim_dir, "timeseries.csv")
    CSV.write(fn, sim.output[:timeseries])

    io_nodes = DataFrame(sim_number=Int[], chemostat_in=Int[], chemostat_out=Int[])
    g = sim.ensemble.ensemble_graph
    for e in edges(g)
        push!(io_nodes, (
            sim_number = sim_number,
            chemostat_in = src(e),
            chemostat_out = dst(e)
        ))
    end
    sim_dir = joinpath(sim.params[:save_name], "$sim_number_str")
    mkpath(sim_dir)
    fn = joinpath(sim_dir, "graph.csv")
    CSV.write(fn, io_nodes)

    rel = relpath(fn, pwd())
    println("Data saved in $rel")

end

#==============================================================================#

function save_checkpoint(sim, tau::Float64)
    
    i = round(tau, digits=3)
    println("   saving at t=$i")
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

    ts = DataFrame()

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
