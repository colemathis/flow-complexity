using JLD2
using CSV
using Pkg
using Arrow

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function create_params_file()

    if isfile("params.jl")
        println("params.jl already exists. Aborting to avoid overwrite.")
        exit(1)
    end

    params_text = """
    nrepeat = 1

    params_template = Dict(
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
    )
    """

    open("params.jl", "w") do io
        write(io, params_text)
    end

    println("Wrote params.jl in ", pwd())

end

#==============================================================================#

function create_slurm_file()

    if isfile("run.slurm")
        println("run.slurm already exists. Aborting to avoid overwrite.")
        exit(1)
    end

    slurm_text = """
    #!/bin/bash

    # Job identification
    #SBATCH --job-name=flow-sim              # descriptive job name
    #SBATCH --array=1-100                    # simulation indices

    # Logging (stdout & stderr go to the same file)
    #SBATCH --output=./data/logs/slurm_%A_%a.log   # %A = array master ID

    # Resources
    #SBATCH --partition=general              # 7‑day queue
    #SBATCH --qos=public                     # default QoS
    #SBATCH --time=24:00:00                  # wall time (hh:mm:ss)
    #SBATCH --mem=4G                         # memory

    #################################
    # SCRIPT BEGINS
    #################################

    module load julia
    mkdir -p ./data/logs
    srun flow launch \${SLURM_ARRAY_TASK_ID}
    """

    open("run.slurm", "w") do io
        write(io, slurm_text)
    end

    println("Wrote run.slurm in ", pwd())

end

#==============================================================================#

function get_relative_path()

    current_dir = pwd()
    milestones_dir = joinpath(dirname(Pkg.project().path), "milestones")
    relative_path = relpath(current_dir, milestones_dir)
    
    return relative_path

end

#==============================================================================#

function write_params_file(params_template, nrepeat)

    array_fn = "./data/params.csv"
    println("Writing parameters for simulation array in $array_fn")
    # Generate all combinations of parameters
    param_keys = collect(keys(params_template))
    values_list = [ isa(params_template[k], AbstractVector) ? params_template[k] : [params_template[k]] for k in param_keys ]
    param_tuples = Iterators.product(values_list...)
    params_array = [ Dict(zip(param_keys, t)) for t in param_tuples ]

    # Duplicate each parameter set nrepeat times so they can be edited independently
    if nrepeat > 1
        params_array = [deepcopy(d) for d in params_array for _ in 1:nrepeat]
    end

    for i in eachindex(params_array)
        params_array[i][:sim_number] = i
    end

    CSV.write(array_fn, params_array)
    println("Done.")
    println("")

end

#==============================================================================#

function extract_sims()

    # old function consolidating JLD2 files into a single CSV file

    directory = "./data/sims/timeseries"
    # nsim = count(isdir, readdir(directory, join=true))
    nsim = count(isfile, readdir(directory, join=true))
    sim_array = []
    params = []
    
    println("Found $nsim simulation files in data folder. Loading data...")
    for i in 1:nsim
        sim_number_str = sim_number_string(i)
        sim_path = joinpath(directory, "$sim_number_str.jld2")
        sim = load(sim_path)["sim"]
        push!(sim_array, sim)
    end
    timeseries = consolidate_array_time_series(sim_array)
    
    println("Loading done. Saving data...")
    
    CSV.write("data/timeseries.csv", timeseries)
    
    # skip saving as jld2 for now
    # @save "data/sim_array.jld2" sim_array

    io_nodes = DataFrame(sim_number=Int[], chemostat_in=Int[], chemostat_out=Int[])
    for i in 1:nsim
        sim_number = i
        g = sim_array[i].ensemble.ensemble_graph
        for e in edges(g)
            push!(io_nodes, (
                sim_number = sim_number,
                chemostat_in = src(e),
                chemostat_out = dst(e)
            ))
        end
    end
    CSV.write("data/graphs.csv", io_nodes)
    
    println("Saving done.")
    println("")
    
end

#==============================================================================#

function extract_sims2()

    # new function consolidating CSV files into a single CSV file
    # and Arrow files for faster loading
    
    sims_dir = "./data/sims"
    sim_folders = filter(f -> isdir(joinpath(sims_dir, f)), readdir(sims_dir))
    
    all_timeseries = DataFrame()
    all_graphs = DataFrame()
    
    for folder in sim_folders
        folder_path = joinpath(sims_dir, folder)
        
        timeseries_path = joinpath(folder_path, "timeseries.csv")
        graph_path = joinpath(folder_path, "graph.csv")
        
        if isfile(timeseries_path)
            ts = DataFrame(CSV.File(timeseries_path))
            ts[!, :sim_number] = fill(parse(Int, folder), nrow(ts))
            append!(all_timeseries, ts)
        end
        
        if isfile(graph_path)
            g = DataFrame(CSV.File(graph_path))
            g[!, :sim_number] = fill(parse(Int, folder), nrow(g))
            append!(all_graphs, g)
        end
    end
    
    CSV.write("data/timeseries.csv", all_timeseries)
    CSV.write("data/graphs.csv", all_graphs)

    Arrow.write("data/timeseries.arrow", all_timeseries; compress=:lz4)
    Arrow.write("data/graphs.arrow", all_graphs; compress=:lz4)
    
    println("Merged timeseries and graphs saved to ./data/")

end

#==============================================================================#

function consolidate_array_time_series(sim_array)

    new_ts = DataFrame()

    nsim = length(sim_array)
    for i in 1:nsim
        sim = sim_array[i]
        append!(new_ts, sim.output[:timeseries])
    end

    return new_ts

end

#==============================================================================#

function launch_simulation(sim; dry_run=false)

    array_fn = "./data/params.csv"
    
    if typeof(sim) == String
        sim = parse(Int, sim)
    end
    
    println("Loading parameters for simulation $sim from file $array_fn")

    df = DataFrame(CSV.File(array_fn))

    params_dict_array = []
    for row in eachrow(df)
        d = Dict()
        for (name, value) in pairs(row)
            if (typeof(value) == String7) | (typeof(value) == String15)
                value = convert(String, value)
            end
            # println("name=$name value=$value")
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

# obsolete ?

# function load_simulation_array()

#     file_path = joinpath("./data", "data.jld2")
#     @load file_path sim_array

#     return sim_array

# end
