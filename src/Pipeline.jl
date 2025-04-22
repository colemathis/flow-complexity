using JLD2
using CSV

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function get_relative_path()

    current_dir = pwd()
    milestones_dir = joinpath(projectdir(), "milestones")
    relative_path = relpath(current_dir, milestones_dir)
    
    return relative_path

end

#==============================================================================#

function write_params_file(params_array)

    array_fn = "./data/params.csv"
    println("Writing parameters for simulation array in $array_fn")
    CSV.write(array_fn, params_array)
    println("Done.")
    println("")

end

#==============================================================================#

function extract()

    directory = "./data/sims"
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

# obsolete ?

# function load_simulation_array()

#     file_path = joinpath("./data", "data.jld2")
#     @load file_path sim_array

#     return sim_array

# end

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
    
end

#==============================================================================#
