using JLD2

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

    array_fn = "./data/array.csv"
    println("Writing parameters for simulation array in $array_fn")
    CSV.write(array_fn, params_array)
    println("Done.")
    println("")

end

#==============================================================================#

function extract()

    directory = "./data/sims"
    nsim = count(isdir, readdir(directory, join=true))
    sim_array = []

    println("Found $nsim simulation files in data folder. Loading data...")
    for i in 1:nsim
        sim_number_str = sim_number_string(i)
        sim_path = joinpath(directory, sim_number_str, "simulation.jld2")
        sim = load(sim_path)["sim"]
        push!(sim_array, sim)
    end
    println("Loading done. Saving data array...")
    
    output_file = joinpath("./data", "data.jld2")
    @save output_file sim_array
    println("Saving done.")
    println("")
    
end

#==============================================================================#

function load_simulation_array()

    file_path = joinpath("./data", "data.jld2")
    @load file_path sim_array

    return sim_array

end

#==============================================================================#

function launch_simulation(sim)

    array_fn = "./data/array.csv"
    
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

    println("Parameters loaded.")
    
    simulation = FlowComplexity.Simulation(; user_params = params_dict)
    RunSimulation(simulation)
    
end

#==============================================================================#

