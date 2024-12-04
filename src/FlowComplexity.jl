module FlowComplexity

using JLD2
using CSV

include("Simulation.jl")

function minimal_addition_chain(n::Int)
    if n == 1
        return [1]
    end

    # Initialize the addition chain with the first element
    chain = [1]

    # A dictionary to store the minimal chain length for each number
    min_chain_length = Dict{Int, Int}()
    min_chain_length[1] = 0

    # A queue for breadth-first search
    queue = [(1, [1])]

    while !isempty(queue)
        (current, current_chain) = popfirst!(queue)

        # Try to extend the chain
        for i in 1:length(current_chain)
            next_value = current + current_chain[i]

            if next_value > n
                continue
            end

            # Check if we found a shorter chain for next_value
            if !haskey(min_chain_length, next_value) || length(current_chain) + 1 < min_chain_length[next_value]
                min_chain_length[next_value] = length(current_chain) + 1
                new_chain = vcat(current_chain, next_value)
                push!(queue, (next_value, new_chain))

                # If we reached n, return the chain
                if next_value == n
                    return new_chain
                end
            end
        end
    end

    return chain
end

function calculate_assembly_index(i)

    chain = minimal_addition_chain(i)
    len = length(chain)-1

    return len

end

function get_relative_path()

    current_dir = pwd()
    milestones_dir = joinpath(projectdir(), "milestones")
    relative_path = relpath(current_dir, milestones_dir)
    
    return relative_path

end

function write_params_file(params_array)

    array_fn = "./data/array.csv"
    CSV.write(array_fn, params_array)

end

function launch_simulation(sim)

    array_fn = "./data/array.csv"
    
    if typeof(sim) == String
        sim = parse(Int, sim)
    end
    
    println("Loading parameter file $array_fn.")
    println("Reading parameters for simulation number $sim.")

    df = DataFrame(CSV.File(array_fn))

    params_dict_array = []
    for row in eachrow(df)
        d = Dict()
        for (name, value) in pairs(row)
            if typeof(value) == String7
                value = convert(String, value)
            end
            # println("name=$name value=$value")
            merge!(d,Dict(name => value))
        end
        push!(params_dict_array,d)
    end

    params_dict = params_dict_array[sim]
    
    simulation = FlowComplexity.Simulation(; user_params = params_dict)
    
    println("Launching simulation.")
    elapsed_time = @elapsed FlowComplexity.RunSimulation(simulation)
    println("Time taken: $elapsed_time seconds")

end

function extract()

    directory = "./data/sims"
    nsim = count(isdir, readdir(directory, join=true))
    sim_array = []
    
    for i in 1:nsim
        sim_number_str = sim_number_string(i)
        sim_path = joinpath(directory, sim_number_str, "simulation.jld2")
        sim = load(sim_path)["sim"]
        push!(sim_array, sim)
    end
    
    output_file = joinpath("./data", "data.jld2")
    @save output_file sim_array
    
end

function load_simulation_array()

    file_path = joinpath("./data", "data.jld2")
    @load file_path sim_array

    return sim_array

end

end
