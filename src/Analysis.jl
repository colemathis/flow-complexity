

#==============================================================================#
# FUNCTIONS
#==============================================================================#

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

#==============================================================================#

function calculate_assembly_index(i)

    chain = minimal_addition_chain(i)
    len = length(chain)-1

    return len

end

#==============================================================================#

# obsolete?

# function calculate_sim_populations(sim)

#     nt = length(sim.output[:timestamps])
#     nchem = sim.params[:N_reactors]

#     # find max specie
#     max_specie = 0
#     for i in 1:nt
#         for j in 1:nchem
#             mols = sim.output[:populations][i][j]
#             if mols != []
#                 m = maximum(mols)
#             else
#                 m = 0
#             end
#             if m > max_specie
#                 max_specie = m
#             end
#         end
#     end
    
#     pop = zeros(Int, nt, nchem, max_specie)

#     for i in 1:nt
#         for j in 1:nchem
#             mols = sim.output[:populations][i][j]
#             species = unique(mols)
#             for s in species
#                 pop[i,j,s] = count(x -> x == s, mols)
#             end
#         end
#     end
    
#     return pop

# end

#==============================================================================#

# obsolete?

# function calculate_array_populations(sim_array)

#     nsim = length(sim_array)
#     sim = sim_array[1]
#     nt = length(sim.output[:timestamps])
#     nchem = sim.params[:N_reactors]

#     # find max specie
#     max_specie = 0
#     for i in 1:nsim
#         for j in 1:nt
#             for k in 1:nchem
#                 mols = sim_array[i].output[:populations][j][k]
#                 if mols != []
#                     m = maximum(mols)
#                 else
#                     m = 0
#                 end
#                 if m > max_specie
#                     max_specie = m
#                 end
#             end
#         end
#     end
    
#     pop = zeros(Int, nsim, nt, nchem, max_specie)

#     for i in 1:nsim
#         for j in 1:nt
#             for k in 1:nchem
#                 mols = sim_array[i].output[:populations][j][k]
#                 species = unique(mols)
#                 for s in species
#                     pop[i,j,k,s] = count(x -> x == s, mols)
#                 end
#             end
#         end
#     end
    
#     return pop
    
# end

#==============================================================================#

# obsolete

# function create_populations_dataframe(sim_array)

#     pop = DataFrame()

#     nsim = length(sim_array)
#     for i in 1:nsim
#         sim = sim_array[i]
#         sim_number = sim.params[:sim_number]
#         nt = length(sim.output[:timestamps])
#         nchem = sim.params[:N_reactors]
#         for j in 1:nt
#             ts = sim.output[:timestamps][j]
#             for k in 1:nchem
#                 mols = sim.output[:populations][j][k]
#                 unique_mols = unique(mols)
#                 for m in unique_mols
#                     f = count(x -> x == m, mols)
#                     r = (sim_number = sim_number,
#                          time = ts,
#                          chemostat_id = k,
#                          integer = m,
#                          frequency = f)
#                     push!(pop, r)
#                 end
#             end
#         end
#     end

#     return pop
    
# end

#==============================================================================#

