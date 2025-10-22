#==============================================================================#
# IMPORTS
#==============================================================================#


#==============================================================================#
# FUNCTIONS
#==============================================================================#

function minimal_addition_chain(n::Int)

    error("uh-oh, minimal_addition_chain is not implemented yet")

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

    error("uh-oh, calculate_assembly_index is not implemented yet")

    chain = minimal_addition_chain(i)
    len = length(chain)-1

    return len

end

#==============================================================================#
# END OF FILE
#==============================================================================#

