using Distributions 
using Random

mutable struct Chemostat
    ## This type contains all the information needed to run the time evolution at one spatial location
    ID::Int64 # Unique identifer
    neighbors::Array{Int64,1} # What other chemostats are you connected to? (Outflowing edges)
    neighbor_flows::Array{Float64,1} # Relative flows to neighbors, should sum to unity 
    molecules::Array{String,1} # list storing all the molecules, non unique, duplicated molecules appear twice
    reaction_probs::Array{Float64,1}  # [constructive, destructive, outflow]
    mass::Int64 # Total mass
    mass_fixed::Bool # Whether the total mass should be fixed 
end

function constructive_rxn(molecules)
    ## Pick two random molecules from an array, 
    ## join them and add the new molecule to the 
    ## vector (removing the original ones)
    shuffle!(molecules) # Shuffle the molecules
    a = pop!(molecules) # take the first one 
    b = pop!(molecules) # and the second one 
    c = a*b # combine them (* in julia concatentates strings)
    push!(molecules, c) # add the new one to the bottom of the list 
    
    return molecules
end

function destructive_rxn(molecules)
    ## Pick a random molecule (of length >1)
    ## split at a random point, add both fragments back to 
    shuffle!(molecules) # Shuffle the molecules
    big_moles = filter(x -> length(x) > 1, molecules) # Ignore monomers
    if length(big_moles) > 0 # If there are any molecules left
        a = pop!(big_moles) # Grab one 
        mole_index = findfirst(x -> x == a,molecules) # pop it out of the molecule list
        a = popat!(molecules, mole_index)
        break_point = sample(1:(length(a)-1)) # find a bond to split it at 
        b = a[1:break_point] # Left side
        c = a[(break_point + 1):end] # right side 
        push!(molecules, b) # push both 
        push!(molecules, c)
    end
    
    return molecules
end