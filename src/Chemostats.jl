using Distributions
using Random
using DrWatson

mutable struct Chemostat
    ## This type contains all the information needed to run the time evolution at one spatial location
    ID::Int64 # Unique identifer within a simulation
    reaction_rate_consts::Vector{Float64}  # [constructive, destructive, outflow]
    molecules::Vector{Int64} # list storing all the molecules, non unique, duplicated molecules appear twice
    mass::Int64 # Total mass
    fixed_mass::Int64  # Target Mass for fixed flow, if 0, mass can vary
    neighbors::Vector{Int64} # What other chemostats are you connected to? (Outflowing edges)
    neighbor_flows::Vector{Float64} # Relative flows to neighbors, should sum to unity
    stabilized_integers::Vector{Int64} # List of integers that are stablized in this reactor
    savename::String # A name that persists between simulations
end

function Chemostat(ID::Int64, reaction_rate_constants::Vector{Float64}; molecules = Vector{Int64}[],
                   mass_fixed=false, neighbors = Vector{Int64}[], neighbor_flows = Vector{Float64}[],
                   stabilized_integers = Vector{Int64}[])
    # Check the savename, maybe you can just load it. 
    this_savename = make_chemostat_name(reaction_rate_constants, stabilized_integers)
    # Calculate the mass
    if length(a) > 0
        mass = sum(molecules)
    else
        mass = 0
    end

    if mass_fixed
        fixed_mass = mass
    else
        fixed_mass = 0
    end
    # Check neighbors 
    if neighbors != []
        n_neighbors = length(neighbors)
        if neighbor_flows ==[]
            neighbor_flows = [1.0/n_neighbors for n in 1:n_neighbors]
        else
            if length(neighbor_flows) != n_neighbors
                error("Number of neighbors and flows are different")
            elseif sum(neighbor_flows) != 1.0
                error("Neighbor Flows don't sum to unit")
            end
        end
    end
    return Chemostat(ID,
                     reaction_rate_constants,
                     molecules,
                     mass,
                     fixed_mass,
                     neighbors,
                     neighbor_flows,
                     stabilized_integers,
                     this_savename)
end

function parse_chemostat_name(name)


    return reaction_rate_constants, stablized_integers
end

function make_chemostat_name(reaction_rate_consts, stabilized_integers)
    name = 
    name = "reaction_rate_consts="*string(reaction_rate_consts)
    name = name*"^stablized_integers="*string(stabilized_integers)
    return name
end
function constructive_rxn(chemostat::Chemostat)
    ## Pick two random molecules from an array,
    ## join them and add the new molecule to the
    ## vector (removing the original ones)
    molecules = chemostat.molecules
    shuffle!(molecules) # Shuffle the molecules
    a = pop!(molecules) # take the first one 
    b = pop!(molecules) # and the second one 
    c = a + b # combine them (add them)
    push!(molecules, c) # add the new one to the bottom of the list 
    chemostat.molecules = molecules
    return chemostat
end

function destructive_rxn(chemostat::Chemostat)
    ## Pick a random molecule (of length >1)
    ## split at a random point, add both fragments back to 
    molecules = chemostat.molecules
    shuffle!(molecules) # Shuffle the molecules
    
    big_moles = filter(x -> x > 1, molecules) # Ignore 1s
    if length(big_moles) > 0 # If there are any molecules left
        a = pop!(big_moles) # Grab one 
        mole_index = findfirst(x -> x == a,molecules) # pop it out of the molecule list
        a = popat!(molecules, mole_index)
        b = sample(1:(a-1)) # find a bond to split it at 
        c = a - b
        push!(molecules, b) # push both 
        push!(molecules, c)
    end
    chemostat.molecules = molecules
    return chemostat
end

function outflow_rxn(chemostat)
    ## Pick two random molecules from an array, 
    ## join them and add the new molecule to the 
    ## vector (removing the original ones)
    molecules = chemostat.molecules
    neighbors = chemostat.neighbors
    neighbor_weights = chemostat.neighbor_flows
    shuffle!(molecules) # Shuffle the molecules
    a = pop!(molecules) # take the first one 
    if neighbors != []
        neighbor = sample(neighbors, Weights(neighbor_weights))
        outflow_direction = Dict(neighbor => a)
    else
        outflow_direction = Dict{Int64,String}()
    end 
    chemostat.molecules = molecules
    return chemostat, outflow_direction
end

function calc_mass(chemostat)
    ## Calculate total mass of molecules
    ## mass is just length 
    molecules = chemostat.molecules
    mass = sum(molecules)
    return mass
end

function calc_propensities(chemostat)
    ## Calculate the propensities for each reaction time
    
    # Constructive
    # depends on n*n-1 (number of both pairwise interactions)
    n = length(chemostat.molecules)
    Cp = chemostat.reaction_rate_consts[1] * n * (n -1)
    
    # Destructive 
    # Depends on number of molecules 
    Cd = chemostat.reaction_rate_consts[2] * n

    # Outflow 
    # Depends on number of molecules 
    Co = chemostat.reaction_rate_consts[3] * n

    propensities = [Cp, Cd, Co]

    return propensities # Constructive, destructive, outflow
end