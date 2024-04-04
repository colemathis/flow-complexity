"""

##########################################################################

CHEMOSTAT



##########################################################################

"""

using Distributions
using Random
using DrWatson

"""

##########################################################################
CHEMOSTAT STRUCT

This type contains all the information needed to run the time evolution at one spatial location.
##########################################################################

    ID                      (int)           = (see Simulation struct)
    reaction_rate_consts    (float vector)  = (see Simulation struct)
    molecules               (int vector)    = list storing all the molecules, non unique, duplicated molecules appear twice
    mass                    (int)           = (see Simulation struct)
    fixed_mass              (int)           = target mass for fixed flow, if 0 mass can vary
    neighbors               (int vector)    = what other chemostats are you connected to? (outflowing edges)
    neighbor_flows          (float vector)  = relative flows to neighbors, should sum to unity
    stabilized_integers     (int vector)    = (see Simulation struct) #p4 to be implemented

##########################################################################

"""

#p3: fix duplicate definition, related to using include instead of import/using
#p3: see https://docs.julialang.org/en/v1/manual/modules/

mutable struct Chemostat

    ID                      ::Int64
    reaction_rate_consts    ::Vector{Float64}
    molecules               ::Vector{Int64}
    mass                    ::Int64
    fixed_mass              ::Int64
    neighbors               ::Vector{Int64}
    neighbor_flows          ::Vector{Float64}
    stabilized_integers     ::Vector{Int64}

end

"""

##########################################################################
GENERATOR FOR CHEMOSTAT STRUCT

Parametrizes the Chemostat struct.
##########################################################################

Input:

    ID                          (int)           = (see Simulation struct)
    reaction_rate_constants     (float vector)  = (see Simulation struct)
    molecules                   (int vector)    = (see Chemostat struct)
    mass_fixed                  (boolean)       = (see Chemostat struct)
    neighbor_flows              (float vector)  = (see Chemostat struct)
    neighbor_flows              (float vector)  = (see Chemostat struct)
    stabilized_integers         (int vector)    = (see Simulation struct)

Output:
    
    Chemostat                   (Chemostat)     = resulting chemostat

##########################################################################

"""

#p3: replaced semi-colon in the following...
function Chemostat(ID                       ::Int64,
                   reaction_rate_constants  ::Vector{Float64}; 
                   molecules                = Vector{Int64}[], #p3 why a semi-colon here?
                   mass_fixed               = false, 
                   neighbors                = Vector{Int64}[], 
                   neighbor_flows           = Vector{Float64}[],
                   stabilized_integers      = Vector{Int64}[]
                   )

#p3: with a comma, to get the fully mixed code to work, but this breaks the distributed one !
# function Chemostat(ID::Int64, reaction_rate_constants::Vector{Float64}, molecules = Vector{Int64}[],
#                     mass_fixed=false, neighbors = Vector{Int64}[], neighbor_flows = Vector{Float64}[],
#                     stabilized_integers = Vector{Int64}[])

    # Check the savename, maybe you can just load it. 
    #p3 what

    # Calculate the total mass from the list of molecules
    if length(molecules) > 0
        mass = sum(molecules)
    else
        mass = 0
    end

    # verify if the mass is fixed
    if mass_fixed
        fixed_mass = mass
    else
        fixed_mass = 0
    end

    # if there’s a list of neighbors, define the flow to be equally distributed (if no flow is defined)
    # or check if flows sum to unity
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

    # return the Chemostat struct
    return Chemostat(ID,
                     reaction_rate_constants,
                     molecules,
                     mass,
                     fixed_mass,
                     neighbors,
                     neighbor_flows,
                     stabilized_integers)
end;

"""

##########################################################################
PERFORM ONE CONSTRUCTIVE REACTION.

Pick two random molecules, join them and remove the original ones.
##########################################################################

Input:

    chemostat   (Chemostat) = chemostat before the reaction occurs

Output:
    
    chemostat   (Chemostat) = chemostat after the reaction occurred

##########################################################################

"""

function constructive_rxn(chemostat::Chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)
    shuffle!(molecules)

    # combine two molecules
    a = pop!(molecules)
    b = pop!(molecules)
    c = a + b

    # add the new one to the end of the list
    push!(molecules, c)

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat

end

"""

##########################################################################
PERFORM ONE DESTRUCTIVE REACTION.

Pick a random molecule, split it and add back fragments to the list.
##########################################################################

Input:

    chemostat   (Chemostat) = chemostat before the reaction occurs

Output:
    
    chemostat   (Chemostat) = chemostat after the reaction occurred

##########################################################################

"""

function destructive_rxn(chemostat::Chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)
    shuffle!(molecules) # Shuffle the molecules
    
    # ignore molecules of length 1
    big_moles = filter(x -> x > 1, molecules)

    if length(big_moles) > 0

        # pick one molecule and remove it from the list
        a = pop!(big_moles)
        mole_index = findfirst(x -> x == a,molecules)
        a = popat!(molecules, mole_index)

        # split it at a random point
        b = sample(1:(a-1))
        c = a - b

        # add back both to the list
        push!(molecules, b)
        push!(molecules, c)

    end

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat

end

"""

##########################################################################
PERFORM ONE OUTFLOW REACTION.

Pick two random molecules, join them, add the new molecule to the list 
and remove the original one. #p1: to complete
##########################################################################

Input:

    chemostat           (Chemostat) = chemostat before the reaction occurs

Output:
    
    chemostat           (Chemostat) = chemostat after the reaction occurred
    outflow_direction   (dict)      = #p1 to complete

##########################################################################

"""

function outflow_rxn(chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)

    # get the list of neighbors and their inflow to this chemostat
    neighbors = chemostat.neighbors
    neighbor_weights = chemostat.neighbor_flows

    # shuffle the molecules and pick one
    shuffle!(molecules)
    a = pop!(molecules)

    # pick a neighbor and ??? #p1: not sure what this is about
    if neighbors != []
        neighbor = sample(neighbors, Weights(neighbor_weights))
        outflow_direction = Dict(neighbor => a)
    else
        outflow_direction = Dict{Int64,String}()
    end 

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat, outflow_direction

end

"""

##########################################################################
CALCULATE THE TOTAL MASS OF MOLECULES

Calculate the total mass of molecules (=length).
##########################################################################

Input:

    chemostat   (Chemostat) = chemostat

Output:
    
    mass        (int)       = total mass

##########################################################################

"""

function calc_mass(chemostat)

    # get the list of molecules
    molecules = chemostat.molecules

    # sum their length
    mass = sum(molecules)

    return mass

end

"""

##########################################################################
CALCULATE THE NUMBER OF ONES.

Calculate the number of molecules of length 1.
##########################################################################

Input:

    chemostat   (Chemostat)     = chemostat

Output:
    
    (implicit)  (int)           = number of molecules of length 1

##########################################################################

"""

function calc_ones(chemostat)

    # get the list of molecules of length 1
    molecules = chemostat.molecules
    all_ones = [m for m in molecules if m == 1]

    # return the total number of molecules of length 1
    return sum(all_ones)

end

"""

##########################################################################
CALCULATE PROPENSITIES

Calculate the propensities to react for each type of reaction 
(constructive, destructive, outflow).
##########################################################################

Input:

    chemostat       (Chemostat)     = chemostat

Output:
    
    propensities    (float vect)    = propensities

##########################################################################

"""

function calc_propensities(chemostat)
    
    # Constructive - depends on n*n-1 (pairwise reaction)
    n = length(chemostat.molecules)
    Cp = chemostat.reaction_rate_consts[1] * n * (n -1)
    
    # Destructive - depends on n
    Cd = chemostat.reaction_rate_consts[2] * n

    # Outflow - depends on n
    Co = chemostat.reaction_rate_consts[3] * n

    propensities = [Cp, Cd, Co]

    return propensities
    
end