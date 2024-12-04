using Distributions
using Random
using DrWatson

#==============================================================================#
# DATA TYPES
#==============================================================================#

mutable struct Chemostat

    ID::Int64
    reaction_rate_consts::Vector{Float64}
    molecules::Vector{Int64}
    mass::Int64
    fixed_mass::Int64
    neighbors::Vector{Int64}
    neighbor_flows::Vector{Float64}

end

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function Chemostat(ID::Int64,
                   reaction_rate_constants::Vector{Float64};
                   molecules=Vector{Int64}[], #p3 why a semi-colon here?
                   mass_fixed=false,
                   neighbors=Vector{Int64}[],
                   neighbor_flows=Vector{Float64}[]
)

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
        if neighbor_flows == []
            neighbor_flows = [1.0 / n_neighbors for n in 1:n_neighbors]
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
        neighbor_flows)

end;

#==============================================================================#

function calc_mass(chemostat)

    molecules = chemostat.molecules
    mass = sum(molecules)

    return mass

end

#==============================================================================#

function calc_ones(chemostat)

    # get the list of molecules of length 1
    molecules = chemostat.molecules
    all_ones = [m for m in molecules if m == 1]

    # return the total number of molecules of length 1
    return sum(all_ones)

end

#==============================================================================#

function calc_propensities(chemostat)

    # Constructive - depends on n*n-1 (pairwise reaction)
    n = length(chemostat.molecules)
    Cp = chemostat.reaction_rate_consts[1] * n * (n - 1)

    # Destructive - depends on n
    sum_ones = calc_ones(chemostat)
    Cd = chemostat.reaction_rate_consts[2] * (n - sum_ones)

    # Outflow - depends on n
    Co = chemostat.reaction_rate_consts[3] * n

    propensities = [Cp, Cd, Co]

    return propensities

end
