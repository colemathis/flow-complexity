include("Chemostats.jl")
incldue("TimeEvolve.jl")

using Graphs

mutable struct Ensemble
    ## This type contains all the information needed to run the time evolution across spatial locations
    reactor_ids::Array{Int64,1} # What are the unique reactor ids
    ensemble_graph::DiGraph # How are the reactors connected (using Graph type from Graphs.jls)
    inflow_ids::Array{Int64,1} # Which reactors are sources? 
end



function make_line_reactors(n, rxn_rates, inflow_mass)



    return line_ensemble
end
