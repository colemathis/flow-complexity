include("Chemostats.jl")

using Graphs

mutable struct Ensemble
    ## This type contains all the information needed to run the time evolution across spatial locations
    reactor_ids::Array{Int64,1} # What are the unique reactor ids
    reactors::Array{Chemostat,1}
    ensemble_graph::DiGraph # How are the reactors connected (using Graph type from Graphs.jls)
    inflow_ids::Array{Int64,1} # Which reactors are sources? 
end

function calc_next_rxn_times(all_propensities, tau)
    ### Calculate the next reaction time for each set of propensities
    # Update time 
    n = length(all_propensities)
    total_p = sum.(all_propensities)
    print(total_p)
    tau_steps = -log.(rand(n))./total_p
    reactor_steps = collect(zip(tau_steps .+ tau, 1:n))
    reactor_steps = sort(reactor_steps, by = x->x[1])
    println(reactor_steps)
    return reactor_steps

end


function make_line_reactors(n, rxn_rates, inflow_mass, initial_mass)

    ensemble_graph = path_digraph(n)
    reactors = Array{Chemostat,1}()
    
    for i in 1:n
        if i ==1
            inflow = inflow_mass 
            molecules = repeat([1], initial_mass)
            mass= initial_mass
        else
            inflow = 0
            molecules = Array{Int64,1}()
            mass = 0
        end
        neighbor = neighbors(ensemble_graph, i)
        neighbor_flow = [1.0]
        this_reactor = Chemostat(i, neighbor, neighbor_flow, molecules, rxn_rates, mass, inflow) 
        push!(reactors, this_reactor)
    end

    line_ensemble = Ensemble(collect(1:n), reactors, ensemble_graph, [1])

    return line_ensemble
end
