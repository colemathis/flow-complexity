#==============================================================================#
# IMPORTS
#==============================================================================#

import Graphs
import StatsBase

#==============================================================================#
# DATA TYPES
#==============================================================================#

mutable struct Ensemble

    reactor_ids     ::Array{Int64,1}
    reactors        ::Array{Chemostat,1}
    ensemble_graph  ::Graphs.DiGraph
    inflow_ids      ::Array{Int64,1}

end

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function Ensemble(N_reactors        ::Int64, 
                  graph_type        ::String, 
                  N_sources         ::Int64, 
                  randomize_edges   ::Bool,
                  mass              ::Int64;
                  chemostat_list    = Vector{Dict{String,Any}}[]
                  )

    if graph_type == "ER"

        ensemble_graph = Graphs.erdos_renyi(N_reactors, 4*N_reactors, is_directed =true)
        inflow_ids, ensemble_graph = find_inflow_nodes(ensemble_graph, N_sources)
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    elseif graph_type == "BA"

        ensemble_graph = Graphs.barabasi_albert(N_reactors, 2, is_directed=true) # Check why is this 2?
        inflow_ids, ensemble_graph  = find_inflow_nodes(ensemble_graph, N_sources)
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    elseif graph_type == "regular"

        ensemble_graph = Graphs.random_regular_digraph(N_reactors, 4)
        inflow_ids, ensemble_graph = find_inflow_nodes(ensemble_graph, N_sources)
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)
    
    elseif graph_type == "lattice"

        (round(sqrt(N_reactors)))^2 == N_reactors || error("lattice graph requires N_reactors to be a perfect square")

        ensemble_graph = lattice_digraph(N_reactors)
        inflow_ids, ensemble_graph  = find_inflow_nodes(ensemble_graph, N_sources)
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    elseif graph_type == "lattice-2way"

        (round(sqrt(N_reactors)))^2 == N_reactors || error("lattice graph requires N_reactors to be a perfect square")
        ensemble_graph = lattice_digraph_bidirectionnal(N_reactors)

        # optionally randomize edges while preserving degree distribution
        randomize_edges && (ensemble_graph = randomize_graph!(ensemble_graph, 10 * Graphs.ne(ensemble_graph)))

        # inflow_ids, ensemble_graph  = find_inflow_nodes(ensemble_graph, N_sources)
        inflow_ids = [1]
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    elseif graph_type == "line"
        # enforce a single source for a line graph
        N_sources == 1 || error("line graph requires exactly one source node")

        ensemble_graph = Graphs.path_digraph(N_reactors)
        inflow_ids = [1]
        chemostat_specs = StatsBase.sample(chemostat_list, N_reactors)
        chemostats = chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    end

    # return the Ensemble
    return Ensemble(collect(1:N_reactors),
                    chemostats,
                    ensemble_graph,
                    inflow_ids,
                    )

end

#==============================================================================#

function lattice_digraph(N)

    # construct a simple digraph with N vertices (and 0 edges)
    L = Graphs.SimpleDiGraph(N)

    # link the vertices together (add edges) to form the lattice
    n = sqrt(N)
    rows = []
    current = 0
    for i in 1:n
        for j in 0:(n-1)
            current_node = i + n*j  
            right_neighbor = current_node + 1
            down_neighbor = current_node + n 
            if i < n
                Graphs.add_edge!(L,(current_node, right_neighbor))
            end
            Graphs.add_edge!(L, (current_node, down_neighbor))
        end
    end

    # return the graph
    return L 

end

#==============================================================================#

function lattice_digraph_bidirectionnal(N)

    # construct a simple digraph with N vertices (and 0 edges)
    L = Graphs.SimpleDiGraph(N)

    # link the vertices together (add edges) to form the lattice
    n = sqrt(N)
    rows = []
    current = 0
    for i in 1:n
        for j in 0:(n-1)
            current_node = i + n*j  
            right_neighbor = current_node + 1
            down_neighbor = current_node + n 
            if i < n
                Graphs.add_edge!(L,(current_node, right_neighbor))
                Graphs.add_edge!(L,(right_neighbor, current_node))
            end
            Graphs.add_edge!(L, (current_node, down_neighbor))
            Graphs.add_edge!(L, (down_neighbor, current_node))
        end
    end

    # return the graph
    return L 

end

#==============================================================================#

function find_inflow_nodes(graph, n_sources)

    error("uh-oh, find_inflow_nodes is not implemented yet")

    # Find the nodes that have no in-edges
    total_edges = length(Graphs.edges(graph))

    # get the in-degree for each node
    in_degrees = Graphs.indegree(graph, Graphs.vertices(graph))
    # determine which node has no in-edge
    no_in_nodes = [i for i in 1:length(in_degrees) if in_degrees[i] == 0]
    # sort nodes by in-degree
    sorted_nodes = sort!(collect(Graphs.vertices(graph)), by = x->Graphs.indegree(graph,x))
    # determine how many nodes have no in-edge
    num_no_in = length(no_in_nodes)

    #p1: not too sure I understand what the rest is about
    if num_no_in >= n_sources
        # You're good just pick a sample of those with 0 in in_degrees
        no_in_nodes = [i for i in 1:length(in_degrees) if in_degrees[i] == 0]
        source_list = sort!(sample(no_in_nodes, n_sources, replace=false))

    else
        # You're going to need to modify the graph to make it work
        source_list = Int64[]
        # Remove the number of nodes you need 
        nodes_to_remove = sorted_nodes[1:n_sources]
        Graphs.rem_vertices!(graph, nodes_to_remove)
        # Get the remaining nodes 
        other_nodes = collect(Graphs.vertices(graph))
        # Make the remaining graph a completely connected graph
        components = Graphs.connected_components(graph)
        while length(components) > 1
            Graphs.add_edge!(graph, rand(components[1]), rand(components[2]))
            components = Graphs.connected_components(graph)
        end
        # Now we're going to connect the source nodes to random other nodes
        random_connections = sample(Graphs.vertices(graph), n_sources)
        for i in 1:n_sources
            Graphs.add_vertex!(graph)
            new_node = length(Graphs.vertices(graph))
            Graphs.add_edge!(graph, new_node, random_connections[i])
            push!(source_list, new_node)
        end

    end

    # Tidy up by making sure you've got the right number of edges
    new_edge_count = length(Graphs.edges(graph))
    total_edges - new_edge_count
    while total_edges - new_edge_count > 0
        source_node = sample(Graphs.vertices(graph))
        destination_node = sample(other_nodes)
        Graphs.add_edge!(graph, source_node, destination_node)
        new_edge_count = length(Graphs.edges(graph))
    end

    return source_list, graph

end

#==============================================================================#

function chemostats_from_specs(ensemble_graph, chemostat_specs, inflow_ids, mass)

    # Make the chemostats with the right neighbors 

    chemostat_list = Chemostat[]
    reactors = collect(Graphs.vertices(ensemble_graph))

    for r in reactors

        specs = StatsBase.sample(chemostat_specs)

        reaction_rate_constants = specs["reaction_rate_constants"]

        # Return a list of outneighbors for node r
        # i.e., the nodes that vertices from r point to
        these_neigbors = Graphs.neighbors(ensemble_graph, r)

        if r in inflow_ids

            molecules = repeat([1], mass)
            mass_fixed= true

            this_chemostat = Chemostat(r, reaction_rate_constants,
                                        molecules = molecules,
                                        mass_fixed = mass_fixed,
                                        neighbors = these_neigbors)
            push!(chemostat_list, this_chemostat)

        else

            this_chemostat = Chemostat(r, reaction_rate_constants, 
                                       neighbors = these_neigbors)
            push!(chemostat_list, this_chemostat)
        end

    end

    return chemostat_list

end

#==============================================================================#

function gen_chemostat_spec_list(N_chemostats   ::Int64,
                                 forward_rate   ::Float64,
                                 backward_rate  ::Float64,
                                 outflow        ::Float64
                                 )

    # dictionary template for chemostat specs
    single_dict = Dict("reaction_rate_constants" => [forward_rate, backward_rate, outflow])

    # create the array of chemostat specs
    spec_list = [single_dict]

    # fill the array with identical copies of the template
    for n in 1:(N_chemostats-1)
        new_spec = copy(single_dict)
        push!(spec_list, new_spec)
    end

    # return the list of chemostat specs
    return spec_list 
    
end

#==============================================================================#

function randomize_graph!(g::Graphs.SimpleDiGraph, nswap::Int)
    # randomize_directed_preserving_degree

    edges = collect(Graphs.edges(g))
    swaps = 0
    tries = 0

    while swaps < nswap && tries < 100 * nswap
        tries += 1

        e1, e2 = rand(edges, 2)
        u1, v1 = Graphs.src(e1), Graphs.dst(e1)
        u2, v2 = Graphs.src(e2), Graphs.dst(e2)

        # Skip if any overlap between the four nodes
        if length(Set([u1, v1, u2, v2])) < 4
            continue
        end

        # Check that mirrored edges exist
        if !(Graphs.has_edge(g, v1, u1) && Graphs.has_edge(g, v2, u2))
            error("Mirrored edges missing: cannot swap ($u1->$v1) and ($u2->$v2) without ($v1->$u1) and ($v2->$u2)")
        end

        # Proposed new edges: (u1→v2) and (u2→v1) and their mirrors (v2→u1) and (v1→u2)
        if Graphs.has_edge(g, u1, v2) || Graphs.has_edge(g, u2, v1) || Graphs.has_edge(g, v2, u1) || Graphs.has_edge(g, v1, u2) ||
           u1 == v2 || u2 == v1 || v2 == u1 || v1 == u2
            continue
        end

        # Perform the swap for both edges and their mirrors
        Graphs.rem_edge!(g, u1, v1)
        Graphs.rem_edge!(g, v1, u1)
        Graphs.rem_edge!(g, u2, v2)
        Graphs.rem_edge!(g, v2, u2)

        Graphs.add_edge!(g, u1, v2)
        Graphs.add_edge!(g, v2, u1)
        Graphs.add_edge!(g, u2, v1)
        Graphs.add_edge!(g, v1, u2)

        # Update edge list and count
        edges = collect(Graphs.edges(g))
        swaps += 1
    end

    return g
end

#==============================================================================#
# END OF FILE
#==============================================================================#
