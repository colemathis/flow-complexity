using DrWatson
using JLD2
using FileIO
using Dates
using DataFrames
using Graphs
using GraphIO

function save_data(timeseries_data, run_parameters, reactors)
    # Save all the data in a new directory with an updated simulation number
    # Save Formats:
    #   - timeseries_data as a bson 
    #   - run parameters as a json
    #   - graph as an edgelist (.txt)
    sim_number = string(get_sim_number())
    save(datadir("sims", sim_number, "timeseries.bson"), timeseries_data)
    save(datadir("sims", sim_number, "parameters.csv"), DataFrame(run_parameters))
    edge_list = generate_edge_list(reactors)
    save(datadir("sims", sim_number, "graph.csv"), edge_list)
    # fname = datadir("sims", savename(parameters, "bson", connector = "^"))
    # save(fname, data)
end

function generate_edge_list(reactors)
    # Generate an edge list as a DataFrame for easy saving
    graph = reactors.ensemble_graph
    sources = []
    destinations = []
    source_inflow = []
    these_edges = edges(graph)
    for e in these_edges
        push!(sources, src(e))
        push!(destinations, dst(e))
        if src(e) in reactors.inflow_ids
            push!(source_inflow, true)
        else
            push!(source_inflow, false)
        end
    end 
    edge_dict = Dict("sources" => sources,
                     "destinations" => destinations,
                     "source_inflow" => source_inflow)
    edge_df = DataFrame(edge_dict)
    return edge_df
end

function get_sim_number()
    # Find the maximum number used as a directory in the sims folder. Add one to that, and return
    # it. This will be the new simulation number
    all_sim_numbers = [s for s in (tryparse.(Int,readdir(datadir("sims")))) if s !== nothing ]
    sim_number = maximum(all_sim_numbers) + 1
    return sim_number
end

