using DrWatson
using FileIO
using Dates
using DataFrames
using Graphs
using GraphIO
using BSON

function save_data(sim, run_parameters)
    # Save all the data in a new directory with an updated simulation number
    # Save Formats: TODO: Update this shit
    #   - timeseries_data as a csv
    #   - run parameters as a csv
    #   - graph as an edgelist (.csv)
    #   - Simulation object as .BSON 

    sim_number = string(sim.sim_number)
    
    # Save time series
    time_series_df = convert_timeseries_to_tidy_df(sim.time_evolution)
    save(datadir("sims", sim_number, "timeseries.csv"), time_series_df)
    # Save Run parameters
    save(datadir("sims", sim_number, "parameters.csv"), DataFrame(run_parameters))
    # Save reactor graph
    reactors = sim.ensemble
    edge_list = generate_edge_list(reactors)
    save(datadir("sims", sim_number, "graph.csv"), edge_list)
    # Save Simulation object 
    save(datadir("sims", sim_number, "simulation.bson"), Dict(:sim =>sim))
    println("Data Saved")
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

function convert_timeseries_to_tidy_df(timeseries)
    recorded_vars = [k for k in keys(timeseries[1])]
    times = [t for t in keys(timeseries[1][recorded_vars[1]])]
    reactors = collect(keys(timeseries))
    data = []
    for r in reactors
        for t in times
            for var in recorded_vars
                if var == :complete_timeseries
                    if timeseries[r][:complete_timeseries][t] != []
                        time_counts = countmap(timeseries[r][:complete_timeseries][t])
                        for (v,c) in time_counts
                            push!(data, Dict("reactor"=> r,"time"=> t, "variable"=>string(v), "value"=> c))
                        end
                    end
                else
                    push!(data, Dict("reactor"=> r,"time" => t, "variable" => String(var), "value"=> get(timeseries[r][var],t,0) ))
                end
            end
        end
    end

    tidy_df = DataFrame(time = map(x -> x["time"], data),
                        reactor = map(x -> x["reactor"], data),
                        variable = map(x -> x["variable"], data),
                        value = map(x -> x["value"], data))

    return tidy_df
end