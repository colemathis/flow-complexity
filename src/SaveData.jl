"""

##########################################################################

SAVEDATA



##########################################################################

"""

using DrWatson
using FileIO
using Dates
using DataFrames
using Graphs
using GraphIO
using BSON

"""

##########################################################################
SAVE DATA

Save all the data in a new directory with an updated simulation number.
##########################################################################

Input:

    sim             (Simulation)    = simulation to save
    run_parameters  (dict)          = parameters to save

Output:

    (none)

##########################################################################

"""

# #p4: Save Formats
#   - timeseries_data as a csv
#   - run parameters as a csv
#   - graph as an edgelist (.csv)
#   - Simulation object as .BSON 

function save_data(sim, 
                   run_parameters
                   )

    # get the sim number
    sim_number = string(sim.sim_number)
    
    if sim.save_time_series
        # #p1: what’s the conversion below?
        time_series_df = convert_timeseries_to_tidy_df(sim.time_evolution)
        # save the time series
        # save(datadir("sims", sim_number, "timeseries.csv"), time_series_df)
        # save(datadir(sim.save_name, "timeseries.csv"), time_series_df)
        fn = joinpath(sim.save_name, "timeseries.csv")
        save(datadir(fn, time_series_df)
    end

    if sim.save_parameters
        # save parameters
        # save(datadir("sims", sim_number, "parameters.csv"), DataFrame(run_parameters))
        # save(datadir(sim.save_name, "parameters.csv"), DataFrame(run_parameters))
        fn = joinpath(sim.save_name, "parameters.csv")
        save(fn, DataFrame(run_parameters))
    end

    if sim.save_graph
        # save reactor graph
        reactors = sim.ensemble
        edge_list = generate_edge_list(reactors)
        # save(datadir("sims", sim_number, "graph.csv"), edge_list)
        # save(datadir(sim.save_name, "graph.csv"), edge_list)
        fn = joinpath(sim.save_name, "graph.csv")
        save(fn, edge_list)
    end

    if sim.save_simulation
        # save simulation object 
        # save(datadir("sims", sim_number, "simulation.bson"), Dict(:sim =>sim))
        # save(datadir(sim.save_name, "simulation.bson"), Dict(:sim =>sim))
        # save(datadir(sim.save_name, "simulation.jld2"), Dict("sim" =>sim))
        fn = joinpath(sim.save_name, "simulation.jld2")
        save(fn, Dict("sim" =>sim))
    end

    println("Data Saved")

end

"""

##########################################################################
GENERATE EDGE LIST AS DATAFRAME

Generate an edge list as a DataFrame for easy saving.
##########################################################################

Input:

    reactors    (Ensemble)      = ensemble to convert

Output:

    edge_df     (DataFrame)     = dataframe containing edges

##########################################################################

"""

function generate_edge_list(reactors)

    # get the graph and edges of the ensemble
    graph = reactors.ensemble_graph
    these_edges = edges(graph)

    # define the named vectors of the dataframe
    sources = []
    destinations = []
    source_inflow = []

    # fill them
    for e in these_edges
        push!(sources, src(e))
        push!(destinations, dst(e))
        if src(e) in reactors.inflow_ids
            push!(source_inflow, true)
        else
            push!(source_inflow, false)
        end
    end 

    # create the dataframe
    edge_dict = Dict("sources" => sources,
                     "destinations" => destinations,
                     "source_inflow" => source_inflow)
    edge_df = DataFrame(edge_dict)

    # return it
    return edge_df

end

"""

##########################################################################
GET SIM NUMBER

Find the maximum number used as a directory in the sims folder. Add one to that, and return
it. This will be the new simulation number.
##########################################################################

Input:

    (none)

Output:

    sim_number  (int)   = next sim number available

##########################################################################

"""

function get_sim_number()

    # read the sim numbers in the "sims" directory
    all_sim_numbers = [s for s in (tryparse.(Int,readdir(datadir("sims")))) if s !== nothing ]

    # add one to the highest sim number
    sim_number = maximum(all_sim_numbers) + 1

    # and return that
    return sim_number

end

"""

##########################################################################
CONVERT TIME SERIES TO DATAFRAME

Convert the time series as a DataFrame for easy saving.
##########################################################################

Input:

    timeseries  (dict)          = time series data

Output:

    tidy_df     (DataFrame)     = dataframe containing the time series

##########################################################################

"""

function convert_timeseries_to_tidy_df(timeseries)

    # get the recorded variables of the timeseries dict
    recorded_vars = [k for k in keys(timeseries[1])]

    # get the times
    times = [t for t in keys(timeseries[1][recorded_vars[1]])]

    # get the reactors
    reactors = collect(keys(timeseries))

    # create a data array and fill, for every reactor and time frame, the variable and its value
    data = []
    for r in reactors
        for t in times
            for var in recorded_vars
                if var == :complete_timeseries #p3 not too sure what’s the difference here
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

    # wrap this up in a dataframe
    tidy_df = DataFrame(time = map(x -> x["time"], data),
                        reactor = map(x -> x["reactor"], data),
                        variable = map(x -> x["variable"], data),
                        value = map(x -> x["value"], data))

    # and return it
    return tidy_df
    
end