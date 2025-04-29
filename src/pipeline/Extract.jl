#==============================================================================#
# IMPORTS
#==============================================================================#

import Distributions

#==============================================================================#
# FUNCTION
#==============================================================================#

function extract_sims()

    # new function consolidating CSV files into a single CSV file
    # and Arrow files for faster loading
    
    sims_dir = "./data/sims"
    sim_folders = filter(f -> isdir(joinpath(sims_dir, f)), readdir(sims_dir))
    
    all_timeseries = DataFrames.DataFrame()
    all_graphs = DataFrames.DataFrame()
    
    for folder in sim_folders
        folder_path = joinpath(sims_dir, folder)
        
        timeseries_path = joinpath(folder_path, "timeseries.csv")
        graph_path = joinpath(folder_path, "graph.csv")
        
        if isfile(timeseries_path)
            ts = DataFrames.DataFrame(CSV.File(timeseries_path))
            ts[!, :sim_number] = fill(parse(Int, folder), DataFrames.nrow(ts))
            append!(all_timeseries, ts)
        end
        
        if isfile(graph_path)
            g = DataFrames.DataFrame(CSV.File(graph_path))
            g[!, :sim_number] = fill(parse(Int, folder), DataFrames.nrow(g))
            append!(all_graphs, g)
        end
    end
    
    CSV.write("data/timeseries.csv", all_timeseries)
    CSV.write("data/graphs.csv", all_graphs)

    Arrow.write("data/timeseries.arrow", all_timeseries; compress=:lz4)
    Arrow.write("data/graphs.arrow", all_graphs; compress=:lz4)
    
    println("Merged timeseries and graphs saved to ./data/")

end

#==============================================================================#
# END OF FILE
#==============================================================================#