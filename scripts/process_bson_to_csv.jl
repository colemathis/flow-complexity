using DataFrames
using FileIO
using Glob
using DelimitedFiles

using Pathways
## July 29th processing initial bson files to csv for plotting

function initial_bsons_to_csv(outputname::String, directory::String)
    all_bsons = glob(directory*"/*.bson")
    complete_df = DataFrame()

    for bfile in all_bsons
        parameters = parse_fname(bfile)
        exp_df = bson_to_tidy_df(bfile)

        for (param, value) in parameters
           exp_df[!,param] = repeat([value], nrow(exp_df)) 
        end
        complete_df = vcat(complete_df, exp_df)
    end
    writedlm(outputname, Iterators.flatten(([names(complete_df)], eachrow(complete_df))), ',')
end

function parse_fname(bson_fname)
    basename = split(bson_fname, "\\")[end]
    string_params = split(basename, "_")[1:end-1]
    repeat = Int64(parse(Float64, string_params[1]))
    iterations = Int64(parse(Float64, string_params[2]))
    outflow = parse(Float64, string_params[3])
    epsilion = parse(Float64, string_params[4])

    parameters = Dict(:repeat => repeat, :iterations => iterations, :outflow => outflow, :epsilion => epsilion)
    return parameters
end

function bson_to_tidy_df(bfile)

    csvfile = split(bfile, ".bson")[1] * ".csv"
    data_dict = load(bfile)
    recorded_vars = [k for k in keys(data_dict)]
    times = [t for t in keys(data_dict[recorded_vars[1]])]
    data = []
    for t in times
        for var in recorded_vars
            if var == :complete_timeseries
                time_counts = countmap(data_dict[:complete_timeseries][t])
                for (v,c) in time_counts
                    push!(data, Dict("time"=> t, "variable"=>string(v), "value"=> c))
                end
            else
                push!(data, Dict("time" => t, "variable" => String(var), "value"=> data_dict[var][t] ))
            end
        end
    end

    tidy_df = DataFrame(time = map(x -> x["time"], data),
                        variable = map(x -> x["variable"], data),
                        value = map(x -> x["value"], data))

    writedlm(csvfile, Iterators.flatten(([names(tidy_df)], eachrow(tidy_df))), ',')
end

