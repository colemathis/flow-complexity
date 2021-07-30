using DataFrames
using FileIO
using Glob
using DelimitedFiles

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
    data_dict = load(bfile)
    recorded_vars = [k for k in keys(data_dict)]
    times = [t for t in keys(data_dict[recorded_vars[1]])]

    tidy_df = DataFrame(time = times)
    for var in recorded_vars
        tidy_df[!, var] = map(x-> data_dict[var][x], tidy_df[!,:time])
    end

    return tidy_df
end

