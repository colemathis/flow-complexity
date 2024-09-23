# Usage: launch.jl <CSV file> <sim>

using DrWatson
@quickactivate :FlowComplexity

using DataFrames

function main(args)
    array_fn, sim = args[1], parse(Int, args[2])
    
    println("Loading parameter file $array_fn.")
    params_df = FlowComplexity.GetParamsFromCSV(array_fn)
    
    println("Reading parameters for simulation number $sim.")
    params_dict = FlowComplexity.DF_row_to_dict(params_df, sim)
    
    println("Creating simulation object.")
    simulation = FlowComplexity.Simulation(; params_dict...)
    
    println("Launching simulation.")
    FlowComplexity.RunSimulation(simulation)
end

main(ARGS)