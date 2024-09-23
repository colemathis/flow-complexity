"""
##########################################################################

launch.jl

Usage:
    launch.jl <CSV file> <sim>

where:
    <CSV file>  : path of the CSV file relative to data/
    <sim>       : simulation number

##########################################################################
"""

"""
##########################################################################
IMPORTS
##########################################################################
"""

using DrWatson
@quickactivate :FlowComplexity

using DataFrames

"""
##########################################################################
FUNCTIONS
##########################################################################
"""

function print_usage()
    println("""
    Usage:
        launch.jl <CSV file> <sim>

    where:
        <CSV file>  : path of the CSV file relative to data/
        <sim>       : simulation number
    """)
end

"""
##########################################################################
MAIN
##########################################################################
"""

# Check the number of command line arguments
if length(ARGS) < 2
    print_usage()
    exit(1)
end

# load the parameters from the CSV, e.g. "sims/MS00_MWE/array1.csv"
array_fn = ARGS[1]
println("Loading parameter file $array_fn.")
params_df = FlowComplexity.GetParamsFromCSV(array_fn)

# get a param dictionary for simulation number <sim>
sim = parse(Int, ARGS[2])
println("Reading parameters for simulation number $sim.")
params_dict = FlowComplexity.DF_row_to_dict(params_df, sim)

# create simulation object
println("Creating simulation object.")
s = FlowComplexity.Simulation(; params_dict...)

# launch it
println("Launching simulation.")
FlowComplexity.RunSimulation(s)
