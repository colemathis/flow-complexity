using FlowComplexity
const FC = FlowComplexity

# relative_path = FC.get_relative_path()

# current_directory = pwd()
# fn = joinpath(current_directory, "params.jl")
# include(fn)

include("params.jl")
mkpath("data")
FC.write_params_file(params_array)
