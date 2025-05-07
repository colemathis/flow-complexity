#==============================================================================#
# IMPORTS
#==============================================================================#

import CSV
import Arrow
import DataFrames
import Pkg

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function create_params_jl_file()

    if isfile("params.jl")
        println("params.jl already exists. Aborting to avoid overwrite.")
        exit(1)
    end

    params_text = """
    nrepeat = 1

    params_template = Dict(
        # simulation parameters
        :save_interval    	=> 1,
        :method         	=> "tau-leaping",
        :dt             	=> 1e-3,
        :random_seed    	=> "random",

        # physical parameters
        :total_time     	=> 100,
        :initial_mass       => 0,

        # topological parameters
        :graph_type         => "lattice-2way",
        :randomize_edges    => false,
        :N_reactors     	=> 25,

        # reaction rates
        :inflow_mols        => 1e3,
        :forward_rate   	=> 1e-3,
        :diffusion_rate   	=> 1e-3,
        :outflow_rate       => 1e-3,
    )
    """

    open("params.jl", "w") do io
        write(io, params_text)
    end

    println("Wrote params.jl in ", pwd())

end

#==============================================================================#

function create_params_csv_file()

    include(joinpath(pwd(), "params.jl")) # creates variables params_template, nrepeat

    # Generate all combinations of parameters
    param_keys = collect(keys(params_template))
    values_list = [ isa(params_template[k], AbstractVector) ? params_template[k] : [params_template[k]] for k in param_keys ]
    param_tuples = Iterators.product(values_list...)
    params_array = [ Dict(zip(param_keys, t)) for t in param_tuples ]

    # Duplicate each parameter set nrepeat times so they can be edited independently
    if nrepeat > 1
        params_array = [deepcopy(d) for d in params_array for _ in 1:nrepeat]
    end

    for i in eachindex(params_array)
        params_array[i][:sim_number] = i
    end

    for i in eachindex(params_array)
        if params_array[i][:outflow_rate] == "equal-to-diffusion-rate"
            params_array[i][:outflow_rate] = params_array[i][:diffusion_rate]
        end
    end

    array_fn = "./data/params.csv"
    println("Writing parameters for simulation array in $array_fn")
    mkpath("data")
    CSV.write(array_fn, params_array)
    println("Done.")
    println("")

end

#==============================================================================#
# END OF FILE
#==============================================================================#



