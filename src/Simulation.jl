include("Pipeline.jl")
include("Chemostats.jl")
include("Ensemble.jl")
include("EvolveStochastic.jl")
include("EvolveTauleaping.jl")
include("Analysis.jl")

using CSV, DataFrames

# Set the backward rate to a constant
const backward_rate = 1.0

#==============================================================================#
# DATA TYPES
#==============================================================================#

mutable struct Simulation

    params   ::Dict{Symbol, Any}
    output   ::Dict{Symbol, Any}
    ensemble ::Ensemble

end

#==============================================================================#

function Simulation(; user_params = Dict{Symbol, Any}())

    params = apply_params_default_values()
    params = apply_params_user_values(params, user_params)
    
    output = Dict{Symbol, Any}()
    ensemble = generate_ensemble_from_specs(params)
    
    return Simulation(params, output, ensemble)
    
end

#==============================================================================#
# FUNCTIONS
#==============================================================================#

function apply_params_default_values()

    default_values = Dict(
        :mass           => 1000,
        :graph_type     => "line",
        :N_reactors     => 1,
        :forward_rate   => 1e-3,
        :outflow_rate   => 0.0,
        :total_time     => 100,
        :output_time    => 1.0,
        :N_inflow       => 1,
        :sim_number     => 1,
        :method         => "exact",
        :method         => "tau-leaping",
        :save_name      => get_relative_path(),
        :random_seed    => 1,
    )
    
    return default_values

end

#==============================================================================#

function apply_params_user_values(previous_params, user_params)

    final_params = copy(previous_params)
    
    # overwrite with the ones provided by the user
    for (user_key, user_value) in user_params
        final_params[user_key] = user_value
    end
    
    return final_params

end

#==============================================================================#

function generate_ensemble_from_specs(params)

    # generate a list of chemostat specifications which can be passed to the Ensemble constructor
    chemostat_list = gen_chemostat_spec_list(params[:N_reactors],
                                             params[:forward_rate],
                                             backward_rate,
                                             params[:outflow_rate])

    # create an ensemble of chemostats with the list of specifications
    ensemble = Ensemble(params[:N_reactors],
                             params[:graph_type],
                             params[:N_inflow],
                             params[:mass],
                             chemostat_list=chemostat_list)

    return ensemble

end

#==============================================================================#

function RunSimulation(sim)

    sim_number = sim.params[:sim_number]

    sim_number_string = lpad(sim_number,6,"0")
    sim.params[:save_name] = projectdir("milestones", sim.params[:save_name], "data", "sims", sim_number_string)
    
    sim.output[:timestamps] = []
    sim.output[:populations] = []

    elapsed_time = @elapsed begin
        if sim.params[:method] == "exact"
            println("Launching simulation using exact algorithm...")
            evolve_distributed_exact(sim)
        elseif sim.params[:method] == "tau-leaping"
            println("Launching simulation using tau-leaping algorithm...")
            evolve_distributed_tau_leaping(sim)
        else
            println("error: method unknown")
            exit()
        end
    end

    elapsed_time = round(elapsed_time, digits = 2)
    println("Sim Completed. Time taken: $elapsed_time seconds.")

    save_data(sim)
    println("")
    
end

#==============================================================================#

function sim_number_string(sim_number)

    sim_number_string = lpad(sim_number, 6, "0")

    return sim_number_string

end

#==============================================================================#

function save_data(sim)

    sim_number = string(sim.params[:sim_number])    
    fn = joinpath(sim.params[:save_name], "simulation.jld2")
    save(fn, Dict("sim" => sim))

    rel = relpath(fn, pwd())
    println("Data saved in $rel")

end

#==============================================================================#

function save_checkpoint(sim, tau::Float64)
    
    i = round(tau, digits=3)
    # println("   saving at t=$i")
    print(".")
    
    # record the current time
    push!(sim.output[:timestamps], i)
    
    # loop over reactors then and add vector of molecules
    push!(sim.output[:populations], [])
    for i in 1:sim.params[:N_reactors]
        copy_of_molecules = copy(sim.ensemble.reactors[i].molecules)
        push!(sim.output[:populations][end], copy_of_molecules)
    end

    # the resulting ndim array will have dimensions:
    # 1 = time id (NOT time stamp)
    # 2 = reactor id
    # 3 = molecule vector
    #
    # so if time id is "t", reactor id is "i" and we want to access molecule "m",
    # sim.populations[t, i, m]
    
end

#==============================================================================#

