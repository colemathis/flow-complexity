# July 13 2022 
# Exploring parameters of well mixed model, line chain, and lattice 

using DrWatson

@quickactivate

include("../src/Simulation.jl")

function logunif(min, max)
    scale = log10(max) - log10(min)
    r = 10^(scale*rand() + log10(min))
    return r 
end


function run_lattice_reactions(N_reactors_list, f_rate_list, o_rate_list)
    for N in N_reactors_list
        for f in f_rate_list
            for o in o_rate_list
                this_sim = Simulation(1000, "lattice", N ,f, o, sim_notes = "Lattice Parameter Sweep July 13 2022")
                RunSimulation(this_sim)
            end
        end
    end
end

function run_line_reactions(N_reactors_list, f_rate_list, o_rate_list)
    for N in N_reactors_list
        for f in f_rate_list
            for o in o_rate_list
                this_sim = Simulation(1000, "line", N ,f, o, sim_notes = "Line Parameter Sweep July 13 2022")
                RunSimulation(this_sim)
            end
        end
    end
end

function run_mixed_reactions(f_rate_list, o_rate_list)
    for f in f_rate_list
        for o in o_rate_list
            this_sim = Simulation(1000, "line", 1,f, o, sim_notes = "Line Parameter Sweep July 13 2022")
            RunSimulation(this_sim)
        end
    end
end

function run_all_topologies()
    f_rate_list = [0.0005, 0.001, 0.005, 0.01, 0.05]
    o_rate_list = [1.0, 2.0, 5.0, 10.0]
    N_reactor_list = [4, 9, 16]
    println("Running Mixed Reactions")
    run_mixed_reactions(f_rate_list, o_rate_list)
    println("Running Line Reactions")
    run_line_reactions(N_reactor_list, f_rate_list, o_rate_list)
    println("Running Lattice Reaction")
    run_lattice_reactions(N_reactor_list, f_rate_list, o_rate_list)
end
