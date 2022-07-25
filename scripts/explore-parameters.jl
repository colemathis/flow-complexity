# July 13 2022 
# Exploring parameters of well mixed model, line chain, and lattice 

using DrWatson
using Distributed
@quickactivate

@everywhere include("../src/Simulation.jl")

function logunif(min, max)
    scale = log10(max) - log10(min)
    r = 10^(scale*rand() + log10(min))
    return r 
end


function run_lattice_reactions(N_reactors_list, f_rate_list, o_rate_list, current_sim)
    @sync @distributed for N in N_reactors_list
        i = 0
        for f in f_rate_list
            for o in o_rate_list
                i += 1
                sim_number  = current_sim + N*1000 + i
                this_sim = Simulation(1000, "lattice", N ,f, o, sim_number = sim_number, notes = "Lattice Parameter Sweep July 14 2022")
                RunSimulation(this_sim)
            end
        end
    end
end

function run_line_reactions(N_reactors_list, f_rate_list, o_rate_list, current_sim)
    @sync @distributed for N in N_reactors_list
        i = 0
        for f in f_rate_list
            for o in o_rate_list
                i += 1
                sim_number  = current_sim + N*1000 + i
                this_sim = Simulation(1000, "line", N ,f, o, sim_number = sim_number, notes = "Line Parameter Sweep July 14 2022")
                RunSimulation(this_sim)
            end
        end
    end
end

function run_mixed_reactions(f_rate_list, o_rate_list)
    for f in f_rate_list
        for o in o_rate_list
            this_sim = Simulation(1000, "line", 1,f, o, notes = "Mixed Parameter Sweep July 14 2022")
            RunSimulation(this_sim)
        end
    end
end

function run_all_topologies()
    f_rate_list = [0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
    o_rate_list = [5.0, 10.0, 50.0, 100.0]
    N_reactor_list = [4, 9, 16, 25]
    current_sim = 26040
    println("Running Line Reactions")
    run_line_reactions(N_reactor_list, f_rate_list, o_rate_list, current_sim)
    current_sim = 26040 + 97
    println("Running Lattice Reaction")
    run_lattice_reactions(N_reactor_list, f_rate_list, o_rate_list, current_sim)
    current_sim = 26040 + (97*2)
    println("Running Mixed Reactions")
    run_mixed_reactions(f_rate_list, o_rate_list)
end
