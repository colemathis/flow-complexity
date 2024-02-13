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
    # compute parameter permutation
    params = [(N,f,o) for N in N_reactors_list, f in f_rate_list, o in o_rate_list]
    # get the number of permutations/simulations
    nparams = length(params)
    # define vectors that will hold the parameters & index
    Ns = zeros(Int, nparams)
    fs = zeros(Float64, nparams)
    os = zeros(Float64, nparams)
    inds = zeros(Int, nparams)
    # pre-calculate the parameters for each simulation
    for i in 1:nparams
        Ns[i] = params[i][1]
        fs[i] = params[i][2]
        os[i] = params[i][3]
        inds[i] = i
    end
    # loop over all parameter combination
    @sync @distributed for i in 1:nparams
        N = Ns[i]
        f = fs[i]
        o = os[i]
        ind = inds[i]
        println("N=$N f=$f o=$o i=$i")
        sim_number  = current_sim + N*1000 + i
        this_sim = Simulation(1000, "lattice", N ,f, o, sim_number = sim_number, notes = "Lattice Parameter Sweep July 14 2022")
        RunSimulation(this_sim)
    end
end

function run_line_reactions(N_reactors_list, f_rate_list, o_rate_list, current_sim)
    # compute parameter permutation
    params = [(N,f,o) for N in N_reactors_list, f in f_rate_list, o in o_rate_list]
    # get the number of permutations/simulations
    nparams = length(params)
    # define vectors that will hold the parameters & index
    Ns = zeros(Int, nparams)
    fs = zeros(Float64, nparams)
    os = zeros(Float64, nparams)
    inds = zeros(Int, nparams)
    # pre-calculate the parameters for each simulation
    for i in 1:nparams
        Ns[i] = params[i][1]
        fs[i] = params[i][2]
        os[i] = params[i][3]
        inds[i] = i
    end
    # loop over all parameter combination
    @sync @distributed for i in 1:nparams
        N = Ns[i]
        f = fs[i]
        o = os[i]
        ind = inds[i]
        println("N=$N f=$f o=$o i=$i")
        sim_number  = current_sim + N*1000 + ind
        this_sim = Simulation(1000, "line", N ,f, o, sim_number = sim_number, notes = "Line Parameter Sweep July 14 2022")
        RunSimulation(this_sim)
    end
end

function run_mixed_reactions(f_rate_list, o_rate_list)
    # compute parameter permutation
    params = [(f,o) for f in f_rate_list, o in o_rate_list]
    # get the number of permutations/simulations
    nparams = length(params)
    # define vectors that will hold the parameters & index
    fs = zeros(Float64, nparams)
    os = zeros(Float64, nparams)
    # pre-calculate the parameters for each simulation
    for i in 1:nparams
        fs[i] = params[i][1]
        os[i] = params[i][2]
    end
    # loop over all parameter combination
    @sync @distributed for i in 1:nparams
        f = fs[i]
        o = os[i]
        println("f=$f o=$o")
        this_sim = Simulation(1000, "line", 1,f, o, notes = "Mixed Parameter Sweep July 14 2022")
        RunSimulation(this_sim)
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

run_all_topologies()
