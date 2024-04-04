"""

##########################################################################

EXPLORE-TEST



##########################################################################

"""

using DrWatson
using Distributed
@quickactivate

@everywhere include("../src/Simulation.jl")

"""

##########################################################################
LOG UNIFORM SAMPLE
#p3 (seems unused)

Generates a random number from a log-uniform distribution between the 
specified minimum and maximum values.
##########################################################################

Input:

    min (float) = minimum
    max (float) = maximum

Output:

    r   (float) = random number

##########################################################################

"""

function logunif(min, max)

    # get the distribution min/max
    scale = log10(max) - log10(min)

    # pick a random number
    r = 10^(scale*rand() + log10(min))

    return r 
end

"""

##########################################################################
RUN LATTICE REACTIONS

Takes lists of parameters, computes all permutations and run these sims
in parallel.
##########################################################################

Input:

    N_reactors_list (array) = number of reactors for each simulation
    f_rate_list     (array) = forward reaction rates
    o_rate_list     (array) = outflow reaction rates
    current_sim     (int)   = next simulation number to be used

Output:

    (none)

##########################################################################

"""

function run_lattice_reactions(N_reactors_list, 
                               f_rate_list, 
                               o_rate_list, 
                               current_sim
                               )

    # compute parameter permutation
    params = [(N,f,o) for N in N_reactors_list, f in f_rate_list, o in o_rate_list]

    # get the number of permutations/simulations
    nparams = length(params)

    # define vectors that will hold the parameters & index
    Ns = zeros(Int, nparams)
    fs = zeros(Float64, nparams)
    os = zeros(Float64, nparams)
    inds = zeros(Int, nparams)      #p3: unused

    # pre-calculate the parameters for each simulation
    for i in 1:nparams
        Ns[i] = params[i][1]
        fs[i] = params[i][2]
        os[i] = params[i][3]
        inds[i] = i                 #p3: unused
    end

    # loop over all parameter combination and run the simulations in parallel
    @sync @distributed for i in 1:nparams

        # get the parameters for this simulation
        N = Ns[i]
        f = fs[i]
        o = os[i]
        ind = inds[i]               #p3: unused

        println("N=$N f=$f o=$o i=$i")

        # calculate the sim number based on the # of reactors and parameter permutation index
        sim_number  = current_sim + N*1000 + i

        # run the simulation
        this_sim = Simulation(1000, "lattice", N ,f, o, sim_number = sim_number, notes = "Lattice Parameter Sweep July 14 2022")
        RunSimulation(this_sim)

    end
end

"""

##########################################################################
RUN LINE REACTIONS

Takes lists of parameters, computes all permutations and run these sims
in parallel.
##########################################################################

Input:

    N_reactors_list (array) = number of reactors for each simulation
    f_rate_list     (array) = forward reaction rates
    o_rate_list     (array) = outflow reaction rates
    current_sim     (int)   = next simulation number to be used

Output:

    (none)

##########################################################################

"""

function run_line_reactions(N_reactors_list, 
                            f_rate_list, 
                            o_rate_list, 
                            current_sim
                            )

    # compute parameter permutation
    params = [(N,f,o) for N in N_reactors_list, f in f_rate_list, o in o_rate_list]

    # get the number of permutations/simulations
    nparams = length(params)

    # define vectors that will hold the parameters & index
    Ns = zeros(Int, nparams)
    fs = zeros(Float64, nparams)
    os = zeros(Float64, nparams)
    inds = zeros(Int, nparams)      #p3: unused

    # pre-calculate the parameters for each simulation
    for i in 1:nparams
        Ns[i] = params[i][1]
        fs[i] = params[i][2]
        os[i] = params[i][3]
        inds[i] = i                 #p3: unused
    end

    # loop over all parameter combination and run the simulations in parallel
    @sync @distributed for i in 1:nparams

        # get the parameters for this simulation
        N = Ns[i]
        f = fs[i]
        o = os[i]
        ind = inds[i]               #p3: unused

        println("N=$N f=$f o=$o i=$i")
        
        # calculate the sim number based on the # of reactors and parameter permutation index
        sim_number  = current_sim + N*1000 + ind    #p3 replace with i

        # run the simulation
        this_sim = Simulation(1000, "line", N ,f, o, sim_number = sim_number, notes = "Line Parameter Sweep July 14 2022")
        RunSimulation(this_sim)

    end
end

"""

##########################################################################
RUN MIXED REACTIONS

Takes lists of parameters, computes all permutations and run these sims
in parallel.
##########################################################################

Input:

    f_rate_list     (array) = forward reaction rates
    o_rate_list     (array) = outflow reaction rates

Output:

    (none)

##########################################################################

"""

function run_mixed_reactions(f_rate_list, 
                             o_rate_list
                             )

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
    
    # loop over all parameter combination and run the simulations in parallel
    @sync @distributed for i in 1:nparams

        # get the parameters for this simulation
        f = fs[i]
        o = os[i]

        println("f=$f o=$o")

        # run the simulation
        this_sim = Simulation(1000, "line", 1,f, o, notes = "Mixed Parameter Sweep July 14 2022")
        RunSimulation(this_sim)

    end
end

"""

##########################################################################
RUN ALL TOPOLOGIES

Defines lists of parameters and runs line, lattice and mixed reactions
exploring all permutations of these parameters.
##########################################################################

Input:

    (none)

Output:

    (none)

##########################################################################

"""

function run_all_topologies()

    # define the parameters we’ll explore
    f_rate_list = [0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
    o_rate_list = [5.0, 10.0, 50.0, 100.0]
    N_reactor_list = [4, 9, 16, 25]

    # define next sim number we’ll be using and run the sim
    current_sim = 26040
    println("Running Line Reactions")
    run_line_reactions(N_reactor_list, f_rate_list, o_rate_list, current_sim)

    # define next sim number we’ll be using and run the sim
    current_sim = 26040 + 97
    println("Running Lattice Reaction")
    run_lattice_reactions(N_reactor_list, f_rate_list, o_rate_list, current_sim)

    # define next sim number we’ll be using and run the sim
    current_sim = 26040 + (97*2)
    println("Running Mixed Reactions")
    run_mixed_reactions(f_rate_list, o_rate_list)

end

run_all_topologies()
