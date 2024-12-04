using StatsBase

#==============================================================================#
# FUNCTION
#==============================================================================#

function evolve_distributed_exact(sim)

    Ensemble = sim.ensemble

    tau = 0.0
    checkpoint = 0.0

    all_propensities = calc_propensities.(Ensemble.reactors)
    next_taus = calc_next_rxn_times(all_propensities, tau)

    while tau < sim.total_time
        
        next_rxn = popfirst!(next_taus)
        next_reactor = next_rxn[2]
        tau = next_rxn[1]

        next_reactor_propensities = all_propensities[next_reactor]
        rxn = sample(["construction", "degradation", "outflow"], Weights(next_reactor_propensities))
        this_chemostat = Ensemble.reactors[next_reactor]

        if rxn == "construction"
            this_chemostat = constructive_rxn(this_chemostat)
        elseif rxn == "degradation"
            this_chemostat = destructive_rxn(this_chemostat)
        elseif rxn == "outflow"
            this_chemostat, outflow_direction = outflow_rxn(this_chemostat)
            if collect(keys(outflow_direction)) != []
                outflow_target = collect(keys(outflow_direction))[1]
                outflow_molecule = outflow_direction[outflow_target]
                push!(Ensemble.reactors[outflow_target].molecules, outflow_molecule)
                #println(Ensemble.reactors[outflow_target].molecules)
                # Update propensitie for outflow target based on diffusion
                all_propensities[outflow_target] = calc_propensities(Ensemble.reactors[outflow_target])
                next_taus = [nt for nt in next_taus if nt[2] != outflow_target]
                outflow_total_p = sum(all_propensities[outflow_target])
                outflow_target_next_tau = tau - log(rand()) / outflow_total_p
                push!(next_taus, (outflow_target_next_tau, outflow_target))
            end
        end

        # Update Chemostat and calculate next reaction time
        Ensemble.reactors[next_reactor] = this_chemostat # i don’t know what this line does
        all_propensities[next_reactor] = calc_propensities(Ensemble.reactors[next_reactor])
        this_chemostat_total_p = sum(all_propensities[next_reactor])
        chemostat_next_tau = tau - log(rand()) / this_chemostat_total_p
        # println(-log(rand()) / this_chemostat_total_p)
        push!(next_taus, (chemostat_next_tau, next_reactor))
        sort!(next_taus, by=x -> x[1])

        keep_ones_fixed(Ensemble)

        if tau > checkpoint
            save_checkpoint(sim, tau)
            checkpoint += sim.output_time
        end

    end

    return

end

#==============================================================================#

function constructive_rxn(chemostat::Chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)
    shuffle!(molecules)

    # combine two molecules
    a = pop!(molecules)
    b = pop!(molecules)
    c = a + b

    # add the new one to the end of the list
    push!(molecules, c)

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat

end

#==============================================================================#

function destructive_rxn(chemostat::Chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)
    shuffle!(molecules) # Shuffle the molecules

    # ignore molecules of length 1
    big_moles = filter(x -> x > 1, molecules)

    if length(big_moles) > 0

        # pick one molecule and remove it from the list
        a = pop!(big_moles)
        mole_index = findfirst(x -> x == a, molecules)
        a = popat!(molecules, mole_index)

        # split it at a random point
        b = sample(1:(a-1))
        c = a - b

        # add back both to the list
        push!(molecules, b)
        push!(molecules, c)

    end

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat

end

#==============================================================================#

function outflow_rxn(chemostat)

    # create temporary molecules array and shuffle it
    molecules = copy(chemostat.molecules)

    # get the list of neighbors and their inflow to this chemostat
    neighbors = chemostat.neighbors
    neighbor_weights = chemostat.neighbor_flows

    # shuffle the molecules and pick one
    shuffle!(molecules)
    a = pop!(molecules)

    # pick a neighbor and ??? #p1: not sure what this is about
    if neighbors != []
        neighbor = sample(neighbors, Weights(neighbor_weights))
        outflow_direction = Dict(neighbor => a)
    else
        outflow_direction = Dict{Int64,String}()
    end

    # copy back the list of molecules
    chemostat.molecules = copy(molecules)

    return chemostat, outflow_direction

end

#==============================================================================#

function evolve_distributed_tau_leaping(sim)

    Ensemble = sim.ensemble

    tau = 0.0
    checkpoint = 0.0

    dt = 0.01
    t = 0

    while tau <= sim.params[:total_time]
        tau = round(tau, digits=5) # workaround for numerical errors
    
        all_propensities = calc_propensities.(Ensemble.reactors)

        # Apply constructive/destructive reactions
        for id in Ensemble.reactor_ids

            this_chemostat = Ensemble.reactors[id]
            this_chemostat_propensities = all_propensities[id]

            A_f = this_chemostat_propensities[1]
            A_b = this_chemostat_propensities[2]

            # Calculate the number of reactions in the time interval tau
            n_f = rand(Poisson(A_f * dt))
            n_b = rand(Poisson(A_b * dt))

            apply_tauleap_constructive_rxn(this_chemostat, tau, n_f)
            apply_tauleap_destructive_rxn(this_chemostat, tau, n_b)
                            
        end

        # # Apply diffusion
        for id in Ensemble.reactor_ids            

            this_chemostat = Ensemble.reactors[id]
            this_chemostat_propensities = all_propensities[id]
    
            A_d = this_chemostat_propensities[3]
            
            # Calculate the number of reactions in the time interval tau
            n_d = rand(Poisson(A_d * dt))

            apply_tauleap_diffusion_rxn(Ensemble, this_chemostat, tau, n_d)
        end                
        
        keep_ones_fixed(Ensemble)

        if tau >= checkpoint
            save_checkpoint(sim, tau)
            checkpoint += sim.params[:output_time]
        end

        tau += dt
    
    end


    return

end

#==============================================================================#

function apply_tauleap_constructive_rxn(this_chemostat, tau, n_f)

    skipped = 0

    n_molecules = length(this_chemostat.molecules) # number of molecules
    nmax_reactions = max(n_molecules - 1, 0)       # number of constructive reactions that can happen (excluding negative numbers)
    if n_f <= nmax_reactions
        # perform number of reactions intended
        n_reactions = n_f
    else
        # perform reduced count of reactions and warn
        n_reactions = nmax_reactions
        skipped = n_f - n_reactions
        println("        warning: skipped $skipped out of $n_f constructive reactions at time $tau")
    end

    for i in 1:n_reactions
        a = pick_molecule_at_random(this_chemostat.molecules)
        b = pick_molecule_at_random(this_chemostat.molecules)        
        c = a + b        
        insert_molecule_at_random(this_chemostat.molecules, c)        
    end

    return skipped

end

#==============================================================================#

function apply_tauleap_destructive_rxn(this_chemostat, tau, n_b)

    skipped = 0

    molecules_gtr_one = filter(x -> x > 1, this_chemostat.molecules)

    nmax_reactions = sum(molecules_gtr_one .-1) # number of destructive reactions that can happen

    if n_b <= nmax_reactions
        # perform number of reactions intended
        n_reactions = n_b
    else
        # perform reduced count of reactions and warn
        n_reactions = nmax_reactions
        skipped = n_b - n_reactions
        println("        warning: skipped $skipped out of $n_b destructive reactions at time $tau")
    end
    
    for i in 1:n_reactions

        # pick a molecule >1
        a = pick_molecule_at_random(molecules_gtr_one)

        # find it in the list of molecules
        mol_idx = findfirst(x -> x == a, this_chemostat.molecules)
        a = popat!(this_chemostat.molecules, mol_idx)
        
        # split it
        b = sample(1:(a-1)) # sample.(UnitRange.(1, (new_mols[(end-n_b):end] .-1)))
        c = a - b

        # add them back to the list
        if b > 1
            insert_molecule_at_random(molecules_gtr_one, b)
        end
        if c > 1
            insert_molecule_at_random(molecules_gtr_one, c)
        end
        insert_molecule_at_random(this_chemostat.molecules, b)
        insert_molecule_at_random(this_chemostat.molecules, c)

    end

    return skipped

end

#==============================================================================#

function apply_tauleap_diffusion_rxn(Ensemble, this_chemostat, tau, n_d)
    
    skipped = 0

    n_molecules = length(this_chemostat.molecules) # number of molecules
    nmax_reactions = n_molecules                   # maximum number of diffusion reactions
    if n_d <= nmax_reactions
        # perform intended number of reactions
        n_reactions = n_d
    else
        # perform reduced number of reactions and warn
        n_reactions = nmax_reactions
        skipped = n_d - n_reactions
        println("        warning: skipped $skipped out of $n_d diffusion reactions at time $tau")
    end

    for i in 1:n_reactions
        a = pick_molecule_at_random(this_chemostat.molecules)

        # get the list of neighbors and their inflow to this chemostat
        neighbors = this_chemostat.neighbors
        neighbor_weights = this_chemostat.neighbor_flows
    
        if neighbors != []
            neighbor = sample(neighbors, Weights(neighbor_weights))
            insert_molecule_at_random(Ensemble.reactors[neighbor].molecules, a)
        end

    end

    return skipped
    
end

#==============================================================================#

function pick_molecule_at_random(molecules)

    n = length(molecules)
    random_index = rand(1:n)
    a = popat!(molecules, random_index)

    return a
    
end

#==============================================================================#

function insert_molecule_at_random(molecules, m)

    n = length(molecules)
    random_index = rand(1:n+1)
    insert!(molecules, random_index, m)

end

#==============================================================================#

function save_checkpoint(sim, tau::Float64)
    
    i = round(tau, digits=3)
    println("   saving at t=$i")

    # record the current time
    push!(sim.output[:timestamps], i)
    
    # loop over reactors then and add vector of molecules
    push!(sim.output[:populations], [])
    for i in sim.params[:N_reactors]
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

function keep_ones_fixed(ensemble)
    for chemostat in ensemble.reactors
        if chemostat.fixed_mass != 0
            # Update the mass of the chemostat
            chemostat.mass = calc_mass(chemostat)
            # Calculate the current number of ones
            current_ones = calc_ones(chemostat)
            # Determine the difference from the fixed mass
            delta_ones = chemostat.fixed_mass - current_ones

            if delta_ones > 0
                # Add the necessary number of ones
                append!(chemostat.molecules, ones(Int, delta_ones))
            end
        end
    end
end

