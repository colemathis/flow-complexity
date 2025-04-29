#==============================================================================#
# IMPORTS
#==============================================================================#


#==============================================================================#
# FUNCTIONS
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

        push!(next_taus, (chemostat_next_tau, next_reactor))
        sort!(next_taus, by=x -> x[1])

        keep_ones_fixed(Ensemble)

        if tau > checkpoint
            save_checkpoint(sim, tau)
            checkpoint += sim.params[:save_interval]
        end

    end

    return

end

#==============================================================================#

function calc_next_rxn_times(all_propensities, tau)

    # how many different propensities do we have?
    n = length(all_propensities)

    # sum these up
    total_p = sum.(all_propensities)
    tau_steps = -log.(rand(n))./total_p
    reactor_steps = collect(zip(tau_steps .+ tau, 1:n))
    reactor_steps = sort(reactor_steps, by = x->x[1])

    return reactor_steps

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
# END OF FILE
#==============================================================================#
