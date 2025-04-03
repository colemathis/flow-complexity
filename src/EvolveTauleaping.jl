using StatsBase
using Dates

#==============================================================================#
# FUNCTION
#==============================================================================#

function evolve_distributed_tau_leaping(sim; dry_run=false)

    Ensemble = sim.ensemble

    tau = 0.0
    checkpoint = 0.0

    dt = sim.params[:dt]
    
    t = 0

    total_constructive_rxn = 0
    total_destructive_rxn = 0
    total_diffusion_rxn = 0

    skipped_constructive_rxn = 0
    skipped_destructive_rxn = 0
    skipped_diffusion_rxn = 0

    if dry_run == true
        # record current time
        dry_start = Dates.now()
    end

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

            total_constructive_rxn += n_f
            total_destructive_rxn += n_b

            n_f_skipped = apply_tauleap_constructive_rxn(this_chemostat, tau, n_f)
            n_b_skipped = apply_tauleap_destructive_rxn(this_chemostat, tau, n_b)

            skipped_constructive_rxn += n_f_skipped
            skipped_destructive_rxn += n_b_skipped
                            
        end

        # # Apply diffusion
        for id in Ensemble.reactor_ids            

            this_chemostat = Ensemble.reactors[id]
            this_chemostat_propensities = all_propensities[id]
    
            A_d = this_chemostat_propensities[3]
            
            # Calculate the number of reactions in the time interval tau
            n_d = rand(Poisson(A_d * dt))

            total_diffusion_rxn += n_d

            n_d_skipped = apply_tauleap_diffusion_rxn(Ensemble, this_chemostat, tau, n_d)
            skipped_diffusion_rxn += n_d_skipped
            
        end                

        if sim.params[:mass] > 0
            keep_ones_fixed(Ensemble)
        end

        if sim.params[:inflow_mols] > 0
            for id in Ensemble.inflow_ids
                this_chemostat = Ensemble.reactors[id]
                nmol = sim.params[:inflow_mols] * dt
                for i in 1:nmol
                    insert_molecule_at_random(this_chemostat.molecules, 1)
                end
            end
        end

        if tau >= checkpoint
            save_checkpoint(sim, tau)
            checkpoint += sim.params[:output_time]
        end

        tau += dt

        if dry_run == true && tau >= sim.params[:total_time] * 0.1
            total_time = Dates.now() - dry_start
            total_time = Dates.value(total_time) / 1000
            forecasted_time = (sim.params[:total_time] / tau) * (total_time)
            println("")
            println("Dry Run Completed. Time taken: $(round(total_time, digits=2)) seconds. Forecasted time: $(round(forecasted_time, digits=2)) seconds.")
            exit()
        end
    
    end

    println("")

    pc = 100 * skipped_constructive_rxn / total_constructive_rxn
    if isnan(pc) pc = 0 end
    pc = round(pc, digits = 0)
    println("Performed $total_constructive_rxn constructive reactions, skipped $skipped_constructive_rxn ($pc %)")

    pc = 100 * skipped_destructive_rxn / total_destructive_rxn
    if isnan(pc) pc = 0 end
    pc = round(pc, digits = 0)
    println("Performed $total_destructive_rxn destructive reactions, skipped $skipped_destructive_rxn ($pc %)")

    pc = 100 * skipped_diffusion_rxn / total_diffusion_rxn
    if isnan(pc) pc = 0 end
    pc = round(pc, digits = 0)
    println("Performed $total_diffusion_rxn diffusion reactions, skipped $skipped_diffusion_rxn ($pc %)")

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
        # println("        warning: skipped $skipped out of $n_f constructive reactions at time $tau")
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
        # println("        warning: skipped $skipped out of $n_b destructive reactions at time $tau")
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
        # println("        warning: skipped $skipped out of $n_d diffusion reactions at time $tau")
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

