using StatsBase
using Statistics
include("Chemostats.jl")
include("Ensemble.jl")

function evolve_well_mixed(chemostat::Chemostat, tau_max::Float64, output_freq::Float64, outputs::Array{Symbol,1})
    ## Time evolution of a single well mixed chemostat

    # Figure out what we need to record
    evolution_outputs = Dict{Symbol,Any}()
    if :complete_timeseries in outputs
        complete_ts = Dict{Int64, Array{Int64,1}}()
    end

    if :molecule_count in outputs
        mole_count = Dict{Int64, Int64}()
    end

    if :average_length in outputs
        ave_lengths = Dict{Int64, Float64}()
        var_lengths = Dict{Int64, Float64}()
    end


    propensities = calc_propensities(chemostat)
    tau = 0.0
    checkpoint = 0.0
    while tau < tau_max
        # Calculate Reaction Propensities
        
        # Pick reaction 
        rxn = sample(["construction", "degradation", "outflow"],Weights(propensities))

        # Execute Reaction 
        if rxn == "construction"
            chemostat = constructive_rxn(chemostat)
        elseif rxn == "degradation"
            chemostat = destructive_rxn(chemostat)
        elseif rxn == "outflow"
            chemostat, outflow_direction = outflow_rxn(chemostat)
        end
        
        # Check Mass 
        chemostat.mass = calc_mass(chemostat)
        if  chemostat.mass != chemostat.mass_fixed && chemostat.mass_fixed != 0
            delta_mass = chemostat.mass_fixed - chemostat.mass
            new_moles = repeat([1], delta_mass)
            chemostat.molecules = vcat(chemostat.molecules, new_moles)
        end
        
        # Update propensities
        propensities = calc_propensities(chemostat)

        # Update time 
        total_p = sum(propensities)
        tau_step = -log(rand())/total_p
        tau += tau_step

        # Record outputs
        if tau > checkpoint
            i = round(tau, digits =3)
            checkpoint += output_freq
            if :complete_timeseries in outputs
                this_ts = Dict(i => chemostat.molecules)
                complete_ts = merge(complete_ts, this_ts)
            end
        
            if :molecule_count in outputs
                this_count = Dict(i => length(chemostat.molecules))
                mole_count = merge(mole_count, this_count)
            end
        
            if :average_length in outputs
                these_lengths = chemostat.molecules
                this_ave = Dict(i => mean(these_lengths))
                ave_lengths = merge(ave_lengths, this_ave)

                this_var = Dict(i => var(these_lengths))
                var_lengths = merge(var_lengths, this_var)
            end
            
        end

        if :complete_timeseries in outputs
            evolution_outputs[:complete_timeseries] = complete_ts
        end

        if :molecule_count in outputs
            evolution_outputs[:molecule_count] = mole_count 
        end

        if :average_length in outputs
            evolution_outputs[:average_lengths] = ave_lengths 
            evolution_outputs[:variance_lengths] = var_lengths
        end
    end
    return evolution_outputs
end



function evolve_distributed(Ensemble::Ensemble, tau_max::Float64, output_freq::Float64, outputs::Array{Symbol,1})
    ## Time evolution of a distributed system 

    # # Figure out what we need to record
    evolution_outputs = Dict{Int64,Any}()
    for id in Ensemble.reactor_ids
        this_reactor_data = Dict{Symbol,Any}()
        for out in outputs
            this_reactor_data[out] = Dict{Any,Any}()
        end
        evolution_outputs[id] = this_reactor_data
    end
    
    tau = 0.0
    checkpoint = 0.0

    all_propensities = calc_propensities.(Ensemble.reactors)
    println(all_propensities)
    next_taus = calc_next_rxn_times(all_propensities, tau)

    while tau < tau_max
        # Determine next reaction 
        next_rxn = popfirst!(next_taus)
        next_reactor = next_rxn[2]
        tau = next_rxn[1]
        #println(tau)
        
        next_reactor_propensities = all_propensities[next_reactor] 
        # Pick reaction 
        rxn = sample(["construction", "degradation", "outflow"], Weights(next_reactor_propensities))
        #println(next_reactor, " \t", rxn)
        this_chemostat = Ensemble.reactors[next_reactor]
        # Execute Reaction and update propensities
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
                outflow_target_next_tau = tau - log(rand())/outflow_total_p
                push!(next_taus, (outflow_target_next_tau, outflow_target))
            end
        end

        # Update Chemostat and calculate next reaction time
        Ensemble.reactors[next_reactor] = this_chemostat
        all_propensities[next_reactor] = calc_propensities(Ensemble.reactors[next_reactor])
        this_chemostat_total_p = sum(all_propensities[next_reactor])
        chemostat_next_tau = tau - log(rand())/this_chemostat_total_p
        push!(next_taus, (chemostat_next_tau, next_reactor))
        sort!(next_taus, by= x->x[1]) 
        # Check Mass
        for chemostat in Ensemble.reactors
            if chemostat.mass_fixed != 0
                chemostat.mass = calc_mass(chemostat)
                if  chemostat.mass != chemostat.mass_fixed
                    delta_mass = chemostat.mass_fixed - chemostat.mass
                    new_moles = repeat([1], delta_mass)
                    chemostat.molecules = vcat(chemostat.molecules, new_moles)
                end
            end
        end

        # Record outputs
        if tau > checkpoint
            i = round(tau, digits =3)
            checkpoint += output_freq
            println(i)
            for id in Ensemble.reactor_ids
                this_reactor_data = evolution_outputs[id]
                if :complete_timeseries in outputs
                    this_ts = Dict(i => Ensemble.reactors[id].molecules)
                    this_reactor_data[:complete_timeseries] = merge(this_reactor_data[:complete_timeseries], this_ts)
                end
            
                if :molecule_count in outputs
                    this_count = Dict(i => length(Ensemble.reactors[id].molecules))
                    this_reactor_data[:molecule_count] = merge(this_reactor_data[:molecule_count], this_count)
                end
            
                if :average_length in outputs
                    these_lengths = Ensemble.reactors[id].molecules
                    this_ave = Dict(i => mean(these_lengths))
                    this_reactor_data[:average_length] = merge(this_reactor_data[:average_length], this_ave)
                end
                
                if :var_length in outputs
                    these_lengths = Ensemble.reactors[id].molecules
                    this_var = Dict(i => var(these_lengths))
                    this_reactor_data[:var_length] = merge(this_reactor_data[:var_lengths], this_var)
                end
            end
        end

        # if :complete_timeseries in outputs
        #     evolution_outputs[:complete_timeseries] = complete_ts
        # end

        # if :molecule_count in outputs
        #     evolution_outputs[:molecule_count] = mole_count 
        # end

        # if :average_length in outputs
        #     evolution_outputs[:average_lengths] = ave_lengths 
        #     evolution_outputs[:variance_lengths] = var_lengths
        # end
    end
    return evolution_outputs
end