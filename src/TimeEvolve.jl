using StatsBase
using Statistics
include("Chemostats.jl")

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