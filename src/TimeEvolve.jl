"""

##########################################################################

TIMEEVOLVE



##########################################################################

"""

using StatsBase
using Statistics

# include("Chemostats.jl")
# include("Ensemble.jl")
#p3 to remove?

"""

##########################################################################
EVOLVE WELL MIXED #p1 broken !

Time evolution of a single well mixed chemostat.
##########################################################################

Input:

    chemostat       (Chemostat) = chemostat struct
    tau_max         (float)     = total simulation time
    output_freq     (float)     = interval at which data is recorded
    outputs         (array)     = variables to record #p1 what’s "symbol" here?

Output:

    final_output    (dict)      = simulation output

##########################################################################

"""

function evolve_well_mixed(chemostat::Chemostat,
    tau_max::Float64,
    output_freq::Float64,
    outputs::Array{Symbol,1}
)

    # create simulation output dictionary
    evolution_outputs = Dict{Symbol,Any}()

    # create dictionary for recorded variables
    if :complete_timeseries in outputs
        complete_ts = Dict{Int64,Array{Int64,1}}()
    end

    if :molecule_count in outputs
        mole_count = Dict{Int64,Int64}()
    end

    if :average_length in outputs
        ave_lengths = Dict{Int64,Float64}()
        var_lengths = Dict{Int64,Float64}()
    end

    # calculate the propensities for the chemostat
    propensities = calc_propensities(chemostat)

    # set the time to zero, and the saving checkpoint to zero too
    tau = 0.0
    checkpoint = 0.0

    # evolve the simulation
    while tau < tau_max

        # pick a reaction
        rxn = sample(["construction", "degradation", "outflow"], Weights(propensities))

        if rxn == "construction"
            chemostat = constructive_rxn(chemostat)
        elseif rxn == "degradation"
            chemostat = destructive_rxn(chemostat)
        elseif rxn == "outflow"
            chemostat, outflow_direction = outflow_rxn(chemostat)
        end

        #p1 commenting this out for now
        # calculate the total mass of the molecules
        chemostat.mass = calc_mass(chemostat)
        # if the mass is different from the fixed mass and the fixed mass has been defined,
        # then add mass 1 molecules to the chemostat
        if chemostat.fixed_mass != chemostat.mass
            delta_mass = chemostat.fixed_mass - chemostat.mass
            new_moles = repeat([1], delta_mass)
            chemostat.molecules = vcat(chemostat.molecules, new_moles)
        end

        # update propensities
        propensities = calc_propensities(chemostat)

        # update time 
        total_p = sum(propensities)

        # calculate the time step and define the new time #p1 what is this formula? gillespie? to check
        tau_step = -log(rand()) / total_p
        tau += tau_step

        # record the output
        if tau > checkpoint

            # record the time
            i = round(tau, digits=3)

            # set the next checkpoint
            checkpoint += output_freq

            # record the other variables
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

        # p2 is this duplicated (outputs vs evolution_outputs) because of performance issues?
        # integrate the data into the output dictionary 
        if :complete_timeseries in outputs
            evolution_outputs[:complete_timeseries] = complete_ts
        end

        if :molecule_count in outputs
            evolution_outputs[:molecule_count] = mole_count
        end

        if :average_length in outputs
            evolution_outputs[:average_lengths] = ave_lengths
            evolution_outputs[:var_lengths] = var_lengths
        end
    end

    # get everything together
    final_output = Dict(1 => evolution_outputs)

    # and return the data
    return final_output
end

"""

##########################################################################
EVOLVE DISTRIBUTED SYSTEM

Time evolution of a distributed system, i.e. an ensemble of chemostats
structured in a specific topology (graph).
##########################################################################

Input:

    Ensemble            (Ensemble)      = ensemble to evolve
    tau_max             (float)         = total simulation time
    output_freq         (float)         = interval at which data is recorded
    outputs             (symbol array)  = variables to record #p1 symbol?
    seed                (int)           = random seed

Output:

    evolution_outputs   (dict)          = simulation data

##########################################################################

"""

function evolve_distributed(Ensemble,
    tau_max::Float64,
    output_freq::Float64,
    outputs::Array{Symbol,1},
    seed::Int64=1337
)

    # set random seed
    Random.seed!(seed)

    # create a dictionary holding the data of all reactors
    evolution_outputs = Dict{Int64,Any}()

    # for each reactor, create a dictionary holding the variables and data
    for id in Ensemble.reactor_ids
        this_reactor_data = Dict{Symbol,Any}()
        for out in outputs
            this_reactor_data[out] = Dict{Any,Any}()
        end
        evolution_outputs[id] = this_reactor_data
    end

    # set the time to zero
    tau = 0.0

    # set the next frame framed to zero too
    checkpoint = 0.0

    # calculate reaction propensities
    all_propensities = calc_propensities.(Ensemble.reactors)
    # returns a vector of vector containing the propensities for each reactor

    # use these propensities to calculate the next reaction times
    next_taus = calc_next_rxn_times(all_propensities, tau)

    switch = false #p3 not sure about this one

    # main simulation loop
    while tau < tau_max

        # determine next reaction (and its reactor) from next reaction times
        next_rxn = popfirst!(next_taus)
        next_reactor = next_rxn[2]
        tau = next_rxn[1]

        # determine the propensities for the reactor who will react next
        next_reactor_propensities = all_propensities[next_reactor]

        # pick a reaction from the propensities
        rxn = sample(["construction", "degradation", "outflow"], Weights(next_reactor_propensities))

        # if switch & (next_reactor != 1) #p3 not sure what’s this
        #     println(next_reactor, " \t", rxn)
        # end

        # get the chemostat for this next reactor
        this_chemostat = Ensemble.reactors[next_reactor]

        # execute reaction and update propensities
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
        push!(next_taus, (chemostat_next_tau, next_reactor))
        sort!(next_taus, by=x -> x[1])
        # Check Total 1s
        for chemostat in Ensemble.reactors
            if chemostat.fixed_mass != 0
                chemostat.mass = calc_mass(chemostat)
                chemo_ones = calc_ones(chemostat)
                if chemo_ones != chemostat.fixed_mass
                    delta_ones = chemostat.fixed_mass - chemo_ones
                    if delta_ones > 0
                        new_moles = repeat([1], delta_ones)
                        chemostat.molecules = vcat(chemostat.molecules, new_moles)
                    end
                end
            end
        end

        # record outputs
        if tau > checkpoint

            # record the time
            i = round(tau, digits=3)

            # set the next record checkpoint
            checkpoint += output_freq

            # for the current time frame, record the data for every reactor in the ensemble
            for id in Ensemble.reactor_ids

                # p1: how is this_reactor_data put back in evolution_outputs?
                this_reactor_data = evolution_outputs[id]
                # println(countmap(Ensemble.reactors[id].molecules))

                # add the time series data (molecules) to the data dictionary
                if :complete_timeseries in outputs
                    this_ts = Dict(i => Ensemble.reactors[id].molecules)
                    this_reactor_data[:complete_timeseries] = merge(this_reactor_data[:complete_timeseries], this_ts)
                end

                # add the molecule count to the data dictionary
                if :molecule_count in outputs
                    this_count = Dict(i => length(Ensemble.reactors[id].molecules))
                    this_reactor_data[:molecule_count] = merge(this_reactor_data[:molecule_count], this_count)
                end

                # add the average of the molecules length to the data dictionary
                if :average_length in outputs
                    these_lengths = Ensemble.reactors[id].molecules
                    this_ave = Dict(i => mean(these_lengths))
                    this_reactor_data[:average_length] = merge(this_reactor_data[:average_length], this_ave)
                end

                # add the variance of the molecule length to the data dictionary
                if :var_length in outputs
                    these_lengths = Ensemble.reactors[id].molecules
                    this_var = Dict(i => var(these_lengths))
                    this_reactor_data[:var_length] = merge(this_reactor_data[:var_length], this_var)
                end
            end
        end

    end

    # return the simulation data
    return evolution_outputs

end
