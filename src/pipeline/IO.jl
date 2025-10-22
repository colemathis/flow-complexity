#==============================================================================#
# IMPORTS
#==============================================================================#

using JLD2

#==============================================================================#
# FUNCTION
#==============================================================================#

function save_adjacency_matrix(sim, fn)
    """
    Save the adjacency matrix to a file.
    """

    # save the adjacency matrix
    A = sim.params[:A]
    @save fn A

end

#==============================================================================#

function load_adjacency_matrix!(sim, fn)
    """
    Load the adjacency matrix from a file.
    """

    # load the adjacency matrix
    @load fn A
    sim.params[:A] = A

end

#==============================================================================#

function save_data(sim)
    """
    Save the simulation data to a file.
    """

    # convert sim number to string
    sim_number = sim.params[:sim_number]
    sim_number_str = sim_number_string(sim_number)

    # determine save path and create the directory if needed
    sim_dir = joinpath("./data/sims/", "$sim_number_str")
    mkpath(sim_dir)

    # save timeseries data
    fn = joinpath(sim_dir, "timeseries.csv")
    CSV.write(fn, sim.output[:timeseries])

    # save graph data
    io_nodes = get_io_nodes(sim)
    fn = joinpath(sim_dir, "graph.csv")
    CSV.write(fn, io_nodes)

    # save metadata
    fn = joinpath(sim_dir, "meta.csv")
    meta = get_metadata(sim)
    CSV.write(fn, meta)

    println("Data saved in $sim_dir")

end

#==============================================================================#

function get_io_nodes(sim)
    """
    Get a DataFrame of the simulation graph.
    """

    sim_number = sim.params[:sim_number]

    # create DataFrame of input/output nodes
    io_nodes = DataFrames.DataFrame(sim_number=Int[], chemostat_in=Int[], chemostat_out=Int[])

    # extract edges from the ensemble graph
    g = sim.ensemble.ensemble_graph
    for e in Graphs.edges(g)
        push!(io_nodes, (
            sim_number = sim_number,
            chemostat_in = Graphs.src(e),
            chemostat_out = Graphs.dst(e)
        ))
    end

    return io_nodes

end

#==============================================================================#

function get_metadata(sim)
    """
    Get a DataFrame of the simulation metadata.
    """

    meta = DataFrames.DataFrame(
        sim_number = sim.params[:sim_number],
        total_time = sim.output[:total_time],
        total_constructive_rxn = get(sim.output, :total_constructive_rxn, missing),
        skipped_constructive_rxn = get(sim.output, :skipped_constructive_rxn, missing),
        total_destructive_rxn = get(sim.output, :total_destructive_rxn, missing),
        skipped_destructive_rxn = get(sim.output, :skipped_destructive_rxn, missing),
        total_diffusion_rxn = get(sim.output, :total_diffusion_rxn, missing),
        skipped_diffusion_rxn = get(sim.output, :skipped_diffusion_rxn, missing),
        total_outflow_rxn = get(sim.output, :total_outflow_rxn, missing),
        skipped_outflow_rxn = get(sim.output, :skipped_outflow_rxn, missing),
    )

    return meta

end

#==============================================================================#
# END OF FILE
#==============================================================================#