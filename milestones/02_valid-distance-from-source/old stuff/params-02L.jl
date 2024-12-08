params_template = DataFrame(
    mass = 1000,
    graph_type = "lattice-2way",
    N_reactors = 25,
    forward_rate = 1e-3,
    outflow_rate = 0.0,
    total_time = 100,
    sim_number = 1,
    save_time_series = true,
    save_parameters = true,
    save_graph = true,
    save_simulation = true,
    save_directory = relative_path,
    notes = "none"
)

# Number of simulations
n = 100

# Customize parameters
vals = exp10.(LinRange(-6,1,n))
