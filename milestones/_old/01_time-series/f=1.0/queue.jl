using DrWatson
@quickactivate :FlowComplexity

using DataFrames

# Define directories and paths
current_dir = pwd()
milestones_dir = joinpath(projectdir(), "milestones")
relative_path = relpath(current_dir, milestones_dir)

# Define a configuration template
params_template = DataFrame(
    mass = 1000,
    graph_type = "line",
    N_reactors = 1,
    forward_rate = 1.0,
    outflow_rate = 5.0,
    total_time = 10,
    sim_number = 1,
    save_time_series = true,
    save_parameters = true,
    save_graph = true,
    save_simulation = true,
    save_directory = relative_path,
    notes = "none"
)

# Number of simulations
n = 10

# Generate DataFrame with n rows based on the template
params_df = vcat([params_template for _ in 1:n]...)

# Update the sim_number for each row
params_df.sim_number .= 1:n

# Save the DataFrame to a CSV file
array_fn = "./data/array.csv"
save(array_fn, params_df)
