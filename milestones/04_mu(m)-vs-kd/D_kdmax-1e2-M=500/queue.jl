using DrWatson
@quickactivate :FlowComplexity

using DataFrames

# Define directories and paths
current_dir = pwd()
milestones_dir = joinpath(projectdir(), "milestones")
relative_path = relpath(current_dir, milestones_dir)

# Define a configuration template
params_template = DataFrame(
    mass = 500,
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

# Generate DataFrame with n rows based on the template
params_df = vcat([params_template for _ in 1:n]...)

# Update the sim_number for each row
params_df.sim_number .= 1:n

# Customize parameters
vals = exp10.(LinRange(-6,2,n))
for i in 1:n
    params_df[i, :outflow_rate] = vals[i]
end

# Save the DataFrame to a CSV file
array_fn = "./data/array.csv"
save(array_fn, params_df)
