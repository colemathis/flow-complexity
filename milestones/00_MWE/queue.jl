using DrWatson
@quickactivate :FlowComplexity

using DataFrames

current_dir = pwd()
milestones_dir = joinpath(projectdir(), "milestones")
relative_path = relpath(current_dir, milestones_dir)

# define a configuration template
params_df = DataFrame(
    mass = 1000,
    graph_type = "line",
    N_reactors = 1,
    forward_rate = 1.0e-2,
    outflow_rate = 5.0,
    total_time = 1.0,
    sim_number = 1,
    save_time_series = true,
    save_parameters = true,
    save_graph = true,
    save_simulation = true,
    save_directory = relative_path,
    notes = "dummy task"
)

# dummy example: copy the template n times and change the sim number
n = 100
first_row = params_df[1, :]
new_rows = [first_row for _ in 1:n]
params_df = DataFrame(new_rows)

for i in 1:n
    params_df[i, :sim_number] = i
end

array_fn = "./data/array.csv"
save(array_fn, params_df) ;
