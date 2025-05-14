using DrWatson
@quickactivate :FlowComplexity

using JLD2

# Define the directory containing the simulations
directory = "./data/sims"

# Get the number of subfolders (i.e., simulations)
nsim = count(isdir, readdir(directory, join=true))

# Initialize an array to hold the simulation objects
sim_array = Vector{FlowComplexity.Simulation}(undef, nsim)

# Load each simulation and store it in the array
for i in 1:nsim
    sim_number_string = lpad(i, 6, "0")
    sim_path = joinpath(directory, sim_number_string, "simulation.jld2")
    sim_array[i] = load(sim_path)["sim"]
end

# Save the array of simulations to a file
output_file = joinpath("./data", "data.jld2")
@save output_file sim_array