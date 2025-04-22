nsim = 100
params_template = Dict(
    :mass		=> 0,
    :inflow_mols        => 0,
    :graph_type         => "lattice-2way",
    :randomize_edges    => true,
    :N_reactors     	=> 25,
    :forward_rate   	=> 1e-3,
    :outflow_rate   	=> 0,
    :total_time     	=> 10000,
    :output_time    	=> 100.0,
    :N_inflow       	=> 1,
    :sim_number     	=> 1,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> 1,
)

params_array = [copy(params_template) for _ in 1:nsim]

n = sqrt(nsim)
n = Int(n)

log_inflow_min = 3
log_inflow_max = 6
log_inflow = LinRange(log_inflow_min, log_inflow_max, n)
inflows = exp10.(log_inflow)

log_diffusion_min = -3
log_diffusion_max = -1
log_diffusion = LinRange(log_diffusion_min, log_diffusion_max, n)
diffusions = exp10.(log_diffusion)

# vals = exp10.(LinRange(3,6,nsim))
# vals = LinRange(1, 100, 10)
# vals = round.(Int, vals)*10

using Distributions

for i in 1:n
    for j in 1:n

        ns = (i-1)*n + j
        inflow = inflows[i]
        diffusion = diffusions[j]

        params_array[ns][:sim_number] = ns
        params_array[ns][:inflow_mols] = inflow
        params_array[ns][:outflow_rate] = diffusion

    end
end
