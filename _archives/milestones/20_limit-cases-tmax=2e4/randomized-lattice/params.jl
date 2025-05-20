# nsim = 40
params_template = Dict(
    :mass		=> 0,
    :inflow_mols        => 0,
    :graph_type         => "lattice-2way",
    :randomize_edges    => true,
    :N_reactors     	=> 25,
    :forward_rate   	=> 1e-3,
    :outflow_rate   	=> 0,
    :total_time     	=> 20000,
    :output_time    	=> 200.0,
    :N_inflow       	=> 1,
    :sim_number     	=> 1,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",
)

# params_array = [copy(params_template) for _ in 1:nsim]

nseries = 4
nrepeat = 10

params_array = [copy(params_template) for _ in 1:nseries]

log_inflow_min = 3
log_inflow_max = 6
inflow_min = exp10(log_inflow_min)
inflow_max = exp10(log_inflow_max)

log_diffusion_min = -3
log_diffusion_max = -1
diffusion_min = exp10(log_diffusion_min)
diffusion_max = exp10(log_diffusion_max)

params_array[1][:inflow_mols] = inflow_min
params_array[1][:outflow_rate] = diffusion_min

params_array[2][:inflow_mols] = inflow_min
params_array[2][:outflow_rate] = diffusion_max

params_array[3][:inflow_mols] = inflow_max
params_array[3][:outflow_rate] = diffusion_min

params_array[4][:inflow_mols] = inflow_max
params_array[4][:outflow_rate] = diffusion_max

# duplicate the params for nrepeat
params_array = vcat((copy(p) for _ in 1:nrepeat for p in params_array)...)

nsim = nseries * nrepeat

for i in 1:nsim
    params_array[i][:sim_number] = i
end