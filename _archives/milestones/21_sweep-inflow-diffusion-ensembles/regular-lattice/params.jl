# nsim = 40
params_template = Dict(
    :mass		=> 0,
    :inflow_mols        => 0,
    :graph_type         => "lattice-2way",
    :randomize_edges    => false,
    :N_reactors     	=> 25,
    :forward_rate   	=> 1e-3,
    :outflow_rate   	=> 0,
    :total_time     	=> 20000,
    :output_time    	=> 20000.0,
    :N_inflow       	=> 1,
    :sim_number     	=> 1,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",
)

# params_array = [copy(params_template) for _ in 1:nsim]

nseries = 100
nrepeat = 10

params_array = [copy(params_template) for _ in 1:nseries]

n = 10

log_inflow_min = 3
log_inflow_max = 6
log_inflow = LinRange(log_inflow_min, log_inflow_max, n)
inflows = exp10.(log_inflow)

log_diffusion_min = -3
log_diffusion_max = -1
log_diffusion = LinRange(log_diffusion_min, log_diffusion_max, n)
diffusions = exp10.(log_diffusion)


for i in 1:n
    for j in 1:n

        ns = (i-1)*n + j
        inflow = inflows[i]
        diffusion = diffusions[j]

        params_array[ns][:inflow_mols] = inflow
        params_array[ns][:outflow_rate] = diffusion

    end
end

# duplicate the params for nrepeat
params_array = vcat((copy(p) for _ in 1:nrepeat for p in params_array)...)

nsim = nseries * nrepeat

for i in 1:nsim
    params_array[i][:sim_number] = i
end