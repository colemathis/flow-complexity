nsim = 100
params_template = Dict(
    :mass           => 0,
    :inflow_mols    => 0,
    :graph_type     => "lattice-2way",
    :N_reactors     => 25,
    :forward_rate   => 1e-3,
    :outflow_rate   => 0,
    :total_time     => 100,
    :output_time    => 1.0,
    :N_inflow       => 1,
    :sim_number     => 1,
    :method         => "tau-leaping",
    :dt             => 1e-3,
    :random_seed    => 1,
)

params_array = [copy(params_template) for _ in 1:nsim]

log_inflow_min = 3
log_inflow_max = 6

log_diffusion_min = -3
log_diffusion_max = -1

# vals = exp10.(LinRange(3,6,nsim))
# vals = LinRange(1, 100, 10)
# vals = round.(Int, vals)*10

using Distributions

for i in 1:nsim

    log_inflow = rand(Uniform(log_inflow_min, log_inflow_max))
    log_diffusion = rand(Uniform(log_diffusion_min, log_diffusion_max))
    
    inflow = exp10(log_inflow)
    diffusion = exp10(log_diffusion)

    params_array[i][:sim_number] = i
    params_array[i][:inflow_mols] = inflow
    params_array[i][:outflow_rate] = diffusion
    
end
