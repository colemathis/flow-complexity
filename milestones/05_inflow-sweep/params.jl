nsim = 10
params_template = Dict(
    :mass           => 0,
    :inflow_mols    => 0,
    :graph_type     => "lattice-2way",
    :N_reactors     => 25,
    :forward_rate   => 1e-3,
    :outflow_rate   => 1e-3,
    :total_time     => 100,
    :output_time    => 1.0,
    :N_inflow       => 1,
    :sim_number     => 1,
    :method         => "tau-leaping",
    :dt             => 1e-3,
    :random_seed    => 1,
)

params_array = [copy(params_template) for _ in 1:nsim]

vals = exp10.(LinRange(3,6,nsim))
# vals = LinRange(1, 100, 10)
# vals = round.(Int, vals)*10

for i in 1:nsim
    params_array[i][:sim_number] = i
    params_array[i][:inflow_mols] = vals[i]
end