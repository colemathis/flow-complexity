nsim = 1
params_template = Dict(
    :mass           => 1000,
    :graph_type     => "line",
    :N_reactors     => 1,
    :forward_rate   => 1e-3,
    :outflow_rate   => 0.0,
    :total_time     => 100,
    :output_time    => 1.0,
    :N_inflow       => 1,
    :sim_number     => 1,
    :method         => "tau-leaping",
    :dt             => 0.01,
    :random_seed    => 1,
)

params_array = [copy(params_template) for _ in 1:nsim]

# vals = exp10.(LinRange(-6,2,n))
# vals = LinRange(1, 100, 10)
# vals = round.(Int, vals)*10

for i in 1:nsim
    params_array[i][:sim_number] = i
    # params_array[i][:outflow_rate] = vals[i]
end