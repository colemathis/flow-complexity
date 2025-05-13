nrepeat = 1

inflow_rates = exp10.(LinRange(2,4,4))
diffusion_rates = exp10.(LinRange(-1,3,100))

params_template = Dict(
    # simulation parameters
    :save_interval    	=> 1e2,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",

    # physical parameters
    :total_time     	=> 1e4,
    :initial_mass       => 0,

    # topological parameters
    :graph_type         => "lattice-2way",
    :randomize_edges    => false,
    :N_reactors     	=> 25,

    # reaction rates
    :inflow_mols        => inflow_rates,
    :forward_rate   	=> 1e-3,
    :diffusion_rate   	=> diffusion_rates,
    :outflow_rate       => "equal-to-diffusion-rate",
)