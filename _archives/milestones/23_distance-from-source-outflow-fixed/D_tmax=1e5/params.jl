nrepeat = 1

nsim = 100
diffusion_rates = exp10.(LinRange(-6,1,nsim))

params_template = Dict(
    # simulation parameters
    :save_interval    	=> 1e3,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",

    # physical parameters
    :total_time     	=> 1e5,
    :initial_mass       => 0,

    # topological parameters
    :graph_type         => "lattice-2way",
    :randomize_edges    => false,
    :N_reactors     	=> 25,

    # reaction rates
    :inflow_mols        => 1e3,
    :forward_rate   	=> 1e-3,
    :diffusion_rate   	=> diffusion_rates,
    :outflow_rate       => "equal-to-diffusion-rate",
)
