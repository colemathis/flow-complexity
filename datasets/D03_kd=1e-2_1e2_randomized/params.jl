nrepeat = 1

params_template = Dict(
    # simulation parameters
    :save_interval    	=> 100,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",

    # physical parameters
    :total_time     	=> 1e4,
    :initial_mass       => 0,

    # topological parameters
    :graph_type         => "lattice-2way",
    :randomize_edges    => true,
    :N_reactors     	=> 25,

    # reaction rates
    :inflow_mols        => 1e4,
    :forward_rate   	=> 1e-3,
    :diffusion_rate   	=> exp10.(LinRange(-2,2,100)),
    :outflow_rate       => "equal-to-diffusion-rate",
)
