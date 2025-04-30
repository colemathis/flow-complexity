nrepeat = 1

params_template = Dict(
    # simulation parameters
    :save_interval    	=> 1,
    :method         	=> "tau-leaping",
    :dt             	=> 1e-3,
    :random_seed    	=> "random",

    # physical parameters
    :total_time     	=> 100,
    :initial_mass       => 0,

    # topological parameters
    :graph_type         => "lattice-2way",
    :randomize_edges    => false,
    :N_reactors     	=> 25,

    # reaction rates
    :inflow_mols        => 1e3,
    :forward_rate   	=> 1e-3,
    :diffusion_rate   	=> 1e-1,
)
