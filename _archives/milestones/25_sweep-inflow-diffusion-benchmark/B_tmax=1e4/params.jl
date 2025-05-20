nrepeat = 1

inflow_mols = exp10.(LinRange(2,6,10))
diffusion_rates = exp10.(LinRange(-6,3,10))

params_template = Dict(
    # simulation parameters
    :save_interval    	=> 10,
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
    :inflow_mols        => inflow_mols,
    :forward_rate   	=> 1e-3,
    :diffusion_rate   	=> diffusion_rates,
    :outflow_rate       => "equal-to-diffusion-rate",
)
