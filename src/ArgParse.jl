#==============================================================================#
# IMPORTS
#==============================================================================#


#==============================================================================#
# FUNCTIONS
#==============================================================================#

function print_help()
    println("""

Usage: flow <command> [options]

Commands:
  params                 Create params.jl in the current folder
  queue                  Create a queue of jobs to be run in data/params.csv
  slurm                  Create a slurm script to run the simulations
  dry <sim_number>       Run the first 10% iterations and prints calculation time
  launch <sim_number>    Launch simulation sim_number
  extract                Extract data in data/sims

  test                   Run commands "params", "queue", "launch" and "extract" in sequence
                         for testing purposes

""")
end

#==============================================================================#

function parse_arguments()
    if isempty(ARGS)
        print_help()
        return
    end

    cmd, rest = ARGS[1], ARGS[2:end]
    commands = Dict(
        "params"  => ()      -> create_params_jl_file(),
        "queue"   => ()      -> create_params_csv_file(),
        "slurm"   => ()      -> create_slurm_file(),
        "dry"     => (a...)  -> length(a) == 1 ? launch_simulation(a[1], dry_run=true) : print_help(),
        "launch"  => (a...)  -> length(a) == 1 ? launch_simulation(a[1]) : print_help(),
        "extract" => ()      -> extract_sims(),
        "test"    => ()      -> begin
            create_params_jl_file()
            create_params_csv_file()
            launch_simulation("1")
            extract_sims()
        end,
    )

    if haskey(commands, cmd)
        commands[cmd](rest...)
    else
        print_help()
    end
end

#==============================================================================#
# END OF FILE
#==============================================================================#