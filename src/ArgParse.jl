#==============================================================================#
# IMPORTS
#==============================================================================#


#==============================================================================#
# FUNCTIONS
#==============================================================================#

function print_help()
    """
    Print the help message.
    """

    println("""

Usage: flow <command> [options]

Commands:
  params                 Create params.jl in the current folder
  queue                  Create a queue of jobs to be run in data/params.csv
  slurm                  Create a slurm script to run the simulations
  dry <sim_number>       Run the first 10% iterations and prints calculation time
  launch <sim_number>    Launch simulation sim_number

  test-params            Run commands "params", "queue", and "launch 1"
  test-queue             Run commands "queue", and "launch 1"

""")

end

#==============================================================================#

function parse_arguments()
    """
    Parse command line arguments.
    """

    if isempty(ARGS)
        print_help()
        return
    end

    # Extract command and remaining arguments
    cmd, rest = ARGS[1], ARGS[2:end]
    commands = Dict(
        "params"         => ()      -> create_params_jl_file(),
        "queue"          => ()      -> create_params_csv_file(),
        "slurm"          => ()      -> create_slurm_file(),
        "dry"            => (a...)  -> length(a) == 1 ? launch_simulation(a[1], dry_run=true) : print_help(),
        "launch"         => (a...)  -> length(a) == 1 ? launch_simulation(a[1]) : print_help(),
        "test-params"    => ()      -> begin
            create_params_jl_file()
            create_params_csv_file()
            launch_simulation("1")
        end,
        "test-queue"    => ()      -> begin
        create_params_csv_file()
        launch_simulation("1")
    end,
    )

    # ...and launch the corresponding function
    if haskey(commands, cmd)
        commands[cmd](rest...)
    else
        print_help()
    end

end

#==============================================================================#
# END OF FILE
#==============================================================================#