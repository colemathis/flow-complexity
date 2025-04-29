module FlowComplexity

include("Simulation.jl")

function queue()
    include(joinpath(pwd(), "params.jl"))
    mkpath("data")
    write_params_file(params_template, nrepeat)
end

function dry(sim_number)
    launch_simulation(sim_number, dry_run=true)
end

function launch(sim_number)
    launch_simulation(sim_number)
end

function extract()
    extract_sims2()
end

function print_help()
    println("""
Usage: flow <command> [options]
Commands:
  queue                  Create a queue of jobs to be run in data/params.csv
  dry <sim_number>       Run the first 10% iterations of simulation sim_number
                         and prints calculation time estimates
  launch <sim_number>    Launch simulation sim_number
  extract                Extract data in data/sims and save it in data/sim_array.jld2
Examples:
  flow queue
  flow dry 99
  flow launch 99
  flow extract
""")
end

# --- Command dispatcher ---
function main()
    if isempty(ARGS)
        print_help()
        return
    end

    cmd, rest = ARGS[1], ARGS[2:end]
    commands = Dict(
        "queue"   => ()      -> queue(),
        "dry"     => (a...)  -> length(a) == 1 ? dry(a[1]) : print_help(),
        "launch"  => (a...)  -> length(a) == 1 ? launch(a[1]) : print_help(),
        "extract" => ()      -> extract()
    )

    if haskey(commands, cmd)
        commands[cmd](rest...)
    else
        print_help()
    end
end

main()

end