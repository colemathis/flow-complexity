module FlowComplexity

include("Simulation.jl")

function parse_args()
    if length(ARGS) == 0
        print_help()
        return
    end

    cmd = ARGS[1]
    if cmd == "queue"
        println("Julia CWD: ", pwd())
        queue()
    elseif cmd == "dry"
        if length(ARGS) < 2
            println("Usage: flow dry <sim_number>")
            return
        end
        sim_number = ARGS[2]
        dry(sim_number)
    elseif cmd == "launch"
        if length(ARGS) < 2
            println("Usage: flow launch <sim_number>")
            return
        end
        sim_number = ARGS[2]
        launch(sim_number)
    elseif cmd == "extract"
        extract()
    else
        print_help()
    end
end

function queue()
    include(joinpath(pwd(), "params.jl"))
    mkpath("data")
    write_params_file(params_array)
end

function dry(sim_number::String)
    launch_simulation(sim_number, dry_run=true)
end

function launch(sim_number::String)
    launch_simulation(sim_number)
end

function extract()
    extract_sims()
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

parse_args()

end
