# FlowComplexity

Investigates how the topology of networks of chemostats influence the complexity of the chemistry.

## Content

- [Installation](#installation)
- [Interactive usage](#interactive-usage)
  - [Loading the module](#loading-the-module)
  - [Running a basic simulation](#running-a-basic-simulation)
- [Usage from the command line](#usage-from-the-command-line)
  - [Creating parameter file](#creating-parameter-file)
  - [Creating queue file](#creating-queue-file)
  - [Running a simulation](#running-a-simulation)

# Installation

1. Install Julia, preferably version `1.11.2` or similar.

2. Clone repo and install requirements.

```shell
git clone https://github.com/colemathis/flow-complexity
cd flow-complexity
julia --project=.
```

then, from the Julia prompt,

```julia
julia> using Pkg
julia> Pkg.instantiate()
```

3. Add folder `flow-complexity` to PATH

e.g., by adding this line to `.bashrc` or `.zshrc`

```shell
export PATH="$PATH:[path to repo]/flow-complexity"
```

# Usage

## Workflow and command line usage

The code is designed to work both in command-line (shell) mode and in interactive mode (e.g., Jupyter Notebook).

You can launch the code using the `flow` command line, followed by one of these commands:

- `params` creates a parameter file named `params.jl` in the current directory
- `queue` creates a csv file named `data/params.csv` containing the parameters for an ensemble of simulations to be executed
- `slurm` creates a SLURM script file named `run.slurm` in the current directory
- `launch <sim number>` launches simulation number `<sim_number>` from the queue file

## Usage in interactive mode

In interactive mode, load Julia using

```bash
cd flow-complexity
julia --project=.
```

then load the `FlowComplexity` module with

```julia
julia> include("src/FlowComplexity.jl")
```

A `Simulation` object can be created using

```julia
s = FlowComplexity.Simulation();
```

Then it can be launched with either

```julia
julia> FlowComplexity.run_simulation(s);
```

or

```julia
julia> s.run()
```

A simulation from an array can be launched using

```julia
julia> s.from_sim_array(1);
```

Data is automatically saved in `./data/<sim_number>`. It can be accessed either from the object itself,

```julia
using Plots

df = s.output[:timeseries]

plt = plot()
for i in 1:10
    d = df[(df.chemostat_id .== 1) .& (df.integer .== i), :]
    plot!(plt, d.time, d.frequency, label="int $i")
end

display(plt)
```

or from the files written to disk

```julia
using CSV, DataFrames, Plots

df = CSV.read("./data/sims/000001/timeseries.csv", DataFrame)

plt = plot()
for i in 1:10
    d = df[(df.chemostat_id .== 1) .& (df.integer .== i), :]
    plot!(plt, d.time, d.frequency, label="int $i")
end

display(plt)
```

## Additional commands

The `A` matrix can be saved/loaded using

```julia
julia> s.save_R_matrix("R-matrix.jld2")
julia> s.load_R_matrix("R-matrix.jld2")
```