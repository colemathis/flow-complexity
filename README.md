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

# Running simulations using the SLURM scheduler

To execute an ensemble of simulations simultaneously on the SLURM scheduler on ASU’s high performance computing cluster, the first step is to create a SLURM job file using

```shell
flow slurm
```

This creates a file named `run.slurm` in the current directory.

You can edit the file to adjust the parameters. For example, if your `params.jl` defines an array of 100 simulations, the following line will let SLURM know to execute simulations 1 to 100:

```shell
#SBATCH --array=1-100                    # simulation indices
```

Other parameters that you can customize include wall time (the time interval after which SLURM will kill a simulation that has not completed) and memory (the amount of RAM memory allocated to each simulation):

```shell
#SBATCH --time=24:00:00                  # wall time (hh:mm:ss)
#SBATCH --mem=4G                         # memory
```

You should not have to customize other parameters, or modify the rest of the script.

After you have made sure the parameters in `slurm.run` are set correctly, you can send your job description to SLURM using

```shell
sbatch slurm.run
```

This will add your list of simulations to the queue. You can confirm the job has been added using the command `myjobs`.

The priority with which they will be run depends on your "fairshare score" — the more you use the cluster, the longer it will take before your jobs are executed.

Refer to the [page on Sol](https://github.com/mathis-group/wiki/wiki/HPC-Cluster) in the wiki for more information about the cluster.