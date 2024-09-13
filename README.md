# flow-complexity

This repository contains source code, notebooks, and maybe scripts to analyze how spatial localization and heterogenaity influences the assembly of complex objects


Simulations are based on the type `Chemostat` in the `Chemostat.jl` file. 

## Installation

### Install Julia

See [this link](https://julialang.org/downloads/). As of today, the latest Julia version is v1.10.0.

### Cloning this project

Clone this repository using

```
git clone https://github.com/ELIFE-ASU/flow-complexity
```

then install the dependencies using

```
cd flow-complexity
julia --project=.
```

type the `]` key to enter pkg mode, then `instantiate`.

Test using:

```
include("scripts/explore-parameters.jl")
run_all_topologies()
```

## Basic Usage

### Minimal Working Example

[to be completed]

### Test timing

Use `test_timing.jl`. This will write output to the file `timing_results.csv` located in the `data` directory (make sure the file exists).

### Explore parameters

The simulations will be save under `data/sims`.

To put together the output data, use `get-parameter-sweep-stats.r` which reads everything under `data/sims` and writes to something like `2022_07_22_parameter_sweep_stats.csv`.

[to be completed]

## Advanced Usage

...

### Running with multiprocessing enabled

For example, to run `explore-test.jl` from the `scripts` folder on 220 CPU,

```
julia -p 220 --project=.. explore-test.jl
```

# Analysis

[to be completed]

## Dr Watson

Saving and File management is handled via DrWatson, make sure to run `@quickactivate` before running sims.
