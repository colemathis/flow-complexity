# flow-complexity

This repository contains source code, notebooks, and maybe scripts to analyze how spatial localization and heterogenaity influences the assembly of complex objects


Simulations are based on the type `Chemostat` in the `Chemostat.jl` file. 

## Installation

1. Use `juliaup status` to see which julia version is installed, then remove/rename as needed $HOME/.julia*
1. Install julia using `juliaup add 1.10.2`
1. From the repo folder, open julia using `julia --project=.`
1. Install required packages with `using Pkg; Pkg.instantiate()` (this takes some time)
1. Exit julia then open it again (without any argument)
1. Install DrWatson using `Pkg.add(name="DrWatson", version="2.17.0")` then pin this version to the home environment using `Pkg.pin(name="DrWatson", version="2.17.0")`.

If stuck in dependency hell, update packages using `Pkg.update()` then commit the updated `Manifest.toml`

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

