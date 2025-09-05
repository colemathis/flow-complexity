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



# test

```shell
[1] alexandre@MacBook-Pro  3 - research/quasi-metabolism/quasi-metabolisms.jl/datasets/00_test   main ± 
$ quasi launch 1
Loading parameters for simulation 1 from file ./data/params.csv
Parameters loaded.

   l = 1000
   random_seed = random
   N = 1000
   μ = 1.0e-5
   T = 100
   k = 50
   s = 0.05
   R = <sparse matrix>
   sim_number = 1

Generated random seed: 8117905094621693033
Sim Completed. Time taken: 2.72 seconds.
Data saved in ./data/sims/000001


[1] alexandre@MacBook-Pro  3 - research/quasi-metabolism/quasi-metabolisms.jl/datasets/00_test   main ± 
$ cd ..

[1] alexandre@MacBook-Pro  work/3 - research/quasi-metabolism/quasi-metabolisms.jl/datasets   main ± 
$ cd ..

[1] alexandre@MacBook-Pro  Documents/work/3 - research/quasi-metabolism/quasi-metabolisms.jl   main ± 
$ ll
drwxr-xr-x alexandre staff  96 B  Thu Sep  4 15:31:26 2025  _i-don’t-know-what-that-is
drwxr-xr-x alexandre staff 128 B  Thu Sep  4 15:20:02 2025  analysis
drwxr-xr-x alexandre staff 192 B  Thu Sep  4 15:20:08 2025  datasets
drwxr-xr-x alexandre staff 160 B  Thu Sep  4 15:20:14 2025  research
drwxr-xr-x alexandre staff 160 B  Thu Sep  4 15:20:20 2025  scripts
drwxr-xr-x alexandre staff 352 B  Thu Sep  4 15:56:18 2025  src
.rw-r--r-- alexandre staff 1.0 KB Mon Sep  1 12:29:57 2025  LICENSE
.rw-r--r-- alexandre staff 121 KB Tue Sep  2 19:35:57 2025  Manifest.toml
.rw-r--r-- alexandre staff 1.0 KB Tue Sep  2 19:35:56 2025  Project.toml
.rwxr-xr-x alexandre staff 303 B  Mon Sep  1 13:28:31 2025  quasi
.rw-r--r-- alexandre staff 2.7 KB Wed Sep  3 11:24:30 2025  README.md
```

