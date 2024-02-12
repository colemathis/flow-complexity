# flow-complexity

This repository contains source code, notebooks, and maybe scripts to analyze how spatial localization and heterogenaity influences the assembly of complex objects


Simulations are based on the type `Chemostat` in the `Chemostat.jl` file. 

# Installation

## Install Julia

See [this link](https://julialang.org/downloads/). As of today, the latest Julia version is v1.10.0).

## Cloning this project

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

# Usage

## Running with multiprocessing enabled

Run Julia using

```
julia --project=. -p [num]
```

where `[num]` represents the number of Julia processes to be launched.

## Dr Watson

Saving and File management is handled via DrWatson, make sure to run `@quickactivate` before running sims.
