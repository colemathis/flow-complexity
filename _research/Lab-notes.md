# May 12th 2022

I'm starting a notebook because I keep forgetting what I'm working on. Hopefully this will completement the git logs

Okay today I was able to update the `Ensemble` constructor function, and write a function in `explore-parameters.jl` that will generate an Ensemble and run the dynamics. I need to check for bugs and work on the visulation side now. 

I need to make sure to use `@quickactivate` with DrWatson otherwise saving does not work as expected. 

# June 24th 2022

Okay I should focus on this project a bit more. My goal today is to make a plan of attack. This should include a few different things. Mostly issues on git to guide future development. 

Main goals:
- Parameter sweep of well mixed case:
    - Plot Max MA, MA distribution fit, integer number against (forward, backward rate)
- Parameter sweep of spacially distributed cases 
    - Plot same things against diffusion rates
    - Start with Line Graphs, then lattice, then random.
    - Compute MA of Graphs. (pick forward/backward rate from sweep above and vary diffusion)

To do those we'll need a coherent data management system, and an efficient analysis pipeline.

# July 6th 

Taking a look after a couple days away. The Simulation class is implemented now. I think the main thing to do is convert the .bson files to a .csv (somewhat efficiently) and then write some R scripts to do the analysis we want.

The save system is working, now its time to write the analysis pipeline. 

# July 12th 

Analysis pipeline mostly sorted out. There's likely some bugs lurking about and I don't know what will happen as scale but we'll see. The big problem now is mostly a problem for the future, once specific compounds are stabilized by different reactors it will be difficult to keep track of all the information in coherent ways. The same problem exists for simulations with different rate constants in different reactors

Need to check the assumptions of the simulation.

# July 13th

Lots of systems in place now. Just debugging the simulations before running the parameter sweep. Something that's distrubing is the degree to which the reactors without inflow seem to be completely static for all time. This strikes me as a serious bug. I can't seem to find the source. 

# Sept 9, 2024

- Seems like Cole has fixed the function call in `test_timing.jl`, I have lanched several simulations varying the parameter set and it executes normally.
- I am now trying to get the parameter explorations to work (`test_timing.jl` and `explore-test.jl`), ideally using multiple CPU cores. The original file coded by Cole (`explore-test.jl`) was still looping sequentially over the different parameters which is why I have made a new version (`explore-test.jl`) that computes in advance the different parameter permutations and lanches the simulations using these pre-computed parameters. Seems to use the full extent of the multiprocessing module as seen from my CPU usage in e.g. `htop`.

# Sept 11, 2024

- Did a bunch of tests in the past few days---`test_timing.jl` and `explore-test.jl` are now working on Sol, the GCloud cluster, etc.

# Sept 12, 2024

- start coding a minimal working example: plotting the population time series

