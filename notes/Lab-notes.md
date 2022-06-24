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