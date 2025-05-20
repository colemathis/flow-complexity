require(tidyverse)

source("analyze-timeseries.r")

all.sim.paths <- Sys.glob(file.path("../data/sims/*"))

all.sim.parameters <- data.frame()
for (p in all.sim.paths){
    sim.num = as.numeric(str_split(p, "/")[[1]][4])
    print(sim.num)
    if (sim.num > 20){ # Totally ad-hoc
        this.sim.params <- read.csv(file.path(p, "parameters.csv"))
        AI.stats <- get_AI_stats_from_sim(sim.num)
        fit.stats <- get_exp_fit_from_sim(sim.num)

        this.sim.params$alpha = fit.stats$alpha
        this.sim.params$alpha.plus = fit.stats$alpha.plus
        this.sim.params$alpha.minus = fit.stats$alpha.minus
        this.sim.params$r2 = fit.stats$r2

        this.sim.params$max.AI <- AI.stats$max.AI
        this.sim.params$median.AI <- AI.stats$median.AI
        this.sim.params$weighted.mean.AI <- AI.stats$weighted.mean.AI

        all.sim.parameters <- rbind(all.sim.parameters, this.sim.params)
    }

    

}

write.csv(all.sim.parameters, "../data/2022_07_22_parameter_sweep_stats.csv")