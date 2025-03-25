############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(dplyr)
library(ggplot2)
library(latex2exp)

############################
# DATA PROCESSING
############################

selected_sim <- 70  
params <- read.csv("data/params.csv")
file_path <- sprintf("data/single-timeseries/single_%d.csv", selected_sim)

processed_data <- if (file.exists(file_path)) {
    read.csv(file_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number == selected_sim, integer %in% 1:10) %>%
        {write.csv(., file_path, row.names = FALSE); .}
}

grid_size <- sqrt(params$N_reactors[params$sim_number == selected_sim])
sim_params <- params %>% filter(sim_number == selected_sim)

############################
# CREATE MULTIPANEL PLOT
############################

p <- ggplot(processed_data, aes(x = time, y = frequency, color = factor(integer))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ chemostat_id, ncol = grid_size, nrow = grid_size, scales = "free_y",
               labeller = labeller(chemostat_id = function(x) paste("chemostat #", x))) +
    theme_bw() +
    theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))

plot_title <- TeX(sprintf("Simulation %d: $log_{10}(I)=%.2f$, $log_{10}(k_d)=%.2f$", 
                           selected_sim, log10(sim_params$inflow_mols), log10(sim_params$outflow_rate)))

p <- p + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))

ggsave("figs/single-timeseries.pdf", plot = p, width = 8, height = 8)
