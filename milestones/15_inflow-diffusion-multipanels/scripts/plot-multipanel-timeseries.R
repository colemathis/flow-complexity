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

selected_sims <- 1:100
params <- read.csv("data/params.csv")
processed_data_path <- "data/multipanel-timeseries/processed_data.csv"

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, integer %in% 1:10) %>%
        group_by(sim_number, time, integer) %>%
        summarise(frequency = sum(frequency), .groups = "drop") %>%
        {write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

merged_data <- merge(params, processed_data, by = "sim_number") %>%
    mutate(inflow_mols = factor(inflow_mols, levels = rev(sort(unique(inflow_mols)))))

p <- ggplot(merged_data, aes(x = time, y = frequency, color = factor(integer))) +
    geom_line(alpha = 0.7) +
    facet_grid(rows = vars(inflow_mols), cols = vars(outflow_rate), scales = "free_y",
               labeller = labeller(.default = function(x) sprintf("%.2f", log10(as.numeric(x))))) +
    theme_bw() +
    theme(axis.text = element_text(size = 6), legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = TeX("$log_{10} (k_d)$"), y = TeX("$log_{10} (I)$"),
         title = "Time Series of Simulations", color = "Integer")

ggsave("figs/multipanel-timeseries.pdf", plot = p, width = 8, height = 8)
