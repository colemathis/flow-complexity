############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)

############################
# DATA PROCESSING
############################

nrows <- 10
selected_sims <- 1:(nrows^2)

params <- read.csv("data/params.csv")
processed_data_path <- "data/multipanel-histograms/processed_data.csv"

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, integer %in% 1:10) %>%
        filter(time == max(time)) %>%
        group_by(sim_number, time, integer) %>%
        reframe(frequency = sum(frequency)) %>%
        ungroup() %>%
        {write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

p <- ggplot(processed_data, aes(x = factor(integer), y = frequency)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
    scale_y_log10() +
    labs(title = "Simulation Results", x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ sim_number, labeller = labeller(sim_number = function(x) {
        params %>% filter(sim_number == x) %>%
        reframe(label = sprintf("I=%.2e \n kd=%.2e", inflow_mols, outflow_rate)) %>%
        pull(label)
    })) +
    geom_text(aes(label = sim_number), x = Inf, y = Inf, hjust = 1.5, vjust = 1.5, size = 10, color = "grey", alpha = 0.5)

ggsave("figs/multipanel-histograms.pdf", plot = p, width = 20, height = 20)
