############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(reshape2)

############################
# DATA PROCESSING
############################

params <- read.csv("data/params.csv")

nrows <- 10
nsims <- nrows^2
selected_sims <- 1:nsims
max_time <- params$total_time[params$sim_number == 1]

processed_data_path <- sprintf("data/heatmap-avg-integer/processed_data.csv")

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, time == max_time) %>%
        group_by(sim_number) %>%
        summarise(mean_integer = sum(integer * frequency) / sum(frequency), .groups = "drop") %>%
        {write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

heatmap_matrix <- matrix(processed_data$mean_integer, nrow = nrows, ncol = nrows)

heatmap_plot <- ggplot(melt(heatmap_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(x = TeX("$\\log_{10} k_d$"), y = TeX("$\\log I$"), fill = "Mean Integer") +
    theme_minimal() +
    scale_x_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$outflow_rate[seq(1, nrows, by = 1)]))) +
    scale_y_reverse(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$inflow_mols[seq(1, nsims, by = nrows)]))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(aes(label = rep(processed_data$sim_number)))

print(heatmap_plot)

ggsave("figs/heatmap.pdf", plot = heatmap_plot, width = 7, height = 6)
