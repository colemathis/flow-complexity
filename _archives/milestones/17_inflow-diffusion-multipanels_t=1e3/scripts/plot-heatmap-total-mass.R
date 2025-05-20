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

processed_data_path <- sprintf("data/heatmap-total-mass/processed_data.csv")

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, time == max_time) %>%
        group_by(sim_number) %>%
        summarise(total_mass = sum(integer * frequency), .groups = "drop") %>%
        {dir.create(dirname(processed_data_path), recursive = TRUE, showWarnings = FALSE); 
        write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

heatmap_matrix <- matrix(processed_data$total_mass, nrow = nrows, ncol = nrows)

heatmap_plot <- ggplot(melt(heatmap_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "black") +
    labs(x = TeX("$log_{10} (k_d)$"), y = TeX("$log_{10} (I)$"), fill = "Total Mass") +
    theme_bw() +
    scale_x_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$outflow_rate[seq(1, nrows, by = 1)])), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$inflow_mols[seq(1, nsims, by = nrows)])), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

# print(heatmap_plot)

# ggsave("figs/heatmap-total-mass.pdf", plot = heatmap_plot, width = 8, height = 7)

# Linear scale
heatmap_plot_linear <- heatmap_plot +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Total Mass Across Simulations \n (Linear Scale Colormap)")

ggsave("figs/heatmap-total-mass-linear.pdf", plot = heatmap_plot_linear, width = 8, height = 7)

# Log scale
heatmap_plot_log <- heatmap_plot +
    scale_fill_distiller(palette = "Spectral", trans = "log", labels = scales::scientific) +
    labs(title = "Total Mass Across Simulations \n (Log Scale Colormap)")

ggsave("figs/heatmap-total-mass-log.pdf", plot = heatmap_plot_log, width = 8, height = 7)
