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
# dir.create(dirname(processed_data_path), recursive = TRUE, showWarnings = FALSE)

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, time == max_time) %>%
        group_by(sim_number) %>%
        summarise(mean_integer = sum(integer * frequency) / sum(frequency), .groups = "drop") %>%
        {dir.create(dirname(processed_data_path), recursive = TRUE, showWarnings = FALSE); 
        write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

heatmap_matrix <- matrix(processed_data$mean_integer, nrow = nrows, ncol = nrows)

heatmap_plot <- ggplot(melt(heatmap_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "black") +
    # scale_fill_gradientn(colors = c("blue", "white", "red")) +
    labs(x = TeX("$log_{10} (k_d)$"), y = TeX("$log_{10} (I)$"), fill = "Mean Integer") +
    theme_bw() +
    scale_x_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$outflow_rate[seq(1, nrows, by = 1)])), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$inflow_mols[seq(1, nsims, by = nrows)])), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

# Linear scale
heatmap_plot_linear <- heatmap_plot +
    scale_fill_distiller(palette = "Spectral") +
    labs(title = "Integer Average Across Simulations \n (Linear Scale Colormap)")

# print(heatmap_plot_linear)
ggsave("figs/heatmap-avg-integer-linear.pdf", plot = heatmap_plot_linear, width = 8, height = 7)

# Log scale
heatmap_plot_log <- heatmap_plot +
    scale_fill_distiller(palette = "Spectral", trans = "log", labels = scales::scientific) +
    labs(title = "Integer Average Across Simulations \n (Log Scale Colormap)")

# print(heatmap_plot_log)
ggsave("figs/heatmap-avg-integer-log.pdf", plot = heatmap_plot_log, width = 8, height = 7)

# print(heatmap_plot)
