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
        {dir.create(dirname(processed_data_path), recursive = TRUE, showWarnings = FALSE); 
        write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

merged_data <- merge(params, processed_data, by = "sim_number") %>%
    mutate(inflow_mols = factor(inflow_mols, levels = rev(sort(unique(inflow_mols)))))

p <- ggplot(merged_data, aes(x = time, y = frequency, color = factor(integer))) +
    geom_line(alpha = 0.7) +
    facet_grid(rows = vars(inflow_mols), cols = vars(outflow_rate), scales = "fixed",
               labeller = labeller(.default = function(x) sprintf("%.2f", log10(as.numeric(x))))) +
    scale_y_log10(labels = label_scientific()) +
    theme_bw() +
    theme(legend.position = "none", axis.text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(color = "white"), 
          strip.background.x = element_rect(fill = "blue"), strip.background.y = element_rect(fill = "red")) +
    labs(x = TeX("time"), y = TeX("frequency"),
         title = ("Time Series of Simulations \n (blue = diffusion, red = inflow)"))

ggsave("figs/multipanel-timeseries.pdf", plot = p, width = 8, height = 8)
