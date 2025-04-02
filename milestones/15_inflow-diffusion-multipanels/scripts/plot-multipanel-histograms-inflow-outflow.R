############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(scales)
library(MASS)
library(purrr)
library(tidyr)

############################
# DATA PROCESSING
############################

nrows <- 10
selected_sims <- 1:(nrows^2)

params <- read.csv("data/params.csv")
processed_data_path <- "data/multipanel-histograms-inflow-outflow/processed_data_200int.csv"

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, integer %in% 1:50) %>%
        filter(time == max(time)) %>%
        # group_by(sim_number, time, integer) %>%
        # reframe(frequency = sum(frequency)) %>%
        # ungroup() %>%
        {write.csv(., processed_data_path, row.names = FALSE); .}
}

############################
# CREATE MULTIPANEL PLOT
############################

# get the parameters and sim data together
merged_data <- merge(params, processed_data, by = "sim_number") %>%
    mutate(inflow_mols = factor(inflow_mols, levels = rev(sort(unique(inflow_mols)))))

# store the logarithms of the integer and frequency
merged_data <- merged_data %>%
    mutate(integer = as.numeric(integer),
           frequency = as.numeric(frequency),
           log_integer = log10(integer),
           log_frequency = log10(frequency)) %>%
    filter(is.finite(log_integer), is.finite(log_frequency))

# fit a power-law to each combination of inflow and outflow
fits <- merged_data %>%
    group_by(inflow_mols, outflow_rate) %>%
    nest() %>%
    mutate(fit = map(data, ~ lm(log_frequency ~ log_integer, data = .x))) %>%
    transmute(inflow_mols, outflow_rate, 
              intercept = map_dbl(fit, ~ coef(.x)[1]), 
              slope = map_dbl(fit, ~ coef(.x)[2]))

# calculate the fitted values from the coefficients
fit_lines <- merged_data %>%
    ungroup() %>% 
    distinct(inflow_mols, outflow_rate, log_integer) %>% 
    inner_join(fits, by = c("inflow_mols", "outflow_rate")) %>%
    mutate(log_integer = as.numeric(log_integer), 
           fitted_y = intercept + slope * log_integer)

p <- ggplot(merged_data, aes(x = log_integer, y = log_frequency)) +
    geom_point(data = merged_data %>% filter(chemostat_id == 1), 
               shape = 4, color = "darkgreen", size = 0.5) +
    geom_point(data = merged_data %>% filter(chemostat_id == 25), 
               shape = 4, color = "darkorange", size = 0.5) +
    # geom_line(data = fit_lines, aes(x = log_integer, y = fitted_y, group = interaction(inflow_mols, outflow_rate)), 
    #           color = "red", linewidth = 0.5) +
    facet_grid(rows = vars(inflow_mols), cols = vars(outflow_rate), scales = "fixed",
                labeller = labeller(.default = function(x) sprintf("%.2f", log10(as.numeric(x))))) +
    theme_bw() +
    theme(legend.position = "none", axis.text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(color = "white"), 
          strip.background.x = element_rect(fill = "blue"), strip.background.y = element_rect(fill = "red")) +
    labs(x = TeX("$log_{10}$ (integer value)"), y = TeX("$log_{10} (frequency)$"),
         title = ("Integer Distributions Inflow/Outflow with Power-Law Fit \n (blue = diffusion, red = inflow)"))

ggsave("figs/multipanel-histograms-inflow-outflow.pdf", plot = p, width = 8, height = 8)
