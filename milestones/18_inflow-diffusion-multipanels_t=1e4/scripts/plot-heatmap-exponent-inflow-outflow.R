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
library(reshape2)

############################
# DATA PROCESSING
############################

nrows <- 10
selected_sims <- 1:(nrows^2)
nsims <- nrows^2

params <- read.csv("data/params.csv")
processed_data_path <- "data/heatmap-exponent-inflow-outflow/processed_data.csv"

processed_data <- if (file.exists(processed_data_path)) {
    read.csv(processed_data_path)
} else {
    read.csv("data/timeseries.csv") %>%
        filter(sim_number %in% selected_sims, integer %in% 1:50) %>%
        filter(time == max(time)) %>%
        # group_by(sim_number, time, integer) %>%
        # reframe(frequency = sum(frequency)) %>%
        # ungroup() %>%
        {dir.create(dirname(processed_data_path), recursive = TRUE, showWarnings = FALSE); 
        write.csv(., processed_data_path, row.names = FALSE); .}
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
fits_inflow <- merged_data %>% filter(chemostat_id == 1) %>% 
    group_by(inflow_mols, outflow_rate) %>%
    nest() %>%
    mutate(fit = map(data, ~ lm(log_frequency ~ log_integer, data = .x))) %>%
    transmute(inflow_mols, outflow_rate, 
              intercept = map_dbl(fit, ~ coef(.x)[1]), 
              slope = map_dbl(fit, ~ coef(.x)[2]))

# calculate the fitted values from the coefficients
# fit_lines_inflow <- merged_data %>%
#     ungroup() %>% 
#     distinct(inflow_mols, outflow_rate, log_integer) %>% 
#     inner_join(fits_inflow, by = c("inflow_mols", "outflow_rate")) %>%
#     mutate(log_integer = as.numeric(log_integer), 
#            fitted_y = intercept + slope * log_integer)

# fit a power-law to each combination of inflow and outflow
fits_outflow <- merged_data %>% filter(chemostat_id == 25) %>% 
    group_by(inflow_mols, outflow_rate) %>%
    nest() %>%
    mutate(fit = map(data, ~ lm(log_frequency ~ log_integer, data = .x))) %>%
    transmute(inflow_mols, outflow_rate, 
              intercept = map_dbl(fit, ~ coef(.x)[1]), 
              slope = map_dbl(fit, ~ coef(.x)[2]))

# calculate the fitted values from the coefficients
# fit_lines_outflow <- merged_data %>%
#     ungroup() %>% 
#     distinct(inflow_mols, outflow_rate, log_integer) %>% 
#     inner_join(fits_outflow, by = c("inflow_mols", "outflow_rate")) %>%
#     mutate(log_integer = as.numeric(log_integer), 
#            fitted_y = intercept + slope * log_integer)

############################
# HEATMAP OF POWER-LAW EXPONENTS
############################

# slope of power-law fits (i.e., k in x^{-k}, so they’re positive)
k_inflow <- -fits_inflow$slope
k_outflow <- -fits_outflow$slope

heatmap_matrix <- matrix(k_outflow-k_inflow, nrow = nrows, ncol = nrows)
# heatmap_matrix <- -heatmap_matrix

heatmap_plot <- ggplot(melt(heatmap_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "black") +
    labs(x = TeX("$log_{10} (k_d)$"), y = TeX("$log_{10} (I)$"), fill = TeX("$Delta k$")) +
    theme_bw() +
    scale_x_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$outflow_rate[seq(1, nrows, by = 1)])), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:nrows, labels = sprintf("%.2f", log10(params$inflow_mols[seq(1, nsims, by = nrows)])), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

# print(heatmap_plot)

# Linear scale
heatmap_plot_linear <- heatmap_plot +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    labs(title = TeX("Difference between exponents of power law fits ($k_{outflow}-k_{inflow}$)"))

ggsave("figs/heatmap-exponents-inflow-outflow-linear.pdf", plot = heatmap_plot_linear, width = 8, height = 7)

# Log scale
# heatmap_plot_log <- heatmap_plot +
#     scale_fill_gradientn(colors = c("blue", "white", "red"), trans = "log", labels = scales::scientific) +
#     labs(title = "(Exponent of Outflow Power-law Fit) - (Exponent of Inflow Power-law Fit)")

# ggsave("figs/heatmap-exponents-inflow-outflow-log.pdf", plot = heatmap_plot_log, width = 8, height = 7)
