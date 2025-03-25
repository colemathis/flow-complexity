############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(grid)

############################
# DATA PROCESSING
############################

# Set the simulation number and read the parameters
selected_sim <- 70                                                          # select the simulation number
params <- read.csv("data/params.csv")                                       # read the params file

# Process the data, or read the processed data if it exists
file_path <- sprintf("data/single-timeseries/single_%d.csv", selected_sim)  # path to processed data
if (file.exists(file_path)) {                                               # check if the processed data exists
    processed_data <- read.csv(file_path)                                    # read the processed data
} else {                                                                    # if the processed data does not exist
    data <- read.csv("data/timeseries.csv")                                 # read the timeseries data
    processed_data <- data %>%                                               # filter the data
        filter(sim_number == selected_sim, integer %in% 1:10)
    write.csv(processed_data, file_path, row.names = FALSE)                  # save the processed data
}

############################
# PLOT THE FIGURE
############################

# Get basic parameters
n_chemostat <- params$N_reactors[params$sim_number == selected_sim]                             # number of chemostats
nrows <- sqrt(n_chemostat)                                                                      # number of rows of the figure
sim_params <- params %>% filter(sim_number == selected_sim)                                     # get the simulation parameters

# Plot the time series as a multipanel figure
p <- ggplot(processed_data, aes(x = time, y = frequency, color = factor(integer))) +             # plot the time series
    geom_line(alpha = 0.7) +                                                                    # add lines
    facet_wrap(~ chemostat_id, ncol = nrows, nrow = nrows, scales = "free_y",                   # create grid plot
                labeller = labeller(chemostat_id = function(x) paste("chemostat #", x))) +      # add plot titles
    xlab(NULL) +                                                                                # remove x-axis label
    ylab(NULL) +                                                                                # remove y-axis label
    theme_bw() +                                                                                # set theme
    theme(legend.position = "none",                                                             # remove legend
            panel.border = element_rect(color = "black", fill = NA))                            # add border around plots

# Add title
t1 <- sprintf("simulation %d", selected_sim)                                                    # add simulation number
t2 <- sprintf("$log_{10}(I)=%.2f$", log10(sim_params$inflow_mols))                              # add inflow rate
t3 <- sprintf("$log_{10}(k_d)=%.2f$", log10(sim_params$outflow_rate))                           # add diffusion coefficient
title <- paste0(t1, ": ", t2, ", ", t3)                                                         # put them together
title <- TeX(title)                                                                             # convert to LaTeX
p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))                         # add to plot and center

# Save the figure
# print(p)
ggsave("figs/single-timeseries.pdf", plot = p, width = 8, height = 8)
