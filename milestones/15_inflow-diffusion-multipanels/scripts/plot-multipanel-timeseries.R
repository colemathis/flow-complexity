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
nsims <- nrows^2
selected_sims <- 1:nsims
params <- read.csv("data/params.csv")

# Process the data, or read the processed data if it exists
file_path <- sprintf("data/multipanel-timeseries/processed_data.csv")       # path to processed data
if (file.exists(file_path)) {                                               # check if the processed data exists
    processed_data <- read.csv(file_path)                                   # read the processed data
} else {                                                                    # if the processed data does not exist
    data <- read.csv("data/timeseries.csv")                                 # read the timeseries data
    # Filter the dataset based on the selected simulations and integer range
    processed_data <- data %>%
        filter(sim_number %in% selected_sims, integer %in% 1:10)
    # Aggregate data: sum frequency over `chemostat_id`
    processed_data <- processed_data %>%
        group_by(sim_number, time, integer) %>%
        summarise(frequency = sum(frequency), .groups = "drop")
    write.csv(processed_data, file_path, row.names = FALSE)                  # save the processed data
}

############################
# PLOT THE FIGURE
############################

merged_data <- merge(params, processed_data, by = "sim_number")

merged_data$inflow_mols <- factor(merged_data$inflow_mols, 
                                    levels = sort(unique(merged_data$inflow_mols), 
                                    decreasing = TRUE))

# Create the plot using facet_wrap
p <- ggplot(merged_data, aes(x = time, y = frequency, color = factor(integer))) +
    geom_line(alpha = 0.7) +

    facet_grid(
        rows = vars(inflow_mols),
        cols = vars(outflow_rate),
        scales = "free_y", 
        labeller = labeller(
                inflow_mols = function(x) sprintf("%.2f", log10(as.numeric(x))),
                outflow_rate = function(x) sprintf("%.2f", log10(as.numeric(x)))
        )
    ) +
    theme_bw() +
    theme(
        strip.text = element_text(color = "blue")
    ) +
    labs(x = TeX("$log_{10} (k_d)$"), y = TeX("$log_{10} (I)$")) +
    theme(
        axis.title.x = element_text(color = "blue"),
        axis.title.y = element_text(color = "blue")
    ) +
    labs(title = "Time Series of Simulations", color = "Integer") +
    theme(
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)
    ) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.border = element_rect(color = "black", fill = NA))

# Save the figure
ggsave("figs/multipanel-timeseries.pdf", plot = p, width = 8, height = 8)







    # Use facet_wrap to create a multi-panel plot
    # facet_wrap(~ sim_number, scales = "free_y", 
    #     labeller = labeller(sim_number = function(x) {
    #     sim_params <- params %>% filter(sim_number == x)
    #     sprintf("I=%.2e \n k_d=%.2e", sim_params$inflow_mols, sim_params$outflow_rate)
    #     })
    # ) +
    # theme(
    #     strip.background = element_blank(),
    #     strip.placement = "outside",
    #     strip.text = element_text(size = 10)
    # ) +
    # labs(x = "Inflow Mols (I)", y = "Outflow Rate (k_d)") +

        # Use facet_wrap to create a multi-panel plot
    # facet_grid(rows = vars(sim_params$inflow_mols), cols = vars(sim_params$outflow_rate), scales = "free_y") +
    # theme(
    #     strip.background = element_blank(),
    #     strip.placement = "outside",
    #     strip.text = element_text(size = 10)
    # ) +
    # labs(x = "Time", y = "Frequency") +

    # Annotate each panel with the sim_number
    # geom_text(data = summarised_data %>% distinct(sim_number), 
    #           aes(x = mean(range(summarised_data$time)), y = mean(range(summarised_data$frequency)), label = sim_number), 
    #           color = "black", size = 10, hjust = -0.5, vjust = 1.5)