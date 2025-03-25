library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)

# Load the dataset
data <- read_csv("data/timeseries.csv")

# Select 25 random `sim_number`
set.seed(123)  # Ensure reproducibility
selected_sims <- sample(unique(data$sim_number), 25)

# Filter dataset for selected simulations
data_filtered <- data %>%
  filter(sim_number %in% selected_sims)

# Aggregate data: average frequency over `chemostat_id`
data_agg <- data_filtered %>%
  group_by(sim_number, time, integer) %>%
  summarise(frequency = mean(frequency), .groups = 'drop')

# Filter for integers 1 to 10
data_agg <- data_agg %>%
  filter(integer %in% 1:10)

# Create the plot
plot_list <- list()

for (i in 1:25) {
  sim <- selected_sims[i]
  p <- ggplot(data_agg %>% filter(sim_number == sim), aes(x = time, y = frequency, color = factor(integer))) +
    geom_line(alpha = 0.7) +
    labs(title = paste("Simulation", sim), x = "Time", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "none")
  plot_list[[i]] <- p
}

# Arrange plots in a 5x5 grid
pdf("figs/timeseries.pdf", width = 15, height = 15)
do.call(grid.arrange, c(plot_list, ncol = 5, nrow = 5))
dev.off()

