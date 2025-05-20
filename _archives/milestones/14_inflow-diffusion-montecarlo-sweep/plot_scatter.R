library(ggplot2)
library(dplyr)

# Read the data
params <- read.csv("./data/params.csv")
timeseries <- read.csv("./data/timeseries.csv")
assembly_data <- read.csv("./data/Assembly-10000.csv")

# Filter for the highest value of "time"
max_time <- max(timeseries$time)
filtered_timeseries <- timeseries %>% filter(time == max_time)

# Compute the average integer value for each sim_number
avg_integer_values <- filtered_timeseries %>%
  group_by(sim_number) %>%
  summarize(avg_integer = mean(integer * frequency, na.rm = TRUE))

# Merge with params data
params <- params %>% left_join(avg_integer_values, by = "sim_number")

# Look up the assembly index based on avg_integer matching integer in assembly_data
params <- params %>% left_join(assembly_data, by = c("avg_integer" = "integer"))

# Create the first scatter plot with avg_integer
p1 <- ggplot(params, aes(x = inflow_mols, y = outflow_rate, color = avg_integer)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Inflow (mols)", y = "Outflow Rate", color = "Avg Integer") +
  theme_minimal()

# Save the first plot
ggsave("figs/scatter_integers.pdf", p1, width = 8, height = 6)

# Create the second scatter plot with assembly index
p2 <- ggplot(params, aes(x = inflow_mols, y = outflow_rate, color = assemblyindex)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Inflow (mols)", y = "Outflow Rate", color = "Assembly Index") +
  theme_minimal()

# Save the second plot
ggsave("figs/scatter_AI.pdf", p2, width = 8, height = 6)

# Display the plots
print(p1)
print(p2)
