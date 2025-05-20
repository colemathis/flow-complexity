# Load libraries
library(readr)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

# Define simulation number
sim <- 1

# Read and prepare graph
g <- read_csv("data/graphs.csv") %>%
  filter(sim_number == sim) %>%
  dplyr::select(chemostat_in, chemostat_out) %>%
  rename(from = chemostat_in, to = chemostat_out) %>%
  graph_from_data_frame(directed = TRUE) %>%
  as_tbl_graph()

# Create and save plot
p <- ggraph(g, layout = "stress") +
  geom_edge_link(arrow = arrow(length = unit(4, "mm")), end_cap = circle(3, "mm")) +
  geom_node_point(size = 4, color = "steelblue") +
  geom_node_text(aes(label = name, color = "red"), repel = TRUE) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(paste0("figs/graph_", sim, ".pdf"), plot = p, width = 6, height = 6, dpi = 300)
