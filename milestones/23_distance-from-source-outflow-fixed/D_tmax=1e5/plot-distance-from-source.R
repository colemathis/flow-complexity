############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(grid)
library(igraph)      # for graph‑based distance calculations

############################
# PARAMETERS
############################

TITLE        <- "Distance from source"
ID           <- "distance-from-source"
USE_CACHE   <- TRUE

DATA_DIR  <- "data"
CACHE_DIR <- paste0("cache/", ID)
FIGS_FILE <- ID
FIGS_DIR  <- paste0("figs")

TIMESERIES_ARROW <- file.path(DATA_DIR, "timeseries.arrow")
PARAMS_CSV       <- file.path(DATA_DIR, "params.csv")
GRAPHS_CSV       <- file.path(DATA_DIR, "graphs.csv")

ASSEMBLY_CSV <- "Assembly-10000.csv"

############################
# FUNCTIONS
############################

load_processed_data <- function() {
    cache_path <- file.path(CACHE_DIR, sprintf("distance-from-source.csv"))

    if (file.exists(cache_path) && USE_CACHE) {
        read_csv(cache_path, show_col_types = FALSE)
    } else {
        data <- open_dataset(TIMESERIES_ARROW, format = "arrow") %>%
            filter(sim_number %in% c(30, 60, 100)) %>%
            collect() %>%
            filter(time == max(time))

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}
  
append_distance_data <- function(ts,
                                 graphs_csv = GRAPHS_CSV) {
    library(igraph)  # ensure igraph is available inside non‑interactive calls
    
    graphs <- read.csv(graphs_csv)
    sim_numbers <- unique(ts$sim_number)
    
    dist_tbl <- lapply(sim_numbers, function(s) {
        g_edges <- graphs %>%
            filter(sim_number == s) %>%
            select(chemostat_in, chemostat_out)
        
        # skip if graph is missing
        if (nrow(g_edges) == 0) return(NULL)
        
        # build an **undirected** graph so distance is symmetric
        g <- graph_from_data_frame(g_edges, directed = FALSE)
        
        # ensure vertex "1" (the source chemostat) exists
        if (!"1" %in% V(g)$name) return(NULL)
        
        d_vec <- distances(g, v = "1")[1, ]
        
        tibble(
            sim_number   = s,
            chemostat_id = as.integer(names(d_vec)),
            distance     = as.numeric(d_vec)
        )
    }) %>% bind_rows()
    
    ts %>% left_join(dist_tbl, by = c("sim_number", "chemostat_id"))
}


join_assembly_index <- function(ts,
                                assembly_csv = ASSEMBLY_CSV,
                                missing_value = 17) {
    ai <- read.csv(assembly_csv)
    ts %>%
        left_join(ai, by = "integer") %>%
        mutate(assemblyindex = ifelse(is.na(assemblyindex), missing_value, assemblyindex))
}

calculate_weighted_stats <- function(ts) {
    ts %>%
        group_by(sim_number, chemostat_id) %>%
        summarize(
            mean_ai = sum(assemblyindex * frequency) / sum(frequency),
            sd_ai   = sqrt(sum(frequency * (assemblyindex - mean(assemblyindex))^2) / sum(frequency)),
            .groups = "drop"
        )
}

create_main_plot <- function(ts) {
  p <- ggplot(ts, aes(x = distance, y = mean_ai, color = factor(sim_number))) +
    geom_point(alpha = 0.40) +
    geom_smooth(data = ts %>% filter(sim_number == 30),
                method = "loess", span = 1.0, se = FALSE) +
    geom_smooth(data = ts %>% filter(sim_number == 60),
                method = "loess", span = 0.75, se = FALSE) +
    geom_smooth(data = ts %>% filter(sim_number == 100),
                method = "loess", span = 4.20, se = FALSE) +
    labs(
      x = TeX("Distance from source $d$"),
      y = "Mean Assembly Index",
      color = "sim_number"
    ) +
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 6.0, ymax = 8.5,
             color = "grey", fill = NA, linetype = "dashed") +
    annotate("text", x = 6, y = 12, label = TeX("$k_d = 10^{-4}$"),
             hjust = 0, size = 3, color = scales::hue_pal()(3)[1]) +
    annotate("text", x = 6, y = 11, label = TeX("$k_d = 10^{-2}$"),
             hjust = 0, size = 3, color = scales::hue_pal()(3)[2]) +
    annotate("text", x = 6, y = 10, label = TeX("$k_d = 10^{1}$"),
             hjust = 0, size = 3, color = scales::hue_pal()(3)[3]) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p
}

create_inset_plot <- function(ts) {
  # create inset plot using data from a specific chemostat and distance
  ts_inset <- ts %>% filter(sim_number == 60, distance == 4)
  p_inset <- ggplot(ts_inset, aes(x = integer, y = frequency, color = factor(chemostat_id))) +
    scale_y_log10() +
    geom_density(aes(y = ..scaled..), adjust = 2, alpha = 1.0) +
    scale_color_manual(values = c("#0000FF", "#2222FF", "#4444FF", "#6666FF", "#9999FF")) +
    labs(
      x = "Integer",
      y = "Density",
      color = "chemostat_id"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  p_inset
}

combine_plots <- function(main_plot, inset_plot,
                          xmin = 0.75, xmax = 6.5,
                          ymin = 1.25, ymax = 5.5) {
  inset_grob <- ggplotGrob(inset_plot)
  combined <- main_plot +
    annotation_custom(
      grob = rectGrob(
        gp = gpar(fill = "white", col = "black", lwd = 1)
      ),
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    ) +
    annotation_custom(
      grob = inset_grob,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  combined
}

############################
# MAIN SCRIPT
############################

# Load and process the data
ts <- load_processed_data()

ts_inset <- ts
ts_inset <- append_distance_data(ts_inset)

ts <- join_assembly_index(ts)
ts <- calculate_weighted_stats(ts)
ts <- append_distance_data(ts)

# Create the main and inset plots
p_main  <- create_main_plot(ts)
p_inset <- create_inset_plot(ts_inset)

# Combine plots into a final output plot
combined_plot <- combine_plots(p_main, p_inset)

# display the combined plot

options(vsc.dev.args = list(width = 80, height = 70, res=300, units = "mm"))
print(combined_plot)

out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
ggsave(filename = out_file, plot = combined_plot, width = 80, height = 70, units = "mm", create.dir = TRUE)