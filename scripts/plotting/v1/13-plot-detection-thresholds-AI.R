############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(grid)
library(igraph)      # for graph‑based distance calculations
library(arrow)      # for reading Arrow files

############################
# PARAMETERS
############################

TITLE        <- "Detection thresholds (AI)"
ID           <- "detection-thresholds-ai"
USE_CACHE   <- TRUE

DATA_DIR  <- "../23_distance-from-source-outflow-fixed/D_tmax=1e5/data"
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
    cache_path <- file.path(CACHE_DIR, sprintf("%s.csv", ID))

    if (file.exists(cache_path) && USE_CACHE) {
        read_csv(cache_path, show_col_types = FALSE)
    } else {
        data <- open_dataset(TIMESERIES_ARROW, format = "arrow") %>%
            filter(time == 1e5) %>%
            collect()

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}

join_assembly_index <- function(ts,
                                assembly_csv = ASSEMBLY_CSV,
                                missing_value = 17) {
    ai <- read.csv(assembly_csv)
    ts %>%
        left_join(ai, by = "integer") %>%
        mutate(assemblyindex = ifelse(is.na(assemblyindex), missing_value, assemblyindex))
}

attach_diffusion_rate <- function(ts, params) {
  ts %>%
    left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")
}

calculate_detection_threshold <- function(ts, threshold) {
    # Compute total number of molecules in each chemostat
    ts_totals <- ts %>%
        group_by(diffusion_rate, chemostat_id) %>%
        summarise(total_molecules = sum(frequency), .groups = "drop")

    # Attach totals and compute concentration of each integer
    ts_conc <- ts %>%
        left_join(ts_totals,
                  by = c("diffusion_rate", "chemostat_id")) %>%
        mutate(concentration = frequency / total_molecules)

    # For each chemostat, find the highest integer whose concentration
    # exceeds the supplied threshold
    detectable <- ts_conc %>%
        filter(concentration > threshold) %>%
        group_by(diffusion_rate, chemostat_id) %>%
        summarise(highest_AI = max(assemblyindex), .groups = "drop")

    # For each diffusion rate, choose the highest of those integers
    detection_thresholds <- detectable %>%
        group_by(diffusion_rate) %>%
        summarise(max_detected_AI = max(highest_AI, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(max_detected_AI = ifelse(
            is.infinite(max_detected_AI), NA_integer_, max_detected_AI))

    # Return a dataframe with one row per diffusion_rate
    detection_thresholds
}

create_main_plot <- function(detection_thresholds) {

  p_main <- ggplot(detection_thresholds,
              aes(
                  x = diffusion_rate, 
                  y = max_detected_AI, 
                  color = factor(threshold)
              )) +
    geom_point(size = 0.5, alpha = 0.25) +
    geom_smooth(method = "loess", span = 0.25, se = FALSE) +
    scale_x_log10() +
    # scale_y_log10() +
    labs(
      x = TeX("Diffusion coefficient"),
      y = "Highest AI detected",
      color = "Threshold"
    ) +
  theme(
    legend.justification=c(0,0), 
    legend.position=c(0.05,0.55),
    legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
    legend.key.size = unit(0.15, "cm"),    # Decrease key size
    legend.text = element_text(size = 8), # Decrease text size
    legend.title = element_text(size = 8) # Decrease title size
  )
  
  p_main
}

############################
# MAIN SCRIPT
############################

ts <- load_processed_data()
params <- read_csv(PARAMS_CSV, show_col_types = FALSE)
ts <- join_assembly_index(ts)
ts <- attach_diffusion_rate(ts, params)

detection_thresholds <- list(
    # `1e-7` = calculate_detection_threshold(ts, 1e-7),
    # `1e-6` = calculate_detection_threshold(ts, 1e-6),
    # `1e-5` = calculate_detection_threshold(ts, 1e-5),
    `1e-4` = calculate_detection_threshold(ts, 1e-4),
    `1e-3` = calculate_detection_threshold(ts, 1e-3),
    `1e-2` = calculate_detection_threshold(ts, 1e-2),
    `1e-1` = calculate_detection_threshold(ts, 1e-1)
) %>% 
  bind_rows(.id = "threshold")
  
p_main  <- create_main_plot(detection_thresholds)

options(vsc.dev.args = list(width = 80, height = 70, res=300, units = "mm"))
print(p_main)

out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
ggsave(filename = out_file, plot = p_main, width = 80, height = 70, units = "mm", create.dir = TRUE)