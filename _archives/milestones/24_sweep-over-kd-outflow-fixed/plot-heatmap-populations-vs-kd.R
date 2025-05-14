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

TITLE        <- "Frequency vs Diffusion Rate"
ID           <- "populations-heatmap-vs-kd"
USE_CACHE    <- TRUE
TIME         <- 1e5
OUTFLOW_ONLY <- TRUE

DATA_DIR  <- "../23_distance-from-source-outflow-fixed/D_tmax=1e5/data"
CACHE_DIR <- paste0("cache/", ID)
FIGS_FILE <- ID
if (OUTFLOW_ONLY) FIGS_FILE <- paste0(ID, "-outflow-only")
FIGS_DIR  <- paste0("figs")

TIMESERIES_ARROW <- file.path(DATA_DIR, "timeseries.arrow")
PARAMS_CSV       <- file.path(DATA_DIR, "params.csv")
GRAPHS_CSV       <- file.path(DATA_DIR, "graphs.csv")

# ASSEMBLY_CSV <- "Assembly-10000.csv"

############################
# FUNCTIONS
############################

load_processed_data <- function() {
    cache_path <- file.path(CACHE_DIR, sprintf("%s.csv", ID))

    if (file.exists(cache_path) && USE_CACHE) {
        read_csv(cache_path, show_col_types = FALSE)
    } else {
        data <- open_dataset(TIMESERIES_ARROW, format = "arrow") %>%
            filter(time == TIME) %>%
            collect()

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}

pre_process_data <- function(ts, params) {

  # Attach diffusion_rate from params to each row of the time‑series data
  ts <- ts %>%
    left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

  ts <- ts %>%
    filter(integer > 1)

  if (OUTFLOW_ONLY) {
    ts <- ts %>%
      filter(chemostat_id == 25)
  }

  # Bin integer and aggregate frequency for heatmap
  bin_edges <- unique(round(exp(seq(log(1), log(1e7), length.out = 51))))
  heatmap_data <- ts %>%
    mutate(integer_bin = cut(
      integer,
      breaks = bin_edges,
      include.lowest = TRUE,
      right = TRUE
    )) %>%
    group_by(diffusion_rate, integer_bin) %>%
    summarise(frequency = sum(frequency), .groups = "drop")

  # Complete the dataset: for non-existent frequencies assign 0
  heatmap_data <- heatmap_data %>%
    complete(diffusion_rate, integer_bin, fill = list(frequency = 0))

  heatmap_data
}

create_main_plot <- function(ts, params) {

  p_main <- ggplot(ts, aes(x = diffusion_rate, y = integer_bin, fill = frequency)) +
    geom_tile() +
    scale_x_log10() +
    scale_fill_viridis_c(name = "Freq.", na.value = "grey") +
    labs(
      x = TeX("Diffusion coefficient"),
      y = "Integer",
    ) +
    scale_y_discrete(
      breaks = levels(ts$integer_bin)[seq(1, length(levels(ts$integer_bin)), length.out = 8)],
      labels = c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7))
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0.95, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
      legend.key.size = unit(0.5, "lines"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )

  if (OUTFLOW_ONLY) {
    p_main <- p_main +
      labs(title = paste(TITLE, "(outflow only)"))
  } else {
    p_main <- p_main +
      labs(title = paste(TITLE, "(whole system)"))
  }

  # add a white line with slope ~ -1
  p_main <- p_main +
    geom_segment(
      aes(x = 1e-6, xend = 1e1, y = levels(ts$integer_bin)[length(levels(ts$integer_bin))], yend = levels(ts$integer_bin)[1]),
      color = "white",
      inherit.aes = FALSE
  )

  p_main
}

############################
# MAIN SCRIPT
############################

ts <- load_processed_data()
params <- read_csv(PARAMS_CSV, show_col_types = FALSE)
ts <- pre_process_data(ts, params)
p_main  <- create_main_plot(ts, params)

width <- 120
height <- 70

options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
print(p_main)

out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
ggsave(filename = out_file, plot = p_main, width = width, height = height, units = "mm", create.dir = TRUE)