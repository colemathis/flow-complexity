################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(grid)
library(igraph)      # for graph‑based distance calculations
library(arrow)      # for reading Arrow files

################################################################################
# PARAMETERS
################################################################################

TITLE        <- "Standard deviation across sweep over diffusion, multiple inflows"
ID           <- "populations-sd-vs-kd-multiple-inflows"
USE_CACHE    <- TRUE
TIME         <- 1e4
N_REACTORS   <- 25
SAVEFIGS     <- TRUE

DATA_DIR     <- "data"
CACHE_DIR    <- file.path("cache", ID)
CACHE_PATH   <- file.path(CACHE_DIR, paste0(ID, ".csv"))
FIGS_FILE    <- ID
FIGS_DIR     <- file.path("figs")

PARAMS_CSV   <- file.path(DATA_DIR, "params.csv")
GRAPHS_CSV   <- file.path(DATA_DIR, "graphs.csv")

################################################################################
# FUNCTIONS
################################################################################

load_params <- function() {
  read_csv(PARAMS_CSV, show_col_types = FALSE)
}

#==============================================================================#

load_and_process_data <- function() {

  timeseries_files <- list.files(DATA_DIR, pattern = "timeseries.csv", recursive = TRUE, full.names = TRUE)

  ts_all <- map_dfr(timeseries_files, function(file) {
    ts <- read_csv(file, show_col_types = FALSE, progress = FALSE)
    process_data(ts)
  }, .progress = TRUE)

  return(ts_all)
}

#==============================================================================#

process_data <- function(ts) {

  ts <- ts %>%
    filter(time == TIME) %>%
    filter(integer == 2)

  ts <- ts %>%
    left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
    left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

  ts %>%
    group_by(diffusion_rate, inflow_mols, integer) %>%
    summarize(
      mean_frequency = sum(frequency, na.rm = TRUE) / N_REACTORS,
      # sd_frequency = sd(frequency),
      # compute the standard deviation with using sd()
      sd_frequency = sqrt(sum((frequency - mean_frequency)^2) / (N_REACTORS - 1)),
      .groups = "drop"
    )

}

#==============================================================================#

save_cached_data <- function(ts_all) {
  dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
  write_csv(ts_all, CACHE_PATH)
}

#==============================================================================#

load_cached_data <- function() {
  read_csv(CACHE_PATH, show_col_types = FALSE)
}

#==============================================================================#

plot_figure <- function(ts) {

  p_main <- ggplot(ts,
              aes(
                  x = diffusion_rate, 
                  y = sd_frequency, 
                  color = factor(log10(inflow_mols))
              )) +
    geom_point(size = 0.5, alpha = 0.25) +
    geom_smooth(method = "loess", span = 0.5, se = FALSE, size = 0.5) +
    scale_x_log10() +
    scale_color_manual(values = c("darkblue", "blue", "yellow", "red")) +
    # scale_color_viridis_d(option = "A") + # Use Viridis color scale for cold-to-hot effect
    labs(
      x = TeX("Diffusion coefficient"),
      y = "Std. deviation (whole systems)",
      color = TeX("Log(I)"),
      title = TITLE,
    ) +
  theme(
    legend.justification=c(0,1), 
    legend.position=c(0.05,0.85),
    legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
    legend.key.size = unit(0.15, "cm"),    # Decrease key size
    # legend.text = element_text(size = 6), # Decrease text size
    # legend.title = element_text(size = 6) # Decrease title size
    plot.title = element_text(size = 6)
  )

}

#==============================================================================#

save_figure <- function(p_main) {

  height <- 70
  width <- 80

  options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
  print(p_main)

  if (SAVEFIGS) {
    out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
    ggsave(filename = out_file, plot = p_main, width = width, height = height, units = "mm", create.dir = TRUE)
  }

}

################################################################################
# MAIN SCRIPT
################################################################################

if (!file.exists(CACHE_PATH) || !USE_CACHE) {
  params <- load_params()
  ts <- load_and_process_data()
  save_cached_data(ts)
} else {
  ts <- load_cached_data()
}

p <- plot_figure(ts)
save_figure(p)

