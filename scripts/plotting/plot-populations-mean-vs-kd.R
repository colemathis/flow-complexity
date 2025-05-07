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

TITLE        <- "Populations across sweep over diffusion"
ID           <- "populations-mean-vs-kd"
USE_CACHE   <- TRUE

DATA_DIR  <- "../23_distance-from-source-outflow-fixed/D_tmax=1e5/data"
CACHE_DIR <- paste0("cache/", ID)
FIGS_FILE <- ID
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
            filter(time == 1e5) %>%
            collect()

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}

create_main_plot <- function(ts, params) {
  # Attach diffusion_rate from params to each row of the time‑series data
  ts <- ts %>%
    left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

  ts <- ts %>% 
  complete(
    chemostat_id,         # every chemostat
    integer,              # every integer that already exists in ts
    diffusion_rate,       # every diffusion_rate that already exists in ts
    fill = list(frequency = 0)   # value to use for the rows you create
  )

  p_main <- ggplot(ts %>% filter(integer >= 2, integer <= 10) %>%
            group_by(integer, diffusion_rate) %>%
            summarize(
              mean_frequency = mean(frequency),
              sd_frequency = sd(frequency),
              .groups = "drop"
            ),
              aes(
                  x = diffusion_rate, 
                  y = mean_frequency, 
                  color = factor(integer)
              )) +
    geom_point(size = 0.5, alpha = 0.25) +
    geom_smooth(method = "loess", span = 0.25, se = FALSE) +
    scale_x_log10() +
    labs(
      x = TeX("Diffusion coefficient"),
      y = "Avg. freq. (whole systems)",
      color = "Integer"
    ) +
  theme(
    legend.justification=c(1,1), 
    legend.position=c(0.20,1.00),
    legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
    legend.key.size = unit(0.15, "cm"),    # Decrease key size
    legend.text = element_text(size = 6), # Decrease text size
    legend.title = element_text(size = 6) # Decrease title size
  )
  
  p_main
}

############################
# MAIN SCRIPT
############################

ts <- load_processed_data()
params <- read_csv(PARAMS_CSV, show_col_types = FALSE)
p_main  <- create_main_plot(ts, params)

options(vsc.dev.args = list(width = 80, height = 70, res=300, units = "mm"))
print(p_main)

out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
ggsave(filename = out_file, plot = p_main, width = 80, height = 70, units = "mm", create.dir = TRUE)