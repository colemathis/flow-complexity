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

TITLE        <- "# molecules and mass across sweep over diffusion"
ID           <- "molecules-and-mass-vs-kd"
TIME        <- 1e5
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
            filter(time == TIME) %>%
            collect()

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}

calculate_number_of_molecules <- function(ts) {
    ts %>%
        group_by(sim_number) %>%
        summarize(
            total_molecules = sum(frequency),
            total_mass = sum(integer * frequency),
            .groups = "drop"
        )
}

create_main_plot <- function(ts, params) {
  # Attach diffusion_rate from params to each row of the time‑series data
  ts <- ts %>%
    left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")
  
  scaling = 4e3

  p_main <- ggplot(ts) +
    geom_point(aes(
        x = diffusion_rate, 
        y = total_molecules, 
        color = "Molecules"
    ), size = 0.5, alpha = 0.25) +
    geom_point(aes(
        x = diffusion_rate, 
        y = total_mass/scaling, 
        color = "Mass"
    ), size = 0.5, alpha = 0.25) +
    scale_y_continuous(
        name = "Molecules",
        sec.axis = sec_axis(~ . * scaling, name = "Mass")
    ) +
    scale_color_manual(
        name = "Legend",
        values = c("Molecules" = "blue", "Mass" = "red")
    ) +
    scale_x_log10() +

    labs(
      x = TeX("Diffusion coefficient")
    ) +
    theme(
      legend.justification = c(0, 0), 
      legend.position = c(0.05, 0.05),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key.size = unit(0.15, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_blank(),
      axis.title.y.right = element_text(color = "red"), # Optional: Customize right axis label color
      axis.title.y.left = element_text(color = "blue")  # Optional: Customize left axis label color
    )
  
  p_main
}

############################
# MAIN SCRIPT
############################

ts <- load_processed_data()
params <- read_csv(PARAMS_CSV, show_col_types = FALSE)
ts <- calculate_number_of_molecules(ts)
p_main  <- create_main_plot(ts, params)

width <- 120
height <- 70

options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
print(p_main)

out_file <- file.path(FIGS_DIR, sprintf("%s.pdf", FIGS_FILE))
ggsave(filename = out_file, plot = p_main, width = width, height = height, units = "mm", create.dir = TRUE)