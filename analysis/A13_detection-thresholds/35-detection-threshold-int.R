################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "35-detection-threshold-int"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

DATA_DIR      			<<- "../../datasets/D04_kd=1e-6_1e1/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
CACHE_FILE       		<<- paste0(ID, ".csv")
FIGS_FILE        		<<- paste0(ID, ".pdf")

CACHE_PATH   			<<- file.path(CACHE_DIR, CACHE_FILE)
PARAMS_PATH  			<<- file.path(DATA_DIR, PARAMS_FILE)

################################################################################
# FUNCTIONS
################################################################################

setwd_to_script_path <- function() {

	# Case 1: Rscript
	args <- commandArgs(trailingOnly = FALSE)
	file_arg <- grep("^--file=", args, value = TRUE)
	if (length(file_arg) > 0) {
		fn = normalizePath(sub("^--file=", "", file_arg))
		setwd(dirname(fn))
		return(invisible())
	}

	# Case 2: VSCode
	if (!is.null(sys.frames()[[1]]$ofile)) {
		fn = normalizePath(sys.frames()[[1]]$ofile)
		setwd(dirname(fn))
		return(invisible())
	}

	# Error
	stop("Cannot determine script path.")

}

#==============================================================================#

load_params <- function() {

	read_csv(PARAMS_PATH, show_col_types = FALSE)

}

#==============================================================================#

load_and_process_time_series <- function() {

	timeseries_files <- list.files(DATA_DIR, pattern = TIMESERIES_FILES, recursive = TRUE, full.names = TRUE)

	ts_all <- map_dfr(timeseries_files, function(file) {
		if (file.info(file)$size == 0) return(NULL)
		ts <- read_csv(file, show_col_types = FALSE, progress = FALSE)
		process_data(ts)
	}, .progress = TRUE)

	# Save the processed data to a cache file
	dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
	write_csv(ts_all, CACHE_PATH)

	return(ts_all)

}

#==============================================================================#

process_data <- function(ts) {

	MAX_TIME   <- params$total_time[1]

	ts <- ts %>%
		filter(time == MAX_TIME)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

    # ts <- ts %>%
    #     group_by(diffusion_rate) %>%
    #     summarize(
    #         total_molecules = sum(frequency),
    #         total_mass = sum(integer * frequency),
    #         .groups = "drop"
    #     )

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
        summarise(highest_integer = max(integer), .groups = "drop")

    # For each diffusion rate, choose the highest of those integers
    detection_thresholds <- detectable %>%
        group_by(diffusion_rate) %>%
        summarise(max_detected_integer = max(highest_integer, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(max_detected_integer = ifelse(
            is.infinite(max_detected_integer), NA_integer_, max_detected_integer))

    # Return a dataframe with one row per diffusion_rate
    detection_thresholds
}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

# plot_figure <- function(ts) {
plot_figure <- function(detection_thresholds) {

  p <- ggplot(detection_thresholds,
              aes(
                  x = diffusion_rate, 
                  y = max_detected_integer, 
                  color = factor(log10(as.numeric(threshold)))
              )) +
    geom_point(size = 0.5, alpha = 0.25) +
	geom_line(stat = "smooth", method = "loess", span = 0.25, se = FALSE, size = 1.00, alpha = 0.75) +
    scale_x_log10(labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))) +
    scale_y_log10(labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))) +
	coord_cartesian(xlim = c(1e-4, 1e1)) +
    labs(
      x = TeX("Diffusion coefficient $k_d$"),
      y = "Highest integer",
      color = TeX("$log[X]$")
    ) +
	theme_minimal(base_size = 11) +
  theme(
    legend.justification=c(1,1), 
    legend.position=c(0.95,0.95),
    legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
    legend.key.size = unit(0.35, "cm"),    # Decrease key size
    # legend.text = element_text(size = 8), # Decrease text size
    # legend.title = element_text(size = 8) # Decrease title size
  )

	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	height	<- 60
	width	<- 80

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
		print(p)
	}

	if (SAVE_FIGS) {
		out_file <- file.path(FIGS_DIR, FIGS_FILE)
		ggsave(filename = out_file, plot = p, width = width, height = height, units = "mm", create.dir = TRUE)
	}

}

################################################################################
# MAIN SCRIPT
################################################################################

# Set-up environment
setwd_to_script_path()

# Load and process data
params <- load_params()

if (file.exists(CACHE_PATH) && USE_CACHE) {
	data   <- load_cached_data()
} else {
	data   <- load_and_process_time_series()
}

data <- data %>% filter(diffusion_rate > 1e-4)

detection_thresholds <- list(
    # `1e-7` = calculate_detection_threshold(ts, 1e-7),
    # `1e-6` = calculate_detection_threshold(ts, 1e-6),
    # `1e-5` = calculate_detection_threshold(ts, 1e-5),
    `1e-4` = calculate_detection_threshold(data, 1e-4),
    `1e-3` = calculate_detection_threshold(data, 1e-3),
    `1e-2` = calculate_detection_threshold(data, 1e-2),
    `1e-1` = calculate_detection_threshold(data, 1e-1)
) %>% 
  bind_rows(.id = "threshold")

# Plot the figure
p <- plot_figure(detection_thresholds)