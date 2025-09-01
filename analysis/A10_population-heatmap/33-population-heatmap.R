################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "33-population-heatmap"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

DATA_DIR      			<<- "../../datasets/D04_kd=1e-6_1e1/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
CACHE_FILE       		<<- paste0(ID)
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
	# write_csv(ts_all, CACHE_PATH)
	write_rds(ts_all, CACHE_PATH)

	return(ts_all)

}

#==============================================================================#

process_data <- function(ts) {

	MAX_TIME   <- params$total_time[1]
	N_REACTORS <- params$N_reactors[1]

	ts <- ts %>%
		filter(time == MAX_TIME, integer > 1)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

	# Bin integer and aggregate frequency for heatmap
	bin_edges <- unique(round(exp(seq(log(1), log(1e7), length.out = 51))))
	ts <- ts %>%
		mutate(integer_bin = cut(
		integer,
		breaks = bin_edges,
		include.lowest = TRUE,
		right = TRUE
		)) %>%
		group_by(diffusion_rate, integer_bin) %>%
		summarise(frequency = sum(frequency), .groups = "drop")

	# Complete the dataset: for non-existent frequencies assign 0
	ts <- ts %>%
		complete(diffusion_rate, integer_bin, fill = list(frequency = 0))

}

#==============================================================================#

load_cached_data <- function() {

	# read_csv(CACHE_PATH, show_col_types = FALSE)
	read_rds(CACHE_PATH)

}

#==============================================================================#

plot_figure <- function(ts) {

  p <- ggplot(ts, aes(x = diffusion_rate, y = integer_bin, fill = frequency)) +
    geom_tile() +
    scale_x_log10(
		labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
	) +
    scale_fill_viridis_c(name = "Freq.", na.value = "grey") +
    labs(
      x = TeX("Diffusion coefficient $k_d$"),
      y = "Integer",
    ) +
    scale_y_discrete(
      breaks = levels(ts$integer_bin)[seq(1, length(levels(ts$integer_bin)), length.out = 8)],
      labels = c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7))
    ) +
	coord_cartesian(xlim = c(1e-4, 1e1)) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = c(0.95, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
      legend.key.size = unit(0.5, "lines"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
	panel.grid = element_blank()
    )

#   if (OUTFLOW_ONLY) {
#     p_main <- p_main +
#       labs(title = paste(TITLE, "(outflow only)"))
#   } else {
#     p_main <- p_main +
#       labs(title = paste(TITLE, "(whole system)"))
#   }

  # add a white line with slope ~ -1
#   p <- p +
#     geom_segment(
#       aes(x = 1e-6, xend = 1e1, y = levels(ts$integer_bin)[length(levels(ts$integer_bin))], yend = levels(ts$integer_bin)[1]),
#       color = "white",
#       inherit.aes = FALSE
#   )


	height	<- 50
	width	<- 80

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
		print(p)
	}

	if (SAVE_FIGS) {
		out_file <- file.path(FIGS_DIR, FIGS_FILE)
		ggsave(filename = out_file, plot = p, width = width, height = height, units = "mm", create.dir = TRUE)
	}

	return(p)

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

# Plot the figure
p <- plot_figure(data)
saveRDS(p, file = file.path(CACHE_DIR, paste0(ID, ".rds")))