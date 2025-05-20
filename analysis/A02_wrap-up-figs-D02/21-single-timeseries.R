################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "21-single-timeseries"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SELECTED_SIM            <<- 100

DATA_DIR      			<<- "../../datasets/D02_kd=1e-2_1e2_lattice/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
META_FILES				<<- "meta.csv"
CACHE_FILE       		<<- paste0(ID, "_", SELECTED_SIM, ".csv")
FIGS_FILE        		<<- paste0(ID, "_", SELECTED_SIM, ".pdf")

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

load_and_process_sim_data <- function() {

	timeseries_files <- list.files(DATA_DIR, pattern = TIMESERIES_FILES, recursive = TRUE, full.names = TRUE)
	meta_files       <- list.files(DATA_DIR, pattern = META_FILES, recursive = TRUE, full.names = TRUE)

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

	ts <- ts %>%
		filter(sim_number == SELECTED_SIM) %>%
		filter(integer <= 10)

	# complete the data with zeros where appropriate
	ts <- ts %>%
		complete(sim_number, chemostat_id, time, integer, fill = list(frequency = 0))

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	N_REACTORS <- params$N_reactors[SELECTED_SIM]
	side <- ceiling(sqrt(N_REACTORS))

	p <- ggplot(ts, aes(x = time, y = frequency, color = factor(integer))) +
		geom_line(linewidth = 0.2) +
		facet_wrap(
			~ chemostat_id,
			ncol = side,
			scales = "fixed",
			labeller = labeller(
				chemostat_id = function(x) paste("chemostat #", x)
			)
		) +
		labs(x = "Time", y = "Frequency") +
		scale_y_log10(
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		) +
		theme_bw(base_size = 8) +
		theme(
			legend.position = "none",
			panel.border = element_rect(color = "black", fill = NA),
			plot.title = element_text(hjust = 0.5)
		) +
		ggtitle(
			TeX(sprintf(
				"%s Simulation %d: $log_{10}(I)=%.2f$, $log_{10}(k_d)=%.2f$",
				"First ten populations.",
				SELECTED_SIM,
				log10(params$inflow_mols[SELECTED_SIM]),
				log10(params$diffusion_rate[SELECTED_SIM])
			))
		)

	height <- 200
	width  <- 200

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res = 300, units = "mm"))
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
	data   <- load_and_process_sim_data()
}

# Plot the figure
plot_figure(data)