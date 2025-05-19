################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "24-multipanel-pdf-inflow-outflow"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

DATA_DIR      			<<- "../../datasets/D03_kd=1e-2_1e2_randomized/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
META_FILES				<<- "meta.csv"
CACHE_FILE       		<<- paste0(ID, ".csv")
FIGS_FILE        		<<- paste0(ID, ".pdf")
FITS_FILE       		<<- paste0(ID, "_powerlaw_fits.csv")

CACHE_PATH   			<<- file.path(CACHE_DIR, CACHE_FILE)
PARAMS_PATH  			<<- file.path(DATA_DIR, PARAMS_FILE)
FITS_PATH       		<<- file.path(CACHE_DIR, FITS_FILE)

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

	MAX_TIME   <- params$total_time[1]
	N_REACTORS <- params$N_reactors[1]

	ts <- ts %>%
		filter(time == MAX_TIME)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

	MAX_INTEGER <<- 1000
	NBINS       <- 20

	ts <- ts %>%
		filter(integer <= MAX_INTEGER) %>%
		filter(chemostat_id %in% c(1, N_REACTORS))

	ts <- ts %>%
		complete(sim_number, chemostat_id, time, integer, fill = list(frequency = 0))

	# Log‑spaced bins
	breaks <- unique(round(exp(seq(log(1), log(MAX_INTEGER), length.out = NBINS + 1))))

	ts <- ts %>%
		mutate(
			bin       = cut(integer, breaks = breaks, include.lowest = TRUE, right = TRUE),
			bin_upper = breaks[-1][as.integer(bin)],
			bin_width = diff(breaks)[as.integer(bin)]
		) %>%
		group_by(sim_number, chemostat_id, bin) %>%
		summarise(
			freq_per_unit = sum(frequency) / first(bin_width),
			bin_upper     = first(bin_upper),
			.groups       = "drop"
		)

	

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#
# Compute power-law fits for each chemostat in each sim
compute_power_law_fits <- function(ts) {

	# Keep only finite, positive values that can be log‑transformed
	fits <- ts %>%
		filter(
			is.finite(freq_per_unit), is.finite(bin_upper),
			freq_per_unit > 0, bin_upper > 0
		) %>%
		group_by(sim_number, chemostat_id) %>%
		summarise(
			slope     = if (n() >= 2)
				coef(lm(log10(freq_per_unit) ~ log10(bin_upper)))[2]
				else NA_real_,
			intercept = if (n() >= 2)
				coef(lm(log10(freq_per_unit) ~ log10(bin_upper)))[1]
				else NA_real_,
			.groups   = "drop"
		)

	# Save the exponents for later reuse
	dir.create(dirname(FITS_PATH), recursive = TRUE, showWarnings = FALSE)
	write_csv(fits, FITS_PATH)

	return(fits)
}

plot_figure <- function(ts, fits) {

	N_SIM <- nrow(params)
	side <- ceiling(sqrt(N_SIM))

	p <- ggplot(ts, aes(x = bin_upper, y = freq_per_unit, color = factor(chemostat_id))) +
		geom_point(alpha = 0.3, size = 1) +
		scale_color_manual(
			values = c("1" = "darkgreen", "25" = "red"),
			breaks = c("1", "25"),
			labels = c("1" = "Chemostat 1", "25" = "Chemostat 25")
		) +
		facet_wrap(
			~ sim_number,
			ncol = side,
			scales = "fixed",
			labeller = as_labeller(function(sim) {
				kd <- params$diffusion_rate[match(sim, params$sim_number)]
				sprintf("%s: log(kd)=%.2f", sim, log10(kd))
			})
		) +
		labs(x = "Integer (bin upper boundary)", y = "Frequency") +
		scale_x_log10(
			breaks = 10^(0:ceiling(log10(MAX_INTEGER))),
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		) +
		scale_y_log10(
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		) +
		theme_bw(base_size = 8) +
		theme(
			legend.position = "none",
			panel.border = element_rect(color = "black", fill = NA),
			plot.title = element_text(hjust = 0.5)
		)

	# ---------------------------------------------------------------------------
	# Add power‑law fit lines (dashed) for each chemostat in each panel
	line_data <<- ts %>%
		group_by(sim_number, chemostat_id) %>%
		summarise(
			x_start = min(bin_upper),
			x_end   = max(bin_upper),
			.groups = "drop"
		) %>%
		left_join(fits, by = c("sim_number", "chemostat_id")) %>%
		mutate(
			y_start = 10^(intercept) * (x_start^slope),
			y_end   = 10^(intercept) * (x_end^slope)
		)

	p <- p +
		geom_abline(
			data = fits %>% filter(!is.na(slope), !is.na(intercept)),
			aes(intercept = intercept,
			    slope     = slope,
			    color     = factor(chemostat_id)),
			size     = 0.5,
			inherit.aes = FALSE
		)

	height <- 300
	width  <- 300

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

# Fit power‑law exponents and plot the figure
fits <- compute_power_law_fits(data)
p    <- plot_figure(data, fits)