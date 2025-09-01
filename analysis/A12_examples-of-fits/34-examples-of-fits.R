################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "34-examples-of-fits"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

MAX_INTEGER <<- 1000
NBINS       <- 20

DATA_DIR      			<<- "../../datasets/D02_kd=1e-2_1e2_lattice/data"
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

#==============================================================================#

plot_figure <- function(ts, fits) {

	SELECTED_SIMS <- c(1, 80)

	ts <- ts %>% filter(sim_number %in% SELECTED_SIMS)
	fits <- fits %>% filter(sim_number %in% SELECTED_SIMS)

	N_SIM <- nrow(params)
	side <- ceiling(sqrt(N_SIM))

	p <- ggplot(ts, aes(x = bin_upper, y = freq_per_unit, color = factor(chemostat_id))) +
		geom_point(alpha = 0.3, size = 1) +
		scale_color_manual(
			values = c("1" = "darkorange", "25" = "magenta"),
			breaks = c("1", "25"),
			labels = c("1" = "Inflow", "25" = "Outflow")
		) +
		facet_grid(
			rows = vars(sim_number),
			labeller = as_labeller(function(sim) {
				sim
			})
		) +
		geom_text(
			data = ts %>%
				distinct(sim_number) %>%
				mutate(
					x = Inf,
					y = Inf,
					label = sprintf("log[10](k[d]) == %.2f",
									 log10(params$diffusion_rate[match(sim_number, params$sim_number)]))
				),
			aes(x = x, y = y, label = label),
			inherit.aes = FALSE,
			hjust = 1.1, vjust = 1.5,
			size = 3,
			parse = TRUE
		) +
		labs(x = "Integers (binned)", y = "Frequency", color = NULL) +
		scale_x_log10(
			breaks = 10^(0:ceiling(log10(MAX_INTEGER))),
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		) +
		scale_y_log10(
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		) +
		theme_minimal(base_size = 11) +
		theme(
			legend.position = c(0.02, 0.54),
			legend.justification = c("left", "bottom"),
			legend.key.size = unit(0.5, "lines"),
			legend.text = element_text(size = 9),
			legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
			panel.border = element_rect(color = "black", fill = NA),
			plot.title = element_text(hjust = 0.5),
			strip.text.y = element_blank()  # Hide facet labels on the side
		)

	# Hide legend for all but the first panel (top panel)
	# This requires patchwork or cowplot for true per-panel legends, but as a workaround:
	# Only show legend if sim_number is the first in SELECTED_SIMS
	if (!is.null(SELECTED_SIMS) && length(SELECTED_SIMS) > 1) {
		p <- p + guides(color = guide_legend(override.aes = list(alpha = 1)))
	}

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
				size     = 1.00,
				alpha    = 0.75,
				inherit.aes = FALSE
				# show.legend = FALSE
			)

	# Add annotation for chemostat_id == 25
	p <- p +
		geom_text(
			data = fits %>% filter(chemostat_id == 25, !is.na(slope), !is.na(intercept)),
			aes(
				x = Inf, y = Inf,
				# label = sprintf("f(x) == %.2f * x^{-%.2f}", 10^intercept, abs(slope))
				label = paste0(TeX(sprintf("$\\alpha = %.2f$", abs(slope))))
			),
			color = "magenta",
			inherit.aes = FALSE,
			hjust = 1.2, vjust = 4.0,
			size = 3,
			parse = TRUE
		)


	height <- 60
	width  <- 80

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