################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(patchwork)   # for combining ggplots

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "24b-plot-pdf-inflow-outflow-exponent"
USE_CACHE  				<<- FALSE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

FITS_PATH               <<- file.path("cache/24-multipanel-pdf-inflow-outflow/24-multipanel-pdf-inflow-outflow_powerlaw_fits.csv")

DATA_DIR      			<<- "../../datasets/D03_kd=1e-2_1e2_randomized/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
META_FILES				<<- "meta.csv"
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

load_and_process_sim_data <- function() {

	ts <- read_csv(FITS_PATH, show_col_types = FALSE)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

}

#==============================================================================#

# process_data <- function(ts) {

# 	MAX_TIME   <- params$total_time[1]
# 	N_REACTORS <- params$N_reactors[1]

# 	ts <- ts %>%
# 		filter(time == MAX_TIME) %>%
# 		filter(integer == 2)

# 	ts <- ts %>%
# 		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
# 		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

# 	ts %>%
# 		group_by(diffusion_rate, inflow_mols, integer) %>%
# 		summarize(
# 			mean_frequency = sum(frequency, na.rm = TRUE) / N_REACTORS,
# 			sd_frequency = sqrt(sum((frequency - mean_frequency)^2) / (N_REACTORS - 1)),
# 			.groups = "drop"
# 		)

# }

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	#--- Prepare data ---------------------------------------------------------#
	ts <- ts %>%
		mutate(neg_slope = -slope)

	# Panel A: Absolute slope by chemostat
	ts <- ts %>%
		mutate(chemostat_label = factor(chemostat_id, levels = c(1, 25), labels = c("inflow", "outflow")))

	p_abs <- ggplot(ts, aes(x = diffusion_rate, y = neg_slope, colour = chemostat_label)) +
		geom_point(alpha = 0.3, size = 1.3) +
		geom_smooth(method = "loess", se = FALSE, span = 0.5) +
		scale_x_log10(labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))) +
		coord_cartesian(ylim = c(0, 8)) +
		scale_colour_manual(values = c("inflow" = "darkgreen", "outflow" = "red"),
											 labels = c("inflow" = "1: inflow", "outflow" = "25: outflow")) +
		labs(x = TeX("$k_d$"),
				 y = TeX("$k$"),
				 colour = "Chemostats") +
		theme_bw() +
		theme(
			legend.position = c(0.02, 0.98),
			legend.justification = c("left", "top"),
			legend.background = element_rect(fill = alpha("white", 0.9), colour = "black"),
			legend.title = element_text(face = "bold")
		)

	# Panel B: Difference in absolute slope between chemostat 1 and 25
	diff_ts <- ts %>%
		filter(chemostat_id %in% c(1, 25)) %>%
		select(sim_number, diffusion_rate, chemostat_id, neg_slope) %>%
		pivot_wider(names_from = chemostat_id, values_from = neg_slope, names_prefix = "ch") %>%
		mutate(delta_abs = ch1 - ch25)

	p_diff <- ggplot(diff_ts, aes(x = diffusion_rate, y = delta_abs)) +
		geom_hline(yintercept = 0, linetype = "dashed") +
		geom_point(size = 1.3, alpha = 0.3, colour = "blue") +
		geom_smooth(method = "loess", se = FALSE, span = 0.5, colour = "blue") +
		scale_x_log10(labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))) +
		coord_cartesian(ylim = c(0, 2.5)) +
		labs(x = TeX("$k_d$"),
				 y = TeX("$\\Delta k$")) +
		annotate("text",
						 x = min(diff_ts$diffusion_rate, na.rm = TRUE),
						 y = max(diff_ts$delta_abs, na.rm = TRUE),
						 label = TeX("$\\Delta k = k_{inflow} - k_{outflow}$"),
						 hjust = 0, vjust = 1, size = 4) +
		theme_bw()

	#--- Combine panels -------------------------------------------------------#
	p <- p_abs / p_diff + plot_layout(heights = c(2, 1))

	#--- Output ---------------------------------------------------------------#
	height <- 120   # mm
	width  <- 120    # mm

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width,
									height = height,
									res = 300,
									units = "mm"))
		print(p)
	}

	if (SAVE_FIGS) {
		out_file <- file.path(FIGS_DIR, FIGS_FILE)
		ggsave(filename = out_file,
				plot = p,
				width = width,
				height = height,
				units = "mm",
				create.dir = TRUE)
	}

	invisible(p)
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
p <- plot_figure(data)