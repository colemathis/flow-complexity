################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(patchwork)   # for combining ggplots
library(grid)        # for rasterGrob and gradient background

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "24c-plot-pdf-outflow-exponent-lattice-vs-randomized"
USE_CACHE  				<<- FALSE
PRINT_FIGS 				<<- FALSE
SAVE_FIGS 			 	<<- TRUE

FITS_PATH_LATTICE       <<- file.path("../A02_wrap-up-figs-D02/cache/24-multipanel-pdf-inflow-outflow/24-multipanel-pdf-inflow-outflow_powerlaw_fits.csv")
FITS_PATH_RANDOMIZED    <<- file.path("../A03_wrap-up-figs-D03/cache/24-multipanel-pdf-inflow-outflow/24-multipanel-pdf-inflow-outflow_powerlaw_fits.csv")

DATA_DIR      			<<- "../../datasets/D02_kd=1e-2_1e2_lattice/data"
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

	ts1 <- read_csv(FITS_PATH_LATTICE, show_col_types = FALSE) %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

	ts2 <- read_csv(FITS_PATH_RANDOMIZED, show_col_types = FALSE) %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

	# join ts1 and ts2 by setting column "topology" to "lattice" or "randomized"
	ts <- ts1 %>%
		mutate(topology = "lattice") %>%
		bind_rows(ts2 %>% mutate(topology = "randomized"))

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	#--- Prepare data ---------------------------------------------------------#
	ts <- ts %>%
		mutate(pos_slope = slope) %>%
		filter(chemostat_id == 25) %>%
		filter(diffusion_rate > -1)

	# Create a left‑to‑right grey→white gradient for the plot background
	gradient <- rasterGrob(
		matrix(colorRampPalette(c("grey60", "white"))(256), ncol = 256),
		width = unit(1, "npc"), height = unit(1, "npc"),
		interpolate = TRUE
	)

	# Shared colour scale
	colour_vals <- c("lattice" = "darkblue", "randomized" = "darkred")

	#--- Main plot (points only, no spline) -----------------------------------#
	p_abs <- ggplot(ts, aes(x = diffusion_rate, y = pos_slope, colour = topology)) +
		# Shade only up to diffusion_rate = 10 so that the plot is pure white beyond
		annotation_custom(gradient, xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf) +
		geom_point(alpha = 0.3, size = 0.5) +
		scale_x_log10(labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))) +
		labs(x = TeX("Diffusion coefficient $k_d$"),
				 y = TeX("Powel-Law Exponent $- \\alpha$"),
				 colour = "Topology") +
		coord_cartesian(ylim = c(-9, -1.0)) +
		scale_color_manual(values = colour_vals, name = "Topology") +
		theme_minimal(base_size = 11) +
		theme(legend.position = "none") +
		theme(panel.grid = element_blank()) +
		annotate("text", x = min(ts$diffusion_rate, na.rm = TRUE), y = -8.8, 
				 label = "Heterogeneous", hjust = 0, vjust = 1, size = 3.5) +
		annotate("text", x = max(ts$diffusion_rate, na.rm = TRUE), y = -8.8, 
				 label = "Well-mixed", hjust = 1, vjust = 1, size = 3.5)

	p_abs <- p_abs + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	#--- Inset: zoomed spline over x=20..100, y=-5..-1 -----------------------#
	p_inset <- ggplot(ts, aes(x = diffusion_rate, y = pos_slope, colour = topology)) +
		geom_point(alpha = 0.2, size = 0.3) +
		geom_line(data = ~ filter(.x, diffusion_rate >= 10 & diffusion_rate <= 100),
				  stat = "smooth", method = "loess", se = FALSE, span = 1.250,
				  linewidth = 1.00, alpha = 0.75) +
		scale_x_log10() +
		coord_cartesian(xlim = c(20, 100), ylim = c(-5, -1.75)) +
		scale_color_manual(values = colour_vals) +
		theme_minimal(base_size = 7) +
		theme(
			panel.grid   = element_blank(),
			axis.title   = element_blank(),
			panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
			plot.background = element_rect(fill = "white", colour = "black", linewidth = 1.0),
			plot.margin  = margin(5, 5, 5, 5),
			legend.position = "inside",
			legend.position.inside = c(0.02, 0.02),
			legend.justification = c("left", "bottom"),
			legend.background = element_blank(),
			legend.title = element_blank(),
			legend.key.size = unit(3, "mm"),
			legend.text = element_text(size = 9),
			legend.margin = margin(0, 0, 0, 0)
		)
	p_inset_grob <- ggplotGrob(p_inset)

	#--- Zoomed-region outline + connector lines + inset ----------------------#
	# Zoom rect in data coords matching the inset range
	# Inset grob placed at: x = log10(0.01)..log10(1.78), y = -8.5..-1.8
	p <- p_abs +
		annotate("rect", xmin = 20, xmax = 100, ymin = -6.75, ymax = -2,
				 fill = NA, colour = "grey40", linetype = "dashed", linewidth = 0.3) +
		# top-left corner of rect → top-left corner of inset
		annotate("segment", x = 20, xend = 0.01, y = -2, yend = -1.8,
				 colour = "grey40", linetype = "dashed", linewidth = 0.3) +
		# bottom-left corner of rect → bottom-left corner of inset
		annotate("segment", x = 20, xend = 1.5, y = -6.75, yend = -8.5,
				 colour = "grey40", linetype = "dashed", linewidth = 0.3) +
		annotation_custom(p_inset_grob,
						  xmin = log10(0.01), xmax = log10(1.78),
						  ymin = -8.5, ymax = -1.8)

	#--- Output ---------------------------------------------------------------#
	height <- 60   # mm
	width  <- 80    # mm

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