################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(cowplot)
library(patchwork)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "combine"
# USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SOURCE_FIG_1            <<- "cache/17-heatmap-skipped-constructive-rxn/17-heatmap-skipped-constructive-rxn.rds"
SOURCE_FIG_2            <<- "cache/18-heatmap-skipped-destructive-rxn/18-heatmap-skipped-destructive-rxn.rds"	
SOURCE_FIG_3            <<- "cache/19-heatmap-skipped-diffusion-rxn/19-heatmap-skipped-diffusion-rxn.rds"	
SOURCE_FIG_4            <<- "cache/20-heatmap-skipped-outflow-rxn/20-heatmap-skipped-outflow-rxn.rds"	

# SELECTED_SIMS 			<<- c(30, 58, 86)

# DATA_DIR      			<<- "../../_archives/milestones/23_distance-from-source-outflow-fixed/D_tmax=1e5/data"
# CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

# PARAMS_FILE 			<<- "params.csv"
# TIMESERIES_FILES		<<- "timeseries.csv"
# META_FILES				<<- "meta.csv"
# CACHE_FILE       		<<- paste0(ID, ".csv")
FIGS_FILE        		<<- paste0(ID, ".pdf")

# CACHE_PATH   			<<- file.path(CACHE_DIR, CACHE_FILE)
# PARAMS_PATH  			<<- file.path(DATA_DIR, PARAMS_FILE)
# ASSEMBLY_PATH           <<- "../../datasets/assembly-1e5.csv"

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

plot_figure <- function(ts) {

	p1 <<- readRDS(SOURCE_FIG_1)
	p2 <<- readRDS(SOURCE_FIG_2)
	p3 <<- readRDS(SOURCE_FIG_3)
	p4 <<- readRDS(SOURCE_FIG_4)

	p1 <- p1 +
	guides(fill = guide_colourbar(title = "Skipped\nreactions (%)\n")) +
	theme(
		legend.title      = element_text(size = 5, hjust = 0.5),
		legend.title.align = 0.5
	)
	# p1 <- p1 + theme(legend.position = "none")
	p2 <- p2 + theme(legend.position = "none")
	p3 <- p3 + theme(legend.position = "none")
	p4 <- p4 + theme(legend.position = "none")

	p1 <- p1 + theme(plot.title = element_text(size = 6, hjust = 0.5))
	p2 <- p2 + theme(plot.title = element_text(size = 6, hjust = 0.5))
	p3 <- p3 + theme(plot.title = element_text(size = 6, hjust = 0.5))
	p4 <- p4 + theme(plot.title = element_text(size = 6, hjust = 0.5))

	main_plot <- (
	((p1 + p2) / (p3 + p4)) |   # add a legend slot
	guide_area()                # (placeholder to collect guides)
	) +
	plot_layout(
		guides = "collect",
		widths = c(8, 1, 0.4)       # 5 % of the width for the legend
	) &
	theme(
		legend.text = element_text(size = 4),
		legend.key.height = unit(5, "mm")  # keep your existing tweak
	)

	theme(
		legend.title       = element_text(size = 5, hjust = 0.5),
		legend.title.align = 0.5
	)

	main_plot <- main_plot &
	theme(
		panel.spacing      = unit(0, "mm"),   # no gaps between panels
		legend.box.spacing = unit(0, "mm"),   # no gap before the legend
		legend.margin      = margin(0, 0, 0, 0),
		plot.margin        = margin(0, 0, 0, -20)  # kill outer white frame
	)

	p <- ggdraw() +
	draw_plot(main_plot)

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

	# Return the plot object
	return(p)

}

################################################################################
# MAIN SCRIPT
################################################################################

# Set-up environment
setwd_to_script_path()

# Load and process data
# params <- load_params()

# if (file.exists(CACHE_PATH) && USE_CACHE) {
# 	data   <- load_cached_data()
# } else {
# 	data   <- load_and_process_sim_data()
# }

# Plot the figure
p <- plot_figure(data)