################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "27-combine-25-with-26-inset"
# USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SOURCE_FIG_1            <<- "25-total-assembly-3sims.rds"
SOURCE_FIG_2            <<- "26-total-assembly-inflow-outflow.rds"

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

    main_plot <- readRDS(SOURCE_FIG_1)
    inset_plot <- readRDS(SOURCE_FIG_2)

    main_plot <- main_plot + 
        coord_cartesian(ylim = c(1e2, 1e6)) +
        theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.direction = "horizontal")

    inset_plot <- inset_plot +
        labs(y = TeX("$\\Delta$Assembly")) +
        scale_x_continuous(breaks = c(1, 1e5), labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))) +
        scale_y_log10(breaks = c(1e1, 1e5), labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))) +
        theme(
            plot.background = element_rect(fill = "transparent", color = NA)
        )

	p <- ggdraw() +
		draw_plot(main_plot) +
		draw_plot(inset_plot, x = 0.25, y = 0.15, width = 0.7, height = 0.45)

	height <- 80
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