################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)
library(cowplot)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "31b-combine-30-31"
# USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SOURCE_FIG_1            <<- "cache/30-mean-vs-kd/30-mean-vs-kd.rds"
SOURCE_FIG_2            <<- "cache/31-std-vs-kd/31-std-vs-kd.rds"

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
    #     coord_cartesian(ylim = c(1e2, 1e6)) +
        # theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1)
		# ) +
		theme(legend.position = c(0.02, 0.98)) +
		theme(legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +
		theme(legend.key.size = unit(0.25, "lines")) +
		theme(legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25))
		

    inset_plot <- inset_plot +
        # labs(y = TeX("$\\Delta$Assembly")) +
        # scale_x_continuous(breaks = c(1, 1e5), labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))) +
        # scale_y_log10(breaks = c(1e1, 1e5), labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))) +
        theme(
			legend.position = "none",
            plot.background = element_rect(fill = "transparent", color = NA)
			) +
			theme(
				axis.title = element_text(color = "black"),
				axis.text = element_text(color = "black")
        )
		inset_plot <- inset_plot +
			theme(
			axis.title = element_text(size = 6),
			axis.text = element_text(size = 6)
			)
		inset_plot <- inset_plot +
			theme(
				panel.grid = element_line(color = "grey90"),
				panel.border = element_rect(color = "black", fill = NA, size = 0.25),
				# panel.grid.minor = element_blank(),
				# panel.grid.major.x = element_blank()
			)
		inset_plot <- inset_plot +
			labs(x = TeX("$k_d$"), y = "SD")
		inset_plot <- inset_plot +
			scale_x_log10(breaks = c(1e-3, 1e-1, 1e1), labels = c(TeX("$10^{-3}$"), "", TeX("$10^{1}$")))
		inset_plot <- inset_plot +
			scale_y_continuous(breaks = c(0, 20, 40, 60, 80), position = "right", labels = c("0", "", "", "", "80"))
		inset_plot <- inset_plot +
			theme(
				panel.background = element_rect(fill = "grey95", color = NA)
			)

	p <- ggdraw() +
		draw_plot(main_plot) +
		draw_plot(inset_plot, x = 0.31, y = 0.47, width = 0.37, height = 0.50)
		# theme(plot.background = element_rect(fill = "transparent", color = NA))

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