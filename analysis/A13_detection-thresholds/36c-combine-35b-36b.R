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

ID         				<<- "36c-combine-35b-36b"
# USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SOURCE_FIG_1            <<- "cache/35b-detection-threshold-int-outflow/35b-detection-threshold-int-outflow.rds"
SOURCE_FIG_2            <<- "cache/36b-detection-threshold-AI-outflow/36b-detection-threshold-AI-outflow.rds"

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
        theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1)
		) +
		theme(legend.text = element_text(size = 6), legend.title = element_text(size = 7)) +
		theme(legend.key.size = unit(0.4, "lines"))

	# main_plot <- main_plot +
	# 	# theme(legend.position = "bottom") +
	# 	guides(color = guide_legend(ncol = 2))
		
	main_plot <- main_plot +
		coord_cartesian(ylim = c(1e0, 1e7))

	# # Ensure all geom_line layers in inset_plot use a thin line
	# main_plot$layers <- lapply(main_plot$layers, function(layer) {
	# if (inherits(layer$geom, "GeomLine")) {
	# 	# ggplot2 & cowplot 3.4+ prefer `linewidth`; older versions use `size`
	# 	# layer$aes_params$linewidth  <- 0.75
	# 	layer$aes_params$alpha  <- 0.5
	# 	# layer$geom_params$linewidth <- 0.25
	# }
	# layer
	# })

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
				panel.border = element_rect(color = "black", fill = NA, size = 0.25)
			)
		inset_plot <- inset_plot +
			theme(
				panel.grid = element_line(color = "grey90"),	
				panel.background = element_rect(fill = "grey95", color = NA)
			)
		inset_plot <- inset_plot +
			labs(x = TeX("$k_d$"), y = "AI")
		
		# Ensure all geom_line layers in inset_plot use a thin line
		inset_plot$layers <- lapply(inset_plot$layers, function(layer) {
		if (inherits(layer$geom, "GeomLine")) {
			# ggplot2 & cowplot 3.4+ prefer `linewidth`; older versions use `size`
			layer$aes_params$linewidth  <- 0.75
			# layer$aes_params$alpha  <- 0.5
			# layer$geom_params$linewidth <- 0.25
		}
		layer
		})

		# Ensure all geom_line and geom_point layers in inset_plot are thin
		inset_plot$layers <- lapply(inset_plot$layers, function(layer) {
		if (inherits(layer$geom, "GeomPoint")) {
			# Small points
			layer$aes_params$size       <- 0.1
			# layer$geom_params$size      <- 0.25
		}
		layer
		})

	p <- ggdraw() +
		draw_plot(main_plot) +
		draw_plot(inset_plot, x = 0.60, y = 0.47, width = 0.38, height = 0.50)
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