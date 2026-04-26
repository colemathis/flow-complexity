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

ID         				<<- "heatmaps"
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SELECTED_SIMS 			<<- c(30, 58, 86)

CACHE_PATH   			<<- "cache/33-distance-from-source/33-distance-from-source.csv"
FIGS_DIR      			<<- "figs"
FIGS_FILE        		<<- paste0(ID, ".pdf")

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

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	# Fix diffusion_rate labels
	ts <- ts %>%
		mutate(diffusion_rate = case_when(
			sim_number == 30 ~ 1e-4,
			sim_number == 58 ~ 1e-2,
			sim_number == 86 ~ 1e0,
			TRUE ~ diffusion_rate
		))

	# Derive grid coordinates from chemostat_id
	ts <- ts %>%
		mutate(
			col_x = (chemostat_id - 1) %% 5 + 1,
			row_y = (chemostat_id - 1) %/% 5 + 1
		)

	# Facet label
	ts <- ts %>%
		mutate(kd_label = factor(
			log10(diffusion_rate),
			levels = c(-4, -2, 0),
			labels = c(
				"log[10](k[d])==-4",
				"log[10](k[d])==-2",
				"log[10](k[d])==0"
			)
		))

	highlight_blue <- ts %>% filter(sim_number %in% c(30, 58), chemostat_id %in% c(5, 21))
	highlight_red  <- ts %>% filter(sim_number %in% c(30, 58), chemostat_id %in% c(9, 13, 17))

	p <- ggplot(ts, aes(x = col_x, y = row_y, fill = mean_ai)) +
		geom_tile(color = "white", linewidth = 0.5) +
		geom_tile(data = highlight_blue, fill = NA, color = "blue", linewidth = 1.2, width = 0.90, height = 0.90) +
		geom_tile(data = highlight_red,  fill = NA, color = "red",  linewidth = 1.2) +
		geom_text(aes(label = chemostat_id), color = "white", size = 4) +
		facet_wrap(~ kd_label, nrow = 1, labeller = label_parsed) +
		scale_y_reverse(expand = c(0, 0)) +
		scale_x_continuous(expand = c(0, 0)) +
		scale_fill_viridis_c(
			option = "turbo",
			name = "Mean AI"
		) +
		coord_fixed() +
		labs(x = NULL, y = NULL) +
		theme_minimal(base_size = 14) +
		theme(
			panel.border = element_blank(),
			legend.position = "left",
			legend.title = element_text(size = 12),
			legend.text = element_text(size = 11),
			strip.text = element_text(size = 13),
			axis.text = element_blank(),
			axis.ticks = element_blank()
		)

	height <- 60
	width  <- 160

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res = 300, units = "mm"))
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

# Load data from existing cache
data <- load_cached_data()

# Plot the figure
p <- plot_figure(data)
