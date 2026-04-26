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

ID         				<<- "25-total-assembly-3sims"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

SELECTED_SIMS 			<<- c(30, 58, 86)

DATA_DIR      			<<- "../../_archives/milestones/23_distance-from-source-outflow-fixed/D_tmax=1e5/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
META_FILES				<<- "meta.csv"
CACHE_FILE       		<<- paste0(ID, ".csv")
FIGS_FILE        		<<- paste0(ID, ".pdf")

CACHE_PATH   			<<- file.path(CACHE_DIR, CACHE_FILE)
PARAMS_PATH  			<<- file.path(DATA_DIR, PARAMS_FILE)
ASSEMBLY_PATH           <<- "../../datasets/assembly-1e5.csv"

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

	assembly_indices <<- read.csv(ASSEMBLY_PATH)

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
		filter(sim_number %in% SELECTED_SIMS)

	# add assembly indices
	missing_value = 21
    ts <- ts %>%
        left_join(assembly_indices, by = "integer") %>%
        mutate(assemblyindex = ifelse(is.na(assemblyindex), missing_value, assemblyindex))

	# merge chemostats
	ts <- ts %>%
		group_by(sim_number, time, integer) %>%          # chemostat_id is absent here
		summarise(
		frequency     = sum(frequency, na.rm = TRUE),  # summed across chemostats
		assemblyindex = first(assemblyindex),          # keep its (unique) value
		.groups       = "drop"
		)

	# calculate total Assembly
    ts <- ts %>%
        group_by(sim_number, time) %>%
        filter(!is.na(frequency), !is.na(assemblyindex)) %>%
        mutate(
            weight       = (frequency - 1) / sum(frequency),
            contribution = exp(assemblyindex) * weight
        ) %>%
        summarise(
            total_assembly    = sum(contribution, na.rm = TRUE),
            max_assemblyindex = max(assemblyindex, na.rm = TRUE),
            .groups = "drop"
        )

	# add diffusion rate
	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	# Main plot: total assembly
	ts_main <- ts # %>% mutate(total_assembly = ifelse(time == 0, 10^1, total_assembly))

	log1p10_trans <- scales::trans_new(
		name      = "log1p10",
		transform = function(x) log10(1 + x),
		inverse   = function(x) 10^x - 1
	)

	main_plot <- ts_main %>%
		ggplot(aes(x = time, y = total_assembly, color = factor(diffusion_rate))) +
		geom_point(size = 0.5, alpha = 0.25) +
		geom_line(stat = "smooth", method = "loess", span = 0.10, se = FALSE, size = 0.6, alpha = 0.75) +
		# scale_y_log10(
		#	labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
		scale_y_continuous(
			trans   = log1p10_trans,
			breaks  = 10^(0:6),
			labels  = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
		) +
		scale_color_viridis_d(
			option = "turbo",
			begin = 0.6,
			end = 1.0,
			direction = -1,
			name = TeX("$\\log(k_d)$"),
			labels = function(x) sprintf("%.0f", log10(as.numeric(x)))
		) +
		labs(
			x = TeX("Time $t$"),
			y = TeX("Assembly $A$")
		) +
		scale_x_continuous(
			breaks = seq(0, 1e5, by = 2e4),
			labels = function(x) ifelse(x %in% c(0, 1e5), c("0", TeX("$10^5$")), "")
		) +
		coord_cartesian(ylim = c(1e2, 1e6)) +
		theme_minimal(base_size = 11) +
		theme(
			panel.grid.minor = element_blank(),
      		legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
			legend.title = element_text(size = 6),
			legend.text = element_text(size = 6),
			legend.direction = "horizontal"
		) +
        theme(legend.position = c(0.03, 0.98), legend.justification = c(0, 1))
	main_plot <- main_plot + guides(color = guide_legend(keyheight = unit(0.5, "lines"), keywidth = unit(1, "lines"), default.unit = "lines"))
	main_plot <- main_plot + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	# Inset plot: max assembly index over time
	inset_plot <- ts %>%
		ggplot(aes(x = time, y = max_assemblyindex, color = factor(diffusion_rate))) +
		geom_point(size = 0.25, alpha = 0.25) +
		geom_line(stat = "smooth", method = "loess", span = 0.05, se = FALSE, size = 0.3, alpha = 0.75) +
		#stat_summary_bin(fun = median, geom = "line", bins = 100, linewidth = 0.5, alpha = 0.75) +
		scale_color_viridis_d(
			option = "turbo",
			begin = 0.6,
			end = 1.0,
			direction = -1
		) +
		labs(
			# x = TeX("$t$"),
			x = "",
			y = TeX("Highest AI")
		) +
		scale_x_continuous(
			breaks = seq(0, 1e5, by = 2e4),
			labels = c("0", "", "", "", "", TeX("$10^5$"))
		) +
		coord_cartesian(ylim = c(0, 25)) +
		theme_bw(base_size = 8) +
		theme(
			panel.grid.minor = element_blank(),
			panel.grid = element_line(color = "grey90"),
			panel.background = element_rect(fill = "grey95", color = NA),
			legend.position = "none",
			plot.background = element_rect(fill = "transparent", color = NA)
		)

	p <- ggdraw() +
		draw_plot(main_plot) +
		draw_plot(inset_plot, x = 0.35, y = 0.20, width = 0.60, height = 0.48)

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

	return(p)

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
saveRDS(p, file = file.path(paste0(ID, ".rds")))

# ################################################################################
# # TEMPORARY: compare three smoothing methods
# ################################################################################

# make_panel <- function(ts, smooth_layer, title) {
# 	ts %>%
# 		ggplot(aes(x = time, y = total_assembly, color = factor(diffusion_rate))) +
# 		geom_point(size = 0.5, alpha = 0.25) +
# 		smooth_layer +
# 		scale_y_log10(
# 			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
# 		) +
# 		scale_color_viridis_d(option = "turbo", begin = 0.6, end = 1.0, direction = -1,
# 			name = TeX("$\\log(k_d)$"),
# 			labels = function(x) sprintf("%.0f", log10(as.numeric(x)))) +
# 		labs(x = TeX("Time $t$"), y = TeX("Assembly $A$"), title = title) +
# 		scale_x_continuous(
# 			breaks = seq(0, 1e5, by = 2e4),
# 			labels = function(x) ifelse(x %in% c(0, 1e5), c("0", TeX("$10^5$")), "")
# 		) +
# 		coord_cartesian(ylim = c(1e2, 1e6)) +
# 		theme_minimal(base_size = 9) +
# 		theme(
# 			panel.grid.minor  = element_blank(),
# 			panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
# 			legend.position   = "none",
# 			plot.title        = element_text(size = 8, hjust = 0.5)
# 		)
# }

# p1 <- make_panel(data,
# 	geom_line(stat = "smooth", method = "loess", span = 0.10, se = FALSE, linewidth = 0.75, alpha = 0.75),
# 	"LOESS")

# p2 <- make_panel(data,
# 	geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 0.75, alpha = 0.75),
# 	"GAM (cubic spline)")

# p3 <- make_panel(data,
# 	stat_summary_bin(fun = median, geom = "line", bins = 40, linewidth = 0.75, alpha = 0.75),
# 	"Running median (40 bins)")

# p_compare <- plot_grid(p1, p2, p3, nrow = 1)
# options(vsc.dev.args = list(width = 240, height = 70, res = 300, units = "mm"))
# print(p_compare)