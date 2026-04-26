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
        filter(!is.na(frequency), !is.na(assemblyindex), assemblyindex > 0) %>%
        mutate(
            weight       = (frequency - 1) / sum(frequency),
            contribution = log(assemblyindex) * weight
        ) %>%
        summarise(
            total_assembly = sum(contribution, na.rm = TRUE),
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

	# log1p10_trans <- scales::trans_new(
	# 	name      = "log1p10",
	# 	transform = function(x) log10(1 + x),
	# 	inverse   = function(x) 10^x - 1
	# )

	# prepend t=0, total_assembly=0 for each sim
	ts <- bind_rows(
		ts %>% distinct(diffusion_rate) %>% mutate(time = 0, total_assembly = 10^(-1)),
		ts
	) %>% arrange(diffusion_rate, time)

	p <- ts %>%
		ggplot(aes(x = time, y = total_assembly, color = factor(diffusion_rate))) +
		geom_point(size = 0.5, alpha = 0.25) +
		geom_line(stat = "smooth", method = "loess", span = 0.05, se = FALSE, size = 0.6, alpha = 0.75) +
		scale_y_log10(
			# trans   = log1p10_trans,
			# breaks = c(1, 2),
			# labels = c(TeX("$10^0$"), TeX("$0.2 \\times 10^1$"))
		) +
		# scale_color_discrete(
		scale_color_viridis_d(
			option = "turbo",   # "magma" or "plasma" give a dark-to-bright progression
			begin = 0.6,        # start slightly lighter than pure black
			end = 1.0,          # stop before pure white
			direction = -1,      # 1 = dark→light
			name = TeX("$\\log(k_d)$"),
			labels = function(x) sprintf("%.0f", log10(as.numeric(x)))
		) +
		labs(
			x = TeX("$t$"),
			y = TeX("$A$")
		) +
		scale_x_continuous(
			breaks = seq(0, 1e5, by = 2e4),
			labels = function(x) ifelse(x %in% c(0, 1e5), c("0", TeX("$10^5$")), "")
		) +
		coord_cartesian(ylim = c(9e-1, 2)) +
		theme_minimal(base_size = 11) +
		theme(
			panel.grid.minor = element_blank(),
			# panel.grid.major.x = element_blank(),
			# legend.position = c(0.85, 0.3),
			# legend.background = element_rect(fill = "white", color = "black"),
      		legend.background = element_rect(fill = "grey95", color = NA),	
			legend.title = element_text(size = 6),
			legend.text = element_text(size = 6),
			legend.direction = "horizontal"
		) +
		# theme(
		# 	panel.grid.major = element_blank(),
		# 	panel.grid.minor = element_blank()
		# )
        theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
		p <- p + guides(color = guide_legend(keyheight = unit(0.5, "lines"), keywidth = unit(1, "lines"), default.unit = "lines"))

	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

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