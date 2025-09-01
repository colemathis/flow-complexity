################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(viridis)
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "30-mean-vs-kd"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- FALSE

DATA_DIR      			<<- "../../datasets/D04_kd=1e-6_1e1/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
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

load_and_process_time_series <- function() {

	timeseries_files <- list.files(DATA_DIR, pattern = TIMESERIES_FILES, recursive = TRUE, full.names = TRUE)

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
		filter(time == MAX_TIME) %>%
		filter(integer <= 10)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

	ts %>%
		group_by(diffusion_rate, integer) %>%
		summarize(
			mean_frequency = sum(frequency, na.rm = TRUE) / N_REACTORS,
			sd_frequency = sqrt(sum((frequency - mean_frequency)^2) / (N_REACTORS - 1)),
			.groups = "drop"
		)

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	ts <- ts %>%
		filter(integer %in% c(1,2,3,5,10))

	p <- ggplot(ts %>% filter(diffusion_rate > 1e-4),
				aes(
					x = diffusion_rate, 
					y = mean_frequency, 
					color = factor(integer)
				))

	p <- p + geom_point(size = 0.5, alpha = 0.25)
	p <- p + geom_line(stat="smooth", method = "loess", span = 0.35, se = FALSE, size = 1, alpha = 0.75)

	p <- p + coord_cartesian(xlim = c(1e-4, 1e1))

	p <- p + scale_x_log10(
		labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
	)
	# p <- p + scale_color_manual(values = c("darkblue", "blue", "yellow", "red"))

	p <- p + scale_color_viridis_d(
	option = "plasma",   # "magma" or "plasma" give a dark-to-bright progression
	begin = 0.1,        # start slightly lighter than pure black
	end = 0.7,          # stop before pure white
	direction = -1,      # 1 = dark→light
	name = "Integer"
	)

	p <- p + labs(
		x = TeX("Diffusion coefficient $k_d$"),
		y = "Average frequency",
		color = TeX("Integer"),
		# title = "Populations of 2’s, sweep over diffusion with multiple inflows",
		# caption = ID
	)

	p <- p + theme_minimal(base_size = 11)

	p <- p + theme(
		legend.justification=c(0,1), 
		legend.position=c(0.05,0.98),
		legend.background = element_rect(fill = "grey95", color = NA), # Optional: Customize legend background
		legend.key.size = unit(0.15, "cm"),    # Decrease key size
		legend.text = element_text(size = 6), # Decrease text size
		legend.title = element_text(size = 7), # Decrease title size
	)

	# show only one decimal in the legend labels
	# p <- p + scale_color_discrete(labels = function(x) {
	# 	x <- as.numeric(x)
	# 	x <- round(x, 1)
	# 	return(x)
	# })

	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	height	<- 60
	width	<- 80

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
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
	data   <- load_and_process_time_series()
}

# Plot the figure
p <- plot_figure(data)
saveRDS(p, file = file.path(CACHE_DIR, paste0(ID, ".rds")))