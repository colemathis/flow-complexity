################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "29-population-histogram"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

DATA_DIR      			<<- "../../datasets/D01_inflow=1e3_1e4_kd=1e-1_1e3/data"
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
	# N_REACTORS <- params$N_reactors[1]

	ts <- ts %>%
		filter(time == MAX_TIME)

	missing_value = 21
    ts <- ts %>%
        left_join(assembly_indices, by = "integer") %>%
        mutate(assemblyindex = ifelse(is.na(assemblyindex), missing_value, assemblyindex))

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	ts <- ts %>%
		group_by(sim_number, integer, assemblyindex) %>%
		summarise(frequency = sum(frequency), .groups = "drop")

	ts_long <- ts %>%
		filter(sim_number == 90, integer < 5000) %>%
		pivot_longer(cols = c(integer, assemblyindex), names_to = "type", values_to = "xval")

	ts_long <- ts_long %>%
		filter(type == "integer")

	# print(ts_long, n = 30)

	# 2) Log-spaced breaks
	min_pos <- min(ts_long$xval[ts_long$xval > 0])
	max_val <- max(ts_long$xval)
	breaks  <- 10^seq(floor(log10(min_pos)), ceiling(log10(max_val)), length.out = 25)

	# 3) Log-binned histogram with width normalization
	hist_data <- ts_long %>%
	filter(xval > 0) %>%
	mutate(bin = cut(xval, breaks = breaks, include_lowest = TRUE, right = TRUE)) %>%
	mutate(bin_idx = as.integer(bin),
			bin_min = breaks[bin_idx],
			bin_max = breaks[bin_idx + 1]) %>%
	group_by(type, bin_min, bin_max) %>%
	summarise(count = sum(frequency), .groups = "drop") %>%
	mutate(width = bin_max - bin_min,
			density = count / width)   # normalized by bin width (per unit x)

	# ts_long <- ts_long %>%
	# 	filter(xval >= 0, xval <= 3000)

	# ts_long <- ts_long %>%
	# 	group_by(type, xval) %>%
	# 	summarise(frequency = sum(frequency), .groups = "drop")

	# print(ts_long,n=30)
	# print(hist_data, n = 30)

	width_factor <- 0.8  # 1.0 = full bin width, <1 narrows the bar

	hist_data_narrow <- hist_data %>%
	mutate(
		bin_mid = 10^((log10(bin_min) + log10(bin_max)) / 2),
		log_half_width = (log10(bin_max) - log10(bin_min)) / 2,
		bin_min_narrow = 10^(log10(bin_mid) - width_factor * log_half_width),
		bin_max_narrow = 10^(log10(bin_mid) + width_factor * log_half_width)
	)

	# p <- ggplot(ts_long, aes(x = xval, y = frequency, fill = type)) +
	# p <- ggplot(hist_data, aes(xmin = bin_min, xmax = bin_max, ymin = 0, ymax = density, fill = type)) +
	p <- ggplot(hist_data_narrow,
		aes(xmin = bin_min_narrow, xmax = bin_max_narrow,
			ymin = 0, ymax = density, fill = type)) +
		# geom_col(
		# 	data = filter(ts_long, type == "integer"),
		# 	position = "identity", alpha = 0.5, width = 0.8
		# ) +
		# geom_line() +
		geom_rect(alpha = 0.7) +
		# geom_col() + 
		# geom_col(
		# 	data = filter(ts_long, type == "assemblyindex"),
		# 	position = "identity", alpha = 0.5, width = 0.8
		# ) +
		labs(
			x = TeX("Integer $z$"),
			y = "Density",
			# title = "Histogram of Frequency vs Integer and Assembly Index",
			# fill = "Type",
			# caption = ID
		) +
		scale_fill_manual(
			values = c("integer" = "#1f77b4", "assemblyindex" = "#ff7f0e"),
			labels = c("assemblyindex" = "Assembly Index", "integer" = "Integer")
		) +
		# scale_x_log10(
		# 	labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
		# ) +
		scale_x_log10(
		breaks = c(1, 10, 100, 1000),
		labels = scales::trans_format(
			"log10",
			function(x) TeX(sprintf("$10^{%g}$", x))
		)
		) +
		scale_y_log10(
			labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
		) +
		theme_minimal(base_size = 11) +
		# theme(
		# 	legend.position = c(0.85, 0.98),
		# 	legend.justification = c("right", "top"),
		# 	# legend.background = element_rect(fill = "grey95", color = NA),
		# 	# legend.text = element_text(size = 6) # Make legend text smaller
		# ) +
		# hide the legend
		theme(legend.position = "none") +
		guides(fill = guide_legend(title = NULL))

	# Add vertical dashed line at y = ymax with rotated label "(truncated)"
	ymax <- max(ts_long$xval, na.rm = TRUE)
	# p <- p +
	# 	geom_vline(xintercept = ymax, linetype = "dashed", color = "grey55") +
	# 	annotate("text", x = ymax, y = max(ts_long$frequency, na.rm = TRUE), 
	# 			 label = "(truncated)", angle = 90, vjust = -0.5, hjust = 1, size = 3)

	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	p <- p +
			theme(legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25))

	# p <- ggplot(ts,
	# 			aes(
	# 				x = diffusion_rate, 
	# 				y = mean_frequency, 
	# 				color = factor(log10(inflow_mols))
	# 			))

	# p <- p + geom_point(size = 0.5, alpha = 0.25)
	# p <- p + geom_smooth(method = "loess", span = 0.5, se = FALSE, size = 0.5)

	# p <- p + scale_x_log10()
	# p <- p + scale_color_manual(values = c("darkblue", "blue", "yellow", "red"))

	# p <- p + labs(
	# 	x = TeX("Diffusion coefficient"),
	# 	y = "Avg. freq. of 2’s (whole sys.)",
	# 	color = TeX("Log(I)"),
	# 	title = "Populations of 2’s, sweep over diffusion with multiple inflows",
	# 	caption = ID
	# )

	# p <- p + theme(
	# 	legend.justification=c(0,0), 
	# 	legend.position=c(0.05,0.05),
	# 	legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
	# 	legend.key.size = unit(0.15, "cm"),    # Decrease key size
	# 	legend.text = element_text(size = 6), # Decrease text size
	# 	legend.title = element_text(size = 6), # Decrease title size
	# 	plot.title = element_text(size = 6),
	# 	plot.caption = element_text(size = 6, color = "grey50")
	# )

	# # show only one decimal in the legend labels
	# p <- p + scale_color_discrete(labels = function(x) {
	# 	x <- as.numeric(x)
	# 	x <- round(x, 1)
	# 	return(x)
	# })

	height	<- 70
	width	<- 80

	if (PRINT_FIGS) {
		options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
		print(p)
	}

	if (SAVE_FIGS) {
		out_file <- file.path(FIGS_DIR, FIGS_FILE)
		ggsave(filename = out_file, plot = p, width = width, height = height, units = "mm", create.dir = TRUE)
	}

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