################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "32-molecules-and-mass"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

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
		filter(time == MAX_TIME)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

    ts <- ts %>%
        group_by(diffusion_rate) %>%
        summarize(
            total_molecules = sum(frequency),
            total_mass = sum(integer * frequency),
            .groups = "drop"
        )

}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

scaling = 4e3

  ts <- ts %>% filter(diffusion_rate > 1e-4)

  p <- ggplot(ts) +
    geom_point(aes(
        x = diffusion_rate, 
        y = total_molecules, 
        color = "Molecules"
    ), size = 0.5, alpha = 0.25) +
	geom_line(aes(
		x = diffusion_rate, 
		y = total_molecules, 
		color = "Molecules"
	), stat="smooth", method = "loess", span = 0.35, se = FALSE, size = 0.75, alpha = 0.75) +
    geom_point(aes(
        x = diffusion_rate, 
        y = total_mass/scaling, 
        color = "Mass"
    ), size = 0.5, alpha = 0.25) +
	geom_line(aes(
		x = diffusion_rate,
		y = total_mass/scaling,
		color = "Mass"
	), stat="smooth", method = "loess", span = 0.35, se = FALSE, size = 0.75, alpha = 0.75) +
    scale_y_continuous(
        name = TeX("Molecules$\\times 10^{3}$"),
		# breaks = c(0, 25000),
		breaks = seq(0, 25000, by = 5000),
		# labels = c("0", "25000"),
		labels = c("0", "5", "10", "15", "20", "25"),
		sec.axis = sec_axis(
			~ . * scaling,
			name = TeX("Mass$\\times 10^{7}$"),
			# breaks = c(0, 1e8),
			breaks = seq(0, 1e8, by = 2.5e7),
			labels = c("0", "2.5", "5.0", "7.5", "10.0"),
			# labels = c("0", TeX("$10^{8}$"))
		),
    ) +
    scale_color_manual(
        name = "Legend",
        values = c("Molecules" = "#3366CC", "Mass" = "#DC3912")
    ) +
    scale_x_log10(
		labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
	) +
    labs(
      x = TeX("Diffusion coefficient $k_d$"),
    ) +
	coord_cartesian(xlim = c(1e-4, 1e1)) +
	theme_minimal(base_size = 11) +
    theme(
      legend.justification = c(0, 0), 
      legend.position = c(0.02, 0.02),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25), # Optional: Customize legend background
      legend.key.size = unit(0.25, "cm"),
      legend.text = element_text(size = 5),
      legend.title = element_blank(),
	# axis.title.y.left = element_text(color = "#3366CC"), # Match "Molecules" curve color (default ggplot2 blue)
	# axis.title.y.right = element_text(color = "#DC3912"), # Match "Mass" curve color (default ggplot2 red)
      axis.title.y.right = element_text(color = "#DC3912"), # Optional: Customize right axis label color
      axis.title.y.left = element_text(color = "#3366CC"),  # Optional: Customize left axis label color
      axis.text.y = element_text(color = "#3366CC"),      # Left‑axis tick labels
      axis.text.y.right = element_text(color = "#DC3912"), # Right‑axis tick labels
    )

	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))

	# p <- ggplot(ts %>% filter(diffusion_rate > 1e-4),
	# 			aes(
	# 				x = diffusion_rate, 
	# 				y = mean_frequency, 
	# 				color = factor(integer)
	# 			))

	# p <- p + geom_point(size = 0.5, alpha = 0.25)
	# p <- p + geom_line(stat="smooth", method = "loess", span = 0.35, se = FALSE, size = 1, alpha = 0.75)

	# p <- p + coord_cartesian(xlim = c(1e-4, 1e1))

	# p <- p + scale_x_log10(
	# 	labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%f}$", x)))
	# )
	# # p <- p + scale_color_manual(values = c("darkblue", "blue", "yellow", "red"))

	# p <- p + labs(
	# 	x = TeX("Diffusion coefficient $k_d$"),
	# 	y = "Average frequency",
	# 	color = TeX("Integer"),
	# 	# title = "Populations of 2’s, sweep over diffusion with multiple inflows",
	# 	# caption = ID
	# )

	# p <- p + theme_minimal(base_size = 11)

	# p <- p + theme(
	# 	legend.justification=c(0,1), 
	# 	legend.position=c(0.05,0.98),
	# 	legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
	# 	legend.key.size = unit(0.15, "cm"),    # Decrease key size
	# 	legend.text = element_text(size = 8), # Decrease text size
	# 	legend.title = element_text(size = 8), # Decrease title size
	# )

	# # show only one decimal in the legend labels
	# p <- p + scale_color_discrete(labels = function(x) {
	# 	x <- as.numeric(x)
	# 	x <- round(x, 1)
	# 	return(x)
	# })

	height	<- 50
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