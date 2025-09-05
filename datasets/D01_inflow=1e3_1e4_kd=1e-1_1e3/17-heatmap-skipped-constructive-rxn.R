################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

ID         				<<- "17-heatmap-skipped-constructive-rxn"
USE_CACHE  				<<- TRUE
PRINT_FIGS 				<<- TRUE
SAVE_FIGS 			 	<<- TRUE

DATA_DIR      			<<- "../../datasets/D01_inflow=1e3_1e4_kd=1e-1_1e3/data"
CACHE_DIR     			<<- file.path("cache", ID)
FIGS_DIR      			<<- "figs"

PARAMS_FILE 			<<- "params.csv"
TIMESERIES_FILES		<<- "timeseries.csv"
META_FILES				<<- "meta.csv"
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

load_and_process_sim_data <- function() {

	timeseries_files <- list.files(DATA_DIR, pattern = TIMESERIES_FILES, recursive = TRUE, full.names = TRUE)
	meta_files       <- list.files(DATA_DIR, pattern = META_FILES, recursive = TRUE, full.names = TRUE)

	ts_all <- map_dfr(meta_files, function(file) {
		if (file.info(file)$size == 0) return(NULL)
		ts <- read_csv(file, show_col_types = FALSE, progress = FALSE)
		# process_data(ts)
	}, .progress = TRUE)

	# Save the processed data to a cache file
	dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
	write_csv(ts_all, CACHE_PATH)

	return(ts_all)

}

#==============================================================================#

# process_data <- function(ts) {

# 	MAX_TIME   <- params$total_time[1]
# 	N_REACTORS <- params$N_reactors[1]

# 	ts <- ts %>%
# 		filter(time == MAX_TIME) %>%
# 		filter(integer == 2)

# 	ts <- ts %>%
# 		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
# 		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

# 	ts %>%
# 		group_by(diffusion_rate, inflow_mols, integer) %>%
# 		summarize(
# 			mean_frequency = sum(frequency, na.rm = TRUE) / N_REACTORS,
# 			sd_frequency = sqrt(sum((frequency - mean_frequency)^2) / (N_REACTORS - 1)),
# 			.groups = "drop"
# 		)

# }

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

	# Calculate the size of the square heatmap
	n <- nrow(ts)
	side <- ceiling(sqrt(n))

	# Pad the data if necessary to make a perfect square
	if (n < side^2) {
		pad <- side^2 - n
		ts <- bind_rows(ts, tibble(
			x = rep(NA, pad),
			y = rep(NA, pad),
			total_time = rep(NA, pad),
			sim_number = rep(NA, pad)
		))
	}

	# Assign x and y positions for square layout: left-to-right, then top-to-bottom,
	# with (1,1) in the upper left corner (y decreases downward)
	ts$x <- rep(1:side, times = side)[1:nrow(ts)]
	ts$y <- rep(side:1, each = side)[1:nrow(ts)]

	p <- ggplot(ts, aes(x = x, y = y, fill = 100*skipped_constructive_rxn/total_constructive_rxn)) +
		geom_tile(color = "white", na.rm = TRUE) +
		# geom_text(aes(label = sim_number), size = 2, na.rm = TRUE) +
		scale_fill_viridis_c(
			name = "% skipped",
			option = "C",
      		limits = c(0, 100),                       # force 0 → 1 range
			na.value = "grey90",
			guide = guide_colorbar(
				barwidth = unit(3, "mm"),
				barheight = unit(20, "mm"),
				title.position = "top",
				title.hjust = 0.5
			)
		) +
		labs(
			title = "Constructive",
			# caption = ID
		) +
		theme_minimal(base_size = 8) +
		theme(
			panel.grid = element_blank(),
			plot.title = element_text(size = 10, hjust = 0.5),
			plot.caption = element_text(size = 7, color = "grey50"),
			legend.position = "right",
			legend.title = element_text(size = 7),
			legend.text = element_text(size = 6),
			legend.key.height = unit(10, "mm"),
			legend.key.width = unit(2, "mm"),
			axis.text = element_blank(),
			axis.ticks = element_blank(),
			axis.title = element_blank()
		) +
		coord_fixed()

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
saveRDS(p, file = file.path(CACHE_DIR, paste0(ID, ".rds")))