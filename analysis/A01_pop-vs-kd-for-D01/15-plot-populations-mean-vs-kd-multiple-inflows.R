################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # includes dplyr, ggplot2, etc.
library(latex2exp)

################################################################################
# MAIN FUNCTION
################################################################################

main <- function() {

  ID         <<- "populations-mean-vs-kd-multiple-inflows"
  USE_CACHE  <<- FALSE
  PRINT_FIGS <<- TRUE
  SAVE_FIGS  <<- FALSE

  DATA_DIR      <<- "data"
  CACHE_DIR     <<- file.path("cache", ID)
  FIGS_DIR      <<- "figs"

  TIMESERIES_FILES <<- "timeseries.csv"
  CACHE_FILE       <<- paste0(ID, ".csv")
  FIGS_FILE        <<- paste0(ID, ".pdf")

  CACHE_PATH   <<- file.path(CACHE_DIR, CACHE_FILE)
  PARAMS_PATH  <<- file.path(DATA_DIR, "params.csv")

  # Set-up environment
  setwd_to_script_path()

  # Load and process data
  if (file.exists(CACHE_PATH) && USE_CACHE) {
    data   <- load_cached_data()
  } else {
    params <- load_params()
    data   <- load_and_process_time_series()
    save_cached_data(data)
  }

  # Plot the figure
  p <- plot_figure(data)
  save_figure(p)

}

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
		ts <- read_csv(file, show_col_types = FALSE, progress = FALSE)
		process_data(ts)
	}, .progress = TRUE)

	return(ts_all)
}

#==============================================================================#

process_data <- function(ts) {

	MAX_TIME   <- params$total_time[1]
	N_REACTORS <- params$N_reactors[1]

	ts <- ts %>%
		filter(time == MAX_TIME) %>%
		filter(integer == 2)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number") %>%
		left_join(params %>% select(sim_number, inflow_mols), by = "sim_number")

	ts %>%
		group_by(diffusion_rate, inflow_mols, integer) %>%
		summarize(
			mean_frequency = sum(frequency, na.rm = TRUE) / N_REACTORS,
			sd_frequency = sqrt(sum((frequency - mean_frequency)^2) / (N_REACTORS - 1)),
			.groups = "drop"
		)

}

#==============================================================================#

save_cached_data <- function(ts_all) {
	dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
	write_csv(ts_all, CACHE_PATH)
}

#==============================================================================#

load_cached_data <- function() {
	read_csv(CACHE_PATH, show_col_types = FALSE)
}

#==============================================================================#

plot_figure <- function(ts) {

	p <- ggplot(ts,
							aes(
									x = diffusion_rate, 
									y = mean_frequency, 
									color = factor(log10(inflow_mols))
							))

  p <- p + geom_point(size = 0.5, alpha = 0.25)
  p <- p + geom_smooth(method = "loess", span = 0.5, se = FALSE, size = 0.5)

  p <- p + scale_x_log10()
  p <- p + scale_color_manual(values = c("darkblue", "blue", "yellow", "red"))

  p <- p + labs(
    x = TeX("Diffusion coefficient"),
    y = "Avg. freq. (whole systems)",
    color = TeX("Log(I)"),
    title = "Populations means across sweep over diffusion, multiple inflows",
  )

  p <- p + labs(caption = ID) +
    theme(
      plot.caption = element_text(size = 6, color = "grey50")
    )

	p <- p + theme(
		legend.justification=c(0,1), 
		legend.position=c(0.05,0.95),
		legend.background = element_rect(fill = "white", color = "black"), # Optional: Customize legend background
		legend.key.size = unit(0.15, "cm"),    # Decrease key size
		# legend.text = element_text(size = 6), # Decrease text size
		# legend.title = element_text(size = 6) # Decrease title size
		plot.title = element_text(size = 6)
	)

}

#==============================================================================#

save_figure <- function(p_main) {

	height <- 70
	width <- 80

  if (PRINT_FIGS) {
    options(vsc.dev.args = list(width = width, height = height, res=300, units = "mm"))
    print(p_main)
  }

	if (SAVE_FIGS) {
		out_file <- file.path(FIGS_DIR, FIGS_FILE)
		ggsave(filename = out_file, plot = p_main, width = width, height = height, units = "mm", create.dir = TRUE)
	}

}

################################################################################
# MAIN SCRIPT
################################################################################

main()