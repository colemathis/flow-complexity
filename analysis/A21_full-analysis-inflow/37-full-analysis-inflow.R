################################################################################
# IMPORTS
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)
library(latex2exp)
library(cowplot)

################################################################################
# PARAMETERS
################################################################################

ID            <<- "37-full-analysis-inflow"
USE_CACHE     <<- TRUE
PRINT_FIGS    <<- TRUE
SAVE_FIGS     <<- TRUE

DATA_DIR      <<- "../../datasets/D01_inflow=1e3_1e4_kd=1e-1_1e3/data"
CACHE_DIR     <<- file.path("cache", ID)
FIGS_DIR      <<- "figs"

PARAMS_FILE   <<- "params.csv"
TIMESERIES_FILES <<- "timeseries.csv"
CACHE_FILE    <<- paste0(ID)
FIGS_FILE     <<- paste0(ID, ".pdf")

CACHE_PATH    <<- file.path(CACHE_DIR, CACHE_FILE)
PARAMS_PATH   <<- file.path(DATA_DIR, PARAMS_FILE)

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

	dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
	write_rds(ts_all, CACHE_PATH)

	return(ts_all)

}

#==============================================================================#

process_data <- function(ts) {

	MAX_TIME <- params$total_time[1]

	ts <- ts %>%
		filter(time == MAX_TIME, integer > 1)

	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate, inflow_mols), by = "sim_number")

	bin_edges <- unique(round(exp(seq(log(1), log(1e7), length.out = 51))))
	ts <- ts %>%
		mutate(integer_bin = cut(
			integer,
			breaks = bin_edges,
			include.lowest = TRUE,
			right = TRUE
		)) %>%
		group_by(inflow_mols, diffusion_rate, integer_bin) %>%
		summarise(frequency = sum(frequency), .groups = "drop")

	ts <- ts %>%
		complete(inflow_mols, diffusion_rate, integer_bin, fill = list(frequency = 0))

}

#==============================================================================#

load_cached_data <- function() {

	read_rds(CACHE_PATH)

}

#==============================================================================#

make_panel <- function(inflow_val, ts, show_legend = FALSE, show_yaxis = FALSE) {

	ts_panel <- ts %>%
		filter(inflow_mols == inflow_val)

	# Derive all bin levels from the same formula used in process_data, independent
	# of what values appear in this panel's data (avoids white tiles where
	# complete() in process_data left gaps for bins with zero observations).
	bin_edges    <- unique(round(exp(seq(log(1), log(1e7), length.out = 51))))
	midpts       <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2
	all_levels   <- levels(cut(midpts, breaks = bin_edges, include.lowest = TRUE, right = TRUE))
	lower_bounds <- as.numeric(gsub("[[(]([^,]+),.*", "\\1", all_levels))
	valid_levels <- all_levels[lower_bounds < 1e5]
	n_valid      <- length(valid_levels)
	break_idx    <- round(seq(1, n_valid, length.out = 6))

	# Explicit cross-join: guaranteed to cover all (kd × bin) combinations regardless
	# of factor vs character type issues with complete()
	ts_panel <- expand_grid(
		diffusion_rate = sort(unique(params$diffusion_rate)),
		integer_bin    = valid_levels
	) %>%
		left_join(
			ts_panel %>%
				mutate(integer_bin = as.character(integer_bin)) %>%
				select(diffusion_rate, integer_bin, frequency),
			by = c("diffusion_rate", "integer_bin")
		) %>%
		mutate(
			frequency = replace_na(frequency, 0),
			frequency = na_if(frequency, 0)
		)

	p <- ggplot(ts_panel, aes(x = log10(diffusion_rate), y = integer_bin, fill = frequency)) +
		geom_tile() +
		scale_x_continuous(
			breaks = -1:3,
			labels = TeX(sprintf("$10^{%d}$", -1:3))
		) +
		scale_fill_viridis_c(name = TeX("Copy number $n$ "), na.value = "grey35") +
		labs(
			x = NULL,
			y = NULL,
			title = TeX(sprintf("$\\log(I) =$%.1f", log10(inflow_val)))
		) +
		scale_y_discrete(
			limits = valid_levels,
			breaks = valid_levels[break_idx],
			labels = c(expression(10^0), expression(10^1), expression(10^2),
			           expression(10^3), expression(10^4), expression(10^5))
		) +
		coord_cartesian(xlim = c(-1, 3)) +
		theme_minimal(base_size = 11) +
		theme(
			plot.title    = element_text(size = 10, hjust = 0.5),
			panel.grid    = element_blank(),
			plot.margin   = margin(t = 2, r = 1, b = 2, l = 2, unit = "mm")
		)

	if (show_legend) {
		p <- p + theme(
			legend.position      = c(0.90, 0.95),
			legend.justification = c("right", "top"),
			legend.background    = element_rect(fill = alpha("white", 0.95), color = NA),
			legend.key.size      = unit(0.5, "lines"),
			legend.text          = element_text(size = 5),
			legend.title         = element_text(size = 5)
		)
	} else {
		p <- p + theme(legend.position = "none")
	}

	if (!show_yaxis) {
		p <- p + theme(
			axis.text.y   = element_blank(),
			axis.ticks.y  = element_blank(),
			plot.margin   = margin(t = 2, r = 1, b = 2, l = 0, unit = "mm")
		)
	}

	return(p)

}

#==============================================================================#

plot_figure <- function(ts) {

	inflow_vals <- sort(unique(ts$inflow_mols))
	panels <- list(
		make_panel(inflow_vals[1], ts, show_legend = TRUE,  show_yaxis = TRUE),
		make_panel(inflow_vals[2], ts, show_legend = FALSE, show_yaxis = FALSE),
		make_panel(inflow_vals[3], ts, show_legend = FALSE, show_yaxis = FALSE),
		make_panel(inflow_vals[4], ts, show_legend = FALSE, show_yaxis = FALSE)
	)

	p_grid <- plot_grid(plotlist = panels, nrow = 1, align = "h", axis = "tb",
	                    rel_widths = c(1.20, 1, 1, 1))

	# Shared axis labels
	x_label <- ggdraw() + draw_label(TeX("Diffusion coefficient $k_d$"), size = 11)
	y_label <- ggdraw() + draw_label(TeX("Integer $z$"), size = 11, angle = 90)

	p_with_xlab <- plot_grid(p_grid, x_label, ncol = 1, rel_heights = c(1, 0.09))
	p <- plot_grid(y_label, p_with_xlab, nrow = 1, rel_widths = c(0.03, 1))

	height <- 55
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

setwd_to_script_path()

params <- load_params()

if (file.exists(CACHE_PATH) && USE_CACHE) {
	data <- load_cached_data()
} else {
	data <- load_and_process_time_series()
}

p <- plot_figure(data)
saveRDS(p, file = file.path(CACHE_DIR, paste0(ID, ".rds")))
