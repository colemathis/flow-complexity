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

ID         				<<- "33-distance-from-source"
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
GRAPHS_PATH 			<<- file.path(DATA_DIR, "graphs.csv")
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

	ts_all <- append_distance_data(ts_all, graphs_csv = GRAPHS_PATH)

	# Save the processed data to a cache file
	dir.create(dirname(CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
	write_csv(ts_all, CACHE_PATH)

	return(ts_all)

}

#==============================================================================#

process_data <- function(ts) {

	ts <- ts %>%
		filter(sim_number %in% SELECTED_SIMS)

	# if ts is empty, return early
	if (nrow(ts) == 0) {
		return(ts)
	}

	# add assembly indices
	missing_value = 21
    ts <- ts %>%
        left_join(assembly_indices, by = "integer") %>%
        mutate(assemblyindex = ifelse(is.na(assemblyindex), missing_value, assemblyindex))

    ts <- ts %>%
        group_by(sim_number, chemostat_id) %>%
        summarize(
            mean_ai = sum(assemblyindex * frequency) / sum(frequency),
            sd_ai   = sqrt(sum(frequency * (assemblyindex - mean(assemblyindex))^2) / sum(frequency)),
            .groups = "drop"
        )

	# add diffusion rate
	ts <- ts %>%
		left_join(params %>% select(sim_number, diffusion_rate), by = "sim_number")

}

append_distance_data <- function(ts,
                                 graphs_csv = GRAPHS_CSV) {
    library(igraph)  # ensure igraph is available inside non‑interactive calls
    
    graphs <- read.csv(graphs_csv)
    sim_numbers <- unique(ts$sim_number)
    
    dist_tbl <- lapply(sim_numbers, function(s) {
        g_edges <- graphs %>%
            filter(sim_number == s) %>%
            select(chemostat_in, chemostat_out)
        
        # skip if graph is missing
        if (nrow(g_edges) == 0) return(NULL)
        
        # build an **undirected** graph so distance is symmetric
        g <- graph_from_data_frame(g_edges, directed = FALSE)
        
        # ensure vertex "1" (the source chemostat) exists
        if (!"1" %in% V(g)$name) return(NULL)
        
        d_vec <- distances(g, v = "1")[1, ]
        
        tibble(
            sim_number   = s,
            chemostat_id = as.integer(names(d_vec)),
            distance     = as.numeric(d_vec)
        )
    }) %>% bind_rows()
    
    ts %>% left_join(dist_tbl, by = c("sim_number", "chemostat_id"))
}

#==============================================================================#

load_cached_data <- function() {

	read_csv(CACHE_PATH, show_col_types = FALSE)

}

#==============================================================================#

plot_figure <- function(ts) {

  # replace diffusion_rate in ts for sim_number 30 to -4
ts <- ts %>%
	mutate(diffusion_rate = case_when(
		sim_number == 30 ~ 1e-4,
		sim_number == 58 ~ 1e-2,
		sim_number == 86 ~ 1e0,
		TRUE ~ diffusion_rate
	))

library(scico)
cols <- scico(3, palette = "navia", end = 0.5)

p <- ggplot(ts, aes(x = distance, y = mean_ai, color = factor(log10(diffusion_rate)))) +
	geom_point(alpha = 0.40) +
	geom_line(data = ts %>% filter(sim_number == 30),
						stat = "smooth", method = "loess", span = 1.0, se = FALSE, size = 0.75, alpha = 0.75) +
	geom_line(data = ts %>% filter(sim_number == 58),
						stat = "smooth", method = "loess", span = 0.75, se = FALSE, size = 0.75, alpha = 0.75) +
	geom_line(data = ts %>% filter(sim_number == 86),
						stat = "smooth", method = "loess", span = 4.20, se = FALSE, size = 0.75, alpha = 0.75) +
	# scale_color_scico_d(palette = "oleron") +
	# scale_colour_manual(values = cols) +
	scale_color_viridis_d(
		option = "turbo",   # "magma" or "plasma" give a dark-to-bright progression
		begin = 0.6,        # start slightly lighter than pure black
		end = 1.0,          # stop before pure white
		direction = -1,      # 1 = dark→light
		name = TeX("$\\log(k_d)$"),
		# labels = function(x) sprintf("%.0f", log10(as.numeric(x)))
	) +
    labs(
      x = TeX("Distance from source $d$"),
      y = "Mean Assembly Index",
      color = TeX("$\\log_{10}(k_d)$")
    ) +
	coord_cartesian(ylim = c(1.5, 12.5)) +
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 6.0, ymax = 9.0,
             color = "darkgreen", fill = NA, 
			 size = 0.5, linetype = "32") +
    theme_minimal(base_size = 11)
	p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
	p <- p + theme(
		legend.position = c(0.02, 0.90),
		legend.justification = "left",
		legend.background = element_rect(fill = "grey95", color = NA),
		legend.text = element_text(size = 8),
		legend.title = element_text(size = 8),
		legend.direction = "horizontal"
	)
	p <- p + guides(color = guide_legend(keyheight = unit(0.5, "lines"), keywidth = unit(0.45, "lines"), default.unit = "lines"))

	p <- p +
			theme(legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25))

	# Add dist.pdf as an inset in the top-right corner
	dist_file <- file.path("dist.pdf")
	# if (file.exists(dist_file)) {
		# dist_plot <- cowplot::ggdraw() + cowplot::draw_image(dist_file, scale = 0.6)
	dist_plot <- cowplot::ggdraw() + cowplot::draw_image("dist.png", scale = 0.6)

	p <- cowplot::ggdraw() +
		cowplot::draw_plot(p) +
		cowplot::draw_plot(dist_plot, x = 0.59, y = 0.54, width = 0.52, height = 0.52)
	# } else {
	# 	warning("File 'dist.pdf' not found. Skipping inset.")
	# }

	# p <- ts %>%
	# 	ggplot(aes(x = time, y = total_assembly, color = factor(diffusion_rate))) +
	# 	geom_point(size = 0.5, alpha = 0.25) +
	# 	geom_line(stat = "smooth", method = "loess", span = 0.30, se = FALSE, size = 0.75, alpha = 0.75) +
	# 	scale_y_log10(
	# 		labels = scales::trans_format("log10", function(x) TeX(sprintf("$10^{%d}$", x)))
	# 	) +
	# 	scale_color_discrete(
	# 		name = TeX("$\\log(k_d)$"),
	# 		labels = function(x) sprintf("%.0f", log10(as.numeric(x)))
	# 	) +
	# 	labs(
	# 		x = TeX("$t$"),
	# 		y = TeX("$A$")
	# 	) +
	# 	scale_x_continuous(
	# 		breaks = seq(0, 1e5, by = 2e4),
	# 		labels = function(x) ifelse(x %in% c(0, 1e5), c("0", TeX("$10^5$")), "")
	# 	) +
	# 	coord_cartesian(ylim = c(1e2, 1e6)) +
	# 	theme_minimal(base_size = 11) +
	# 	theme(
	# 		panel.grid.minor = element_blank(),
	# 		# panel.grid.major.x = element_blank(),
	# 		# legend.position = c(0.85, 0.3),
	# 		# legend.background = element_rect(fill = "white", color = "black"),
    #   		legend.background = element_rect(fill = "grey95", color = NA),	
	# 		legend.title = element_text(size = 6),
	# 		legend.text = element_text(size = 6),
	# 		legend.direction = "horizontal"
	# 	) +
	# 	# theme(
	# 	# 	panel.grid.major = element_blank(),
	# 	# 	panel.grid.minor = element_blank()
	# 	# )
    #     theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
	# 	p <- p + guides(color = guide_legend(keyheight = unit(0.5, "lines"), keywidth = unit(1, "lines"), default.unit = "lines"))


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