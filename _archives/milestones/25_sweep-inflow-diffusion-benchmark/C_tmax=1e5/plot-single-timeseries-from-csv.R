############################
# IMPORTS
############################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)   # brings in ggplot2, dplyr, readr, scales, …
library(latex2exp)
library(arrow)

############################
# PARAMETERS
############################

TITLE        <- "First ten populations."
ID           <- "single-timeseries"
# selected_sim <- 25
USE_CACHE    <- FALSE

DATA_DIR  <- "data"
CACHE_DIR <- paste0("cache/", ID)
FIGS_FILE <- paste0(ID, "_")
FIGS_DIR  <- paste0("figs")

# TIMESERIES_ARROW <- file.path(DATA_DIR, "timeseries.arrow")

PARAMS_CSV       <- file.path(DATA_DIR, "params.csv")

############################
# FUNCTIONS
############################

load_processed_data <- function(sim_id) {
    cache_path <- file.path(CACHE_DIR, sprintf("sim_%d.csv", sim_id))

    if (file.exists(cache_path) && USE_CACHE) {
        read_csv(cache_path, show_col_types = FALSE)
    } else {
        TIMESERIES_CSV   <- file.path(DATA_DIR, "sims", sprintf("%06d", sim_id), "timeseries.csv")
        data <- tryCatch({
            open_dataset(TIMESERIES_CSV, format = "csv") %>%
                filter(sim_number == sim_id, integer %in% 1:10) %>%
                collect()
        }, error = function(e) {
            message(sprintf("Error loading dataset for sim_id %d: %s", sim_id, e$message))
            tibble()  # Return an empty dataset
        })

        dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
        write_csv(data, cache_path)
        data
    }
}

add_blind_data <- function(data, sim_id, grid_size) {
    blind <- expand.grid(
        sim_number   = sim_id,
        time         = 0:1,
        integer      = 0,
        chemostat_id = 1:(grid_size^2),
        frequency    = NA_real_
    )
    bind_rows(data, blind)
}

plot_timeseries <- function(data, grid_size, sim_params, sim_id) {
    ggplot(data, aes(x = time, y = frequency, color = factor(integer))) +
        geom_line(alpha = 0.7) +
        facet_wrap(
            ~ chemostat_id,
            ncol    = grid_size,
            nrow    = grid_size,
            scales  = "fixed",
            labeller = labeller(
                chemostat_id = function(x) paste("chemostat #", x)
            )
        ) +
        scale_y_log10(labels = scales::label_scientific()) +
        theme_bw() +
        theme(
            legend.position = "none",
            panel.border    = element_rect(color = "black", fill = NA)
        ) +
        ggtitle(
            TeX(sprintf(
                "%s Simulation %d: $log_{10}(I)=%.2f$, $log_{10}(k_d)=%.2f$",
                TITLE,
                sim_id,
                log10(sim_params$inflow_mols),
                log10(sim_params$diffusion_rate)
            ))
        ) +
        theme(plot.title = element_text(hjust = 0.5))
}

############################
# MAIN SCRIPT
############################

do_one_sim <- function(selected_sim) {

    params       <- read_csv(PARAMS_CSV, show_col_types = FALSE)

    grid_size   <- sqrt(params$N_reactors[params$sim_number == selected_sim])
    sim_params  <- params %>% filter(sim_number == selected_sim)

    processed_data <- load_processed_data(selected_sim)

    if (nrow(processed_data) == 0) {
        cat(sprintf("No data for sim_id %d\n", selected_sim))
        return()
    }

    processed_data <- processed_data %>% mutate(frequency = na_if(frequency, 0))
    processed_data <- add_blind_data(processed_data, selected_sim, grid_size)

    plot <- plot_timeseries(processed_data, grid_size, sim_params, selected_sim)

    # options(vsc.dev.args = list(width = 8, height = 8, res=300, units = "in"))
    # print(plot)

    out_file <- file.path(FIGS_DIR, sprintf("%s%d.pdf", FIGS_FILE, selected_sim))
    ggsave(out_file, plot = plot, width = 8, height = 8, create.dir = TRUE)

}

for (selected_sim in 1:100) {

    cat(sprintf("Processing sim_id %d\n", selected_sim))
    do_one_sim(selected_sim)

}

