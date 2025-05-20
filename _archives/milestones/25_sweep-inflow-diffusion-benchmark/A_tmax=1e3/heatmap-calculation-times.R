#!/usr/bin/env Rscript
# --- prerequisites -----------------------------------------------------------
library(tidyverse)    # supplies purrr, dplyr, tidyr, ggplot2
library(stringr)      # string helpers
# --- collect & parse log files ----------------------------------------------
log_files <- list.files("data/logs", pattern = "\\.log$", full.names = TRUE)

log_df <- map_dfr(log_files, function(path) {
  sim_id <- str_extract(basename(path), "(?<=_)[0-9]+(?=\\.log$)") %>% as.integer()
  
  lines  <- readLines(path, warn = FALSE)
  hit    <- str_subset(lines, "^Sim Completed\\.")
  
  secs <- if (length(hit)) str_extract(hit[1], "[0-9]+\\.?[0-9]*") %>% as.numeric() else NA_real_
  
  tibble(simulation_id = sim_id,
         calculation_time = secs / 3600)         # hours
})

# --- ensure full 1–100 grid --------------------------------------------------
heat_df <- tibble(simulation_id = 1:100) %>%                # guaranteed grid
  left_join(log_df, by = "simulation_id") %>%
  mutate(
    row = ceiling(simulation_id / 10),                      # 1-based rows
    col = ((simulation_id - 1) %% 10) + 1                   # 1-based cols
  )

# --- 10 × 10 heat-map --------------------------------------------------------
plot <- ggplot(heat_df, aes(col, row, fill = calculation_time)) +
    geom_tile(colour = "grey70") +
    coord_fixed() +
    scale_x_continuous(breaks = 1:10) +
    scale_y_reverse(breaks = 1:10) +                          # row 1 at top
    scale_fill_viridis_c(na.value = "white", name = "Hours", breaks = seq(1, 25, by = 4), limits = c(0, 25)) +
    labs(title = "Simulation Calculation Time (hours) — tmax=1e3",
             x = NULL, y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

# options(vsc.dev.args = list(width = 8, height = 8, res=300, units = "in"))
# print(plot)

out_file <- file.path("figs", sprintf("heatmap-calculation-times.pdf"))
ggsave(out_file, plot = plot, width = 8, height = 8, create.dir = TRUE)