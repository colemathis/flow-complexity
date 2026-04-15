################################################################################
# Residence Time Analysis — Plotting
#
# Reads CSV output from residence_time.jl and produces two SI figures:
#   1. fpt_comparison.pdf        — Mean FPT per randomized realization vs lattice
#   2. calibration_fpt_vs_topology.pdf — Gap decomposition: total vs FPT-explained
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

SAVE_FIGS <- TRUE
FIGS_DIR  <- "figs"

################################################################################
# LOAD DATA
################################################################################

fpt     <- read_csv("data/fpt_summary.csv", show_col_types = FALSE)
cal     <- read_csv("data/calibration.csv", show_col_types = FALSE)
meta    <- read_csv("data/meta.csv", show_col_types = FALSE)

################################################################################
# PLOT 1 — FPT Comparison (sorted scatter + reference lines)
################################################################################

fpt_rand <- fpt %>%
  filter(topology != "lattice") %>%
  filter(!is.na(mean_fpt)) %>%
  arrange(mean_fpt) %>%
  mutate(rank = row_number())

lattice_mean <- meta$lattice_mean_fpt
rand_grand_mean <- meta$randomized_mean_fpt

p1 <- ggplot(fpt_rand, aes(x = rank, y = mean_fpt)) +
  geom_point(colour = "#E06050", size = 1.5, alpha = 0.7) +
  geom_hline(aes(yintercept = lattice_mean, linetype = "Lattice mean FPT"),
             colour = "#4682B4", linewidth = 0.8) +
  geom_hline(aes(yintercept = rand_grand_mean, linetype = "Randomized mean"),
             colour = "#E06050", linewidth = 0.8) +
  scale_linetype_manual(
    name = NULL,
    values = c("Lattice mean FPT" = "dashed", "Randomized mean" = "dotted")
  ) +
  labs(
    x = "Realization (sorted by mean FPT)",
    y = "Mean FPT (steps)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = "grey70")
  )

if (SAVE_FIGS) {
  ggsave(
    file.path(FIGS_DIR, "fpt_comparison.pdf"),
    plot = p1, width = 80, height = 60, units = "mm"
  )
  cat("Saved fpt_comparison.pdf\n")
}

################################################################################
# PLOT 2 — Calibration: gap decomposition
################################################################################

shift_label <- sprintf("\u00d7%.2f shift", meta$fpt_shift_factor)

p2 <- ggplot(cal, aes(x = log_kd)) +
  geom_ribbon(aes(ymin = 0, ymax = total_gap_smooth, fill = "Total topology gap"),
              alpha = 0.25) +
  geom_line(aes(y = total_gap_smooth, colour = "Total topology gap"),
            linewidth = 0.6) +
  geom_ribbon(aes(ymin = 0, ymax = fpt_gap_smooth, fill = "FPT-explained gap"),
              alpha = 0.35) +
  geom_line(aes(y = fpt_gap_smooth, colour = "FPT-explained gap"),
            linewidth = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
  scale_colour_manual(
    name = NULL,
    values = c("Total topology gap" = "#E06050", "FPT-explained gap" = "#4682B4"),
    labels = c(
      "Total topology gap" = "Total topology gap",
      "FPT-explained gap"  = paste0("FPT-explained gap (", shift_label, ")")
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Total topology gap" = "#E06050", "FPT-explained gap" = "#4682B4"),
    labels = c(
      "Total topology gap" = "Total topology gap",
      "FPT-explained gap"  = paste0("FPT-explained gap (", shift_label, ")")
    )
  ) +
  labs(
    x = TeX("$\\log_{10}(k_d)$"),
    y = TeX("$\\Delta\\alpha$ (lattice $-$ randomized)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = "grey70"),
    legend.text = element_text(size = 7)
  )

if (SAVE_FIGS) {
  ggsave(
    file.path(FIGS_DIR, "calibration_fpt_vs_topology.pdf"),
    plot = p2, width = 80, height = 60, units = "mm"
  )
  cat("Saved calibration_fpt_vs_topology.pdf\n")
}
