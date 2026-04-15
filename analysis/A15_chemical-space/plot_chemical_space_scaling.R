################################################################################
# Chemical Space Scaling — Integer vs. Molecular Assembly Space
#
# Plots achievable assembly index (AI) as a function of chemical space size
# for integer chemistry (used in the simulations) vs. real molecular chemistry
# (from GDB-17 empirical data, Ruddigkeit et al. 2012).
#
# Key result: for any given space size, integer chemistry reaches higher AI
# than molecular chemistry, because molecular space grows as ~alpha^k
# (alpha ≈ 5.3) vs. 2^k for integers. This directly shows that our integer
# model overestimates achievable AI.
################################################################################

options(conflicts.policy = list(warn.conflicts = FALSE))

library(tidyverse)
library(latex2exp)

################################################################################
# PARAMETERS
################################################################################

SAVE_FIGS <- TRUE
FIGS_DIR  <- "figs"
K_MAX     <- 20   # extend predictions beyond GDB data (GDB goes to k=17)

################################################################################
# LOAD DATA
################################################################################

gdb <- read_csv("data/gdb_molecule_counts.csv", show_col_types = FALSE)

################################################################################
# COMPUTE CUMULATIVE COUNTS AND INTEGER SPACE
################################################################################

df <- gdb %>%
  arrange(hac) %>%
  mutate(
    n_mol_cumul = cumsum(n_molecules),
    k           = hac,               # AI ≈ heavy-atom count (conservative)
    n_int       = 2^k                # integer-space upper bound
  )

# Fit log-linear model to cumulative molecular counts: log10(N) = a + b*k
fit <- lm(log10(n_mol_cumul) ~ k, data = df)
alpha_log10 <- coef(fit)[["k"]]      # slope on log10 scale
alpha       <- 10^alpha_log10         # base of exponential growth
intercept   <- coef(fit)[["(Intercept)"]]

cat(sprintf("Fit: log10(N_mol) = %.3f + %.3f * k\n", intercept, alpha_log10))
cat(sprintf("Molecular growth base alpha = %.2f (cf. integer base = 2)\n", alpha))
cat(sprintf("Ratio base alpha/2 = %.2f  =>  R(k) ~ %.2f^k\n", alpha / 2, alpha / 2))

# Extend to prediction range
df_pred <- tibble(
  k = seq(1, K_MAX),
  log10_n_int = k * log10(2),
  log10_n_mol = intercept + alpha_log10 * k
)

# Define both curves over a shared log10(N) range for the flipped-axis plot
logN_max <- max(df_pred$log10_n_mol)  # molecular curve extends further
logN_seq <- seq(0, logN_max, length.out = 200)

df_int_curve <- tibble(
  log10_N = logN_seq,
  k = sqrt(2) * log10_N / log10(2)    # k = sqrt(2) * log2(N)  [paper upper bound]
)

df_mol_curve <- tibble(
  log10_N = seq(log10(df$n_mol_cumul[1]), logN_max, length.out = 200),
  k = (log10_N - intercept) / alpha_log10  # invert the fit
)

################################################################################
# PLOT — Achievable AI vs. space size (flipped axes)
#
# x-axis: log10(N) — size of chemical space explored
# y-axis: k — maximum assembly index reachable
#
# Reading: for any given space size, integer chemistry reaches higher AI.
################################################################################

# Aesthetics
col_int <- "#4682B4"   # blue
col_mol <- "#E06050"   # red/coral
col_fill <- "grey85"   # ribbon fill

# For the annotation: at a fixed log10(N) value, what AI does each curve reach?
annot_logN <- 10
annot_k_int <- sqrt(2) * annot_logN / log10(2)      # k = sqrt(2) * log2(N)
annot_k_mol <- (annot_logN - intercept) / alpha_log10 # k from the fit
annot_delta <- round(annot_k_int - annot_k_mol)

brace_mid <- mean(c(annot_k_int, annot_k_mol))

# Equation annotations along curves — compute visual tilt angles
slope_int_data <- sqrt(2) / log10(2)       # dk/d(log10 N) for integer curve
slope_mol_data <- 1 / alpha_log10          # dk/d(log10 N) for molecular curve

x_data_range <- logN_max
y_data_range <- max(df_int_curve$k) * 1.02  # approx upper y limit
aspect_ratio <- 60 / 80                    # plot height / width (mm)

angle_int <- atan(slope_int_data * aspect_ratio * x_data_range / y_data_range) * 180 / pi
angle_mol <- atan(slope_mol_data * aspect_ratio * x_data_range / y_data_range) * 180 / pi

# Center each annotation along the curve
eq_x       <- logN_max / 2
eq_k_int   <- sqrt(2) * eq_x / log10(2)
eq_k_mol   <- (eq_x - intercept) / alpha_log10
eq_offset  <- 2.0  # vertical offset above/below curve

# GDB data range boundary
gdb_logN_max <- log10(sum(gdb$n_molecules))  # log10(166.4B) ≈ 11.2

# Ribbon: shade the gap between curves
df_ribbon <- tibble(
  log10_N = logN_seq,
  k_int = sqrt(2) * log10_N / log10(2),
  k_mol = (log10_N - intercept) / alpha_log10
)

int_label <- TeX("Integer")
mol_label <- TeX("Molecular")

p <- ggplot() +
  # 1. Shaded gap between curves
  geom_ribbon(
    data = df_ribbon,
    aes(x = log10_N, ymin = k_mol, ymax = k_int),
    fill = col_fill, alpha = 0.5
  ) +
  # Integer space line — mapped to legend
  geom_line(
    data = df_int_curve,
    aes(x = log10_N, y = k, colour = "int", linetype = "int"),
    linewidth = 0.5
  ) +
  # Molecular space line — mapped to legend
  geom_line(
    data = df_mol_curve,
    aes(x = log10_N, y = k, colour = "mol", linetype = "mol"),
    linewidth = 0.5
  ) +
  # GDB empirical data points
  geom_point(
    data = df,
    aes(x = log10(n_mol_cumul), y = k),
    colour = col_mol, size = 1.2, alpha = 0.8
  ) +
  # 2. AT detection threshold at AI = 15
  geom_hline(yintercept = 15, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
  annotate(
    "text", x = 0.3, y = 17.5,
    label = "AI = 15",
    size = 3, hjust = 0, colour = "grey40"
  ) +
  # 5. Mark GDB data range
  geom_vline(xintercept = gdb_logN_max, linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  annotate(
    "text", x = gdb_logN_max - 0.2, y = 2.5,
    label = "GDB-17", size = 3, colour = "grey50", hjust = 1, fontface = "italic"
  ) +
  # Delta annotation
  annotate(
    "segment",
    x = annot_logN, xend = annot_logN,
    y = annot_k_mol, yend = annot_k_int,
    colour = "black", linewidth = 1.0
  ) +
  # Label box centered on line
  annotate(
    "rect",
    xmin = annot_logN - 0.8, xmax = annot_logN + 0.8,
    ymin = brace_mid - 1.8, ymax = brace_mid + 1.8,
    fill = "grey92", colour = NA
  ) +
  annotate(
    "text",
    x = annot_logN,
    y = brace_mid,
    label = TeX(sprintf("$\\delta A \\approx %d$", annot_delta)),
    size = 3, colour = "black", fontface = "bold"
  ) +
  # Equation labels along curves
  annotate(
    "text", x = eq_x - 1, y = eq_k_int + eq_offset - 5,
    label = TeX("$k \\approx \\sqrt{2}\\,\\log_2 N$"),
    angle = angle_int + 3, size = 3, colour = col_int, vjust = 0
  ) +
  annotate(
    "text", x = eq_x, y = eq_k_mol - eq_offset,
    label = TeX("$k \\approx 1.38 \\times \\log_{10} N$"),
    angle = angle_mol, size = 3, colour = col_mol, vjust = 1
  ) +
  # Scales — legend with equations
  scale_colour_manual(
    name = NULL,
    values = c("int" = col_int, "mol" = col_mol),
    labels = c("int" = int_label, "mol" = mol_label)
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("int" = "solid", "mol" = "dashed"),
    labels = c("int" = int_label, "mol" = mol_label)
  ) +
  scale_y_continuous(
    breaks = seq(0, 60, by = 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = TeX("Reachable chemical space $\\log_{10} (N)$"),
    y = "Assembly index A"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = alpha("white", 0.9), colour = "grey70"),
  )

if (SAVE_FIGS) {
  dir.create(FIGS_DIR, showWarnings = FALSE)
  ggsave(
    file.path(FIGS_DIR, "chemical-space-scaling.pdf"),
    plot = p, width = 80, height = 60, units = "mm"
  )
  cat("Saved chemical-space-scaling.pdf\n")
}

################################################################################
# PRINT SUMMARY TABLE
################################################################################

summary_tbl <- df %>%
  transmute(
    k,
    log10_N_int = round(log10(n_int), 2),
    log10_N_mol = round(log10(n_mol_cumul), 2),
    log10_R     = round(log10(n_mol_cumul) - log10(n_int), 2)
  )

cat("\n=== Summary: Integer vs. Molecular Space ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat(sprintf("\nAnalytical approximation: R(k) ≈ (%.1f/2)^k ≈ %.1f^k\n", alpha, alpha / 2))
