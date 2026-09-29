# ==============================================================================
# FIGURE 1 — CLIMATIC CONTEXT AND SAMPLING COVERAGE
# ==============================================================================

source("00_config.R")

age_model <- load_age_model()
samples <- readRDS(file.path(DERIVED_DIR, "sample_info.rds")) |>
  dplyr::arrange(master_depth_cm, library_id)

age_to_depth <- function(age_ka) {
  stats::approx(
    x = age_model$age_calBP / 1000,
    y = age_model$depth_cm,
    xout = age_ka,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

min_depth <- min(age_model$depth_cm, na.rm = TRUE)
max_depth <- max(age_model$depth_cm, na.rm = TRUE)
boundary_ages_ka <- c(11.7, 12.9, 14.7)
boundary_depths <- age_to_depth(boundary_ages_ka)

intervals <- data.frame(
  label = c("Early\nHolocene", "Younger\nDryas", "Bølling–\nAllerød", "Older Late\nGlacial"),
  top = c(min_depth, boundary_depths[1], boundary_depths[2], boundary_depths[3]),
  bottom = c(boundary_depths[1], boundary_depths[2], boundary_depths[3], max_depth),
  stringsAsFactors = FALSE
)

interval_fills <- c("#dce6f1", "#eee5c4", "#dcebd5", "#e6dceb")
depth_ticks <- seq(0, 600, by = 100)

max_modeled_age_bp <- max(age_model$age_calBP, na.rm = TRUE)
plateau_depth <- min(
  age_model$depth_cm[abs(age_model$age_calBP - max_modeled_age_bp) < 1e-8],
  na.rm = TRUE
)

age_axis_depths <- c(min_depth, boundary_depths[1], boundary_depths[2], boundary_depths[3])
age_axis_labels <- c("9.5", "11.7", "12.9", "14.7")

samples$x_offset <- 0
split_idx <- split(seq_len(nrow(samples)), samples$master_depth_cm)
for (idx in split_idx) {
  if (length(idx) > 1) {
    samples$x_offset[idx] <- seq(-0.025, 0.025, length.out = length(idx))
  }
}

make_figure <- function() {
  par(mar = c(3.0, 5.2, 5.3, 7.2), xpd = NA)

  plot(
    NA, NA,
    xlim = c(0.3, 4.7),
    ylim = c(max_depth, min_depth),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )

  climate_left <- 1.10
  climate_right <- 2.45
  sample_left <- 3.20
  sample_right <- 3.68
  sample_center <- (sample_left + sample_right) / 2

  for (i in seq_len(nrow(intervals))) {
    rect(
      climate_left,
      intervals$bottom[i],
      climate_right,
      intervals$top[i],
      col = interval_fills[i],
      border = "grey35",
      lwd = 1.1
    )
    text(
      (climate_left + climate_right) / 2,
      (intervals$top[i] + intervals$bottom[i]) / 2,
      labels = intervals$label[i],
      cex = 0.95,
      font = 2
    )
  }

  text(
    (climate_left + climate_right) / 2,
    min_depth - 15,
    labels = "Climatic interval",
    cex = 1.0,
    font = 2
  )

  rect(
    sample_left,
    max_depth,
    sample_right,
    min_depth,
    col = "white",
    border = "grey35",
    lwd = 1.0
  )

  segments(
    sample_center - 0.12 + samples$x_offset,
    samples$master_depth_cm,
    sample_center + 0.12 + samples$x_offset,
    samples$master_depth_cm,
    col = "grey25",
    lwd = 0.9
  )

  text(
    sample_center,
    min_depth - 34,
    labels = paste0("sedaDNA\nsamples\n(n = ", nrow(samples), ")"),
    cex = 0.92,
    font = 2
  )

  for (d in boundary_depths) {
    segments(
      climate_left - 0.12,
      d,
      sample_right + 0.12,
      d,
      col = "grey70",
      lty = 2,
      lwd = 0.9
    )
  }

  axis(side = 2, at = depth_ticks, labels = depth_ticks, las = 1, cex.axis = 0.95)
  mtext("Composite depth (cm)", side = 2, line = 3.6, cex = 1.05)

  axis(
    side = 4,
    at = age_axis_depths,
    labels = age_axis_labels,
    las = 1,
    cex.axis = 0.95,
    lwd = 0,
    lwd.ticks = 1
  )

  usr <- par("usr")
  age_axis_x <- usr[2]
  segments(age_axis_x, min_depth, age_axis_x, plateau_depth, col = "grey20", lwd = 1.0)

  bracket_x <- age_axis_x + 0.16
  segments(bracket_x, plateau_depth, bracket_x, max_depth, col = "grey50", lwd = 0.9)
  segments(bracket_x, plateau_depth, bracket_x + 0.07, plateau_depth, col = "grey50", lwd = 0.9)
  segments(bracket_x, max_depth, bracket_x + 0.07, max_depth, col = "grey50", lwd = 0.9)

  text(
    bracket_x + 0.10,
    (plateau_depth + max_depth) / 2,
    labels = "17.2 ka\n(constant)",
    cex = 0.63,
    col = "grey35",
    adj = c(0, 0.5)
  )

  mtext("Working modeled age (ka cal BP)", side = 4, line = 4.7, cex = 1.05)
}

png(
  filename = file.path(PLOT_DIR, "Figure_01_climatic_context.png"),
  width = 1800,
  height = 2600,
  res = 300
)
make_figure()
dev.off()

pdf(
  file = file.path(PLOT_DIR, "Figure_01_climatic_context.pdf"),
  width = 6.0,
  height = 8.5
)
make_figure()
dev.off()
