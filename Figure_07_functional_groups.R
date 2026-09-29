# ==============================================================================
# FIGURE 7 — PLANT FUNCTIONAL-GROUP STRATIGRAPHY
# ==============================================================================

source("00_config.R")
library(tidyverse)

community_long <- readRDS(file.path(DERIVED_DIR, "community_long_primary.rds"))
samples <- readRDS(file.path(DERIVED_DIR, "sample_info.rds"))

functional_long <- community_long |>
  filter(!is.na(plot_group)) |>
  mutate(plot_group = factor(plot_group, levels = FUNCTIONAL_GROUP_LEVELS))

group_counts <- functional_long |>
  group_by(library_id, plot_group) |>
  summarise(group_count = sum(genus_count, na.rm = TRUE), .groups = "drop")

library_totals <- community_long |>
  group_by(library_id) |>
  summarise(
    retained_plant_count_total = max(retained_plant_genus_count_total, na.rm = TRUE),
    .groups = "drop"
  )

group_metrics <- expand_grid(
  library_id = samples$library_id,
  plot_group = FUNCTIONAL_GROUP_LEVELS
) |>
  left_join(
    samples |> select(library_id, core, master_depth_cm),
    by = "library_id"
  ) |>
  left_join(group_counts, by = c("library_id", "plot_group")) |>
  left_join(library_totals, by = "library_id") |>
  mutate(
    group_count = replace_na(group_count, 0),
    relative_value = if_else(
      retained_plant_count_total > 0,
      group_count / retained_plant_count_total,
      NA_real_
    ),
    plot_group = factor(plot_group, levels = FUNCTIONAL_GROUP_LEVELS)
  )

profile_data <- group_metrics |>
  filter(!is.na(relative_value)) |>
  group_by(plot_group, core, master_depth_cm) |>
  summarise(relative_value = mean(relative_value, na.rm = TRUE), .groups = "drop") |>
  arrange(plot_group, core, master_depth_cm)

functional_labels <- c(
  "Trees" = "Trees",
  "Shrubs" = "Shrubs",
  "Forbs" = "Forbs",
  "Graminoids" = "Graminoids",
  "Ferns" = "Ferns",
  "Bryophytes" = "Bryophytes",
  "Submerged macrophytes" = "Submerged\nmacrophytes",
  "Emergent & wetland" = "Emergent &\nwetland"
)

functional_breaks <- function(limits) {
  upper <- max(limits, na.rm = TRUE)
  if (upper <= 0.02) {
    c(0, 0.005)
  } else if (upper <= 0.08) {
    c(0, 0.05)
  } else if (upper <= 0.30) {
    c(0, 0.10)
  } else {
    c(0, 0.50)
  }
}

functional_percent_labels <- function(x) {
  if (any(x > 0 & x < 0.01, na.rm = TRUE)) {
    scales::percent(x, accuracy = 0.1)
  } else {
    scales::percent(x, accuracy = 1)
  }
}

p_functional <- ggplot(
  profile_data,
  aes(
    y = master_depth_cm,
    xmin = 0,
    xmax = relative_value,
    group = interaction(plot_group, core),
    fill = plot_group
  )
) +
  geom_ribbon(alpha = 0.82, orientation = "y", na.rm = TRUE) +
  facet_wrap(
    ~plot_group,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(plot_group = functional_labels)
  ) +
  scale_fill_manual(values = FUNCTIONAL_GROUP_COLORS, drop = FALSE) +
  scale_x_continuous(
    breaks = functional_breaks,
    labels = functional_percent_labels,
    expand = expansion(mult = c(0, 0.04))
  ) +
  depth_scale_with_age() +
  labs(
    title = "Plant functional-group representation through the sediment sequence",
    x = "Share of retained plant reads",
    y = "Composite master depth (cm)",
    fill = "Functional group"
  ) +
  thesis_theme(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 25, margin = margin(b = 12)),
    strip.text = element_text(face = "bold", size = 17, lineheight = 0.95),
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 20, margin = margin(l = 10)),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 17),
    axis.text.y.right = element_text(size = 17),
    panel.spacing.x = grid::unit(1.15, "lines"),
    legend.position = "none",
    plot.margin = margin(t = 12, r = 24, b = 14, l = 14)
  )

save_figure(
  p_functional,
  "Figure_07_functional_groups",
  width = 18,
  height = 7.5
)
