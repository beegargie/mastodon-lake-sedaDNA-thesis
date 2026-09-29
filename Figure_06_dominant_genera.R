# ==============================================================================
# FIGURE 6 — DATASET-WIDE READ REPRESENTATION OF DOMINANT PLANT GENERA
# ==============================================================================

source("00_config.R")
library(tidyverse)

TOP_N <- 9L

community_wide <- readRDS(file.path(DERIVED_DIR, "community_counts_wide_primary.rds"))

mat <- community_wide |>
  select(-library_id) |>
  as.matrix()
storage.mode(mat) <- "numeric"

genus_summary <- tibble(
  genus = colnames(mat),
  authenticated_reads = colSums(mat, na.rm = TRUE)
) |>
  filter(authenticated_reads > 0) |>
  arrange(desc(authenticated_reads))

top_genera <- genus_summary |>
  slice_head(n = TOP_N) |>
  pull(genus)

pie_data <- genus_summary |>
  mutate(category = if_else(genus %in% top_genera, genus, "Other")) |>
  group_by(category) |>
  summarise(authenticated_reads = sum(authenticated_reads, na.rm = TRUE), .groups = "drop") |>
  mutate(
    category = factor(category, levels = c(top_genera, "Other")),
    proportion = authenticated_reads / sum(authenticated_reads),
    percent = scales::percent(proportion, accuracy = 0.1),
    slice_label = if_else(proportion >= 0.04, percent, "")
  ) |>
  arrange(category)

legend_text <- ifelse(
  as.character(pie_data$category) == "Other",
  paste0("'Other (", pie_data$percent, ")'"),
  paste0("italic('", as.character(pie_data$category), "')~' (", pie_data$percent, ")'")
)
legend_labels <- parse(text = legend_text)

p_pie <- ggplot(
  pie_data,
  aes(x = "", y = proportion, fill = category)
) +
  geom_col(width = 1, color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = slice_label),
    position = position_stack(vjust = 0.5),
    size = 5.4,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  scale_fill_brewer(
    palette = "Paired",
    breaks = levels(pie_data$category),
    labels = legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "Dominant authenticated plant genera by read representation",
    fill = NULL
  ) +
  theme_void(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 19, hjust = 0.5, margin = margin(b = 12)),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 18, r = 20, b = 15, l = 20)
  )

save_figure(
  p_pie,
  "Figure_06_dominant_genera",
  width = 11,
  height = 7.5
)
