# ==============================================================================
# FIGURE 4 — LIBRARY SEQUENCING DEPTH AND EUKARYOTIC MAPPING
# ==============================================================================

source("00_config.R")
library(tidyverse)

qc <- readr::read_tsv(
  QC_TABLE_PATH,
  show_col_types = FALSE,
  na = c("", "NA")
)

samples <- readRDS(file.path(DERIVED_DIR, "sample_info.rds"))

required_qc <- c("library_id", "sample_class", "reads_sequenced", "pct_mapped_euk")
missing_qc <- setdiff(required_qc, names(qc))
if (length(missing_qc) > 0) {
  stop("QC table is missing required columns: ", paste(missing_qc, collapse = ", "))
}

lake <- qc |>
  filter(sample_class == "lake_sample") |>
  transmute(
    library_id = as.character(library_id),
    reads_sequenced = as.numeric(reads_sequenced),
    pct_mapped_euk = as.numeric(pct_mapped_euk)
  ) |>
  left_join(
    samples |> select(library_id, master_depth_cm),
    by = "library_id"
  ) |>
  mutate(reads_sequenced_million = reads_sequenced / 1e6) |>
  arrange(master_depth_cm)

if (nrow(lake) != 92) {
  warning("Expected 92 lake libraries but found ", nrow(lake), ".")
}

if (any(is.na(lake$master_depth_cm))) {
  stop("One or more lake libraries are missing composite master depth.")
}

median_reads_million <- median(lake$reads_sequenced_million, na.rm = TRUE)
median_pct_mapped_euk <- median(lake$pct_mapped_euk, na.rm = TRUE)

plot_long <- lake |>
  transmute(
    master_depth_cm,
    `A. Sequencing depth (million reads)` = reads_sequenced_million,
    `B. Eukaryotic mapping (%)` = pct_mapped_euk
  ) |>
  pivot_longer(
    cols = -master_depth_cm,
    names_to = "panel",
    values_to = "value"
  )

reference_lines <- tibble(
  panel = c(
    "A. Sequencing depth (million reads)",
    "B. Eukaryotic mapping (%)"
  ),
  yintercept = c(median_reads_million, median_pct_mapped_euk)
)

median_keys <- plot_long |>
  group_by(panel) |>
  summarise(panel_max = max(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    key_y = panel_max * 1.10,
    label = case_when(
      panel == "A. Sequencing depth (million reads)" ~
        paste0("Median = ", round(median_reads_million, 1), " million reads"),
      panel == "B. Eukaryotic mapping (%)" ~
        paste0("Median = ", round(median_pct_mapped_euk, 1), "%")
    )
  )

p_qc <- ggplot(plot_long, aes(x = master_depth_cm, y = value)) +
  geom_hline(
    data = reference_lines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_point(size = 2.3, alpha = 0.85) +
  geom_segment(
    data = median_keys,
    aes(x = 245, xend = 295, y = key_y, yend = key_y),
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = median_keys,
    aes(x = 305, y = key_y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 4.5,
    inherit.aes = FALSE
  ) +
  facet_wrap(~panel, scales = "free_y", ncol = 1) +
  labs(
    title = "Library-level QC across the Mastodon Lake sediment sequence",
    x = "Composite master depth (cm)",
    y = NULL
  ) +
  scale_x_continuous(
    breaks = seq(0, 650, by = 100),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 17)
  )

save_figure(
  p_qc,
  "Figure_04_library_QC",
  width = 9,
  height = 8
)
