# ==============================================================================
# FIGURE 10 — AUTHENTICATED VERTEBRATE-FAMILY DETECTIONS
# ==============================================================================

source("00_config.R")
library(tidyverse)
library(patchwork)

STRICT_THRESHOLD <- 100L
PERMISSIVE_THRESHOLD <- 10L
PERMISSIVE_TOP_N <- 30L
DETECTION_COLOR <- "#1F4E79"

samples <- readRDS(file.path(DERIVED_DIR, "sample_info.rds")) |>
  transmute(
    library_id = as.character(library_id),
    core,
    master_depth_cm = as.numeric(master_depth_cm)
  ) |>
  distinct(library_id, .keep_all = TRUE) |>
  arrange(master_depth_cm, library_id) |>
  mutate(age_calBP = age_from_depth(master_depth_cm))

family_meta <- readRDS(FAMILY_DATA_PATH)

required_family_columns <- c(
  "library_id",
  "family",
  "rank",
  "A",
  "Zfit",
  "total_n_reads_from_mD",
  "n_alns"
)
missing_family_columns <- setdiff(required_family_columns, names(family_meta))
if (length(missing_family_columns) > 0) {
  stop("Family-level input is missing required columns: ", paste(missing_family_columns, collapse = ", "))
}

coalesce_columns <- function(data, candidates) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0) {
    return(rep("", nrow(data)))
  }

  out <- as.character(data[[found[[1]]]])
  if (length(found) > 1) {
    for (column_name in found[-1]) {
      out <- dplyr::coalesce(out, as.character(data[[column_name]]))
    }
  }
  dplyr::coalesce(out, "")
}

phylum_value <- coalesce_columns(family_meta, c("phylum.x", "phylum.y", "phylum"))
class_value <- coalesce_columns(family_meta, c("class.x", "class.y", "class"))

family_meta <- family_meta |>
  mutate(
    library_id = as.character(library_id),
    family = as.character(family),
    A_value = suppressWarnings(as.numeric(A)),
    Zfit_value = suppressWarnings(as.numeric(Zfit)),
    abundance_reads = suppressWarnings(as.numeric(total_n_reads_from_mD)),
    support_n_alns = suppressWarnings(as.numeric(n_alns)),
    phylum_value = phylum_value,
    class_value = class_value,
    is_vertebrate =
      str_detect(str_to_lower(phylum_value), "chordata|vertebrata") |
      str_detect(
        str_to_lower(class_value),
        paste(
          c("mammalia", "aves", "actinopter", "chondrichthyes", "amphibia", "reptilia"),
          collapse = "|"
        )
      )
  ) |>
  filter(
    library_id %in% samples$library_id,
    rank == "family",
    is_vertebrate,
    is.finite(A_value),
    is.finite(Zfit_value),
    A_value >= DAMAGE_AMPLITUDE_MIN,
    Zfit_value >= ZFIT_MIN,
    !is.na(family),
    family != "",
    is.finite(abundance_reads),
    abundance_reads > 0
  )

vertebrate_family_all <- family_meta |>
  group_by(library_id, family) |>
  summarise(
    vertebrate_abundance_reads = sum(abundance_reads, na.rm = TRUE),
    vertebrate_support_n_alns = sum(support_n_alns, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    samples |> select(library_id, core, master_depth_cm, age_calBP),
    by = "library_id"
  ) |>
  arrange(master_depth_cm, family)

permissive_detections <- vertebrate_family_all |>
  filter(vertebrate_abundance_reads >= PERMISSIVE_THRESHOLD)

strict_detections <- vertebrate_family_all |>
  filter(vertebrate_abundance_reads >= STRICT_THRESHOLD)

write_tsv(
  permissive_detections,
  file.path(TABLE_DIR, "08_vertebrate_family_permissive_detections.tsv")
)
write_tsv(
  strict_detections,
  file.path(TABLE_DIR, "08_vertebrate_family_strict_detections.tsv")
)

# ==============================================================================
# FAMILY SUMMARIES AND DISPLAY ORDER
# ==============================================================================

summarize_families <- function(detections) {
  detections |>
    group_by(family) |>
    summarise(
      retained_authenticated_reads = sum(vertebrate_abundance_reads, na.rm = TRUE),
      n_libraries = n_distinct(library_id),
      shallowest_detection_depth_cm = min(master_depth_cm, na.rm = TRUE),
      deepest_detection_depth_cm = max(master_depth_cm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      shallowest_detection_age_calBP = age_from_depth(shallowest_detection_depth_cm),
      deepest_detection_age_calBP = age_from_depth(deepest_detection_depth_cm)
    )
}

strict_family_summary <- summarize_families(strict_detections)
permissive_family_summary <- summarize_families(permissive_detections)

strict_displayed <- strict_family_summary
permissive_displayed <- permissive_family_summary |>
  arrange(desc(retained_authenticated_reads), desc(n_libraries)) |>
  slice_head(n = PERMISSIVE_TOP_N)

strict_family_order <- strict_displayed |>
  arrange(
    desc(shallowest_detection_depth_cm),
    desc(n_libraries),
    desc(retained_authenticated_reads)
  ) |>
  pull(family)

permissive_family_order <- permissive_displayed |>
  arrange(
    desc(shallowest_detection_depth_cm),
    desc(n_libraries),
    desc(retained_authenticated_reads)
  ) |>
  pull(family)

write_tsv(
  strict_family_summary |> arrange(desc(retained_authenticated_reads)),
  file.path(TABLE_DIR, "Figure_10_strict_family_summary.tsv")
)
write_tsv(
  permissive_family_summary |> arrange(desc(retained_authenticated_reads)),
  file.path(TABLE_DIR, "Figure_10_permissive_family_summary.tsv")
)

sample_plot_positions <- samples |>
  arrange(master_depth_cm, library_id) |>
  mutate(
    previous_depth = lag(master_depth_cm),
    next_depth = lead(master_depth_cm),
    previous_gap = master_depth_cm - previous_depth,
    next_gap = next_depth - master_depth_cm,
    local_gap = case_when(
      is.na(previous_gap) & is.na(next_gap) ~ 5,
      is.na(previous_gap) ~ next_gap,
      is.na(next_gap) ~ previous_gap,
      TRUE ~ pmin(previous_gap, next_gap)
    ),
    tile_height = pmin(pmax(local_gap * 0.80, 1.5), 7)
  ) |>
  select(library_id, core, master_depth_cm, age_calBP, tile_height)

make_heatmap_data <- function(detections, family_order) {
  expand_grid(
    library_id = sample_plot_positions$library_id,
    family = family_order
  ) |>
    left_join(sample_plot_positions, by = "library_id") |>
    left_join(
      detections |>
        filter(family %in% family_order) |>
        transmute(library_id, family, retained = 1L),
      by = c("library_id", "family")
    ) |>
    mutate(
      retained = replace_na(retained, 0L),
      family = factor(family, levels = family_order)
    )
}

strict_heat <- make_heatmap_data(strict_detections, strict_family_order)
permissive_heat <- make_heatmap_data(permissive_detections, permissive_family_order)

make_detection_panel <- function(data, threshold, title, x_text_size) {
  ggplot(
    data,
    aes(
      x = family,
      y = master_depth_cm,
      height = tile_height,
      fill = factor(retained, levels = c(0, 1))
    )
  ) +
    geom_tile(width = 0.94, color = "grey90", linewidth = 0.08) +
    scale_fill_manual(
      values = c("0" = "grey97", "1" = DETECTION_COLOR),
      labels = c("0" = "Not retained", "1" = paste0("≥", threshold, " reads")),
      drop = FALSE
    ) +
    scale_x_discrete(drop = FALSE, guide = guide_axis(angle = 60, n.dodge = 1)) +
    depth_scale_with_age() +
    labs(
      title = title,
      x = "Vertebrate family",
      y = "Composite master depth (cm)",
      fill = "Detection"
    ) +
    thesis_theme(base_size = 22) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 28, face = "bold", margin = margin(b = 12)),
      axis.title.x = element_text(size = 24, margin = margin(t = 14)),
      axis.title.y = element_text(size = 24, margin = margin(r = 12)),
      axis.title.y.right = element_text(size = 24, margin = margin(l = 12)),
      axis.text.x = element_text(size = x_text_size, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 18),
      axis.text.y.right = element_text(size = 18),
      legend.position = "right",
      legend.title = element_text(size = 21, face = "bold"),
      legend.text = element_text(size = 20),
      legend.key.height = grid::unit(1.0, "cm"),
      legend.key.width = grid::unit(1.0, "cm"),
      plot.margin = margin(t = 14, r = 24, b = 80, l = 14)
    )
}

p_strict <- make_detection_panel(
  strict_heat,
  STRICT_THRESHOLD,
  "A. Strict authenticated vertebrate-family detections",
  24
)

p_permissive <- make_detection_panel(
  permissive_heat,
  PERMISSIVE_THRESHOLD,
  "B. Permissive authenticated vertebrate-family detections",
  20
)

combined_vertebrate_plot <- p_strict / p_permissive +
  patchwork::plot_layout(heights = c(1, 1))

save_figure(
  combined_vertebrate_plot,
  "Figure_10_vertebrate_detections",
  width = 20,
  height = 18
)

# ==============================================================================
# SEQUENCING-DEPTH ANALYSIS REPORTED IN THE RESULTS
# ==============================================================================

qc <- readr::read_tsv(
  QC_TABLE_PATH,
  show_col_types = FALSE,
  na = c("", "NA")
) |>
  filter(sample_class == "lake_sample") |>
  transmute(
    library_id = as.character(library_id),
    reads_sequenced = as.numeric(reads_sequenced)
  ) |>
  distinct(library_id, .keep_all = TRUE)

make_library_recovery <- function(detections) {
  recovery <- detections |>
    group_by(library_id) |>
    summarise(
      total_retained_vertebrate_family_reads = sum(vertebrate_abundance_reads, na.rm = TRUE),
      n_retained_vertebrate_families = n_distinct(family),
      .groups = "drop"
    )

  samples |>
    select(library_id, core, master_depth_cm) |>
    left_join(qc, by = "library_id") |>
    left_join(recovery, by = "library_id") |>
    mutate(
      total_retained_vertebrate_family_reads = replace_na(
        total_retained_vertebrate_family_reads,
        0
      ),
      n_retained_vertebrate_families = replace_na(
        n_retained_vertebrate_families,
        0L
      )
    )
}

correlation_summary <- function(data, criterion) {
  cor_reads <- cor.test(
    data$reads_sequenced,
    data$total_retained_vertebrate_family_reads,
    method = "spearman",
    exact = FALSE,
    use = "complete.obs"
  )

  cor_families <- cor.test(
    data$reads_sequenced,
    data$n_retained_vertebrate_families,
    method = "spearman",
    exact = FALSE,
    use = "complete.obs"
  )

  tibble(
    criterion = criterion,
    response = c(
      "Total retained vertebrate-family reads",
      "Number of retained vertebrate families"
    ),
    spearman_rho = c(unname(cor_reads$estimate), unname(cor_families$estimate)),
    p_value = c(cor_reads$p.value, cor_families$p.value),
    n_libraries = c(
      sum(complete.cases(data$reads_sequenced, data$total_retained_vertebrate_family_reads)),
      sum(complete.cases(data$reads_sequenced, data$n_retained_vertebrate_families))
    )
  )
}

permissive_recovery <- make_library_recovery(permissive_detections)
strict_recovery <- make_library_recovery(strict_detections)

sequencing_depth_correlations <- bind_rows(
  correlation_summary(permissive_recovery, "Permissive (>=10 reads)"),
  correlation_summary(strict_recovery, "Strict (>=100 reads)")
)

write_tsv(
  sequencing_depth_correlations,
  file.path(TABLE_DIR, "Figure_10_sequencing_depth_correlations.tsv")
)
