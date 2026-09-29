# ==============================================================================
# FIGURE 5 — RAREFIED PLANT DIVERSITY
# ==============================================================================

source("00_config.R")
library(tidyverse)
library(vegan)
library(mgcv)

set.seed(SET_RANDOM_SEED)

community_wide <- readRDS(file.path(DERIVED_DIR, "community_counts_wide_primary.rds"))
samples <- readRDS(file.path(DERIVED_DIR, "sample_info.rds"))

# ==============================================================================
# RAREFACTION AND DIVERSITY
# ==============================================================================

mat <- as.matrix(
  community_wide[
    ,
    setdiff(names(community_wide), "library_id"),
    drop = FALSE
  ]
)

storage.mode(mat) <- "numeric"
rownames(mat) <- community_wide$library_id
mat <- round(mat)

sample_totals <- rowSums(mat, na.rm = TRUE)
positive_totals <- sample_totals[sample_totals > 0]

if (length(positive_totals) < 2) {
  stop("Too few non-zero libraries for rarefaction.")
}

rarefaction_depth <- floor(
  as.numeric(
    quantile(
      positive_totals,
      probs = 1 - RAREFACTION_RETAIN_SAMPLE_FRACTION,
      names = FALSE,
      type = 1
    )
  )
)
rarefaction_depth <- max(1L, rarefaction_depth)
eligible <- sample_totals >= rarefaction_depth

rare_out <- vegan::rarefy(
  mat[eligible, , drop = FALSE],
  sample = rarefaction_depth,
  se = TRUE
)

if (is.matrix(rare_out) && all(c("S", "se") %in% rownames(rare_out))) {
  rarefied_richness <- as.numeric(rare_out["S", , drop = TRUE])
  rarefied_se <- as.numeric(rare_out["se", , drop = TRUE])
  rarefied_library_ids <- colnames(rare_out)
} else {
  rarefied_richness <- as.numeric(rare_out)
  rarefied_se <- rep(NA_real_, length(rarefied_richness))
  rarefied_library_ids <- rownames(mat[eligible, , drop = FALSE])
}

eligible_mat <- mat[eligible, , drop = FALSE]

random_diversity <- map_dfr(
  seq_len(RAREFACTION_REPEATS),
  function(iteration) {
    rare_mat <- vegan::rrarefy(eligible_mat, sample = rarefaction_depth)
    tibble(
      iteration = iteration,
      library_id = rownames(rare_mat),
      shannon = vegan::diversity(rare_mat, index = "shannon"),
      inverse_simpson = vegan::diversity(rare_mat, index = "invsimpson")
    )
  }
)

random_diversity_summary <- random_diversity |>
  group_by(library_id) |>
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon),
    inverse_simpson_mean = mean(inverse_simpson),
    inverse_simpson_sd = sd(inverse_simpson),
    .groups = "drop"
  )

diversity_summary <- tibble(
  library_id = rownames(mat),
  total_retained_counts = sample_totals,
  observed_richness = rowSums(mat > 0, na.rm = TRUE),
  eligible_for_rarefaction = eligible
) |>
  left_join(
    tibble(
      library_id = rarefied_library_ids,
      expected_rarefied_richness = rarefied_richness,
      rarefied_richness_se = rarefied_se
    ),
    by = "library_id"
  ) |>
  left_join(random_diversity_summary, by = "library_id") |>
  left_join(samples, by = "library_id") |>
  arrange(master_depth_cm, library_id)

write_tsv(
  diversity_summary,
  file.path(TABLE_DIR, "Figure_05_diversity_rarefaction_summary.tsv")
)

# ==============================================================================
# GAM TRENDS AND FINAL FIGURE
# ==============================================================================

diversity_long <- diversity_summary |>
  filter(eligible_for_rarefaction) |>
  select(
    library_id,
    core,
    master_depth_cm,
    shannon_mean,
    inverse_simpson_mean
  ) |>
  pivot_longer(
    cols = c(shannon_mean, inverse_simpson_mean),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      inverse_simpson_mean = "Inverse Simpson diversity",
      shannon_mean = "Shannon diversity"
    ),
    metric = factor(
      metric,
      levels = c("Inverse Simpson diversity", "Shannon diversity")
    )
  ) |>
  filter(is.finite(master_depth_cm), is.finite(value))

fit_gam <- function(data) {
  model <- mgcv::gam(
    value ~ s(master_depth_cm, k = 6),
    data = data,
    method = "REML"
  )

  prediction_depths <- tibble(
    master_depth_cm = seq(
      min(data$master_depth_cm, na.rm = TRUE),
      max(data$master_depth_cm, na.rm = TRUE),
      length.out = 300
    )
  )

  predictions <- predict(model, newdata = prediction_depths, se.fit = TRUE)
  prediction_depths |>
    mutate(
      fit = as.numeric(predictions$fit),
      se = as.numeric(predictions$se.fit)
    )
}

gam_predictions <- diversity_long |>
  group_split(metric) |>
  map_dfr(
    function(df) {
      metric_name <- as.character(unique(df$metric))
      fit_gam(df) |>
        mutate(metric = metric_name)
    }
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("Inverse Simpson diversity", "Shannon diversity")
    )
  ) |>
  arrange(metric, master_depth_cm)

p_diversity <- ggplot(
  diversity_long,
  aes(x = value, y = master_depth_cm)
) +
  geom_path(
    data = gam_predictions,
    aes(x = fit, y = master_depth_cm, group = metric),
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_point(aes(shape = core), size = 2.3, alpha = 0.85) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  depth_scale_with_age() +
  labs(
    title = "Rarefied plant diversity through the sediment sequence",
    x = "Diversity index",
    y = "Composite master depth (cm)",
    shape = "Core"
  ) +
  thesis_theme()

save_figure(
  p_diversity,
  "Figure_05_plant_diversity",
  width = 11,
  height = 6
)
