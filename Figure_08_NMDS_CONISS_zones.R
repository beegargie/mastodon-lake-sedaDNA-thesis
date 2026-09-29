# ==============================================================================
# FIGURE 8 — NMDS WITH SIX CONISS ZONES AND ZONE-CHARACTERISTIC GENERA
# ==============================================================================

source("00_config.R")

library(tidyverse)
library(vegan)
library(rioja)
library(ggrepel)

set.seed(SET_RANDOM_SEED)

# ==============================================================================
# SETTINGS
# ==============================================================================

NMDS_GENERA_PER_ZONE <- 2L

NMDS_ZONE_GENUS_MIN_PREVALENCE <- max(
  2L,
  NMDS_PASSIVE_MIN_PREVALENCE
)

NMDS_ZONE_GENUS_MIN_WITHIN_ZONE_DETECTIONS <- 2L

NMDS_ZONE_HULL_ALPHA <- 0.14

selected_n_zones <- CONISS_N_ZONES

zone_levels <- paste0(
  "Zone ",
  LETTERS[seq_len(selected_n_zones)]
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

community_wide <- readRDS(
  file.path(
    DERIVED_DIR,
    "community_counts_wide_primary.rds"
  )
)

community_long <- readRDS(
  file.path(
    DERIVED_DIR,
    "community_long_primary.rds"
  )
)

samples <- readRDS(
  file.path(
    DERIVED_DIR,
    "sample_info.rds"
  )
)

# ==============================================================================
# BUILD COMMUNITY MATRIX FOR NMDS
# ==============================================================================

mat <- as.matrix(
  community_wide[
    ,
    setdiff(
      names(community_wide),
      "library_id"
    ),
    drop = FALSE
  ]
)

storage.mode(mat) <- "numeric"

rownames(mat) <- community_wide$library_id

prevalence <- colSums(
  mat > 0,
  na.rm = TRUE
)

mat <- mat[
  rowSums(
    mat,
    na.rm = TRUE
  ) > 0,
  
  prevalence >= MIN_SAMPLES_PER_GENUS_FOR_ORDINATION,
  
  drop = FALSE
]

meta <- samples |>
  filter(
    library_id %in% rownames(mat)
  ) |>
  arrange(
    match(
      library_id,
      rownames(mat)
    )
  )

if (
  nrow(mat) < 4 ||
  ncol(mat) < 2
) {
  stop(
    "Too few retained samples or genera for NMDS."
  )
}

# ==============================================================================
# HELLINGER TRANSFORMATION + BRAY-CURTIS NMDS
# ==============================================================================

hellinger <- vegan::decostand(
  mat,
  method = "hellinger"
)

nmds <- vegan::metaMDS(
  hellinger,
  distance = "bray",
  k = NMDS_DIMENSIONS,
  trymax = NMDS_TRYMAX,
  autotransform = FALSE,
  trace = FALSE
)

site_scores_matrix <- vegan::scores(
  nmds,
  display = "sites"
)

site_scores <- as.data.frame(
  site_scores_matrix
) |>
  rownames_to_column(
    "library_id"
  ) |>
  left_join(
    meta,
    by = "library_id"
  ) |>
  arrange(
    master_depth_cm,
    library_id
  )

# ==============================================================================
# RECONSTRUCT SIX-ZONE CONISS SOLUTION
# ==============================================================================

coniss_sample_order <- samples |>
  arrange(
    master_depth_cm,
    library_id
  )

coniss_mat <- as.matrix(
  community_wide[
    ,
    setdiff(
      names(community_wide),
      "library_id"
    ),
    drop = FALSE
  ]
)

storage.mode(coniss_mat) <- "numeric"

rownames(coniss_mat) <- community_wide$library_id

coniss_mat <- coniss_mat[
  coniss_sample_order$library_id,
  ,
  drop = FALSE
]

coniss_prevalence <- colSums(
  coniss_mat > 0,
  na.rm = TRUE
)

coniss_mat <- coniss_mat[
  rowSums(
    coniss_mat,
    na.rm = TRUE
  ) > 0,
  
  coniss_prevalence >= MIN_SAMPLES_PER_GENUS_FOR_ZONATION,
  
  drop = FALSE
]

coniss_sample_order <- coniss_sample_order |>
  filter(
    library_id %in% rownames(coniss_mat)
  ) |>
  arrange(
    match(
      library_id,
      rownames(coniss_mat)
    )
  )

coniss_hellinger <- vegan::decostand(
  coniss_mat,
  method = "hellinger"
)

coniss_distance <- stats::dist(
  coniss_hellinger,
  method = "euclidean"
)

coniss_local <- rioja::chclust(
  coniss_distance,
  method = "coniss"
)

zone_membership <- tibble(
  library_id = rownames(coniss_mat),
  
  zone_raw = stats::cutree(
    coniss_local,
    k = selected_n_zones
  )
) |>
  left_join(
    coniss_sample_order |>
      select(
        library_id,
        master_depth_cm
      ),
    by = "library_id"
  )

# ==============================================================================
# LABEL ZONES FROM DEEPEST TO SHALLOWEST
# ==============================================================================

zone_label_order <- zone_membership |>
  group_by(
    zone_raw
  ) |>
  summarise(
    median_depth_cm = median(
      master_depth_cm,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(
    desc(median_depth_cm)
  ) |>
  mutate(
    zone = paste0(
      "Zone ",
      LETTERS[row_number()]
    )
  )

zone_membership <- zone_membership |>
  left_join(
    zone_label_order,
    by = "zone_raw"
  ) |>
  mutate(
    age_calBP = age_from_depth(
      master_depth_cm
    ),
    
    zone = factor(
      zone,
      levels = zone_levels
    )
  ) |>
  select(
    library_id,
    zone,
    master_depth_cm,
    age_calBP
  )

zone_membership |>
  group_by(
    zone
  ) |>
  summarise(
    shallowest_depth_cm = min(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    deepest_depth_cm = max(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  print()

# ==============================================================================
# ATTACH CONISS ZONES TO NMDS SCORES
# ==============================================================================

site_scores_zoned <- site_scores |>
  inner_join(
    zone_membership |>
      select(
        library_id,
        zone
      ),
    by = "library_id"
  ) |>
  mutate(
    zone = factor(
      zone,
      levels = zone_levels
    )
  )

# ==============================================================================
# ZONE SUMMARIES
# ==============================================================================

zone_summary <- site_scores_zoned |>
  group_by(
    zone
  ) |>
  summarise(
    n_libraries = n(),
    
    minimum_depth_cm = min(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    maximum_depth_cm = max(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    median_depth_cm = median(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    median_age_calBP = age_from_depth(
      median_depth_cm
    ),
    
    centroid_NMDS1 = mean(
      NMDS1,
      na.rm = TRUE
    ),
    
    centroid_NMDS2 = mean(
      NMDS2,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# ==============================================================================
# ZONE COLORS
# ==============================================================================

zone_colors <- c(
  "Zone A" = "#0072B2",
  "Zone B" = "#6A3D9A",
  "Zone C" = "#009E73",
  "Zone D" = "#E69F00",
  "Zone E" = "#D55E00",
  "Zone F" = "#CC79A7"
)

# ==============================================================================
# CONVEX HULLS
# ==============================================================================

hull_data <- site_scores_zoned |>
  group_by(
    zone
  ) |>
  filter(
    n() >= 3
  ) |>
  slice(
    chull(
      NMDS1,
      NMDS2
    )
  ) |>
  ungroup()

# ==============================================================================
# TEMPORAL CENTROID TRAJECTORY
# ==============================================================================

zone_centroids <- zone_summary |>
  transmute(
    zone,
    NMDS1 = centroid_NMDS1,
    NMDS2 = centroid_NMDS2,
    median_depth_cm,
    median_age_calBP
  ) |>
  arrange(
    desc(median_depth_cm)
  )

oldest_centroid <- zone_centroids |>
  slice_head(
    n = 1
  )

youngest_centroid <- zone_centroids |>
  slice_tail(
    n = 1
  )

# ==============================================================================
# PASSIVE GENUS SCORES
# ==============================================================================

species_scores_matrix <- vegan::wascores(
  site_scores_matrix,
  hellinger
)

species_scores_all <- as.data.frame(
  species_scores_matrix
) |>
  rownames_to_column(
    "genus"
  ) |>
  mutate(
    ordination_strength = sqrt(
      NMDS1^2 +
        NMDS2^2
    )
  ) |>
  left_join(
    tibble(
      genus = colnames(mat),
      
      prevalence = colSums(
        mat > 0,
        na.rm = TRUE
      )
    ),
    by = "genus"
  ) |>
  left_join(
    community_long |>
      distinct(
        genus,
        dataset_total_abundance_reads,
        biogeo
      ),
    by = "genus"
  )

# ==============================================================================
# IDENTIFY ZONE-CHARACTERISTIC GENERA
# ==============================================================================

hellinger_long_zoned <- as.data.frame(
  hellinger
) |>
  rownames_to_column(
    "library_id"
  ) |>
  as_tibble() |>
  pivot_longer(
    -library_id,
    names_to = "genus",
    values_to = "hellinger_abundance"
  ) |>
  inner_join(
    site_scores_zoned |>
      select(
        library_id,
        zone
      ),
    by = "library_id"
  )

zone_genus_association <- map_dfr(
  zone_levels,
  
  function(one_zone) {
    
    hellinger_long_zoned |>
      group_by(
        genus
      ) |>
      summarise(
        
        zone_mean_hellinger = mean(
          hellinger_abundance[
            as.character(zone) == one_zone
          ],
          na.rm = TRUE
        ),
        
        outside_mean_hellinger = mean(
          hellinger_abundance[
            as.character(zone) != one_zone
          ],
          na.rm = TRUE
        ),
        
        within_zone_detections = sum(
          hellinger_abundance[
            as.character(zone) == one_zone
          ] > 0,
          na.rm = TRUE
        ),
        
        overall_prevalence = sum(
          hellinger_abundance > 0,
          na.rm = TRUE
        ),
        
        .groups = "drop"
      ) |>
      mutate(
        zone = one_zone,
        
        zone_association =
          zone_mean_hellinger -
          outside_mean_hellinger
      )
  }
) |>
  mutate(
    zone = factor(
      zone,
      levels = zone_levels
    )
  ) |>
  left_join(
    species_scores_all |>
      select(
        genus,
        dataset_total_abundance_reads,
        biogeo
      ),
    by = "genus"
  )

best_zone_per_genus <- zone_genus_association |>
  filter(
    overall_prevalence >=
      NMDS_ZONE_GENUS_MIN_PREVALENCE,
    
    within_zone_detections >=
      NMDS_ZONE_GENUS_MIN_WITHIN_ZONE_DETECTIONS,
    
    is.finite(
      zone_association
    )
  ) |>
  group_by(
    genus
  ) |>
  slice_max(
    order_by = zone_association,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  filter(
    zone_association > 0
  )

selected_zone_genera <- best_zone_per_genus |>
  group_by(
    zone
  ) |>
  arrange(
    desc(zone_association),
    desc(dataset_total_abundance_reads),
    .by_group = TRUE
  ) |>
  slice_head(
    n = NMDS_GENERA_PER_ZONE
  ) |>
  ungroup() |>
  left_join(
    species_scores_all |>
      select(
        genus,
        NMDS1,
        NMDS2,
        ordination_strength,
        prevalence
      ),
    by = "genus"
  ) |>
  arrange(
    zone,
    desc(zone_association)
  )

# ==============================================================================
# PLOT LIMITS
# ==============================================================================

x_values <- c(
  site_scores_zoned$NMDS1,
  selected_zone_genera$NMDS1
)

y_values <- c(
  site_scores_zoned$NMDS2,
  selected_zone_genera$NMDS2
)

x_pad <- 0.12 *
  diff(
    range(
      x_values,
      na.rm = TRUE
    )
  )

y_pad <- 0.12 *
  diff(
    range(
      y_values,
      na.rm = TRUE
    )
  )

if (
  !is.finite(x_pad) ||
  x_pad <= 0
) {
  x_pad <- 0.25
}

if (
  !is.finite(y_pad) ||
  y_pad <= 0
) {
  y_pad <- 0.25
}

# ==============================================================================
# FINAL NMDS FIGURE
# ==============================================================================

p_nmds_taxa <- ggplot() +
  
  geom_polygon(
    data = hull_data,
    aes(
      x = NMDS1,
      y = NMDS2,
      group = zone,
      fill = zone
    ),
    alpha = NMDS_ZONE_HULL_ALPHA,
    color = NA
  ) +
  
  geom_path(
    data = zone_centroids,
    aes(
      x = NMDS1,
      y = NMDS2,
      group = 1
    ),
    color = "grey35",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = zone_centroids,
    aes(
      x = NMDS1,
      y = NMDS2
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    stroke = 0.8,
    fill = "white",
    color = "grey35"
  ) +
  
  geom_text(
    data = oldest_centroid,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = "Oldest"
    ),
    inherit.aes = FALSE,
    nudge_x = 0.08,
    nudge_y = -0.08,
    size = 7,
    fontface = "bold",
    color = "grey25"
  ) +
  
  geom_text(
    data = youngest_centroid,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = "Youngest"
    ),
    inherit.aes = FALSE,
    nudge_x = 0.08,
    nudge_y = 0.08,
    size = 7,
    fontface = "bold",
    color = "grey25"
  ) +
  
  geom_point(
    data = site_scores_zoned,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = zone,
      shape = core
    ),
    size = 4,
    alpha = 0.95
  ) +
  
  ggrepel::geom_text_repel(
    data = selected_zone_genera,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = genus
    ),
    inherit.aes = FALSE,
    size = 5.5,
    fontface = "italic",
    color = "grey15",
    box.padding = 0.4,
    point.padding = 0.1,
    min.segment.length = 0.05,
    segment.color = "grey70",
    segment.size = 0.3,
    max.overlaps = Inf,
    seed = SET_RANDOM_SEED
  ) +
  
  scale_color_manual(
    values = zone_colors,
    breaks = zone_levels,
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = zone_colors,
    breaks = zone_levels,
    drop = FALSE,
    guide = "none"
  ) +
  
  coord_cartesian(
    xlim =
      range(
        x_values,
        na.rm = TRUE
      ) +
      c(
        -x_pad,
        x_pad
      ),
    
    ylim =
      range(
        y_values,
        na.rm = TRUE
      ) +
      c(
        -y_pad,
        y_pad
      ),
    
    clip = "off"
  ) +
  
  labs(
    title = "Plant-community NMDS with zone-characteristic genera",
    
    subtitle = paste0(
      nrow(selected_zone_genera),
      " genera selected by enrichment within the six CONISS zones; stress = ",
      round(
        nmds$stress,
        3
      )
    ),
    
    x = "NMDS1",
    y = "NMDS2",
    color = "CONISS zone",
    shape = "Core"
  ) +
  
  thesis_theme(
    base_size = 20
  ) +
  
  theme(
    
    
    plot.title = element_text(
      size = 25,
      face = "bold",
      margin = margin(
        b = 5
      )
    ),
    
    
    plot.subtitle = element_text(
      size = 19,
      margin = margin(
        b = 12
      )
    ),
    
    
    axis.title.x = element_text(
      size = 20,
      margin = margin(
        t = 10
      )
    ),
    
    axis.title.y = element_text(
      size = 20,
      margin = margin(
        r = 10
      )
    ),
    
    
    axis.text.x = element_text(
      size = 18
    ),
    
    axis.text.y = element_text(
      size = 18
    ),
    
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 17
    ),
    
    legend.key.height = grid::unit(
      0.85,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.85,
      "cm"
    ),
    
    legend.spacing.y = grid::unit(
      0.18,
      "cm"
    ),
    
    
    plot.margin = margin(
      t = 12,
      r = 18,
      b = 14,
      l = 18
    )
  )

# ==============================================================================
# ==============================================================================

# ==============================================================================
# SAVE SUPPORTING TABLES
# ==============================================================================

write_tsv(
  zone_membership |> arrange(desc(master_depth_cm)),
  file.path(TABLE_DIR, "Figure_08_CONISS_zone_membership.tsv")
)

write_tsv(
  zone_summary,
  file.path(TABLE_DIR, "Figure_08_CONISS_zone_summary.tsv")
)

write_tsv(
  zone_genus_association,
  file.path(TABLE_DIR, "Figure_08_zone_genus_associations.tsv")
)

# ==============================================================================
# SAVE
# ==============================================================================

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Figure_08_NMDS_CONISS_zones.png"
  ),
  plot = p_nmds_taxa,
  width = 14,
  height = 10.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Figure_08_NMDS_CONISS_zones.pdf"
  ),
  plot = p_nmds_taxa,
  width = 14,
  height = 10.2,
  bg = "white"
)
