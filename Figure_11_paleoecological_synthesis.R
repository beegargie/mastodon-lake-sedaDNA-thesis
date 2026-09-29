# ==============================================================================
# FIGURE 11 — PALEOECOLOGICAL SYNTHESIS
# ==============================================================================

source("00_config.R")

library(tidyverse)
library(vegan)
library(rioja)

N_ZONES <- CONISS_N_ZONES
MAX_GROUPS_TO_TEST <- CONISS_MAX_GROUPS_TO_TEST

# ==============================================================================
# 1. LOAD DATA
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

if (!exists("age_from_depth")) {
  stop("age_from_depth() must be defined in 00_config.R.")
}

# ==============================================================================
# 2. PREPARE COMMUNITY MATRIX AND FIT CONISS
# ==============================================================================

sample_order <- samples |>
  arrange(
    master_depth_cm,
    library_id
  )

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

mat <- mat[
  sample_order$library_id,
  ,
  drop = FALSE
]

nonempty_samples <- rowSums(
  mat,
  na.rm = TRUE
) > 0

mat <- mat[
  nonempty_samples,
  ,
  drop = FALSE
]

sample_order <- sample_order |>
  filter(
    library_id %in% rownames(mat)
  ) |>
  arrange(
    match(
      library_id,
      rownames(mat)
    )
  )

prevalence <- colSums(
  mat > 0,
  na.rm = TRUE
)

mat <- mat[
  ,
  prevalence >= MIN_SAMPLES_PER_GENUS_FOR_ZONATION,
  drop = FALSE
]

hellinger <- vegan::decostand(
  mat,
  method = "hellinger"
)

hellinger_distance <- stats::dist(
  hellinger,
  method = "euclidean"
)

coniss <- rioja::chclust(
  hellinger_distance,
  method = "coniss"
)

# ==============================================================================
# 2A. SIX-ZONE MEMBERSHIP
# ==============================================================================

zone_membership <- tibble(
  library_id = rownames(mat),
  
  zone_raw = stats::cutree(
    coniss,
    k = N_ZONES
  )
) |>
  
  left_join(
    sample_order |>
      select(
        library_id,
        master_depth_cm
      ),
    by = "library_id"
  ) |>
  
  arrange(
    master_depth_cm,
    library_id
  )

zone_key <- zone_membership |>
  
  group_by(
    zone_raw
  ) |>
  
  summarise(
    median_depth_cm = median(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    minimum_depth_cm = min(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    maximum_depth_cm = max(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  
  arrange(
    desc(
      median_depth_cm
    )
  ) |>
  
  mutate(
    zone = paste0(
      "Zone ",
      LETTERS[row_number()]
    )
  )

zone_membership <- zone_membership |>
  
  left_join(
    zone_key |>
      select(
        zone_raw,
        zone
      ),
    by = "zone_raw"
  )

zone_ranges <- zone_membership |>
  
  group_by(
    zone
  ) |>
  
  summarise(
    min_depth_cm = min(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    max_depth_cm = max(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  
  arrange(
    min_depth_cm
  )

zone_boundaries <- tibble(
  boundary_depth_cm = (
    zone_ranges$max_depth_cm[
      -nrow(zone_ranges)
    ] +
      zone_ranges$min_depth_cm[
        -1
      ]
  ) / 2
)

# ==============================================================================
# 3. AGE-AXIS HELPER
# ==============================================================================

make_age_axis <- function(
    depth_values
) {
  
  depth_range <- range(
    depth_values,
    na.rm = TRUE
  )
  
  axis_depths <- unique(
    c(
      depth_range[[1]],
      200,
      400,
      600
    )
  )
  
  axis_depths <- axis_depths[
    axis_depths >= depth_range[[1]] &
      axis_depths <= depth_range[[2]]
  ]
  
  tibble(
    depth_cm = axis_depths,
    
    age_ka = age_from_depth(
      axis_depths
    ) / 1000
  )
}

# ==============================================================================
# 6. PALEOECOLOGICAL SYNTHESIS
# ==============================================================================

# ==============================================================================
# 6A. TAXA RETAINED FOR THE SYNTHESIS FIGURE
# ==============================================================================

SYNTHESIS_GENERA <- c(
  
  "Picea",
  "Larix",
  "Betula",
  "Populus",
  "Alnus",
  
  "Salix",
  "Dryas",
  
  "Potamogeton",
  "Ceratophyllum",
  "Stuckenia",
  "Utricularia",
  
  "Carex",
  "Schoenoplectus"
)

available_genera <- intersect(
  SYNTHESIS_GENERA,
  
  unique(
    as.character(
      community_long$genus
    )
  )
)

missing_selected_genera <- setdiff(
  SYNTHESIS_GENERA,
  available_genera
)

if (
  length(
    missing_selected_genera
  ) > 0
) {
  
  warning(
    "Selected synthesis genera not found in community_long: ",
    paste(
      missing_selected_genera,
      collapse = ", "
    )
  )
}

# ==============================================================================
# 6B. IDENTIFY AND AGGREGATE FERNS
# ==============================================================================

fern_rows <- community_long |>
  
  mutate(
    plot_group_character =
      as.character(
        plot_group
      ),
    
    plot_group_lower =
      str_to_lower(
        plot_group_character
      )
  ) |>
  
  filter(
    str_detect(
      plot_group_lower,
      "fern"
    )
  )

if (
  nrow(
    fern_rows
  ) == 0
) {
  
  stop(
    "Could not identify fern genera from plot_group."
  )
}

fern_group_name <- fern_rows |>
  
  filter(
    !is.na(
      plot_group_character
    )
  ) |>
  
  distinct(
    plot_group_character
  ) |>
  
  slice_head(
    n = 1
  ) |>
  
  pull(
    plot_group_character
  )

# ==============================================================================
# 6C. FUNCTIONAL-GROUP METADATA
# ==============================================================================

selected_meta <- community_long |>
  
  filter(
    genus %in%
      available_genera
  ) |>
  
  group_by(
    genus
  ) |>
  
  summarise(
    plot_group = {
      
      x <- as.character(
        na.omit(
          plot_group
        )
      )
      
      if (
        length(
          x
        ) > 0
      ) {
        x[[1]]
      } else {
        NA_character_
      }
    },
    
    .groups = "drop"
  )

selected_meta <- bind_rows(
  selected_meta,
  
  tibble(
    genus = "Ferns_total",
    plot_group =
      fern_group_name
  )
)

# ==============================================================================
# 6D. ORDER PANELS BY FUNCTIONAL GROUP
# ==============================================================================

synthesis_order <- c(
  
  "Picea",
  "Larix",
  "Betula",
  "Populus",
  "Alnus",
  
  "Salix",
  "Dryas",
  
  "Ferns_total",
  
  "Potamogeton",
  "Ceratophyllum",
  "Stuckenia",
  "Utricularia",
  
  "Carex",
  "Schoenoplectus"
)

synthesis_order <- synthesis_order[
  synthesis_order %in%
    selected_meta$genus
]

synthesis_display_names <- ifelse(
  synthesis_order ==
    "Ferns_total",
  
  "Ferns (total)",
  
  synthesis_order
)

# ==============================================================================
# 6E. BUILD SELECTED-GENUS DATA
# ==============================================================================

selected_genus_long <- community_long |>
  
  filter(
    genus %in%
      available_genera,
    
    library_id %in%
      rownames(
        mat
      )
  ) |>
  
  transmute(
    library_id,
    
    genus =
      as.character(
        genus
      ),
    
    relative_abundance
  )

fern_long <- fern_rows |>
  
  filter(
    library_id %in%
      rownames(
        mat
      )
  ) |>
  
  group_by(
    library_id
  ) |>
  
  summarise(
    relative_abundance = sum(
      relative_abundance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  
  mutate(
    genus = "Ferns_total"
  ) |>
  
  select(
    library_id,
    genus,
    relative_abundance
  )

synthesis_long <- bind_rows(
  selected_genus_long,
  fern_long
)

# ==============================================================================
# 6F. BUILD STRATIGRAPHIC MATRIX
# ==============================================================================

strat_data <- synthesis_long |>
  
  mutate(
    genus = factor(
      genus,
      levels =
        synthesis_order
    )
  ) |>
  
  pivot_wider(
    names_from =
      genus,
    
    values_from =
      relative_abundance,
    
    values_fill = 0
  ) |>
  
  right_join(
    sample_order |>
      select(
        library_id,
        master_depth_cm
      ),
    
    by = "library_id"
  ) |>
  
  arrange(
    match(
      library_id,
      rownames(
        mat
      )
    )
  )

for (
  taxon_name in synthesis_order
) {
  
  if (
    !taxon_name %in%
    names(
      strat_data
    )
  ) {
    
    strat_data[[taxon_name]] <- 0
  }
}

strat_depth <-
  strat_data$master_depth_cm

strat_matrix <- strat_data |>
  
  select(
    all_of(
      synthesis_order
    )
  ) |>
  
  as.data.frame()

strat_matrix <-
  100 * strat_matrix

# ==============================================================================
# 6G. PANEL-SPECIFIC X-AXIS LIMITS
# ==============================================================================

curve_maxima <- vapply(
  strat_matrix,
  
  function(
    x
  ) {
    
    max(
      x,
      na.rm = TRUE
    )
  },
  
  numeric(
    1
  )
)

choose_panel_max <- function(
    x
) {
  
  if (
    !is.finite(x) ||
    x <= 0
  ) {
    
    return(2)
  }
  
  
  if (x <= 2) {
    
    return(2)
    
  } else if (x <= 5) {
    
    return(5)
    
  } else if (x <= 10) {
    
    return(10)
    
  } else if (x <= 20) {
    
    return(20)
    
  } else if (x <= 40) {
    
    return(40)
    
  } else if (x <= 60) {
    
    return(60)
    
  } else if (x <= 80) {
    
    return(80)
    
  } else {
    
    return(100)
  }
}

panel_maxima <- vapply(
  curve_maxima,
  choose_panel_max,
  numeric(
    1
  )
)

panel_limits <- cbind(
  minimum = rep(
    0,
    length(
      panel_maxima
    )
  ),
  
  maximum =
    panel_maxima
)

rownames(
  panel_limits
) <- names(
  panel_maxima
)

# ==============================================================================
# 6H. FUNCTIONAL-GROUP COLORS
# ==============================================================================

group_lookup <- setNames(
  selected_meta$plot_group,
  selected_meta$genus
)

selected_groups <-
  group_lookup[
    synthesis_order
  ]

selected_colors <- unname(
  FUNCTIONAL_GROUP_COLORS[
    selected_groups
  ]
)

selected_colors[
  is.na(
    selected_colors
  )
] <- "grey55"

# ==============================================================================
# 6I. VERTEBRATE COLORS
# ==============================================================================

SEDA_COLORS <- c(
  "Equidae" =
    "#E69F00",
  
  "Elephantidae" =
    "#1F4E79"
)

FOSSIL_COLORS <- c(
  "Mastodon" =
    "#2E8B57",
  
  "Mammoth" =
    "#7A5195"
)

# ==============================================================================
# 6J. AGE AXIS
# ==============================================================================

strat_age_axis <- make_age_axis(
  strat_depth
)

# ==============================================================================
# 6K. ZONE RANGES
# ==============================================================================

zone_key_synthesis <- zone_membership |>
  
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
    desc(
      median_depth_cm
    )
  ) |>
  
  mutate(
    zone = paste0(
      "Zone ",
      LETTERS[
        row_number()
      ]
    )
  )

zone_membership_synthesis <-
  zone_membership |>
  
  select(
    -any_of(
      "zone"
    )
  ) |>
  
  left_join(
    zone_key_synthesis |>
      select(
        zone_raw,
        zone
      ),
    
    by = "zone_raw"
  )

zone_ranges_plot <-
  zone_membership_synthesis |>
  
  group_by(
    zone
  ) |>
  
  summarise(
    min_depth_cm = min(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    max_depth_cm = max(
      master_depth_cm,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  
  arrange(
    min_depth_cm
  ) |>
  
  mutate(
    midpoint_depth_cm = (
      min_depth_cm +
        max_depth_cm
    ) / 2
  )

zone_boundaries_synthesis <- tibble(
  boundary_depth_cm = (
    zone_ranges_plot$max_depth_cm[
      -nrow(
        zone_ranges_plot
      )
    ] +
      zone_ranges_plot$min_depth_cm[
        -1
      ]
  ) / 2
)

# ==============================================================================
# 6L. INVERT WORKING AGE MODEL
# ==============================================================================

depth_from_age <- function(
    target_age_calBP
) {
  
  lower_depth <- min(
    strat_depth,
    na.rm = TRUE
  )
  
  
  upper_depth <- max(
    strat_depth,
    na.rm = TRUE
  )
  
  
  vapply(
    target_age_calBP,
    
    function(
    target_age
    ) {
      
      stats::optimize(
        function(
    depth_cm
        ) {
          
          abs(
            age_from_depth(
              depth_cm
            ) -
              target_age
          )
        },
    
    interval = c(
      lower_depth,
      upper_depth
    )
    
      )$minimum
    },
    
    numeric(
      1
    )
  )
}

# ==============================================================================
# 6M. CLIMATE TRANSITIONS
# ==============================================================================

climate_events <- tibble(
  event = c(
    "Holocene onset",
    "YD onset",
    "B-A onset"
  ),
  
  age_calBP = c(
    11700,
    12900,
    14700
  )
) |>
  
  mutate(
    depth_cm =
      depth_from_age(
        age_calBP
      )
  )

climate_intervals <- tibble(
  climate = c(
    "Holocene",
    "YD",
    "B-A",
    "Pre-B-A"
  ),
  
  upper_depth = c(
    min(
      strat_depth,
      na.rm = TRUE
    ),
    
    climate_events$depth_cm[
      climate_events$event ==
        "Holocene onset"
    ],
    
    climate_events$depth_cm[
      climate_events$event ==
        "YD onset"
    ],
    
    climate_events$depth_cm[
      climate_events$event ==
        "B-A onset"
    ]
  ),
  
  lower_depth = c(
    climate_events$depth_cm[
      climate_events$event ==
        "Holocene onset"
    ],
    
    climate_events$depth_cm[
      climate_events$event ==
        "YD onset"
    ],
    
    climate_events$depth_cm[
      climate_events$event ==
        "B-A onset"
    ],
    
    max(
      strat_depth,
      na.rm = TRUE
    )
  )
) |>
  
  mutate(
    midpoint_depth = (
      upper_depth +
        lower_depth
    ) / 2
  )

# ==============================================================================
# 6N. REGIONAL PROBOSCIDEAN FOSSIL CONTEXT
# ==============================================================================

fossil_events <- tibble(
  event = c(
    "Mammoth",
    "Mastodon"
  ),
  
  age_old_calBP = c(
    13390,
    12810
  ),
  
  age_young_calBP = c(
    13020,
    12560
  )
) |>
  
  mutate(
    age_mid_calBP = (
      age_old_calBP +
        age_young_calBP
    ) / 2,
    
    depth_old_cm =
      depth_from_age(
        age_old_calBP
      ),
    
    depth_young_cm =
      depth_from_age(
        age_young_calBP
      ),
    
    depth_mid_cm =
      depth_from_age(
        age_mid_calBP
      ),
    
    age_label = paste0(
      sprintf(
        "%.2f",
        age_old_calBP /
          1000
      ),
      
      "–",
      
      sprintf(
        "%.2f",
        age_young_calBP /
          1000
      ),
      
      " ka"
    ),
    
    point_color =
      unname(
        FOSSIL_COLORS[
          event
        ]
      )
  )

# ==============================================================================
# 6O. YOUNGEST PERMISSIVE SEDADNA DETECTIONS
# ==============================================================================

PERMISSIVE_VERTEBRATE_PATH <- file.path(
  TABLE_DIR,
  "08_vertebrate_family_permissive_detections.tsv"
)

if (
  !file.exists(
    PERMISSIVE_VERTEBRATE_PATH
  )
) {
  
  stop(
    "Could not find permissive vertebrate table:\n",
    PERMISSIVE_VERTEBRATE_PATH
  )
}

seda_last_detections <-
  readr::read_tsv(
    PERMISSIVE_VERTEBRATE_PATH,
    show_col_types = FALSE
  ) |>
  
  filter(
    family %in% c(
      "Elephantidae",
      "Equidae"
    )
  ) |>
  
  group_by(
    family
  ) |>
  
  slice_min(
    order_by =
      master_depth_cm,
    
    n = 1,
    
    with_ties = FALSE
  ) |>
  
  ungroup() |>
  
  mutate(
    age_calBP =
      age_from_depth(
        master_depth_cm
      ),
    
    age_label = sprintf(
      "%.2f ka",
      age_calBP /
        1000
    )
  )

# ==============================================================================
# 6P. SYNTHESIS-FIGURE LAYOUT
# ==============================================================================

SYNTHESIS_Y_BOTTOM <- 0.285
SYNTHESIS_Y_TOP <- 0.655

SYNTHESIS_X_LEFT <- 0.118
SYNTHESIS_X_RIGHT <- 0.645

# CONISS

ZONE_FIG_LEFT <- 0.650
ZONE_FIG_RIGHT <- 0.690

CLIMATE_FIG_LEFT <- 0.703
CLIMATE_FIG_RIGHT <- 0.743

VERTEBRATE_FIG_LEFT <- 0.765
VERTEBRATE_FIG_RIGHT <- 0.910

AGE_FIG_LEFT <- 0.940
AGE_FIG_RIGHT <- 0.958

AGE_TITLE_LEFT <- 0.980
AGE_TITLE_RIGHT <- 0.999

# ==============================================================================
# 7. FIGURE ANNOTATION FUNCTIONS
# ==============================================================================

# ==============================================================================
# 7A. MAIN TITLE
# ==============================================================================

add_synthesis_title <- function() {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      SYNTHESIS_X_LEFT,
      0.96,
      0.935,
      0.995
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim = c(
      0,
      1
    )
  )
  
  
  text(
    x = 0,
    
    y = 0.55,
    
    labels =
      "Mastodon Lake paleoecological synthesis",
    
    adj = c(
      0,
      0.5
    ),
    
    cex = 2.25,
    font = 2
  )
}

# ==============================================================================
# 7B. DEPTH-AXIS TITLE
# ==============================================================================

add_depth_axis_title <- function() {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      0.008,
      0.058,
      SYNTHESIS_Y_BOTTOM,
      SYNTHESIS_Y_TOP
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim = c(
      0,
      1
    )
  )
  
  
  text(
    x = 0.52,
    
    y = 0.5,
    
    labels =
      "Composite master depth (cm)",
    
    srt = 90,
    
    cex = 2.00
  )
}

# ==============================================================================
# 7C. PLANT-PANEL X AXES
# ==============================================================================

format_axis_value <- function(
    x
) {
  
  if (
    abs(
      x - round(x)
    ) < 1e-8
  ) {
    
    return(
      as.character(
        round(x)
      )
    )
  }
  
  
  format(
    round(
      x,
      1
    ),
    
    trim = TRUE,
    scientific = FALSE
  )
}

add_plant_x_axes <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  for (
    i in seq_along(
      diagram$figs
    )
  ) {
    
    one_max <-
      panel_limits[
        i,
        "maximum"
      ]
    
    
    par(
      fig =
        diagram$figs[[i]],
      
      mar = c(
        0,
        0,
        0,
        0
      ),
      
      new = TRUE,
      xpd = NA
    )
    
    
    plot.new()
    
    
    plot.window(
      xlim = c(
        0,
        one_max
      ),
      
      ylim =
        diagram$ylim,
      
      xaxs = "i",
      yaxs = "i"
    )
    
    
    axis(
      side = 1,
      
      at = c(
        0,
        one_max
      ),
      
      labels = FALSE,
      
      tcl = -0.30,
      
      lwd = 0,
      lwd.ticks = 1
    )
    
    
    axis(
      side = 1,
      
      at =
        one_max,
      
      labels =
        format_axis_value(
          one_max
        ),
      
      tick = FALSE,
      
      cex.axis = 1.68,
      
      mgp = c(
        3,
        0.65,
        0
      ),
      
      hadj = 1,
      
      lwd = 0
    )
  }
}

# ==============================================================================
# 7D. PLANT-PANEL REFERENCE LINES
# ==============================================================================

add_reference_lines_to_plant_panels <- function(
    diagram
) {
  
  for (
    one_depth in
    zone_boundaries_synthesis$boundary_depth_cm
  ) {
    
    rioja::addZone(
      diagram,
      
      upper =
        one_depth,
      
      col = "red3",
      
      lwd = 1.15
    )
  }
  
  
  for (
    one_depth in
    climate_events$depth_cm
  ) {
    
    rioja::addZone(
      diagram,
      
      upper =
        one_depth,
      
      col = "grey15",
      
      lty = 2,
      
      lwd = 1.35
    )
  }
}

# ==============================================================================
# 7E. REFERENCE LINES THROUGH DENDROGRAM
# ==============================================================================

add_reference_lines_to_dendrogram <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  last_plant_right <-
    diagram$figs[[length(diagram$figs)]][2]
  
  
  dendrogram_left <-
    last_plant_right + 0.002
  
  
  if (
    dendrogram_left >=
    SYNTHESIS_X_RIGHT
  ) {
    
    return(
      invisible(
        NULL
      )
    )
  }
  
  
  par(
    fig = c(
      dendrogram_left,
      SYNTHESIS_X_RIGHT,
      diagram$box[[3]],
      diagram$box[[4]]
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = FALSE
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim =
      diagram$ylim,
    
    yaxs = "i"
  )
  
  
  abline(
    h =
      zone_boundaries_synthesis$boundary_depth_cm,
    
    col = "red3",
    
    lwd = 1.15
  )
  
  
  abline(
    h =
      climate_events$depth_cm,
    
    col = "grey15",
    
    lty = 2,
    
    lwd = 1.35
  )
}

# ==============================================================================
# 7F. SHARED PLANT X-AXIS TITLE
# ==============================================================================

add_relative_abundance_title <- function() {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      SYNTHESIS_X_LEFT,
      SYNTHESIS_X_RIGHT,
      0.145,
      0.220
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim = c(
      0,
      1
    )
  )
  
  
  text(
    x = 0.5,
    
    y = 0.5,
    
    labels =
      "Relative authenticated plant reads (%)",
    
    cex = 1.95
  )
}

# ==============================================================================
# 7G. CONISS ZONE STRIP
# ==============================================================================

add_zone_strip <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      ZONE_FIG_LEFT,
      ZONE_FIG_RIGHT,
      diagram$box[[3]],
      diagram$box[[4]]
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim =
      diagram$ylim,
    
    yaxs = "i"
  )
  
  
  abline(
    h =
      zone_boundaries_synthesis$boundary_depth_cm,
    
    col = "red3",
    
    lwd = 1.15
  )
  
  
  text(
    x = 0.45,
    
    y =
      zone_ranges_plot$midpoint_depth_cm,
    
    labels =
      zone_ranges_plot$zone,
    
    cex = 1.75,
    
    font = 2
  )
  
  
  mtext(
    "CONISS",
    
    side = 3,
    
    line = 0.35,
    
    cex = 1.62,
    
    font = 2
  )
}

# ==============================================================================
# 7H. CLIMATE STRIP
# ==============================================================================

add_climate_strip <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      CLIMATE_FIG_LEFT,
      CLIMATE_FIG_RIGHT,
      diagram$box[[3]],
      diagram$box[[4]]
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim =
      diagram$ylim,
    
    yaxs = "i"
  )
  
  
  abline(
    h =
      climate_events$depth_cm,
    
    col = "grey15",
    
    lty = 2,
    
    lwd = 1.35
  )
  
  
  text(
    x = 0.53,
    
    y =
      climate_intervals$midpoint_depth,
    
    labels =
      climate_intervals$climate,
    
    cex = 1.75
  )
  
  
  mtext(
    "Climate",
    
    side = 3,
    
    line = 0.35,
    
    cex = 1.62,
    
    font = 2
  )
}

# ==============================================================================
# 7I. VERTEBRATE CHRONOLOGY PANEL
# ==============================================================================

add_vertebrate_panel <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      VERTEBRATE_FIG_LEFT,
      VERTEBRATE_FIG_RIGHT,
      diagram$box[[3]],
      diagram$box[[4]]
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim =
      diagram$ylim,
    
    yaxs = "i"
  )
  
  
  
  segments(
    x0 = 0.53,
    
    y0 = min(
      strat_depth
    ),
    
    x1 = 0.53,
    
    y1 = max(
      strat_depth
    ),
    
    col = "grey90",
    
    lwd = 1
  )
  
  
  abline(
    h =
      climate_events$depth_cm,
    
    col = "grey15",
    
    lty = 2,
    
    lwd = 1.35
  )
  
  
  header_y <- min(
    strat_depth,
    na.rm = TRUE
  ) + 27
  
  
  
  text(
    x = 0.20,
    
    y = header_y,
    
    labels =
      "Fossil LO",
    
    font = 2,
    
    cex = 1.72
  )
  
  
  text(
    x = 0.79,
    
    y = header_y,
    
    labels =
      "sedaDNA",
    
    font = 2,
    
    cex = 1.72
  )
  
  
  mtext(
    "Vertebrate chronology",
    
    side = 3,
    
    line = 1.8,
    
    cex = 1.82,
    
    font = 2
  )
  
  
  # ---------------------------------------------------------------------------
  # REGIONAL FOSSIL LAST OCCURRENCES
  # ---------------------------------------------------------------------------
  
  fossil_x <- 0.06
  fossil_label_x <- 0.13
  
  
  fossil_events_plot <- fossil_events |>
    
    mutate(
      label_depth = case_when(
        event == "Mastodon" ~
          depth_mid_cm - 34,
        
        event == "Mammoth" ~
          depth_mid_cm + 34,
        
        TRUE ~
          depth_mid_cm
      )
    )
  
  
  for (
    i in seq_len(
      nrow(
        fossil_events_plot
      )
    )
  ) {
    
    this_color <-
      fossil_events_plot$point_color[[i]]
    
    
    segments(
      x0 =
        fossil_x,
      
      y0 =
        fossil_events_plot$depth_young_cm[[i]],
      
      x1 =
        fossil_x,
      
      y1 =
        fossil_events_plot$depth_old_cm[[i]],
      
      col =
        this_color,
      
      lwd = 2.0
    )
    
    
    points(
      x =
        fossil_x,
      
      y =
        fossil_events_plot$depth_mid_cm[[i]],
      
      pch = 21,
      
      bg =
        this_color,
      
      col =
        this_color,
      
      lwd = 1.2,
      
      cex = 1.60
    )
    
    
    segments(
      x0 =
        fossil_x,
      
      y0 =
        fossil_events_plot$depth_mid_cm[[i]],
      
      x1 =
        fossil_label_x - 0.012,
      
      y1 =
        fossil_events_plot$label_depth[[i]],
      
      col =
        this_color,
      
      lwd = 1.2
    )
    
    
    text(
      x =
        fossil_label_x,
      
      y =
        fossil_events_plot$label_depth[[i]],
      
      labels = paste0(
        fossil_events_plot$event[[i]],
        "\n",
        fossil_events_plot$age_label[[i]]
      ),
      
      adj = c(
        0,
        0.5
      ),
      
      cex = 1.65
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # SEDADNA DETECTIONS
  # ---------------------------------------------------------------------------
  
  if (
    nrow(
      seda_last_detections
    ) > 0
  ) {
    
    seda_x <- 0.66
    
    
    seda_last_detections_plot <-
      seda_last_detections |>
      
      mutate(
        label_depth = case_when(
          family == "Equidae" ~
            360,
          
          family == "Elephantidae" ~
            525,
          
          TRUE ~
            master_depth_cm + 70
        ),
        
        label_x = case_when(
          family == "Equidae" ~
            0.77,
          
          family == "Elephantidae" ~
            0.77,
          
          TRUE ~
            0.77
        ),
        
        point_color =
          unname(
            SEDA_COLORS[
              family
            ]
          )
      )
    
    
    
    points(
      x = rep(
        seda_x,
        nrow(
          seda_last_detections_plot
        )
      ),
      
      y =
        seda_last_detections_plot$master_depth_cm,
      
      pch = 19,
      
      col =
        seda_last_detections_plot$point_color,
      
      cex = 1.70
    )
    
    
    for (
      i in seq_len(
        nrow(
          seda_last_detections_plot
        )
      )
    ) {
      
      this_family <-
        seda_last_detections_plot$family[[i]]
      
      
      this_color <-
        seda_last_detections_plot$point_color[[i]]
      
      
      this_depth <-
        seda_last_detections_plot$master_depth_cm[[i]]
      
      
      this_label_x <-
        seda_last_detections_plot$label_x[[i]]
      
      
      this_label_depth <-
        seda_last_detections_plot$label_depth[[i]]
      
      
      # -----------------------------------------------------------------------
      # EQUIDAE ONLY:
      # -----------------------------------------------------------------------
      
      if (
        this_family == "Equidae"
      ) {
        
        segments(
          x0 =
            seda_x,
          
          y0 =
            this_depth,
          
          x1 =
            this_label_x - 0.018,
          
          y1 =
            this_label_depth,
          
          col =
            this_color,
          
          lwd = 1.55
        )
      }
      
      
      
      
      text(
        x =
          this_label_x,
        
        y =
          this_label_depth,
        
        labels = paste0(
          this_family,
          "\n",
          seda_last_detections_plot$age_label[[i]]
        ),
        
        adj = c(
          0,
          0.5
        ),
        
        cex = 1.65
      )
    }
  }
}

# ==============================================================================
# 7J. AGE-AXIS NUMBERS
# ==============================================================================

add_age_axis_to_synthesis <- function(
    diagram
) {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      AGE_FIG_LEFT,
      AGE_FIG_RIGHT,
      diagram$box[[3]],
      diagram$box[[4]]
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim =
      diagram$ylim,
    
    yaxs = "i"
  )
  
  
  axis(
    side = 4,
    
    at =
      strat_age_axis$depth_cm,
    
    labels = sprintf(
      "%.1f",
      strat_age_axis$age_ka
    ),
    
    las = 1,
    
    cex.axis = 1.78,
    
    tck = -0.09,
    
    mgp = c(
      1.8,
      0.02,
      0
    ),
    
    lwd = 0,
    lwd.ticks = 1
  )
  
  
  segments(
    x0 = 1,
    
    y0 = min(
      diagram$ylim
    ),
    
    x1 = 1,
    
    y1 = max(
      diagram$ylim
    )
  )
}

# ==============================================================================
# 7K. AGE-AXIS TITLE
# ==============================================================================

add_age_axis_title <- function() {
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      AGE_TITLE_LEFT,
      AGE_TITLE_RIGHT,
      SYNTHESIS_Y_BOTTOM,
      SYNTHESIS_Y_TOP
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim = c(
      0,
      1
    )
  )
  
  
  text(
    x = 0.50,
    
    y = 0.5,
    
    labels =
      "Age (ka cal BP)",
    
    srt = 90,
    
    cex = 1.95
  )
}

# ==============================================================================
# 7L. BOTTOM LEGEND
# ==============================================================================

add_bottom_legend <- function() {
  
  used_groups <- unique(
    selected_groups
  )
  
  
  used_groups <- used_groups[
    !is.na(
      used_groups
    )
  ]
  
  
  used_group_colors <- unname(
    FUNCTIONAL_GROUP_COLORS[
      used_groups
    ]
  )
  
  
  used_group_labels <- used_groups |>
    
    stringr::str_replace_all(
      "_",
      " "
    ) |>
    
    stringr::str_to_sentence()
  
  
  old_par <- par(
    c(
      "fig",
      "mar",
      "usr",
      "xpd"
    )
  )
  
  
  on.exit(
    par(
      fig = old_par$fig,
      mar = old_par$mar,
      usr = old_par$usr,
      xpd = old_par$xpd
    ),
    add = TRUE
  )
  
  
  par(
    fig = c(
      0.025,
      0.950,
      0.002,
      0.145
    ),
    
    mar = c(
      0,
      0,
      0,
      0
    ),
    
    new = TRUE,
    xpd = NA
  )
  
  
  plot.new()
  
  
  plot.window(
    xlim = c(
      0,
      1
    ),
    
    ylim = c(
      0,
      1
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # FUNCTIONAL-GROUP KEY
  # ---------------------------------------------------------------------------
  
  legend_y <- 0.74
  
  
  group_x <- c(
    0.10,   # Trees
    0.24,   # Shrubs
    0.37,   # Ferns
    0.52,   # Submerged macrophytes
    0.79    # Emergent & wetland
  )
  
  
  if (
    length(
      used_group_labels
    ) != length(
      group_x
    )
  ) {
    
    group_x <- seq(
      0.10,
      0.82,
      length.out =
        length(
          used_group_labels
        )
    )
  }
  
  
  for (
    i in seq_along(
      used_group_labels
    )
  ) {
    
    rect(
      xleft =
        group_x[[i]] - 0.018,
      
      ybottom =
        legend_y - 0.025,
      
      xright =
        group_x[[i]] - 0.007,
      
      ytop =
        legend_y + 0.025,
      
      col =
        used_group_colors[[i]],
      
      border = NA
    )
    
    
    text(
      x =
        group_x[[i]],
      
      y =
        legend_y,
      
      labels =
        used_group_labels[[i]],
      
      adj = c(
        0,
        0.5
      ),
      
      cex = 1.95
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # CONISS / CLIMATE KEY
  # ---------------------------------------------------------------------------
  
  legend(
    x = 0.5,
    
    y = 0.18,
    
    legend = c(
      "CONISS boundary",
      "Climatic transition"
    ),
    
    col = c(
      "red3",
      "grey15"
    ),
    
    lty = c(
      1,
      2
    ),
    
    lwd = c(
      1.8,
      1.35
    ),
    
    horiz = TRUE,
    
    xjust = 0.5,
    yjust = 0.5,
    
    bty = "n",
    
    cex = 1.85,
    
    x.intersp = 0.80
  )
}

# ==============================================================================
# 8. DRAW FINAL PALEOECOLOGICAL SYNTHESIS FIGURE
# ==============================================================================

grDevices::png(
  filename = file.path(
    PLOT_DIR,
    "Figure_11_paleoecological_synthesis.png"
  ),
  
  width = 10000,
  
  height = 4900,
  
  res = 300,
  
  pointsize = 18
)

tryCatch(
  {
    
    par(
      oma = c(
        1.5,
        2.0,
        0.8,
        2.5
      ),
      
      cex.axis = 1.75,
      
      xpd = NA
    )
    
    
    diagram <- rioja::strat.plot(
      d =
        strat_matrix,
      
      yvar =
        strat_depth,
      
      y.rev =
        REVERSE_DEPTH_AXIS,
      
      
      # ------------------------------------------------------------------------
      # EQUAL-WIDTH PANELS WITH PANEL-SPECIFIC PERCENT SCALES
      # ------------------------------------------------------------------------
      
      scale.percent = FALSE,
      
      graph.widths = rep(
        1,
        ncol(
          strat_matrix
        )
      ),
      
      minmax =
        panel_limits,
      
      scale.minmax = TRUE,
      
      x.axis = FALSE,
      
      
      # ------------------------------------------------------------------------
      # POSITION
      # ------------------------------------------------------------------------
      
      xLeft =
        SYNTHESIS_X_LEFT,
      
      xRight =
        SYNTHESIS_X_RIGHT,
      
      yBottom =
        SYNTHESIS_Y_BOTTOM,
      
      yTop =
        SYNTHESIS_Y_TOP,
      
      
      # ------------------------------------------------------------------------
      # LABELS
      # ------------------------------------------------------------------------
      
      title = "",
      
      ylabel = "",
      
      ann = FALSE,
      
      x.names =
        synthesis_display_names,
      
      
      # ------------------------------------------------------------------------
      # CURVES
      # ------------------------------------------------------------------------
      
      plot.line = FALSE,
      
      plot.bar = FALSE,
      
      plot.poly = TRUE,
      
      col.poly =
        selected_colors,
      
      col.poly.line =
        NA,
      
      
      # ------------------------------------------------------------------------
      # PANEL SPACING AND TEXT
      # ------------------------------------------------------------------------
      
      xSpace = 0.006,
      
      cex.xlabel = 1.92,
      
      srt.xlabel = 90,
      
      cex.title = 2.05,
      
      cex.ylabel = 1.95,
      
      cex.yaxis = 1.85,
      
      cex.axis = 1.68,
      
      mgp = c(
        2.2,
        0.72,
        0
      ),
      
      
      # ------------------------------------------------------------------------
      # CONISS DENDROGRAM
      # ------------------------------------------------------------------------
      
      clust =
        coniss,
      
      clust.width =
        0.060
    )
    
    
    # --------------------------------------------------------------------------
    # REFERENCE LINES
    # --------------------------------------------------------------------------
    
    add_reference_lines_to_plant_panels(
      diagram
    )
    
    
    add_reference_lines_to_dendrogram(
      diagram
    )
    
    
    # --------------------------------------------------------------------------
    # AXES AND TITLES
    # --------------------------------------------------------------------------
    
    add_plant_x_axes(
      diagram
    )
    
    
    add_synthesis_title()
    
    
    add_depth_axis_title()
    
    
    add_relative_abundance_title()
    
    
    # --------------------------------------------------------------------------
    # RIGHT-HAND CONTEXTUAL PANELS
    # --------------------------------------------------------------------------
    
    add_zone_strip(
      diagram
    )
    
    
    add_climate_strip(
      diagram
    )
    
    
    add_vertebrate_panel(
      diagram
    )
    
    
    add_age_axis_to_synthesis(
      diagram
    )
    
    
    add_age_axis_title()
    
    
    add_bottom_legend()
  },
  
  finally = {
    
    grDevices::dev.off()
  }
)

# ==============================================================================
# 9. SAVE SUPPORTING TABLES
# ==============================================================================

write_tsv(
  zone_membership_synthesis |>
    
    arrange(
      desc(
        master_depth_cm
      )
    ),
  
  file.path(
    TABLE_DIR,
    "Figure_11_CONISS_zone_membership.tsv"
  )
)

write_tsv(
  climate_events,
  
  file.path(
    TABLE_DIR,
    "Figure_11_climate_transitions.tsv"
  )
)

write_tsv(
  fossil_events,
  
  file.path(
    TABLE_DIR,
    "Figure_11_proboscidean_fossil_occurrences.tsv"
  )
)

write_tsv(
  seda_last_detections,
  
  file.path(
    TABLE_DIR,
    "Figure_11_vertebrate_shallowest_detections.tsv"
  )
)
