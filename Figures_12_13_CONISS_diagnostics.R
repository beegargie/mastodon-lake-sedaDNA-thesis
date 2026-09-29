# ==============================================================================
# FIGURES 12–13 — CONISS BROKEN-STICK ASSESSMENT AND DENDROGRAM
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
# 4. BROKEN-STICK ASSESSMENT
# ==============================================================================

ng_to_test <- min(
  MAX_GROUPS_TO_TEST,
  nrow(mat)
)

broken_stick <- rioja::bstick(
  coniss,
  ng = ng_to_test,
  plot = FALSE
) |>
  as_tibble()

first_unsupported <- which(
  broken_stick$dispersion <
    broken_stick$bstick
)

if (
  length(
    first_unsupported
  ) > 0
) {
  
  first_unsupported_n_groups <-
    broken_stick$nGroups[
      first_unsupported[[1]]
    ]
  
  supported_zones <-
    first_unsupported_n_groups - 1L
  
  if (
    supported_zones != N_ZONES
  ) {
    
    warning(
      "The broken-stick stopping rule supports ",
      supported_zones,
      " zones, while N_ZONES is set to ",
      N_ZONES,
      "."
    )
  }
}

# ==============================================================================
# 4A. SAVE BROKEN-STICK FIGURE
# ==============================================================================

grDevices::png(
  filename = file.path(
    PLOT_DIR,
    "Figure_12_CONISS_broken_stick.png"
  ),
  
  width = 3200,
  height = 2250,
  res = 300,
  pointsize = 18
)

tryCatch(
  {
    
    par(
      mar = c(
        6.0,
        7.0,
        4.5,
        2.0
      ) + 0.1,
      
      mgp = c(
        3.8,
        0.9,
        0
      )
    )
    
    
    rioja::bstick(
      coniss,
      ng = ng_to_test,
      plot = TRUE,
      main = "",
      cex = 1.50,
      cex.axis = 1.45
    )
    
    
    legend(
      "topright",
      
      legend = c(
        "Observed CONISS",
        "Broken-stick expectation"
      ),
      
      col = c(
        "black",
        "red"
      ),
      
      lty = 1,
      pch = 1,
      bty = "n",
      cex = 1.40
    )
    
    
    mtext(
      "CONISS broken-stick assessment",
      side = 3,
      line = 2.0,
      font = 2,
      cex = 1.70
    )
  },
  
  finally = {
    grDevices::dev.off()
  }
)

# ==============================================================================
# 5. CONISS DENDROGRAM BY DEPTH WITH SECONDARY AGE AXIS
# ==============================================================================

dendrogram_age_axis <- make_age_axis(
  sample_order$master_depth_cm
)

grDevices::png(
  filename = file.path(
    PLOT_DIR,
    "Figure_13_CONISS_dendrogram.png"
  ),
  
  width = 3200,
  height = 5400,
  res = 300,
  pointsize = 18
)

tryCatch(
  {
    
    par(
      mar = c(
        6.2,
        7.2,
        4.2,
        6.2
      ) + 0.1,
      
      cex.main = 1.45
    )
    
    
    dendrogram_labels <- ifelse(
      seq_along(
        sample_order$master_depth_cm
      ) %% 2 == 1 |
        
        seq_along(
          sample_order$master_depth_cm
        ) == nrow(
          sample_order
        ),
      
      format(
        round(
          sample_order$master_depth_cm,
          1
        ),
        
        trim = TRUE,
        scientific = FALSE
      ),
      
      ""
    )
    
    
    plot(
      coniss,
      
      labels =
        dendrogram_labels,
      
      xvar =
        sample_order$master_depth_cm,
      
      hang = -1,
      horiz = TRUE,
      
      x.rev =
        REVERSE_DEPTH_AXIS,
      
      main =
        "CONISS constrained community clustering",
      
      cex = 1.12
    )
    
    
    mtext(
      "CONISS fusion height",
      side = 1,
      line = 2.8,
      cex = 1.55
    )
    
    
    mtext(
      "Composite master depth (cm)",
      side = 2,
      line = 4.8,
      cex = 1.55
    )
    
    
    axis(
      side = 4,
      
      at =
        dendrogram_age_axis$depth_cm,
      
      labels = sprintf(
        "%.1f",
        dendrogram_age_axis$age_ka
      ),
      
      las = 1,
      cex.axis = 1.48,
      lwd = 0,
      lwd.ticks = 1
    )
    
    
    dendrogram_usr <- par(
      "usr"
    )
    
    
    segments(
      x0 =
        dendrogram_usr[[2]],
      
      y0 = min(
        sample_order$master_depth_cm,
        na.rm = TRUE
      ),
      
      x1 =
        dendrogram_usr[[2]],
      
      y1 = max(
        sample_order$master_depth_cm,
        na.rm = TRUE
      ),
      
      xpd = NA
    )
    
    
    mtext(
      "Age (ka cal BP)",
      side = 4,
      line = 3.8,
      cex = 1.55
    )
  },
  
  finally = {
    grDevices::dev.off()
  }
)
