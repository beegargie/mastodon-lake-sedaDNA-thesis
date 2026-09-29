# ==============================================================================
# FIGURES 14–19 — INDIVIDUAL-LIBRARY ANCIENT-DNA DAMAGE PROFILES
# ==============================================================================

source("00_config.R")

library(tidyverse)

# ==============================================================================
# 1. SETTINGS
# ==============================================================================

A_THRESHOLD <- DAMAGE_AMPLITUDE_MIN
ZFIT_THRESHOLD <- ZFIT_MIN
N_DAMAGE_POSITIONS <- 15

LIBRARIES_PER_SHEET <- 16L
SHEET_NCOL <- 4L

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

genus_path <- file.path(
  DERIVED_DIR,
  "genus_meta_prepared_rows.rds"
)

community_path <- file.path(
  DERIVED_DIR,
  "community_long_primary.rds"
)

metadata_path <- file.path(
  DERIVED_DIR,
  "sample_info.rds"
)

if (!file.exists(genus_path)) {
  stop("Missing required file: ", genus_path)
}

if (!file.exists(community_path)) {
  stop("Missing required file: ", community_path)
}

if (!file.exists(metadata_path)) {
  stop("Missing required file: ", metadata_path)
}

genus_damage <- readRDS(genus_path)
community_long <- readRDS(community_path)
metadata <- readRDS(metadata_path)

# ==============================================================================
# 3. PREPARE AUTHENTICATED PLANT DATA
# ==============================================================================

if (!"genus" %in% names(community_long)) {
  stop("community_long_primary.rds must contain a 'genus' column.")
}

if (!all(c("library_id", "genus", "A", "Zfit") %in% names(genus_damage))) {
  stop(
    "genus_meta_prepared_rows.rds must contain library_id, genus, A, and Zfit."
  )
}

if (!all(c("library_id", "master_depth_cm") %in% names(metadata))) {
  stop(
    "sample_info.rds must contain library_id and master_depth_cm."
  )
}

required_damage_columns <- c(
  paste0("fwK", 0:(N_DAMAGE_POSITIONS - 1)),
  paste0("fwN", 0:(N_DAMAGE_POSITIONS - 1)),
  paste0("fwdx", 0:(N_DAMAGE_POSITIONS - 1)),
  paste0("bwK", 0:(N_DAMAGE_POSITIONS - 1)),
  paste0("bwN", 0:(N_DAMAGE_POSITIONS - 1)),
  paste0("bwdx", 0:(N_DAMAGE_POSITIONS - 1))
)

missing_damage_columns <- setdiff(
  required_damage_columns,
  names(genus_damage)
)

if (length(missing_damage_columns) > 0) {
  stop(
    "Missing required damage columns: ",
    paste(missing_damage_columns, collapse = ", ")
  )
}

retained_plant_genera <- community_long |>
  distinct(genus) |>
  filter(
    !is.na(genus),
    genus != ""
  ) |>
  pull(genus) |>
  as.character() |>
  unique()

genus_damage <- genus_damage |>
  mutate(
    genus = as.character(genus),
    A = suppressWarnings(as.numeric(A)),
    Zfit = suppressWarnings(as.numeric(Zfit)),
    auth_pass = if ("auth_pass_primary" %in% names(genus_damage)) {
      as.logical(auth_pass_primary)
    } else {
      is.finite(A) &
        is.finite(Zfit) &
        A >= A_THRESHOLD &
        Zfit >= ZFIT_THRESHOLD
    }
  ) |>
  filter(
    genus %in% retained_plant_genera,
    auth_pass
  )

if ("is_sediment" %in% names(genus_damage)) {
  
  genus_damage <- genus_damage |>
    filter(is_sediment)
  
} else if ("core" %in% names(metadata)) {
  
  sediment_library_ids <- metadata |>
    filter(core %in% c("R", "L")) |>
    distinct(library_id) |>
    pull(library_id)
  
  genus_damage <- genus_damage |>
    filter(library_id %in% sediment_library_ids)
}

# ==============================================================================
# 4. RESHAPE TERMINAL-DAMAGE DATA
# ==============================================================================

extract_damage_end <- function(
    df,
    direction = c("forward", "reverse")
) {
  
  direction <- match.arg(direction)
  
  if (direction == "forward") {
    
    k_prefix <- "fwK"
    n_prefix <- "fwN"
    fit_prefix <- "fwdx"
    damage_label <- "5′ C→T"
    
  } else {
    
    k_prefix <- "bwK"
    n_prefix <- "bwN"
    fit_prefix <- "bwdx"
    damage_label <- "3′ G→A"
  }
  
  
  map_dfr(
    0:(N_DAMAGE_POSITIONS - 1),
    function(position) {
      
      tibble(
        library_id = df$library_id,
        genus = df$genus,
        position_from_terminus = position,
        
        K = suppressWarnings(
          as.numeric(
            df[[paste0(k_prefix, position)]]
          )
        ),
        
        N = suppressWarnings(
          as.numeric(
            df[[paste0(n_prefix, position)]]
          )
        ),
        
        fitted_frequency = suppressWarnings(
          as.numeric(
            df[[paste0(fit_prefix, position)]]
          )
        ),
        
        damage_end = damage_label
      )
    }
  )
}

damage_long <- bind_rows(
  
  extract_damage_end(
    genus_damage,
    "forward"
  ),
  
  extract_damage_end(
    genus_damage,
    "reverse"
  )
  
) |>
  mutate(
    
    smiley_x = case_when(
      
      damage_end == "5′ C→T" ~
        position_from_terminus,
      
      damage_end == "3′ G→A" ~
        (2 * N_DAMAGE_POSITIONS - 1) -
        position_from_terminus,
      
      TRUE ~
        NA_real_
    )
  )

# ==============================================================================
# 5. POOL DAMAGE SIGNAL WITHIN EACH LIBRARY
# ==============================================================================

library_damage_summary <- damage_long |>
  group_by(
    library_id,
    damage_end,
    position_from_terminus,
    smiley_x
  ) |>
  summarise(
    
    pooled_K = sum(
      K,
      na.rm = TRUE
    ),
    
    pooled_N = sum(
      N,
      na.rm = TRUE
    ),
    
    observed_frequency = if_else(
      pooled_N > 0,
      pooled_K / pooled_N,
      NA_real_
    ),
    
    weighted_fitted_frequency = {
      
      valid <- is.finite(fitted_frequency) &
        is.finite(N) &
        N > 0
      
      if (any(valid)) {
        
        weighted.mean(
          fitted_frequency[valid],
          w = N[valid],
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      }
    },
    
    .groups = "drop"
  )

# ==============================================================================
# 6. ADD DEPTH INFORMATION AND ASSIGN PLOT PAGES
# ==============================================================================

metadata_for_damage <- metadata |>
  select(
    library_id,
    master_depth_cm,
    any_of("core")
  ) |>
  distinct(
    library_id,
    .keep_all = TRUE
  )

library_damage_summary <- library_damage_summary |>
  left_join(
    metadata_for_damage,
    by = "library_id"
  )

library_sheet_info <- library_damage_summary |>
  group_by(library_id) |>
  summarise(
    
    master_depth_cm = first(
      master_depth_cm[
        is.finite(master_depth_cm)
      ],
      default = NA_real_
    ),
    
    .groups = "drop"
  ) |>
  arrange(
    master_depth_cm,
    library_id
  ) |>
  mutate(
    
    library_order = row_number(),
    
    sheet = ceiling(
      library_order /
        LIBRARIES_PER_SHEET
    ),
    
    facet_label = paste0(
      
      library_id,
      
      if_else(
        is.finite(master_depth_cm),
        
        paste0(
          "\n",
          format(
            round(master_depth_cm, 1),
            trim = TRUE
          ),
          " cm"
        ),
        
        ""
      )
    )
  )

library_sheet_data <- library_damage_summary |>
  inner_join(
    library_sheet_info |>
      select(
        library_id,
        library_order,
        sheet,
        facet_label
      ),
    by = "library_id"
  ) |>
  arrange(
    library_order,
    smiley_x
  )

maximum_damage_frequency <- max(
  c(
    0.45,
    library_sheet_data$observed_frequency,
    library_sheet_data$weighted_fitted_frequency
  ),
  na.rm = TRUE
)

library_smiley_ylim <- ceiling(
  maximum_damage_frequency / 0.05
) * 0.05

n_sheets <- max(
  library_sheet_info$sheet,
  na.rm = TRUE
)

# ==============================================================================
# 7. GENERATE GROUPED DAMAGE-PROFILE FIGURES
# ==============================================================================

dir.create(
  PLOT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

walk(
  seq_len(n_sheets),
  function(one_sheet) {
    
    one_sheet_info <- library_sheet_info |>
      filter(
        sheet == one_sheet
      )
    
    one_sheet_data <- library_sheet_data |>
      filter(
        sheet == one_sheet
      ) |>
      mutate(
        facet_label = factor(
          facet_label,
          levels = one_sheet_info$facet_label
        )
      )
    
    
    p_sheet <- ggplot(
      one_sheet_data,
      aes(
        x = smiley_x,
        group = damage_end
      )
    ) +
      
      geom_vline(
        xintercept =
          N_DAMAGE_POSITIONS - 0.5,
        linetype = "dotted",
        linewidth = 0.9
      ) +
      
      geom_line(
        aes(
          y = weighted_fitted_frequency,
          linetype = damage_end
        ),
        linewidth = 1.25,
        na.rm = TRUE
      ) +
      
      geom_point(
        aes(
          y = observed_frequency,
          shape = damage_end
        ),
        size = 2.8,
        stroke = 0.9,
        na.rm = TRUE
      ) +
      
      facet_wrap(
        ~ facet_label,
        ncol = SHEET_NCOL,
        scales = "fixed"
      ) +
      
      scale_x_continuous(
        
        breaks = c(
          0,
          5,
          10,
          14.5,
          19,
          24,
          29
        ),
        
        labels = c(
          "0",
          "5",
          "10",
          "14",
          "10",
          "5",
          "0"
        ),
        
        minor_breaks = NULL,
        
        expand = expansion(
          mult = c(
            0.04,
            0.04
          )
        )
      ) +
      
      coord_cartesian(
        ylim = c(
          0,
          library_smiley_ylim
        )
      ) +
      
      labs(
        x = "Position from read terminus (bp)",
        y = "Substitution frequency",
        linetype = NULL,
        shape = NULL
      ) +
      
      thesis_theme(
        base_size = 24
      ) +
      
      theme(
        
        strip.text = element_text(
          face = "bold",
          size = 21
        ),
        
        axis.text.x = element_text(
          size = 19,
          margin = margin(
            t = 7
          )
        ),
        
        axis.text.y = element_text(
          size = 19
        ),
        
        axis.title.x = element_text(
          face = "bold",
          size = 24,
          margin = margin(
            t = 16
          )
        ),
        
        axis.title.y = element_text(
          face = "bold",
          size = 24,
          margin = margin(
            r = 16
          )
        ),
        
        axis.ticks = element_line(
          linewidth = 0.9
        ),
        
        axis.ticks.length = grid::unit(
          0.22,
          "cm"
        ),
        
        legend.position = "bottom",
        
        legend.text = element_text(
          size = 20
        ),
        
        legend.key.width = grid::unit(
          2.5,
          "cm"
        ),
        
        legend.key.height = grid::unit(
          0.8,
          "cm"
        ),
        
        panel.spacing = grid::unit(
          1.2,
          "lines"
        )
      ) +
      
      guides(
        linetype = guide_legend(
          override.aes = list(
            linewidth = 1.6
          )
        ),
        
        shape = guide_legend(
          override.aes = list(
            size = 4
          )
        )
      )
    
    
    ggsave(
      filename = file.path(
        PLOT_DIR,
        paste0(
          "Figure_",
          sprintf("%02d", 13L + one_sheet),
          "_damage_profiles.png"
        )
      ),
      
      plot = p_sheet,
      width = 20,
      height = 15,
      dpi = 300
    )
  }
)
