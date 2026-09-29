# ==============================================================================
# SHARED SETTINGS AND HELPERS
# ==============================================================================

GENUS_DATA_PATH <- "genus_meta.rds"
FAMILY_DATA_PATH <- "family_meta.rds"
SAMPLE_METADATA_PATH <- "SS0203_sample_metadata.tsv"
FUNCTIONAL_GROUPS_PATH <- "SS0203_plant_functional_groups.tsv"
AGE_MODEL_PATH <- "age_model_MST.csv"
QC_TABLE_PATH <- "SS0203_library_QC_full.tsv"

DERIVED_DIR <- "derived_data"
PLOT_DIR <- "figures"
TABLE_DIR <- "tables"

dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# CHRONOLOGY AND PLOTTING
# ==============================================================================

SHOW_AGE_AXIS <- TRUE
AGE_AXIS_LABEL <- "Age (cal BP)"
BASE_PLOT_FONT_SIZE <- 15
REVERSE_DEPTH_AXIS <- TRUE

load_age_model <- function(path = AGE_MODEL_PATH) {
  if (!file.exists(path)) {
    stop("Age-model file not found: ", path)
  }

  age_model <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("depth_cm", "age_calBP")
  missing <- setdiff(required, names(age_model))

  if (length(missing) > 0) {
    stop("Age model is missing required columns: ", paste(missing, collapse = ", "))
  }

  age_model$depth_cm <- suppressWarnings(as.numeric(age_model$depth_cm))
  age_model$age_calBP <- suppressWarnings(as.numeric(age_model$age_calBP))
  age_model <- age_model[
    is.finite(age_model$depth_cm) & is.finite(age_model$age_calBP),
    ,
    drop = FALSE
  ]
  age_model <- age_model[order(age_model$depth_cm), , drop = FALSE]

  if (nrow(age_model) < 2) {
    stop("Age model contains fewer than two valid depth-age points.")
  }

  age_model
}

age_from_depth <- function(depth_cm, age_model = NULL) {
  if (is.null(age_model)) {
    age_model <- load_age_model()
  }

  stats::approx(
    x = age_model$depth_cm,
    y = age_model$age_calBP,
    xout = depth_cm,
    rule = 2,
    ties = "ordered"
  )$y
}

age_ka_from_depth <- function(depth_cm, age_model = NULL) {
  age_from_depth(depth_cm, age_model) / 1000
}

age_axis_labels <- function(depth_values) {
  ages <- age_from_depth(depth_values)
  ifelse(
    is.finite(ages),
    format(round(ages), big.mark = ",", scientific = FALSE, trim = TRUE),
    ""
  )
}

depth_scale_with_age <- function() {
  secondary_age_axis <- if (SHOW_AGE_AXIS) {
    ggplot2::dup_axis(name = AGE_AXIS_LABEL, labels = age_axis_labels)
  } else {
    ggplot2::waiver()
  }

  if (REVERSE_DEPTH_AXIS) {
    ggplot2::scale_y_reverse(
      sec.axis = secondary_age_axis,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    )
  } else {
    ggplot2::scale_y_continuous(
      sec.axis = secondary_age_axis,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    )
  }
}

thesis_theme <- function(base_size = BASE_PLOT_FONT_SIZE) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 3, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, lineheight = 1.05),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      strip.text = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      plot.margin = ggplot2::margin(10, 16, 10, 10)
    )
}

save_figure <- function(plot, stem, width, height, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )

  invisible(NULL)
}

# ==============================================================================
# PLANT AUTHENTICATION AND COMMUNITY SETTINGS
# ==============================================================================

SUPPORT_COLUMN <- "n_alns"
ABUNDANCE_COLUMN <- "total_n_reads_from_mD"
DAMAGE_AMPLITUDE_MIN <- 0.05
ZFIT_MIN <- 2
MIN_DATASET_ABUNDANCE_READS <- 100

DAMAGE_AMPLITUDE_COLUMN_CANDIDATES <- c(
  "A",
  "damage_amplitude",
  "MAP_damage",
  "MAPfit"
)

ZFIT_COLUMN_CANDIDATES <- c(
  "Zfit",
  "zfit",
  "MAP_significance",
  "significance"
)

FUNCTIONAL_GROUP_LEVELS <- c(
  "Trees",
  "Shrubs",
  "Forbs",
  "Graminoids",
  "Ferns",
  "Bryophytes",
  "Submerged macrophytes",
  "Emergent & wetland"
)

EXPECTED_FUNCTIONAL_GROUP_COUNTS <- c(
  "Trees" = 49L,
  "Shrubs" = 32L,
  "Forbs" = 88L,
  "Graminoids" = 17L,
  "Ferns" = 23L,
  "Bryophytes" = 5L,
  "Submerged macrophytes" = 12L,
  "Emergent & wetland" = 54L
)

FUNCTIONAL_GROUP_COLORS <- c(
  "Trees" = "#2E7D32",
  "Shrubs" = "#8D6E63",
  "Forbs" = "#AB47BC",
  "Graminoids" = "#F9A825",
  "Ferns" = "#66A61E",
  "Bryophytes" = "#00695C",
  "Submerged macrophytes" = "#00ACC1",
  "Emergent & wetland" = "#1E88E5"
)

# ==============================================================================
# DIVERSITY, ORDINATION, AND ZONATION SETTINGS
# ==============================================================================

SET_RANDOM_SEED <- 123
RAREFACTION_RETAIN_SAMPLE_FRACTION <- 0.80
RAREFACTION_REPEATS <- 100

NMDS_DIMENSIONS <- 2
NMDS_TRYMAX <- 200
NMDS_PASSIVE_MIN_PREVALENCE <- 5
MIN_SAMPLES_PER_GENUS_FOR_ORDINATION <- 2
MIN_SAMPLES_PER_GENUS_FOR_ZONATION <- 2

CONISS_N_ZONES <- 6L
CONISS_MAX_GROUPS_TO_TEST <- 12L
