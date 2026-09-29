# ==============================================================================
# PREPARE AUTHENTICATED PLANT COMMUNITY DATA
# ==============================================================================

source("00_config.R")
library(tidyverse)

# ==============================================================================
# HELPERS
# ==============================================================================

resolve_column <- function(data, candidates, label) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  found[[1]]
}

first_non_missing <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) NA else v[[1]]
}

make_exclusive_functional_group <- function(group_value, growth_habit_value, wetland_value) {
  group_clean <- str_to_lower(str_trim(coalesce(group_value, "")))
  habit_clean <- str_to_lower(str_trim(coalesce(growth_habit_value, "")))
  wetland_clean <- str_to_upper(str_trim(coalesce(wetland_value, "")))

  case_when(
    habit_clean %in% c("submerged", "floating") | group_clean == "aquatic" ~
      "Submerged macrophytes",
    habit_clean == "moss" | group_clean == "bryophyte" ~
      "Bryophytes",
    habit_clean == "tree" ~
      "Trees",
    habit_clean %in% c("shrub", "dwarf_shrub", "subshrub", "vine") ~
      "Shrubs",
    habit_clean %in% c("fern", "lycophyte") ~
      "Ferns",
    wetland_clean %in% c("OBL", "FACW") &
      habit_clean %in% c("forb", "graminoid", "horsetail") ~
      "Emergent & wetland",
    habit_clean == "graminoid" | group_clean == "graminoid" ~
      "Graminoids",
    habit_clean %in% c("forb", "horsetail") | group_clean == "forb" ~
      "Forbs",
    group_clean == "woody" ~
      "Shrubs",
    TRUE ~ NA_character_
  )
}

# ==============================================================================
# LOAD AND STANDARDIZE INPUTS
# ==============================================================================

raw <- readRDS(GENUS_DATA_PATH)
metadata <- read_tsv(SAMPLE_METADATA_PATH, show_col_types = FALSE) |>
  mutate(
    age_calBP = age_from_depth(master_depth_cm),
    age_kaBP = age_calBP / 1000
  )

groups <- read_tsv(FUNCTIONAL_GROUPS_PATH, show_col_types = FALSE)

damage_column <- resolve_column(
  raw,
  DAMAGE_AMPLITUDE_COLUMN_CANDIDATES,
  "the damage-amplitude column"
)

zfit_column <- resolve_column(
  raw,
  ZFIT_COLUMN_CANDIDATES,
  "the Zfit column"
)

required_raw <- c(
  "library_id",
  "genus",
  "rank",
  SUPPORT_COLUMN,
  ABUNDANCE_COLUMN,
  damage_column,
  zfit_column,
  "family",
  "accession.version"
)

missing_raw <- setdiff(required_raw, names(raw))
if (length(missing_raw) > 0) {
  stop("Raw genus data are missing required columns: ", paste(missing_raw, collapse = ", "))
}

# ==============================================================================
# FUNCTIONAL-GROUP LOOKUP
# ==============================================================================

functional_lookup <- groups |>
  mutate(use_in_figures = str_to_lower(as.character(use_in_figures))) |>
  filter(use_in_figures == "yes") |>
  transmute(
    genus = as.character(genus),
    classification_family = as.character(family),
    primary_group = as.character(group),
    growth_habit = as.character(growth_habit),
    wetland_nwpl_mw = as.character(wetland_nwpl_mw),
    wetland_evidence = as.character(wetland_evidence),
    biogeo = str_to_lower(as.character(biogeo)),
    classification_source = as.character(source),
    classification_note = as.character(note),
    plot_group = make_exclusive_functional_group(group, growth_habit, wetland_nwpl_mw)
  ) |>
  mutate(plot_group = factor(plot_group, levels = FUNCTIONAL_GROUP_LEVELS))

if (anyDuplicated(functional_lookup$genus) > 0) {
  stop("Duplicated genus names were found in the functional-group lookup.")
}

if (any(is.na(functional_lookup$plot_group))) {
  stop("At least one retained genus could not be assigned to a functional group.")
}

functional_group_check <- functional_lookup |>
  count(plot_group, name = "observed_genera") |>
  complete(
    plot_group = factor(FUNCTIONAL_GROUP_LEVELS, levels = FUNCTIONAL_GROUP_LEVELS),
    fill = list(observed_genera = 0L)
  ) |>
  mutate(
    expected_genera = as.integer(EXPECTED_FUNCTIONAL_GROUP_COUNTS[as.character(plot_group)]),
    matches_expected = observed_genera == expected_genera
  )

if (!all(functional_group_check$matches_expected)) {
  stop("Functional-group counts do not match the expected 280-genus classification.")
}

curated_genera <- functional_lookup$genus

# ==============================================================================
# AUTHENTICATE AND AGGREGATE PLANT READS
# ==============================================================================

metadata_fields_in_raw <- intersect(setdiff(names(metadata), "library_id"), names(raw))

dat <- raw |>
  select(-all_of(metadata_fields_in_raw)) |>
  left_join(metadata, by = "library_id") |>
  mutate(
    genus = as.character(genus),
    rank = as.character(rank),
    support_n_alns = as.numeric(.data[[SUPPORT_COLUMN]]),
    abundance_reads = as.numeric(.data[[ABUNDANCE_COLUMN]]),
    damage_amplitude = as.numeric(.data[[damage_column]]),
    zfit = as.numeric(.data[[zfit_column]]),
    A = damage_amplitude,
    Zfit = zfit,
    is_blank = is.na(core) |
      str_detect(str_to_lower(coalesce(field_sample_id, "")), "blank|control|ntc"),
    is_sediment = core %in% c("R", "L") & !is_blank,
    auth_pass_primary = is.finite(damage_amplitude) &
      is.finite(zfit) &
      damage_amplitude >= DAMAGE_AMPLITUDE_MIN &
      zfit >= ZFIT_MIN,
    in_curated_plant_set = genus %in% curated_genera
  ) |>
  left_join(functional_lookup, by = "genus")

plant_counts <- dat |>
  filter(
    rank == "genus",
    is_sediment,
    in_curated_plant_set,
    auth_pass_primary,
    !is.na(genus),
    genus != "",
    is.finite(abundance_reads),
    abundance_reads >= 0
  ) |>
  group_by(library_id, field_sample_id, core, master_depth_cm, genus) |>
  summarise(
    genus_count = sum(abundance_reads, na.rm = TRUE),
    genus_support_n_alns = sum(support_n_alns, na.rm = TRUE),
    raw_family = first_non_missing(family),
    classification_family = first_non_missing(classification_family),
    primary_group = first_non_missing(primary_group),
    growth_habit = first_non_missing(growth_habit),
    wetland_nwpl_mw = first_non_missing(wetland_nwpl_mw),
    wetland_evidence = first_non_missing(wetland_evidence),
    biogeo = first_non_missing(biogeo),
    plot_group = first_non_missing(as.character(plot_group)),
    .groups = "drop"
  ) |>
  mutate(plot_group = factor(plot_group, levels = FUNCTIONAL_GROUP_LEVELS))

genus_totals <- plant_counts |>
  group_by(genus) |>
  summarise(
    dataset_total_abundance_reads = sum(genus_count, na.rm = TRUE),
    prevalence_libraries = n_distinct(library_id[genus_count > 0]),
    .groups = "drop"
  ) |>
  right_join(functional_lookup, by = "genus") |>
  mutate(
    dataset_total_abundance_reads = replace_na(dataset_total_abundance_reads, 0),
    prevalence_libraries = replace_na(prevalence_libraries, 0L)
  )

if (nrow(genus_totals) != 280) {
  stop("The curated plant lookup does not contain exactly 280 genera.")
}

if (any(genus_totals$dataset_total_abundance_reads < MIN_DATASET_ABUNDANCE_READS)) {
  stop("At least one curated genus has fewer than 100 authenticated reads across the dataset.")
}

sample_info <- metadata |>
  filter(core %in% c("R", "L")) |>
  arrange(master_depth_cm, library_id)

community_long <- expand_grid(
  library_id = sample_info$library_id,
  genus = curated_genera
) |>
  left_join(
    sample_info |>
      select(library_id, field_sample_id, core, master_depth_cm, age_calBP, age_kaBP),
    by = "library_id"
  ) |>
  left_join(
    plant_counts |>
      select(library_id, genus, genus_count, genus_support_n_alns),
    by = c("library_id", "genus")
  ) |>
  left_join(genus_totals, by = "genus") |>
  mutate(
    genus_count = replace_na(genus_count, 0),
    genus_support_n_alns = replace_na(genus_support_n_alns, 0)
  ) |>
  group_by(library_id) |>
  mutate(
    retained_plant_genus_count_total = sum(genus_count, na.rm = TRUE),
    relative_abundance = if_else(
      retained_plant_genus_count_total > 0,
      genus_count / retained_plant_genus_count_total,
      NA_real_
    )
  ) |>
  ungroup() |>
  mutate(plot_group = factor(plot_group, levels = FUNCTIONAL_GROUP_LEVELS))

community_counts_wide <- community_long |>
  select(library_id, genus, genus_count) |>
  pivot_wider(names_from = genus, values_from = genus_count, values_fill = 0)

# ==============================================================================
# SAVE DERIVED DATA
# ==============================================================================

saveRDS(dat, file.path(DERIVED_DIR, "genus_meta_prepared_rows.rds"))
saveRDS(sample_info, file.path(DERIVED_DIR, "sample_info.rds"))
saveRDS(community_long, file.path(DERIVED_DIR, "community_long_primary.rds"))
saveRDS(community_counts_wide, file.path(DERIVED_DIR, "community_counts_wide_primary.rds"))

write_tsv(sample_info, file.path(TABLE_DIR, "sample_info_with_age.tsv"))
write_tsv(
  genus_totals |> arrange(desc(dataset_total_abundance_reads)),
  file.path(TABLE_DIR, "plant_genus_dataset_totals.tsv")
)
write_tsv(functional_group_check, file.path(TABLE_DIR, "functional_group_count_validation.tsv"))
