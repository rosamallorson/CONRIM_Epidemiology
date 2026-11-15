library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(broom)
library(ggplot2)
library(janitor)
library(busoga)

source("scripts/output_helpers.R")

# Loads and cleans all data for Forest_plots.R and Barplot_NMR.R (both source this file).
# Run from the CONRIM project folder (paths below are relative to it).
# Produces: analysis_df, analysis_wash, clean_wash_vars

# Load, clean and merge the survey data with the busoga package.
# Fresh exports are downloaded from ODK Central (needs ODK_EMAIL / ODK_PASSWORD in the
# environment); set options(conrim.download = FALSE) to reuse the CSVs in JinjaFiles/
# under BUSOGA_ROOT. The result is cached for the session in `busoga_prep`.
if (!exists("busoga_prep")) {
  busoga_prep <- prepare_survey_data(
    busoga_paths(),
    download = getOption("conrim.download", TRUE)
  )
  odk_logout()
}

# One row per report (short form, with district/subcounty) and best 28-day follow-up per participant
short_raw <- busoga_prep$survey %>% janitor::clean_names() %>%
  rename(survey_participant_district = district, survey_participant_subcounty = subcounty)
fu_raw    <- busoga_prep$followup %>% janitor::clean_names() %>%
  rename(infant_alive = fu_infant_alive, infant_dob = fu_infant_dob, infant_dod = fu_infant_dod)

# Standardize helper for "alive" coding (handles "yes"/"no")
is_alive <- function(x) {
  x_clean <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_clean %in% c("yes") ~ TRUE,
    x_clean %in% c("no") ~ FALSE,
    TRUE ~ NA
  )
 }

# Prepare short form (First Visit)
sf <- short_raw %>%
  transmute(
    survey_participant_number = trimws(as.character(survey_participant_number)),
    survey_survey_date        = as.Date(survey_survey_date),
    infant_dob_sf             = as.Date(infant_dob),
    infant_alive_sf           = is_alive(infant_alive),
    infant_dod_sf             = as.Date(infant_dod),
    infant_cod_sf             = infant_causeofdeath,
    infant_cod_other_sf       = infant_causeofdeath_other,
    infant_sex                = factor(infant_sex),
    infant_cordcare, infant_cordcare_other, environment_eat_drink,
    environment_practices, environment_practices_other, environment_childpractices, environment_statepractices, environment_statepractices_other,
    environment_watersrc_rainy, environment_watersrc_rainy_other,
    environment_watersrc_dry, environment_watersrc_dry_other,
    environment_handwash_when, environment_handwash, environment_latrine, environment_animals, environment_animals_other,
    survey_participant_district, survey_participant_subcounty, pregnancy_folic_uptake, introduction_hh_twins
  )

# Prepare 28-day follow-up (FU)
fu <- fu_raw %>%
  transmute(
    survey_participant_number = trimws(as.character(survey_participant_number)),
    address_age_at_followup    = as.numeric(address_age_at_followup),
    infant_dob_fu             = as.Date(infant_dob),
    infant_alive_fu           = is_alive(infant_alive),
    infant_dod_fu             = as.Date(infant_dod),
    infant_cod_fu             = infant_causeofdeath,
    infant_cod_other_fu       = infant_causeofdeath_other
  )

# Deduplicate Short Form (1 row per participant)
sf_clean <- sf %>%
  filter(!is.na(survey_participant_number), survey_participant_number != "") %>%
  arrange(survey_participant_number, is.na(infant_alive_sf), survey_survey_date) %>%
  distinct(survey_participant_number, .keep_all = TRUE)

# Deduplicate Follow-Up (1 row per participant, prioritizing known vital status & latest visit)
fu_clean <- fu %>%
  filter(!is.na(survey_participant_number), survey_participant_number != "") %>%
  arrange(survey_participant_number, is.na(infant_alive_fu), desc(address_age_at_followup)) %>%
  distinct(survey_participant_number, .keep_all = TRUE)

# Join the deduplicated tables
merged_data <- sf_clean %>%
  left_join(fu_clean, by = "survey_participant_number", relationship = "one-to-one") %>%
  mutate(
    dob = coalesce(infant_dob_sf, infant_dob_fu),
    dod = coalesce(infant_dod_sf, infant_dod_fu),
    age_at_death_days = as.numeric(dod - dob),
    dob_mismatch = !is.na(infant_dob_sf) & !is.na(infant_dob_fu) & (infant_dob_sf != infant_dob_fu),

# Audit status across First Visit -> Follow-Up
    survival_trajectory = case_when(
      infant_alive_sf == TRUE  & infant_alive_fu == TRUE  ~ "Valid: Alive -> Alive",
      infant_alive_sf == TRUE  & infant_alive_fu == FALSE ~ "Valid: Alive -> Died at FU",
      infant_alive_sf == FALSE & is.na(infant_alive_fu)   ~ "Valid: Died at First Visit (No FU)",
      infant_alive_sf == FALSE & infant_alive_fu == FALSE ~ "Valid: Died at First Visit -> Confirmed Dead",
      infant_alive_sf == TRUE  & is.na(infant_alive_fu)   ~ "Pending/Lost to FU: Alive at First Visit, Missing FU",
      infant_alive_sf == FALSE & infant_alive_fu == TRUE  ~ "ERROR: Dead at First Visit -> Alive at FU",
      is.na(infant_alive_sf)                              ~ "ERROR: Missing First Visit Status",
      TRUE                                                ~ "Other / Review"
    ),
    date_error = case_when(
      !is.na(dod) & !is.na(dob) & dod < dob ~ "ERROR: Death date before DOB",
      !is.na(dod) & !is.na(survey_survey_date) & infant_alive_sf == TRUE & dod < survey_survey_date ~ "ERROR: Alive at survey but DOD is earlier",
      !is.na(age_at_death_days) & age_at_death_days < 0 ~ "ERROR: Negative age at death",
      TRUE ~ "OK"
    )
  ) 
# Save/view records that don't make sense for data cleaning log
inconsistent_records <- merged_data %>%
  filter(str_starts(survival_trajectory, "ERROR") | date_error != "OK" | dob_mismatch) %>%
  select(
    survey_participant_number, survey_survey_date,
    infant_dob_sf, infant_dob_fu,
    infant_alive_sf, infant_alive_fu,
    infant_dod_sf, infant_dod_fu,
    survival_trajectory, date_error
  )

save_table(inconsistent_records, "inconsistent_records")
if (interactive()) View(inconsistent_records)

# Filter out impossible trajectories/dates and define the 28-day outcome:
# sb_death_28d = 1 for a stillbirth or a death by day 28, 0 if alive past day 28
# (stillbirths are INCLUDED here; nmr_wash below is the live-births-only population).
analysis_df <- merged_data %>%
  filter(!str_starts(survival_trajectory, "ERROR"), date_error == "OK") %>%
  mutate(
    died_any = (infant_alive_sf == FALSE) | (!is.na(infant_alive_fu) & infant_alive_fu == FALSE),
    stillbirth = as.integer(infant_alive_sf == FALSE & infant_cod_sf %in% "Stillbirth"),
    sb_death_28d = case_when(
      # Died within 28 days of birth (0 to 28 days)
      died_any & !is.na(age_at_death_days) & age_at_death_days >= 0 & age_at_death_days <= 28 ~ 1L,
      # Died at first visit or 28-day FU without an exact DOD recorded
      died_any & is.na(age_at_death_days) & (is.na(address_age_at_followup) | address_age_at_followup <= 28) ~ 1L,
      # Died after day 28 (post-neonatal death -> survived neonatal period)
      died_any & !is.na(age_at_death_days) & age_at_death_days > 28 ~ 0L,
      # Confirmed alive at 28-day follow-up
      infant_alive_fu == TRUE ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(sb_death_28d))

# Print overall rates (per 1,000): stillbirth-or-death per birth, and the true NMR per live birth
analysis_df %>%
  summarise(
    births                       = n(),
    stillbirths                  = sum(stillbirth),
    neonatal_deaths              = sum(sb_death_28d) - stillbirths,
    stillbirth_or_death_per_1000 = round(sum(sb_death_28d) / births * 1000, 1),
    nmr_per_1000                 = round(neonatal_deaths / (births - stillbirths) * 1000, 1)
  ) %>%
  print()

# 1. Inspect top raw categories across all 7 WASH variables (including eat_drink)
raw_wash_vars <- c(
  "environment_eat_drink",
  "environment_watersrc_rainy",
  "environment_watersrc_dry",
  "environment_handwash",
  "environment_handwash_when",
  "environment_latrine",
  "environment_animals"
)

lapply(analysis_df[raw_wash_vars], function(x) head(sort(table(x), decreasing = TRUE), 8))

# Explicit water-source mapping (matched on lower-case labels with underscores removed,
# so rainy "Tap_water" and dry "Tapwater" both map to "tapwater").
# Harvested rainwater is counted as protected (assumption: JMP "improved" source).
protected_sources   <- c("borehole", "tapwater", "protectedspring", "protectedwell", "harvestedrainwater")
unprotected_sources <- c("shallowwell", "unprotectedwell")

classify_water <- function(x) {
  key <- gsub("_", "", tolower(trimws(x)))
  case_when(
    key %in% protected_sources   ~ "Yes",
    key %in% unprotected_sources ~ "No",
    TRUE                         ~ NA_character_
  )
}

# STEP 2: Clean Multi-Select & Text WASH Variables Using Exact ODK Labels
# (a function so the same recoding applies to the mortality and stillbirth populations)
add_wash_vars <- function(df) {
  df %>%
  mutate(
    # 1. Latrine (Yes vs. No; blanks -> NA)
    has_latrine = case_when(
      tolower(trimws(environment_latrine)) == "yes" ~ "Yes",
      tolower(trimws(environment_latrine)) == "no"  ~ "No",
      TRUE ~ NA_character_
    ),

    # 2. Eat/Drink given (Any reported item vs. Blank/None)
    any_eat_drink = case_when(
      !is.na(environment_eat_drink) & trimws(environment_eat_drink) != "" ~ "Yes",
      TRUE ~ "No"
    ),

    # 3. Animals present (None_4rm_above -> No; Poultry/Goats/Cows -> Yes; blanks -> NA)
    animals_present = case_when(
      is.na(environment_animals) | trimws(environment_animals) == "" ~ NA_character_,
      str_detect(tolower(environment_animals), "none_4rm_above") ~ "No",
      TRUE ~ "Yes"
    ),

    # 4. Handwashing agent (Soap/Sanitiser vs. Water only; blanks -> NA)
    handwash_soap = case_when(
      is.na(environment_handwash) | trimws(environment_handwash) == "" ~ NA_character_,
      str_detect(tolower(environment_handwash), "soap|sanitiser") ~ "Yes",
      TRUE ~ "No"
    ),

    # 5. Handwashing at all 3 critical times (before_eating, before_handling, before_feeding_baby)
    handwash_all_3_times = case_when(
      is.na(environment_handwash_when) | trimws(environment_handwash_when) == "" ~ NA_character_,
      str_detect(tolower(environment_handwash_when), "before_eating") &
        str_detect(tolower(environment_handwash_when), "before_handling") &
        str_detect(tolower(environment_handwash_when), "before_feeding_baby") ~ "Yes",
      TRUE ~ "No"
    ),

    # 6. Rainy season water source (Protected vs. Unprotected; blanks, "_specify" & unmapped -> NA)
    protected_water_rainy = classify_water(environment_watersrc_rainy),

    # 7. Dry season water source (Protected vs. Unprotected; blanks, "_specify" & unmapped -> NA)
    protected_water_dry = classify_water(environment_watersrc_dry)
  )
}

analysis_wash <- add_wash_vars(analysis_df)

clean_wash_vars <- c(
  "any_eat_drink",
  "has_latrine",
  "animals_present",
  "handwash_soap",
  "handwash_all_3_times",
  "protected_water_rainy",
  "protected_water_dry"
)

# Verify clean counts before plotting
lapply(analysis_wash[clean_wash_vars], table, useNA = "ifany")

# Live births only (stillbirths excluded): population for the true neonatal mortality rate.
# neonatal_death_28d = 1 if died by day 28, 0 if alive past day 28.
nmr_wash <- analysis_wash %>%
  filter(stillbirth == 0L) %>%
  mutate(neonatal_death_28d = sb_death_28d)

# Multi-category WASH variables (all individual categories kept), used by the
# detailed odds-ratio forest plots. Blanks/_specify -> NA.
add_wash_cats <- function(df) {
  df %>%
    mutate(
      # Rainy water source
      rainy_water = case_when(
        is.na(environment_watersrc_rainy) | trimws(environment_watersrc_rainy) %in% c("", "_specify") ~ NA_character_,
        TRUE ~ trimws(as.character(environment_watersrc_rainy))
      ),
      # Dry water source
      dry_water = case_when(
        is.na(environment_watersrc_dry) | trimws(environment_watersrc_dry) %in% c("", "_specify") ~ NA_character_,
        TRUE ~ trimws(as.character(environment_watersrc_dry))
      ),
      # Handwashing agent (keep soap, water, sanitiser separately)
      handwash_agent = case_when(
        is.na(environment_handwash) | trimws(environment_handwash) == "" ~ NA_character_,
        TRUE ~ trimws(as.character(environment_handwash))
      ),
      # Latrine (Yes vs. No)
      latrine_status = case_when(
        tolower(trimws(environment_latrine)) == "yes" ~ "Yes",
        tolower(trimws(environment_latrine)) == "no"  ~ "No",
        TRUE ~ NA_character_
      ),
      # Eat/Drink categories (standardize Milk & Breastmilk spellings)
      eat_drink_cat = case_when(
        is.na(environment_eat_drink) | trimws(environment_eat_drink) == "" ~ "None / Blank",
        str_detect(tolower(environment_eat_drink), "breast") ~ "Breastmilk",
        str_detect(tolower(environment_eat_drink), "milk")   ~ "Milk (Non-breast)",
        str_detect(tolower(environment_eat_drink), "bcg")    ~ "BCG",
        TRUE ~ "Other Item"
      ),
      # Multi-select Animals -> individual animal indicators (Yes vs. No)
      animal_poultry = if_else(is.na(environment_animals) | trimws(environment_animals) == "", NA_character_,
                               if_else(str_detect(environment_animals, "Poultry"), "Yes", "No")),
      animal_goats   = if_else(is.na(environment_animals) | trimws(environment_animals) == "", NA_character_,
                               if_else(str_detect(environment_animals, "Goats"), "Yes", "No")),
      animal_cows    = if_else(is.na(environment_animals) | trimws(environment_animals) == "", NA_character_,
                               if_else(str_detect(environment_animals, "Cows"), "Yes", "No")),

      # Multi-select Handwash When -> individual timing indicators (Yes vs. No)
      hw_before_eating   = if_else(is.na(environment_handwash_when) | trimws(environment_handwash_when) == "", NA_character_,
                                   if_else(str_detect(environment_handwash_when, "before_eating"), "Yes", "No")),
      hw_before_handling = if_else(is.na(environment_handwash_when) | trimws(environment_handwash_when) == "", NA_character_,
                                   if_else(str_detect(environment_handwash_when, "before_handling"), "Yes", "No")),
      hw_before_feeding  = if_else(is.na(environment_handwash_when) | trimws(environment_handwash_when) == "", NA_character_,
                                   if_else(str_detect(environment_handwash_when, "before_feeding_baby"), "Yes", "No"))
    )
}

analysis_wash_cats <- add_wash_cats(analysis_df)
nmr_wash_cats      <- add_wash_cats(nmr_wash)

multi_cat_vars <- c(
  "eat_drink_cat",
  "rainy_water",
  "dry_water",
  "handwash_agent",
  "latrine_status",
  "animal_poultry",
  "animal_goats",
  "animal_cows",
  "hw_before_eating",
  "hw_before_handling",
  "hw_before_feeding"
)

# All births with a known, consistent first-visit status, for stillbirth analyses.
# Unlike analysis_df this keeps live births with no 28-day follow-up, because
# stillbirth is known at the first visit. Stillbirth = recorded dead with cause "Stillbirth".
births_wash <- merged_data %>%
  filter(!str_starts(survival_trajectory, "ERROR"), date_error == "OK") %>%
  mutate(stillbirth = as.integer(infant_alive_sf == FALSE & infant_cod_sf %in% "Stillbirth")) %>%
  add_wash_vars()
