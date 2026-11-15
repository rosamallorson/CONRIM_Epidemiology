# Univariate logistic regression forest plots (WASH factors and infant sex).
# Outcome: sb_death_28d = stillbirth or death by day 28 (1) vs. alive past day 28 (0).
# "Events" in the plot text = stillbirths + deaths by day 28.
# For the true NMR outcome (stillbirths excluded) with Firth regression, see Forest_OR_firth.R.
# Run from the CONRIM project folder.
source("scripts/output_helpers.R")


# Load and clean the data once per session (skipped if already in memory;
# run rm(analysis_wash) or source("scripts/Load_data.R") to force a reload)
if (!exists("analysis_wash_cats") || !"sb_death_28d" %in% names(analysis_wash)) source("scripts/Load_data.R")
source("scripts/forest_helpers.R")

# --- PLOT 1: Univariate Forest Plot for WASH ---
wash_or <- fit_univariate_or(analysis_wash, clean_wash_vars, "sb_death_28d")

p_wash <- forest_table_plot(
  wash_or,
  title    = "Univariate Association of WASH Factors with Stillbirth or 28-Day Death",
  subtitle = "Odds Ratios (95% Wald CI) on log scale; dashed line = null effect (OR = 1.0)",
  point_size = 3.2
)

show_plot(p_wash, "Forest_plots_wash", width = 12, height = 8)

# --- PLOT 2: Univariate Forest Plot for Infant Sex ---
sex_or <- fit_univariate_or(analysis_wash, "infant_sex", "sb_death_28d", ref = c(infant_sex = "Female"))

p_sex <- forest_table_plot(
  sex_or,
  title       = "Univariate Association of Infant Sex with Stillbirth or 28-Day Death",
  subtitle    = "Odds Ratio (95% Wald CI) on log scale",
  point_color = "#E67E22",
  point_size  = 4,
  xlim        = c(0.2, 5)
)

show_plot(p_sex, "Forest_plots_sex", width = 12, height = 8)

# --- PLOT 3: Univariate Forest Plot for WASH (Multi-Category) ---
# Models use the largest category (or "No" for binary variables) as the reference group
wash_or_detailed <- fit_univariate_or(analysis_wash_cats, multi_cat_vars, "sb_death_28d")

p_wash_detailed <- forest_table_plot(
  wash_or_detailed,
  title     = "Univariate WASH Associations with Stillbirth or 28-Day Death (All Categories)",
  subtitle  = "Odds Ratios (95% Wald CI) on log scale vs. largest reference category",
  base_size = 10.5
)

show_plot(p_wash_detailed, "Forest_plots_wash_detailed", width = 12, height = 10)

# --- PLOT 4: Grouped by WASH domain ---
wash_or_grouped <- group_wash_domains(wash_or_detailed)

p_wash_grouped <- forest_table_plot(
  wash_or_grouped,
  label_col = "clean_label",
  facet_col = "domain",
  title     = "Univariate WASH Associations with Stillbirth or 28-Day Death",
  subtitle  = "Grouped by WASH domain; Odds Ratios (95% Wald CI) on log scale"
)

show_plot(p_wash_grouped, "Forest_plots_wash_grouped", width = 12, height = 11)
