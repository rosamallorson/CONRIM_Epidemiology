# Stillbirth rate (per 1,000 total births): forest plot and Firth odds ratios.
# Run from the CONRIM project folder.
#
# Stillbirth = infant recorded dead at the first visit with cause "Stillbirth".
# Total births = live births + stillbirths among reports with a known, consistent
# first-visit status (births_wash, built in Load_data.R).
#
# None of the stillbirths has WASH answers (the WASH questions appear to be skipped
# for them), so the strata here are the characteristics recorded for every report.
source("scripts/output_helpers.R")


library(logistf)

# Load and clean the data once per session (skipped if already in memory;
# run rm(births_wash) or source("scripts/Load_data.R") to force a reload)
if (!exists("births_wash") || !"survey_participant_subcounty" %in% names(births_wash)) source("scripts/Load_data.R")

stratum_labels <- c(
  overall    = "Overall (all births)",
  infant_sex = "Infant sex",
  district   = "District",
  subcounty  = "Subcounty",
  folic      = "Folic acid in pregnancy",
  plurality  = "Babies in pregnancy"
)

# Wilson score interval for a proportion, returned per 1,000
wilson_per_1000 <- function(events, n, z = 1.96) {
  p <- events / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  list(low = (centre - half) * 1000, high = (centre + half) * 1000)
}

sb_data <- births_wash %>%
  transmute(
    stillbirth,
    overall    = "All births",
    infant_sex = as.character(infant_sex),
    district   = survey_participant_district,
    subcounty  = survey_participant_subcounty,
    folic      = pregnancy_folic_uptake,
    plurality  = case_when(
      introduction_hh_twins == "One"           ~ "One baby",
      introduction_hh_twins == "More_than_one" ~ "More than one baby",
      TRUE                                     ~ NA_character_
    )
  )

sb_long <- sb_data %>%
  pivot_longer(-stillbirth, names_to = "stratum", values_to = "category") %>%
  filter(!is.na(category), category != "")

sb_rows <- sb_long %>%
  group_by(stratum, category) %>%
  summarise(births = n(), stillbirths = sum(stillbirth), .groups = "drop") %>%
  mutate(
    rate    = stillbirths / births * 1000,
    ci_low  = wilson_per_1000(stillbirths, births)$low,
    ci_high = wilson_per_1000(stillbirths, births)$high,
    domain  = factor(stratum_labels[stratum], levels = stratum_labels)
  )

overall_sb <- sum(births_wash$stillbirth) / nrow(births_wash) * 1000

# Monospaced y-axis text carries label, rate (95% CI) and stillbirths/births columns
sb_text <- sb_rows %>%
  mutate(
    row_text  = category,
    rate_text = sprintf("%.1f (%.1f-%.1f)", rate, ci_low, ci_high),
    n_text    = sprintf("%d/%d", stillbirths, births)
  )

hdr <- c("Group", "Rate (95% CI)", "Stillbirths/births")
w <- c(max(nchar(c(sb_text$row_text, hdr[1]))),
       max(nchar(c(sb_text$rate_text, hdr[2]))),
       max(nchar(c(sb_text$n_text, hdr[3]))))
fmt_row <- function(a, b, c) {
  paste(str_pad(a, w[1], "right"), str_pad(b, w[2], "right"), str_pad(c, w[3], "right"), sep = "  ")
}

sb_text <- sb_text %>% mutate(row_label = fmt_row(row_text, rate_text, n_text))
header_row <- tibble::tibble(
  row_label = fmt_row(hdr[1], hdr[2], hdr[3]),
  domain    = factor(stratum_labels[1], levels = stratum_labels)
)

# Highest rates at the top of each panel
sb_plot_df <- bind_rows(sb_text, header_row) %>%
  mutate(row_label = factor(row_label, levels = c(
    unique(row_label[!is.na(rate)][order(rate[!is.na(rate)])]),
    header_row$row_label
  )))

p_sb_forest <- ggplot(sb_plot_df, aes(x = rate, y = row_label)) +
  geom_vline(xintercept = overall_sb, linetype = "dashed", color = "#D9534F", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, color = "#2C3E50", na.rm = TRUE) +
  geom_point(size = 2.8, color = "#1F78B4", na.rm = TRUE) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Stillbirth Rate by Infant and Household Characteristics",
    subtitle = sprintf("Stillbirths per 1,000 total births with 95%% Wilson CIs; dashed line = overall rate (%.1f)", overall_sb),
    x        = "Stillbirth rate (per 1,000 total births)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(family = "mono", hjust = 0, size = rel(0.85)),
    strip.text.y     = element_text(angle = 0, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "#F0F4F8", color = NA),
    panel.spacing.y  = unit(0.6, "lines")
  )

show_plot(p_sb_forest, "Forest_stillbirth_sb_forest", width = 13, height = 9)

# Univariate Firth logistic regression (penalized; suited to 44 events), reference = largest group
sb_or <- map_dfr(setdiff(names(stratum_labels), "overall"), function(v) {
  d <- sb_long %>% filter(stratum == v)
  ref_cat <- names(sort(table(d$category), decreasing = TRUE))[1]
  d$category <- relevel(factor(d$category), ref = ref_cat)
  fit <- logistf(stillbirth ~ category, data = d)
  tibble(
    stratum   = stratum_labels[[v]],
    level     = sub("^category", "", names(fit$coefficients))[-1],
    reference = ref_cat,
    or        = exp(fit$coefficients[-1]),
    ci_low    = exp(fit$ci.lower[-1]),
    ci_high   = exp(fit$ci.upper[-1]),
    p_value   = fit$prob[-1]
  )
})

sb_or %>%
  mutate(across(c(or, ci_low, ci_high), ~ round(.x, 2)), p_value = signif(p_value, 2)) %>%
  as.data.frame() %>%
  save_table("stillbirth_firth_or") %>%
  print()
