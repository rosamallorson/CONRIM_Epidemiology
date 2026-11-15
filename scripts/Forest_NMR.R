# Forest plot of the true 28-day neonatal mortality rate (NMR, per 1,000 live births):
# stillbirths are excluded, by WASH exposure, overall and by infant sex, with Wilson 95% CIs.
# Run from the CONRIM project folder.
source("scripts/output_helpers.R")


# Load and clean the data once per session (skipped if already in memory;
# run rm(nmr_wash) or source("scripts/Load_data.R") to force a reload)
if (!exists("nmr_wash")) source("scripts/Load_data.R")

wash_labels <- c(
  any_eat_drink         = "Anything given to eat/drink",
  has_latrine           = "Household latrine",
  animals_present       = "Animals in household",
  handwash_soap         = "Handwashing with soap/sanitiser",
  handwash_all_3_times  = "Handwashing at all 3 key times",
  protected_water_rainy = "Protected water (rainy season)",
  protected_water_dry   = "Protected water (dry season)"
)

# Wilson score interval for a proportion, returned per 1,000
wilson_per_1000 <- function(deaths, n, z = 1.96) {
  p <- deaths / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  list(low = (centre - half) * 1000, high = (centre + half) * 1000)
}

nmr_long <- nmr_wash %>%
  pivot_longer(cols = all_of(clean_wash_vars), names_to = "wash_var", values_to = "category") %>%
  filter(!is.na(category), category != "") %>%
  select(wash_var, category, infant_sex, neonatal_death_28d)

# Overall (all live births, regardless of exposure) as its own panel
overall_long <- nmr_wash %>%
  transmute(wash_var = "overall", category = "Overall", infant_sex, neonatal_death_28d)

wash_labels <- c(overall = "Overall (all live births)", wash_labels)

# One row per exposure level for all infants ("All") and for each sex
nmr_rows <- bind_rows(
  nmr_long %>% mutate(sex = "All"),
  nmr_long %>% filter(!is.na(infant_sex), infant_sex != "") %>% mutate(sex = as.character(infant_sex)),
  overall_long %>% mutate(sex = "All"),
  overall_long %>% filter(!is.na(infant_sex), infant_sex != "") %>% mutate(sex = as.character(infant_sex))
) %>%
  group_by(wash_var, category, sex) %>%
  summarise(n_births = n(), deaths = sum(neonatal_death_28d), .groups = "drop") %>%
  mutate(
    nmr     = deaths / n_births * 1000,
    ci_low  = wilson_per_1000(deaths, n_births)$low,
    ci_high = wilson_per_1000(deaths, n_births)$high,
    domain  = factor(wash_labels[wash_var], levels = wash_labels),
    sex     = factor(sex, levels = c("All", "Female", "Male"))
  )

overall_nmr <- sum(nmr_wash$neonatal_death_28d) / nrow(nmr_wash) * 1000

# Monospaced y-axis text carries label, NMR (95% CI) and deaths/N columns
nmr_text <- nmr_rows %>%
  mutate(
    row_text = paste0(category, ": ", sex),
    nmr_text = sprintf("%.1f (%.1f-%.1f)", nmr, ci_low, ci_high),
    n_text   = sprintf("%d/%d", deaths, n_births),
    ord      = (category != "No") * 10 + as.integer(sex)
  )

hdr <- c("Group", "NMR (95% CI)", "Deaths/N")
w <- c(max(nchar(c(nmr_text$row_text, hdr[1]))),
       max(nchar(c(nmr_text$nmr_text, hdr[2]))),
       max(nchar(c(nmr_text$n_text, hdr[3]))))
fmt_row <- function(a, b, c) {
  paste(str_pad(a, w[1], "right"), str_pad(b, w[2], "right"), str_pad(c, w[3], "right"), sep = "  ")
}

nmr_text <- nmr_text %>% mutate(row_label = fmt_row(row_text, nmr_text, n_text))
header_row <- tibble::tibble(
  row_label = fmt_row(hdr[1], hdr[2], hdr[3]),
  domain    = factor(wash_labels[1], levels = wash_labels)
)

nmr_plot_df <- bind_rows(nmr_text, header_row) %>%
  mutate(row_label = factor(row_label, levels = c(
    unique(row_label[!is.na(ord)][order(-ord[!is.na(ord)])]),
    header_row$row_label
  )))

p_nmr_forest <- ggplot(nmr_plot_df, aes(x = nmr, y = row_label, color = sex)) +
  geom_vline(xintercept = overall_nmr, linetype = "dashed", color = "#D9534F", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, alpha = n_births < 50), height = 0.25, na.rm = TRUE) +
  geom_point(aes(alpha = n_births < 50), size = 2.8, na.rm = TRUE) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_discrete(na.translate = FALSE) +
  # Faded points flag groups with fewer than 50 infants
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), labels = c("50 or more", "Fewer than 50"),
                     breaks = c(FALSE, TRUE), na.translate = FALSE) +
  labs(
    title    = "28-Day Neonatal Mortality Rate by WASH Exposure and Infant Sex",
    subtitle = sprintf("Deaths by day 28 per 1,000 live births (stillbirths excluded) with 95%% Wilson CIs; dashed line = overall NMR (%.1f)", overall_nmr),
    x        = "Neonatal mortality rate (per 1,000 live births)",
    y        = NULL,
    color    = "Infants",
    alpha    = "Group size"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(family = "mono", hjust = 0, size = rel(0.85)),
    strip.text.y     = element_text(angle = 0, face = "bold", hjust = 0),
    strip.background = element_rect(fill = "#F0F4F8", color = NA),
    panel.spacing.y  = unit(0.6, "lines")
  )

show_plot(p_nmr_forest, "Forest_NMR_nmr_forest", width = 13, height = 9)
