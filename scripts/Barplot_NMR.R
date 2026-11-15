# Stratified true 28-day neonatal mortality rate (NMR, per 1,000 live births) bar plot
# (WASH x infant sex). Stillbirths are excluded.
# Run from the CONRIM project folder.
source("scripts/output_helpers.R")


# Load and clean the data once per session (skipped if already in memory;
# run rm(nmr_wash) or source("scripts/Load_data.R") to force a reload)
if (!exists("nmr_wash")) source("scripts/Load_data.R")

strat_nmr <- nmr_wash %>%
  filter(!is.na(infant_sex), infant_sex != "") %>%
  pivot_longer(cols = all_of(clean_wash_vars), names_to = "wash_domain", values_to = "category") %>%
  filter(!is.na(category), category != "") %>%
  group_by(wash_domain, category, infant_sex) %>%
  summarise(
    n_births = n(),
    deaths   = sum(neonatal_death_28d),
    nmr      = (deaths / n_births) * 1000,
    # Wilson score interval (better than Wald for small n / rare events)
    p        = deaths / n_births,
    wilson_c = (p + 1.96^2 / (2 * n_births)) / (1 + 1.96^2 / n_births),
    wilson_h = 1.96 * sqrt(p * (1 - p) / n_births + 1.96^2 / (4 * n_births^2)) / (1 + 1.96^2 / n_births),
    ci_low   = (wilson_c - wilson_h) * 1000,
    ci_high  = (wilson_c + wilson_h) * 1000,
    .groups  = "drop"
  ) %>%
  select(-p, -wilson_c, -wilson_h)

overall_nmr <- sum(nmr_wash$neonatal_death_28d) / nrow(nmr_wash) * 1000

p_strat_nmr <- ggplot(strat_nmr, aes(x = category, y = nmr, fill = infant_sex)) +
  geom_hline(yintercept = overall_nmr, linetype = "dashed", color = "#D9534F", linewidth = 0.7) +
  geom_col(aes(alpha = n_births < 50, group = infant_sex), position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_text(
    aes(y = ci_high, label = paste0("n=", n_births), group = infant_sex),
    position = position_dodge(width = 0.8),
    vjust = -0.4, size = 2.8
  ) +
  facet_wrap(~ wash_domain, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  # Faded bars flag groups with fewer than 50 infants
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), labels = c("50 or more", "Fewer than 50")) +
  labs(
    title    = "28-Day Neonatal Mortality Rate Stratified by WASH Exposure and Sex",
    subtitle = sprintf("Deaths by day 28 per 1,000 live births (stillbirths excluded) with 95%% Wilson confidence intervals; n = live births in group; dashed line = overall NMR (%.1f)", overall_nmr),
    x        = NULL,
    y        = "Neonatal Mortality Rate (per 1,000 live births)",
    fill     = "Infant Sex",
    alpha    = "Group size"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

show_plot(p_strat_nmr, "Barplot_NMR_strat_nmr", width = 14, height = 10)
