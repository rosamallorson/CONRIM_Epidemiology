# Bar plot of the combined outcome "stillbirth or death by day 28" (WASH x infant sex).
# This is NOT the neonatal mortality rate: stillbirths are counted as events.
# For the true NMR (live births only) see Barplot_NMR.R.
# Run from the CONRIM project folder.
source("scripts/output_helpers.R")


# Load and clean the data once per session (skipped if already in memory;
# run rm(analysis_wash) or source("scripts/Load_data.R") to force a reload)
if (!exists("analysis_wash") || !"sb_death_28d" %in% names(analysis_wash)) source("scripts/Load_data.R")

strat_sbd <- analysis_wash %>%
  filter(!is.na(infant_sex), infant_sex != "") %>%
  pivot_longer(cols = all_of(clean_wash_vars), names_to = "wash_domain", values_to = "category") %>%
  filter(!is.na(category), category != "") %>%
  group_by(wash_domain, category, infant_sex) %>%
  summarise(
    n_births = n(),
    events   = sum(sb_death_28d),
    rate     = (events / n_births) * 1000,
    # Wilson score interval (better than Wald for small n / rare events)
    p        = events / n_births,
    wilson_c = (p + 1.96^2 / (2 * n_births)) / (1 + 1.96^2 / n_births),
    wilson_h = 1.96 * sqrt(p * (1 - p) / n_births + 1.96^2 / (4 * n_births^2)) / (1 + 1.96^2 / n_births),
    ci_low   = (wilson_c - wilson_h) * 1000,
    ci_high  = (wilson_c + wilson_h) * 1000,
    .groups  = "drop"
  ) %>%
  select(-p, -wilson_c, -wilson_h)

overall_sbd <- sum(analysis_wash$sb_death_28d) / nrow(analysis_wash) * 1000

p_strat_sbd <- ggplot(strat_sbd, aes(x = category, y = rate, fill = infant_sex)) +
  geom_hline(yintercept = overall_sbd, linetype = "dashed", color = "#D9534F", linewidth = 0.7) +
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
    title    = "Stillbirth or 28-Day Death Stratified by WASH Exposure and Sex",
    subtitle = sprintf("Stillbirths + deaths by day 28 per 1,000 births with 95%% Wilson confidence intervals; n = infants in group; dashed line = overall rate (%.1f)", overall_sbd),
    x        = NULL,
    y        = "Stillbirths and deaths by day 28 (per 1,000 births)",
    fill     = "Infant Sex",
    alpha    = "Group size"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

show_plot(p_strat_sbd, "Barplot_stillbirth_or_death_strat_sbd", width = 14, height = 10)
