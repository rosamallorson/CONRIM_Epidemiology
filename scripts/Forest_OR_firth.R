# Odds-ratio forest plots with Firth penalized logistic regression (profile-likelihood CIs).
# Compares two outcomes:
#   - true NMR outcome: death by day 28 among live births (stillbirths excluded), nmr_wash
#   - stillbirth or death by day 28 (stillbirths included), analysis_wash
# Run from the CONRIM project folder. Needs the logistf package.
source("scripts/output_helpers.R")


library(logistf)

# Load and clean the data once per session (skipped if already in memory;
# run rm(nmr_wash_cats) or source("scripts/Load_data.R") to force a reload)
if (!exists("nmr_wash_cats") || !"sb_death_28d" %in% names(analysis_wash)) source("scripts/Load_data.R")
source("scripts/forest_helpers.R")

sex_ref  <- c(infant_sex = "Female")
bin_vars <- c(clean_wash_vars, "infant_sex")

# Firth fits for both outcomes
or_sbd_bin <- fit_univariate_or(analysis_wash,     bin_vars,       "sb_death_28d",       "firth", ref = sex_ref)
or_nmr_bin <- fit_univariate_or(nmr_wash,          bin_vars,       "neonatal_death_28d", "firth", ref = sex_ref)
or_sbd_cat <- fit_univariate_or(analysis_wash_cats, multi_cat_vars, "sb_death_28d",       "firth")
or_nmr_cat <- fit_univariate_or(nmr_wash_cats,      multi_cat_vars, "neonatal_death_28d", "firth")

# --- True NMR outcome: forest plots ---
p_nmr_or <- forest_table_plot(
  or_nmr_bin,
  title        = "WASH Factors and Infant Sex vs. Neonatal Death by Day 28 (Stillbirths Excluded)",
  subtitle     = "Univariate odds ratios (95% Firth penalized-likelihood CI) on log scale",
  point_size   = 3.2,
  events_label = "Deaths"
)

show_plot(p_nmr_or, "Forest_OR_firth_nmr_or", width = 12, height = 8)

p_nmr_or_grouped <- forest_table_plot(
  group_wash_domains(or_nmr_cat),
  label_col    = "clean_label",
  facet_col    = "domain",
  title        = "WASH Associations with Neonatal Death by Day 28 (Stillbirths Excluded)",
  subtitle     = "Grouped by WASH domain; univariate odds ratios (95% Firth CI) on log scale vs. largest category",
  events_label = "Deaths"
)

show_plot(p_nmr_or_grouped, "Forest_OR_firth_nmr_or_grouped", width = 12, height = 11)

# --- Comparison of the two outcomes ---
sbd_name <- "Stillbirth or death by day 28"
nmr_name <- "Neonatal death (stillbirths excluded)"

compare_or <- function(sbd, nmr, label_col = "label", facet_col = NULL) {
  keep_cols <- c("variable", "level", label_col, facet_col)
  one <- function(d, suffix) {
    d %>%
      drop_unestimable() %>%
      select(all_of(keep_cols), or_est, conf.low, conf.high, p.value, deaths, n, ref_deaths, ref_n) %>%
      rename_with(~ paste0(.x, "_", suffix), c(or_est, conf.low, conf.high, p.value, deaths, n, ref_deaths, ref_n))
  }
  inner_join(one(nmr, "nmr"), one(sbd, "sbd")[, setdiff(names(one(sbd, "sbd")), c(label_col, facet_col))],
             by = c("variable", "level"))
}

cmp_bin <- compare_or(or_sbd_bin, or_nmr_bin)
cmp_cat <- compare_or(group_wash_domains(or_sbd_cat), group_wash_domains(or_nmr_cat),
                      label_col = "clean_label", facet_col = "domain")

# Side-by-side forest plot: both outcomes per row, with both ORs in the left columns
compare_forest <- function(cmp, label_col = "label", facet_col = NULL, title, subtitle, base_size = 11) {
  cmp <- cmp %>%
    mutate(
      row_text = .data[[label_col]],
      nmr_text = sprintf("%.2f (%.2f-%.2f)", or_est_nmr, conf.low_nmr, conf.high_nmr),
      sbd_text = sprintf("%.2f (%.2f-%.2f)", or_est_sbd, conf.low_sbd, conf.high_sbd)
    )
  hdr <- c("Exposure", "Neonatal death OR", "Stillbirth/death OR")
  w <- c(max(nchar(c(cmp$row_text, hdr[1]))), max(nchar(c(cmp$nmr_text, hdr[2]))), max(nchar(c(cmp$sbd_text, hdr[3]))))
  fmt <- function(a, b, c) paste(str_pad(a, w[1], "right"), str_pad(b, w[2], "right"), str_pad(c, w[3], "right"), sep = "  ")
  cmp <- cmp %>% mutate(row_label = fmt(row_text, nmr_text, sbd_text))
  header_row <- tibble::tibble(row_label = fmt(hdr[1], hdr[2], hdr[3]))

  long <- bind_rows(
    cmp %>% transmute(row_label, !!!if (!is.null(facet_col)) rlang::syms(facet_col), outcome = nmr_name,
                      or = or_est_nmr, lo = conf.low_nmr, hi = conf.high_nmr),
    cmp %>% transmute(row_label, !!!if (!is.null(facet_col)) rlang::syms(facet_col), outcome = sbd_name,
                      or = or_est_sbd, lo = conf.low_sbd, hi = conf.high_sbd)
  )

  if (!is.null(facet_col)) {
    fac_levels <- levels(droplevels(cmp[[facet_col]]))
    header_row[[facet_col]] <- factor(fac_levels[1], levels = fac_levels)
    long[[facet_col]] <- droplevels(long[[facet_col]])
  }

  ord <- cmp$row_label[order(cmp$or_est_nmr)]
  long <- bind_rows(long, header_row) %>%
    mutate(
      row_label = factor(row_label, levels = c(ord, header_row$row_label)),
      outcome   = factor(outcome, levels = c(nmr_name, sbd_name))
    )

  m <- max(long$hi, 1 / long$lo, na.rm = TRUE)

  g <- ggplot(long, aes(x = or, y = row_label, color = outcome)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#D9534F", linewidth = 0.8) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.4,
                  position = position_dodge(width = 0.6, orientation = "y"), na.rm = TRUE) +
    geom_point(size = 2.2, position = position_dodge(width = 0.6, orientation = "y"), na.rm = TRUE) +
    scale_x_log10(limits = c(1 / m, m)) +
    scale_color_discrete(na.translate = FALSE) +
    labs(title = title, subtitle = subtitle, x = "Odds Ratio (95% Firth CI, log scale)", y = NULL, color = "Outcome") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(family = "mono", hjust = 0, size = rel(0.85)),
      legend.position  = "bottom"
    )

  if (!is.null(facet_col)) {
    g <- g +
      facet_grid(reformulate(".", facet_col), scales = "free_y", space = "free_y") +
      theme(
        strip.text.y     = element_text(angle = 0, face = "bold", hjust = 0),
        strip.background = element_rect(fill = "#F0F4F8", color = NA),
        panel.spacing.y  = unit(0.6, "lines")
      )
  }
  g
}

p_or_compare <- compare_forest(
  cmp_bin,
  title    = "Binary WASH Factors and Infant Sex: Neonatal Death vs. Stillbirth or Death",
  subtitle = "Univariate odds ratios (95% Firth CI) on log scale; stillbirths have no WASH answers"
)

show_plot(p_or_compare, "Forest_OR_firth_or_compare", width = 14, height = 10)

p_or_compare_grouped <- compare_forest(
  cmp_cat,
  label_col = "clean_label",
  facet_col = "domain",
  title     = "WASH Categories: Neonatal Death vs. Stillbirth or Death",
  subtitle  = "Grouped by WASH domain; univariate odds ratios (95% Firth CI) on log scale"
)

show_plot(p_or_compare_grouped, "Forest_OR_firth_or_compare_grouped", width = 14, height = 12)

# Comparison table (binary factors and categories), rows where the two outcomes differ most first
or_comparison <- bind_rows(
  cmp_bin %>% mutate(group = "Binary factors"),
  cmp_cat %>% mutate(group = "Categories", label = clean_label)
) %>%
  transmute(
    group, label,
    nmr_deaths_n = sprintf("%d/%d vs %d/%d", deaths_nmr, n_nmr, ref_deaths_nmr, ref_n_nmr),
    nmr_or       = sprintf("%.2f (%.2f-%.2f)", or_est_nmr, conf.low_nmr, conf.high_nmr),
    nmr_p        = signif(p.value_nmr, 2),
    sbd_events_n = sprintf("%d/%d vs %d/%d", deaths_sbd, n_sbd, ref_deaths_sbd, ref_n_sbd),
    sbd_or       = sprintf("%.2f (%.2f-%.2f)", or_est_sbd, conf.low_sbd, conf.high_sbd),
    sbd_p        = signif(p.value_sbd, 2),
    log_or_diff  = abs(log(or_est_nmr) - log(or_est_sbd))
  ) %>%
  arrange(desc(log_or_diff))

options(width = 250)
or_comparison %>% select(-log_or_diff) %>% as.data.frame() %>% save_table("firth_or_comparison_nmr_vs_sbd") %>% print()
