# Shared helpers for the odds-ratio forest plots
# (sourced by Forest_plots.R and Forest_OR_firth.R).
# Needs the packages loaded by Load_data.R; method = "firth" also needs logistf.

# Attach events / N for each exposure level and its reference group.
# df_sub holds the outcome as `y` and the exposure as `exposure`.
add_group_counts <- function(res, df_sub, ref_level) {
  counts <- df_sub %>%
    group_by(exposure = as.character(exposure)) %>%
    summarise(deaths = sum(y), n = n(), .groups = "drop")
  ref <- filter(counts, exposure == ref_level)

  res %>%
    left_join(counts, by = c("level" = "exposure")) %>%
    mutate(ref_deaths = ref$deaths, ref_n = ref$n)
}

# Univariate logistic regression of `outcome` on each variable in `vars`.
# Reference group: `ref[[var]]` if given, else "No", else the largest group.
# Groups with fewer than `min_n` rows are dropped.
# method = "wald" (glm, Wald CI) or "firth" (penalized likelihood, profile CI).
fit_univariate_or <- function(data, vars, outcome, method = c("wald", "firth"),
                              ref = NULL, min_n = 10) {
  method <- match.arg(method)

  map_dfr(vars, function(v) {
    df_sub <- data %>%
      select(y = all_of(outcome), exposure = all_of(v)) %>%
      filter(!is.na(exposure), exposure != "") %>%
      group_by(exposure) %>%
      filter(n() >= min_n) %>%
      ungroup()

    if (n_distinct(df_sub$exposure) < 2) return(NULL)

    counts  <- sort(table(df_sub$exposure), decreasing = TRUE)
    ref_cat <- if (!is.null(ref) && v %in% names(ref)) {
      ref[[v]]
    } else if ("No" %in% names(counts)) {
      "No"
    } else {
      names(counts)[1]
    }
    df_sub$exposure <- relevel(factor(df_sub$exposure), ref = ref_cat)

    if (method == "wald") {
      res <- glm(y ~ exposure, data = df_sub, family = binomial()) %>%
        tidy() %>%
        filter(term != "(Intercept)") %>%
        transmute(
          level     = str_remove(term, "^exposure"),
          or_est    = exp(estimate),
          conf.low  = exp(estimate - 1.96 * std.error),
          conf.high = exp(estimate + 1.96 * std.error),
          p.value
        )
    } else {
      fit <- logistf::logistf(y ~ exposure, data = df_sub)
      res <- tibble(
        level     = sub("^exposure", "", names(fit$coefficients))[-1],
        or_est    = exp(fit$coefficients[-1]),
        conf.low  = exp(fit$ci.lower[-1]),
        conf.high = exp(fit$ci.upper[-1]),
        p.value   = unname(fit$prob[-1])
      )
    }

    res %>%
      mutate(variable = v, label = paste0(v, ": ", level, " (vs. ", ref_cat, ")")) %>%
      add_group_counts(df_sub, ref_cat)
  })
}

# Drop estimates whose CI is unestimable (separation / empty cells -> 0 to Inf)
drop_unestimable <- function(res, max_ci_ratio = 1000) {
  keep <- is.finite(res$conf.low) & is.finite(res$conf.high) &
    res$conf.low > 0 & (res$conf.high / res$conf.low) < max_ci_ratio
  if (any(!keep)) {
    message("Excluded (unestimable CI): ", paste(res$label[!keep], collapse = "; "))
  }
  res[keep, ]
}

# Add WASH domain (factor) and short row labels to multi-category results
group_wash_domains <- function(res) {
  res %>%
    mutate(
      domain = case_when(
        variable == "eat_drink_cat" ~ "Pre-lacteal / Eat & Drink",
        variable == "rainy_water" ~ "Water Source (Rainy)",
        variable == "dry_water" ~ "Water Source (Dry)",
        variable == "handwash_agent" ~ "Handwashing Agent",
        str_starts(variable, "hw_before_") ~ "Handwashing Timing",
        variable == "latrine_status" ~ "Sanitation / Latrine",
        str_starts(variable, "animal_") ~ "Household Animals",
        TRUE ~ variable
      ),
      domain = factor(domain, levels = c(
        "Water Source (Rainy)",
        "Water Source (Dry)",
        "Sanitation / Latrine",
        "Handwashing Agent",
        "Handwashing Timing",
        "Household Animals",
        "Pre-lacteal / Eat & Drink"
      )),
      clean_label = case_when(
        variable == "animal_poultry"     ~ "Poultry: Yes (vs. No)",
        variable == "animal_goats"       ~ "Goats: Yes (vs. No)",
        variable == "animal_cows"        ~ "Cows: Yes (vs. No)",
        variable == "hw_before_eating"   ~ "Before eating: Yes (vs. No)",
        variable == "hw_before_handling" ~ "Before handling: Yes (vs. No)",
        variable == "hw_before_feeding"  ~ "Before feeding baby: Yes (vs. No)",
        TRUE ~ str_remove(label, "^[^:]+: ")
      )
    )
}

# Forest plot with label, OR (95% CI) and events/N columns on the left.
# Columns are built into monospaced y-axis text so they stay aligned (incl. facets).
forest_table_plot <- function(res, label_col = "label", facet_col = NULL,
                              title, subtitle, point_color = "#1F78B4",
                              point_size = 2.8, base_size = 11, xlim = NULL,
                              events_label = "Events") {
  res <- drop_unestimable(res) %>%
    mutate(
      row_text = .data[[label_col]],
      or_text  = sprintf("%.2f (%.2f-%.2f)", or_est, conf.low, conf.high),
      n_text   = sprintf("%d/%d vs %d/%d", deaths, n, ref_deaths, ref_n)
    )

  hdr <- c("Exposure", "OR (95% CI)", paste0(events_label, "/N (group vs ref)"))
  w <- c(max(nchar(c(res$row_text, hdr[1]))),
         max(nchar(c(res$or_text, hdr[2]))),
         max(nchar(c(res$n_text, hdr[3]))))
  fmt <- function(a, b, c) {
    paste(str_pad(a, w[1], "right"), str_pad(b, w[2], "right"), str_pad(c, w[3], "right"), sep = "  ")
  }

  res <- res %>% mutate(row_label = fmt(row_text, or_text, n_text))
  header_row <- tibble::tibble(row_label = fmt(hdr[1], hdr[2], hdr[3]))

  if (!is.null(facet_col)) {
    fac_levels <- levels(droplevels(res[[facet_col]]))
    header_row[[facet_col]] <- factor(fac_levels[1], levels = fac_levels)
    res[[facet_col]] <- droplevels(res[[facet_col]])
  }

  plot_df <- bind_rows(res, header_row) %>%
    mutate(row_label = factor(row_label, levels = c(row_label[!is.na(or_est)][order(or_est[!is.na(or_est)])],
                                                   header_row$row_label)))

  # Symmetric limits on the log scale so OR = 1 sits in the middle
  if (is.null(xlim)) {
    m <- max(res$conf.high, 1 / res$conf.low)
    xlim <- c(1 / m, m)
  }

  g <- ggplot(plot_df, aes(x = or_est, y = row_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#D9534F", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25,
                   color = "#2C3E50", na.rm = TRUE) +
    geom_point(size = point_size, color = point_color, na.rm = TRUE) +
    scale_x_log10(limits = xlim) +
    labs(title = title, subtitle = subtitle, x = "Odds Ratio (95% CI, log scale)", y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(family = "mono", hjust = 0, size = rel(0.85))
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
