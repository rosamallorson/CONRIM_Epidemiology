# CONRIM_Epidemiology

Stillbirth, 28-day neonatal mortality and household water, sanitation and hygiene (WASH) in Busoga (Jinja District, Uganda), from the CONRIM community surveillance data.

## Data

[scripts/Load_data.R](scripts/Load_data.R) uses the `busoga` package (`prepare_survey_data()`) to download the CONRIMSHORTFRM (first visit) and CONRIMFOLLOW (28-day follow-up) forms from ODK Central, then derives the analysis data sets.

- Needs `ODK_EMAIL`, `ODK_PASSWORD` and `BUSOGA_ROOT` (the `busoga` project folder) in the environment, for example in `~/.Renviron`.
- To reuse the cached CSVs instead of downloading, run `options(conrim.download = FALSE)` before sourcing.
- Reports with inconsistent vital status or dates are dropped and logged in `output/tables/inconsistent_records.csv`.

## Running

All R scripts are in `scripts/`. Run them from the project root (this folder), for example `source("scripts/Forest_NMR.R")`. Every script sources `scripts/Load_data.R` if the data are not already in memory.

| Script (in `scripts/`) | Output |
|---|---|
| `Forest_stillbirth.R` | Stillbirth rate by stratum, with Firth odds ratios |
| `Forest_NMR.R` | Neonatal mortality rate (NMR) by stratum and sex |
| `Forest_stillbirth_or_death.R` | Stillbirth-or-28-day-death rate by stratum and sex |
| `Forest_plots.R` | Univariate WASH odds ratios (Wald) |
| `Forest_OR_firth.R` | Firth odds ratios, NMR outcome vs stillbirth-or-death outcome |
| `Barplot_NMR.R`, `Barplot_stillbirth_or_death.R` | Rates by WASH exposure and infant sex |
| `WASH_stillbirth_NMR_report.qmd` | HTML and PPTX report |

Render the report with a UTF-8 locale, or en dashes can print as `<U+2013>`:

```
LC_ALL=en_US.UTF-8 quarto render WASH_stillbirth_NMR_report.qmd --output-dir output/report
```

## Outputs

Everything is written (relative to the project root) to `output/`: `plots/` (PNG), `tables/` (CSV), `report/` (HTML and PPTX) and `logs/`. Shared helpers are in [scripts/output_helpers.R](scripts/output_helpers.R), which also forces a UTF-8 locale. [scripts/forest_helpers.R](scripts/forest_helpers.R) holds the odds-ratio plotting functions.

## Results from the latest run (fresh ODK download, 15 Nov 2026)

All estimates are unadjusted and univariate, with many comparisons and no multiplicity correction. Treat them as hypothesis-generating, not causal.

**Rates**

| Measure | Value |
|---|---|
| Births (stillbirths and live births with consistent data) | 4,906 |
| Stillbirths | 44 |
| Neonatal deaths (live births, by day 28) | 104 |
| Stillbirth or death by day 28, per 1,000 births | 30.2 |
| NMR, per 1,000 live births | 21.4 |

**Stillbirth** (Firth OR, 44 events)
- Higher odds with no folic acid in pregnancy (OR 2.48, 95% CI 1.18-4.80) and with a multiple pregnancy (OR 2.62, 1.02-5.74).
- No clear difference by infant sex or district.
- Jinja South Division had lower odds than the Jinja North Division reference (OR 0.24, 0.03-0.94). This is one of 8 subcounty contrasts with small counts.

**WASH and neonatal death** (Firth OR; the live-birth NMR outcome and the stillbirth-or-death outcome give identical WASH estimates, because no stillbirth has WASH answers)

| Exposure | OR (95% CI) |
|---|---|
| Household latrine (vs none) | 0.20 (0.10-0.44) |
| Soap or sanitiser for handwashing (vs water only) | 0.24 (0.14-0.45) |
| Handwashing before feeding baby | 0.37 (0.22-0.61) |
| Handwashing at all 3 critical times | 0.49 (0.27-0.82) |
| Protected dry-season water source (vs unprotected) | 0.24 (0.12-0.54) |
| Dry-season shallow well (vs borehole) | 6.06 (2.46-13.24) |
| Protected rainy-season water source | 0.85 (0.29-4.12) |
| Animals in the household | 0.73 (0.47-1.13) |

- Some groups are very small. The "no latrine" group has 111 infants, "water only" handwashing has 241, and several water-source categories have fewer than 100. Wide intervals and a few estimates with zero events reflect this.
- The reference groups for the multi-category water source models are the largest category (tap water for rainy season, borehole for dry season).
- Some categories had unestimable intervals and were dropped from the Wald plots (for example sanitiser, BCG and breastmilk as pre-lacteal items).

## Known limitations

- Stillbirths have no WASH answers (the questions appear to be skipped), so WASH analyses cover live births only.
- Infants without a 28-day follow-up and without a recorded death are excluded from the 28-day outcome.
- Harvested rainwater is counted as a protected source (assumed improved).
- `geom_errorbarh()` is deprecated in ggplot2 4.0 and prints a warning. It does not affect results.
