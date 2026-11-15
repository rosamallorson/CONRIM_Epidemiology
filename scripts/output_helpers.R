# Output folders and helpers shared by all analysis scripts.
# Plots go to output/plots/, tables to output/tables/.

# Force a UTF-8 locale so characters such as the en dash are not written as "<U+2013>"
# (Rscript/quarto often start in the "C" locale).
for (loc in c("en_US.UTF-8", "C.UTF-8")) {
  if (l10n_info()$`UTF-8` || nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
}

output_dir <- "output"
plots_dir  <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Save a ggplot to output/plots/<name>.png then print it.
show_plot <- function(p, name, width = 12, height = 8) {
  ggplot2::ggsave(file.path(plots_dir, paste0(name, ".png")), p,
                  width = width, height = height, dpi = 300, bg = "white")
  print(p)
}

# Write a table to output/tables/<name>.csv (UTF-8) and return it invisibly.
save_table <- function(df, name) {
  utils::write.csv(df, file.path(tables_dir, paste0(name, ".csv")), row.names = FALSE,
                   fileEncoding = "UTF-8")
  invisible(df)
}
