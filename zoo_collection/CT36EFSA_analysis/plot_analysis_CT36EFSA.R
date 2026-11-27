# plot_analysis_CT36EFSA.R -- figures of report Chapter 2 (Figs 2.2-2.19): sensitivity of
# the Add-my-Pet 0-variate predictions to the scaled functional response f.
# (Same computations as the former dataAnalysis.Rmd, saved as PNG files.)
#
# Inputs (this folder):
#   variation_f_food - <sp>.csv  f, data, prd (prediction at f), prd_f1 (prediction at f = 1),
#                                symbol, units, description
#   ref_values - <sp>.csv        symbol, units, description, min, max (range of observed values)
# Outputs (default: this folder):
#   <sp>_prediction.png   predictions vs % of food reduction, data range (green band),
#                         AmP data value (grey line)      -> Figs 2.2, 2.4, ..., 2.18
#   <sp>_elasticity.png   % change of each variable vs % of food reduction, linear fit
#                         through 0 (red, slope = elasticity) and quadratic fit (grey)
#                                                         -> Figs 2.3, 2.5, ..., 2.19
#   elasticity_CT36EFSA.csv  slopes and R2 of both fits, all species
# Usage: Rscript plot_analysis_CT36EFSA.R [outdir]   (run from any folder)

suppressPackageStartupMessages({ library(ggplot2) })

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else here
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# species in report order (Figs 2.2/2.3 = Perdix, ..., 2.18/2.19 = Podarcis)
species <- c("Perdix_perdix", "Motacilla_flava", "Alauda_arvensis", "Turdus_philomelos",
             "Sorex_araneus", "Crocidura_russula", "Ichthyosaura_alpestris",
             "Rana_temporaria", "Podarcis_muralis")

theme_set(theme_minimal(base_size = 11) + theme(strip.text = element_text(size = 9)))
wrap <- label_wrap_gen(width = 25)

fits <- list()
for (sp in species) {
  v <- read.csv(file.path(here, sprintf("variation_f_food - %s.csv", sp)))
  r <- read.csv(file.path(here, sprintf("ref_values - %s.csv", sp)))
  v$prd <- as.numeric(v$prd); v$data <- as.numeric(v$data)
  lev <- unique(paste0(v$description, ", ", v$symbol))   # panel order = order in the input file
  # min-max of the data, matched on the symbol only (the descriptions differ between files)
  v <- merge(v, r[, c("symbol", "min", "max")], by = "symbol", all.x = TRUE, sort = FALSE)
  v$panel <- paste0(v$description, ", ", v$symbol)
  v$panel <- factor(v$panel, levels = lev)
  v$red <- (1 - v$f) * 100                     # % of food reduction
  v$eff <- (1 - v$prd / v$prd_f1) * 100        # % of reduction in the variable
  n_panel <- nlevels(v$panel)
  width <- if (n_panel > 10) 12 else 10

  # ---- predictions vs data range
  g1 <- ggplot(v) +
    geom_ribbon(aes(red, ymin = min, ymax = max), alpha = 0.5, fill = "#C8E897") +
    geom_line(aes(red, data), colour = "grey30") +
    geom_point(aes(red, prd), na.rm = TRUE) +
    scale_y_continuous(limits = c(0, NA)) +
    facet_wrap(~panel, labeller = wrap, nrow = 2, scales = "free") +
    labs(x = "Percentage of food reduction",
         y = "Prediction (black points)\nData min-max (green band)")
  ggsave(file.path(outdir, sprintf("%s_prediction.png", sp)), g1,
         width = width, height = 5, dpi = 200, bg = "white")

  # ---- elasticity: linear fit through 0 (slope) and quadratic fit through 0
  fit <- do.call(rbind, lapply(split(v[is.finite(v$eff), ], ~panel, drop = TRUE), function(d) {
    lin <- lm(eff ~ 0 + red, data = d)
    quad <- lm(eff ~ 0 + red + I(red^2), data = d)
    data.frame(species = sp, panel = d$panel[1], symbol = d$symbol[1], n = nrow(d),
               elasticity = unname(coef(lin)[1]), R2_lin = summary(lin)$r.squared,
               quad_lin = unname(coef(quad)[1]), quad_sq = unname(coef(quad)[2]),
               R2_quad = summary(quad)$r.squared)
  }))
  fits[[sp]] <- fit
  g2 <- ggplot(v[is.finite(v$eff), ], aes(red, eff)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ 0 + x + I(x^2), se = TRUE,
                colour = "grey20", fill = "lightgray") +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = TRUE,
                colour = "#aa1122", fill = "lightgray") +
    geom_text(data = fit, aes(x = -Inf, y = Inf, label = sprintf("elasticity = %.3f", ifelse(abs(elasticity) < 5e-4, 0, elasticity))),
               hjust = -0.05, vjust = 1.3, size = 3, colour = "#aa1122") +
    facet_wrap(~panel, labeller = wrap, nrow = 2, scales = "free") +
    labs(x = "Percentage of food reduction", y = "Percentage of reduction in variable")
  ggsave(file.path(outdir, sprintf("%s_elasticity.png", sp)), g2,
         width = width, height = 5, dpi = 200, bg = "white")
  cat("  saved", sp, "\n")
}
fits <- do.call(rbind, fits)
fits$panel <- as.character(fits$panel)
write.csv(fits, file.path(outdir, "elasticity_CT36EFSA.csv"), row.names = FALSE)
