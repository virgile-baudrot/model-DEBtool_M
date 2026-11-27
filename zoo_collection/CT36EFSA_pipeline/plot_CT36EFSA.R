# plot_CT36EFSA.R -- figures of the CT36 EFSA report from the CSV files written by
# run_octave_CT36EFSA.m.
# Scenario "acute" (default) draws:
#   <sp>_growth_combined.png        top part of Figs 3.3-3.11 (acute, n-day window)
#                                   (Podarcis: <sp>_combined_dynamics.png + <sp>_combined_fecundity.png;
#                                    Rana: also <sp>_LdL_combined.png)
#   <sp>_comparison_variables.png   bottom part of Figs 3.3-3.11 (chronic outputs, from the outputs_*.csv)
#   <sp>_growth_combined_relative.png  same as growth_combined, in % of the factor-1 curve
# Scenario "full" draws the same growth figures from traj_*_full.csv, named <sp>_full_*.png
# (APPENDIX C, Figs C.1-C.12); no chronic-output figure (it does not depend on the scenario).
# Usage: Rscript plot_CT36EFSA.R <species> [outroot] [acute|full]

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
sp <- args[1]
outroot <- if (length(args) >= 2) args[2] else
  file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))),
            "..", "CT36EFSA_results")
scenario <- if (length(args) >= 3) args[3] else "acute"
stopifnot(scenario %in% c("acute", "full"))
is_full <- scenario == "full"
sfx <- if (is_full) "_full" else ""
dir <- file.path(outroot, sp)
sp_lab <- gsub("_", " ", sp)

modes <- switch(sp,
  Sorex_araneus = , Crocidura_russula = c("f_lactation", "f_juvenile", "f_adult"),
  Podarcis_muralis = c("f_embryo", "f_juvenile", "f_adult"),
  c("f_nest", "f_juvenile", "f_adult"))

# ---------- variable descriptions (from master_run.m)
desc <- c(
  ab = "age at birth", tx = "time since birth at fledging", tp = "time since birth at puberty",
  tR = "time since birth at 1st brood", am = "life span", Wwb = "wet weight at birth",
  Wwi = "ultimate wet weight", Wwim = "ultimate wet weight for males", Ww0 = "initial wet weight",
  Ri = "maximum reproductive rate", tg = "gestation time", Li = "ultimate total length",
  Wwx = "wet weight at weaning", Wwp = "wet weight at puberty",
  tj = "time since birth at metamorphosis", Lj = "SVL at metamorphosis",
  Lp = "SVL at puberty (females)", Lpm = "SVL at puberty (males)", Lb = "SVL at birth")
if (sp %in% c("Sorex_araneus", "Crocidura_russula")) desc["tx"] <- "time since birth at weaning"
if (sp %in% c("Rana_temporaria", "Ichthyosaura_alpestris", "Podarcis_muralis")) {
  desc[c("Li", "Lim", "Wwi", "Wwim")] <- c("ultimate SVL (females)", "ultimate SVL (males)",
                                           "ultimate wet weight (females)", "ultimate wet weight (males)")
}
if (sp == "Podarcis_muralis") desc["Ww0"] <- "wet weight of eggs"

theme_set(theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5)))

# ---------- growth curves (acute scenario)
dat <- read.csv(file.path(dir, sprintf("data_%s.csv", sp)))
traj <- do.call(rbind, lapply(modes, function(m) {
  d <- read.csv(file.path(dir, sprintf("traj_%s_%s%s.csv", sp, m, sfx))); d$mode <- m; d
}))

f_scale <- scale_colour_viridis_c(option = "turbo", direction = -1, limits = c(0.5, 1),
                                  breaks = seq(0.5, 1, 0.1), name = NULL)

growth_panel <- function(m, panel, xlab, ylab, title) {
  d <- traj[traj$mode == m & traj$panel == panel, ]
  d <- d[order(d$f_value, decreasing = TRUE), ]            # low f drawn on top
  obs <- dat[dat$panel == panel, ]
  ggplot(d, aes(x, y, group = f_value, colour = f_value)) +
    geom_line(data = d[abs(d$f_value - 1) > 1e-6, ], linewidth = 0.7, alpha = 0.7) +
    geom_line(data = d[abs(d$f_value - 1) < 1e-6, ], linewidth = 1.3) +
    { if (nrow(obs)) geom_point(data = obs, aes(x, y, shape = source), inherit.aes = FALSE,
                                colour = "black", size = 1.6) } +
    scale_shape_manual(values = c(data = 16, tL_K = 16, tL_D = 1), guide = "none") +
    f_scale + expand_limits(y = 0) + labs(x = xlab, y = ylab, title = title) +
    guides(colour = guide_colourbar(title = "reduction\nfactor", barheight = unit(5, "cm")))
}
# same panel, y = proportional change relative to the reference curve (factor 1) at the same time, in %
growth_panel_rel <- function(m, panel, ylab, title, xlab = "Time since birth (d)") {
  d <- traj[traj$mode == m & traj$panel == panel, ]
  ref <- d[abs(d$f_value - 1) < 1e-6, ]
  ref <- ref[order(ref$x), ]
  # reference value at the same x (time, or SVL for the weight-SVL panels), interpolated on the reference curve
  d$y <- 100 * (d$y / approx(ref$x, ref$y, d$x, ties = mean)$y - 1)
  d <- d[is.finite(d$y), ]
  ggplot(d, aes(x, y, group = f_value, colour = f_value)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_line(data = d[abs(d$f_value - 1) > 1e-6, ], linewidth = 0.7, alpha = 0.7) +
    geom_line(data = d[abs(d$f_value - 1) < 1e-6, ], linewidth = 1.3) +
    f_scale + labs(x = xlab, y = ylab, title = title) +
    guides(colour = guide_colourbar(title = "reduction\nfactor", barheight = unit(5, "cm")))
}
save <- function(p, name, w, h) {
  if (is_full) name <- sub(sp, paste0(sp, "_full"), name, fixed = TRUE)
  ggsave(file.path(dir, name), p, width = w, height = h, dpi = 200, bg = "white")
  cat("  saved", file.path(dir, name), "\n")
}
var_title <- function(m) sprintf("%s - %s variation", sp_lab, gsub("_", " ", m))

if ("tW" %in% traj$panel) {                                         # birds and shrews
  ps <- lapply(modes, function(m) growth_panel(m, "tW", "Time since birth (d)", "Wet weight (g)", var_title(m)))
  save(wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect"), sprintf("%s_growth_combined.png", sp), 14, 4)
  if (!is_full) save(ps[[2]], sprintf("%s_f_juvenile.png", sp), 5.5, 4)            # single panel (Fig 3.2 for Alauda)
  ps <- lapply(modes, function(m) growth_panel_rel(m, "tW", "Change in wet weight vs factor 1 (%)", var_title(m)))
  if (!is_full) save(wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect"), sprintf("%s_growth_combined_relative.png", sp), 14, 4)
} else if (sp == "Podarcis_muralis") {
  dyn <- list(c("tL_f", "Time since birth (d)", "SVL (cm)", "SVL vs Time (Female)"),
              c("tL_m", "Time since birth (d)", "SVL (cm)", "SVL vs Time (Male)"),
              c("LW_f", "SVL (cm)", "Wet weight (g)", "Weight vs SVL (Female)"),
              c("LW_m", "SVL (cm)", "Wet weight (g)", "Weight vs SVL (Male)"))
  ps <- list()
  for (m in modes[2:3]) for (q in dyn)
    ps[[length(ps) + 1]] <- growth_panel(m, q[1], q[2], q[3], sprintf("%s - %s", q[4], gsub("_", " ", m))) +
      theme(plot.title = element_text(size = 9))
  save(wrap_plots(ps, nrow = 2) + plot_layout(guides = "collect"), sprintf("%s_combined_dynamics.png", sp), 16, 6.5)
  ps <- list()
  for (m in modes[2:3]) for (q in dyn) {
    ylab <- if (q[1] %in% c("tL_f", "tL_m")) "Change in SVL vs factor 1 (%)" else "Change in weight at same SVL (%)"
    ps[[length(ps) + 1]] <- growth_panel_rel(m, q[1], ylab, sprintf("%s - %s", q[4], gsub("_", " ", m)), xlab = q[2]) +
      theme(plot.title = element_text(size = 9))
  }
  if (!is_full) save(wrap_plots(ps, nrow = 2) + plot_layout(guides = "collect"), sprintf("%s_combined_dynamics_relative.png", sp), 16, 6.5)
  ps <- lapply(modes, function(m) growth_panel(m, "LN", "SVL (cm)", "Yearly fecundity (#)",
                                               sprintf("Fecundity vs SVL (Female) - %s", gsub("_", " ", m))))
  save(wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect"), sprintf("%s_combined_fecundity.png", sp), 14, 4)
} else {                                                             # amphibians: female row, male row
  ps <- c(lapply(modes, function(m) growth_panel(m, "tL_f", "Time since birth (d)", "SVL (cm)",
                                                 sprintf("%s - female", gsub("_", " ", m)))),
          lapply(modes, function(m) growth_panel(m, "tL_m", "Time since birth (d)", "SVL (cm)",
                                                 sprintf("%s - male", gsub("_", " ", m)))))
  save(wrap_plots(ps, nrow = 2) + plot_layout(guides = "collect") +
         plot_annotation(title = sp_lab, theme = theme(plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5))),
       sprintf("%s_growth_combined.png", sp), 14, 7.5)
  ps <- c(lapply(modes, function(m) growth_panel_rel(m, "tL_f", "Change in SVL vs factor 1 (%)",
                                                     sprintf("%s - female", gsub("_", " ", m)))),
          lapply(modes, function(m) growth_panel_rel(m, "tL_m", "Change in SVL vs factor 1 (%)",
                                                     sprintf("%s - male", gsub("_", " ", m)))))
  if (!is_full) save(wrap_plots(ps, nrow = 2) + plot_layout(guides = "collect") +
         plot_annotation(title = sp_lab, theme = theme(plot.title = element_text(face = "bold.italic", size = 15, hjust = 0.5))),
       sprintf("%s_growth_combined_relative.png", sp), 14, 7.5)
  if ("LdL" %in% traj$panel) {
    ps <- lapply(modes, function(m) growth_panel(m, "LdL", "SVL (cm)", "Change in SVL (cm/d)",
                                                 sprintf("%s LdL - %s", sp_lab, gsub("_", " ", m))))
    save(wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect"), sprintf("%s_LdL_combined.png", sp), 14, 4)
  }
}

# ---------- chronic outputs vs reduction factor (from outputs_*.csv)
if (is_full) quit(save = "no")
out <- do.call(rbind, lapply(modes, function(m) {
  d <- read.csv(file.path(dir, sprintf("outputs_%s_%s.csv", sp, m)))
  vars <- setdiff(names(d), "f_value")
  data.frame(f_value = rep(d$f_value, length(vars)), var = rep(vars, each = nrow(d)),
             value = unlist(d[vars]), mode = m)
}))
vars <- unique(out$var)
out$var <- factor(out$var, levels = vars,
                  labels = sprintf("%s\n(%s)", ifelse(is.na(desc[vars]), vars, desc[vars]), vars))
out$mode <- factor(out$mode, levels = modes)
cols <- setNames(c("#D95319", "#77AC30", "#0072BD"), modes)       # MATLAB orange, green, blue
shps <- setNames(c(16, 15, 17), modes)
p <- ggplot(out, aes(f_value, value, colour = mode, shape = mode)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(size = 1.4, na.rm = TRUE) +
  facet_wrap(~var, nrow = 2, scales = "free_y") + expand_limits(y = 0) +
  scale_colour_manual(values = cols, name = NULL) + scale_shape_manual(values = shps, name = NULL) +
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "Reduction factor applied to the baseline functional response", y = NULL) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9), axis.title.x = element_text(face = "bold"))
save(p, sprintf("%s_comparison_variables.png", sp), 2.6 * ceiling(length(vars) / 2) + 1, 5.5)
