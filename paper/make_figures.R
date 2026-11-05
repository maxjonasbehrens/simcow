# Figs 2 and 3 of the paper, from the simulation results in paper/results/.
# Fig 2: cross-validated J_g (= mean weight of subgroup g) vs shift, all external subgroups;
#        individual weights for one grid cell, recomputed with the seeds of run_simulation.R.
# Fig 3: RMSE per run by scenario, similarity, and ESS ratio, all 32,400 runs.
# Usage (from any directory): Rscript paper/make_figures.R
# Needs the replayed data sets in paper/output/cache/ (written by run_simulation.R) for Fig 2B.
# Writes paper/output/figures/Fig2 and Fig3 (.pdf, .tif). Fig 1: make_fig1.R; Fig 4: make_fig4.R.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(gridExtra) })
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "..", "R", "create_cow.R"))
figdir <- file.path(here, "output", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
res <- file.path(here, "results")
pdf(NULL)   # ggplotGrob() opens a device; keep it from writing Rplots.pdf to the working directory

samples <- c(local = "#1baf7a", global = "#eb6834", D = "#2a78d6", IW = "#eda100", mixed = "#6250d6")
sample_labels <- c(local = "local", global = "global", D = "weighted", IW = "covariate-shift IW", mixed = "mixed")
scenario_labels <- c(AY = "outcome", AX = "covariate", AH = "covariate + outcome")
# PLOS figure files: Arial, Times, or Symbol at 8-12 pt; base_size 10 keeps theme_bw's
# smallest text (axis, legend, strips: rel 0.8) at 8 pt. Width 7.5 in, TIFF 600 dpi, LZW, RGB.
theme_paper <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom")
save_fig <- function(p, file, width, height) {
  ggsave(paste0(file, ".pdf"), p, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(file, ".tif"), p, width = width, height = height, units = "in", dpi = 600,
         device = ragg::agg_tiff, compression = "lzw")
}

# ---------------- Fig 2 ----------------
sub <- read.csv(file.path(res, "sim_subgroups.csv")) |>
  mutate(scenario = factor(scenario_labels[scenario], scenario_labels))
jg <- sub |> group_by(scenario, n_target, shift) |>
  summarise(med = median(J_cf), q25 = quantile(J_cf, .25), q75 = quantile(J_cf, .75), .groups = "drop")
p2a <- ggplot(jg, aes(shift, med, color = factor(n_target), fill = factor(n_target))) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(. ~ scenario) +
  scale_color_manual(values = c("10" = "#6baed6", "20" = "#2171b5", "30" = "#08306b"),
                     aesthetics = c("color", "fill"), name = expression(n[target])) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Shift of external subgroup", y = expression("Subgroup weight" ~ italic(J)[g])) +
  theme_paper + theme(legend.position = "top")

# Individual weights for one grid cell: 9 external subgroups (shift 1/3 to 3), n_g = n_target = 20
cache <- file.path(here, "output", "cache", "replay_kept.rds")
if (!file.exists(cache)) stop("Replayed data sets missing: run paper/run_simulation.R first (it writes ", cache, ")")
kept <- readRDS(cache)
param_grid <- expand.grid(scenario = c("AX","AY","AH"), n_cat = c(2,4,8,10), n_train = c(10,20,30),
                          n_target = c(10,20,30), to = c(3,2,1), num_env_predictors = c(3),
                          stabilization_method = c("none","truncate","stabilize","calibrate"),
                          stringsAsFactors = FALSE)
covs <- c("X_global_1", paste0("X_env_", 1:3))
cell <- which(sapply(kept, function(x) with(param_grid[x$i, ], n_cat == 10 && to == 3 && n_train == 20 && n_target == 20)))
w_cell <- lapply(cell, function(k) {
  x <- kept[[k]]; cw <- create_cow(x$tr, "A", 1, covs, outcome = "Y", seed = 1e6 + k)   # seed as in run_simulation.R
  ext <- x$tr$A != "1"
  list(w = data.frame(scenario = param_grid$scenario[x$i], group = as.character(x$tr$A[ext]), w = cw$weight[ext]),
       info = cbind(row = x$i, sim = x$s, group = as.integer(attr(cw, "subgroups")$cohort), J = attr(cw, "subgroups")$J))
})
w_ind <- do.call(rbind, lapply(w_cell, `[[`, "w"))
info_cell <- as.data.frame(do.call(rbind, lapply(w_cell, `[[`, "info")))
chk <- merge(info_cell, sub, by = c("row", "sim", "group"))
stopifnot(nrow(chk) == nrow(info_cell), max(abs(chk$J - chk$J_cf)) < 1e-12)
cat("Fig 2 check: recomputed J_g of", nrow(chk), "subgroups equals paper/results/sim_subgroups.csv\n")
w_ind <- w_ind |> mutate(shift = (as.numeric(group) - 1) / 9 * 3,
                         scenario = factor(scenario_labels[scenario], scenario_labels))
y_max <- w_ind |> group_by(scenario, shift) |> summarise(u = boxplot.stats(w)$stats[5], .groups = "drop") |> pull(u) |> max()
p2b <- ggplot(w_ind, aes(shift, w, group = shift)) +
  geom_boxplot(fill = samples[["D"]], alpha = 0.35, color = "grey25", outlier.shape = NA, width = 0.22, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6, fill = "white", color = "grey10") +
  facet_grid(. ~ scenario) +
  scale_x_continuous(breaks = 0:3) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = "Shift of external subgroup", y = expression("Individual weight" ~ w[j])) +
  theme_paper
g2 <- rbind(ggplotGrob(p2a + labs(tag = "A")), ggplotGrob(p2b + labs(tag = "B")), size = "first")
save_fig(g2, file.path(figdir, "Fig2"), 7.5, 5.6)

# ---------------- Fig 3 ----------------
# Bins: ESS ratio of the weighted sample, rounded; a ">= 5" bin keeps every run.
runs <- read.csv(file.path(res, "sim_runs.csv"))
long <- runs |>
  mutate(ess_bin = factor(pmin(pmax(round(ess_D), 1), 5), 1:5)) |>
  select(row, sim, scenario, to, ess_bin, all_of(paste0("rmse_", names(samples)))) |>
  pivot_longer(starts_with("rmse_"), names_to = "sample", names_prefix = "rmse_", values_to = "RMSE") |>
  filter(!is.na(RMSE)) |>
  mutate(sample = factor(sample, names(samples)),
         scenario = factor(scenario_labels[scenario], scenario_labels),
         similarity = factor(c("similar (shift up to 1)", "medium similarity (shift up to 2)", "dissimilar (shift up to 3)")[to],
                             c("dissimilar (shift up to 3)", "medium similarity (shift up to 2)", "similar (shift up to 1)")))
cat("Fig 3:", n_distinct(paste(long$row, long$sim)), "runs;", sum(is.na(runs$rmse_local)), "runs without local fit;",
    "fewest runs per box:", min(count(long, scenario, similarity, ess_bin, sample)$n), "\n")
whisk <- long |> group_by(scenario, similarity, ess_bin, sample) |>
  summarise(lo = boxplot.stats(RMSE)$stats[1], hi = boxplot.stats(RMSE)$stats[5], .groups = "drop")
p3 <- ggplot(long, aes(ess_bin, RMSE, fill = sample)) +
  geom_boxplot(outlier.shape = NA, color = "grey25", linewidth = 0.25, width = 0.8, position = position_dodge(0.85)) +
  facet_grid(scenario ~ similarity) +
  scale_fill_manual(values = samples, labels = sample_labels, name = "Sample") +
  scale_x_discrete(labels = expression("" <= 1, "2", "3", "4", "" >= 5)) +
  coord_cartesian(ylim = c(min(whisk$lo), max(whisk$hi))) +
  labs(x = expression("ESS ratio" ~ ESS[weighted] / n[0] ~ "(rounded)"), y = "RMSE on the target test data") +
  theme_paper
save_fig(p3, file.path(figdir, "Fig3"), 7.5, 6)
cat("Written:", file.path(figdir, c("Fig2", "Fig3")), "(.pdf, .tif)\n")
