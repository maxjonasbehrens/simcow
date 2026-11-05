# Fig 4 of the paper (clinical application): per-patient mean absolute error in Gy and weighted
# dose CDFs for the five SENIOR target centers. The figure shows patient-level values, so it is
# written to <outdir> (outside this repository), next to the output of run_senior.R.
# Usage: Rscript paper/make_fig4.R <data.csv> <outdir>   (run paper/run_senior.R first)
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(gridExtra) })
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "R", "read_clinical.R"))
source(file.path(here, "..", "R", "create_cow.R"))
args <- commandArgs(TRUE)
if (length(args) != 2) stop("Usage: Rscript paper/make_fig4.R <data.csv> <outdir>")
outdir <- outdir_outside_repo(args[2], file.path(here, ".."))
pdf(NULL)   # ggplotGrob() opens a device; keep it from writing Rplots.pdf to the working directory

samples <- c(local = "#1baf7a", global = "#eb6834", weighted = "#2a78d6")
theme_paper <- theme_bw(base_size = 10, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom")

senior <- read_clinical(args[1])
covs <- c("BeginRadioTherapie", "CCI", "Hb", "Leukozyten")
centers <- c(2, 3, 6, 7, 11)
tmp <- as.data.frame(senior |> filter(Zentrum %in% centers))
metrics <- read.csv(file.path(outdir, "senior_metrics.csv"))
center_lab <- setNames(paste0("Center ", metrics$center, " (n = ", metrics$n0, ")")[metrics$model == "weighted"], centers)

obs <- do.call(rbind, lapply(centers, function(t) read.csv(file.path(outdir, paste0("senior_obs_errors_center", t, ".csv"))))) |>
  arrange(match(center, centers), model != "weighted", model, ID) |>   # row order of the paper's run, so the jitter matches
  mutate(sample = factor(model, names(samples)), center = factor(center_lab[as.character(center)], center_lab))
set.seed(1)
obs <- obs |> mutate(x = as.numeric(sample) + runif(n(), -0.25, 0.25))
meds <- obs |> group_by(center, sample) |> summarise(m = median(mae_gy), .groups = "drop")
n_lab <- metrics |>
  mutate(sample = factor(model, names(samples)), center = factor(center_lab[as.character(center)], center_lab),
         label = case_when(model == "local" ~ paste0("n = ", n0), model == "global" ~ paste0("n = ", nrow(tmp)),
                           TRUE ~ paste0("ESS = ", round(ESS_median))))
p4a <- ggplot(obs, aes(x, mae_gy)) +
  geom_line(aes(group = ID), color = "grey70", alpha = 0.4, linewidth = 0.2) +
  geom_point(aes(color = sample), size = 1.1) +
  geom_segment(data = meds, aes(x = as.numeric(sample) - 0.3, xend = as.numeric(sample) + 0.3, y = m, yend = m),
               color = "grey10", linewidth = 0.6) +
  geom_text(data = n_lab, aes(as.numeric(sample), 0.12, label = label), size = 8 / .pt, family = "Arial", color = "grey30") +
  facet_grid(center ~ ., switch = "y") +
  scale_color_manual(values = samples, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = names(samples)) +
  scale_y_log10(breaks = c(0.3, 1, 3, 10, 30), limits = c(0.09, 70)) +
  labs(x = NULL, y = "Mean absolute error per patient (Gy, log scale)") +
  theme_paper + theme(strip.placement = "outside")

# Weighted dose CDFs with the weights estimated on all patients (descriptive)
cdf <- do.call(rbind, lapply(centers, function(t) {
  w <- create_cow(tmp, "Zentrum", t, covs, outcome = "DosisEQDW", seed = 3e6 + 1000 * t)$weight
  dose <- sqrt(tmp$DosisEQDW)
  step <- function(x, w, s) { o <- order(x); data.frame(dose = x[o], cdf = cumsum(w[o]) / sum(w), sample = s) }
  rbind(step(dose[tmp$Zentrum == t], rep(1, sum(tmp$Zentrum == t)), "local"),
        step(dose, rep(1, nrow(tmp)), "global"),
        step(dose, w, "weighted")) |> mutate(center = t)
})) |> mutate(sample = factor(sample, names(samples)), center = factor(center_lab[as.character(center)], center_lab))
p4b <- ggplot(cdf, aes(dose, cdf, color = sample, linetype = sample)) +
  geom_step(linewidth = 0.6) +
  facet_grid(center ~ .) +
  scale_color_manual(values = samples, name = NULL) +
  scale_linetype_manual(values = c(local = "dashed", global = "dotted", weighted = "solid"), name = NULL) +
  scale_y_continuous(position = "right", breaks = c(0, 0.5, 1)) +
  labs(x = "Radiotherapy dose (Gy)", y = "(weighted) CDF") +
  theme_paper + theme(strip.text = element_blank())
g4a <- ggplotGrob(p4a)
g4a$widths[g4a$layout$l[grepl("^panel", g4a$layout$name)][1]] <- unit(2, "null")   # error panels twice as wide as CDFs
g4 <- cbind(g4a, ggplotGrob(p4b), size = "max")
ggsave(file.path(outdir, "Fig4.pdf"), g4, width = 7.5, height = 8.5, units = "in", device = cairo_pdf)
ggsave(file.path(outdir, "Fig4.tif"), g4, width = 7.5, height = 8.5, units = "in", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw")
cat("Written:", file.path(outdir, "Fig4"), "(.pdf, .tif)\n")
