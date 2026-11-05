# Fig 1 of the paper: the three types of distribution shift in the simulation as DAGs.
# Symbols in Times italic, scenario labels in Arial, 8-12 pt, 7.5 in wide, TIFF 600 dpi LZW RGB
# (and a PDF for viewing).
# Usage (from any directory): Rscript paper/make_fig1.R
# Writes paper/output/figures/Fig1 (.pdf, .tif).
suppressPackageStartupMessages(library(ggplot2))
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
figdir <- file.path(here, "output", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

r <- 0.09  # node radius (data units; coord_fixed keeps circles round)
nodes <- data.frame(id = c("S", "H", "Y", "Xs", "Xg"),
                    x = c(0.12, 0.42, 0.82, 0.42, 0.82), y = c(0.88, 0.58, 0.58, 0.12, 0.12),
                    lab = c("italic(A)", "italic(H)", "italic(Y)", "italic(X)[plain(s)]", "italic(X)[plain(g)]"))
common <- data.frame(from = c("H", "H", "Xs", "Xg"), to = c("Y", "Xs", "Y", "Y"))
scen <- c("covariate", "outcome", "covariate + outcome")
edges <- rbind(cbind(common, scenario = rep(scen, each = 4)),
               data.frame(from = "S", to = c("Xs", "Y", "H"), scenario = scen))
edges <- within(edges, {
  x0 <- nodes$x[match(from, nodes$id)]; y0 <- nodes$y[match(from, nodes$id)]
  x1 <- nodes$x[match(to, nodes$id)];   y1 <- nodes$y[match(to, nodes$id)]
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2)
  xs <- x0 + r * (x1 - x0) / d; ys <- y0 + r * (y1 - y0) / d       # start on the source circle
  xe <- x1 - r * (x1 - x0) / d; ye <- y1 - r * (y1 - y0) / d       # end on the target circle
})
circ <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
  a <- seq(0, 2 * pi, length.out = 100)
  data.frame(id = nodes$id[i], x = nodes$x[i] + r * cos(a), y = nodes$y[i] + r * sin(a))
}))
rep_scen <- function(d) do.call(rbind, lapply(scen, function(s) cbind(d, scenario = s)))
lv <- function(d) transform(d, scenario = factor(scenario, scen))

p <- ggplot() +
  geom_segment(data = lv(edges), aes(xs, ys, xend = xe, yend = ye),
               arrow = arrow(length = unit(0.07, "in"), type = "closed"), linewidth = 0.45) +
  geom_path(data = lv(rep_scen(circ)), aes(x, y, group = id), linewidth = 0.4) +
  # math symbols in Times italic (PLOS allows Arial, Times, Symbol); 12 pt keeps the plotmath
  # subscripts (scaled by 0.7) above the 8 pt minimum
  geom_text(data = lv(rep_scen(nodes)), aes(x, y, label = lab), parse = TRUE, family = "Times New Roman", size = 12 / .pt) +
  facet_wrap(~scenario, nrow = 1, strip.position = "bottom") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(strip.text = element_text(size = 10, margin = margin(t = 4)), panel.spacing = unit(0.35, "in"),
        plot.margin = margin(6, 10, 6, 10))

# bg = "white": theme_void has no background, and PLOS TIFFs must have no alpha channel
ggsave(file.path(figdir, "Fig1.pdf"), p, width = 7.5, height = 2.5, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(figdir, "Fig1.tif"), p, width = 7.5, height = 2.5, units = "in", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", bg = "white")
cat("Written:", file.path(figdir, "Fig1"), "(.pdf, .tif)\n")
