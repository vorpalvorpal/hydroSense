#!/usr/bin/env Rscript
# Plot the issue #15 crossover-scale experiment (B.S01, ARA): deterministic
# centre vs the posterior median + 50% band under c = HC5 / HC1 / HC0.1.
# Reads dev/crossover_experiment.qs2 (written by dev/crossover_experiment.R).
#   Rscript dev/crossover_plot.R

suppressMessages({ library(dplyr); library(ggplot2) })
res <- qs2::qs_read("dev/crossover_experiment.qs2")
labs_c <- names(res)                              # HC5, HC1, HC0.1
pal <- c(HC5 = "#377eb8", HC1 = "#ff7f00", "HC0.1" = "#e41a1c")

START <- min(res[[1]]$band$date); END <- max(res[[1]]$band$date)
xlim  <- ggplot2::scale_x_date(limits = c(START, END))
thm   <- ggplot2::theme_minimal(base_size = 10)
ylab  <- "% species affected (ARA)"
det   <- res[[1]]$det                             # ~c-independent reference

## Panel A: the three medians + deterministic, overlaid.
meds <- dplyr::bind_rows(lapply(labs_c, function(nm)
  dplyr::transmute(res[[nm]]$band, date, median, crossover = nm)))
meds$crossover <- factor(meds$crossover, levels = labs_c)
pA <- ggplot() +
  geom_line(data = det, aes(date, mspaf), colour = "grey45",
            linetype = "21", linewidth = 0.5) +
  geom_line(data = meds, aes(date, median, colour = crossover), linewidth = 0.6) +
  scale_colour_manual(values = pal, name = "crossover c") +
  xlim + labs(title = "A. Posterior median by crossover (grey dashed = deterministic centre)",
              x = NULL, y = ylab) + thm + theme(legend.position = "top")

## Panels B-D: one per crossover, band + median + deterministic.
panel <- function(nm, letter) {
  b <- res[[nm]]$band
  ggplot(b, aes(date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = pal[[nm]], alpha = 0.25) +
    geom_line(data = det, aes(date, mspaf), colour = "grey45",
              linetype = "21", linewidth = 0.5) +
    geom_line(aes(y = median), colour = pal[[nm]], linewidth = 0.6) +
    xlim + labs(title = sprintf("%s. c = %s: median + 50%% band vs deterministic",
                                letter, nm), x = NULL, y = ylab) + thm
}
pB <- panel("HC5", "B"); pC <- panel("HC1", "C"); pD <- panel("HC0.1", "D")

g <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(pA, pB, pC, pD, ncol = 1)
} else {
  gridExtra::grid.arrange(pA, pB, pC, pD, ncol = 1)
}
ggsave("dev/crossover_compare.png", g, width = 10, height = 13, dpi = 120)
cat("WROTE dev/crossover_compare.png\n")
