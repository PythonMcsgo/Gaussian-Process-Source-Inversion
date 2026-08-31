
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

script_location <- function() {
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  if (length(frames)) {
    return(dirname(normalizePath(tail(frames, 1)[[1]], winslash = "/",
                                  mustWork = FALSE)))
  }
  script_args <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
  if (length(hit)) {
    return(dirname(normalizePath(hit[1], winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
script_dir <- script_location()
root <- normalizePath(get_arg("root", file.path(script_dir, "..")),
                      winslash = "/", mustWork = TRUE)
result_prefix <- get_arg("prefix", "active_single_")
report_dir <- normalizePath(
  get_arg("output", file.path(root, "results", "active_learning_6sources_report")),
  winslash = "/", mustWork = FALSE
)
dirs <- list(
  root = report_dir,
  figures = file.path(report_dir, "figures"),
  tables = file.path(report_dir, "tables"),
  data = file.path(report_dir, "data")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

source_order <- c("spike", "mixture", "dipole", "ring", "edge", "sine")
source_labels <- c(
  spike = "Narrow spike", mixture = "Broad + sharp", dipole = "Dipole",
  ring = "Ring", edge = "Moving edge", sine = "Sine-cosine"
)
method_order <- c(
  "uniform", "random", "max_source_variance", "variance_times_gradient"
)
method_labels <- c(
  uniform = "Uniform", random = "Random",
  max_source_variance = "Maximum variance",
  variance_times_gradient = "Variance x gradient"
)
method_colours <- c(
  uniform = "#440154", random = "#414487",
  max_source_variance = "#22A884",
  variance_times_gradient = "#FDE725"
)

read_one <- function(source, filename) {
  path <- file.path(root, "results", paste0(result_prefix, source),
                    "results", filename)
  if (!file.exists(path)) stop("Missing active-learning result: ", path)
  z <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  z$source <- source
  z$source_label <- unname(source_labels[source])
  z
}

summary_all <- do.call(rbind, lapply(source_order, read_one,
                                    filename = "single_dataset_metrics.csv"))
gain_all <- do.call(rbind, lapply(source_order, read_one,
                                 filename = "uniform_comparison.csv"))
calibration_all <- do.call(rbind, lapply(
  source_order, read_one, filename = "calibration.csv"
))

summary_all <- summary_all[
  order(match(summary_all$source, source_order),
        match(summary_all$method, method_order)),
]
summary_all$source_label <- factor(summary_all$source_label,
                                   levels = unname(source_labels[source_order]))
summary_all$method_label <- factor(
  unname(method_labels[summary_all$method]),
  levels = unname(method_labels[method_order])
)
gain_all$source_label <- factor(gain_all$source_label,
                                levels = unname(source_labels[source_order]))
gain_all$method_label <- factor(
  unname(method_labels[gain_all$method]),
  levels = unname(method_labels[method_order[-1]])
)
calibration_all$source_label <- factor(
  calibration_all$source_label, levels = unname(source_labels[source_order])
)
calibration_all$method_label <- factor(
  unname(method_labels[calibration_all$method]),
  levels = unname(method_labels[method_order])
)

write.csv(summary_all, file.path(dirs$data, "six_source_summary.csv"),
          row.names = FALSE)
write.csv(gain_all, file.path(dirs$data, "six_source_relative_gain.csv"),
          row.names = FALSE)
write.csv(calibration_all,
          file.path(dirs$data, "six_source_calibration.csv"), row.names = FALSE)

theme_report <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

p_rmse <- ggplot(summary_all,
                 aes(method_label, RMSE, fill = method)) +
  geom_col(width = 0.70, colour = "grey20", linewidth = 0.20) +
  facet_wrap(~source_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = method_colours, guide = "none") +
  labs(
    title = "Source recovery at a total budget of 16 observations",
    subtitle = "One fixed noisy dataset per source; all designs share that dataset",
    x = NULL, y = "Source RMSE"
  ) +
  theme_report +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 8))
ggsave(file.path(dirs$figures, "rmse_six_sources.png"), p_rmse,
       width = 11.5, height = 7.4, dpi = 220, bg = "white")

p_gain <- ggplot(gain_all,
                 aes(method_label, relative_RMSE_gain_percent, colour = method)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
  geom_point(size = 2.2) +
  facet_wrap(~source_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = method_colours, guide = "none") +
  labs(
    title = "RMSE gain relative to the uniform design",
    subtitle = "Positive values favour the alternative on the fixed dataset",
    x = NULL, y = "Relative RMSE gain (%)"
  ) +
  theme_report +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
ggsave(file.path(dirs$figures, "relative_gain_six_sources.png"), p_gain,
       width = 11.5, height = 7.4, dpi = 220, bg = "white")

p_calibration <- ggplot(
  calibration_all,
  aes(nominal, empirical, colour = method, group = method)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.75) + geom_point(size = 1.3) +
  facet_wrap(~source_label, ncol = 3) +
  scale_colour_manual(values = method_colours,
                      labels = method_labels, name = "Design") +
  coord_equal(xlim = c(0.5, 1), ylim = c(0.45, 1)) +
  labs(title = "Spatial credible-interval calibration",
       x = "Nominal credibility", y = "Empirical coverage") +
  theme_report
ggsave(file.path(dirs$figures, "calibration_six_sources.png"),
       p_calibration, width = 11.5, height = 7.8, dpi = 220, bg = "white")

# Source catalogue and representative sensor designs. The truth is scaled by
# max(abs(f)) only for this location plot so signed and non-negative sources
# share one colour scale; all numerical scores use the unscaled source.
design_fields <- list(); design_points <- list(); catalogue <- list()
field_id <- point_id <- catalogue_id <- 1
for (source in source_order) {
  path <- file.path(root, "results", paste0(result_prefix, source), "results",
                    "complete_experiment.rds")
  if (!file.exists(path)) stop("Missing complete experiment: ", path)
  fit <- readRDS(path)
  scale <- max(abs(fit$f_truth))
  base <- data.frame(
    fit$truth_grid, value = fit$f_truth / scale,
    source = source, source_label = unname(source_labels[source])
  )
  catalogue[[catalogue_id]] <- base
  catalogue_id <- catalogue_id + 1
  for (method in method_order) {
    key <- paste(method, 16, sep = "::")
    selected <- fit$representative[[key]]$indices
    one <- base
    one$method <- method
    one$method_label <- unname(method_labels[method])
    design_fields[[field_id]] <- one
    field_id <- field_id + 1
    pts <- fit$candidate_grid[selected, c("x", "y"), drop = FALSE]
    pts$source <- source
    pts$source_label <- unname(source_labels[source])
    pts$method <- method
    pts$method_label <- unname(method_labels[method])
    design_points[[point_id]] <- pts
    point_id <- point_id + 1
  }
}
catalogue <- do.call(rbind, catalogue)
design_fields <- do.call(rbind, design_fields)
design_points <- do.call(rbind, design_points)
for (object in c("catalogue", "design_fields", "design_points")) {
  z <- get(object)
  z$source_label <- factor(z$source_label,
                           levels = unname(source_labels[source_order]))
  if ("method_label" %in% names(z)) {
    z$method_label <- factor(z$method_label,
                             levels = unname(method_labels[method_order]))
  }
  assign(object, z)
}

source_panel <- function(fields, title, sensors = NULL) {
  signed <- min(fields$value) < -sqrt(.Machine$double.eps)
  fill_scale <- if (signed) {
    scale_fill_gradient2(
      low = "#440154", mid = "white", high = "#FDE725", midpoint = 0,
      limits = c(-1, 1), oob = scales::squish, name = "Scaled f"
    )
  } else {
    scale_fill_viridis_c(
      option = "D", limits = c(0, 1), oob = scales::squish,
      name = "Scaled f"
    )
  }
  ggplot(fields, aes(x, y, fill = value)) +
    geom_raster() +
    {if (is.null(sensors)) NULL else geom_point(
      data = sensors, inherit.aes = FALSE, aes(x, y),
      shape = 21, fill = "white", colour = "black", stroke = 0.25, size = 1
    )} +
    fill_scale + coord_equal(expand = FALSE) +
    labs(title = title, x = "x", y = "y") + theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
          legend.key.height = grid::unit(0.48, "cm"),
          legend.key.width = grid::unit(0.18, "cm"))
}

source_panels <- lapply(source_order, function(source) {
  source_panel(catalogue[catalogue$source == source, ],
               unname(source_labels[source]))
})
p_sources <- arrangeGrob(
  grobs = source_panels, ncol = 3,
  top = grid::textGrob("Six source families",
                       gp = grid::gpar(fontface = "bold"))
)
ggsave(file.path(dirs$figures, "source_catalogue.png"), p_sources,
       width = 10.4, height = 7.0, dpi = 220, bg = "white")

save_design_figure <- function(sources, filename) {
  panels <- list()
  for (source in sources) {
    for (method in method_order) {
      fields <- design_fields[
        design_fields$source == source & design_fields$method == method, ]
      points <- design_points[
        design_points$source == source & design_points$method == method, ]
      panels[[length(panels) + 1]] <- source_panel(
        fields,
        paste(unname(source_labels[source]), unname(method_labels[method]),
              sep = " - "),
        points
      )
    }
  }
  p <- arrangeGrob(grobs = panels, ncol = 4)
  ggsave(file.path(dirs$figures, filename), p,
         width = 12.2, height = 8.4, dpi = 220, bg = "white")
}
save_design_figure(source_order[1:3], "designs_spike_mixture_dipole.png")
save_design_figure(source_order[4:6], "designs_ring_edge_sine.png")

format_cell <- function(value) sprintf("%.4f", value)
rmse_wide <- data.frame(Source = unname(source_labels[source_order]),
                        check.names = FALSE)
for (method in method_order) {
  z <- summary_all[summary_all$method == method, ]
  z <- z[match(source_order, z$source), ]
  rmse_wide[[unname(method_labels[method])]] <-
    format_cell(z$RMSE)
}
write.csv(rmse_wide, file.path(dirs$data, "rmse_table.csv"), row.names = FALSE)

gain_wide <- data.frame(Source = unname(source_labels[source_order]),
                        check.names = FALSE)
for (method in method_order[-1]) {
  z <- gain_all[gain_all$method == method, ]
  z <- z[match(source_order, z$source), ]
  gain_wide[[unname(method_labels[method])]] <-
    sprintf("%.1f", z$relative_RMSE_gain_percent)
}
write.csv(gain_wide, file.path(dirs$data, "gain_table.csv"), row.names = FALSE)

winner <- do.call(rbind, lapply(source_order, function(source) {
  z <- summary_all[summary_all$source == source, ]
  data.frame(
    Source = unname(source_labels[source]),
    `Lowest RMSE` = as.character(z$method_label[which.min(z$RMSE)]),
    `Lowest calibration MAE` =
      as.character(z$method_label[which.min(z$CalibrationMAE)]),
    `Lowest CRPS` = as.character(z$method_label[which.min(z$CRPS)]),
    check.names = FALSE
  )
}))
write.csv(winner, file.path(dirs$data, "metric_winners.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(dirs$data, "session_info.txt"))
message("Active-learning report assets written to: ", report_dir)
