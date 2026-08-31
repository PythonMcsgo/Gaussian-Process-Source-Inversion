
if (!requireNamespace("ggplot2", quietly=TRUE)) {
  stop("Install ggplot2 first: install.packages('ggplot2')")
}

script_location <- function() {
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  if (length(frames)) {
    return(dirname(normalizePath(tail(frames, 1)[[1]], winslash="/")))
  }
  z <- commandArgs(trailingOnly=FALSE)
  hit <- sub("^--file=", "", z[grepl("^--file=", z)])
  if (length(hit)) dirname(normalizePath(hit[[1]], winslash="/")) else getwd()
}

root <- normalizePath(file.path(script_location(), ".."), winslash="/",
                      mustWork=TRUE)
result_dir <- file.path(root, "results", "main_experiment")
result_fig <- file.path(result_dir, "figures")
chapter_fig <- result_fig
dir.create(result_fig, recursive=TRUE, showWarnings=FALSE)
dir.create(chapter_fig, recursive=TRUE, showWarnings=FALSE)

metrics <- read.csv(file.path(result_dir, "gp_metrics.csv"),
                    check.names=FALSE)
calibration <- read.csv(file.path(result_dir, "gp_calibration.csv"),
                        check.names=FALSE)
metrics <- metrics[metrics$hyper_method == "multistart_eb", ]
calibration <- calibration[calibration$hyper_method == "multistart_eb", ]

kernel_order <- c("rbf", "matern05", "matern15", "matern25")
kernel_labels <- c(rbf="RBF", matern05="Matern 0.5",
                   matern15="Matern 1.5", matern25="Matern 2.5")
source_order <- c("spike", "broad_plus_sharp", "dipole", "ring",
                  "moving_edge", "sine")
source_labels <- c(spike="Narrow spike",
                   broad_plus_sharp="Broad + sharp", dipole="Dipole",
                   ring="Ring", moving_edge="Moving edge",
                   sine="Sine-cosine")
colours <- c("RBF"="#440154", "Matern 0.5"="#3B528B",
             "Matern 1.5"="#21918C", "Matern 2.5"="#5DC863")

metrics$kernel_label <- factor(kernel_labels[metrics$kernel],
                               levels=unname(kernel_labels[kernel_order]))
score_names <- c("RMSE", "CalibrationMAE", "CRPS", "NLPD")
score_labels <- c(RMSE="RMSE", CalibrationMAE="Calibration MAE",
                  CRPS="CRPS", NLPD="NLPD")
score_rows <- do.call(rbind, lapply(score_names, function(metric) {
  z <- stats::aggregate(metrics[[metric]],
                        by=list(kernel=metrics$kernel_label), mean)
  names(z)[2] <- "value"
  z$metric <- score_labels[[metric]]
  z
}))
score_rows$metric <- factor(score_rows$metric,
                            levels=unname(score_labels[score_names]))

p_scores <- ggplot2::ggplot(score_rows,
                            ggplot2::aes(kernel, value, fill=kernel)) +
  ggplot2::geom_col(width=.68) +
  ggplot2::facet_wrap(~metric, nrow=1, scales="free_y") +
  ggplot2::scale_fill_manual(values=colours, guide="none") +
  ggplot2::labs(x=NULL, y="Six-source mean",
                title="Fixed-design reconstruction and probabilistic scores",
                subtitle="Lower values are better for every displayed score") +
  ggplot2::theme_minimal(base_size=11) +
  ggplot2::theme(panel.grid.major.x=ggplot2::element_blank(),
                 panel.grid.minor=ggplot2::element_blank(),
                 panel.background=ggplot2::element_rect(fill="white", colour=NA),
                 plot.background=ggplot2::element_rect(fill="white", colour=NA),
                 strip.text=ggplot2::element_text(face="bold"),
                 axis.text.x=ggplot2::element_text(angle=25, hjust=1))

calibration$kernel_label <- factor(kernel_labels[calibration$kernel],
                                   levels=unname(kernel_labels[kernel_order]))
calibration$source_label <- factor(source_labels[calibration$source],
                                   levels=unname(source_labels[source_order]))

make_calibration_plot <- function(source_key) {
  z <- calibration[calibration$source == source_key, , drop=FALSE]
  ggplot2::ggplot(
    z,
    ggplot2::aes(nominal, empirical, colour=kernel_label,
                 group=kernel_label)) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype=2,
                         colour="#555B63", linewidth=.65) +
    ggplot2::geom_line(linewidth=.95) +
    ggplot2::geom_point(size=2.05) +
    ggplot2::coord_equal(xlim=c(.48, 1.01), ylim=c(.48, 1.01),
                         expand=FALSE) +
    ggplot2::scale_colour_manual(values=colours, name="Kernel") +
    ggplot2::labs(x="Nominal credibility",
                  y="Empirical spatial coverage",
                  title=paste(source_labels[[source_key]],
                              "source: fixed-design calibration")) +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::theme(
      panel.grid.minor=ggplot2::element_blank(),
      panel.background=ggplot2::element_rect(fill="white", colour=NA),
      plot.background=ggplot2::element_rect(fill="white", colour=NA),
      legend.position="bottom",
      legend.box.margin=ggplot2::margin(t=-2),
      plot.title=ggplot2::element_text(hjust=.5, face="bold")
    )
}

save_both <- function(name, plot, width, height) {
  for (directory in c(result_fig, chapter_fig)) {
    ggplot2::ggsave(file.path(directory, name), plot, width=width,
                    height=height, dpi=240, bg="white")
  }
}
save_both("fixed_design_score_overview.png", p_scores, 11.6, 4.6)
for (source_key in source_order) {
  save_both(paste0("fixed_design_calibration_", source_key, ".png"),
            make_calibration_plot(source_key), 7.4, 5.9)
}
message("Fixed-design scores and six separate calibration figures refreshed.")
