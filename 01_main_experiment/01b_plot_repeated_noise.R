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
data_file <- file.path(root, "results", "main_experiment",
                       "repeated_noise_raw_metrics.csv")
if (!file.exists(data_file)) stop("Missing repeated-noise data: ", data_file)

raw <- read.csv(data_file, check.names=FALSE)
required <- c("replication", "eta", "source", "method", "RMSE")
if (!all(required %in% names(raw))) {
  stop("Repeated-noise data do not contain the required columns.")
}

rep_mean <- stats::aggregate(RMSE ~ eta + replication + method, raw,
                             FUN=function(x) mean(x, na.rm=TRUE))
method_order <- c("Tikhonov discrepancy", "RBF", "Matern 0.5",
                  "Matern 1.5", "Matern 2.5")
method_labels <- c("Tikhonov", "RBF", "M0.5", "M1.5", "M2.5")
rep_mean$method <- factor(rep_mean$method, levels=method_order,
                          labels=method_labels)
eta_order <- sort(unique(rep_mean$eta))
eta_labels <- setNames(sprintf("eta == 10^%d", round(log10(eta_order))),
                       format(eta_order, scientific=TRUE))
rep_mean$noise_level <- factor(
  eta_labels[format(rep_mean$eta, scientific=TRUE)],
  levels=unname(eta_labels))

colours <- c("Tikhonov"="#6B7280", "RBF"="#440154",
             "M0.5"="#31688E", "M1.5"="#35B779", "M2.5"="#FDE725")

p <- ggplot2::ggplot(rep_mean,
                     ggplot2::aes(method, RMSE, fill=method)) +
  ggplot2::geom_boxplot(width=.82, linewidth=.58,
                        outlier.shape=21, outlier.size=1.45,
                        outlier.stroke=.45, outlier.colour="#30343B",
                        outlier.fill="white", colour="#30343B", alpha=.92) +
  ggplot2::facet_wrap(~noise_level, nrow=1, scales="free_y",
                      labeller=ggplot2::label_parsed) +
  ggplot2::scale_fill_manual(values=colours, guide="none") +
  ggplot2::labs(x=NULL, y="Six-source mean RMSE",
                title="RMSE across 20 paired-noise repetitions") +
  ggplot2::theme_minimal(base_size=11) +
  ggplot2::theme(panel.grid.major.x=ggplot2::element_blank(),
                 panel.grid.minor=ggplot2::element_blank(),
                 panel.background=ggplot2::element_rect(fill="white",
                                                        colour=NA),
                 plot.background=ggplot2::element_rect(fill="white",
                                                       colour=NA),
                 panel.spacing=grid::unit(1.2, "lines"),
                 strip.text=ggplot2::element_text(face="bold"),
                 plot.title=ggplot2::element_text(hjust=.5, face="bold"),
                 axis.text.x=ggplot2::element_text(face="bold"))

result_figure <- file.path(root, "results", "main_experiment", "figures",
                           "repeated_noise_rmse_boxplot.png")
chapter_figure <- result_figure
dir.create(dirname(result_figure), recursive=TRUE, showWarnings=FALSE)
ggplot2::ggsave(result_figure, p, width=11.8, height=4.35, dpi=240,
                bg="white")
ggplot2::ggsave(chapter_figure, p, width=11.8, height=4.35, dpi=240,
                bg="white")
message("Repeated-noise box plot written to ", chapter_figure)
