# Assemble the six independent ARWM runs into report-ready tables and figures.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Run install.packages('ggplot2').")
  }
})

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", full_args[grepl("^--file=", full_args)])
script_file <- if (length(file_arg)) {
  normalizePath(file_arg[1], winslash = "/", mustWork = FALSE)
} else tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                error = function(e) NA_character_)
script_dir <- if (is.na(script_file) || !length(script_file)) getwd() else
  dirname(script_file)
root <- normalizePath(get_arg("root", file.path(dirname(script_dir), "results")),
                      winslash = "/", mustWork = FALSE)
result_prefix <- get_arg("prefix", "MCMC_results_")
report_dir <- normalizePath(get_arg("output", file.path(root, "MCMC_report")),
                            winslash = "/", mustWork = FALSE)
data_dir <- file.path(report_dir, "data")
figure_dir <- file.path(report_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sources <- c("spike", "mixture", "dipole", "ring", "edge", "sine")
labels <- c(
  spike = "Narrow spike", mixture = "Broad + sharp", dipole = "Dipole",
  ring = "Ring", edge = "Moving edge", sine = "Sine-cosine"
)

read_one <- function(source, filename) {
  read.csv(file.path(root, paste0(result_prefix, source), filename),
           check.names = FALSE)
}
bind_files <- function(filename) {
  do.call(rbind, lapply(sources, read_one, filename = filename))
}

metrics <- bind_files("mcmc_metrics.csv")
summary <- bind_files("mcmc_diagnostic_summary.csv")
detail <- bind_files("mcmc_diagnostics.csv")
chains <- bind_files("mcmc_chain_diagnostics.csv")
stability <- bind_files("mcmc_stability.csv")
calibration <- bind_files("mcmc_calibration.csv")
adaptation <- bind_files("mcmc_adaptation_history.csv")

for (object_name in c("metrics", "summary", "detail", "chains", "stability",
                      "calibration", "adaptation")) {
  object <- get(object_name)
  object$source_label <- unname(labels[object$source])
  assign(object_name, object)
  write.csv(object, file.path(data_dir, paste0(object_name, ".csv")),
            row.names = FALSE)
}

diag_long <- rbind(
  data.frame(source = summary$source_label, diagnostic = "Maximum Rhat",
             value = pmax(summary$max_Rhat_hyperparameters,
                          summary$max_Rhat_downstream)),
  data.frame(source = summary$source_label, diagnostic = "Minimum bulk ESS",
             value = pmin(summary$min_ESS_hyperparameters,
                          summary$min_ESS_downstream))
)
p_diag <- ggplot(diag_long, aes(reorder(source, value), value)) +
  geom_col(fill = "#2C7FB8", width = .72) +
  geom_hline(data = data.frame(diagnostic = "Maximum Rhat", threshold = 1.01),
             aes(yintercept = threshold), linetype = 2, colour = "#D95F0E") +
  facet_wrap(~diagnostic, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(x = NULL, y = NULL,
       title = "Worst-case four-chain convergence diagnostics",
       subtitle = "Rhat includes hyperparameters and four downstream posterior means") +
  theme_minimal(base_size = 11)
ggsave(file.path(figure_dir, "diagnostic_summary.png"), p_diag,
       width = 10, height = 5.3, dpi = 200)

p_accept <- ggplot(chains,
                   aes(factor(chain), sampling_acceptance_rate,
                       group = source_label, colour = source_label)) +
  geom_hline(yintercept = .234, linetype = 2, colour = "grey35") +
  geom_line(linewidth = .65) + geom_point(size = 1.8) +
  facet_wrap(~source_label, ncol = 3) +
  scale_colour_viridis_d(option = "D", guide = "none") +
  labs(x = "Chain", y = "Sampling-phase acceptance rate",
       title = "Frozen-proposal acceptance rates",
       subtitle = "The dashed 0.234 line is an adaptation reference, not a pass/fail target") +
  theme_minimal(base_size = 10)
ggsave(file.path(figure_dir, "acceptance_rates.png"), p_accept,
       width = 10.5, height = 6.5, dpi = 200)

p_stability <- ggplot(stability,
                      aes(retained_draws, mean_abs_change_from_Smax,
                          colour = source_label)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  facet_wrap(~source_label, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(120, 240, 480)) +
  scale_colour_viridis_d(option = "D", guide = "none") +
  labs(x = "Balanced mixture components S",
       y = "Mean absolute change from S=480",
       title = "Monte Carlo stability of the posterior mean") +
  theme_minimal(base_size = 10)
ggsave(file.path(figure_dir, "posterior_mean_stability.png"), p_stability,
       width = 10.5, height = 6.5, dpi = 200)

p_variance <- ggplot(stability,
                     aes(retained_draws, mean_between_variance /
                           mean_total_variance, colour = source_label)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  facet_wrap(~source_label, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(120, 240, 480)) +
  scale_colour_viridis_d(option = "D", guide = "none") +
  labs(x = "Balanced mixture components S",
       y = "Between / total variance",
       title = "Contribution of hyperparameter uncertainty") +
  theme_minimal(base_size = 10)
ggsave(file.path(figure_dir, "variance_decomposition.png"), p_variance,
       width = 10.5, height = 6.5, dpi = 200)

p_cal <- ggplot(calibration, aes(nominal, empirical)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey35") +
  geom_line(colour = "#2C7FB8", linewidth = .7) +
  geom_point(colour = "#D95F0E", size = 1.5) +
  facet_wrap(~source_label, ncol = 3) +
  coord_equal(xlim = c(.48, 1), ylim = c(.48, 1)) +
  labs(x = "Nominal credibility", y = "Empirical coverage",
       title = "Calibration of the four-chain hyperparameter mixture") +
  theme_minimal(base_size = 10)
ggsave(file.path(figure_dir, "calibration_six_sources.png"), p_cal,
       width = 10.5, height = 7, dpi = 200)

message("Report data and figures written to ", report_dir)
