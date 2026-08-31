#!/usr/bin/env Rscript

# Repeated-noise robustness study for the six-source 2D Poisson experiment.
# The same standard-normal noise direction is used across eta within a
# replication.  Consequently, method and noise-level contrasts are paired.

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

profile <- "paper"
MASTER_SEED <- 20161387
repetitions <- 20
workers <- as.integer(get_arg("workers", 4))
force <- as_flag(get_arg("force", "false"))
noise_levels <- c(1e-3, 1e-2, 1e-1)
if (repetitions < 2) stop("At least two repetitions are required.")
if (any(!is.finite(noise_levels)) || any(noise_levels <= 0)) {
  stop("Noise levels must be positive finite numbers.")
}

script_location <- function() {
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  if (length(frames)) {
    return(dirname(normalizePath(tail(frames, 1)[[1]], winslash = "/",
                                  mustWork = FALSE)))
  }
  args_all <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
  if (length(hit)) {
    return(dirname(normalizePath(hit[1], winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
script_dir <- script_location()
root <- dirname(script_dir)
fixed_script <- file.path(script_dir, "fixed_design_and_hyperparameters.R")
out_root <- normalizePath(get_arg(
  "output", file.path(root, "results", "repeated_noise")
), winslash = "/", mustWork = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(out_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
jobs <- expand.grid(replication = seq_len(repetitions), eta = noise_levels,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
jobs$seed <- MASTER_SEED + jobs$replication * 100000L
jobs$tag <- sprintf("eta_%1.0e_rep_%02d", jobs$eta, jobs$replication)
jobs$directory <- file.path(out_root, "jobs", jobs$tag)

run_one <- function(j) {
  row <- jobs[j, ]
  metrics_path <- file.path(row$directory, "results", "metrics.csv")
  tikh_path <- file.path(row$directory, "results", "tikhonov_metrics.csv")
  cal_path <- file.path(row$directory, "results", "calibration_curves.csv")
  if (!force && all(file.exists(c(metrics_path, tikh_path, cal_path)))) {
    return(data.frame(job = row$tag, status = 0L, elapsed = 0))
  }
  dir.create(row$directory, recursive = TRUE, showWarnings = FALSE)
  command <- c(
    shQuote(fixed_script), paste0("--profile=", profile),
    paste0("--seed=", row$seed), paste0("--noise-ratio=", row$eta),
    "--kernels=rbf,matern05,matern15,matern25",
    "--hyper-methods=multistart_eb", "--skip-ensemble=true",
    "--no-plots=true", "--lean-output=true", "--force=true",
    paste0("--output=", shQuote(row$directory))
  )
  log_file <- file.path(row$directory, "run.log")
  start <- proc.time()[3]
  status <- system2(rscript, command, stdout = log_file, stderr = log_file)
  data.frame(job = row$tag, status = as.integer(status),
             elapsed = unname(proc.time()[3] - start))
}

workers <- max(1, min(workers, nrow(jobs)))
message(sprintf("Running %d paired-noise jobs with %d worker(s).", nrow(jobs), workers))
if (workers > 1) {
  cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, c("jobs", "force", "fixed_script", "profile",
                                     "rscript", "run_one"), envir = environment())
  status <- do.call(rbind, parallel::parLapply(cluster, seq_len(nrow(jobs)), run_one))
  parallel::stopCluster(cluster)
  on.exit(NULL, add = FALSE)
} else {
  status <- do.call(rbind, lapply(seq_len(nrow(jobs)), run_one))
}
write.csv(status, file.path(out_root, "job_status.csv"), row.names = FALSE)
if (any(status$status != 0L)) {
  stop("Repeated-noise jobs failed: ", paste(status$job[status$status != 0L], collapse = ", "))
}

read_job <- function(j) {
  row <- jobs[j, ]
  result_dir <- file.path(row$directory, "results")
  gp <- read.csv(file.path(result_dir, "metrics.csv"), check.names = FALSE)
  gp <- gp[gp$hyper_method == "multistart_eb", ]
  gp$method <- gp$label
  gp$replication <- row$replication; gp$eta <- row$eta; gp$seed <- row$seed
  gp <- gp[c("replication", "eta", "seed", "source", "source_label", "method",
             "RMSE", "MAE", "Coverage95", "Width95", "CalibrationMAE", "CRPS", "NLPD")]
  tik <- read.csv(file.path(result_dir, "tikhonov_metrics.csv"), check.names = FALSE)
  tik$method <- "Tikhonov discrepancy"
  tik$replication <- row$replication; tik$eta <- row$eta; tik$seed <- row$seed
  tik$Coverage95 <- tik$Width95 <- tik$CalibrationMAE <- tik$CRPS <- tik$NLPD <- NA_real_
  tik <- tik[names(gp)]
  cal <- read.csv(file.path(result_dir, "calibration_curves.csv"), check.names = FALSE)
  cal <- cal[cal$hyper_method == "multistart_eb", ]
  cal$method <- cal$label; cal$replication <- row$replication
  cal$eta <- row$eta; cal$seed <- row$seed
  list(metrics = rbind(gp, tik), calibration = cal)
}
pieces <- lapply(seq_len(nrow(jobs)), read_job)
raw <- do.call(rbind, lapply(pieces, `[[`, "metrics"))
calibration <- do.call(rbind, lapply(pieces, `[[`, "calibration"))
write.csv(raw, file.path(out_root, "raw_metrics.csv"), row.names = FALSE)
write.csv(calibration, file.path(out_root, "raw_calibration.csv"), row.names = FALSE)
main_result_dir <- file.path(root, "results", "main_experiment")
dir.create(main_result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(raw, file.path(main_result_dir, "repeated_noise_raw_metrics.csv"),
          row.names = FALSE)

metrics <- c("RMSE", "MAE", "Coverage95", "CalibrationMAE", "CRPS", "NLPD")
rep_mean <- aggregate(raw[metrics],
  by = raw[c("eta", "replication", "method")], FUN = function(x) mean(x, na.rm = TRUE))
for (metric in metrics) rep_mean[[metric]][!is.finite(rep_mean[[metric]])] <- NA_real_
summary_rows <- lapply(split(rep_mean, interaction(rep_mean$eta, rep_mean$method, drop = TRUE)), function(z) {
  values <- lapply(metrics, function(metric) {
    x <- z[[metric]][is.finite(z[[metric]])]
    if (!length(x)) return(c(mean = NA, sd = NA, mcse = NA, lower = NA, upper = NA))
    se <- sd(x) / sqrt(length(x)); critical <- qt(.975, df = max(length(x) - 1, 1))
    c(mean = mean(x), sd = sd(x), mcse = se,
      lower = mean(x) - critical * se, upper = mean(x) + critical * se)
  })
  out <- data.frame(eta = z$eta[1], method = z$method[1], repetitions = nrow(z))
  for (i in seq_along(metrics)) for (s in names(values[[i]])) {
    out[[paste(metrics[i], s, sep = "_")]] <- values[[i]][s]
  }
  out
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(summary$eta, summary$RMSE_mean), ]
write.csv(summary, file.path(out_root, "summary_with_mcse.csv"), row.names = FALSE)

# Within-replication, within-source paired RMSE differences relative to
# Matern 0.5.  The inferential unit is a noise replication after averaging
# the six source-specific contrasts.
gp_raw <- raw[raw$method != "Tikhonov discrepancy", ]
reference <- gp_raw[gp_raw$method == "Matern 0.5", c("eta", "replication", "source", "RMSE")]
names(reference)[4] <- "RMSE_reference"
paired <- merge(gp_raw, reference, by = c("eta", "replication", "source"))
paired$delta_RMSE <- paired$RMSE - paired$RMSE_reference
paired_rep <- aggregate(delta_RMSE ~ eta + replication + method, paired, mean)
paired_summary <- do.call(rbind, lapply(split(paired_rep,
  interaction(paired_rep$eta, paired_rep$method, drop = TRUE)), function(z) {
  se <- sd(z$delta_RMSE) / sqrt(nrow(z)); critical <- qt(.975, nrow(z) - 1)
  data.frame(eta = z$eta[1], method = z$method[1], repetitions = nrow(z),
             mean_delta_RMSE = mean(z$delta_RMSE), mcse = se,
             lower95 = mean(z$delta_RMSE) - critical * se,
             upper95 = mean(z$delta_RMSE) + critical * se)
}))
paired_summary <- paired_summary[order(paired_summary$eta, paired_summary$mean_delta_RMSE), ]
write.csv(paired, file.path(out_root, "paired_source_differences.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(out_root, "paired_method_differences.csv"), row.names = FALSE)

# Calibration is averaged across the six truth families inside each
# replication, then summarised across independent noise replications.
cal_rep <- aggregate(empirical ~ eta + replication + method + nominal,
                     calibration, mean)
cal_summary <- do.call(rbind, lapply(split(cal_rep,
  interaction(cal_rep$eta, cal_rep$method, cal_rep$nominal, drop = TRUE)), function(z) {
  data.frame(eta = z$eta[1], method = z$method[1], nominal = z$nominal[1],
             empirical_mean = mean(z$empirical),
             empirical_mcse = sd(z$empirical) / sqrt(nrow(z)),
             empirical_q025 = quantile(z$empirical, .025),
             empirical_q975 = quantile(z$empirical, .975))
}))
write.csv(cal_summary, file.path(out_root, "calibration_summary.csv"), row.names = FALSE)

colours <- c("Tikhonov discrepancy" = "#6B7280", "RBF" = "#440154",
             "Matern 0.5" = "#31688E", "Matern 1.5" = "#35B779",
             "Matern 2.5" = "#FDE725")
display <- c("Tikhonov discrepancy" = "Tikhonov", "RBF" = "RBF",
             "Matern 0.5" = "M0.5", "Matern 1.5" = "M1.5",
             "Matern 2.5" = "M2.5")
etas <- sort(unique(summary$eta)); methods <- names(colours)
png(file.path(figure_dir, "repeated_noise_rmse.png"), 2200, 900, res = 180)
old <- par(mfrow = c(1, length(etas)), mar = c(7, 4.5, 3, .8), las = 1)
for (eta in etas) {
  z <- summary[summary$eta == eta, ]; z <- z[match(methods, z$method), ]
  at <- seq_len(nrow(z)); ylim <- range(z$RMSE_lower, z$RMSE_upper, na.rm = TRUE)
  plot(at, z$RMSE_mean, ylim = ylim, xaxt = "n", xlab = "", ylab = "RMSE",
       pch = 19, col = colours[z$method], main = bquote(eta == .(eta)))
  arrows(at, z$RMSE_lower, at, z$RMSE_upper, angle = 90, code = 3,
         length = .04, col = colours[z$method], lwd = 1.5)
  axis(1, at, labels = display[z$method], las = 1, cex.axis = .82)
  grid(nx = NA, col = "grey90")
}
par(old); dev.off()

png(file.path(figure_dir, "repeated_noise_paired.png"), 2100, 900, res = 180)
old <- par(mfrow = c(1, length(etas)), mar = c(7, 4.5, 3, .8), las = 1)
for (eta in etas) {
  z <- paired_summary[paired_summary$eta == eta & paired_summary$method != "Matern 0.5", ]
  z <- z[order(z$mean_delta_RMSE), ]; at <- seq_len(nrow(z))
  lim <- range(c(z$lower95, z$upper95, 0), finite = TRUE)
  plot(at, z$mean_delta_RMSE, ylim = lim, xaxt = "n", xlab = "",
       ylab = if (eta == etas[1]) expression(Delta*" RMSE versus Matern 0.5") else "",
       pch = 19,
       col = colours[z$method], main = bquote(eta == .(eta)))
  abline(h = 0, lty = 2, col = "grey45")
  arrows(at, z$lower95, at, z$upper95, angle = 90, code = 3,
         length = .04, col = colours[z$method], lwd = 1.5)
  axis(1, at, labels = display[z$method], las = 1, cex.axis = .82)
  grid(nx = NA, col = "grey90")
}
par(old); dev.off()

png(file.path(figure_dir, "repeated_noise_calibration.png"), 2200, 900, res = 180)
old <- par(mfrow = c(1, length(etas)), mar = c(4.4, 4.5, 3, .8), las = 1)
for (eta in etas) {
  z <- cal_summary[cal_summary$eta == eta, ]
  plot(c(.45, 1), c(.45, 1), type = "n", xlab = "Nominal credibility",
       ylab = "Empirical coverage", main = bquote(eta == .(eta)))
  abline(0, 1, lty = 2, col = "grey45"); grid(col = "grey90")
  for (method in setdiff(methods, "Tikhonov discrepancy")) {
    zz <- z[z$method == method, ]; zz <- zz[order(zz$nominal), ]
    lines(zz$nominal, zz$empirical_mean, type = "b", pch = 19, lwd = 1.8,
          col = colours[method])
  }
  if (eta == etas[1]) legend("topleft", setdiff(methods, "Tikhonov discrepancy"),
    col = colours[setdiff(methods, "Tikhonov discrepancy")], lwd = 1.8,
    pch = 19, bty = "n", cex = .75)
}
par(old); dev.off()

fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, digits = d, format = "f"))
escape_tex <- function(x) gsub("%", "\\\\%", x, fixed = TRUE)
table_lines <- c("\\begin{table}[H]", "\\centering",
  paste0("\\caption{Repeated-noise robustness under multi-start empirical Bayes. ",
         "Entries are means across ", repetitions,
         " paired noise replications after averaging the six sources; parentheses give Monte Carlo standard errors.}"),
  "\\label{tab:repeated-noise}", "\\small", "\\begin{tabular}{clrrrr}",
  "\\toprule", "$\\eta$ & Method & RMSE & Coverage95 & Cal.MAE & CRPS\\\\", "\\midrule")
for (eta in etas) {
  z <- summary[summary$eta == eta, ]
  for (i in seq_len(nrow(z))) table_lines <- c(table_lines, sprintf(
    "$10^{%d}$ & %s & %s (%s) & %s (%s) & %s (%s) & %s (%s)\\\\",
    round(log10(eta)), escape_tex(z$method[i]), fmt(z$RMSE_mean[i], 4), fmt(z$RMSE_mcse[i], 4),
    fmt(z$Coverage95_mean[i]), fmt(z$Coverage95_mcse[i]),
    fmt(z$CalibrationMAE_mean[i]), fmt(z$CalibrationMAE_mcse[i]),
    fmt(z$CRPS_mean[i], 4), fmt(z$CRPS_mcse[i], 4)))
  if (eta != tail(etas, 1)) table_lines <- c(table_lines, "\\addlinespace")
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
table_dir <- file.path(out_root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(table_lines, file.path(table_dir, "repeated_noise.tex"))

paired_lines <- c("\\begin{table}[H]", "\\centering",
  "\\caption{Paired RMSE differences relative to Mat\\'ern $\\nu=0.5$. Negative values favour the row method; intervals are 95\\% $t$ intervals over noise replications.}",
  "\\label{tab:paired-rmse}", "\\small", "\\begin{tabular}{clrrr}",
  "\\toprule", "$\\eta$ & Method & Mean difference & Lower 95\\% & Upper 95\\%\\\\", "\\midrule")
for (i in seq_len(nrow(paired_summary))) if (paired_summary$method[i] != "Matern 0.5") {
  z <- paired_summary[i, ]; paired_lines <- c(paired_lines, sprintf(
    "$10^{%d}$ & %s & %s & %s & %s\\\\", round(log10(z$eta)),
    escape_tex(z$method), fmt(z$mean_delta_RMSE, 4), fmt(z$lower95, 4), fmt(z$upper95, 4)))
}
paired_lines <- c(paired_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(paired_lines, file.path(table_dir, "repeated_noise_paired.tex"))

message("Repeated-noise study complete: ", out_root)
