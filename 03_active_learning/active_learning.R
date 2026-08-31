
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

# 1. Configuration --------------------------------------------------------

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

no_plots <- isTRUE(getOption(
  "source.inversion.active.no_plots",
  as_flag(get_arg("no-plots", "false"))
))
source_names <- c("spike", "mixture", "dipole", "ring", "edge", "sine")
source_arg <- getOption("source.inversion.active.source",
                        get_arg("source", "all"))

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

# RStudio one-click entry. With no command-line arguments, this file runs all
# six sources in isolated environments inside the current R session and then
# assembles the combined tables and figures. No child Rscript process is used,
# which avoids Windows/RStudio process-launch and quoting failures.
if (identical(source_arg, "all")) {
  project_root <- normalizePath(dirname(script_dir), winslash = "/",
                                mustWork = TRUE)
  message("Running the formal budget-16 experiment for all six sources.")
  for (src in source_names) {
    message("\n[", match(src, source_names), "/", length(source_names),
            "] Source: ", src)
    old_options <- options(
      source.inversion.active.source = src,
      source.inversion.active.no_plots = no_plots
    )
    tryCatch(
      sys.source(file.path(script_dir, "active_learning.R"),
                 envir = new.env(parent = globalenv())),
      error = function(e) {
        stop("Active-learning experiment failed for source '", src,
             "': ", conditionMessage(e), call. = FALSE)
      },
      finally = options(old_options)
    )
  }

  message("\n[Final] Assembling six-source tables and figures.")
  sys.source(
    file.path(script_dir, "assemble_active_learning_6sources_report.R"),
    envir = new.env(parent = globalenv())
  )
  report_output <- file.path(project_root, "results",
                             "active_learning_6sources_report")
  message("\nAll six sources completed successfully.")
  message("Combined results: ", report_output)
} else {
source_name <- match.arg(source_arg, source_names)

cfg <- list(
  seed = 20161387,
  n_inverse = 15, n_truth = 41, n_modes = 36,
  n_candidate_side = 25, noise_ratio = 1e-2,
  initial_side = 2, budgets = 16,
  n_starts = 4, optim_maxit = 120,
  levels = c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.975, 0.99)
)

cfg$seed <- as.integer(get_arg("seed", cfg$seed))
cfg$noise_ratio <- as.numeric(get_arg("noise-ratio", cfg$noise_ratio))
if (!is.finite(cfg$noise_ratio) || cfg$noise_ratio <= 0) {
  stop("--noise-ratio must be positive.")
}
if (any(sqrt(cfg$budgets) %% 1 != 0) || min(cfg$budgets) <= cfg$initial_side^2) {
  stop("Every budget must be a square number larger than the initial design.")
}

noise_tag <- formatC(cfg$noise_ratio, format = "e", digits = 0)
source_folders <- c(
  spike = "Spike", mixture = "Broad_Plus_Sharp", dipole = "Dipole",
  ring = "Ring", edge = "Moving_Edge", sine = "Sine"
)
result_prefix <- "active_single_"
output_default <- file.path(dirname(script_dir), "results",
                            paste0(result_prefix, source_name))
out_dir <- normalizePath(get_arg("output", output_default), winslash = "/",
                         mustWork = FALSE)
dirs <- list(
  root = out_dir, figures = file.path(out_dir, "figures"),
  results = file.path(out_dir, "results")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

method_labels <- c(
  uniform = "Uniform",
  random = "Random",
  max_source_variance = "Max source variance",
  variance_times_gradient = "Variance x |gradient mean|"
)
method_colours <- c(
  uniform = "#440154", random = "#414487",
  max_source_variance = "#22A884",
  variance_times_gradient = "#FDE725"
)
sequential_methods <- setdiff(names(method_labels), "uniform")

# 2. Numerical helpers ----------------------------------------------------

safe_chol <- function(a, initial = 1e-12, max_tries = 10) {
  a <- (a + t(a)) / 2
  scale <- max(mean(abs(diag(a))), 1e-12)
  for (i in 0:max_tries) {
    jitter <- if (i == 0) 0 else initial * scale * 10^(i - 1)
    r <- try(chol(a + diag(jitter, nrow(a))), silent = TRUE)
    if (!inherits(r, "try-error")) return(r)
  }
  stop("Cholesky factorisation failed.")
}
chol_solve <- function(r, b) backsolve(r, forwardsolve(t(r), b))
rmse <- function(x, truth) sqrt(mean((x - truth)^2))
mae <- function(x, truth) mean(abs(x - truth))

regular_grid <- function(n, lower = NULL, upper = NULL) {
  axis <- if (is.null(lower)) (seq_len(n) - 0.5) / n else
    seq(lower, upper, length.out = n)
  expand.grid(x = axis, y = axis, KEEP.OUT.ATTRS = FALSE)
}

pairwise_sqdist <- function(a, b = a) {
  a <- as.matrix(a[, c("x", "y")]); b <- as.matrix(b[, c("x", "y")])
  pmax(outer(rowSums(a^2), rowSums(b^2), "+") - 2 * tcrossprod(a, b), 0)
}

matern05 <- function(d2, ell) exp(-sqrt(d2) / ell)

build_forward_matrix <- function(obs, source_grid, n_modes) {
  modes <- expand.grid(n = seq_len(n_modes), m = seq_len(n_modes),
                       KEEP.OUT.ATTRS = FALSE)
  beta <- pi^2 * (modes$n^2 + modes$m^2)
  phi <- function(points) {
    sx <- outer(points$x, modes$n, function(x, n) sin(n * pi * x))
    sy <- outer(points$y, modes$m, function(y, m) sin(m * pi * y))
    2 * sx * sy
  }
  weight <- 1 / nrow(source_grid)
  weight * sweep(phi(obs), 2, 1 / beta, "*") %*% t(phi(source_grid))
}

point_keys <- function(points) sprintf("%.10f_%.10f", points$x, points$y)
match_points <- function(points, catalogue) {
  answer <- match(point_keys(points), point_keys(catalogue))
  if (anyNA(answer)) stop("A requested design point is absent from candidate grid.")
  answer
}

# Dense bilinear interpolation matrix.  It is used only for acquisitions and
# trajectory diagnostics; final posterior fields use exact GP cross-covariance.
bilinear_matrix <- function(from_grid, to_grid) {
  xs <- sort(unique(from_grid$x)); ys <- sort(unique(from_grid$y))
  nx <- length(xs); ny <- length(ys)
  P <- matrix(0, nrow(to_grid), nrow(from_grid))
  index <- function(ix, iy) ix + (iy - 1) * nx
  for (r in seq_len(nrow(to_grid))) {
    x <- to_grid$x[r]; y <- to_grid$y[r]
    ix0 <- max(1, min(nx - 1, findInterval(x, xs)))
    iy0 <- max(1, min(ny - 1, findInterval(y, ys)))
    ix1 <- ix0 + 1; iy1 <- iy0 + 1
    tx <- (x - xs[ix0]) / (xs[ix1] - xs[ix0])
    ty <- (y - ys[iy0]) / (ys[iy1] - ys[iy0])
    weights <- c((1 - tx) * (1 - ty), tx * (1 - ty),
                 (1 - tx) * ty, tx * ty)
    columns <- c(index(ix0, iy0), index(ix1, iy0),
                 index(ix0, iy1), index(ix1, iy1))
    P[r, columns] <- weights
  }
  P
}

gradient_magnitude <- function(values, n_side, spacing) {
  z <- matrix(values, nrow = n_side, ncol = n_side)
  gx <- gy <- matrix(0, n_side, n_side)
  gx[2:(n_side - 1), ] <-
    (z[3:n_side, , drop = FALSE] - z[1:(n_side - 2), , drop = FALSE]) /
    (2 * spacing)
  gx[1, ] <- (z[2, ] - z[1, ]) / spacing
  gx[n_side, ] <- (z[n_side, ] - z[n_side - 1, ]) / spacing
  gy[, 2:(n_side - 1)] <-
    (z[, 3:n_side, drop = FALSE] - z[, 1:(n_side - 2), drop = FALSE]) /
    (2 * spacing)
  gy[, 1] <- (z[, 2] - z[, 1]) / spacing
  gy[, n_side] <- (z[, n_side] - z[, n_side - 1]) / spacing
  as.vector(sqrt(gx^2 + gy^2))
}

# 3. Problem and common design grids --------------------------------------

latent_grid <- regular_grid(cfg$n_inverse)
truth_grid <- regular_grid(cfg$n_truth)
candidate_grid <- regular_grid(cfg$n_candidate_side, 0.1, 0.9)
candidate_spacing <- 0.8 / (cfg$n_candidate_side - 1)

initial_axis <- seq(0.1, 0.9, length.out = cfg$initial_side)
initial_grid <- expand.grid(x = initial_axis, y = initial_axis,
                            KEEP.OUT.ATTRS = FALSE)
initial_idx <- match_points(initial_grid, candidate_grid)
uniform_idx <- setNames(lapply(cfg$budgets, function(budget) {
  side <- as.integer(round(sqrt(budget)))
  axis <- seq(0.1, 0.9, length.out = side)
  match_points(
    expand.grid(x = axis, y = axis, KEEP.OUT.ATTRS = FALSE),
    candidate_grid
  )
}), as.character(cfg$budgets))

source_catalogue <- list(
  spike = list(
    label = "Hidden narrow spike",
    fun = function(x, y) {
      exp(-((x - 0.73)^2 + (y - 0.27)^2) / (2 * 0.035^2))
    }
  ),
  mixture = list(
    label = "Broad plus sharp source",
    fun = function(x, y) {
      0.7 * exp(-((x - 0.35)^2 + (y - 0.65)^2) / (2 * 0.18^2)) +
        1.2 * exp(-((x - 0.75)^2 + (y - 0.25)^2) / (2 * 0.035^2))
    }
  ),
  dipole = list(
    label = "Positive-negative dipole",
    fun = function(x, y) {
      exp(-((x - 0.35)^2 + (y - 0.50)^2) / (2 * 0.07^2)) -
        exp(-((x - 0.65)^2 + (y - 0.50)^2) / (2 * 0.07^2))
    }
  ),
  ring = list(
    label = "Ring-shaped source",
    fun = function(x, y) {
      radius <- sqrt((x - 0.5)^2 + (y - 0.5)^2)
      exp(-((radius - 0.22)^2) / (2 * 0.035^2))
    }
  ),
  edge = list(
    label = "Moving-edge source",
    fun = function(x, y) {
      as.numeric(x > 0.45 + 0.15 * sin(4 * pi * y))
    }
  ),
  sine = list(
    label = "Sine-cosine source",
    fun = function(x, y) sin(pi * x) * cos(2 * pi * y)
  )
)

source_problem <- source_catalogue[[source_name]]
source_label <- source_problem$label
f_truth <- source_problem$fun(truth_grid$x, truth_grid$y)
source_is_signed <- min(f_truth) < -sqrt(.Machine$double.eps)

message("Building Green matrices...")
H_candidate_truth <- build_forward_matrix(candidate_grid, truth_grid,
                                            cfg$n_modes)
H_candidate <- build_forward_matrix(candidate_grid, latent_grid, cfg$n_modes)
u_clean_all <- as.vector(H_candidate_truth %*% f_truth)
sigma_noise <- cfg$noise_ratio * max(abs(u_clean_all))

D2_ll <- pairwise_sqdist(latent_grid)
D2_tl <- pairwise_sqdist(truth_grid, latent_grid)
P_candidate <- bilinear_matrix(latent_grid, candidate_grid)
P_truth <- bilinear_matrix(latent_grid, truth_grid)

# 4. Common initial MAP hyperparameters ----------------------------------

fit_initial_map <- function(y_all) {
  H0 <- H_candidate[initial_idx, , drop = FALSE]
  lower <- c(log_ell = log(0.015), log_sigma_f = log(0.02))
  upper <- c(log_ell = log(1.20), log_sigma_f = log(5.00))
  prior_mean <- c(log_ell = log(0.16), log_sigma_f = log(0.65))
  prior_sd <- c(log_ell = 0.80, log_sigma_f = 0.75)
  target <- function(theta) {
    K <- exp(2 * theta["log_sigma_f"]) *
      matern05(D2_ll, exp(theta["log_ell"]))
    S <- H0 %*% K %*% t(H0) + diag(sigma_noise^2, nrow(H0))
    R <- try(safe_chol(S), silent = TRUE)
    if (inherits(R, "try-error")) return(-1e100)
    y <- y_all[initial_idx]
    alpha <- chol_solve(R, y)
    ll <- -0.5 * (sum(y * alpha) + 2 * sum(log(diag(R))) +
                    length(y) * log(2 * pi))
    ll + sum(stats::dnorm(theta, prior_mean, prior_sd, log = TRUE))
  }
  objective <- function(theta) -target(theta)
  starts <- rbind(
    c(log_ell = log(0.16), log_sigma_f = log(0.65)),
    c(log_ell = log(0.05), log_sigma_f = log(0.65)),
    c(log_ell = log(0.40), log_sigma_f = log(0.65)),
    c(log_ell = log(0.12), log_sigma_f = log(1.20))
  )[seq_len(cfg$n_starts), , drop = FALSE]
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    stats::optim(starts[i, ], objective, method = "L-BFGS-B",
                 lower = lower, upper = upper,
                 control = list(maxit = cfg$optim_maxit, factr = 1e7))
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  list(theta = best$par, objective = best$value,
       convergence = best$convergence)
}

kernel_blocks <- function(theta) {
  ell <- exp(theta["log_ell"]); sigma_f2 <- exp(2 * theta["log_sigma_f"])
  list(
    Kll = sigma_f2 * matern05(D2_ll, ell),
    Ktl = sigma_f2 * matern05(D2_tl, ell),
    sigma_f2 = sigma_f2
  )
}

# 5. Sequential posterior and acquisition --------------------------------

initial_state <- function(theta, y_all) {
  blocks <- kernel_blocks(theta); K <- blocks$Kll
  H0 <- H_candidate[initial_idx, , drop = FALSE]
  S <- H0 %*% K %*% t(H0) + diag(sigma_noise^2, nrow(H0))
  R <- safe_chol(S)
  KHo <- K %*% t(H0)
  m <- as.vector(KHo %*% chol_solve(R, y_all[initial_idx]))
  root <- forwardsolve(t(R), t(KHo))
  C <- (K - crossprod(root) + t(K - crossprod(root))) / 2
  B <- C %*% t(H_candidate)
  mean_candidate <- as.vector(P_candidate %*% m)
  PC <- P_candidate %*% C
  var_candidate <- pmax(rowSums(PC * P_candidate), 0)
  obs_var <- pmax(rowSums(H_candidate * t(B)) + sigma_noise^2, 1e-14)
  list(m = m, C = C, B = B, mean_candidate = mean_candidate,
       var_candidate = var_candidate, obs_var = obs_var,
       selected = initial_idx)
}

select_next <- function(state, strategy, random_order, random_position) {
  remaining <- setdiff(seq_len(nrow(candidate_grid)), state$selected)
  if (strategy == "random") {
    while (random_position <= length(random_order) &&
           random_order[random_position] %in% state$selected) {
      random_position <- random_position + 1
    }
    return(list(index = random_order[random_position],
                random_position = random_position + 1))
  }

  score <- switch(
    strategy,
    max_source_variance = state$var_candidate,
    variance_times_gradient = state$var_candidate * gradient_magnitude(
      state$mean_candidate, cfg$n_candidate_side, candidate_spacing
    ),
    stop("Unknown strategy: ", strategy)
  )
  score[state$selected] <- -Inf
  list(index = which.max(score), random_position = random_position)
}

update_state <- function(state, index, observed_value) {
  c_source <- state$B[, index]
  denominator <- max(state$obs_var[index], 1e-14)
  h <- H_candidate[index, ]
  residual <- observed_value - sum(h * state$m)
  cross_observation <- as.vector(H_candidate %*% c_source)
  cross_candidate_source <- as.vector(P_candidate %*% c_source)

  state$m <- state$m + c_source * residual / denominator
  state$C <- state$C - tcrossprod(c_source) / denominator
  state$C <- (state$C + t(state$C)) / 2
  state$B <- state$B - tcrossprod(c_source, cross_observation) / denominator
  state$mean_candidate <- state$mean_candidate +
    cross_candidate_source * residual / denominator
  state$var_candidate <- pmax(
    state$var_candidate - cross_candidate_source^2 / denominator, 0
  )
  state$obs_var <- pmax(
    state$obs_var - cross_observation^2 / denominator, sigma_noise^2
  )
  state$selected <- c(state$selected, index)
  state
}

run_path <- function(strategy, theta, y_all) {
  state <- initial_state(theta, y_all)
  set.seed(cfg$seed + 7919 + match(strategy, sequential_methods))
  random_order <- sample(setdiff(seq_len(nrow(candidate_grid)), initial_idx))
  random_position <- 1
  trajectory <- data.frame(
    n = length(state$selected),
    RMSE = rmse(as.vector(P_truth %*% state$m), f_truth)
  )
  snapshots <- list()
  while (length(state$selected) < max(cfg$budgets)) {
    choice <- select_next(state, strategy, random_order, random_position)
    random_position <- choice$random_position
    state <- update_state(state, choice$index, y_all[choice$index])
    trajectory <- rbind(trajectory, data.frame(
      n = length(state$selected),
      RMSE = rmse(as.vector(P_truth %*% state$m), f_truth)
    ))
    if (length(state$selected) %in% cfg$budgets) {
      snapshots[[as.character(length(state$selected))]] <- state$selected
    }
  }
  list(indices = snapshots, trajectory = trajectory)
}

posterior_truth <- function(theta, indices, y_all) {
  blocks <- kernel_blocks(theta)
  H <- H_candidate[indices, , drop = FALSE]
  Klo <- blocks$Kll %*% t(H)
  Kto <- blocks$Ktl %*% t(H)
  S <- H %*% Klo + diag(sigma_noise^2, nrow(H))
  R <- safe_chol(S)
  mean <- as.vector(Kto %*% chol_solve(R, y_all[indices]))
  root <- forwardsolve(t(R), t(Kto))
  variance <- pmax(blocks$sigma_f2 - colSums(root^2), 0)
  latent_root <- forwardsolve(t(R), t(Klo))
  latent_variance <- pmax(blocks$sigma_f2 - colSums(latent_root^2), 0)
  list(mean = mean, variance = variance, sd = sqrt(variance),
       latent_variance = latent_variance)
}

score_posterior <- function(posterior) {
  sd <- pmax(posterior$sd, 1e-10)
  z <- (f_truth - posterior$mean) / sd
  crps <- sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
  nlpd <- 0.5 * log(2 * pi * sd^2) +
    0.5 * (f_truth - posterior$mean)^2 / sd^2
  coverage <- vapply(cfg$levels, function(level) {
    quantile <- qnorm((1 + level) / 2)
    mean(abs(f_truth - posterior$mean) <= quantile * sd)
  }, numeric(1))
  j95 <- which.min(abs(cfg$levels - 0.95))
  list(
    metrics = c(
      RMSE = rmse(posterior$mean, f_truth),
      MAE = mae(posterior$mean, f_truth),
      Coverage95 = coverage[j95],
      CalibrationMAE = mean(abs(coverage - cfg$levels)),
      CRPS = mean(crps), NLPD = mean(nlpd), MeanSD = mean(sd),
      MeanVariance = mean(posterior$variance),
      LatentIntegratedVariance = mean(posterior$latent_variance)
    ),
    calibration = data.frame(nominal = cfg$levels, empirical = coverage)
  )
}

# 6. One fixed noisy dataset ---------------------------------------------

message("Generating the fixed noisy dataset...")
set.seed(cfg$seed)
noise_all <- rnorm(nrow(candidate_grid), 0, sigma_noise)
y_all <- u_clean_all + noise_all
fit <- fit_initial_map(y_all)
theta <- fit$theta
hyperparameters <- data.frame(
  source = source_name, source_label = source_label,
  ell = exp(theta["log_ell"]), sigma_f = exp(theta["log_sigma_f"]),
  convergence = fit$convergence
)

paths <- lapply(sequential_methods, function(method) {
  run_path(method, theta, y_all)
})
names(paths) <- sequential_methods
trajectories <- do.call(rbind, lapply(sequential_methods, function(method) {
  data.frame(
    source = source_name, source_label = source_label,
    method = method, label = unname(method_labels[method]),
    paths[[method]]$trajectory
  )
}))

metric_rows <- list(); calibration_rows <- list(); representative <- list()
row_id <- calibration_id <- 1
for (budget in cfg$budgets) {
  designs <- c(
    list(uniform = uniform_idx[[as.character(budget)]]),
    lapply(paths, function(path) path$indices[[as.character(budget)]])
  )
  for (method in names(method_labels)) {
    indices <- designs[[method]]
    posterior <- posterior_truth(theta, indices, y_all)
    scores <- score_posterior(posterior)
    one_metrics <- scores$metrics
    metric_rows[[row_id]] <- data.frame(
      source = source_name, source_label = source_label,
      budget = budget, method = method,
      label = unname(method_labels[method]), as.list(one_metrics),
      stringsAsFactors = FALSE
    )
    calibration_rows[[calibration_id]] <- data.frame(
      source = source_name, source_label = source_label,
      budget = budget, method = method,
      label = unname(method_labels[method]), scores$calibration,
      stringsAsFactors = FALSE
    )
    representative[[paste(method, budget, sep = "::")]] <- list(
      indices = indices, posterior = posterior, metrics = one_metrics
    )
    row_id <- row_id + 1
    calibration_id <- calibration_id + 1
  }
}
metrics <- do.call(rbind, metric_rows)
calibration <- do.call(rbind, calibration_rows)

uniform_comparison <- merge(
  metrics[, c("budget", "method", "label", "RMSE")],
  metrics[metrics$method == "uniform", c("budget", "RMSE")],
  by = "budget", suffixes = c("", "_uniform")
)
uniform_comparison <- uniform_comparison[uniform_comparison$method != "uniform", ]
uniform_comparison$RMSE_difference <-
  uniform_comparison$RMSE - uniform_comparison$RMSE_uniform
uniform_comparison$relative_RMSE_gain_percent <-
  100 * (uniform_comparison$RMSE_uniform - uniform_comparison$RMSE) /
  uniform_comparison$RMSE_uniform

# 7. Output tables --------------------------------------------------------

write.csv(metrics, file.path(dirs$results, "single_dataset_metrics.csv"),
          row.names = FALSE)
write.csv(trajectories, file.path(dirs$results, "rmse_trajectories.csv"),
          row.names = FALSE)
write.csv(hyperparameters, file.path(dirs$results, "initial_map_hyperparameters.csv"),
          row.names = FALSE)
write.csv(calibration, file.path(dirs$results, "calibration.csv"),
          row.names = FALSE)
write.csv(uniform_comparison,
          file.path(dirs$results, "uniform_comparison.csv"),
          row.names = FALSE)
write.csv(data.frame(
  parameter = c("source", names(cfg), "sigma_noise", "max_abs_u"),
  value = c(source_name,
            vapply(cfg, function(x) paste(x, collapse = ","), character(1)),
            sigma_noise, max(abs(u_clean_all)))
), file.path(dirs$results, "configuration.csv"), row.names = FALSE)

# 8. Figures --------------------------------------------------------------

field_plot <- function(values, title, option = "D", limits = NULL,
                       sensors = NULL, signed = FALSE) {
  df <- data.frame(truth_grid, value = values)
  fill_scale <- if (signed) {
    scale_fill_gradient2(
      low = "#440154", mid = "white", high = "#FDE725", midpoint = 0,
      limits = limits, oob = scales::squish, name = "value"
    )
  } else {
    scale_fill_viridis_c(option = option, limits = limits,
                         oob = scales::squish, name = "value")
  }
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    {if (is.null(sensors)) NULL else geom_point(
      data = sensors, aes(x, y), inherit.aes = FALSE,
      shape = 21, fill = "white", colour = "black", stroke = 0.25, size = 1
    )} +
    fill_scale +
    coord_equal(expand = FALSE) +
    labs(title = title, x = "x", y = "y") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
          legend.key.height = unit(0.60, "cm"),
          legend.key.width = unit(0.22, "cm"))
}

if (!no_plots) {
  uniform_points <- metrics[metrics$method == "uniform", ]
  p_trace <- ggplot(trajectories, aes(n, RMSE, colour = method)) +
    geom_line(linewidth = 0.75) +
    geom_point(data = uniform_points,
               aes(x = budget, y = RMSE, colour = method),
               inherit.aes = FALSE, size = 2.8, shape = 18) +
    scale_colour_manual(values = method_colours, labels = method_labels) +
    geom_vline(xintercept = cfg$budgets, linetype = 3, colour = "grey60") +
    labs(title = paste0(source_label, ": RMSE versus observation budget"),
         subtitle = "Lines: sequential designs; diamonds: fixed uniform grids",
         x = "Number of observations", y = "RMSE", colour = "Design") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(dirs$figures, "rmse_trajectory.png"), p_trace,
         width = 10.5, height = 7, dpi = 180, bg = "white")

  p_gain <- ggplot(uniform_comparison,
                   aes(label, relative_RMSE_gain_percent, fill = method)) +
    geom_hline(yintercept = 0, colour = "grey35") +
    geom_col(width = 0.68) +
    facet_wrap(~paste0("Budget = ", budget), scales = "free_y") +
    scale_fill_manual(values = method_colours, guide = "none") +
    labs(title = "RMSE gain relative to the uniform design",
         subtitle = "Positive values favour the alternative design",
         x = NULL, y = "Relative RMSE gain (%)") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  ggsave(file.path(dirs$figures, "relative_rmse_gain.png"), p_gain,
         width = 11, height = 6.5, dpi = 180, bg = "white")

  p_cov <- ggplot(metrics, aes(label, Coverage95, fill = method)) +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey35") +
    geom_col(width = 0.68) +
    facet_wrap(~paste0("Budget = ", budget)) +
    scale_fill_manual(values = method_colours, guide = "none") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Empirical spatial 95% coverage",
         x = NULL, y = "Coverage") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  ggsave(file.path(dirs$figures, "coverage95.png"), p_cov,
         width = 11, height = 6.5, dpi = 180, bg = "white")

  p_cal <- ggplot(calibration,
                  aes(nominal, empirical, colour = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                colour = "grey45") +
    geom_line(linewidth = 0.75) + geom_point(size = 1.5) +
    scale_colour_manual(values = method_colours, labels = method_labels) +
    scale_fill_manual(values = method_colours, guide = "none") +
    coord_equal(xlim = c(min(cfg$levels), 1), ylim = c(0, 1)) +
    facet_wrap(~paste0("Budget = ", budget)) +
    labs(title = paste(source_label, "spatial calibration"),
         x = "Nominal credibility", y = "Empirical coverage",
         colour = "Design") +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(dirs$figures, "calibration_curves.png"), p_cal,
         width = 8.5, height = 7, dpi = 180, bg = "white")

  for (budget in cfg$budgets) {
    keys <- paste(names(method_labels), budget, sep = "::")
    posts <- lapply(keys, function(key) representative[[key]]$posterior)
    names(posts) <- names(method_labels)
    mean_values <- c(f_truth, unlist(lapply(posts, `[[`, "mean")))
    mean_limit <- if (source_is_signed) {
      rep(max(abs(mean_values)), 2) * c(-1, 1)
    } else {
      range(c(0, mean_values))
    }
    truth_limit <- if (source_is_signed) {
      rep(max(abs(f_truth)), 2) * c(-1, 1)
    } else {
      range(c(0, f_truth))
    }
    error_limit <- c(0, max(unlist(lapply(posts, function(p) {
      abs(f_truth - p$mean)
    }))))
    sd_limit <- c(0, max(unlist(lapply(posts, `[[`, "sd"))))
    panels <- list()
    for (method in names(method_labels)) {
      item <- representative[[paste(method, budget, sep = "::")]]
      panels[[length(panels) + 1]] <- field_plot(
        item$posterior$mean,
        sprintf("%s mean\nRMSE=%.4f", method_labels[method],
                item$metrics["RMSE"]), "D", mean_limit,
        signed = source_is_signed
      )
      panels[[length(panels) + 1]] <- field_plot(
        abs(f_truth - item$posterior$mean),
        paste(method_labels[method], "|error|"), "A", error_limit
      )
      panels[[length(panels) + 1]] <- field_plot(
        item$posterior$sd, paste(method_labels[method], "posterior SD"),
        "A", sd_limit
      )
    }
    ggsave(
      file.path(dirs$figures, paste0("representative_fields_budget_", budget,
                                     ".png")),
      arrangeGrob(grobs = panels, ncol = 3,
                  top = grid::textGrob(
                     paste(source_label, "fixed dataset, budget", budget),
                    gp = grid::gpar(fontface = "bold")
                  )),
      width = 14, height = 19, dpi = 170, bg = "white"
    )

    design_panels <- lapply(names(method_labels), function(method) {
      item <- representative[[paste(method, budget, sep = "::")]]
      field_plot(
        f_truth, paste(method_labels[method], "sensor locations"),
        "D", truth_limit, candidate_grid[item$indices, , drop = FALSE],
        signed = source_is_signed
      )
    })
    ggsave(
      file.path(dirs$figures, paste0("representative_designs_budget_", budget,
                                     ".png")),
      arrangeGrob(grobs = design_panels, ncol = 2,
                  top = grid::textGrob(
                    paste("Selected observation locations, budget", budget),
                    gp = grid::gpar(fontface = "bold")
                  )),
      width = 13.5, height = 9, dpi = 180, bg = "white"
    )
  }
}

saveRDS(list(
  cfg = cfg, source = source_name, source_label = source_label,
  truth_grid = truth_grid, candidate_grid = candidate_grid,
  f_truth = f_truth, sigma_noise = sigma_noise,
  metrics = metrics,
  trajectories = trajectories, hyperparameters = hyperparameters,
  calibration = calibration, uniform_comparison = uniform_comparison,
  representative = representative
), file.path(dirs$results, "complete_experiment.rds"))
capture.output(sessionInfo(), file = file.path(dirs$results, "session_info.txt"))

message("Completed: ", dirs$root)
print(metrics[order(metrics$budget, metrics$RMSE),
              c("budget", "label", "RMSE", "MAE", "Coverage95",
                "CalibrationMAE", "CRPS")])
message("\nRelative RMSE gain versus uniform:")
print(uniform_comparison[
  order(uniform_comparison$budget,
        -uniform_comparison$relative_RMSE_gain_percent), ])
}
