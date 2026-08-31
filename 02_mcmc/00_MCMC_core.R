suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Run install.packages('ggplot2').")
  }
})

# 1. Configuration ---------------------------------------------------------

cli <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- cli[startsWith(cli, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
as_flag <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

script_location <- function() {
  # When this core is sourced by a source-specific wrapper, the sourced file
  # path is more reliable than commandArgs(--file), which still names the
  # wrapper and may be relative to the wrapper's temporary working directory.
  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  if (length(frame_files)) {
    return(dirname(normalizePath(tail(frame_files, 1)[[1]],
                                  winslash = "/", mustWork = FALSE)))
  }
  full <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", full[grepl("^--file=", full)])
  if (length(hit)) {
    return(dirname(normalizePath(hit[1], winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

MASTER_SEED <- as.integer(get_arg("seed", "20161387"))
force_refit <- as_flag(get_arg("force", "false"))
no_plots <- as_flag(get_arg("no-plots", "false"))

cfg <- list(
  seed = MASTER_SEED,
  n_inverse = 15,
  n_truth = 41,
  n_modes = 36,
  n_obs_side = 7,
  noise_ratio = 1e-2,
  map_starts = 6,
  map_maxit = 180,
  rwmh_iter = 4000,
  rwmh_burn = 1200,
  rwmh_thin = 4,
  n_chains = 4,
  mixture_support = 480,
  mixture_mc_draws = 800,
  stability_sizes = c(120, 240, 480),
  credibility = c(.50, .60, .70, .80, .90, .95, .975, .99)
)

# A source-specific entry script may request a diagnostics-driven extension.
# The spike wrapper does this because its first four-chain run did not meet the
# rank-normalized Rhat < 1.01 criterion for a peak-nearby downstream mean.
if (exists("ARWM_ITER_OVERRIDE", inherits = TRUE)) {
  cfg$rwmh_iter <- as.integer(get("ARWM_ITER_OVERRIDE", inherits = TRUE))
}
if (exists("ARWM_BURN_OVERRIDE", inherits = TRUE)) {
  cfg$rwmh_burn <- as.integer(get("ARWM_BURN_OVERRIDE", inherits = TRUE))
}

cfg$noise_ratio <- as.numeric(get_arg("noise-ratio", cfg$noise_ratio))
cfg$rwmh_iter <- as.integer(get_arg("rwmh-iter", cfg$rwmh_iter))
cfg$rwmh_burn <- as.integer(get_arg("rwmh-burn", cfg$rwmh_burn))
cfg$rwmh_thin <- as.integer(get_arg("rwmh-thin", cfg$rwmh_thin))

if (!is.finite(cfg$seed)) stop("--seed must be a finite integer.")
if (!is.finite(cfg$noise_ratio) || cfg$noise_ratio <= 0) {
  stop("--noise-ratio must be positive, for example 1e-2.")
}
if (cfg$rwmh_iter <= cfg$rwmh_burn || cfg$rwmh_burn < 10L ||
    cfg$rwmh_thin < 1) stop("Invalid RWMH iteration/burn/thin settings.")
retained_per_chain <- length(seq.int(cfg$rwmh_burn + 1, cfg$rwmh_iter,
                                     by = cfg$rwmh_thin))
if (retained_per_chain * cfg$n_chains < max(cfg$stability_sizes)) {
  stop("The four chains must retain at least 480 draws in total.")
}

script_dir <- script_location()
if (!exists("SOURCE_KEY", inherits = TRUE)) {
  stop("Run one of 01_MCMC_spike.R through 06_MCMC_sine.R, not the core directly.")
}
SOURCE_KEY <- get("SOURCE_KEY", inherits = TRUE)
output_dir <- normalizePath(get_arg(
  "output", file.path(dirname(script_dir), "results",
                      paste0("MCMC_results_", SOURCE_KEY))
), winslash = "/", mustWork = FALSE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Stable linear algebra and diagnostics ---------------------------------

regular_grid <- function(n) {
  axis <- (seq_len(n) - 0.5) / n
  expand.grid(x = axis, y = axis, KEEP.OUT.ATTRS = FALSE)
}

pairwise_sqdist <- function(a, b = a) {
  a <- as.matrix(a[, c("x", "y")])
  b <- as.matrix(b[, c("x", "y")])
  pmax(outer(rowSums(a^2), rowSums(b^2), "+") - 2 * tcrossprod(a, b), 0)
}

safe_chol <- function(A, initial = 1e-12, max_tries = 9) {
  A <- (A + t(A)) / 2
  scale <- max(mean(abs(diag(A))), 1e-12)
  for (i in 0:max_tries) {
    jitter <- if (i == 0L) 0 else initial * scale * 10^(i - 1)
    R <- try(chol(A + diag(jitter, nrow(A))), silent = TRUE)
    if (!inherits(R, "try-error")) return(R)
  }
  stop("Cholesky factorisation failed after adaptive jitter.")
}

chol_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))
log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
rmse <- function(x, truth) sqrt(mean((x - truth)^2))
mae <- function(x, truth) mean(abs(x - truth))

simple_ess <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 4 || !is.finite(sd(x)) || sd(x) == 0) return(as.numeric(n))
  rho <- as.numeric(acf(x, plot = FALSE, lag.max = min(100, n - 1),
                        demean = TRUE)$acf)[-1]
  if (!length(rho)) return(as.numeric(n))
  pairs <- rho[seq(1, length(rho), by = 2)]
  even <- rho[seq(2, length(rho), by = 2)]
  if (length(even)) pairs[seq_along(even)] <- pairs[seq_along(even)] + even
  first_bad <- which(pairs <= 0)[1]
  keep <- if (is.na(first_bad)) length(pairs) else first_bad - 1
  if (keep <= 0L) return(as.numeric(n))
  tau <- max(1 + 2 * sum(rho[seq_len(min(2 * keep, length(rho)))]), 1)
  min(n / tau, n)
}

split_chains <- function(chains) {
  chains <- as.matrix(chains)
  half <- floor(nrow(chains) / 2)
  if (half < 4) stop("Each retained chain needs at least eight draws.")
  cbind(chains[seq_len(half), , drop = FALSE],
        chains[(nrow(chains) - half + 1):nrow(chains), , drop = FALSE])
}

basic_rhat <- function(chains) {
  chains <- as.matrix(chains)
  n <- nrow(chains)
  chain_means <- colMeans(chains)
  W <- mean(apply(chains, 2, var))
  B <- n * var(chain_means)
  if (!is.finite(W) || W <= 0) return(1)
  sqrt(((n - 1) * W / n + B / n) / W)
}

rank_transform <- function(x) {
  z <- rank(x, ties.method = "average")
  qnorm((z - 3 / 8) / (length(z) + 1 / 4))
}

# Vehtari et al.'s rank-normalized split-Rhat: the maximum of the ordinary
# rank-normalized and folded rank-normalized split diagnostics.
rank_normalized_rhat <- function(chains) {
  split <- split_chains(chains)
  ranked <- matrix(rank_transform(as.vector(split)), nrow(split), ncol(split))
  folded_raw <- abs(split - median(split))
  folded <- matrix(rank_transform(as.vector(folded_raw)),
                   nrow(split), ncol(split))
  max(basic_rhat(ranked), basic_rhat(folded))
}

# Multi-chain initial-positive-sequence ESS applied to rank-normalized split
# chains. This is the bulk ESS used for the diagnostics and MCSE calculation.
rank_normalized_ess <- function(chains) {
  split <- split_chains(chains)
  z <- matrix(rank_transform(as.vector(split)), nrow(split), ncol(split))
  n <- nrow(z); m <- ncol(z)
  W <- mean(apply(z, 2, var))
  B <- n * var(colMeans(z))
  var_plus <- (n - 1) * W / n + B / n
  if (!is.finite(var_plus) || var_plus <= 0) return(as.numeric(n * m))
  max_lag <- min(n - 1, 100L)
  acov <- vapply(seq_len(m), function(j) {
    as.numeric(acf(z[, j], type = "covariance", plot = FALSE,
                   lag.max = max_lag, demean = TRUE)$acf)
  }, numeric(max_lag + 1))
  mean_acov <- rowMeans(acov)
  rho <- c(1, 1 - (W - mean_acov[-1]) / var_plus)
  pair_start <- seq(1, length(rho) - 1, by = 2)
  pairs <- rho[pair_start] + rho[pair_start + 1]
  first_negative <- which(pairs < 0)[1]
  if (!is.na(first_negative)) pairs <- pairs[seq_len(first_negative - 1)]
  if (!length(pairs)) return(as.numeric(n * m))
  if (length(pairs) > 1) {
    for (k in 2:length(pairs)) pairs[k] <- min(pairs[k], pairs[k - 1])
  }
  tau <- max(-1 + 2 * sum(pairs), 1)
  min(n * m / tau, n * m)
}

balanced_pool <- function(chain_list, total_draws) {
  if (total_draws %% length(chain_list) != 0L) {
    stop("total_draws must be divisible by the number of chains.")
  }
  n_available <- min(vapply(chain_list, function(z) nrow(z$draws), integer(1)))
  n_each <- total_draws %/% length(chain_list)
  if (n_each > n_available) stop("Not enough retained draws in each chain.")
  index <- unique(round(seq(1, n_available, length.out = n_each)))
  selected <- lapply(chain_list, function(z) z$draws[index, , drop = FALSE])
  do.call(rbind, lapply(seq_len(n_each), function(i) {
    do.call(rbind, lapply(selected, function(z) z[i, , drop = FALSE]))
  }))
}

# 3. Six source terms and the Dirichlet Green forward map ------------------

source_catalogue <- function() {
  list(
    spike = list(
      label = "Hidden narrow spike",
      fun = function(x, y) exp(-((x - .73)^2 + (y - .27)^2) /
                                  (2 * .035^2))
    ),
    mixture = list(
      label = "Broad plus sharp source",
      fun = function(x, y) {
        .7 * exp(-((x - .35)^2 + (y - .65)^2) / (2 * .18^2)) +
          1.2 * exp(-((x - .75)^2 + (y - .25)^2) / (2 * .035^2))
      }
    ),
    dipole = list(
      label = "Positive-negative dipole",
      fun = function(x, y) {
        exp(-((x - .35)^2 + (y - .50)^2) / (2 * .07^2)) -
          exp(-((x - .65)^2 + (y - .50)^2) / (2 * .07^2))
      }
    ),
    ring = list(
      label = "Ring-shaped source",
      fun = function(x, y) {
        radius <- sqrt((x - .5)^2 + (y - .5)^2)
        exp(-(radius - .22)^2 / (2 * .035^2))
      }
    ),
    edge = list(
      label = "Moving-edge source",
      fun = function(x, y) as.numeric(x > .45 + .15 * sin(4 * pi * y))
    ),
    sine = list(
      label = "Sine-cosine source",
      fun = function(x, y) sin(pi * x) * cos(2 * pi * y)
    )
  )
}

# Dirichlet eigenfunctions phi_nm=2 sin(n*pi*x) sin(m*pi*y) are orthonormal
# on the unit square and -Delta phi_nm=lambda_nm phi_nm.  Therefore
# G(z,x)=sum_nm phi_nm(z)phi_nm(x)/lambda_nm.  The final division by the
# number of source-grid cells is the midpoint quadrature weight.
build_forward_matrix <- function(obs, source_grid, n_modes) {
  modes <- expand.grid(n = seq_len(n_modes), m = seq_len(n_modes),
                       KEEP.OUT.ATTRS = FALSE)
  lambda <- pi^2 * (modes$n^2 + modes$m^2)
  basis <- function(points) {
    2 * sin(pi * outer(points$x, modes$n)) *
      sin(pi * outer(points$y, modes$m))
  }
  sweep(basis(obs), 2, lambda, "/") %*% t(basis(source_grid)) /
    nrow(source_grid)
}

latent_grid <- regular_grid(cfg$n_inverse)
truth_grid <- regular_grid(cfg$n_truth)
obs_axis <- seq(.1, .9, length.out = cfg$n_obs_side)
obs_grid <- expand.grid(x = obs_axis, y = obs_axis, KEEP.OUT.ATTRS = FALSE)

message("Constructing Green matrices...")
H_inverse <- build_forward_matrix(obs_grid, latent_grid, cfg$n_modes)
H_truth <- build_forward_matrix(obs_grid, truth_grid, cfg$n_modes)

catalogue <- source_catalogue()
if (length(SOURCE_KEY) != 1 || !SOURCE_KEY %in% names(catalogue)) {
  stop("SOURCE_KEY must be one of: ", paste(names(catalogue), collapse = ", "))
}
catalogue <- catalogue[SOURCE_KEY]

make_problem <- function(key, index) {
  item <- catalogue[[key]]
  truth <- item$fun(truth_grid$x, truth_grid$y)
  u_clean <- as.vector(H_truth %*% truth)
  # The stated experimental convention: sigma_n=eta*max_z |u(z)|.
  sigma_noise <- max(cfg$noise_ratio * max(abs(u_clean)), 1e-10)
  set.seed(cfg$seed + index * 1009)
  list(
    key = key,
    label = item$label,
    truth = truth,
    u_clean = u_clean,
    y = u_clean + rnorm(length(u_clean), 0, sigma_noise),
    sigma_noise = sigma_noise
  )
}
problems <- Map(make_problem, names(catalogue), seq_along(catalogue))
names(problems) <- names(catalogue)

# 4. Matern-0.5 source prior and marginal hyperposterior -------------------

D2_ll <- pairwise_sqdist(latent_grid)
D2_tl <- pairwise_sqdist(truth_grid, latent_grid)

matern05_covariance <- function(d2, theta) {
  ell <- exp(unname(theta["log_ell"]))
  sigma_f <- exp(unname(theta["log_sigma_f"]))
  sigma_f^2 * exp(-sqrt(d2) / ell)
}

parameter_info <- function(problem) {
  y_scale <- max(sd(problem$y), sqrt(mean(problem$y^2)), 1e-8)
  noise_lower <- log(max(problem$sigma_noise * .01, 1e-10))
  noise_upper <- log(max(problem$sigma_noise * 20, y_scale, 1e-7))
  init <- c(log_ell = log(.18), log_sigma_f = log(.60),
            log_sigma_n = log(problem$sigma_noise))
  lower <- c(log_ell = log(.015), log_sigma_f = log(.015),
             log_sigma_n = noise_lower)
  upper <- c(log_ell = log(1.20), log_sigma_f = log(5),
             log_sigma_n = noise_upper)
  prior_mean <- init
  prior_sd <- c(log_ell = .90, log_sigma_f = .70, log_sigma_n = .35)
  list(init = init, lower = lower, upper = upper,
       prior_mean = prior_mean, prior_sd = prior_sd)
}

log_marginal <- function(theta, y) {
  ans <- tryCatch({
    K <- matern05_covariance(D2_ll, theta)
    S <- H_inverse %*% K %*% t(H_inverse) +
      diag(exp(2 * theta["log_sigma_n"]), length(y))
    R <- safe_chol(S)
    alpha <- chol_solve(R, y)
    -.5 * sum(y * alpha) - sum(log(diag(R))) - length(y) * log(2 * pi) / 2
  }, error = function(e) -Inf)
  if (is.finite(ans)) ans else -Inf
}

log_prior <- function(theta, info) {
  if (any(theta < info$lower) || any(theta > info$upper) ||
      any(!is.finite(theta))) return(-Inf)
  sum(dnorm(theta, info$prior_mean, info$prior_sd, log = TRUE))
}

log_target_factory <- function(problem, info) {
  function(theta) {
    names(theta) <- names(info$init)
    ll <- log_marginal(theta, problem$y)
    if (!is.finite(ll)) return(-Inf)
    ll + log_prior(theta, info)
  }
}

# Multi-start MAP is used only to initialise both chains at the same sensible
# point.  It is not used as the final Full Bayesian estimate.
map_initialise <- function(problem, info, seed) {
  set.seed(seed)
  starts <- matrix(info$init, cfg$map_starts, length(info$init), byrow = TRUE,
                   dimnames = list(NULL, names(info$init)))
  if (cfg$map_starts >= 2) starts[2, "log_ell"] <- log(.45)
  if (cfg$map_starts >= 3) {
    for (i in 3:cfg$map_starts) {
      starts[i, ] <- info$prior_mean +
        rnorm(length(info$init), 0, .65 * info$prior_sd)
    }
  }
  starts <- pmax(starts, matrix(info$lower + 1e-9, nrow(starts),
                                length(info$lower), byrow = TRUE))
  starts <- pmin(starts, matrix(info$upper - 1e-9, nrow(starts),
                                length(info$upper), byrow = TRUE))
  target <- log_target_factory(problem, info)
  objective <- function(theta) {
    value <- target(theta)
    if (is.finite(value)) -value else 1e100
  }
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    tryCatch(optim(starts[i, ], objective, method = "L-BFGS-B",
                   lower = info$lower, upper = info$upper,
                   control = list(maxit = cfg$map_maxit, factr = 1e7)),
             error = function(e) list(par = starts[i, ], value = Inf,
                                      convergence = 999))
  })
  values <- vapply(fits, function(z) -z$value, numeric(1))
  if (!any(is.finite(values))) stop("All MAP starts failed for ", problem$key)
  best <- fits[[which.max(values)]]
  theta <- best$par
  names(theta) <- names(info$init)
  list(theta = theta, log_target = max(values), starts = starts,
       start_targets = values, convergence = best$convergence)
}

# 5. Adaptive random-walk Metropolis ---------------------------------------

# Adaptation is stopped at burn-in.  Hence all retained iterations use a fixed
# Gaussian proposal and form an ordinary Markov chain with the desired target.
run_rwmh <- function(target, init, info, seed) {
  start_time <- proc.time()[3]
  set.seed(seed)
  p <- length(init)
  chain <- matrix(NA_real_, cfg$rwmh_iter, p,
                  dimnames = list(NULL, names(init)))
  trace <- numeric(cfg$rwmh_iter)
  accepted <- logical(cfg$rwmh_iter)
  current <- init
  current_lp <- target(current)
  proposal_base <- diag(c(.14, .12, .09)^2, p)
  proposal_cov <- proposal_base
  Rprop <- chol(proposal_cov)
  # Conventional optimal-scaling reference for a p-dimensional Gaussian
  # random-walk proposal. Robbins-Monro subsequently adapts this scalar,
  # while the empirical covariance adapts the proposal shape. There is no
  # shrinkage blend with proposal_base in the adaptive covariance.
  log_scale <- log(2.38^2 / p)
  adaptation_count <- 0L
  adaptation_history <- data.frame(
    iteration = integer(), batch = integer(), batch_acceptance = numeric(),
    gain = numeric(), log_scale = numeric(), scale = numeric(),
    min_eigenvalue = numeric(), max_eigenvalue = numeric()
  )
  chain[1, ] <- current
  trace[1] <- current_lp

  for (iteration in 2:cfg$rwmh_iter) {
    proposal <- current + drop(t(Rprop) %*% rnorm(p))
    names(proposal) <- names(init)
    proposal_lp <- if (any(proposal < info$lower | proposal > info$upper)) {
      -Inf
    } else target(proposal)
    log_ratio <- proposal_lp - current_lp
    if (is.finite(log_ratio) && log(runif(1)) < min(0, log_ratio)) {
      current <- proposal
      current_lp <- proposal_lp
      accepted[iteration] <- TRUE
    }
    chain[iteration, ] <- current
    trace[iteration] <- current_lp

    # Two standard warm-up adaptations:
    #   (i) empirical covariance learns posterior shape/correlation;
    #  (ii) Robbins-Monro learns the global random-walk scale.
    # Both stop at the warm-up boundary, after which proposal_cov is frozen.
    if (iteration <= cfg$rwmh_burn && iteration >= 100L &&
        iteration %% 50L == 0L) {
      adaptation_count <- adaptation_count + 1
      window <- (iteration - 49):iteration
      gain <- adaptation_count^(-.6)
      batch_acceptance <- mean(accepted[window])
      log_scale <- pmin(log(100), pmax(log(.01),
        log_scale + gain * (batch_acceptance - .234)))
      empirical <- cov(chain[seq_len(iteration), , drop = FALSE])
      if (any(!is.finite(empirical))) empirical <- proposal_base
      proposal_cov <- exp(log_scale) * (empirical + diag(1e-8, p))
      Rprop <- safe_chol(proposal_cov)
      eigenvalues <- eigen(proposal_cov, symmetric = TRUE,
                           only.values = TRUE)$values
      adaptation_history <- rbind(adaptation_history, data.frame(
        iteration = iteration, batch = adaptation_count,
        batch_acceptance = batch_acceptance, gain = gain,
        log_scale = log_scale, scale = exp(log_scale),
        min_eigenvalue = min(eigenvalues),
        max_eigenvalue = max(eigenvalues)
      ))
    }
  }
  keep <- seq.int(cfg$rwmh_burn + 1, cfg$rwmh_iter, by = cfg$rwmh_thin)
  list(
    sampler = "arwm",
    draws = chain[keep, , drop = FALSE],
    chain = chain,
    log_trace = trace,
    warmup_acceptance_rate = mean(accepted[2:cfg$rwmh_burn]),
    sampling_acceptance_rate = mean(accepted[(cfg$rwmh_burn + 1):
                                               cfg$rwmh_iter]),
    retained_iterations = keep,
    proposal_cov_final = proposal_cov,
    final_scale = exp(log_scale),
    adaptation_history = adaptation_history,
    runtime_seconds = unname(proc.time()[3] - start_time)
  )
}

# 6. Conditional GP and Full Bayesian mixture ------------------------------

conditional_source_posterior <- function(theta, problem) {
  Kll <- matern05_covariance(D2_ll, theta)
  Ktl <- matern05_covariance(D2_tl, theta)
  sigma_n <- exp(theta["log_sigma_n"])
  S <- H_inverse %*% Kll %*% t(H_inverse) +
    diag(sigma_n^2, length(problem$y))
  R <- safe_chol(S)
  cross <- Ktl %*% t(H_inverse)
  alpha <- chol_solve(R, problem$y)
  mean <- as.vector(cross %*% alpha)
  solved <- chol_solve(R, t(cross))
  prior_diag <- rep(exp(2 * theta["log_sigma_f"]), nrow(truth_grid))
  variance <- pmax(prior_diag - rowSums(cross * t(solved)), 0)
  list(mean = mean, variance = variance, sd = sqrt(variance))
}

representative_points <- data.frame(
  point = c("south_west", "centre", "south_east", "north_east"),
  x = c(.25, .50, .75, .75),
  y = c(.25, .50, .25, .75)
)
D2_rl <- pairwise_sqdist(representative_points, latent_grid)

conditional_representative_mean <- function(theta, problem) {
  Kll <- matern05_covariance(D2_ll, theta)
  Krl <- matern05_covariance(D2_rl, theta)
  S <- H_inverse %*% Kll %*% t(H_inverse) +
    diag(exp(2 * theta["log_sigma_n"]), length(problem$y))
  R <- safe_chol(S)
  alpha <- chol_solve(R, problem$y)
  as.vector(Krl %*% t(H_inverse) %*% alpha)
}

diagnose_scalar_chains <- function(chain_matrix, quantity, type,
                                   total_runtime) {
  pooled <- as.vector(chain_matrix)
  ess <- rank_normalized_ess(chain_matrix)
  data.frame(
    quantity = quantity,
    type = type,
    rank_split_Rhat = rank_normalized_rhat(chain_matrix),
    bulk_ESS = ess,
    MCSE_mean = sd(pooled) / sqrt(max(ess, 1)),
    posterior_mean = mean(pooled),
    posterior_sd = sd(pooled),
    ESS_per_second = ess / max(total_runtime, 1e-12)
  )
}

build_formal_diagnostics <- function(chain_list, problem) {
  total_runtime <- sum(vapply(chain_list, function(z) z$runtime_seconds,
                              numeric(1)))
  hyper_rows <- lapply(colnames(chain_list[[1]]$draws), function(parameter) {
    matrix_by_chain <- do.call(cbind, lapply(chain_list, function(z) {
      z$draws[, parameter]
    }))
    diagnose_scalar_chains(matrix_by_chain, parameter, "hyperparameter",
                           total_runtime)
  })

  functional_by_chain <- lapply(chain_list, function(z) {
    t(vapply(seq_len(nrow(z$draws)), function(i) {
      conditional_representative_mean(z$draws[i, ], problem)
    }, numeric(nrow(representative_points))))
  })
  functional_rows <- lapply(seq_len(nrow(representative_points)), function(j) {
    matrix_by_chain <- do.call(cbind, lapply(functional_by_chain,
                                             function(z) z[, j]))
    diagnose_scalar_chains(
      matrix_by_chain,
      paste0("mu_f(", representative_points$x[j], ",",
             representative_points$y[j], ")"),
      "downstream_posterior_mean", total_runtime
    )
  })
  list(
    detail = do.call(rbind, c(hyper_rows, functional_rows)),
    functional_by_chain = functional_by_chain,
    total_runtime = total_runtime
  )
}

# Exact pointwise log predictive density of a finite Gaussian mixture.
mixture_nlpd <- function(component_mean, component_sd, truth) {
  S <- nrow(component_mean)
  density_log <- vapply(seq_along(truth), function(j) {
    sdj <- pmax(component_sd[, j], 1e-10)
    log_sum_exp(dnorm(truth[j], component_mean[, j], sdj, log = TRUE)) - log(S)
  }, numeric(1))
  -mean(density_log)
}

# Empirical CRPS from pointwise posterior draws.  For sorted x_(i),
# (1/(2S^2))*sum_ij |x_i-x_j| = sum_i (2i-S-1)x_(i)/S^2.
empirical_crps <- function(draws, truth) {
  S <- nrow(draws)
  first <- colMeans(abs(sweep(draws, 2, truth, "-")))
  ordered <- apply(draws, 2, sort)
  coeff <- 2 * seq_len(S) - S - 1
  second <- as.vector(crossprod(coeff, ordered)) / S^2
  mean(first - second)
}

integrate_hyperparameters <- function(fit, problem, seed) {
  draws <- fit$draws
  support_index <- unique(round(seq(1, nrow(draws), length.out =
                                    min(cfg$mixture_support, nrow(draws)))))
  support <- draws[support_index, , drop = FALSE]
  component_mean <- matrix(NA_real_, nrow(support), nrow(truth_grid))
  component_var <- matrix(NA_real_, nrow(support), nrow(truth_grid))
  for (s in seq_len(nrow(support))) {
    conditional <- conditional_source_posterior(support[s, ], problem)
    component_mean[s, ] <- conditional$mean
    component_var[s, ] <- conditional$variance
  }
  posterior_mean <- colMeans(component_mean)
  posterior_variance <- colMeans(component_var + component_mean^2) -
    posterior_mean^2
  posterior_sd <- sqrt(pmax(posterior_variance, 0))

  # Draw pointwise marginal samples from the finite mixture.  They are used for
  # quantiles, coverage and CRPS; the reported mean/variance above are analytic.
  set.seed(seed)
  component_choice <- sample(seq_len(nrow(support)), cfg$mixture_mc_draws,
                             replace = TRUE)
  marginal_draws <- matrix(NA_real_, cfg$mixture_mc_draws, nrow(truth_grid))
  for (m in seq_len(cfg$mixture_mc_draws)) {
    s <- component_choice[m]
    marginal_draws[m, ] <- rnorm(nrow(truth_grid), component_mean[s, ],
                                 sqrt(pmax(component_var[s, ], 0)))
  }
  lower95 <- apply(marginal_draws, 2, quantile, probs = .025, names = FALSE)
  upper95 <- apply(marginal_draws, 2, quantile, probs = .975, names = FALSE)

  calibration <- do.call(rbind, lapply(cfg$credibility, function(level) {
    tail <- (1 - level) / 2
    lower <- apply(marginal_draws, 2, quantile, probs = tail, names = FALSE)
    upper <- apply(marginal_draws, 2, quantile, probs = 1 - tail,
                   names = FALSE)
    data.frame(nominal = level,
               empirical = mean(problem$truth >= lower & problem$truth <= upper))
  }))

  metrics <- data.frame(
    RMSE = rmse(posterior_mean, problem$truth),
    MAE = mae(posterior_mean, problem$truth),
    Coverage95 = mean(problem$truth >= lower95 & problem$truth <= upper95),
    Width95 = mean(upper95 - lower95),
    Calibration_MAE = mean(abs(calibration$empirical - calibration$nominal)),
    CRPS = empirical_crps(marginal_draws, problem$truth),
    NLPD = mixture_nlpd(component_mean, sqrt(component_var), problem$truth)
  )
  list(mean = posterior_mean, sd = posterior_sd,
       lower95 = lower95, upper95 = upper95,
       component_mean = component_mean, component_var = component_var,
       marginal_draws = marginal_draws, calibration = calibration,
       metrics = metrics)
}

stability_from_components <- function(posterior, problem) {
  available <- nrow(posterior$component_mean)
  sizes <- cfg$stability_sizes[cfg$stability_sizes <= available]
  if (!length(sizes)) stop("No requested stability size is available.")
  fields <- lapply(sizes, function(S) {
    mu <- posterior$component_mean[seq_len(S), , drop = FALSE]
    vv <- posterior$component_var[seq_len(S), , drop = FALSE]
    mean_field <- colMeans(mu)
    within <- colMeans(vv)
    between <- pmax(colMeans(mu^2) - mean_field^2, 0)
    total <- within + between
    list(S = S, mean = mean_field, within = within, between = between,
         total = total, lower = mean_field - 1.96 * sqrt(total),
         upper = mean_field + 1.96 * sqrt(total))
  })
  reference <- fields[[length(fields)]]
  table <- do.call(rbind, lapply(fields, function(z) data.frame(
    retained_draws = z$S,
    RMSE = rmse(z$mean, problem$truth),
    mean_within_variance = mean(z$within),
    mean_between_variance = mean(z$between),
    mean_total_variance = mean(z$total),
    mean_interval_width = mean(z$upper - z$lower),
    mean_abs_change_from_Smax = mean(abs(z$mean - reference$mean)),
    mean_abs_sd_change_from_Smax =
      mean(abs(sqrt(z$total) - sqrt(reference$total)))
  )))
  list(table = table, fields = fields)
}

# 7. Visualisation ---------------------------------------------------------

field_plot <- function(values, title, option = "D", limits = NULL) {
  dat <- data.frame(truth_grid, value = values)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_fill_viridis_c(option = option, limits = limits,
                                  name = "value") +
    ggplot2::labs(title = title, x = "x", y = "y") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5,
                                                       face = "bold",
                                                       size = 10),
                   panel.grid = ggplot2::element_blank())
  p
}

save_reconstruction_plot <- function(problem, posterior) {
  common <- range(c(problem$truth, posterior$mean))
  plots <- list(
    field_plot(problem$truth, "True source", "D", common),
    field_plot(posterior$mean,
               sprintf("Posterior mean (RMSE=%.4f)",
                       posterior$metrics$RMSE), "D", common),
    field_plot(abs(posterior$mean - problem$truth),
               sprintf("Absolute error (MAE=%.4f)", posterior$metrics$MAE), "A"),
    field_plot(posterior$sd,
               sprintf("Posterior SD (95%% coverage=%.3f)",
                       posterior$metrics$Coverage95), "B")
  )
  stem <- file.path(figure_dir, paste0(problem$key, "_arwm"))
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    grDevices::png(paste0(stem, "_reconstruction.png"), width = 1900,
                   height = 1500, res = 190)
    gridExtra::grid.arrange(grobs = plots, ncol = 2,
                            top = paste(problem$label, "-",
                                        "four-chain ARWM Full Bayes"))
    grDevices::dev.off()
  } else {
    for (j in seq_along(plots)) {
      ggplot2::ggsave(paste0(stem, "_panel", j, ".png"), plots[[j]],
                      width = 5.4, height = 4.4, dpi = 180)
    }
  }
}

save_trace_plot <- function(chain_list, problem) {
  dat <- do.call(rbind, lapply(seq_along(chain_list), function(chain_id) {
    one <- as.data.frame(chain_list[[chain_id]]$chain)
    one$iteration <- seq_len(nrow(one))
    one$chain <- factor(chain_id)
    one
  }))
  long <- reshape(dat,
                  varying = c("log_ell", "log_sigma_f", "log_sigma_n"),
                  v.names = "value", timevar = "parameter",
                  times = c("log ell", "log sigma_f", "log sigma_n"),
                  idvar = c("chain", "iteration"), direction = "long")
  p <- ggplot2::ggplot(long, ggplot2::aes(iteration, value)) +
    ggplot2::geom_line(ggplot2::aes(colour = chain), linewidth = .25,
                       alpha = .85) +
    ggplot2::geom_vline(xintercept = cfg$rwmh_burn,
                        linetype = 2, colour = "black") +
    ggplot2::facet_grid(parameter ~ chain, scales = "free_y") +
    ggplot2::scale_colour_viridis_d(option = "D", guide = "none") +
    ggplot2::labs(title = paste(problem$label, "- four ARWM chains"),
                  subtitle = "Dashed line marks the end of warm-up",
                  x = "iteration", y = NULL) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir,
                            paste0(problem$key, "_arwm_trace.png")),
                  p, width = 11.2, height = 7.0, dpi = 180)
}

save_stability_plot <- function(stability, problem) {
  table <- stability$table
  long <- rbind(
    data.frame(S = table$retained_draws, quantity = "RMSE",
               value = table$RMSE),
    data.frame(S = table$retained_draws, quantity = "Mean within variance",
               value = table$mean_within_variance),
    data.frame(S = table$retained_draws, quantity = "Mean between variance",
               value = table$mean_between_variance),
    data.frame(S = table$retained_draws, quantity = "Mean 95% band width",
               value = table$mean_interval_width)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(S, value)) +
    ggplot2::geom_line(colour = "#2C7FB8", linewidth = .8) +
    ggplot2::geom_point(colour = "#D95F0E", size = 2.1) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = 2) +
    ggplot2::scale_x_continuous(breaks = cfg$stability_sizes) +
    ggplot2::labs(title = paste(problem$label, "- Monte Carlo stability"),
                  x = "Pooled retained draws S", y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(file.path(figure_dir,
                            paste0(problem$key, "_stability.png")),
                  p, width = 8.5, height = 6.0, dpi = 190)
}

save_adaptation_plot <- function(chain_list, problem) {
  dat <- do.call(rbind, lapply(seq_along(chain_list), function(chain_id) {
    one <- chain_list[[chain_id]]$adaptation_history
    one$chain <- factor(chain_id)
    one
  }))
  long <- rbind(
    data.frame(iteration = dat$iteration, chain = dat$chain,
               quantity = "Batch acceptance", value = dat$batch_acceptance),
    data.frame(iteration = dat$iteration, chain = dat$chain,
               quantity = "Global scale s", value = dat$scale),
    data.frame(iteration = dat$iteration, chain = dat$chain,
               quantity = "Largest proposal eigenvalue",
               value = dat$max_eigenvalue),
    data.frame(iteration = dat$iteration, chain = dat$chain,
               quantity = "Smallest proposal eigenvalue",
               value = dat$min_eigenvalue)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(iteration, value, colour = chain)) +
    ggplot2::geom_line(linewidth = .55) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = 2) +
    ggplot2::scale_colour_viridis_d(option = "D", name = "Chain") +
    ggplot2::labs(
      title = paste(problem$label, "- warm-up adaptation"),
      subtitle = "Empirical-covariance shape and Robbins-Monro global scale; both freeze after warm-up",
      x = "Warm-up iteration", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir,
                            paste0(problem$key, "_arwm_adaptation.png")),
                  p, width = 10.5, height = 6.5, dpi = 190)
}

# 8. Complete source-specific four-chain experiment ------------------------

cache_path <- file.path(output_dir, "complete_mcmc_experiment.rds")
if (file.exists(cache_path) && !force_refit) {
  message("Loading existing result: ", cache_path)
  experiment <- readRDS(cache_path)
} else {
  problem <- problems[[1]]
  info <- parameter_info(problem)
  target <- log_target_factory(problem, info)
  message("MAP initialisation: ", problem$key)
  map <- map_initialise(problem, info, cfg$seed + 10007)

  offsets <- rbind(
    c(0, 0, 0),
    c(-.35, .25, -.15),
    c(.35, -.25, .15),
    c(-.20, -.30, .25)
  )
  colnames(offsets) <- names(map$theta)
  chain_list <- vector("list", cfg$n_chains)
  for (chain_id in seq_len(cfg$n_chains)) {
    init <- map$theta + offsets[chain_id, ]
    init <- pmin(pmax(init, info$lower + 1e-8), info$upper - 1e-8)
    names(init) <- names(map$theta)
    message("Running ARWM chain ", chain_id, "/", cfg$n_chains,
            " for ", problem$key, "...")
    chain_list[[chain_id]] <- run_rwmh(
      target, init, info,
      cfg$seed + 10007 + chain_id * 503
    )
    chain_list[[chain_id]]$chain_id <- chain_id
    chain_list[[chain_id]]$initial_theta <- init
  }

  formal <- build_formal_diagnostics(chain_list, problem)
  formal$detail$source <- problem$key
  pooled_draws <- balanced_pool(chain_list, max(cfg$stability_sizes))
  combined_fit <- list(draws = pooled_draws, sampler = "arwm")
  posterior <- integrate_hyperparameters(
    combined_fit, problem, cfg$seed + 50997
  )
  stability <- stability_from_components(posterior, problem)
  stability$table$source <- problem$key

  metrics <- cbind(
    data.frame(source = problem$key, source_label = problem$label,
               sampler = "four-chain ARWM", kernel = "Matern 0.5",
               noise_ratio = cfg$noise_ratio,
               sigma_noise = problem$sigma_noise),
    posterior$metrics
  )
  calibration <- transform(posterior$calibration,
                           source = problem$key, sampler = "ARWM")
  chain_diagnostics <- do.call(rbind, lapply(chain_list, function(z) {
    final_cov <- z$proposal_cov_final
    data.frame(
      source = problem$key,
      chain = z$chain_id,
      retained_draws = nrow(z$draws),
      warmup_acceptance_rate = z$warmup_acceptance_rate,
      sampling_acceptance_rate = z$sampling_acceptance_rate,
      final_scale = z$final_scale,
      final_cov_11 = final_cov[1, 1],
      final_cov_22 = final_cov[2, 2],
      final_cov_33 = final_cov[3, 3],
      final_cov_12 = final_cov[1, 2],
      final_cov_13 = final_cov[1, 3],
      final_cov_23 = final_cov[2, 3],
      runtime_seconds = z$runtime_seconds,
      initial_log_ell = z$initial_theta["log_ell"],
      initial_log_sigma_f = z$initial_theta["log_sigma_f"],
      initial_log_sigma_n = z$initial_theta["log_sigma_n"]
    )
  }))
  hyper_diag <- formal$detail[formal$detail$type == "hyperparameter", ]
  functional_diag <- formal$detail[
    formal$detail$type == "downstream_posterior_mean", ]
  diagnostic_summary <- data.frame(
    source = problem$key,
    chains = cfg$n_chains,
    pooled_retained_draws = nrow(pooled_draws),
    max_Rhat_hyperparameters = max(hyper_diag$rank_split_Rhat),
    min_ESS_hyperparameters = min(hyper_diag$bulk_ESS),
    max_Rhat_downstream = max(functional_diag$rank_split_Rhat),
    min_ESS_downstream = min(functional_diag$bulk_ESS),
    mean_sampling_acceptance = mean(chain_diagnostics$sampling_acceptance_rate),
    min_ESS_per_second = min(formal$detail$ESS_per_second),
    total_sampling_seconds = formal$total_runtime
  )
  draw_frame <- do.call(rbind, lapply(chain_list, function(z) {
    one <- as.data.frame(z$draws)
    one$source <- problem$key
    one$chain <- z$chain_id
    one$retained_draw <- seq_len(nrow(one))
    one
  }))
  adaptation_frame <- do.call(rbind, lapply(chain_list, function(z) {
    one <- z$adaptation_history
    one$source <- problem$key
    one$chain <- z$chain_id
    one
  }))

  experiment <- list(
    cfg = cfg, obs_grid = obs_grid, latent_grid = latent_grid,
    truth_grid = truth_grid, H_inverse = H_inverse, H_truth = H_truth,
    problem = problem, map = map, chains = chain_list,
    posterior = posterior, stability = stability,
    metrics = metrics, calibration = calibration,
    diagnostics = formal$detail,
    diagnostic_summary = diagnostic_summary,
    chain_diagnostics = chain_diagnostics,
    adaptation_history = adaptation_frame,
    hyperparameter_draws = draw_frame
  )
  if (!no_plots) {
    save_reconstruction_plot(problem, posterior)
    save_trace_plot(chain_list, problem)
    save_stability_plot(stability, problem)
    save_adaptation_plot(chain_list, problem)
  }
  saveRDS(experiment, cache_path)
}

write.csv(experiment$metrics, file.path(output_dir, "mcmc_metrics.csv"),
          row.names = FALSE)
write.csv(experiment$calibration, file.path(output_dir, "mcmc_calibration.csv"),
          row.names = FALSE)
write.csv(experiment$diagnostics, file.path(output_dir, "mcmc_diagnostics.csv"),
          row.names = FALSE)
write.csv(experiment$diagnostic_summary,
          file.path(output_dir, "mcmc_diagnostic_summary.csv"), row.names = FALSE)
write.csv(experiment$chain_diagnostics,
          file.path(output_dir, "mcmc_chain_diagnostics.csv"), row.names = FALSE)
write.csv(experiment$adaptation_history,
          file.path(output_dir, "mcmc_adaptation_history.csv"), row.names = FALSE)
write.csv(experiment$stability$table,
          file.path(output_dir, "mcmc_stability.csv"), row.names = FALSE)
write.csv(experiment$hyperparameter_draws,
          file.path(output_dir, "hyperparameter_draws.csv"), row.names = FALSE)
write.csv(data.frame(
  parameter = c("seed", "source", "kernel", "noise_ratio", "n_inverse",
                "n_truth", "n_modes", "n_obs_side", "algorithm", "chains",
                "iterations", "warmup", "thin"),
  value = c(cfg$seed, SOURCE_KEY, "Matern 0.5", cfg$noise_ratio,
            cfg$n_inverse, cfg$n_truth, cfg$n_modes, cfg$n_obs_side,
            "Dual-adaptation Random-Walk Metropolis with frozen sampling proposal",
            cfg$n_chains,
            cfg$rwmh_iter, cfg$rwmh_burn, cfg$rwmh_thin)
), file.path(output_dir, "configuration.csv"), row.names = FALSE)

if (!no_plots) {
  p_cal <- ggplot2::ggplot(
    experiment$calibration,
    ggplot2::aes(nominal, empirical, group = 1)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
                         colour = "grey35") +
    ggplot2::geom_line(linewidth = .75, colour = "#2C7FB8") +
    ggplot2::geom_point(size = 1.8, colour = "#D95F0E") +
    ggplot2::coord_equal(xlim = c(.45, 1), ylim = c(.45, 1)) +
    ggplot2::labs(title = "Pointwise credible-interval calibration",
                  subtitle = "Diagonal: empirical coverage equals nominal credibility",
                  x = "Nominal credibility", y = "Empirical coverage") +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(file.path(figure_dir, "calibration.png"), p_cal,
                  width = 10.5, height = 7.2, dpi = 190)
}

capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
message("Finished. Results are in: ", output_dir)
print(experiment$metrics)
