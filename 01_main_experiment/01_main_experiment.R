
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install ggplot2 first: install.packages('ggplot2')")
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("Install gridExtra first: install.packages('gridExtra')")
})

MASTER_SEED <- 20161387
N_INVERSE <- 15
N_TRUTH <- 41
N_MODES <- 36
N_OBS_SIDE <- 7
NOISE_RATIO <- 1e-2
N_STARTS <- 6
OPTIM_MAXIT <- 180
CREDIBILITY <- c(.50, .60, .70, .80, .90, .95, .975, .99)

script_path <- function() {
  z <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", z[grepl("^--file=", z)])
  if (length(hit)) dirname(normalizePath(hit[1], winslash = "/")) else getwd()
}
ROOT <- normalizePath(file.path(script_path(), ".."), winslash = "/",
                      mustWork = FALSE)
OUT <- file.path(ROOT, "results", "main_experiment")
FIG <- file.path(OUT, "figures")
CHAPTER_FIG <- FIG
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(CHAPTER_FIG, recursive = TRUE, showWarnings = FALSE)

safe_chol <- function(A, initial = 1e-12, tries = 10) {
  scale <- max(mean(diag(A)), 1)
  for (i in 0:tries) {
    jitter <- if (i == 0) 0 else initial * scale * 10^(i - 1)
    R <- tryCatch(chol(A + diag(jitter, nrow(A))), error = function(e) NULL)
    if (!is.null(R)) return(R)
  }
  stop("Cholesky factorisation failed after adaptive jitter.")
}
chol_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))
sqdist <- function(a, b = a) {
  a <- as.matrix(a); b <- as.matrix(b)
  a2 <- rowSums(a^2); b2 <- rowSums(b^2)
  pmax(outer(a2, b2, "+") - 2 * tcrossprod(a, b), 0)
}

make_grid <- function(n, cell_centred = FALSE) {
  axis <- if (cell_centred) (seq_len(n) - .5) / n else seq(0, 1, length.out = n)
  expand.grid(x = axis, y = axis)
}
inverse_grid <- make_grid(N_INVERSE, TRUE)
truth_grid <- make_grid(N_TRUTH, TRUE)
obs_axis <- seq(.1, .9, length.out = N_OBS_SIDE)
obs_grid <- expand.grid(x = obs_axis, y = obs_axis)
inverse_weight <- 1 / N_INVERSE^2
truth_weight <- 1 / N_TRUTH^2

green_matrix <- function(obs, src, weight, n_modes = N_MODES) {
  p <- rep(seq_len(n_modes), each = n_modes)
  q <- rep(seq_len(n_modes), times = n_modes)
  lambda <- pi^2 * (p^2 + q^2)
  phi <- function(z) 2 * sin(pi * outer(z$x, p)) * sin(pi * outer(z$y, q))
  (phi(obs) %*% (t(phi(src)) / lambda)) * weight
}
H_inverse <- green_matrix(obs_grid, inverse_grid, inverse_weight)
H_truth <- green_matrix(obs_grid, truth_grid, truth_weight)

source_functions <- list(
  spike = function(x, y) exp(-((x-.73)^2 + (y-.27)^2)/(2*.035^2)),
  broad_plus_sharp = function(x, y) {
    .7 * exp(-((x-.35)^2 + (y-.65)^2)/(2*.18^2)) +
      1.2 * exp(-((x-.75)^2 + (y-.25)^2)/(2*.035^2))
  },
  dipole = function(x, y) {
    exp(-((x-.35)^2 + (y-.50)^2)/(2*.07^2)) -
      exp(-((x-.65)^2 + (y-.50)^2)/(2*.07^2))
  },
  ring = function(x, y) {
    r <- sqrt((x-.5)^2 + (y-.5)^2)
    exp(-((r-.22)^2)/(2*.035^2))
  },
  moving_edge = function(x, y) as.numeric(x > .45 + .15*sin(4*pi*y)),
  sine = function(x, y) sin(pi*x) * cos(2*pi*y)
)
source_labels <- c(
  spike = "Hidden narrow spike", broad_plus_sharp = "Broad plus sharp",
  dipole = "Positive-negative dipole", ring = "Ring-shaped source",
  moving_edge = "Moving-edge source", sine = "Sine-cosine source"
)

kernel_correlation <- function(name, d2, ell) {
  r <- sqrt(d2)
  switch(name,
    rbf = exp(-d2/(2*ell^2)),
    matern05 = exp(-r/ell),
    matern15 = (1 + sqrt(3)*r/ell) * exp(-sqrt(3)*r/ell),
    matern25 = (1 + sqrt(5)*r/ell + 5*d2/(3*ell^2)) * exp(-sqrt(5)*r/ell),
    stop("Unknown kernel: ", name)
  )
}
kernel_labels <- c(rbf="RBF", matern05="Matern 0.5",
                   matern15="Matern 1.5", matern25="Matern 2.5")
D2_ll <- sqdist(inverse_grid)
D2_tl <- sqdist(truth_grid, inverse_grid)

covariance_blocks <- function(kernel, theta) {
  ell <- exp(theta[1]); sf <- exp(theta[2]); sn <- exp(theta[3])
  Kll <- sf^2 * kernel_correlation(kernel, D2_ll, ell)
  Ktl <- sf^2 * kernel_correlation(kernel, D2_tl, ell)
  S <- H_inverse %*% Kll %*% t(H_inverse) + diag(sn^2, nrow(H_inverse))
  list(Kll=Kll, Ktl=Ktl, S=S, sf=sf, sn=sn)
}

log_marginal <- function(theta, kernel, y) {
  b <- covariance_blocks(kernel, theta)
  R <- tryCatch(safe_chol(b$S), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  -.5 * drop(crossprod(y, chol_solve(R, y))) - sum(log(diag(R))) -
    .5 * length(y) * log(2*pi)
}
log_prior <- function(theta, sigma_noise) {
  sum(dnorm(theta,
    mean=c(log(.18), log(.60), log(sigma_noise)),
    sd=c(.90, .70, .35), log=TRUE))
}

fit_hyperparameters <- function(kernel, y, sigma_noise, method) {
  lower <- c(log(.015), log(.015), log(max(sigma_noise*.01, 1e-10)))
  upper <- c(log(1.50), log(5.00), log(max(sigma_noise*20, sd(y), 1e-7)))
  centre <- c(log(.18), log(.60), log(sigma_noise))
  set.seed(MASTER_SEED + match(kernel, names(kernel_labels))*103 +
             match(method, c("multistart_eb", "map"))*1009)
  starts <- rbind(centre, matrix(runif((N_STARTS-1)*3, lower, upper),
                                 ncol=3, byrow=TRUE))
  objective <- function(theta) {
    value <- log_marginal(theta, kernel, y)
    if (method == "map") value <- value + log_prior(theta, sigma_noise)
    if (!is.finite(value)) 1e100 else -value
  }
  fits <- lapply(seq_len(N_STARTS), function(i) optim(
    starts[i,], objective, method="L-BFGS-B", lower=lower, upper=upper,
    control=list(maxit=OPTIM_MAXIT)))
  values <- vapply(fits, `[[`, numeric(1), "value")
  best <- fits[[which.min(values)]]
  names(best$par) <- c("log_ell", "log_sigma_f", "log_sigma_n")
  list(theta=best$par, objective=best$value, convergence=best$convergence,
       starts=data.frame(start=seq_len(N_STARTS), objective=values))
}

gp_posterior <- function(kernel, theta, y) {
  b <- covariance_blocks(kernel, theta)
  R <- safe_chol(b$S)
  A <- b$Kll %*% t(H_inverse)
  mean_latent <- drop(A %*% chol_solve(R, y))
  mean_truth <- drop(b$Ktl %*% t(H_inverse) %*% chol_solve(R, y))
  cross <- b$Ktl %*% t(H_inverse)
  root <- forwardsolve(t(R), t(cross))
  variance_truth <- pmax(b$sf^2 - colSums(root^2), 0)
  list(mean=mean_truth, sd=sqrt(variance_truth), mean_latent=mean_latent)
}

tikhonov_discrepancy <- function(H, y, sigma_noise, tau=1.05) {
  s <- svd(H); gamma <- s$d^2; yq <- crossprod(s$u, y)
  target <- tau*sqrt(length(y))*sigma_noise
  residual <- function(log_lambda) {
    lambda <- exp(log_lambda)
    sqrt(sum((lambda/(gamma+lambda)*yq)^2))
  }
  grid <- seq(log(max(gamma)*1e-12), log(max(gamma)*1e8), length.out=220)
  values <- vapply(grid, residual, numeric(1)) - target
  hit <- which(values[-length(values)]*values[-1] <= 0)[1]
  if (is.na(hit)) {
    log_lambda <- grid[which.min(abs(values))]
    status <- "nearest"
  } else {
    log_lambda <- uniroot(function(z) residual(z)-target,
                          interval=c(grid[hit], grid[hit+1]),
                          tol=1e-10)$root
    status <- "root"
  }
  lambda <- exp(log_lambda)
  coef <- s$d/(gamma+lambda) * yq
  estimate <- drop(s$v %*% coef)
  list(mean_latent=estimate, lambda=lambda, residual=residual(log_lambda),
       target=target, ratio=residual(log_lambda)/(sqrt(length(y))*sigma_noise),
       status=status)
}

interpolate_regular_grid <- function(values, from_n=N_INVERSE, to_n=N_TRUTH) {
  Z <- matrix(values, nrow=from_n, ncol=from_n)
  x0 <- (seq_len(from_n)-.5)/from_n; x1 <- seq(0,1,length.out=to_n)
  along_x <- vapply(seq_len(from_n), function(j)
    approx(x0, Z[,j], xout=x1, rule=2)$y, numeric(to_n))
  along_y <- t(vapply(seq_len(to_n), function(i)
    approx(x0, along_x[i,], xout=x1, rule=2)$y, numeric(to_n)))
  as.vector(along_y)
}

gaussian_scores <- function(mean, sd, truth) {
  sd <- pmax(sd, 1e-10); z <- (truth-mean)/sd
  calibration <- data.frame(nominal=CREDIBILITY, empirical=vapply(
    CREDIBILITY, function(a) {
      q <- qnorm((1+a)/2); mean(abs(truth-mean) <= q*sd)
    }, numeric(1)))
  crps <- sd * (z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi))
  metrics <- c(RMSE=sqrt(mean((mean-truth)^2)), MAE=mean(abs(mean-truth)),
    Coverage95=mean(abs(truth-mean) <= qnorm(.975)*sd),
    CalibrationMAE=mean(abs(calibration$empirical-calibration$nominal)),
    CRPS=mean(crps), NLPD=mean(.5*log(2*pi*sd^2)+.5*z^2))
  list(metrics=metrics, calibration=calibration)
}

field_plot <- function(grid, value, title, limits=NULL, palette="viridis") {
  d <- data.frame(grid, value=value)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, fill=value)) +
    ggplot2::geom_raster() + ggplot2::coord_equal(expand=FALSE) +
    ggplot2::theme_minimal(base_size=9) +
    ggplot2::theme(panel.grid=ggplot2::element_blank(),
                   plot.title=ggplot2::element_text(hjust=.5, size=9),
                   legend.position="right") + ggplot2::labs(title=title, fill="value")
  if (palette == "magma") p + ggplot2::scale_fill_viridis_c(option="magma", limits=limits)
  else p + ggplot2::scale_fill_viridis_c(option="viridis", limits=limits)
}

all_metrics <- list(); all_calibration <- list(); all_hyper <- list()
all_tikh <- list(); record <- 1L
for (source_name in names(source_functions)) {
  set.seed(MASTER_SEED + match(source_name, names(source_functions))*1009L)
  f_truth <- with(truth_grid, source_functions[[source_name]](x,y))
  u_clean <- drop(H_truth %*% f_truth)
  sigma_noise <- NOISE_RATIO * max(abs(u_clean))
  y <- u_clean + rnorm(length(u_clean), sd=sigma_noise)
  tik <- tikhonov_discrepancy(H_inverse, y, sigma_noise)
  tik_truth <- interpolate_regular_grid(tik$mean_latent)
  all_tikh[[source_name]] <- data.frame(source=source_name, lambda=tik$lambda,
    residual=tik$residual, target=tik$target, ratio=tik$ratio, status=tik$status,
    RMSE=sqrt(mean((tik_truth-f_truth)^2)), MAE=mean(abs(tik_truth-f_truth)))

  for (method in c("multistart_eb", "map")) {
    fits <- list()
    for (kernel in names(kernel_labels)) {
      fit <- fit_hyperparameters(kernel, y, sigma_noise, method)
      post <- gp_posterior(kernel, fit$theta, y)
      score <- gaussian_scores(post$mean, post$sd, f_truth)
      fits[[kernel]] <- list(fit=fit, post=post, score=score)
      all_metrics[[record]] <- data.frame(source=source_name,
        source_label=source_labels[source_name], hyper_method=method,
        kernel=kernel, kernel_label=kernel_labels[kernel],
        as.list(score$metrics), row.names=NULL)
      all_calibration[[record]] <- cbind(data.frame(source=source_name,
        hyper_method=method, kernel=kernel), score$calibration)
      all_hyper[[record]] <- data.frame(source=source_name, hyper_method=method,
        kernel=kernel, parameter=names(fit$theta), log_value=fit$theta,
        natural_value=exp(fit$theta), convergence=fit$convergence)
      record <- record + 1L
    }
    if (method == "multistart_eb") {
      field_limit <- range(f_truth, tik_truth,
        unlist(lapply(fits, function(z) z$post$mean)), finite=TRUE)
      mean_panels <- c(list(field_plot(truth_grid, f_truth, "True source", field_limit),
        field_plot(truth_grid, tik_truth,
          sprintf("Tikhonov mean; RMSE %.3f", all_tikh[[source_name]]$RMSE), field_limit)),
        lapply(names(fits), function(k) field_plot(truth_grid, fits[[k]]$post$mean,
          sprintf("%s mean; RMSE %.3f", kernel_labels[k], fits[[k]]$score$metrics["RMSE"]), field_limit)))
      error_limit <- c(0, max(abs(tik_truth-f_truth),
        unlist(lapply(fits, function(z) abs(z$post$mean-f_truth)))))
      error_panels <- c(list(field_plot(truth_grid, abs(tik_truth-f_truth),
        "Tikhonov absolute error", error_limit, "magma")),
        lapply(names(fits), function(k) field_plot(truth_grid,
          abs(fits[[k]]$post$mean-f_truth), paste(kernel_labels[k], "absolute error"),
          error_limit, "magma")))
      sd_limit <- c(0, max(unlist(lapply(fits, function(z) z$post$sd))))
      sd_panels <- lapply(names(fits), function(k) field_plot(truth_grid,
        fits[[k]]$post$sd, paste(kernel_labels[k], "posterior SD"), sd_limit, "magma"))
      # The three rows contain different numbers of scientifically meaningful
      # panels.  Nesting separate row layouts avoids blank placeholder cells.
      posterior_layout <- gridExtra::arrangeGrob(
        gridExtra::arrangeGrob(grobs=mean_panels,
                               ncol=length(mean_panels)),
        gridExtra::arrangeGrob(grobs=error_panels,
                               ncol=length(error_panels)),
        gridExtra::arrangeGrob(grobs=sd_panels,
                               ncol=length(sd_panels)),
        ncol=1, heights=c(1,1,1))
      grDevices::png(file.path(FIG, paste0(source_name,"_posterior_fields.png")),
                     3000, 1750, res=220)
      grid::grid.draw(posterior_layout)
      grDevices::dev.off()

      if (source_name == "moving_edge") {
        edge_layout <- gridExtra::arrangeGrob(
          gridExtra::arrangeGrob(grobs=error_panels,
                                 ncol=length(error_panels)),
          gridExtra::arrangeGrob(grobs=sd_panels,
                                 ncol=length(sd_panels)),
          ncol=1, heights=c(1,1))
        grDevices::png(file.path(CHAPTER_FIG, "diagnostics_edge.png"),
                       3000, 1250, res=220)
        grid::grid.draw(edge_layout)
        grDevices::dev.off()
      }
    }
  }
}

metrics <- do.call(rbind, all_metrics)
calibration <- do.call(rbind, all_calibration)
hyperparameters <- do.call(rbind, all_hyper)
tikhonov <- do.call(rbind, all_tikh)
write.csv(metrics, file.path(OUT,"gp_metrics.csv"), row.names=FALSE)
write.csv(calibration, file.path(OUT,"gp_calibration.csv"), row.names=FALSE)
write.csv(hyperparameters, file.path(OUT,"hyperparameters.csv"), row.names=FALSE)
write.csv(tikhonov, file.path(OUT,"tikhonov_metrics.csv"), row.names=FALSE)
write.csv(data.frame(parameter=c("master_seed","inverse_grid","truth_grid",
  "green_modes_per_axis","observation_grid","noise_ratio","optim_starts"),
  value=c(MASTER_SEED,N_INVERSE,N_TRUTH,N_MODES,N_OBS_SIDE,NOISE_RATIO,N_STARTS)),
  file.path(OUT,"configuration.csv"), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(OUT,"session_info.txt"))
message("Main experiment complete: ", OUT)
