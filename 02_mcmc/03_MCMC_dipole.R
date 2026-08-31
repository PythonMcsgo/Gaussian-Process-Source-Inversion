# Fully Bayesian inversion of the positive-negative dipole source.
SOURCE_KEY <- "dipole"
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", full_args[grepl("^--file=", full_args)])
wrapper_file <- if (length(file_arg)) file_arg[1] else
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
wrapper_dir <- if (is.null(wrapper_file) || !length(wrapper_file)) getwd() else
  dirname(normalizePath(wrapper_file, winslash = "/", mustWork = FALSE))
source(file.path(wrapper_dir, "00_MCMC_core.R"), chdir = TRUE)
