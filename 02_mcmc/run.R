#!/usr/bin/env Rscript

cli <- commandArgs(trailingOnly = TRUE)
script_location <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  command_files <- sub("^--file=", "", full[grepl("^--file=", full)])
  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  frame_dirs <- if (length(frame_files)) dirname(unlist(frame_files, use.names=FALSE)) else character()
  seed_dirs <- c(dirname(command_files), frame_dirs, getwd())
  candidates <- unique(c(seed_dirs, file.path(seed_dirs, "02_mcmc"),
                         file.path(dirname(seed_dirs), "02_mcmc")))
  candidates <- candidates[dir.exists(candidates)]
  hits <- candidates[vapply(candidates, function(path) {
    all(file.exists(file.path(path, c("00_MCMC_core.R",
                                     "07_assemble_report.R"))))
  }, logical(1))]
  if (!length(hits)) stop("Cannot locate the 02_mcmc directory.")
  normalizePath(hits[1], winslash="/", mustWork=TRUE)
}

here <- script_location()
project_root <- dirname(here)
results_root <- file.path(project_root, "results")
rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")

source_keys <- c("spike", "mixture", "dipole", "ring", "edge", "sine")
jobs <- c(
  "01_MCMC_spike.R", "02_MCMC_broad_plus_sharp.R", "03_MCMC_dipole.R",
  "04_MCMC_ring.R", "05_MCMC_moving_edge.R", "06_MCMC_sine.R"
)
prefix <- "MCMC_results_"
forward <- cli

for (i in seq_along(jobs)) {
  output <- file.path(results_root, paste0(prefix, source_keys[i]))
  args <- c(shQuote(file.path(here, jobs[i])), forward,
            paste0("--output=", shQuote(output)))
  message("Running ", jobs[i])
  status <- system2(rscript, args)
  if (status != 0) stop("MCMC job failed: ", jobs[i])
}

assembler_args <- c(
  shQuote(file.path(here, "07_assemble_report.R")),
  paste0("--root=", shQuote(results_root)),
  paste0("--prefix=", prefix),
  paste0("--output=", shQuote(file.path(results_root, "MCMC_report")))
)
status <- system2(rscript, assembler_args)
if (status != 0) stop("MCMC report assembly failed.")
message("02 MCMC completed.")
