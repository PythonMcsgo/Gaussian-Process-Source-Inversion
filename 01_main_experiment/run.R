script_location <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  command_files <- sub("^--file=", "", full[grepl("^--file=", full)])
  frame_files <- Filter(Negate(is.null),
                        lapply(sys.frames(), function(z) z$ofile))
  frame_dirs <- if (length(frame_files)) {
    dirname(unlist(frame_files, use.names = FALSE))
  } else character()
  seed_dirs <- c(dirname(command_files), frame_dirs, getwd())
  candidates <- unique(c(
    seed_dirs,
    file.path(seed_dirs, "01_main_experiment"),
    file.path(dirname(seed_dirs), "01_main_experiment")
  ))
  candidates <- candidates[dir.exists(candidates)]
  required <- c("01_main_experiment.R", "01a_generate_repeated_noise.R",
                "01b_plot_repeated_noise.R",
                "01c_plot_fixed_design_summaries.R")
  hits <- candidates[vapply(candidates, function(path) {
    all(file.exists(file.path(path, required)))
  }, logical(1))]
  if (!length(hits)) {
    stop("Cannot locate the 01_main_experiment directory. Open its run.R ",
         "from the complete project folder.")
  }
  normalizePath(hits[1L], winslash = "/", mustWork = TRUE)
}

here <- script_location()
run <- function(filename) {
  message("Running ", filename)
  path <- file.path(here, filename)
  if (!file.exists(path)) stop("Missing script: ", path)
  source(path, chdir = TRUE, local = new.env(parent = globalenv()),
         echo = FALSE, encoding = "UTF-8")
  invisible(TRUE)
}

run("01_main_experiment.R")
run("01c_plot_fixed_design_summaries.R")
run("01a_generate_repeated_noise.R")
run("01b_plot_repeated_noise.R")
message("01 main experiment and summary figures completed.")
