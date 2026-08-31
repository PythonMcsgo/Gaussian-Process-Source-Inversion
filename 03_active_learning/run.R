
script_location <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  command_files <- sub("^--file=", "", full[grepl("^--file=", full)])
  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(z) z$ofile))
  frame_dirs <- if (length(frame_files)) dirname(unlist(frame_files, use.names=FALSE)) else character()
  seed_dirs <- c(dirname(command_files), frame_dirs, getwd())
  candidates <- unique(c(seed_dirs, file.path(seed_dirs, "03_active_learning"),
                         file.path(dirname(seed_dirs), "03_active_learning")))
  candidates <- candidates[dir.exists(candidates)]
  hits <- candidates[vapply(candidates, function(path) {
    all(file.exists(file.path(path, c("active_learning.R",
                          "assemble_6sources.R"))))
  }, logical(1))]
  if (!length(hits)) stop("Cannot locate the 03_active_learning directory.")
  normalizePath(hits[1], winslash="/", mustWork=TRUE)
}

here <- script_location()

# Report-only entry point. The six source-specific experiments must already be
# present under ../results/active_single_<source>. This script does not rerun
# any inverse problem; it only rebuilds the combined tables and figures.
sys.source(
  file.path(here, "assemble_6sources.R"),
  envir = new.env(parent = globalenv())
)
message("Six-source active-learning report completed.")
