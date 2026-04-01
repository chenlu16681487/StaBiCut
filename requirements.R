# requirements.R
# StaBiCut repository dependency helper
# This script checks and optionally installs the R packages used by the
# current script-style StaBiCut repository.

core_packages <- c(
  "survival",
  "survminer",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "openxlsx",
  "splines",
  "forcats"
)

optional_packages <- c(
  "mixtools",
  "forestplot",
  "fmsb"
)

all_packages <- unique(c(core_packages, optional_packages))

install_missing_packages <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = repos)
  } else {
    message("No missing packages detected.")
  }
}

check_packages <- function(pkgs) {
  data.frame(
    package = pkgs,
    installed = vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    stringsAsFactors = FALSE
  )
}

export_session_info <- function(file = "sessionInfo.txt") {
  writeLines(capture.output(sessionInfo()), con = file)
  message("sessionInfo() written to: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

# Optional renv bootstrap:
# install.packages("renv")
# renv::init(bare = TRUE)     # run once in the repository root
# renv::snapshot()            # writes renv.lock based on currently used packages

message("=== StaBiCut package check ===")
print(check_packages(core_packages))
message("=== Optional packages ===")
print(check_packages(optional_packages))

# Uncomment to install anything missing:
# install_missing_packages(all_packages)

# Uncomment to export current package/session metadata:
# export_session_info("sessionInfo.txt")
