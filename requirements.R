# requirements.R
# StaBiCut repository dependency helper
# This script checks and optionally installs the R packages used by the
# current script-style StaBiCut repository.

all_packages <- c(
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

message("=== StaBiCut package check ===")
print(check_packages(all_packages))

# Uncomment to install anything missing:
# install_missing_packages(all_packages)

# Uncomment to export current package/session metadata:
# export_session_info("sessionInfo.txt")

