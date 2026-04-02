############################################################
# StaBiCut v2 — CRC/TCGA example workflow
# File: examples/main_crc_tcga.R
#
# Purpose:
# This script provides a publication-oriented but still
# user-runnable example for applying StaBiCut v2 to a
# predefined candidate-gene set in a TCGA-style CRC cohort.
#
# Recommended use:
# 1) Adapt the paths below to the local repository and data.
# 2) Load a prepared expression matrix and matched clinical data.
# 3) Run one full StaBiCut analysis on the candidate genes.
# 4) Optionally run cross-seed robustness analysis.
# 5) Export result tables and summary figures.
#
# Scope of this example:
# - Single-run StaBiCut analysis.
# - Composite-score summary plot.
# - Expected-direction summary table.
# - Optional multi-seed stability assessment.
#
# Notes:
# - This script is designed for external users running the
#   repository examples/ workflow.
# - It does not reconstruct representative-seed figure panels;
#   those steps are handled in examples/representative_seed_rebuild.R.
############################################################

# ==========================================================
# 0. User-defined paths and analysis configuration
# ==========================================================

# Root directory containing the StaBiCut source scripts.
project_dir <- "."

# Directory containing the prepared input data object.
data_dir <- "./data"

# Output directory for this example run.
out_dir <- file.path("./output", "crc_tcga_example")

# Input RData expected to contain:
#   - mrna_expr_tpm : numeric matrix, genes x samples
#   - clinicalSE    : data.frame with matched clinical variables
input_rdata <- file.path(data_dir, "TCGA-CRC-sur.Rdata")

# StaBiCut source files.
mod_file   <- file.path(project_dir, "R", "modules_StaBiCut_v2.R")
run_file   <- file.path(project_dir, "R", "run_StaBiCut_v2.R")
panel_file <- file.path(project_dir, "R", "Panel_helper_StaBiCut_v2.R")

# Core analysis parameters.
analysis_seed <- 1L
n_boot <- 1000L
minprop <- 0.25
force_direction <- FALSE

# Optional cross-seed robustness block.
run_multiseed <- TRUE
stability_seeds <- 1:20

# Candidate genes evaluated in the CRC application.
geneset <- c(
  "SPOCK2", "PYCR1", "CA4", "CES1", "ABCB1", "ZG16",
  "TNXB", "HMCN2", "MEP1A", "SLC37A2", "CHGB"
)

# Gene-specific directional priors.
gene_prior_table <- data.frame(
  Gene = geneset,
  expected_dir = c(rep("adverse_high", 2), rep("protective_high", 9)),
  stringsAsFactors = FALSE
)

# ==========================================================
# 1. Load required packages and StaBiCut source files
# ==========================================================

# Core packages required by the currently visible StaBiCut v2
# runner/modules/helper scripts in this repository snapshot.
required_pkgs <- c(
  "survival", "dplyr", "survminer", "ggplot2",
  "patchwork", "openxlsx", "splines", "tidyr", "forcats"
)

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

source(mod_file, encoding = "UTF-8")
source(run_file, encoding = "UTF-8")
source(panel_file, encoding = "UTF-8")

message("StaBiCut v2 source files loaded successfully.")

# ==========================================================
# 2. Load input data
# ==========================================================

load(input_rdata)
stopifnot(exists("mrna_expr_tpm"), exists("clinicalSE"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(out_dir, "single_run_panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "stability_runs"), recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# 3. Run one full StaBiCut analysis
# ==========================================================

# This step performs the core StaBiCut workflow: expression
# preprocessing, tumor-only survival analysis, deterministic
# direction-aware cutoff scanning, local bootstrap re-support,
# and integrated five-component scoring.
set.seed(analysis_seed)

res_main <- run_batch_sur_cutpoint_analysis_v2(
  exprset = mrna_expr_tpm,
  geneset = geneset,
  clin = clinicalSE,
  gene_prior_table = gene_prior_table,
  force_direction = force_direction,
  n_boot = n_boot,
  minprop = minprop,
  save_boot_rds = FALSE,
  save_plots = TRUE,
  plot_dir = file.path(out_dir, "single_run_panels"),
  seed = analysis_seed
)

results_df <- res_main$results_df
print(results_df)

write.csv(
  results_df,
  file = file.path(out_dir, "StaBiCut_crc_tcga_results.csv"),
  row.names = FALSE
)

saveRDS(
  res_main,
  file = file.path(out_dir, "StaBiCut_crc_tcga_single_run.rds")
)

# ==========================================================
# 4. Export compact single-run summaries
# ==========================================================

# Composite-score overview for the current run.
plot_cutoff_composite_summary_v2(
  results_df = results_df,
  out_pdf = file.path(out_dir, "StaBiCut_crc_tcga_composite_summary.pdf"),
  top_n = NULL,
  order_by = "composite_score"
)

# Effective expected-direction summary used by the cutoff scan.
expected_dir_table <- get_expected_dir_table_v2(
  exprset = mrna_expr_tpm,
  clin = clinicalSE,
  geneset = geneset,
  gene_prior_table = gene_prior_table,
  force_direction = force_direction
)

write.csv(
  expected_dir_table,
  file = file.path(out_dir, "StaBiCut_crc_tcga_expected_dir_table.csv"),
  row.names = FALSE
)

# ==========================================================
# 5. Optional cross-seed stability analysis
# ==========================================================

# This block is intended for users who want to evaluate whether
# cross-gene prioritization remains stable across repeated seeds.
if (isTRUE(run_multiseed)) {
  all_runs <- run_multiseed_stability_v2(
    seeds = stability_seeds,
    n_boot = n_boot,
    exprset = mrna_expr_tpm,
    geneset = geneset,
    clin = clinicalSE,
    gene_prior_table = gene_prior_table,
    force_direction = force_direction,
    minprop = minprop,
    save_each_seed_rds = TRUE,
    out_dir = file.path(out_dir, "stability_runs")
  )

  save(
    all_runs,
    file = file.path(out_dir, "stability_runs", "all_runs.RData")
  )

  all_runs_slim <- lapply(seq_along(all_runs), function(i) {
    x <- all_runs[[i]]
    if (is.null(x$results_df) || nrow(x$results_df) == 0) return(NULL)
    df <- x$results_df
    df$seed <- x$seed
    df
  })
  all_runs_slim <- Filter(Negate(is.null), all_runs_slim)

  save(
    all_runs_slim,
    file = file.path(out_dir, "stability_runs", "all_runs_slim.RData")
  )

  export_stability_excel_v2(
    all_runs = all_runs_slim,
    stab_summary = NULL,
    out_file = file.path(out_dir, "stability_runs", "StaBiCut_stability_summary.xlsx"),
    top_ks = c(1, 2, 3, 5),
    overwrite = TRUE
  )

  plot_stability_dualpanel_v2(
    all_runs = all_runs_slim,
    value_col = "composite_score",
    top_ks = c(1, 2, 3, 5),
    out_pdf = file.path(out_dir, "stability_runs", "stability_dualpanel.pdf"),
    width = 11,
    height = 5.5
  )

  plot_multiseed_composite_summary_v2(
    all_runs = all_runs,
    out_pdf = file.path(out_dir, "stability_runs", "multiseed_composite_mean_summary.pdf"),
    top_n = NULL,
    order_by = "composite_mean",
    p_col = "Cutoff_P",
    p_agg = "median",
    use_adjusted_p = FALSE
  )
}

# ==========================================================
# 6. Optional reproducibility record
# ==========================================================
save_session_info <- TRUE

if (isTRUE(save_session_info)) {
  writeLines(
    capture.output(sessionInfo()),
    con = file.path(out_dir, "sessionInfo_main_crc_tcga.txt")
  )
}
message("StaBiCut CRC/TCGA example workflow completed successfully.")
