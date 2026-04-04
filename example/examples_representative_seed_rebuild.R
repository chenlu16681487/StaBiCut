############################################################
# StaBiCut v2 — Representative-seed reconstruction example
# File: examples/representative_seed_rebuild.R
#
# Purpose:
# This script reconstructs representative-seed outputs from a
# previously completed multi-seed StaBiCut run.
#
# Recommended use:
# 1) Run examples/main_crc_tcga.R first with run_multiseed=TRUE.
# 2) Load the saved all_runs / all_runs_slim objects.
# 3) Select a representative seed for a target gene.
# 4) Rebuild the single-gene main panel and multi-gene sheets.
#
# Scope of this example:
# - Representative-seed selection.
# - Re-loading the full result object for the chosen seed.
# - Rebuilding a target-gene panel.
# - Rebuilding ranked multi-gene summary sheets.
############################################################

# ==========================================================
# 0. User-defined paths and reconstruction settings
# ==========================================================

project_dir <- "."
data_dir <- "./example_data"
run_dir <- file.path("./output", "crc_tcga_example")
stability_dir <- file.path(run_dir, "stability_runs")

expr_rdata <- file.path(data_dir, "exprset_11DEGs_tcga_crc.Rdata")
clin_rdata <- file.path(data_dir, "clinicalSE_11DEGs_tcga_crc.Rdata")
all_runs_rdata <- file.path(stability_dir, "all_runs.RData")
all_runs_slim_rdata <- file.path(stability_dir, "all_runs_slim.RData")

core_file   <- file.path(project_dir, "scripts", "modules_core_StaBiCut_v2.R")
plot_file   <- file.path(project_dir, "scripts", "modules_plot_single_StaBiCut_v2.R")
stab_file   <- file.path(project_dir, "scripts", "modules_stability_summary_StaBiCut_v2.R")
seed_file   <- file.path(project_dir, "scripts", "modules_seed_selection_StaBiCut_v2.R")
helper_file <- file.path(project_dir, "scripts", "multi-gene panel_helper_StaBiCut_v2.R")
run_file    <- file.path(project_dir, "scripts", "run_StaBiCut_v2.R")

gene_x <- "ZG16"
force_direction <- FALSE

# Same candidate-gene set and prior table used in the CRC run.
geneset <- c(
  "SPOCK2", "PYCR1", "CA4", "CES1", "ABCB1", "ZG16",
  "TNXB", "HMCN2", "MEP1A", "SLC37A2", "CHGB"
)

gene_prior_table <- data.frame(
  Gene = geneset,
  expected_dir = c(rep("adverse_high", 2), rep("protective_high", 9)),
  stringsAsFactors = FALSE
)

# ==========================================================
# 1. Load required packages and StaBiCut source files
# ==========================================================

required_pkgs <- c(
  "survival", "dplyr", "survminer", "ggplot2",
  "patchwork", "openxlsx", "splines", "tidyr", "forcats"
)

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

source(core_file, encoding = "UTF-8")
source(plot_file, encoding = "UTF-8")
source(stab_file, encoding = "UTF-8")
source(seed_file, encoding = "UTF-8")
source(helper_file, encoding = "UTF-8")
source(run_file, encoding = "UTF-8")

message("StaBiCut v2 source files loaded successfully.")

# ==========================================================
# 2. Load input data and saved multi-seed results
# ==========================================================

load(expr_rdata)
load(clin_rdata)
stopifnot(exists("mrna_expr_tpm"), exists("clinicalSE"))

load(all_runs_rdata)
load(all_runs_slim_rdata)
stopifnot(exists("all_runs"), exists("all_runs_slim"))

dir.create(file.path(run_dir, "representative_seed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "representative_seed", "single_gene_panel"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "representative_seed", "multi_gene_sheets"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "representative_seed", "cutscan_exports"), recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# 3. Recompute expected-direction summary
# ==========================================================

expected_dir_table <- get_expected_dir_table_v2(
  exprset = mrna_expr_tpm,
  clin = clinicalSE,
  geneset = geneset,
  gene_prior_table = gene_prior_table,
  force_direction = force_direction
)

write.csv(
  expected_dir_table,
  file = file.path(run_dir, "representative_seed", "expected_dir_table.csv"),
  row.names = FALSE
)

# ==========================================================
# 4. Select a representative seed for the target gene
# ==========================================================

# Representative-seed selection is performed on the slim
# multi-seed result object to identify a run that is close to
# the cross-seed central tendency rather than an extreme run.
rep_excel_res <- export_rep_seed_excel_v2(
  all_runs_slim = all_runs_slim,
  gene_x = gene_x,
  out_file = file.path(run_dir, "representative_seed", paste0("representative_seed_", gene_x, ".xlsx")),
  overwrite = TRUE
)

best_seed <- rep_excel_res$representative_seed
message("Representative seed selected for ", gene_x, ": ", best_seed)

# ==========================================================
# 5. Load the full result object for the representative seed
# ==========================================================

rep_panel_obj <- load_rep_seed_full_result_v2(
  representative_seed = best_seed,
  gene_x = gene_x,
  rds_dir = stability_dir
)

expected_dir_x <- expected_dir_table$expected_dir[
  match(rep_panel_obj$gene_x, expected_dir_table$Gene)
]

# ==========================================================
# 6. Rebuild the target-gene main panel
# ==========================================================
# cutoff-stability plotting: A. cut scan、B. expression density、C. bootstrap density、D. KM plot
panel_main_res <- plot_gene_panel_main_v2(
  df = rep_panel_obj$gene_df,
  gene = rep_panel_obj$gene_x,
  best_cut = rep_panel_obj$best_cut,
  hr = rep_panel_obj$gene_result_row$HR[1],
  ci_low = rep_panel_obj$gene_result_row$CI_low[1],
  ci_high = rep_panel_obj$gene_result_row$CI_high[1],
  cutoff_p = rep_panel_obj$gene_result_row$Cutoff_P[1],
  boot = rep_panel_obj$boot_vec,
  scan_df = rep_panel_obj$scan_df,
  expected_dir = expected_dir_x,
  plot_dir = file.path(run_dir, "representative_seed", "single_gene_panel"),
  width = 12,
  height = 10,
  xps_aspect = 0.625,
  km_aspect = 1
)

# Optional spline-extended panel for the same target gene.
plot_gene_panel_with_spline_v2(
  pan = panel_main_res,
  df = rep_panel_obj$gene_df,
  gene = rep_panel_obj$gene_x,
  plot_dir = file.path(run_dir, "representative_seed", "single_gene_panel"),
  xps_aspect = 0.625,
  width = 15,
  height = 18
)

# ==========================================================
# 7. Export one-seed cut-scan PDFs for all candidate genes
# ==========================================================

export_one_seed_cutscan_pdfs_v2(
  all_runs = all_runs,
  seed_index = best_seed,
  geneset = geneset,
  expected_dir_table = expected_dir_table,
  out_dir = file.path(run_dir, "representative_seed", "cutscan_exports")
)

# ==========================================================
# 8. Rebuild ranked multi-gene sheets from the representative seed
# ==========================================================

# The gene display order was determined from the across-seed mean ranking, 
# whereas the actual panel plots were reconstructed from the representative seed.
res_mean_plot <- plot_multiseed_composite_summary_v2(
  all_runs = all_runs,
  out_pdf = file.path(run_dir, "representative_seed", "multiseed_composite_mean_summary.pdf"),
  top_n = NULL,
  order_by = "composite_mean",
  p_col = "Cutoff_P",
  p_agg = "median",
  use_adjusted_p = FALSE
)

#Plot the remaining 10 DEGs
full_result <- rep_panel_obj$full_result
genes_ranked <- as.character(res_mean_plot$summary_table$Gene)

genes_sheet_1 <- genes_ranked[c(1, 3:6)]
genes_sheet_2 <- genes_ranked[7:11]

plot_multi_gene_sheet_main_v2(
  gene_order = genes_sheet_1,
  results_df = full_result$results_df,
  df_cache = full_result$df_cache,
  boot_cache = full_result$boot_cache,
  scan_cache = full_result$scan_cache,
  expected_dir_table = expected_dir_table,
  out_pdf = file.path(run_dir, "representative_seed", "multi_gene_sheets", "stability_5genes_main_1.pdf"),
  width = 25.5,
  height = 24,
  widths = c(1.15, 1.15, 1.15, 0.72),
  km_aspect = 1,
  xps_aspect = 0.625,
  show_row_tags = TRUE
)

plot_multi_gene_sheet_main_v2(
  gene_order = genes_sheet_2,
  results_df = full_result$results_df,
  df_cache = full_result$df_cache,
  boot_cache = full_result$boot_cache,
  scan_cache = full_result$scan_cache,
  expected_dir_table = expected_dir_table,
  out_pdf = file.path(run_dir, "representative_seed", "multi_gene_sheets", "stability_5genes_main_2.pdf"),
  width = 25.5,
  height = 24,
  widths = c(1.15, 1.15, 1.15, 0.72),
  km_aspect = 1,
  xps_aspect = 0.625,
  show_row_tags = TRUE
)

# ==========================================================
# 9. Optional reproducibility record
# ==========================================================
save_session_info <- TRUE

if (isTRUE(save_session_info)) {
  writeLines(
    capture.output(sessionInfo()),
    con = file.path(out_dir, "sessionInfo_main_crc_tcga.txt")
  )
}
message("Representative-seed reconstruction completed successfully.")
