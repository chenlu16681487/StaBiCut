############################################################
# StaBiCut v2 — Multi-seed stability utilities
# File: R/modules_stability_summary_StaBiCut_v2.R
#
# Scope:
#   Utilities for repeated StaBiCut runs across random seeds and for
#   summarizing, visualizing, and exporting cross-seed stability results.
#
# Included in this file:
#   - repeated multi-seed runner
#   - summary tables across seeds
#   - Excel export of stability results
#   - column dictionary for results_df
#   - stability heatmap
#   - stability dual-panel plot
#
# Dependencies assumed to be already sourced:
#   - modules_core_StaBiCut_v2.R
#   - run_StaBiCut_v2.R (for run_multiseed_stability_v2)
#
# Notes:
#   - This script does not call library() globally.
############################################################

# ==========================================================
# 8. Multi-seed stability utilities
# ==========================================================

#' Run the full StaBiCut pipeline across multiple random seeds
#'
#' @param seeds Integer vector of random seeds.
#' @param n_boot Integer scalar giving the number of bootstrap replicates per run.
#' @param save_each_seed_rds Logical; if TRUE, save each seed-level result object
#'   as an RDS file.
#' @param out_dir Character scalar giving the output directory.
#' @param ... Additional arguments passed to
#'   run_batch_sur_cutpoint_analysis_v2().
#'
#' @return Named list of per-seed StaBiCut run objects, one element per seed.
#' @export
run_multiseed_stability_v2 <- function(
    seeds = 1:20,
    n_boot = 500,
    save_each_seed_rds = TRUE,
    out_dir = "./stability_runs",
    ...
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  all_runs <- vector("list", length(seeds))
  names(all_runs) <- paste0("seed_", seeds)

  for (i in seq_along(seeds)) {
    sd <- seeds[i]
    set.seed(sd)

    res <- run_batch_sur_cutpoint_analysis_v2(
      n_boot = n_boot,
      save_plots = FALSE,
      plot_dir = out_dir,
      ...
    )

    res$seed <- sd
    all_runs[[i]] <- res

    if (save_each_seed_rds) {
      saveRDS(res, file = file.path(out_dir, paste0("results_seed_", sd, ".rds")))
    }
  }

  all_runs
}

#' Summarize cross-seed stability for one target gene and all genes
#'
#' @param all_runs List of per-seed results tables or objects coercible via
#'   dplyr::bind_rows(), typically seed-level `results_df` tables.
#' @param target_gene Character scalar giving the focal gene for TopK frequency
#'   evaluation.
#' @param top_k Integer scalar giving the TopK threshold used for the target
#'   gene stability summary.
#'
#' @return Named list containing:
#'   `target_topk_rate`, `score_summary`, `cut_summary`, and `top_tbl`.
#' @export
summarize_stability_v2 <- function(all_runs, target_gene = "ZG16", top_k = 2) {
  df_all <- dplyr::bind_rows(all_runs)

# 1) Within-seed ranking and target-gene TopK frequency
  top_tbl <- df_all  %>%
    dplyr::group_by(seed) %>%
    dplyr::arrange(dplyr::desc(composite_score), .by_group = TRUE)  %>%
    dplyr::mutate(rank_in_seed = dplyr::row_number())  %>%
    dplyr::ungroup()

  freq_target_topk <- mean(top_tbl$Gene == target_gene & top_tbl$rank_in_seed <= top_k, na.rm = TRUE)

# 2) Cross-seed mean and standard deviation of composite scores by gene
  score_summary <- df_all %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      composite_mean = mean(composite_score, na.rm = TRUE),
      composite_sd = sd(composite_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(composite_mean))

# 3) Cross-seed summary of selected cutoffs: median, IQR, and range
  cut_summary <- df_all %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      cutoff_median = median(Cutoff, na.rm = TRUE),
      cutoff_IQR = IQR(Cutoff, na.rm = TRUE),
      cutoff_range = diff(range(Cutoff, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(cutoff_range))

  list(
    target_topk_rate = freq_target_topk,
    score_summary = score_summary,
    cut_summary = cut_summary,
    top_tbl = top_tbl
  )
}

#' Export cross-seed stability summaries to an Excel workbook
#'
#' @param all_runs List of per-seed results tables, typically `all_runs_slim`.
#' @param out_file Character scalar giving the output workbook path.
#' @param top_ks Integer vector of TopK thresholds to summarize.
#' @param overwrite Logical; if TRUE, overwrite an existing workbook.
#'
#' @return A named list containing the exported summary tables used to build the
#'   Excel workbook, returned invisibly.
#' @export
export_stability_excel_v2 <- function(
    all_runs,
    out_file = "StaBiCut_Stability_Summary.xlsx",
    top_ks = c(2, 3, 5),
    overwrite = TRUE
) {
  stopifnot(length(all_runs) > 0)

  df_all <- dplyr::bind_rows(all_runs)
  if (!"seed" %in% colnames(df_all)) stop("Column 'seed' not found in all_runs.")
  if (!"Gene" %in% colnames(df_all)) stop("Column 'Gene' not found in all_runs.")
  if (!"composite_score" %in% colnames(df_all)) stop("Column 'composite_score' not found in all_runs.")

   # ---- 1) Combine all seed-level results into one long-form table ----
   df_all <- df_all %>%
                mutate(
                  seed = as.integer(seed),
                  Gene = as.character(Gene))

  # ---- 2) Cross-seed summary of composite scores and component scores ----
  score_mean_sd <- df_all %>%
                group_by(Gene) %>%
                summarise(
                  n_seeds = n_distinct(seed),
                  composite_mean = mean(composite_score, na.rm = TRUE),
                  composite_sd   = sd(composite_score, na.rm = TRUE),
                  .groups = "drop"
                ) %>%
                arrange(desc(composite_mean))

  # Optional: also append mean ± SD for individual scoring components
  score_components <- intersect(
    c("bootstrap_score", "hazard_score", "tn_score", "density_score", "balance_score"),
    colnames(df_all)
  )
 if (length(score_components) > 0) {
                comp_tbl <- df_all %>%
                  group_by(Gene) %>%
                  summarise(across(all_of(score_components),
                                   list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))),
                            .groups = "drop")
                score_mean_sd <- left_join(score_mean_sd, comp_tbl, by = "Gene")
              }

  # ---- 3) Cross-seed summary of selected cutoffs: median, IQR, and range ----
  cutoff_median_iqr <- df_all %>%
                group_by(Gene) %>%
                summarise(
                  n_seeds = n_distinct(seed),
                  cutoff_median = median(Cutoff, na.rm = TRUE),
                  cutoff_IQR    = IQR(Cutoff, na.rm = TRUE),
                  cutoff_range  = diff(range(Cutoff, na.rm = TRUE)),
                  .groups = "drop"
                ) %>%
                arrange(desc(cutoff_range))
  
  # ---- 4) Within-seed ranking and gene-level TopK frequency summary ----
  rank_tbl <- df_all %>%
                group_by(seed) %>%
                arrange(desc(composite_score), .by_group = TRUE) %>%
                mutate(rank_in_seed = row_number()) %>%
                ungroup()

  topk_frequency <- lapply(top_ks, function(k) {
                rank_tbl %>%
                  group_by(Gene) %>%
                  summarise(
                    n_seeds = n_distinct(seed),
                    freq = mean(rank_in_seed <= k, na.rm = TRUE),
                    .groups = "drop"
                  ) %>%
                  mutate(top_k = k)
              }) %>%
                bind_rows() %>%
                tidyr::pivot_wider(
                  names_from = top_k,
                  values_from = freq,
                  names_prefix = "top",
                  values_fill = 0
                ) %>%
                arrange(desc(top2), desc(top3), desc(top5))
  
   # ---- 5) Normalize and integrate column dictionaries for results and stability summaries ----
  .normalize_dictionary_df <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.data.frame(x, stringsAsFactors = FALSE)

    if (all(c("Column name", "Category", "Definition") %in% colnames(x))) {
      return(x[, c("Column name", "Category", "Definition"), drop = FALSE])
    }
   
   cn_std <- trimws(colnames(x))

    map_name <- c(
      "Column name" = "Column name",
      "Column.name" = "Column name",
      "column_name" = "Column name",
      "column" = "Column name",
      "name" = "Column name",
      "Category" = "Category",
      "category" = "Category",
      "Definition" = "Definition",
      "definition" = "Definition",
      "Description" = "Definition",
      "description" = "Definition"
    )

    cn_new <- ifelse(cn_std %in% names(map_name), map_name[cn_std], cn_std)
    colnames(x) <- cn_new

    if (!"Column name" %in% colnames(x)) x[["Column name"]] <- NA_character_
    if (!"Category" %in% colnames(x)) x[["Category"]] <- NA_character_
    if (!"Definition" %in% colnames(x)) x[["Definition"]] <- NA_character_
    x[, c("Column name", "Category", "Definition"), drop = FALSE]
  }

  # results_df dictionary
  dict_results <- NULL
  if (exists("build_column_dictionary_resultsdf_v2", mode = "function")) {
                dict_results <- tryCatch(
                  build_column_dictionary_resultsdf_v2(),
                  error = function(e) NULL
                )
              }
              
              if (is.null(dict_results)) {
                dict_results <- data.frame(
                  `Column name` = colnames(df_all),
                  Category = "results_df",
                  Definition = "See StaBiCut v2 output definition.",
                  stringsAsFactors = FALSE
                )
              }

  dict_results <- .normalize_dictionary_df(dict_results)

  # stability summary dictionary
  dict_stability <- data.frame(
    `Column name` = c(
      "seed", "rank_in_seed", "n_seeds", "composite_mean", "composite_sd",
      paste0(score_components, "_mean"), paste0(score_components, "_sd"),
      "cutoff_median", "cutoff_IQR", "cutoff_range", paste0("top", top_ks)
    ),
    Category = c(
      rep("Stability runs", 2),
      rep("Stability summary", 1 + 2 + 2 * length(score_components) + 3 + length(top_ks))
    ),
    Definition = c(
      "Random seed used for one full pipeline run.",
      "Within-seed rank by descending composite_score (1 = best).",
      "Number of seeds included for the gene.",
      "Mean composite_score across seeds.",
      "Standard deviation of composite_score across seeds.",
      if (length(score_components) > 0) paste0("Mean of ", score_components, " across seeds.") else character(0),
      if (length(score_components) > 0) paste0("Standard deviation of ", score_components, " across seeds.") else character(0),
      "Median of per-seed optimal Cutoff across seeds.",
      "Interquartile range of the selected cutoff across seeds.",
      "Range of per-seed optimal Cutoff across seeds (max - min).",
      paste0("Frequency of being ranked within top ", top_ks, " across seeds (0–1).")
    ),
    stringsAsFactors = FALSE
  )

  dict_stability <- .normalize_dictionary_df(dict_stability)

  column_dictionary <- dplyr::bind_rows(dict_results, dict_stability)
              
  if (!"Column name" %in% colnames(column_dictionary)) {
       stop("Column dictionary construction failed: 'Column name' column is missing after binding.")
     }
              
  column_dictionary <- column_dictionary %>%
       dplyr::distinct(`Column name`, .keep_all = TRUE)

  # ---- 6) Create and write the Excel workbook ----
  wb <- openxlsx::createWorkbook()
   # Sheet1
  openxlsx::addWorksheet(wb, "Seed_results_long")
  openxlsx::writeDataTable(wb, "Seed_results_long", df_all)
  # Sheet2
  openxlsx::addWorksheet(wb, "Score_mean_sd")
  openxlsx::writeDataTable(wb, "Score_mean_sd", score_mean_sd)
  # Sheet3
  openxlsx::addWorksheet(wb, "Cutoff_median_IQR")
  openxlsx::writeDataTable(wb, "Cutoff_median_IQR", cutoff_median_iqr)
  # Sheet4
  openxlsx::addWorksheet(wb, "TopK_frequency")
  openxlsx::writeDataTable(wb, "TopK_frequency", topk_frequency)
  # Sheet5
  openxlsx::addWorksheet(wb, "Column_dictionary")
  openxlsx::writeDataTable(wb, "Column_dictionary", column_dictionary)

  openxlsx::saveWorkbook(wb, out_file, overwrite = overwrite)

  invisible(list(
    out_file = out_file,
    df_all = df_all,
    score_mean_sd = score_mean_sd,
    cutoff_median_iqr = cutoff_median_iqr,
    topk_frequency = topk_frequency,
    column_dictionary = column_dictionary
  ))
}

#' Column dictionary for the StaBiCut per-gene results table
#'
#' @return A data.frame documenting the main columns in `results_df`, including
#'   column name, category, and definition.
#' @export
build_column_dictionary_resultsdf_v2 <- function() {
  data.frame(
    `Column name` = c(
      "Gene", "HR", "CI_low", "CI_high", "P",
      "Cutoff", "Cutoff_P", "p_spline",
      "Bootstrap_Median", "Bootstrap_IQR", "Bootstrap_CI_low", "Bootstrap_CI_high", "Bootstrap_CI_range",
      "IQR_rel", "CI_rel", "bootstrap_score",
      "expression_range", "quantile_rank", "density_score",
      "min_group_prop", "low_group_n", "high_group_n",
      "TN_log2FC", "direction_concordance", "hazard_dir_consistency",
      "P_adj", "Cutoff_P_adj",
      "hazard_score", "tn_score", "balance_score",
      "composite_score", "Rank"
    ),
    Category = c(
      rep("Gene annotation", 1),
      rep("Survival statistics", 4),
      rep("Cutoff detection", 3),
      rep("Bootstrap stability", 5),
      rep("Bootstrap stability (relative)", 3),
      rep("Cutoff position / distribution", 3),
      rep("Group robustness", 3),
      rep("Tumor-normal trend", 2),
      rep("Bootstrap hazard consistency", 1),
      rep("Multiple testing", 2),
      rep("Scoring components", 3),
      rep("Integrated score", 2)
    ),
                    Definition = c(
                  "Gene symbol.",
                  "Hazard ratio (HR) for high vs. low groups (Cox PH model with group defined by Cutoff).",
                  "Lower bound of 95% CI of HR.",
                  "Upper bound of 95% CI of HR.",
                  "Wald test P-value from Cox PH model.",
                  
                  "Optimal expression threshold selected by deterministic cutpoint scan (direction constrained when expected_dir is available; minprop applied).",
                  "Log-rank P-value (survdiff) comparing survival between groups defined by Cutoff.",
                  "P-value from likelihood ratio test comparing Cox null vs. ns(expr, df=3) spline model (global nonlinear support).",
                  
                  "Median cutoff across bootstrap resamples.",
                  "Interquartile range (IQR) of bootstrap cutoff distribution.",
                  "Lower bound of 95% CI (2.5th percentile) of bootstrap cutoff distribution.",
                  "Upper bound of 95% CI (97.5th percentile) of bootstrap cutoff distribution.",
                  "Bootstrap_CI_high - Bootstrap_CI_low.",
                  
                  "Bootstrap_IQR normalized by tumor expression IQR (IQR_rel = boot_IQR / IQR(expr)).",
                  "Bootstrap CI width normalized by tumor expression IQR (CI_rel = (CI_high-CI_low)/IQR(expr)).",
                  "Stability score mapped from IQR_rel and CI_rel: exp(-k1*IQR_rel) * exp(-k2*CI_rel). Higher = more stable.",
                  
                  "Robust expression range in tumor samples (95th–5th percentile).",
                  "Empirical CDF quantile position of Cutoff in tumor expression distribution (0–1).",
                  "Cutoff plausibility score: 1 - |quantile_rank - 0.5| (higher = closer to central density).",
                  
                  "Minimum proportion of samples in either group at Cutoff (min(low, high)/total).",
                  "Number of samples in low group (expr ≤ Cutoff).",
                  "Number of samples in high group (expr > Cutoff).",
                  
                  "Median(tumor) - median(normal) on log2(TPM+1) scale.",
                  "Sign agreement between TN_log2FC and log(HR): +1 concordant, -1 discordant, NA otherwise.",
                  "Proportion of bootstrap cutoffs whose inferred HR direction matches the reference HR direction (sign(log(HR))).",
                  
                  "BH-adjusted P-value for Cox Wald P.",
                  "BH-adjusted P-value for Cutoff_P (log-rank).",
                  
                  "Alias of hazard_dir_consistency used in scoring (0–1).",
                  "Tumor–normal concordance score used in scoring (1 if direction_concordance==1; 0 if -1).",
                  "Group balance score used in scoring: min(1, max(0, min_group_prop*2)).",
                  
                  "Weighted composite score combining bootstrap_score, hazard_score, tn_score, density_score, balance_score.",
                  "Rank by descending composite_score (ties use rank(..., ties.method='min'))."
                ),
    stringsAsFactors = FALSE
  )
}

#' Plot a stability heatmap across seeds
#'
#' @param all_runs List of per-seed results tables containing at least `Gene`,
#'   `seed`, and the selected `value_col`.
#' @param value_col Character scalar giving the column to visualize.
#' @param out_pdf Optional character scalar giving the output PDF path.
#' @param width,height Numeric scalars passed to `ggplot2::ggsave()`.
#'
#' @return A ggplot object, returned invisibly.
#' @export
plot_stability_heatmap_v2 <- function(all_runs, value_col = "composite_score", out_pdf = NULL, width = 8, height = 5) {
              df_all <- dplyr::bind_rows(all_runs)

              if (!"seed" %in% colnames(df_all)) {
                stop("No column 'seed' found in all_runs. Ensure run_multiseed_stability_v2 adds res$seed.")
              }
              if (!value_col %in% colnames(df_all)) {
                stop("value_col not found: ", value_col)
              }
              
              df_all <- df_all %>%
                mutate(
                  seed = as.factor(seed),
                  Gene = as.character(Gene)
                )
              
              gene_order <- df_all %>%
                group_by(Gene) %>%
                summarise(mu = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
                arrange(desc(mu)) %>%
                pull(Gene)
              
              df_all$Gene <- factor(df_all$Gene, levels = rev(gene_order))
              df_all$seed <- factor(df_all$seed, levels = as.character(sort(unique(df_all$seed))))
              
              p <- ggplot(df_all, aes(x = seed, y = Gene, fill = .data[[value_col]])) +
                geom_tile(color = "white", linewidth = 0.3) +
                scale_fill_viridis_c(option = "B", direction = -1, name = value_col) +
                theme_minimal(base_size = 13) +
                theme(
                  axis.title = element_blank(),
                  panel.grid = element_blank(),
                  legend.position = "right"
                ) +
                labs(title = "Stability heatmap across random seeds")
              
              if (!is.null(out_pdf)) {
                ggsave(out_pdf, p, width = width, height = height)
              }
              invisible(p)
            }


#' Plot a dual-panel summary of cross-seed stability
#'
#' @param all_runs List of per-seed results tables containing at least `Gene`,
#'   `seed`, and the selected `value_col`.
#' @param value_col Character scalar giving the column visualized in the heatmap
#'   panel.
#' @param top_ks Integer vector of TopK thresholds used in the ranking-stability
#'   panel.
#' @param out_pdf Character scalar giving the output PDF path.
#' @param width,height Numeric scalars passed to `ggplot2::ggsave()`.
#'
#' @return A named list containing the assembled plot and the
#'   underlying summary tables, returned invisibly.
#' @export
plot_stability_dualpanel_v2 <- function(
    all_runs,
    value_col = "composite_score",
    top_ks = c(1, 2, 3, 5),
    out_pdf = "stability_dualpanel_v2.pdf",
    width = 11,
    height = 5.5
) {
              # ----------------------------
              # 1) bind runs
              # ----------------------------
              df_all <- dplyr::bind_rows(all_runs)
              
              stopifnot("seed" %in% colnames(df_all))
              stopifnot("Gene" %in% colnames(df_all))
              stopifnot(value_col %in% colnames(df_all))
              
              df_all <- df_all %>%
                mutate(
                  seed = as.integer(seed),
                  Gene = as.character(Gene)
                )
              
              # ----------------------------
              # 2) gene order by mean composite score
              #    highest score shown at top
              # ----------------------------
              gene_order <- df_all %>%
                group_by(Gene) %>%
                summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
                arrange(desc(mean_value)) %>%
                pull(Gene)
              
              gene_levels <- rev(gene_order)  # ggplot y-axis: last level appears at top
              
              df_all <- df_all %>%
                mutate(Gene = factor(Gene, levels = gene_levels))
              
              # ----------------------------
              # 3) Panel A: stability heatmap
              # ----------------------------
              pA <- ggplot(df_all, aes(x = factor(seed), y = Gene, fill = .data[[value_col]])) +
                geom_tile(color = "white", linewidth = 0.35) +
                scale_fill_gradientn(
                  colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5"),
                  limits = c(min(df_all[[value_col]], na.rm = TRUE),
                             max(df_all[[value_col]], na.rm = TRUE)),
                  name = "Composite\nscore"
                ) +
                labs(
                  x = "Random seed",
                  y = NULL,
                  title = "A  Stability heatmap across seeds"
                ) +
                theme_bw(base_size = 11) +
                theme(
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(face = "italic"),
                  plot.title = element_text(face = "bold", size = 12),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 9)
                )
              
              # ----------------------------
              # 4) compute TopK frequency
              # ----------------------------
              rank_tbl <- df_all %>%
                group_by(seed) %>%
                arrange(desc(.data[[value_col]]), .by_group = TRUE) %>%
                mutate(rank_in_seed = row_number()) %>%
                ungroup()
              
              topk_frequency <- lapply(top_ks, function(k) {
                rank_tbl %>%
                  group_by(Gene) %>%
                  summarise(
                    n_seeds = n_distinct(seed),
                    freq = mean(rank_in_seed <= k, na.rm = TRUE),
                    .groups = "drop"
                  ) %>%
                  mutate(top_k = paste0("top", k))
              }) %>%
                bind_rows()
              
              topk_frequency <- topk_frequency %>%
                mutate(
                  Gene = factor(as.character(Gene), levels = gene_levels),
                  top_k = factor(top_k, levels = paste0("top", top_ks))
                )
              
              # ----------------------------
              # 5) Panel B: TopK mini-heatmap
              # ----------------------------
              pB <- ggplot(topk_frequency, aes(x = top_k, y = Gene, fill = freq)) +
                geom_tile(color = "white", linewidth = 0.5) +
                geom_text(aes(label = sprintf("%.2f", freq)), size = 3) +
                scale_fill_gradientn(
                  colours = c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d"),
                  limits = c(0, 1),
                  name = "TopK\nfrequency"
                ) +
                labs(
                  x = NULL,
                  y = NULL,
                  title = "B  Ranking stability summary"
                ) +
                theme_bw(base_size = 11) +
                theme(
                  panel.grid = element_blank(),
                  axis.text.x = element_text(face = "bold"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(face = "bold", size = 12),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 9)
                )
              
              # ----------------------------
              # 6) combine
              # ----------------------------
              p <- pA + pB +
                plot_layout(widths = c(2.3, 1.1), guides = "collect") &
                theme(legend.position = "right")
              
              ggsave(out_pdf, p, width = width, height = height)
              return(invisible(list(
                plot = p,
                heatmap_data = df_all,
                topk_data = topk_frequency
              )))
            }
