############################################################
# StaBiCut v2 — Main runner
# File: R/run_StaBiCut_v2.R
#
# This script defines the main pipeline function:
#   run_batch_sur_cutpoint_analysis_v2()
#
# Expected sourcing order in a script or example:
#   source("R/modules_core_StaBiCut_v2.R")
#   source("R/modules_plot_single_StaBiCut_v2.R")
#   source("R/Multi-gene panel_helper_StaBiCut_v2.R")      # optional
#   source("R/modules_stability_summary_StaBiCut_v2.R") # optional
#   source("R/modules_seed_selection_StaBiCut_v2.R")    # optional
#   source("R/run_StaBiCut_v2.R")
#
# Design notes:
#   - This file contains the primary per-gene analysis loop and composite
#     scoring logic.
#   - It assumes the analytical modules have already been sourced.
#   - It should not call rm(list = ls()), library(), or load().
############################################################

#' Run StaBiCut v2 for a set of candidate genes
#'
#' @param exprset Numeric expression matrix on TPM scale with genes in rows
#'   and samples in columns.
#' @param geneset Character vector of gene symbols to evaluate.
#' @param clin Data frame of clinical annotations, with rows matched to the
#'   expression-matrix columns by sample identifier.
#' @param gene_prior_table Optional data frame with columns Gene and
#'   expected direction.
#' @param force_direction Logical; if `TRUE`, enables a stress-testing
#'   direction-inference branch when no valid gene-specific prior is available.
#'   In this branch, expected direction is assigned from the sign of the
#'   tumor-normal shift (`TN_log2FC`) alone.
#' @param n_boot Integer number of bootstrap iterations per gene.
#' @param adjust_method Multiple-testing correction method passed to
#'   stats::p.adjust(). Default: "BH".
#' @param minprop Minimum allowed group proportion during cutoff scanning.
#' @param score_threshold Optional composite-score threshold used when
#'   generating an additional top-gene summary plot.
#' @param save_boot_rds Logical; if TRUE, saves retained bootstrap cutoffs for
#'   each gene.
#' @param boot_dir Output directory for bootstrap RDS files. Defaults to a
#'   subdirectory of plot_dir.
#' @param seed Optional integer recorded in saved bootstrap RDS files.
#' @param save_plots Logical; if TRUE, writes single-gene and summary plots.
#' @param plot_dir Output directory for figures.
#'
#' @return Named list with elements:
#'   - results_df: per-gene summary table
#'   - df_cache: per-gene analysis tables used for modeling and plotting
#'   - boot_cache: per-gene bootstrap objects
#'   - scan_cache: per-gene cutoff-scan tables
#'   - exprset_full: full log2(TPM + 1) data frame after low-expression filtering
#'   - exprset_tumor: tumor-only log2(TPM + 1) data frame used for survival analysis
#' @export
run_batch_sur_cutpoint_analysis_v2 <- function(
    exprset, geneset, clin,
    gene_prior_table = NULL,
    force_direction = FALSE,
    n_boot = 500,
    adjust_method = "BH",
    minprop = 0.25,
    score_threshold = 0.6,
    save_boot_rds = FALSE,
    boot_dir = NULL,
    seed = NA_integer_,
    save_plots = TRUE,
    plot_dir = "./output"
) {
  # ---------- I/O setup and output directories ----------
  if (is.null(boot_dir)) boot_dir <- file.path(plot_dir, "boot_rds")
  if (!dir.exists(plot_dir) && save_plots) dir.create(plot_dir, recursive = TRUE)
  
  # ---------- Caches for downstream plotting and reconstruction ----------
  df_cache   <- list()
  boot_cache <- list()
  scan_cache <- list()
  
  results <- list()
  
  # ---------- 1) Expression preprocessing ----------
  # log2(TPM+1) transformation
  exprset <- log2(exprset + 1)
  
  # remove genes with very low mean expression to reduce noise
  exprset <- exprset[rowMeans(exprset) > 0.5, , drop = FALSE]
  exprset_full <- exprset
  
  # ---------- 2) Tumor-only survival subset and endpoint construction ----------
  is_tumor <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
  exprset_tumor <- exprset[, is_tumor, drop = FALSE]
  clin_tumor <- clin[is_tumor, , drop = FALSE]
  
  # build survival endpoints (months)
  clin_tumor$time <- ifelse(!is.na(clin_tumor$days_to_death),
                            clin_tumor$days_to_death,
                            clin_tumor$days_to_last_follow_up) / 30
  clin_tumor$event <- ifelse(clin_tumor$vital_status %in% c("Dead", "dead", "DEAD"), 1, 0)
  clin_tumor <- clin_tumor[!is.na(clin_tumor$time) & !is.na(clin_tumor$event), , drop = FALSE]
  
  # enforce matched sample IDs
  sample_ids <- intersect(colnames(exprset_tumor), rownames(clin_tumor))
  exprset_tumor <- exprset_tumor[, sample_ids, drop = FALSE]
  clin_tumor <- clin_tumor[sample_ids, , drop = FALSE]
  
  message("Number of valid tumor samples: ", nrow(clin_tumor))
  
  # ---------- 3) Per-gene StaBiCut analysis ----------
  for (gene in geneset) {
    if (!gene %in% rownames(exprset_tumor)) {
      message("Gene not found in tumor expression matrix: ", gene)
      next
    }
    
    expr <- as.numeric(exprset_tumor[gene, ])
    
    # minimal per-gene analysis table
    df <- data.frame(
      expr  = expr,
      time  = as.numeric(clin_tumor$time),
      event = as.numeric(clin_tumor$event)
    )
    
    # winsorize tumor expression to reduce the influence of extreme outliers
    df$expr <- winsorize(df$expr, lower = 0.05, upper = 0.95)
    expr_vec <- df$expr
    
    # tumor/normal vectors for direction prior and TN trend
    expr_full_gene <- as.numeric(exprset_full[gene, ])
    sample_type_full <- as.numeric(substr(colnames(exprset_full), 14, 15))
    tumor_expr  <- expr_full_gene[sample_type_full < 10]
    normal_expr <- expr_full_gene[sample_type_full >= 10]
    
    expected_dir <- infer_expected_direction(
      gene = gene,
      df = df,
      tumor_values = tumor_expr,
      normal_values = normal_expr,
      gene_prior_table = gene_prior_table,
      force_direction = force_direction
    )
    
    scan_res <- scan_cutpoints_v2(
      df = df,
      expected_dir = expected_dir,
      minprop = minprop,
      grid = "quantile",
      min_events_per_group = 5
    )
    best_cut <- scan_res$best_cut
    if (!is.finite(best_cut)) {
      message("No valid cutoff for ", gene, " (", scan_res$reason, ")")
      next
    }
    
    # define groups at the selected cutoff
    df$group <- ifelse(df$expr > best_cut, "High", "Low")
    df$group <- factor(df$group, levels = c("Low", "High"))
    
    # Dichotomized Cox PH model at the selected cutoff
    fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ group, data = df), error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)
    
    hr      <- as.numeric(s$coefficients[,"exp(coef)"])[1]
    p_val   <- as.numeric(s$coefficients[,"Pr(>|z|)"])[1]
    ci_low  <- as.numeric(s$conf.int[,"lower .95"])[1]
    ci_high <- as.numeric(s$conf.int[,"upper .95"])[1]
    
    # log-rank p-value at cutoff
    surv_test <- survival::survdiff(survival::Surv(time, event) ~ group, data = df)
    cutoff_p <- 1 - stats::pchisq(surv_test$chisq, df = 1)
    
    # local bootstrap re-support of the selected cutoff under the direction constraint
    boot <- bootstrap_cutoffs_v2(
      df = df,
      best_cut = best_cut,
      expected_dir = expected_dir,
      n_boot = n_boot,
      minprop = minprop,
      grid = "quantile",
      n_grid = 120,
      restrict_window = TRUE,
      window_iqr_mult = 0.5,
      min_valid = 50,
      min_events_per_group = 5
    )
    
    # ---------- cache per-gene objects (for later multi-gene panels) ----------
    df_cache[[gene]]   <- df
    boot_cache[[gene]] <- boot
    scan_cache[[gene]] <- scan_res$scan_df
    
    # optional: save boot RDS per gene
    if (save_boot_rds) {
      if (!dir.exists(boot_dir)) dir.create(boot_dir, recursive = TRUE)
      saveRDS(
        list(
          gene = gene,
          seed = seed,
          best_cut = best_cut,
          expected_dir = expected_dir,
          n_boot = n_boot,
          minprop = minprop,
          window_iqr_mult = 0.5,
          boot_cuts = boot$boot_cuts
        ),
        file = file.path(boot_dir, paste0(gene, "_bootcuts_n", n_boot, ".rds"))
      )
    }
    
    # compute the bootstrap stability score from IQR_rel and CI_rel, normalized to the tumor-expression IQR
    k1 <- 3; k2 <- 2
    expr_iqr <- stats::IQR(df$expr, na.rm = TRUE)
    if (!is.finite(expr_iqr) || expr_iqr <= 0 ||
        !is.finite(boot$iqr) || !is.finite(boot$ci_low) || !is.finite(boot$ci_high)) {
      IQR_rel <- NA_real_; CI_rel <- NA_real_; boot_score_rel <- NA_real_
    } else {
      IQR_rel <- boot$iqr / expr_iqr
      CI_rel  <- (boot$ci_high - boot$ci_low) / expr_iqr
      boot_score_rel <- exp(-k1 * IQR_rel) * exp(-k2 * CI_rel)
    }
    
    # post-selection cutoff-position, group-balance, and TN-concordance metrics
    cutoff_metrics <- analyze_cutoff_position_v2(
      expr_tumor = df$expr,
      cutoff = best_cut,
      tumor_values = tumor_expr,
      normal_values = normal_expr,
      HR = hr
    )
    
    # hazard direction consistency across bootstrap cutoffs
    hazard_dir_consistency <- compute_hazard_dir_consistency_v2(
      df = df,
      bootstrap_cutoffs = boot$boot_cuts,
      hr_ref = hr
    )
    
    # optionally generate and save the single-gene summary panels
    if (save_plots) {
      pan <-plot_gene_panel_main_v2(
        df = df, gene = gene,
        best_cut = best_cut,
        hr = hr, ci_low = ci_low, ci_high = ci_high,
        cutoff_p = cutoff_p,boot = boot,
        scan_df = scan_res$scan_df,
        expected_dir = expected_dir,
        plot_dir = plot_dir,
        width       = 12,
        height      = 10,
        xps_aspect  = 0.625,
        km_aspect   = 1)
      
      panel_res <- plot_gene_panel_with_spline_v2(
        pan = pan,df = df, gene = gene,
        xps_aspect = 0.625,
        plot_dir = plot_dir,
        width       = 15,
        height      = 18)
      
      p_spline <- panel_res$p_spline
    } else {
      p_spline <- evaluate_spline_support_v2(df = df)$p_spline
    }
    
    results[[gene]] <- data.frame(
      Gene = gene,
      HR = hr, CI_low = ci_low, CI_high = ci_high, P = p_val,
      Cutoff = best_cut, Cutoff_P = cutoff_p,
      p_spline = p_spline,
      Bootstrap_Median = boot$median,
      Bootstrap_IQR = boot$iqr,
      Bootstrap_CI_low = boot$ci_low,
      Bootstrap_CI_high = boot$ci_high,
      Bootstrap_CI_range = boot$ci_high - boot$ci_low,
      IQR_rel = IQR_rel, CI_rel = CI_rel,
      bootstrap_score = boot_score_rel,
      expression_range = cutoff_metrics$expression_range,
      quantile_rank = cutoff_metrics$quantile_rank,
      density_score = cutoff_metrics$density_score,
      min_group_prop = cutoff_metrics$min_group_prop,
      low_group_n = cutoff_metrics$low_group_n,
      high_group_n = cutoff_metrics$high_group_n,
      TN_log2FC = cutoff_metrics$TN_log2FC,
      direction_concordance = cutoff_metrics$direction_concordance,
      hazard_dir_consistency = hazard_dir_consistency,
      stringsAsFactors = FALSE
    )
  }
  
  # ---------- 4) Assemble results, adjust P values, and compute composite ranking ----------
  results_df <- do.call(rbind, Filter(is.data.frame, results))
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No genes produced valid cutoffs. Consider relaxing min_events_per_group/minprop.")
    return(list(
      results_df = results_df,
      df_cache = df_cache,
      boot_cache = boot_cache,
      scan_cache = scan_cache,
      exprset_full = exprset_full,
      exprset_tumor = exprset_tumor
    ))
  }
  
  # multiple testing correction
  results_df$P_adj <- stats::p.adjust(as.numeric(results_df$P), method = adjust_method)
  results_df$Cutoff_P_adj <- stats::p.adjust(as.numeric(results_df$Cutoff_P), method = adjust_method)
  
  # scoring components (kept identical to your current implementation)
  results_df$bootstrap_score <- as.numeric(results_df$bootstrap_score)
  results_df$hazard_score <- as.numeric(results_df$hazard_dir_consistency)
  results_df$tn_score <- ifelse(results_df$direction_concordance == 1, 1,
                                ifelse(results_df$direction_concordance == -1, 0, NA_real_))
  results_df$density_score <- as.numeric(results_df$density_score)
  results_df$balance_score <- pmin(1, pmax(0, as.numeric(results_df$min_group_prop) * 2))
  
  w_boot <- 0.40; w_den <- 0.20; w_bal <- 0.10; w_haz <- 0.15; w_tn <- 0.15
  results_df$composite_score <- with(results_df,
                                     w_boot * bootstrap_score +
                                       w_haz  * hazard_score +
                                       w_tn   * tn_score +
                                       w_den  * density_score +
                                       w_bal  * balance_score)
  
  results_df$Rank <- rank(-results_df$composite_score, ties.method = "min")
  results_df <- results_df[order(results_df$Rank, results_df$P_adj), ]
  
  # ---------- 5) Optional summary plotting ----------
  if (save_plots && nrow(results_df) > 0) {
    
    plot_cutoff_composite_summary_v2(results_df,out_pdf = file.path(plot_dir, "test_composite_score_plot.pdf"),
                                     top_n = NULL,order_by ="composite_score")
    
    top_genes <- results_df[results_df$composite_score >= score_threshold &
                              is.finite(results_df$composite_score), , drop = FALSE]
    if (nrow(top_genes) > 0) {
      dir.create(file.path(plot_dir, "top_only"), showWarnings = FALSE, recursive = TRUE)
        
      plot_cutoff_composite_summary_v2(results_df = top_genes,out_pdf = file.path(plot_dir, "top_only", "test_top_genes_composite_score_plot.pdf"),
                                       top_n = NULL,order_by = "composite_score")
      }

  }
  
  # ---------- Return ----------
  list(
    results_df = results_df,
    df_cache = df_cache,
    boot_cache = boot_cache,
    scan_cache = scan_cache,
    exprset_full = exprset_full,
    exprset_tumor = exprset_tumor
  )
}