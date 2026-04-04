############################################################
# StaBiCut v2 — Expected-direction tables and representative-seed utilities
# File: R/modules_seed_selection_StaBiCut_v2.R
#
# Scope:
#   Utilities for extracting the expected-direction layer actually used in the
#   analysis, exporting per-seed diagnostic plots, summarizing cross-seed mean
#   composite patterns, and selecting / loading representative seed-level runs
#   for figure assembly.
#
# Included in this file:
#   - expected-direction table builder
#   - per-seed cut-scan export helper
#   - cross-seed mean composite summary plot
#   - representative-seed selection
#   - loader for full representative-seed result objects
#   - Excel export for representative-seed ranking metrics
#
# Dependencies assumed to be already sourced:
#   - modules_core_StaBiCut_v2.R
#   - modules_plot_single_StaBiCut_v2.R
#
# Notes:
#   - This script does not call library() globally.
############################################################

# ==========================================================
# 9. Utilities for prior extraction and representative runs
# ==========================================================

#' Extract the expected-direction layer actually used for a gene set
#'
#' This helper mirrors the preprocessing steps used in the main runner so that
#' the exported direction table remains aligned with the analysis path actually
#' applied before cutoff scanning.
#'
#' @param exprset Numeric TPM matrix with genes in rows and samples in columns.
#' @param clin Clinical data frame.
#' @param geneset Character vector of target genes.
#' @param gene_prior_table Optional prior table containing gene-specific
#'   expected directions.
#' @param force_direction Logical; whether to use the stress-testing branch of
#'   direction inference.
#' @param winsor_lower Lower winsorization quantile.
#' @param winsor_upper Upper winsorization quantile.
#'
#' @return Data frame containing `Gene`, `expected_dir`, `TN_log2FC`, and
#'   `cox_beta`.
#' @export
get_expected_dir_table_v2 <- function(exprset,clin, geneset,gene_prior_table = NULL,
                                                  force_direction = FALSE,winsor_lower = 0.05,
                                                  winsor_upper = 0.95) {
              exprset <- log2(exprset + 1)
              
              exprset <- exprset[rowMeans(exprset) > 0.5, , drop = FALSE]
              exprset_full <- exprset
              
              is_tumor <- as.numeric(substr(colnames(exprset), 14, 15)) < 10
              exprset_tumor <- exprset[, is_tumor, drop = FALSE]
              clin_tumor <- clin[is_tumor, , drop = FALSE]
              
              clin_tumor$time <- ifelse(!is.na(clin_tumor$days_to_death),
                                        clin_tumor$days_to_death,
                                        clin_tumor$days_to_last_follow_up) / 30
              clin_tumor$event <- ifelse(clin_tumor$vital_status %in% c("Dead", "dead", "DEAD"), 1, 0)
              clin_tumor <- clin_tumor[!is.na(clin_tumor$time) & !is.na(clin_tumor$event), , drop = FALSE]
              
              sample_ids <- intersect(colnames(exprset_tumor), rownames(clin_tumor))
              exprset_tumor <- exprset_tumor[, sample_ids, drop = FALSE]
              clin_tumor <- clin_tumor[sample_ids, , drop = FALSE]
              
              out <- lapply(geneset, function(gene) {
                if (!gene %in% rownames(exprset_tumor) || !gene %in% rownames(exprset_full)) {
                  return(data.frame(
                    Gene = gene,
                    expected_dir = NA_character_,
                    TN_log2FC = NA_real_,
                    cox_beta = NA_real_,
                    stringsAsFactors = FALSE
                  ))
                }
                
                expr <- as.numeric(exprset_tumor[gene, ])
                df <- data.frame(
                  expr  = expr,
                  time  = as.numeric(clin_tumor$time),
                  event = as.numeric(clin_tumor$event)
                )
                
                df$expr <- winsorize(df$expr, lower = winsor_lower, upper = winsor_upper)
                
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
                
                TN_log2FC <- if (length(tumor_expr) > 0 && length(normal_expr) > 0) {
                  median(tumor_expr, na.rm = TRUE) - median(normal_expr, na.rm = TRUE)
                } else {
                  NA_real_
                }
                
                fit_cont <- tryCatch(
                  survival::coxph(survival::Surv(time, event) ~ expr, data = df),
                  error = function(e) NULL
                )
                beta <- if (!is.null(fit_cont)) as.numeric(stats::coef(fit_cont)[1]) else NA_real_
                
                data.frame(
                  Gene = gene,
                  expected_dir = if (is.null(expected_dir)) NA_character_ else expected_dir,
                  TN_log2FC = TN_log2FC,
                  cox_beta = beta,
                  stringsAsFactors = FALSE
                )
              })
              
              do.call(rbind, out)
            }

#' Export cutoff-scan PDFs for all genes from one seed-level result
#'
#' @param all_runs List of full seed-level StaBiCut result objects.
#' @param seed_index Integer index within `all_runs`.
#' @param geneset Character vector of genes to export.
#' @param expected_dir_table Data frame containing `Gene` and `expected_dir`.
#' @param out_dir Character scalar giving the output directory.
#'
#' @return A named list containing export metadata, returned invisibly. 
#' @export
export_one_seed_cutscan_pdfs_v2 <- function(
    all_runs,
    seed_index = 1,
    geneset,
    expected_dir_table,
    out_dir = "./one_seed_cutscan"
)  {
              if (is.null(all_runs) || length(all_runs) == 0) {
                stop("all_runs is empty or NULL.")
              }
              
              if (seed_index < 1 || seed_index > length(all_runs)) {
                stop("seed_index is out of range. length(all_runs) = ", length(all_runs))
              }
              
              if (is.null(expected_dir_table) || !is.data.frame(expected_dir_table)) {
                stop("expected_dir_table must be a data.frame.")
              }
              
              if (!all(c("Gene", "expected_dir") %in% colnames(expected_dir_table))) {
                stop("expected_dir_table must contain columns: Gene, expected_dir")
              }
              
              # extract results of the specific seed
              res_seed <- all_runs[[seed_index]]
              
              seed_id <- if (!is.null(res_seed$seed)) {
                res_seed$seed
              } else {
                seed_index
              }
              
              seed_out_dir <- file.path(out_dir, paste0("seed_", seed_id))
              if (!dir.exists(seed_out_dir)) dir.create(seed_out_dir, recursive = TRUE)
              
              export_log <- vector("list", length(geneset))
              
              for (i in seq_along(geneset)) {
                gene <- geneset[i]
                
                message("Processing seed ", seed_id, " | gene: ", gene)
                
                scan_df_i <- res_seed$scan_cache[[gene]]
                
                if (is.null(scan_df_i) || !is.data.frame(scan_df_i) || nrow(scan_df_i) == 0) {
                  export_log[[i]] <- data.frame(
                    seed = seed_id,
                    Gene = gene,
                    status = "skip_scan_df_missing",
                    best_cut = NA_real_,
                    expected_dir = NA_character_,
                    out_pdf = NA_character_,
                    stringsAsFactors = FALSE
                  )
                  next
                }
                
                ridx <- match(gene, res_seed$results_df$Gene)
                if (is.na(ridx)) {
                  export_log[[i]] <- data.frame(
                    seed = seed_id,
                    Gene = gene,
                    status = "skip_gene_not_in_results_df",
                    best_cut = NA_real_,
                    expected_dir = NA_character_,
                    out_pdf = NA_character_,
                    stringsAsFactors = FALSE
                  )
                  next
                }
                
                best_cut_i <- res_seed$results_df$Cutoff[ridx]
                
                # extract expected_dir only
                didx <- match(gene, expected_dir_table$Gene)
                expected_dir_i <- if (!is.na(didx)) expected_dir_table$expected_dir[didx] else NA_character_
                
                if (is.na(expected_dir_i) || length(expected_dir_i) == 0) {
                  expected_dir_i <- NULL
                }
                
                out_pdf_i <- file.path(seed_out_dir, paste0(gene, "_seed", seed_id, "_cut_scan_wald.pdf"))
                
                tryCatch({
                  plot_cut_scan_curve_v2(
                    scan_df = scan_df_i,
                    gene = gene,
                    best_cut = best_cut_i,
                    expected_dir = expected_dir_i,
                    out_pdf = out_pdf_i
                  )
                  
                  export_log[[i]] <- data.frame(
                    seed = seed_id,
                    Gene = gene,
                    status = "ok",
                    best_cut = best_cut_i,
                    expected_dir = if (is.null(expected_dir_i)) NA_character_ else expected_dir_i,
                    out_pdf = out_pdf_i,
                    stringsAsFactors = FALSE
                  )
                  
                }, error = function(e) {
                  export_log[[i]] <<- data.frame(
                    seed = seed_id,
                    Gene = gene,
                    status = paste0("error: ", conditionMessage(e)),
                    best_cut = best_cut_i,
                    expected_dir = if (is.null(expected_dir_i)) NA_character_ else expected_dir_i,
                    out_pdf = out_pdf_i,
                    stringsAsFactors = FALSE
                  )
                })
              }
              
              export_log_df <- do.call(rbind, export_log)
              
              write.csv(
                export_log_df,
                file = file.path(seed_out_dir, paste0("seed_", seed_id, "_cutscan_export_log.csv")),
                row.names = FALSE
              )
              
              return(invisible(list(
                seed_id = seed_id,
                out_dir = seed_out_dir,
                export_log = export_log_df
              )))
            }


#' Plot the mean composite summary across multiple seeds
#'
#' @param all_runs List of full seed-level StaBiCut result objects.
#' @param out_pdf Character scalar giving the output PDF path.
#' @param top_n Optional integer giving the number of genes to display.
#' @param order_by Character string specifying the ordering rule.
#' @param p_col Character scalar giving the P-value column to summarize.
#' @param p_agg Character string specifying the aggregation rule for P values.
#' @param use_adjusted_p Logical; whether adjusted P values should be used for
#'   significance annotation.
#' @param star_size Numeric scalar giving the size of significance stars.
#' @param bar_fill Fill color for the mean composite-score bars.
#'
#' @return A named list containing the assembled plot and its main
#'   summary tables, retuned invisibly. 
#' @export
plot_multiseed_composite_summary_v2 <- function(
    all_runs,
    out_pdf = "test_multiseed_composite_mean_summary.pdf",
    top_n = NULL,
    order_by = c("composite_mean", "Rank_mean"),
    p_col = c("Cutoff_P", "P"),
    p_agg = c("median", "mean"),
    use_adjusted_p = FALSE,
    star_size = 5.2,
    bar_fill = "grey35"
) {
  order_by <- match.arg(order_by)
  p_col <- match.arg(p_col)
  p_agg <- match.arg(p_agg)

  df_list <- lapply(seq_along(all_runs), function(i) {
    x <- all_runs[[i]]
    if (is.null(x$results_df) || !is.data.frame(x$results_df)) return(NULL)
    df <- x$results_df
    df$seed <- if (!is.null(x$seed)) x$seed else i
    df
  })
  df_all <- dplyr::bind_rows(df_list)
  if (nrow(df_all) == 0) stop("No valid results_df found inside all_runs.")

  need_cols <- c("Gene", "composite_score", "bootstrap_score", "hazard_score", "tn_score", "density_score", "balance_score")
  miss <- setdiff(need_cols, colnames(df_all))
  if (length(miss) > 0) stop("df_all is missing required columns: ", paste(miss, collapse = ", "))
  if (!p_col %in% colnames(df_all)) stop("Requested significance column not found: ", p_col)

  p_use_col <- p_col
  if (isTRUE(use_adjusted_p)) {
    p_adj_col <- paste0(p_col, "_adj")
    if (p_adj_col %in% colnames(df_all)) p_use_col <- p_adj_col
  }

  df_all <- df_all %>%
                mutate(
                  Gene = as.character(Gene),
                  seed = as.integer(seed),
                  composite_score = as.numeric(composite_score),
                  bootstrap_score = as.numeric(bootstrap_score),
                  hazard_score = as.numeric(hazard_score),
                  tn_score = as.numeric(tn_score),
                  density_score = as.numeric(density_score),
                  balance_score = as.numeric(balance_score),
                  p_use = as.numeric(.data[[p_use_col]])
                ) %>%
                mutate(
                  tn_score = ifelse(is.na(tn_score), 0, tn_score)
                )

  p_to_star <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  }

 # summarize
  score_summary <- df_all %>%
                group_by(Gene) %>%
                summarise(
                  n_seeds = n_distinct(seed),
                  composite_mean = mean(composite_score, na.rm = TRUE),
                  composite_sd   = sd(composite_score, na.rm = TRUE),
                  
                  bootstrap_mean = mean(bootstrap_score, na.rm = TRUE),
                  hazard_mean    = mean(hazard_score, na.rm = TRUE),
                  tn_mean        = mean(tn_score, na.rm = TRUE),
                  density_mean   = mean(density_score, na.rm = TRUE),
                  balance_mean   = mean(balance_score, na.rm = TRUE),
                  
                  p_summary = if (p_agg == "median") median(p_use, na.rm = TRUE) else mean(p_use, na.rm = TRUE),
                  sig_freq_0.05 = mean(p_use < 0.05, na.rm = TRUE),
                  Rank_mean = mean(Rank, na.rm = TRUE),
                  .groups = "drop"
                ) %>%
                mutate(
                  sig_star = p_to_star(p_summary)
                )

   # order
  if (order_by == "Rank_mean" && any(is.finite(score_summary$Rank_mean))) {
                score_summary <- score_summary %>% arrange(Rank_mean, desc(composite_mean))
              } else {
                score_summary <- score_summary %>% arrange(desc(composite_mean))
              }
              
              if (!is.null(top_n)) {
                score_summary <- score_summary %>% slice_head(n = top_n)
              }
              
              gene_levels <- rev(score_summary$Gene)
              score_summary$Gene <- factor(score_summary$Gene, levels = gene_levels)

  # heatmap 
  heat_df <- score_summary %>%
                select(Gene, bootstrap_mean, density_mean, hazard_mean, tn_mean, balance_mean) %>%
                pivot_longer(
                  cols = -Gene,
                  names_to = "Dimension",
                  values_to = "Support"
                )
              
              heat_df$Dimension <- factor(
                heat_df$Dimension,
                levels = c("bootstrap_mean", "density_mean", "hazard_mean", "tn_mean", "balance_mean"),
                labels = c(
                  "Bootstrap score\n(0.40)",
                  "Cutoff plausibility score\n(0.20)",
                  "Hazard score\n(0.15)",
                  "TN score\n(0.15)",
                  "Balance score\n(0.10)"
                )
              )
              
  # bar
  xmax_data <- max(
                score_summary$composite_mean + ifelse(is.na(score_summary$composite_sd), 0, score_summary$composite_sd),
                na.rm = TRUE
              )
              
              xmax <- max(1.08, xmax_data + 0.12)
              
              score_summary <- score_summary %>%
                mutate(
                  label_x = pmin(composite_mean + composite_sd + 0.025, xmax - 0.14),
                  star_x  = xmax - 0.17
                )
              
              p_bar <- ggplot(score_summary, aes(x = composite_mean, y = Gene)) +
                geom_col(width = 0.72, fill = bar_fill) +
                geom_errorbar(
                  aes(
                    xmin = pmax(composite_mean - composite_sd, 0),
                    xmax = pmin(composite_mean + composite_sd, xmax)
                  ),
                  width = 0.18,
                  linewidth = 0.55
                ) +
                geom_text(aes(x = label_x, label = sprintf("%.2f", composite_mean)), hjust = 0, size = 3.8) +
                geom_text(aes(x = star_x, label = sig_star), hjust = 0, vjust = 0.5, size = star_size, fontface = "bold") +
                scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
                labs(x = "Mean composite score (0–1)", y = NULL, title = "Composite") +
                theme_bw(base_size = 12) +
                theme(
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  axis.text.y = element_text(face = "bold"),
                  axis.title.x = element_text(margin = margin(t = -48)),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5)
                )

  p_heat <- ggplot(heat_df, aes(x = Dimension, y = Gene, fill = Support)) +
                geom_tile(color = "white", linewidth = 0.8) +
                geom_text(aes(label = sprintf("%.2f", Support)), size = 3.2) +
                scale_fill_gradient(low = "#FFF7F3",
                                    high =  "#980043" ,
                                    limits = c(0, 1),name = "Support (0–1)") +
                labs(x = "Dimensions", y = NULL, title = "StaBiCut mean composite stability score") +
                theme_bw(base_size = 12) +
                theme(
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid = element_blank()
                )
              
              p <- p_bar + p_heat + patchwork::plot_layout(widths = c(1.05, 1.55))
              
              ggsave(
                out_pdf,
                p,
                width = 11,
                height = max(4.5, 0.45 * nrow(score_summary) + 1.8)
              )
              
              invisible(list(
                plot = p,
                summary_table = score_summary,
                heat_df = heat_df,
                df_all = df_all
              ))
}

#' Select a representative seed for visualization
#'
#' @param all_runs_slim List of per-seed slim result tables.
#' @param gene_x Character scalar giving the anchor gene used for
#'   representative-seed selection.
#' @param rank_weight Numeric weight for rank closeness.
#' @param score_weight Numeric weight for score closeness.
#' @param global_weight Numeric weight for global consistency.
#' @param cor_method Character string specifying the correlation method.
#' @param verbose Logical; if TRUE, print diagnostic messages.
#'
#' @return A named list containing the selected seed and the
#'   associated summary tables used for representative-seed ranking, returned invisibly.
#' @export
select_representative_seed_v2 <- function(
    all_runs_slim,
    gene_x = "ZG16",
    rank_weight = 0.35,
    score_weight = 0.45,
    global_weight = 0.20,
    cor_method = "spearman",
    verbose = TRUE
) {
  df_all <- dplyr::bind_rows(all_runs_slim)

  required_cols <- c("Gene", "seed", "composite_score")
  miss <- setdiff(required_cols, colnames(df_all))
  if (length(miss) > 0) stop("all_runs_slim is missing required columns: ", paste(miss, collapse = ", "))

  df_all <- df_all %>%
                mutate(
                  Gene = as.character(Gene),
                  seed = as.integer(seed),
                  composite_score = as.numeric(composite_score)
                )

  if (!gene_x %in% df_all$Gene) {
                stop("gene_x = '", gene_x, "' is not in all_runs_slim.")
              }

              # across-seed mean summary
              mean_tbl <- df_all %>%
                group_by(Gene) %>%
                summarise(
                  composite_mean = mean(composite_score, na.rm = TRUE),
                  composite_sd   = sd(composite_score, na.rm = TRUE),
                  .groups = "drop"
                ) %>%
                arrange(desc(composite_mean)) %>%
                mutate(mean_rank = row_number())
              
              gene_mean <- mean_tbl %>%
                filter(Gene == gene_x) %>%
                pull(composite_mean)
              
              gene_mean_rank <- mean_tbl %>%
                filter(Gene == gene_x) %>%
                pull(mean_rank)
              
              mean_rank_vec <- mean_tbl$mean_rank
              names(mean_rank_vec) <- mean_tbl$Gene
              
              mean_score_vec <- mean_tbl$composite_mean
              names(mean_score_vec) <- mean_tbl$Gene
              
              # per-seed ranking
              rank_tbl <- df_all %>%
                group_by(seed) %>%
                arrange(desc(composite_score), .by_group = TRUE) %>%
                mutate(rank_in_seed = row_number()) %>%
                ungroup()
              
              gene_tbl <- rank_tbl %>%
                filter(Gene == gene_x) %>%
                transmute(
                  seed,
                  Gene = Gene,
                  composite_score,
                  rank_in_seed,
                  distance_to_mean = abs(composite_score - gene_mean)
                )
              
              # modal rank of GeneX
              rank_mode_tbl <- gene_tbl %>%
                count(rank_in_seed, sort = TRUE)
              
              modal_rank <- rank_mode_tbl$rank_in_seed[1]
              
              gene_tbl <- gene_tbl %>%
                mutate(
                  distance_rank_to_mode = abs(rank_in_seed - modal_rank),
                  distance_rank_to_meanrank = abs(rank_in_seed - gene_mean_rank)
                )
              
              # global consistency of each seed
              #  correlation with mean composite score
              #  correlation with mean ranking
              seed_list <- split(rank_tbl, rank_tbl$seed)
              
              global_tbl <- lapply(seed_list, function(df_seed) {
                seed_id <- unique(df_seed$seed)
                
                seed_scores <- df_seed$composite_score
                names(seed_scores) <- df_seed$Gene
                
                seed_ranks <- df_seed$rank_in_seed
                names(seed_ranks) <- df_seed$Gene
                
                common_genes <- intersect(names(mean_score_vec), names(seed_scores))
                
                score_cor <- suppressWarnings(
                  cor(seed_scores[common_genes],
                      mean_score_vec[common_genes],
                      method = cor_method,
                      use = "pairwise.complete.obs")
                )
                
                rank_cor <- suppressWarnings(
                  cor(-seed_ranks[common_genes], -mean_rank_vec[common_genes],
                         method = cor_method,
                         use = "pairwise.complete.obs")
                )
                
                mean_abs_rank_diff <- mean(
                  abs(seed_ranks[common_genes] - mean_rank_vec[common_genes]),
                  na.rm = TRUE
                )
                
                data.frame(
                  seed = seed_id,
                  score_cor_with_mean = score_cor,
                  rank_cor_with_mean  = rank_cor,
                  mean_abs_rank_diff  = mean_abs_rank_diff,
                  stringsAsFactors = FALSE
                )
              }) %>%
                bind_rows()
              
              # merge and compute representative score
              out_tbl <- gene_tbl %>%
                left_join(global_tbl, by = "seed")
              
              # normalize to 0~1 for further weighting
              .rescale01 <- function(x, reverse = FALSE) {
                if (all(is.na(x))) return(rep(NA_real_, length(x)))
                rng <- range(x, na.rm = TRUE)
                if (diff(rng) == 0) {
                  z <- rep(0, length(x))
                } else {
                  z <- (x - rng[1]) / diff(rng)
                }
                if (reverse) z <- 1 - z
                z
              }
              
              out_tbl <- out_tbl %>%
                mutate(
                  score_closeness_norm = .rescale01(distance_to_mean, reverse = TRUE),
                  rank_closeness_norm  = .rescale01(distance_rank_to_mode, reverse = TRUE),
                  global_consistency_norm = rowMeans(
                    cbind(
                      .rescale01(score_cor_with_mean, reverse = FALSE),
                      .rescale01(rank_cor_with_mean, reverse = FALSE),
                      .rescale01(mean_abs_rank_diff, reverse = TRUE)
                    ),
                    na.rm = TRUE
                  )
                ) %>%
                mutate(
                  representative_score =
                    score_weight  * score_closeness_norm +
                    rank_weight   * rank_closeness_norm +
                    global_weight * global_consistency_norm
                ) %>%
                arrange(desc(representative_score), distance_to_mean, distance_rank_to_mode)
              
              best_seed <- out_tbl$seed[1]
              
              # output
              if (isTRUE(verbose)) {
                cat("\n==============================\n")
                cat("Representative seed selection\n")
                cat("==============================\n")
                cat("GeneX                :", gene_x, "\n")
                cat("Across-seed mean     :", round(gene_mean, 4), "\n")
                cat("Across-seed mean rank:", gene_mean_rank, "\n")
                cat("Modal rank of GeneX  :", modal_rank, "\n")
                cat("Selected seed        :", best_seed, "\n\n")
                
                cat("Top candidate seeds:\n")
                print(
                  out_tbl %>%
                    select(
                      seed, composite_score, distance_to_mean,
                      rank_in_seed, distance_rank_to_mode,
                      score_cor_with_mean, rank_cor_with_mean,
                      mean_abs_rank_diff, representative_score
                    ) %>%
                    head(10),
                  row.names = FALSE
                )
                cat("\n")
              }
              
              invisible(list(
                representative_seed = best_seed,
                gene_x = gene_x,
                gene_mean = gene_mean,
                gene_mean_rank = gene_mean_rank,
                modal_rank = modal_rank,
                seed_summary = out_tbl,
                mean_summary = mean_tbl
              ))}

#' Load one representative seed as a full result object
#' 
#' @param representative_seed Integer seed identifier.
#' @param gene_x Character scalar giving the focal gene.
#' @param rds_dir Character scalar giving the directory containing per-seed RDS
#'   files.
#' @param file_prefix Character scalar giving the filename prefix for per-seed
#'   RDS files.
#'
#' @return Named list containing the full representative-seed result object and
#'   gene-level extracts needed for downstream panel assembly.
#' @export
load_rep_seed_full_result_v2 <- function(
    representative_seed,
    gene_x = "ZG16",
    rds_dir = "./stability_runs",
    file_prefix = "results_seed_"
) {
  rds_file <- file.path(rds_dir, paste0(file_prefix, representative_seed, ".rds"))
  if (!file.exists(rds_file)) stop("RDS file not found: ", rds_file)

  res <- readRDS(rds_file)
  if (is.data.frame(res)) {
    warning("The loaded RDS appears to contain only results_df rather than a full result object.")
    return(list(
      representative_seed = representative_seed,
      gene_x = gene_x,
      rds_file = rds_file,
      full_result = res,
      gene_result_row = res %>% filter(Gene == gene_x)
      ))
  }
  if (is.null(res$results_df)) stop("The RDS file does not contain results_df.")

   gene_row <- res$results_df %>%
                filter(Gene == gene_x)

  if (nrow(gene_row) == 0) stop("Gene ", gene_x, " is not present in results_df for this seed.")

  out <- list(
                representative_seed = representative_seed,
                gene_x = gene_x,
                rds_file = rds_file,
                full_result = res,
                gene_result_row = gene_row,
                
                scan_df = if (!is.null(res$scan_cache) && !is.null(res$scan_cache[[gene_x]])) res$scan_cache[[gene_x]] else NULL,
                boot_vec = if (!is.null(res$boot_cache) && !is.null(res$boot_cache[[gene_x]])) res$boot_cache[[gene_x]] else NULL,
                gene_df = if (!is.null(res$df_cache)   && !is.null(res$df_cache[[gene_x]]))   res$df_cache[[gene_x]]   else NULL,
                
                best_cut = gene_row$Cutoff[1],
                expected_dir = gene_row$expected_dir[1],
                composite_score = gene_row$composite_score[1],
                rank = gene_row$Rank[1]
              )
              
              return(out)
}

#' Export representative-seed selection results to Excel
#'
#' @param all_runs_slim List of slim per-seed result tables.
#' @param gene_x Character scalar giving the focal gene.
#' @param out_file Character scalar giving the output workbook path.
#' @param overwrite Logical; if TRUE, overwrite an existing workbook.
#' @param rank_weight Numeric scalar giving the weight assigned to rank closeness
#'   when computing the integrated representative-seed score.
#' @param score_weight Numeric scalar giving the weight assigned to composite-score
#'   closeness to the across-seed mean when computing the integrated representative-seed score.
#' @param global_weight Numeric scalar giving the weight assigned to global ranking
#'   consistency when computing the integrated representative-seed score.
#'
#' @return A named list containing the selected seed and the output
#'   workbook path, returned invisibly.
#' @export
export_rep_seed_excel_v2 <- function(all_runs_slim,
                                                 gene_x = "ZG16",
                                                 out_file = "representative_seed_summary.xlsx",
                                                 overwrite = TRUE,
                                                 rank_weight = 0.35,
                                                 score_weight = 0.45,
                                                 global_weight = 0.20) {

              # select representative seed
              rep_res <- select_representative_seed_v2(
                all_runs_slim = all_runs_slim,
                gene_x = gene_x,
                rank_weight = rank_weight,
                score_weight = score_weight,
                global_weight = global_weight,
                verbose = FALSE
              )
              
              best_seed <- rep_res$representative_seed
              
              # selected seed results
              one_seed_df <- dplyr::bind_rows(all_runs_slim) %>%
                dplyr::filter(seed == best_seed)
              
              # representative_seed sheet
              representative_seed_tbl <- data.frame(
                Gene = gene_x,
                representative_seed = best_seed,
                gene_mean = rep_res$gene_mean,
                gene_mean_rank = rep_res$gene_mean_rank,
                modal_rank = rep_res$modal_rank,
                stringsAsFactors = FALSE
              )
              
              # Column dictionary
              column_dictionary <- data.frame(
                `Column name` = c(
                  # representative_seed
                  "Gene",
                  "representative_seed",
                  "gene_mean",
                  "gene_mean_rank",
                  "modal_rank",
                  
                  # seed_summary
                  "seed",
                  "composite_score",
                  "distance_to_mean",
                  "rank_in_seed",
                  "distance_rank_to_mode",
                  "distance_rank_to_meanrank",
                  "score_cor_with_mean",
                  "rank_cor_with_mean",
                  "mean_abs_rank_diff",
                  "score_closeness_norm",
                  "rank_closeness_norm",
                  "global_consistency_norm",
                  "representative_score",
                  
                  # mean_summary
                  "composite_mean",
                  "composite_sd"
                ),
                Category = c(
                  rep("representative_seed", 5),
                  rep("seed_summary", 13),
                  rep("mean_summary", 2)
                ),
                Definition = c(
                  # representative_seed
                  "Candidate gene used to identify the most representative seed for downstream figure display.",
                  "Recommended seed for figure-level visualization, selected computationally to best match the across-seed central tendency of the specified gene.",
                  "Across-seed mean composite score of the specified gene.",
                  "Rank of the specified gene when genes are ordered by across-seed mean composite score.",
                  "Most frequently observed within-seed rank of the specified gene across all evaluated seeds.",
                  
                  # seed_summary
                  "Random seed used for one full StaBiCut run.",
                  "Composite score of the specified gene in the corresponding seed-specific run.",
                  "Absolute distance between the specified gene's seed-specific composite score and its across-seed mean composite score.",
                  "Within-seed rank of the specified gene based on descending composite score.",
                  "Absolute distance between the specified gene's seed-specific rank and its modal rank across all seeds.",
                  "Absolute distance between the specified gene's seed-specific rank and its across-seed mean rank.",
                  "Correlation between the full gene-level composite-score vector in the current seed and the across-seed mean composite-score vector.",
                  "Correlation between the full gene-level ranking in the current seed and the across-seed mean ranking.",
                  "Mean absolute difference between seed-specific gene ranks and across-seed mean ranks, summarizing global ranking deviation.",
                  "Rescaled similarity score reflecting how closely the specified gene's seed-specific composite score matches its across-seed mean score (higher = closer).",
                  "Rescaled similarity score reflecting how closely the specified gene's seed-specific rank matches its modal rank (higher = closer).",
                  "Rescaled summary score reflecting how closely the overall gene ranking pattern in the current seed matches the across-seed mean ranking structure.",
                  "Integrated representativeness score used to select the recommended seed, combining gene-specific score closeness, rank closeness, and global ranking consistency.",
                  
                  # mean_summary
                  "Across-seed mean composite score for each gene.",
                  "Across-seed standard deviation of the composite score for each gene."
                ),
                stringsAsFactors = FALSE,
                check.names = FALSE
              )
              
              if (!"Column name" %in% colnames(column_dictionary) && "Column.name" %in% colnames(column_dictionary)) {
                colnames(column_dictionary)[colnames(column_dictionary) == "Column.name"] <- "Column name"
              }
              
              column_dictionary <- column_dictionary %>%
                dplyr::distinct(`Column name`, .keep_all = TRUE)
              
              # write workbook
              wb <- openxlsx::createWorkbook()
              
              openxlsx::addWorksheet(wb, "representative_seed")
              openxlsx::writeDataTable(wb, "representative_seed", representative_seed_tbl)
              
              openxlsx::addWorksheet(wb, "seed_summary")
              openxlsx::writeDataTable(wb, "seed_summary", rep_res$seed_summary)
              
              openxlsx::addWorksheet(wb, "mean_summary")
              openxlsx::writeDataTable(wb, "mean_summary", rep_res$mean_summary)
              
              openxlsx::addWorksheet(wb, "selected_seed_results")
              openxlsx::writeDataTable(wb, "selected_seed_results", one_seed_df)
              
              openxlsx::addWorksheet(wb, "Column_dictionary")
              openxlsx::writeDataTable(wb, "Column_dictionary", column_dictionary)
              
              openxlsx::saveWorkbook(wb, out_file, overwrite = overwrite)
              
              invisible(list(
                out_file = out_file,
                representative_seed = best_seed,
                representative_seed_tbl = representative_seed_tbl,
                seed_summary = rep_res$seed_summary,
                mean_summary = rep_res$mean_summary,
                selected_seed_results = one_seed_df,
                column_dictionary = column_dictionary
              ))
            }
            