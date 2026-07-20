############################################################
# StaBiCut v2 — Single-gene plotting modules
# File: R/modules_plot_single_StaBiCut_v2.R
#
# Scope:
#   Plot constructors and single-gene panel assemblers used by the main
#   StaBiCut runner. These functions visualize the evidence supporting a
#   selected cutoff, but they do not modify analytical results.
#
# Included in this file:
#   - tumor-expression density plot
#   - bootstrap cutoff-density plot
#   - spline-effect plot
#   - Kaplan–Meier plot
#   - single-gene four-panel and five-panel assemblers
#   - composite summary plot for one run
#   - cut-scan diagnostic plot
#
# Dependencies assumed to be already sourced:
#   - modules_core_StaBiCut_v2.R
#
# Notes:
#   - patchwork, ggplot2, survminer, and survival must be installed.
#   - This script does not call library() globally.
############################################################

# ==========================================================
# 7. Single-gene plotting utilities
# ==========================================================

#' Tumor-expression density plot with the selected cutoff
#'
#' @param expr_vec Numeric expression vector used for density estimation.
#' @param gene Character scalar giving the gene label used in the plot title and
#'   x-axis label.
#' @param cutoff Numeric scalar giving the selected cutoff to be shown as a
#'   vertical reference line.
#' @return A ggplot object.
#' @export
build_expr_density_plot_v2 <- function(expr_vec, gene, cutoff) {
  expr_df <- data.frame(expr = as.numeric(expr_vec))
  cutoff <- as.numeric(cutoff)[1]

  ggplot2::ggplot(expr_df, ggplot2::aes(x = expr)) +
    ggplot2::geom_density(fill = "#00c2d1", alpha = 0.25) +
    ggplot2::geom_vline(xintercept = cutoff, color = "#ea4846", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = paste0("Expression distribution of ", gene),
      subtitle = paste0("Cutoff = ", signif(cutoff, 3)),
      x = paste(gene, "expression (log2 TPM)"),
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}

#' Density-based summary of bootstrap cutoff structure
#'
#' This helper provides auxiliary, density-based diagnostics of the retained
#' bootstrap cutoff distribution. These summaries are not used in the current
#' StaBiCut v2 stability score or composite ranking, but are retained for
#' diagnostic plotting, structural inspection, and potential future extensions.
#'
#' @param bootstrap_cutoffs Numeric vector of retained bootstrap cutoffs.
#' @param min_n Integer scalar giving the minimum number of cutoffs required for
#'   structure evaluation.
#' @param bw_mult Numeric scalar giving the multiplier applied to the
#'   `stats::density()` bandwidth estimate when defining the neighborhood around
#'   the dominant mode.
#'
#' @return A named list containing auxiliary diagnostics of bootstrap-cutoff
#'   structure: the number of valid cutoffs, bootstrap IQR, density-peak
#'   contrast, dominant-mode mass, dominant-mode location, and a status message.
#' @export
evaluate_bootstrap_structure_v2 <- function(bootstrap_cutoffs, min_n = 30, bw_mult = 1.0) {
  x <- bootstrap_cutoffs[is.finite(bootstrap_cutoffs)]
  if (length(x) < min_n) {
    return(list(
      n = length(x),
      boot_iqr = NA_real_,
      peak_contrast = NA_real_,
      main_mode_mass = NA_real_,
      main_mode = NA_real_,
      reason = "Too few valid bootstrap cutoffs"
    ))
  }

  dens <- stats::density(x)
  y <- dens$y

  peaks <- which(diff(sign(diff(y))) == -2) + 1
  if (length(peaks) == 0) {
    return(list(
      n = length(x),
      boot_iqr = stats::IQR(x),
      peak_contrast = NA_real_,
      main_mode_mass = NA_real_,
      main_mode = NA_real_,
      reason = "No detectable density peak"
    ))
  }

  ord <- order(y[peaks], decreasing = TRUE)
  main_peak <- peaks[ord[1]]
  second_height <- if (length(ord) >= 2) y[peaks[ord[2]]] else 0
  peak_contrast <- second_height / y[main_peak]

  bw <- dens$bw * bw_mult
  main_mode <- dens$x[main_peak]
  main_mode_mass <- mean(abs(x - main_mode) <= bw)

  list(
    n = length(x),
    boot_iqr = stats::IQR(x),
    peak_contrast = as.numeric(peak_contrast),
    main_mode_mass = as.numeric(main_mode_mass),
    main_mode = as.numeric(main_mode),
    reason = "OK"
  )
}

#' Visualize the retained bootstrap-cutoff distribution
#'
#' This plotting helper displays the density of retained bootstrap cutoffs and
#' annotates the figure with externally computed stability summaries.
#' It does not compute bootstrap-structure diagnostics itself.
#'
#' @param bootstrap_cutoffs Numeric vector of retained bootstrap cutoffs.
#' @param cutoff Numeric scalar giving the reference cutoff shown as a vertical
#'   dashed line.
#' @param gene Character scalar giving the gene label used in the plot title.
#' @param IQR_rel Numeric scalar giving the relative bootstrap IQR to be shown
#'   in the plot subtitle.
#' @param CI_rel Numeric scalar giving the relative bootstrap interval width to
#'   be shown in the plot subtitle.
#' @param bootstrap_score Numeric scalar giving the bootstrap stability score to
#'   be shown in the plot subtitle.
#'
#' @return A ggplot object. If fewer than 10 finite bootstrap cutoffs are
#'   available, an empty placeholder plot with a diagnostic subtitle is
#'   returned.
#' @export
build_bootstrap_density_plot_clean_v2 <- function(
    bootstrap_cutoffs,
    cutoff,
    gene,
    IQR_rel = NA_real_,
    CI_rel = NA_real_,
    bootstrap_score = NA_real_
) {
  x <- bootstrap_cutoffs[is.finite(bootstrap_cutoffs)]
  if (length(x) < 10) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(
          title = paste("Bootstrap cutoff distribution –", gene),
          subtitle = "Too few valid bootstrap cutoffs"
        )
    )
  }

  df <- data.frame(cutoff = x)

  subtitle_txt <- paste0(
    "IQR_rel = ", ifelse(is.finite(IQR_rel), sprintf("%.3f", IQR_rel), "NA"),
    "; CI_rel = ", ifelse(is.finite(CI_rel), sprintf("%.3f", CI_rel), "NA"),
    "; bootstrap score = ", ifelse(is.finite(bootstrap_score), sprintf("%.3f", bootstrap_score), "NA")
  )

  ggplot2::ggplot(df, ggplot2::aes(x = cutoff)) +
    ggplot2::geom_density(na.rm = TRUE, fill = "#69b3a2", alpha = 0.40) +
    ggplot2::geom_vline(xintercept = cutoff, color = "#F46D43", linetype = "dashed", linewidth = 1.0) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste("Bootstrap cutoff distribution –", gene),
      subtitle = subtitle_txt,
      x = "Cutoff",
      y = "Density"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

#' Build a spline-effect plot from a Cox model with natural splines
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   `expr` is modeled using `splines::ns(expr, df = 3)` in the Cox model.
#' @param gene Character scalar giving the gene label used in plot annotation.
#'
#' @return A named list with components `p_spline` and `plot`, where `p_spline`
#'   is the likelihood-ratio P value for spline support and `plot` is the
#'   spline-effect ggplot. If spline fitting or term extraction fails, a
#'   placeholder plot is returned.
#' @export
build_spline_plot_v2 <- function(df, gene) {
  p_spline <- NA_real_

  fit_null <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ 1, data = df),
    error = function(e) NULL
  )
  fit_spl <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ splines::ns(expr, df = 3), data = df),
    error = function(e) NULL
  )

  if (!is.null(fit_null) && !is.null(fit_spl)) {
    an <- tryCatch(stats::anova(fit_null, fit_spl), error = function(e) NULL)
    if (!is.null(an) && nrow(an) >= 2 && all(c("Chisq", "Df") %in% colnames(an))) {
      chisq_val <- an$Chisq[2]
      df_val <- an$Df[2]
      if (is.finite(chisq_val) && is.finite(df_val) && df_val > 0) {
        p_spline <- stats::pchisq(chisq_val, df = df_val, lower.tail = FALSE)
      }
    }
  }

  if (is.null(fit_spl)) {
    return(list(
      p_spline = p_spline,
      plot = ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::labs(title = paste("Spline effect of", gene), subtitle = "Spline fit failed")
    ))
  }

  term <- suppressWarnings(
    tryCatch(stats::termplot(fit_spl, se = TRUE, plot = FALSE), error = function(e) NULL)
  )
  if (is.null(term) || !"expr" %in% names(term)) {
    return(list(
      p_spline = p_spline,
      plot = ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::labs(title = paste("Spline effect of", gene), subtitle = "termplot failed")
    ))
  }

  spline_df <- data.frame(
    expr = term$expr$x,
    fit = term$expr$y,
    upper = term$expr$y + 2 * term$expr$se,
    lower = term$expr$y - 2 * term$expr$se
  )

  xr <- range(spline_df$expr, na.rm = TRUE)
  yr <- range(c(spline_df$lower, spline_df$upper), na.rm = TRUE)
  x_annot <- xr[1] + 0.03 * diff(xr)
  y_annot <- yr[2] - 0.05 * diff(yr)

  p <- ggplot2::ggplot(spline_df, ggplot2::aes(x = expr, y = fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "#acb1f1") +
    ggplot2::geom_line(linewidth = 1.0, color = "#2b6ed4") +
    ggplot2::geom_rug(data = df, ggplot2::aes(x = expr), inherit.aes = FALSE, sides = "b", alpha = 0.25) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::annotate(
      "text",
      x = x_annot,
      y = y_annot,
      label = paste0(
        "Spline effect of ", gene,
        "   (P = ", ifelse(is.finite(p_spline), signif(p_spline, 3), "NA"), ")"
      ),
      hjust = 0,
      vjust = 1,
      size = 4.8
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4))
    ) +
    ggplot2::labs(
      x = paste0(gene, " expression (log2 TPM)"),
      y = "Partial log HR"
    )

  list(p_spline = p_spline, plot = p)
}

#' Build a Kaplan-Meier plot at the selected cutoff
#'
#' @param df data.frame containing `group`, `time`, and `event`, where `group`
#'   is the cutoff-defined Low/High classification used for survival fitting.
#' @param gene Character scalar retained for interface consistency.
#' @param hr Numeric scalar giving the hazard ratio shown in the plot
#'   annotation.
#' @param ci_low Numeric scalar giving the lower 95% confidence bound shown in
#'   the plot annotation.
#' @param ci_high Numeric scalar giving the upper 95% confidence bound shown in
#'   the plot annotation.
#' @param pval Numeric scalar giving the displayed P value.
#' @param pal Named character vector giving the colors for `Low` and `High`.
#'
#' @return A ggplot object.
#' @export
build_survival_plot_v2 <- function(df, gene, hr, ci_low, ci_high, pval,
                                          pal = c("Low" = "#c1e153", "High" = "#f383b6")) {
  df$group <- factor(df$group, levels = c("Low", "High"))
  surv_fit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  p <- survminer::ggsurvplot(
    surv_fit, data = df,
    pval = FALSE, conf.int = TRUE,
    palette = c(pal["Low"], pal["High"]),
    legend.title = "Group", legend.labs = c("Low", "High"),
    ggtheme = ggplot2::theme_bw(base_size = 14) + ggplot2::theme(legend.position = "top"),
    censor.shape = 124, censor.size = 3)

  ci_text <- paste0(sprintf("%.3f", hr), " [", sprintf("%.3f", ci_low), ", ", sprintf("%.3f", ci_high), "]")
  xmin <- min(df$time, na.rm = TRUE)
  xmax <- max(df$time, na.rm = TRUE)
  x_pos <- xmin + 0.02 * (xmax - xmin)

  p$plot + ggplot2::annotate(
      "text", x = x_pos, y = 0.08,
      label = paste0("P = ", sprintf("%.3f", pval), "\nHR = ", ci_text),
      size = 4.5, color = "black", hjust = 0, vjust = 0, lineheight = 0.95) +
    ggplot2::xlab("Month") +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "on") +
    ggplot2::theme(text = ggplot2::element_text(size = 14)) +
    ggplot2::theme(panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
}

# Internal helper: apply consistent plot margins
.set_panel_margin_v2 <- function(p, t = 4, r = 4, b = 4, l = 4) {
  p + ggplot2::theme(plot.margin = ggplot2::margin(t, r, b, l))
}

# Internal helper: add a panel tag
.add_tag_v2 <- function(p, tag) {
  p +
    ggplot2::labs(tag = tag) +
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 15, face = "bold"),
      plot.tag.position = c(0.05, 1.00)
    )
}

# Internal helper: add a panel border
.set_panel_border_v2 <- function(p, lw = 0.6) {
  p + ggplot2::theme(
    panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = lw)
  )
}

#' Build and save the four-panel single-gene summary
#'
#' The returned object contains the four primary panels used in the StaBiCut v2
#' single-gene view:
#'   A) cutoff scan
#'   B) expression density with the selected cutoff
#'   C) bootstrap-cutoff density
#'   D) Kaplan-Meier plot at the selected cutoff
#'
#' @param df data.frame with columns `expr`, `time`, `event`, and `group`.
#'   `group` should represent the cutoff-defined Low/High classification.
#' @param gene Character scalar giving the gene label used in titles and the
#'   output filename.
#' @param best_cut Numeric scalar giving the selected cutoff.
#' @param hr,ci_low,ci_high,cutoff_p Numeric scalars used for Kaplan-Meier plot
#'   annotation.
#' @param boot Named list returned by `bootstrap_cutoffs_v2()`.
#' @param scan_df data.frame returned by `scan_cutpoints_v2()`.
#' @param expected_dir Character scalar or `NULL`, used to indicate direction-
#'   compatible candidates in the cutoff-scan panel.
#' @param plot_dir Character scalar giving the output directory for the saved
#'   PDF.
#' @param width,height Numeric scalars passed to `ggplot2::ggsave()`.
#' @param xps_aspect Numeric scalar giving the aspect ratio for the first three
#'   panels.
#' @param km_aspect Numeric scalar giving the aspect ratio for the Kaplan-Meier
#'   panel.
#'
#' @return An invisible named list containing the assembled panel (`panel`), the
#'   four component plots (`pA`, `pB`, `pC`, `pD`), and the saved PDF path
#'   (`out_pdf`).
#' @export
plot_gene_panel_main_v2 <- function(
    df, gene, best_cut, hr, ci_low, ci_high,
    cutoff_p, boot, scan_df, expected_dir = NULL,
    plot_dir, width = 12, height = 10,
    xps_aspect = 0.625,
    km_aspect = 1
) {
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  # A: cut scan (Wald ~ cut)
  if (!is.null(scan_df) && is.data.frame(scan_df) && nrow(scan_df) > 0) {
    scan_df$pass_dir <- TRUE
    if (!is.null(expected_dir)) {
      if (expected_dir == "protective_high") scan_df$pass_dir <- scan_df$hr < 1
      if (expected_dir == "adverse_high")    scan_df$pass_dir <- scan_df$hr > 1
    }

    pA <- ggplot2::ggplot(scan_df, ggplot2::aes(x = cut, y = wald)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(ggplot2::aes(alpha = pass_dir), size = 1.6) +
      ggplot2::geom_vline(xintercept = best_cut, linetype = "dashed", color = "#ea4846", linewidth = 1) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste0(gene, " cutoff scan (Wald ~ cut)"),
        subtitle = paste0(
          "best_cut = ", signif(best_cut, 4),
          if (!is.null(expected_dir)) paste0(" | dir = ", expected_dir) else ""
        ),
        x = "Cutoff",
        y = "Wald statistic"
      ) +
      ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15), guide = "none") +
      ggplot2::theme(aspect.ratio = xps_aspect) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 6, b = 2)),
        plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 4))
      )

    pA <- .set_panel_margin_v2(pA, 4, 14, 4, 4)
    pA <- .add_tag_v2(pA, "A")
  } else {
    pA <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0(gene, " cutoff scan"), subtitle = "scan_df unavailable")
    pA <- .set_panel_margin_v2(pA, 4, 4, 4, 4)
    pA <- .add_tag_v2(pA, "A")
  }

  # B: expression density
  pB <- build_expr_density_plot_v2(df$expr, gene, best_cut) +
    ggplot2::theme(aspect.ratio = xps_aspect, legend.position = "none") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 6, b = 2)),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 4))
    )
  pB <- .set_panel_margin_v2(pB, 4, 4, 4, 14)
  pB <- .add_tag_v2(pB, "B")

 # C: bootstrap density
  expr_iqr <- stats::IQR(df$expr, na.rm = TRUE)
  if (!is.finite(expr_iqr) || expr_iqr <= 0 ||
      !is.finite(boot$iqr) || !is.finite(boot$ci_low) || !is.finite(boot$ci_high)) {
    IQR_rel <- NA_real_
    CI_rel <- NA_real_
    boot_score_rel <- NA_real_
  } else {
      k1 <- 3
      k2 <- 2
      IQR_rel <- boot$iqr / expr_iqr
      CI_rel  <- (boot$ci_high - boot$ci_low) / expr_iqr
      boot_score_rel <- exp(-k1 * IQR_rel) * exp(-k2 * CI_rel)
  }

  pC <- build_bootstrap_density_plot_clean_v2(
    bootstrap_cutoffs = boot$boot_cuts,
    cutoff = best_cut, gene = gene,
    IQR_rel = IQR_rel, CI_rel = CI_rel,
    bootstrap_score = boot_score_rel) +
    ggplot2::theme(aspect.ratio = xps_aspect, legend.position = "none") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 6, b = 2)),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 4))
    )
  pC <- .set_panel_margin_v2(pC, 4, 14, 4, 4)
  pC <- .add_tag_v2(pC, "C")

  # D: Kaplan-Meier plot
  pD <- build_survival_plot_v2(
    df = df, gene = gene, hr = hr,
    ci_low = ci_low, ci_high = ci_high,
    pval = cutoff_p) +
    ggplot2::theme(aspect.ratio = km_aspect, legend.position = "top") +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "on") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(margin = ggplot2::margin(t = 6, b = 2)),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(b = 4))
    )
  pD <- .set_panel_margin_v2(pD, 4, 4, 4, 14)
  pD <- .add_tag_v2(pD, "D")
  pD <- .set_panel_border_v2(pD)

 # force consistent aspect / legend behavior
  pA <- pA + ggplot2::theme(aspect.ratio = xps_aspect, legend.position = "none")
  pB <- pB + ggplot2::theme(aspect.ratio = xps_aspect, legend.position = "none")
  pC <- pC + ggplot2::theme(aspect.ratio = xps_aspect, legend.position = "none")
  pD <- pD + ggplot2::theme(aspect.ratio = km_aspect, legend.position = "top")

 # Simple 2x2 layout
  panel <- (pA | pB) / (pC | pD)
  out_pdf <- file.path(plot_dir, paste0(gene, "_panel_main_v2.pdf"))
  ggplot2::ggsave(filename = out_pdf, plot = panel, width = width, height = height, limitsize = FALSE)

  invisible(list(panel = panel, pA = pA, pB = pB, pC = pC, pD = pD, out_pdf = out_pdf))
}

#' Extend the single-gene summary with a spline-effect panel and save a
#' five-panel PDF
#'
#' @param pan Named list returned by `plot_gene_panel_main_v2()`, containing `pA`, `pB`, `pC`, and `pD`.
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#' @param gene Character scalar giving the gene label used in the spline panel
#'   and output filename.
#' @param plot_dir Character scalar giving the output directory.
#' @param xps_aspect Numeric scalar giving the aspect ratio used for the spline
#'   panel.
#' @param width,height Numeric scalars passed to `ggplot2::ggsave()`.
#'
#' @return An invisible named list containing the assembled five-panel figure
#'   (`panel`), the spline P value (`p_spline`), and the saved PDF path
#'   (`out_pdf`).
#' @export
plot_gene_panel_with_spline_v2 <- function(
    pan, df, gene, plot_dir,
    xps_aspect = 0.625,
    width = 15, height = 18) {
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

 # pan from cutoff-stability plotting: including A: cut scan (Wald ~ cut), B: expression density, C: bootstrap density, D: Kaplan-Meier plot
 # E: spline effect
  spl <- build_spline_plot_v2(df, gene)
  pE <- spl$plot + ggplot2::theme(legend.position = "none", aspect.ratio = xps_aspect)
  pE <- .set_panel_margin_v2(pE, 2, 4, 4, 4)
  pE <- .add_tag_v2(pE, "E")
 
  # Simple 2x2 layout
  design <- "AB\nCD\nEE\n"
  panel <- patchwork::wrap_plots(A = pan$pA, B = pan$pB, C = pan$pC, D = pan$pD, E = pE, design = design)

  out_pdf <- file.path(plot_dir, paste0(gene, "_panel_with_spline.pdf"))
  ggplot2::ggsave(out_pdf, panel, width = width, height = height, limitsize = FALSE)

  invisible(list(panel = panel, p_spline = spl$p_spline, out_pdf = out_pdf))
}

#' Plot the within-run composite score summary across genes
#'
#' This function visualizes the per-gene composite score and its component
#' support values for a single StaBiCut run.
#'
#' @param results_df Per-gene results table containing composite and component
#'   scores.
#' @param out_pdf Character scalar giving the output PDF path.
#' @param top_n Optional integer giving the number of top genes to display.
#' @param order_by Character string specifying the ordering rule.
#'
#' @return A patchwork plot, returned invisibly.
#' @export
plot_cutoff_composite_summary_v2 <- function(
    results_df,
    out_pdf = "test_composite_score_plot.pdf",
    top_n = NULL,
    order_by = c("composite_score", "Rank")
) {
  order_by <- match.arg(order_by)
  stopifnot(is.data.frame(results_df))

  need_cols <- c("Gene", "bootstrap_score", "hazard_score", "tn_score",
                 "density_score", "balance_score", "composite_score")
  miss <- setdiff(need_cols, colnames(results_df))
  if (length(miss) > 0) {
    stop("results_df is missing required columns: ", paste(miss, collapse = ", "))
  }

        df <- results_df %>%
                dplyr::select(
                  Gene, composite_score,
                  bootstrap_score, hazard_score, tn_score, density_score, balance_score,
                  dplyr::any_of("Rank")
                ) %>%
                dplyr::mutate(
                  across(c(composite_score, bootstrap_score, hazard_score,
                           tn_score, density_score, balance_score), as.numeric)
                )
              
        df <- df %>%
                mutate(
                  tn_score = ifelse(is.na(tn_score), 0, tn_score),
                  composite_score = ifelse(is.na(composite_score), 0, composite_score)
                )

  if (order_by == "Rank" && "Rank" %in% colnames(df)) {
              df <- df %>% arrange(Rank, desc(composite_score))
              } else {
                df <- df %>% arrange(desc(composite_score))
              }

  if (!is.null(top_n)) {
                df <- df %>% slice_head(n = top_n)
              }

  gene_levels <- rev(df$Gene)
  df$Gene <- factor(df$Gene, levels = gene_levels)

  heat_df <- df %>%
    dplyr::select(Gene, bootstrap_score, hazard_score, tn_score, density_score, balance_score) %>%
                pivot_longer(
                  cols = -Gene,
                  names_to = "Dimension",
                  values_to = "Support")

  heat_df$Dimension <- factor(
    heat_df$Dimension,
    levels = c("bootstrap_score", "hazard_score", "tn_score", "density_score", "balance_score")
  )

 # Left: composite bar
  p_bar <- ggplot2::ggplot(df, ggplot2::aes(x = composite_score, y = Gene)) +
    ggplot2::geom_col(width = 0.72, fill = "grey35") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", composite_score)), hjust = -0.15, size = 3.8) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(1.02, max(df$composite_score, na.rm = TRUE) + 0.08)),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Composite score (0–1)", y = NULL, title = "Composite") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  #Right: five dimensional heatmap
  p_heat <- ggplot2::ggplot(heat_df, ggplot2::aes(x = Dimension, y = Gene, fill = Support)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Support)), size = 3.2) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 1), name = "Support (0–1)") +
    ggplot2::labs(x = "Dimensions", y = NULL, title = "StaBiCut composite stability score") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  p <- p_bar + p_heat + patchwork::plot_layout(widths = c(1.0, 1.5))
  ggplot2::ggsave(out_pdf, p, width = 11, height = max(4.5, 0.45 * nrow(df) + 1.8))
  invisible(p)
}

#' Plot a diagnostic cutoff-scan curve
#'
#' @param scan_df data.frame returned from `scan_cutpoints_v2()`.
#' @param gene Character scalar giving the gene label used in the plot title.
#' @param best_cut Numeric scalar giving the selected cutoff.
#' @param expected_dir Character scalar or `NULL`, used to indicate
#'   direction-compatible candidates in the plot.
#' @param out_pdf Character scalar giving the output PDF path.
#'
#' @return Invisibly returns the ggplot object, or `NULL` if `scan_df` is
#'   missing or empty.
#' @export
plot_cut_scan_curve_v2 <- function(scan_df, gene, best_cut, expected_dir = NULL, out_pdf) {
  if (is.null(scan_df) || nrow(scan_df) == 0) return(invisible(NULL))
  scan_df$pass_dir <- TRUE
  if (!is.null(expected_dir)) {
    if (expected_dir == "protective_high") scan_df$pass_dir <- scan_df$hr < 1
    if (expected_dir == "adverse_high")    scan_df$pass_dir <- scan_df$hr > 1
  }

  p <- ggplot2::ggplot(scan_df, ggplot2::aes(x = cut, y = wald)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(alpha = pass_dir), size = 1.6) +
    ggplot2::geom_vline(xintercept = best_cut, linetype = "dashed", color = "#ea4846", linewidth = 1) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0(gene, " cutoff scan (Wald ~ cut)"),
      subtitle = paste0("best_cut=", signif(best_cut, 4), if (!is.null(expected_dir)) paste0(" | dir=", expected_dir) else ""),
      x = "cutoff",
      y = "Wald statistic"
    ) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15), guide = "none")

  ggplot2::ggsave(out_pdf, p, width = 8, height = 4.8)
  invisible(p)
}

