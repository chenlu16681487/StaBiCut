############################################################
# StaBiCut v2 — Core analytical modules
# File: R/modules_core_StaBiCut_v2.R
#
# Scope:
#   Core utilities for preprocessing, direction-prior inference,
#   deterministic cutpoint scanning, fixed-cutoff Cox modeling,
#   direction-consistent local bootstrap re-support, and post-selection
#   cutoff-position / direction-consistency metrics.
#
# Included in this file:
#   - low-level helpers (.q, .clamp, .make_group)
#   - winsorization
#   - expected-direction inference
#   - deterministic cutpoint scan
#   - fixed-cutoff Cox model wrapper
#   - local bootstrap stability assessment
#   - post-selection plausibility / concordance metrics
#   - spline-support testing
#
# Intentionally NOT included:
#   - single-gene plotting
#   - multi-gene panel layout
#   - multi-seed summary plots / exports
#   - representative-seed selection utilities
#
# Notes:
#   - This script is designed to be sourced before run_StaBiCut_v2.R.
#   - It does not call library(), rm(list = ls()), or load().
############################################################

# ==========================================================
# 0. Internal helpers
# ==========================================================

#' Safe scalar quantile
#'
#' @param x Numeric vector.
#' @param p Numeric scalar in [0, 1].
#'
#' @return Numeric scalar.
.q <- function(x, p) {
  stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7)
}

#' Clamp values to a closed interval
#'
#' @param x Numeric vector or scalar.
#' @param lo Lower bound.
#' @param hi Upper bound.
#'
#' @return Numeric vector or scalar with values truncated to [lo, hi].
.clamp <- function(x, lo, hi) {
  pmin(hi, pmax(lo, x))
}

#' Construct a reproducible binary group factor from expression and cutoff
#'
#' @param expr Numeric expression vector.
#' @param cut Numeric scalar cutoff.
#'
#' @return Factor with fixed levels c("Low", "High").
.make_group <- function(expr, cut) {
  factor(ifelse(expr > cut, "High", "Low"), levels = c("Low", "High"))
}

# ==========================================================
# 1. Preprocessing utilities
# ==========================================================

#' Winsorize a numeric vector
#'
#' This helper limits extreme values to the specified lower and upper quantile
#' thresholds. In StaBiCut v2 it is used to reduce the influence of extreme
#' tumor-expression outliers on cutoff scanning and downstream survival fits.
#'
#' @param x Numeric vector.
#' @param lower Lower winsorization quantile. Default is 0.05.
#' @param upper Upper winsorization quantile. Default is 0.95.
#'
#' @return Numeric vector with extreme values truncated.
#' @export
winsorize <- function(x, lower = 0.05, upper = 0.95) {
  qnt <- stats::quantile(x, probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  x[x < qnt[1]] <- qnt[1]
  x[x > qnt[2]] <- qnt[2]
  x
}

# ==========================================================
# 2. Direction-prior inference
# ==========================================================

#' Infer the expected hazard direction for high expression
#'
#' Direction labels follow the StaBiCut v2 convention:
#'   - "protective_high": higher expression is expected to associate with
#'     lower hazard (HR < 1 for High vs Low)
#'   - "adverse_high": higher expression is expected to associate with
#'     higher hazard (HR > 1 for High vs Low)
#'   - NULL: no direction constraint is imposed during cutoff scanning
#'
#' Priority of evidence:
#'   1) Gene-specific prior table (if available)
#'   2) Fallback inference from tumor-normal shift (TN_log2FC) and continuous Cox beta,
#'      but only when both are available and concordant
#'   3) Stress-testing mode based on tumor-normal shift alone
#'
#' @param gene Character scalar; gene symbol.
#' @param df Data frame with required columns: expr, time, event.
#' @param tumor_values Numeric vector of full tumor expression values.
#' @param normal_values Numeric vector of full normal expression values.
#' @param gene_prior_table Optional data frame containing columns Gene and
#'   expected_dir.
#' @param force_direction Logical; if TRUE, use the tumor-normal shift (TN_log2FC) alone
#'   as a stress-testing heuristic when no explicit gene-specific prior is used.
#' @param min_abs_beta Numeric; Optional minimum absolute continuous Cox coefficient
#'   required before fallback inference is accepted.
#' @param allow_tumor_only_direction_heuristic Logical; if TRUE, a weak fallback based on
#'   tumor median versus tumor lower quartile is allowed when normal samples are
#'   unavailable.
#'
#' @return Character scalar ("protective_high" or "adverse_high") or NULL.
#' @export
infer_expected_direction <- function(
    gene,
    df,
    tumor_values,
    normal_values,
    gene_prior_table = NULL,
    force_direction = FALSE,
    min_abs_beta = 0,
    allow_tumor_only_direction_heuristic = FALSE
) {
  # --- A) gene-specific prior table ---
  if (!is.null(gene_prior_table) &&
      all(c("Gene", "expected_dir") %in% colnames(gene_prior_table)) &&
      gene %in% gene_prior_table$Gene) {
    dir0 <- gene_prior_table$expected_dir[match(gene, gene_prior_table$Gene)]
    if (dir0 %in% c("protective_high", "adverse_high")) 
      return(dir0)
     }

  has_normal <- !is.null(normal_values) && length(normal_values) > 0 && any(!is.na(normal_values))
  has_tumor  <- !is.null(tumor_values)  && length(tumor_values)  > 0 && any(!is.na(tumor_values))

  TN_log2FC <- NA_real_
  if (has_tumor && has_normal) {
    TN_log2FC <- stats::median(tumor_values, na.rm = TRUE) -
      stats::median(normal_values, na.rm = TRUE)
  } else if (allow_tumor_only_direction_heuristic && has_tumor) {
    # Weak tumor-only fallback used only when normal samples are unavailable.
    # This is a heuristic proxy and should not be interpreted as a true
    # tumor–normal expression shift.
    TN_log2FC <- stats::median(tumor_values, na.rm = TRUE) - .q(tumor_values, 0.25)
  }

  beta <- NA_real_
  fit_cont <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ expr, data = df),
    error = function(e) NULL
  )
  if (!is.null(fit_cont)) {
    beta <- as.numeric(stats::coef(fit_cont)[1])
  }

  # --- B) infer only if TN and beta are concordant ---
  if (!force_direction) {
    if (!is.na(TN_log2FC) && !is.na(beta) && abs(beta) >= min_abs_beta) {
      if (TN_log2FC < 0 && beta < 0) return("protective_high")
      if (TN_log2FC > 0 && beta > 0) return("adverse_high")
    }
    return(NULL)
  }

  # --- C) stress-testing mode ---
  if (!is.na(TN_log2FC)) {
    if (TN_log2FC < 0) return("protective_high")
    if (TN_log2FC > 0) return("adverse_high")
  }

  NULL
}

# ==========================================================
# 3. Deterministic cutoff scan (direction constrained)
# ==========================================================

#' Deterministically evaluate candidate cutoffs within the admissible
#' central expression range
#'
#' Candidate cutoffs are generated within the central tumor-expression range
#' defined by `minprop`. Each admissible cutoff is used to dichotomize samples
#' into Low/High groups and is then evaluated by univariable Cox proportional
#' hazards regression. Candidates failing group-size, event-count, or model-
#' stability criteria are discarded. When an expected direction is available,
#' only direction-compatible candidates are retained for final ranking.
#'
#' Under the default hierarchical ranking rule used in the current StaBiCut v2
#' runner, candidates are ordered by:
#'   1) decreasing scan score (Wald statistic by default),
#'   2) increasing P value,
#'   3) a final deterministic ordering term when required.
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   `expr` is the winsorized tumor-expression vector used for cutoff scanning;
#'   `time` is follow-up time; `event` is the binary event indicator.
#' @param expected_dir character scalar or `NULL`. Allowed directional labels
#'   are `"protective_high"` and `"adverse_high"`. If `NULL`, no direction
#'   constraint is applied during candidate filtering.
#' @param minprop numeric scalar in `(0, 0.5)`. Minimum allowable sample
#'   fraction in either dichotomized group.
#' @param grid character string specifying candidate generation mode.
#'   Must be one of `"unique"` or `"quantile"`.
#' @param n_grid integer scalar. Number of quantile-grid points used to propose
#'   candidate cutoffs when `grid = "quantile"`.
#' @param tie_method character string specifying the deterministic candidate-
#'   ranking rule. Must be one of `"max_wald_then_min_p"`, `"max_wald"`, or
#'   `"min_p"`.
#' @param min_events_per_group integer scalar. Minimum number of events required
#'   in each dichotomized group for a candidate cutoff to remain admissible.
#' @param se_max numeric scalar. Upper bound on the Cox standard error used as a
#'   practical filter against unstable or separation-like solutions.
#' @param coef_max numeric scalar. Upper bound on the absolute Cox coefficient
#'   used as a practical filter against unstable or separation-like solutions.
#' @param ref_cut numeric scalar or `NULL`. Optional reference cutoff used for
#'   sticky scoring when distance penalization is enabled.
#' @param lambda non-negative numeric scalar. Penalty strength for distance from
#'   `ref_cut`. A distance-based penalty is applied only when `lambda > 0` and
#'   `ref_cut` is finite.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{best_cut}{Numeric scalar giving the selected cutoff, or `NA_real_`
#'       if no valid candidate remains.}
#'     \item{expected_dir}{Character scalar or `NULL` indicating the direction
#'       constraint carried into the scan.}
#'     \item{scan_df}{A data.frame of admissible candidate cutoffs and their
#'       scan statistics, or `NULL` if no valid candidates remain.}
#'     \item{reason}{Character scalar describing the scan outcome.}
#'     \item{n_candidates}{Integer scalar giving the number of generated
#'       candidate cutoffs entering the scan before final admissibility
#'       filtering.}
#'   }
#' @export
scan_cutpoints_v2 <- function(
    df,
    expected_dir = NULL,
    minprop = 0.25,
    grid = c("unique", "quantile"),
    n_grid = 200,
    tie_method = c("max_wald_then_min_p", "max_wald", "min_p"),
    min_events_per_group = 5,
    se_max = 3,
    coef_max = 5,
    ref_cut = NULL,
    lambda = 0.0
) {
  stopifnot(all(c("expr", "time", "event") %in% colnames(df)))

  grid <- match.arg(grid)
  tie_method <- match.arg(tie_method)

  lo <- .q(df$expr, minprop)
  hi <- .q(df$expr, 1 - minprop)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    return(list(
      best_cut = NA_real_,
      expected_dir = expected_dir,
      scan_df = NULL,
      reason = "Invalid cut-search range (lo >= hi)",
      n_candidates = 0
    ))
  }

  if (grid == "unique") {
    cuts <- sort(unique(df$expr))
  } else {
    probs <- seq(minprop, 1 - minprop, length.out = n_grid)
    cuts <- sort(unique(stats::quantile(df$expr, probs = probs, na.rm = TRUE, names = FALSE)))
  }

  cuts <- cuts[cuts > lo & cuts < hi]
  if (length(cuts) < 10) {
    return(list(
      best_cut = NA_real_,
      expected_dir = expected_dir,
      scan_df = NULL,
      reason = "Too few candidate cutoffs after interval filtering",
      n_candidates = length(cuts)
    ))
  }

  out <- lapply(cuts, function(cut) {
    g <- .make_group(df$expr, cut)
    tab <- table(g)
    if (any(tab < 2)) return(NULL)

    n_low  <- as.integer(tab["Low"])
    n_high <- as.integer(tab["High"])
    n_total <- n_low + n_high
    prop_low  <- n_low / n_total
    prop_high <- n_high / n_total
    if (min(prop_low, prop_high) < minprop) return(NULL)

    ev_low  <- sum(df$event[g == "Low"]  == 1, na.rm = TRUE)
    ev_high <- sum(df$event[g == "High"] == 1, na.rm = TRUE)
    if (min(ev_low, ev_high) < min_events_per_group) return(NULL)

    dd <- data.frame(time = df$time, event = df$event, g = g)
    fit <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ g, data = dd),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    s <- summary(fit)
    coef1 <- as.numeric(stats::coef(fit)[1])
    se1   <- as.numeric(s$coefficients[, "se(coef)"][1])

   # --- built-in separation / edge-peak filtering ---
    if (!is.finite(coef1) || !is.finite(se1)) return(NULL)
    if (se1 > se_max || abs(coef1) > coef_max) return(NULL)

    hr <- as.numeric(s$coefficients[, "exp(coef)"][1])
    wald <- as.numeric(s$wald["test"])
    p <- as.numeric(s$coefficients[, "Pr(>|z|)"][1])
    if (!is.finite(hr) || !is.finite(wald) || !is.finite(p)) return(NULL)

    data.frame(
      cut = cut,
      hr = hr,
      wald = wald,
      p = p,
      coef = coef1,
      se = se1,
      n_low = n_low,
      n_high = n_high,
      ev_low = ev_low,
      ev_high = ev_high,
      stringsAsFactors = FALSE
    )
  })

  scan_df <- do.call(rbind, out)
  if (is.null(scan_df) || nrow(scan_df) == 0) {
    return(list(
      best_cut = NA_real_,
      expected_dir = expected_dir,
      scan_df = NULL,
      reason = "No valid cutoffs after Cox fitting and feasibility filters",
      n_candidates = length(cuts)
    ))
  }

 # direction constraint
  if (!is.null(expected_dir)) {
    if (expected_dir == "protective_high") scan_df <- scan_df[scan_df$hr < 1, , drop = FALSE]
    if (expected_dir == "adverse_high")    scan_df <- scan_df[scan_df$hr > 1, , drop = FALSE]
  }
  if (nrow(scan_df) == 0) {
    return(list(
      best_cut = NA_real_,
      expected_dir = expected_dir,
      scan_df = NULL,
      reason = "All candidate cutoffs were removed by the direction constraint",
      n_candidates = length(cuts)
    ))
  }

  # Optional distance-to-reference penalty for candidate scoring (useful for local bootstrap support)
  scan_df$score <- scan_df$wald
  if (is.finite(lambda) && lambda > 0 && is.finite(ref_cut)) {
    iqr0 <- stats::IQR(df$expr, na.rm = TRUE)
    if (is.finite(iqr0) && iqr0 > 0) {
      scan_df$score <- scan_df$wald - lambda * abs(scan_df$cut - ref_cut) / iqr0
    }
  }

  ord <- switch(
    tie_method,
    "max_wald_then_min_p" = order(-scan_df$score, scan_df$p, abs(scan_df$hr - 1)),
    "max_wald"            = order(-scan_df$score),
    "min_p"               = order(scan_df$p, -scan_df$score)
  )

  list(
    best_cut = as.numeric(scan_df$cut[ord[1]]),
    expected_dir = expected_dir,
    scan_df = scan_df,
    reason = "OK",
    n_candidates = length(cuts)
  )
}

# ==========================================================
# 4. Cox model at selected cutoff
# ==========================================================

#' Fit a cutoff-defined dichotomized Cox model
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   `expr` is dichotomized at `cutoff` to define the Low/High groups used in
#'   the Cox proportional hazards model.
#' @param cutoff Numeric scalar used to define the dichotomized groups.
#'
#' @return A named list containing the cutoff, hazard ratio, 95% confidence
#'   interval, Wald P value, group sizes, event counts, and the fitted Cox model
#'   object. If model fitting fails, the same output structure is returned with
#'   `NA` summary statistics and `fit = NULL`.
#' @export
fit_cutoff_cox_v2 <- function(df, cutoff) {
  g <- .make_group(df$expr, cutoff)
  tab <- table(g)
  ev_low  <- sum(df$event[g == "Low"]  == 1, na.rm = TRUE)
  ev_high <- sum(df$event[g == "High"] == 1, na.rm = TRUE)

  dd <- data.frame(time = df$time, event = df$event, g = g)
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ g, data = dd),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(
      cutoff = cutoff,
      hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p = NA_real_,
      n_low = as.integer(tab["Low"]),
      n_high = as.integer(tab["High"]),
      events_low = ev_low,
      events_high = ev_high,
      fit = NULL
    ))
  }

  s <- summary(fit)
  list(
    cutoff = as.numeric(cutoff),
    hr = as.numeric(s$coefficients[, "exp(coef)"][1]),
    ci_low = as.numeric(s$conf.int[, "lower .95"][1]),
    ci_high = as.numeric(s$conf.int[, "upper .95"][1]),
    p = as.numeric(s$coefficients[, "Pr(>|z|)"][1]),
    n_low = as.integer(tab["Low"]),
    n_high = as.integer(tab["High"]),
    events_low = ev_low,
    events_high = ev_high,
    fit = fit
  )
}

# ==========================================================
# 5. Direction-consistent local bootstrap support
# ==========================================================

#' Bootstrap local support for the selected cutoff under an optional
#' direction constraint
#'
#' This procedure is intentionally local. When `restrict_window = TRUE`,
#' candidate cutoffs in each bootstrap replicate are proposed only from a
#' neighborhood centered on the originally selected cutoff. The aim is to assess
#' local reproducibility of the selected threshold rather than to repeatedly
#' re-optimize a new global cutoff in each resampled dataset.
#'
#' Importantly, local restriction affects only candidate generation. Once a
#' bootstrap-selected cutoff is proposed, its hazard direction is re-evaluated
#' on the full bootstrap sample before retention.
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   Bootstrap resampling is performed row-wise on this table.
#' @param best_cut Numeric scalar giving the cutoff selected from the primary
#'   scan.
#' @param expected_dir Character scalar or `NULL`. If provided, the same
#'   direction constraint is applied during bootstrap rescanning and rechecked
#'   on the full bootstrap sample before a cutoff is retained.
#' @param n_boot Integer scalar giving the number of bootstrap replicates.
#' @param minprop Numeric scalar in `(0, 0.5)` giving the minimum allowable
#'   sample fraction in either group.
#' @param grid Character string specifying candidate generation mode:
#'   `"quantile"` or `"unique"`.
#' @param n_grid Integer scalar giving the number of candidate grid points used
#'   in each bootstrap rescan.
#' @param restrict_window Logical; if `TRUE`, candidate generation is restricted
#'   to a local window around `best_cut`.
#' @param window_iqr_mult Numeric scalar giving the width multiplier applied to
#'   `IQR(df$expr)` when constructing the local candidate window.
#' @param min_valid Integer scalar giving the minimum number of retained
#'   bootstrap cutoffs required before summary statistics are considered
#'   interpretable.
#' @param min_events_per_group Integer scalar giving the minimum required number
#'   of events in each dichotomized group.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{boot_cuts}{Numeric vector of retained bootstrap cutoffs.}
#'     \item{n_valid}{Integer scalar giving the number of retained bootstrap
#'       cutoffs.}
#'     \item{median}{Numeric scalar giving the median retained cutoff, or
#'       `NA_real_` if `n_valid < min_valid`.}
#'     \item{iqr}{Numeric scalar giving the interquartile range of retained
#'       cutoffs, or `NA_real_` if `n_valid < min_valid`.}
#'     \item{ci_low}{Numeric scalar giving the 2.5th percentile of retained
#'       cutoffs, or `NA_real_` if `n_valid < min_valid`.}
#'     \item{ci_high}{Numeric scalar giving the 97.5th percentile of retained
#'       cutoffs, or `NA_real_` if `n_valid < min_valid`.}
#'     \item{mode_hint}{Numeric scalar giving the location of the highest peak in the
#'      kernel density estimate of the retained bootstrap cutoffs, or `NA_real_` if
#'      `n_valid < min_valid`.}
#       mode_hint is a descriptive, density-based summary of the retained bootstrap
#       cutoff distribution. It is not used in the current stability score or
#       composite ranking, but may be useful for diagnostic plotting.
#'     \item{reason}{Character scalar describing the bootstrap outcome.}
#'   }
#' @export
bootstrap_cutoffs_v2 <- function(
    df,
    best_cut,
    expected_dir = NULL,
    n_boot = 500,
    minprop = 0.25,
    grid = c("quantile", "unique"),
    n_grid = 150,
    restrict_window = TRUE,
    window_iqr_mult = 0.5,
    min_valid = 50,
    min_events_per_group = 5
) {
  grid <- match.arg(grid)

 # define window
  lo_win <- -Inf
  hi_win <- Inf
  if (restrict_window && is.finite(best_cut)) {
    win <- stats::IQR(df$expr, na.rm = TRUE) * window_iqr_mult
    if (is.finite(win) && win > 0) {
      lo_win <- best_cut - win
      hi_win <- best_cut + win
    }
  }

  boot_cuts <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(nrow(df), replace = TRUE)
    dfb <- df[idx, , drop = FALSE]
    dfb2 <- dfb

    if (is.finite(lo_win) || is.finite(hi_win)) {
      dfb2$._win_ok <- dfb2$expr > lo_win & dfb2$expr < hi_win
    } else {
      dfb2$._win_ok <- TRUE
    }

    df_cand <- dfb2[dfb2$._win_ok, c("expr", "time", "event"), drop = FALSE]
    if (nrow(df_cand) < 50) next

    scan <- scan_cutpoints_v2(
      df = df_cand,
      expected_dir = expected_dir,
      minprop = minprop,
      grid = grid,
      n_grid = n_grid,
      tie_method = "max_wald_then_min_p",
      min_events_per_group = min_events_per_group
    )

    if (is.finite(scan$best_cut)) {
      fitb <- fit_cutoff_cox_v2(dfb2, scan$best_cut)
      if (!is.na(fitb$hr)) {
        if (is.null(expected_dir) ||
            (expected_dir == "protective_high" && fitb$hr < 1) ||
            (expected_dir == "adverse_high"    && fitb$hr > 1)) {
          boot_cuts[b] <- scan$best_cut
        }
      }
    }
  }

  boot_cuts <- boot_cuts[is.finite(boot_cuts)]
  n_valid <- length(boot_cuts)

  if (n_valid < min_valid) {
    return(list(
      boot_cuts = boot_cuts,
      n_valid = n_valid,
      median = NA_real_,
      iqr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      mode_hint = NA_real_,
      reason = paste0("Too few valid bootstrap cutoffs (", n_valid, " < ", min_valid, ")")
    ))
  }

  dens <- stats::density(boot_cuts)
  mode_hint <- dens$x[which.max(dens$y)]
  ci <- stats::quantile(boot_cuts, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  list(
    boot_cuts = boot_cuts,
    n_valid = n_valid,
    median = stats::median(boot_cuts, na.rm = TRUE),
    iqr = stats::IQR(boot_cuts, na.rm = TRUE),
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    mode_hint = as.numeric(mode_hint),
    reason = "OK"
  )
}

# ==========================================================
# 6. Post-selection distributional and directional metrics
# ==========================================================

#' Quantile-based expression range
#' @param x Numeric vector.
#' @param lo Lower quantile. Default is 0.05.
#' @param hi Upper quantile. Default is 0.95.
#' @return Numeric scalar giving the difference between the `hi` and `lo`
#'   quantiles of `x`.
#' @export
calc_expression_range <- function(x, lo = 0.05, hi = 0.95) {
  q <- stats::quantile(x, probs = c(lo, hi), na.rm = TRUE, names = FALSE)
  as.numeric(q[2] - q[1])
}

#' Empirical quantile rank of a cutoff within the tumor-expression distribution
#' @param expr_tumor Numeric tumor-expression vector.
#' @param cutoff Numeric scalar cutoff.
#' @return Numeric scalar in [0, 1].
#' @export
calc_quantile_rank   <- function(expr_tumor, cutoff) {
    ec <- stats::ecdf(as.numeric(expr_tumor))
    as.numeric(ec(cutoff))
  }

#' Centrality-based cutoff plausibility score
#'
#' This quantity is stored historically as `density_score` in the current
#' StaBiCut v2 code path, but mathematically it reflects ECDF-based centrality
#' rather than a kernel-density estimate.
#' @param quantile_rank Numeric scalar in [0, 1].
#' @return Numeric scalar in [0, 1].
#' @export
calc_density_score <- function(quantile_rank) {
  1 - abs(quantile_rank - 0.5)
}

#' Group balance at the selected cutoff
#' @param expr_tumor Numeric tumor-expression vector.
#' @param cutoff Numeric scalar cutoff.
#' @return List with low_n, high_n, and min_prop.
#' @export
calc_group_balance <- function(expr_tumor, cutoff) {
  grp <- ifelse(expr_tumor > cutoff, "High", "Low")
  tab <- table(grp)
  n_low <- as.integer(tab["Low"])
  n_high <- as.integer(tab["High"])
  n_total <- n_low + n_high
  min_prop <- min(n_low, n_high) / n_total
  list(low_n = n_low, high_n = n_high, min_prop = as.numeric(min_prop))
}

#' Tumor-normal median shift on log2(TPM+1) scale
#' @param tumor_values Numeric tumor-expression vector.
#' @param normal_values Numeric normal-expression vector.
#' @return Numeric scalar.
#' @export
calc_TN_log2FC <- function(tumor_values, normal_values) {
  stats::median(tumor_values, na.rm = TRUE) - stats::median(normal_values, na.rm = TRUE)
}

#' Directional concordance between tumor-normal shift and cutoff-based log(HR)
#' @param TN_log2FC Numeric scalar.
#' @param HR Numeric scalar hazard ratio.
#' @return Numeric scalar: `+1` for concordant directions, `-1` for discordant
#'   directions, and `NA_real_` otherwise.
#' @export
calc_direction_concordance <- function(TN_log2FC, HR) {
  if (is.na(TN_log2FC) || is.na(HR) || HR <= 0) return(NA_real_)
  s1 <- sign(TN_log2FC)
  s2 <- sign(log(HR))
  if (s1 == 0 || s2 == 0) return(NA_real_)
  as.numeric(s1 * s2)
}

#' Evaluate post-selection distributional and directional metrics for a selected cutoff
#'
#' @param expr_tumor Numeric tumor-expression vector used to evaluate cutoff
#'   position, group balance, and expression-range metrics.
#' @param cutoff Numeric scalar giving the selected cutoff.
#' @param tumor_values Numeric tumor-expression vector used for the tumor-normal
#'   median-shift calculation.
#' @param normal_values Numeric normal-expression vector used for the tumor-normal
#'   median-shift calculation.
#' @param HR Numeric scalar giving the hazard ratio from the cutoff-defined Cox
#'   model.
#'
#' @return A named list containing cutoff-position metrics (`quantile_rank`,
#'   `density_score`, and `expression_range`), group-balance metrics
#'   (`min_group_prop`, `low_group_n`, and `high_group_n`), the tumor-normal
#'   median shift (`TN_log2FC`), and the directional concordance score.
#' @export
analyze_cutoff_position_v2 <- function(expr_tumor, cutoff, tumor_values, normal_values, HR) {
  expr_tumor <- as.numeric(expr_tumor)
  quantile_rank <- calc_quantile_rank(expr_tumor, cutoff)
  density_score <- calc_density_score(quantile_rank)
  bal <- calc_group_balance(expr_tumor, cutoff)
  TN_log2FC <- calc_TN_log2FC(tumor_values, normal_values)
  direction_concordance <- calc_direction_concordance(TN_log2FC, HR)
  expression_range <- calc_expression_range(expr_tumor, lo = 0.05, hi = 0.95)

  list(
    quantile_rank = quantile_rank,
    density_score = density_score,
    min_group_prop = bal$min_prop,
    low_group_n = bal$low_n,
    high_group_n = bal$high_n,
    expression_range = expression_range,
    TN_log2FC = TN_log2FC,
    direction_concordance = direction_concordance
  )
}

#' Assess hazard-direction consistency of retained bootstrap cutoffs on the
#' original dataset
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   Retained bootstrap cutoffs are re-evaluated on this original analysis
#'   table.
#' @param bootstrap_cutoffs Numeric vector of retained bootstrap cutoffs.
#' @param hr_ref Numeric scalar giving the reference hazard ratio from the
#'   selected cutoff on the original dataset.
#'
#' @return Numeric scalar in `[0, 1]` giving the proportion of retained
#'   bootstrap cutoffs whose re-estimated hazard direction matches the
#'   reference direction on the original dataset, or `NA_real_` if direction
#'   consistency cannot be evaluated reliably.
#' @export
compute_hazard_dir_consistency_v2 <- function(df, bootstrap_cutoffs, hr_ref) {
  x <- bootstrap_cutoffs[is.finite(bootstrap_cutoffs)]
  if (length(x) < 30) return(NA_real_)
  if (is.na(hr_ref) || hr_ref <= 0) return(NA_real_)

  ref_dir <- sign(log(hr_ref))
  if (is.na(ref_dir) || ref_dir == 0) return(NA_real_)

  dir_vec <- vapply(x, function(cut) {
    g2 <- .make_group(df$expr, cut)
    dd <- data.frame(time = df$time, event = df$event, g2 = g2)
    m <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ g2, data = dd),
      error = function(e) NULL
    )
    if (is.null(m)) return(NA_real_)
    hr <- exp(stats::coef(m))[1]
    if (!is.finite(hr) || hr <= 0) return(NA_real_)
    sign(log(hr))
  }, numeric(1))

  dir_vec <- dir_vec[is.finite(dir_vec)]
  if (length(dir_vec) < 30) return(NA_real_)
  mean(dir_vec == ref_dir)
}

#' Evaluate spline support using a likelihood-ratio comparison
#'
#' @param df data.frame with numeric columns `expr`, `time`, and `event`.
#'   `expr` is modeled as a natural spline with 3 degrees of freedom in the Cox
#'   model.
#' @param gene Optional gene label retained only for compatibility.
#' @param plot_dir Unused placeholder retained only for compatibility with
#'   earlier code.
#'
#' @return A named list with components `p_spline` and `spline_support`, where
#'   `p_spline` is the likelihood-ratio P value comparing the null Cox model
#'   with the spline-based Cox model, and `spline_support` is a logical flag
#'   indicating whether `p_spline < 0.05`.
#' @export
evaluate_spline_support_v2 <- function(df, gene = NA_character_, plot_dir = NULL) {
  p_spline <- NA_real_
  spline_support <- FALSE

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
        spline_support <- is.finite(p_spline) && p_spline < 0.05
      }
    }
  }

  list(p_spline = p_spline, spline_support = spline_support)
}

