############################################################
# StaBiCut v2 — Meta biological consistency layer
# File: StaBiCut_v2_Meta_biological_consistency_20260731.R
#
# Purpose:
#   Build an external multi-cohort tumor-normal (TN) biological-consistency
#   layer for the 11 StaBiCut temporal candidate genes.
#
# Design:
#   - GSE39582 remains tumor-only survival/cutoff validation.
#   - Independent GEO tumor-normal cohorts are NOT merged into one expression
#     matrix. Each cohort is analyzed internally first to avoid cohort/platform
#     confounding.
#   - Per-cohort TN effects are meta-analyzed across cohorts to generate a
#     Meta-TN reference layer that can replace the TCGA-reference TN sensitivity
#     layer in external reporting.
#
# Expected local input directory:
#   E:/OpenCode/ZG16_validation/Meta_data
# containing files such as:
#   GSE21510_series_matrix.txt.gz
#   GSE32323_series_matrix.txt.gz
#   GSE41258_series_matrix.txt.gz
#
# Optional input:
#   E:/OpenCode/ZG16_validation/GSE39582_StaBiCut_external_validation/
#     GSE39582_StaBiCut_external_results.csv
#   If present, this script appends a Meta-TN layer to the GSE39582 external
#   StaBiCut result table and recalculates the five-layer external score.
############################################################

rm(list = ls())

## =========================
## 0. Paths and configuration
## =========================
base_dir <- "E:/OpenCode/ZG16_validation"
meta_dir <- file.path(base_dir, "New_meta")
out_dir  <- file.path(base_dir, "TCGA_meta0731")
pdf_dir  <- file.path(out_dir, "figures")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Optional GSE39582 external result table from v7 runner
external_result_file <- file.path(
  base_dir,
  "GSE39582_StaBiCut_external_validation",
  "GSE39582_StaBiCut_external_results.csv"
)

## Candidate genes and expected hazard directions used in StaBiCut v2.
## If expected_dir_table_11DEGs.csv exists, it is used; otherwise this fallback
## is used to keep the script runnable.
geneset <- c(
  "SPOCK2", "PYCR1", "CA4", "CES1", "ABCB1", "ZG16",
  "TNXB", "HMCN2", "MEP1A", "SLC37A2", "CHGB"
)

fallback_expected_dir <- data.frame(
  Gene = geneset,
  expected_dir = c(
    "adverse_high", "adverse_high", "protective_high", "protective_high",
    "protective_high", "protective_high", "protective_high", "protective_high",
    "protective_high", "protective_high", "protective_high"
  ),
  stringsAsFactors = FALSE
)

prior_candidates <- c(
  file.path(base_dir, "ZG16_validation", "expected_dir_table_11DEGs.csv"),
  file.path(base_dir, "expected_dir_table_11DEGs.csv"),
  file.path(base_dir, "scripts", "expected_dir_table_11DEGs.csv"),
  file.path(getwd(), "expected_dir_table_11DEGs.csv")
)
prior_candidates <- prior_candidates[file.exists(prior_candidates)]
if (length(prior_candidates) > 0) {
  expected_dir_table <- read.csv(prior_candidates[1], check.names = FALSE, stringsAsFactors = FALSE)
  expected_dir_table <- expected_dir_table[match(geneset, expected_dir_table$Gene), , drop = FALSE]
  expected_dir_table <- expected_dir_table[, intersect(c("Gene", "expected_dir", "TN_log2FC"), colnames(expected_dir_table)), drop = FALSE]
  message("Loaded expected direction table: ", prior_candidates[1])
} else {
  expected_dir_table <- fallback_expected_dir
  message("expected_dir_table_11DEGs.csv not found; using embedded fallback expected directions.")
}

## Dataset configuration.
## Notes:
##   - GSE21510/GSE32323 are GPL570; GSE41258 is GPL96.
##   - The primary Meta-TN model contains exactly three non-overlapping cohorts.
##   - include_regex/exclude_regex are intentionally conservative. All sample
##     assignments are exported for manual review.
dataset_config <- data.frame(
  GSE = c("GSE21510", "GSE32323", "GSE41258"),
  platform = c("GPL570", "GPL570", "GPL96"),
  file = file.path(
    meta_dir,
    paste0(c("GSE21510", "GSE32323", "GSE41258"), "_series_matrix.txt.gz")
  ),
  expected_tumor = c(19L, 17L, 182L),
  expected_normal = c(25L, 17L, 53L),
  ## Primary colorectal tumor terms.
  tumor_regex = c(
    "(^|[,;: _-])(cancer|tumou?r|carcinoma|adenocarcinoma|crc)([,;: _-]|$)",
    "(^|[,;: _-])(cancer|tumou?r|carcinoma|adenocarcinoma|crc)([,;: _-]|$)",
    "primary[ _-]*tumou?r|primary[ _-]*cancer|colon[ _-]*adenocarcinoma|colon[ _-]*carcinoma"
  ),
  ## Normal colorectal mucosa terms.
  normal_regex = c(
    "normal|non[- ]?cancer|non[- ]?tumou?r|mucosa",
    "normal|non[- ]?cancer|non[- ]?tumou?r|mucosa",
    "normal[ _-]*(colon|colorectal|mucosa)|corresponding[ _-]*normal[ _-]*mucosa"
  ),
  ## Exclude non-primary disease states from the TN layer.
  exclude_regex = c(
    "cell[ _-]*line|aza|5[- ]?aza|metastasis|liver|lung|polyp|adenoma|xenograft|\\blcm\\b|laser[ _-]*capture|microdissect",
    "cell[ _-]*line|aza|5[- ]?aza|metastasis|liver|lung|polyp|adenoma|xenograft",
    "metastasis|normal[ _-]*liver|normal[ _-]*lung|polyp|adenoma|microadenoma|cell[ _-]*line|breast|prostate|xenograft|_ez"
  ),
  stringsAsFactors = FALSE
)

stopifnot(identical(dataset_config$GSE, c("GSE21510", "GSE32323", "GSE41258")))

## GSE21510 contains separately normalized LCM and homogenized specimens.
## Restrict the Meta-TN comparison to the homogenized tumor/normal subset.
exclude_lcm_in_GSE21510 <- TRUE

## =========================
## 1. Packages
## =========================
cran_pkgs <- c("dplyr", "tidyr", "stringr", "readr", "ggplot2", "openxlsx", "patchwork")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("AnnotationDbi", "hgu133plus2.db", "hgu133a.db")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

has_metafor <- requireNamespace("metafor", quietly = TRUE)
if (!has_metafor) {
  message("Package 'metafor' not found. The script will use fixed-effect inverse-variance meta-analysis as fallback.")
}

## =========================
## 2. GEO series matrix parser
## =========================
is_gzip_file <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && all(as.integer(magic) == c(31L, 139L))
}

read_text_lines_auto <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE) || is_gzip_file(path)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE, encoding = "UTF-8")
}

clean_geo_field <- function(x) {
  x <- gsub('^"|"$', "", x)
  x <- gsub('\\\\"', '"', x)
  x
}

read_gse_series_matrix <- function(file) {
  x <- read_text_lines_auto(file)
  i1 <- grep("^!series_matrix_table_begin", x)
  i2 <- grep("^!series_matrix_table_end", x)
  if (length(i1) != 1 || length(i2) != 1 || i2 <= i1) {
    stop("Cannot find valid series matrix table boundary in: ", file)
  }

  expr_txt <- x[(i1 + 1):(i2 - 1)]
  expr <- utils::read.delim(
    text = paste(expr_txt, collapse = "\n"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    sep = "\t"
  )
  probe_col <- colnames(expr)[1]
  rownames(expr) <- as.character(expr[[probe_col]])
  expr[[probe_col]] <- NULL
  expr <- as.matrix(expr)
  suppressWarnings(mode(expr) <- "numeric")

  meta_lines <- x[seq_len(i1 - 1)]
  sample_lines <- meta_lines[grepl("^!Sample_", meta_lines)]

  parse_sample_line <- function(pattern) {
    ln <- sample_lines[grepl(pattern, sample_lines)]
    if (length(ln) == 0) return(NULL)
    vals <- strsplit(ln[1], "\t", fixed = TRUE)[[1]]
    vals <- clean_geo_field(vals)
    vals[-1]
  }

  geo_accession <- parse_sample_line("^!Sample_geo_accession")
  title <- parse_sample_line("^!Sample_title")
  source_name <- parse_sample_line("^!Sample_source_name_ch1")

  if (!is.null(geo_accession) && length(geo_accession) == ncol(expr)) {
    colnames(expr) <- geo_accession
  }

  clin <- data.frame(
    sample = colnames(expr),
    geo_accession = if (!is.null(geo_accession) && length(geo_accession) == ncol(expr)) geo_accession else colnames(expr),
    title = if (!is.null(title) && length(title) == ncol(expr)) title else NA_character_,
    source_name = if (!is.null(source_name) && length(source_name) == ncol(expr)) source_name else NA_character_,
    stringsAsFactors = FALSE
  )

  char_lines <- sample_lines[grepl("^!Sample_characteristics_ch1", sample_lines)]
  if (length(char_lines) > 0) {
    char_mat <- lapply(char_lines, function(ln) {
      vals <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      vals <- clean_geo_field(vals)
      vals[-1]
    })

    for (k in seq_along(char_mat)) {
      vals <- char_mat[[k]]
      if (length(vals) != nrow(clin)) next
      key <- trimws(sub(":.*$", "", vals))
      val <- trimws(sub("^[^:]*:\\s*", "", vals))
      if (length(unique(key[!is.na(key) & nzchar(key)])) == 1) {
        nm <- make.names(unique(key[!is.na(key) & nzchar(key)]), unique = TRUE)
      } else {
        nm <- paste0("characteristics_", k)
      }
      clin[[nm]] <- val
    }
  }

  rownames(clin) <- clin$sample
  list(expr_probe = expr, clin_raw = clin)
}

## =========================
## 3. Probe annotation and sample classification
## =========================
get_probe_map <- function(probe_ids, platform) {
  if (platform == "GPL570") {
    pkg <- "hgu133plus2.db"
    keytype <- "PROBEID"
  } else if (platform == "GPL96") {
    pkg <- "hgu133a.db"
    keytype <- "PROBEID"
  } else {
    stop("Unsupported platform: ", platform)
  }

  ann_db <- get(pkg)
  AnnotationDbi::select(
    ann_db,
    keys = probe_ids,
    keytype = keytype,
    columns = c("SYMBOL", "GENENAME", "ENTREZID")
  ) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::distinct(PROBEID, SYMBOL, .keep_all = TRUE)
}

collapse_probe_to_gene <- function(expr_probe, probe_map, geneset = NULL) {
  common <- intersect(rownames(expr_probe), probe_map$PROBEID)
  expr_probe2 <- expr_probe[common, , drop = FALSE]
  map2 <- probe_map[match(common, probe_map$PROBEID), , drop = FALSE]

  probe_median <- apply(expr_probe2, 1, median, na.rm = TRUE)
  map2$probe_median <- probe_median[map2$PROBEID]

  keep_map <- map2 %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::arrange(dplyr::desc(probe_median), PROBEID, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  expr_gene <- expr_probe2[keep_map$PROBEID, , drop = FALSE]
  rownames(expr_gene) <- keep_map$SYMBOL
  if (!is.null(geneset)) {
    expr_gene <- expr_gene[intersect(geneset, rownames(expr_gene)), , drop = FALSE]
  }
  list(expr_gene = expr_gene, selected_probe_map = keep_map)
}

classify_samples <- function(clin, cfg_row) {
  ## Cohort-specific inclusion rules, aligned 1:1 with the external
  ## reproducibility forest meta-analysis (ZG16, 5 cohorts):
  ##   - GSE21510 / GSE32323: homogenized bulk tissues only (title-based).
  ##     LCM specimens are excluded because the LCM and homogenized batches
  ##     were normalized separately.
  ##   - GSE41258: official analysis set only ("included in analysis = Yes"),
  ##     tissue == "primary tumor" (Tumor) or "normal colon" (Normal).
  if (cfg_row$GSE %in% c("GSE21510", "GSE32323")) {
    title_l <- tolower(clin$title)
    class <- ifelse(grepl("cancer, homogenized", title_l, fixed = TRUE), "Tumor",
             ifelse(grepl("normal, homogenized", title_l, fixed = TRUE), "Normal",
                    "Other_or_excluded"))
    return(data.frame(
      sample = rownames(clin),
      class = class,
      review_text = clin$title,
      stringsAsFactors = FALSE
    ))
  }

  if (cfg_row$GSE == "GSE41258") {
    ## The first characteristics block mixes "tissue:" and "cell line:" keys,
    ## so the parser stores it under a generic name; locate the tissue column by
    ## content, preferring characteristics_* columns over titles.
    is_char_col <- grepl("^characteristics_", colnames(clin))
    find_col <- function(clin, pattern, prefer_char = TRUE) {
      hit <- vapply(clin, function(z) {
        is.character(z) && any(grepl(pattern, tolower(z), fixed = FALSE), na.rm = TRUE)
      }, logical(1))
      ord <- if (prefer_char) order(is_char_col, decreasing = TRUE) else seq_along(hit)
      nm <- names(hit)[ord][which(hit[ord])[1]]
      nm
    }
    tissue_col <- find_col(clin, "primary tumor|normal colon", prefer_char = TRUE)
    incl_col   <- if (!is.null(clin$included.in.analysis)) "included.in.analysis" else find_col(clin, "^included in analysis", prefer_char = FALSE)
    tissue <- if (!is.null(tissue_col)) tolower(clin[[tissue_col]]) else rep("", nrow(clin))
    incl   <- if (!is.null(incl_col)) tolower(clin[[incl_col]]) else rep("", nrow(clin))
    class <- ifelse(incl == "yes" & tissue == "primary tumor", "Tumor",
             ifelse(incl == "yes" & tissue == "normal colon", "Normal",
                    "Other_or_excluded"))
    return(data.frame(
      sample = rownames(clin),
      class = class,
      review_text = paste0(
        "tissue: ", clin[[if (!is.null(tissue_col)) tissue_col else colnames(clin)[1]]],
        " ; included in analysis: ", clin[[if (!is.null(incl_col)) incl_col else colnames(clin)[1]]]
      ),
      stringsAsFactors = FALSE
    ))
  }

  text_cols <- vapply(clin, function(z) is.character(z) || is.factor(z), logical(1))
  sample_text <- apply(clin[, text_cols, drop = FALSE], 1, function(z) paste(z, collapse = " ; "))
  sample_text_l <- tolower(sample_text)

  is_excl <- stringr::str_detect(sample_text_l, stringr::regex(cfg_row$exclude_regex, ignore_case = TRUE))
  is_tumor <- stringr::str_detect(sample_text_l, stringr::regex(cfg_row$tumor_regex, ignore_case = TRUE))
  is_normal <- stringr::str_detect(sample_text_l, stringr::regex(cfg_row$normal_regex, ignore_case = TRUE))

  if (cfg_row$GSE == "GSE21510" && isTRUE(exclude_lcm_in_GSE21510)) {
    is_excl <- is_excl | stringr::str_detect(sample_text_l, stringr::regex("\\bLCM\\b|laser microdissection", ignore_case = TRUE))
  }

  class <- rep("Other_or_excluded", length(sample_text_l))
  class[!is_excl & is_normal & !is_tumor] <- "Normal"
  class[!is_excl & is_tumor & !is_normal] <- "Tumor"
  ## If both regexes hit, use simple priority rules.
  class[!is_excl & is_tumor & is_normal & stringr::str_detect(sample_text_l, "normal")] <- "Normal"
  class[!is_excl & is_tumor & is_normal & !stringr::str_detect(sample_text_l, "normal")] <- "Tumor"

  data.frame(
    sample = rownames(clin),
    class = class,
    review_text = sample_text,
    stringsAsFactors = FALSE
  )
}

## =========================
## 4. Per-cohort TN effects
## =========================
hedges_g_effect <- function(x_tumor, x_normal) {
  x_tumor <- x_tumor[is.finite(x_tumor)]
  x_normal <- x_normal[is.finite(x_normal)]
  n1 <- length(x_tumor); n0 <- length(x_normal)
  if (n1 < 2 || n0 < 2) {
    return(c(g = NA_real_, se_g = NA_real_, d = NA_real_))
  }
  m1 <- mean(x_tumor); m0 <- mean(x_normal)
  s1 <- stats::sd(x_tumor); s0 <- stats::sd(x_normal)
  sp <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
  if (!is.finite(sp) || sp <= 0) return(c(g = NA_real_, se_g = NA_real_, d = NA_real_))
  d <- (m1 - m0) / sp
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- J * d
  var_g <- (n1 + n0) / (n1 * n0) + (g^2) / (2 * (n1 + n0 - 2))
  c(g = g, se_g = sqrt(var_g), d = d)
}

per_gene_effects_one_dataset <- function(expr_gene, sample_class, GSE, platform, geneset) {
  tumor_samples <- sample_class$sample[sample_class$class == "Tumor"]
  normal_samples <- sample_class$sample[sample_class$class == "Normal"]
  tumor_samples <- intersect(tumor_samples, colnames(expr_gene))
  normal_samples <- intersect(normal_samples, colnames(expr_gene))

  out <- lapply(geneset, function(gene) {
    if (!gene %in% rownames(expr_gene)) {
      return(data.frame(
        GSE = GSE, platform = platform, Gene = gene,
        n_tumor = length(tumor_samples), n_normal = length(normal_samples),
        mean_tumor = NA_real_, mean_normal = NA_real_,
        median_tumor = NA_real_, median_normal = NA_real_,
        logFC_mean = NA_real_, logFC_median = NA_real_, se_logFC = NA_real_,
        Hedges_g = NA_real_, se_Hedges_g = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    xt <- as.numeric(expr_gene[gene, tumor_samples])
    xn <- as.numeric(expr_gene[gene, normal_samples])
    xt <- xt[is.finite(xt)]; xn <- xn[is.finite(xn)]
    n1 <- length(xt); n0 <- length(xn)
    hg <- hedges_g_effect(xt, xn)
    se_lfc <- if (n1 >= 2 && n0 >= 2) sqrt(stats::var(xt) / n1 + stats::var(xn) / n0) else NA_real_
    data.frame(
      GSE = GSE, platform = platform, Gene = gene,
      n_tumor = n1, n_normal = n0,
      mean_tumor = ifelse(n1 > 0, mean(xt), NA_real_),
      mean_normal = ifelse(n0 > 0, mean(xn), NA_real_),
      median_tumor = ifelse(n1 > 0, median(xt), NA_real_),
      median_normal = ifelse(n0 > 0, median(xn), NA_real_),
      logFC_mean = ifelse(n1 > 0 && n0 > 0, mean(xt) - mean(xn), NA_real_),
      logFC_median = ifelse(n1 > 0 && n0 > 0, median(xt) - median(xn), NA_real_),
      se_logFC = se_lfc,
      Hedges_g = as.numeric(hg["g"]),
      se_Hedges_g = as.numeric(hg["se_g"]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

## =========================
## 5. Meta-analysis helpers
## =========================
meta_one_gene <- function(df, effect_col = "Hedges_g", se_col = "se_Hedges_g") {
  dd <- df %>%
    dplyr::filter(is.finite(.data[[effect_col]]), is.finite(.data[[se_col]]), .data[[se_col]] > 0)
  k <- nrow(dd)
  if (k == 0) {
    return(data.frame(k = 0, estimate = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
                      ci_low = NA_real_, ci_high = NA_real_, tau2 = NA_real_, I2 = NA_real_))
  }

  yi <- dd[[effect_col]]
  sei <- dd[[se_col]]
  vi <- sei^2

  if (has_metafor && k >= 2) {
    fit <- metafor::rma.uni(yi = yi, vi = vi, method = "REML")
    data.frame(
      k = k,
      estimate = as.numeric(fit$b[1]),
      se = as.numeric(fit$se),
      z = as.numeric(fit$zval),
      p = as.numeric(fit$pval),
      ci_low = as.numeric(fit$ci.lb),
      ci_high = as.numeric(fit$ci.ub),
      tau2 = as.numeric(fit$tau2),
      I2 = as.numeric(fit$I2),
      stringsAsFactors = FALSE
    )
  } else {
    ## Fixed-effect fallback.
    w <- 1 / vi
    est <- sum(w * yi) / sum(w)
    se <- sqrt(1 / sum(w))
    z <- est / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    data.frame(
      k = k,
      estimate = est,
      se = se,
      z = z,
      p = p,
      ci_low = est - 1.96 * se,
      ci_high = est + 1.96 * se,
      tau2 = NA_real_,
      I2 = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

calc_meta_support_score <- function(estimate, p, consistency_prop,
                                    effect_cap = 2, p_log_cap = 10) {
  effect_component <- pmin(abs(as.numeric(estimate)) / effect_cap, 1)
  p_component <- pmin(-log10(pmax(as.numeric(p), .Machine$double.xmin)) / p_log_cap, 1)
  score <- as.numeric(consistency_prop) * effect_component * p_component
  score[!is.finite(score)] <- NA_real_
  score
}

## =========================
## 6. Main analysis
## =========================
all_effects <- list()
all_counts <- list()
all_sample_reviews <- list()
all_probe_maps <- list()

for (i in seq_len(nrow(dataset_config))) {
  cfg <- dataset_config[i, ]
  message("\nProcessing ", cfg$GSE, " | ", cfg$platform)
  if (!file.exists(cfg$file)) {
    warning("Missing series matrix file: ", cfg$file)
    next
  }

  gse <- read_gse_series_matrix(cfg$file)
  probe_map <- get_probe_map(rownames(gse$expr_probe), cfg$platform)
  collapsed <- collapse_probe_to_gene(gse$expr_probe, probe_map, geneset = geneset)
  sample_class <- classify_samples(gse$clin_raw, cfg)
  sample_class$GSE <- cfg$GSE
  sample_class$platform <- cfg$platform

  observed_tumor <- sum(sample_class$class == "Tumor", na.rm = TRUE)
  observed_normal <- sum(sample_class$class == "Normal", na.rm = TRUE)
  if (observed_tumor != cfg$expected_tumor || observed_normal != cfg$expected_normal) {
    failed_audit <- file.path(out_dir, paste0(cfg$GSE, "_FAILED_sample_class_review.csv"))
    write.csv(sample_class, failed_audit, row.names = FALSE)
    stop(
      cfg$GSE, " sample-count audit failed: observed T=", observed_tumor,
      ", N=", observed_normal, "; expected T=", cfg$expected_tumor,
      ", N=", cfg$expected_normal, ". Review: ", failed_audit
    )
  }

  counts <- sample_class %>%
    dplyr::count(GSE, platform, class, name = "n") %>%
    tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
    dplyr::mutate(
      n_series_total = ncol(gse$expr_probe),
      n_genes_mapped_in_panel = nrow(collapsed$expr_gene)
    )

  eff <- per_gene_effects_one_dataset(
    expr_gene = collapsed$expr_gene,
    sample_class = sample_class,
    GSE = cfg$GSE,
    platform = cfg$platform,
    geneset = geneset
  )

  all_effects[[cfg$GSE]] <- eff
  all_counts[[cfg$GSE]] <- counts
  all_sample_reviews[[cfg$GSE]] <- sample_class
  selected_probe_map <- collapsed$selected_probe_map
  selected_probe_map$GSE <- cfg$GSE
  selected_probe_map$platform <- cfg$platform
  all_probe_maps[[cfg$GSE]] <- selected_probe_map

  write.csv(sample_class, file.path(out_dir, paste0(cfg$GSE, "_sample_class_review.csv")), row.names = FALSE)
  write.csv(eff, file.path(out_dir, paste0(cfg$GSE, "_per_gene_TN_effects.csv")), row.names = FALSE)
}

per_cohort_effects <- dplyr::bind_rows(all_effects)
sample_counts <- dplyr::bind_rows(all_counts)
sample_reviews <- dplyr::bind_rows(all_sample_reviews)
probe_maps <- dplyr::bind_rows(all_probe_maps)

write.csv(per_cohort_effects, file.path(out_dir, "MetaTN_per_cohort_effects_11genes.csv"), row.names = FALSE)
write.csv(sample_counts, file.path(out_dir, "MetaTN_dataset_sample_counts.csv"), row.names = FALSE)
write.csv(sample_reviews, file.path(out_dir, "MetaTN_all_sample_class_review.csv"), row.names = FALSE)
write.csv(probe_maps, file.path(out_dir, "MetaTN_selected_probe_per_gene_by_dataset.csv"), row.names = FALSE)

## Meta-analyze Hedges' g and mean logFC separately.
meta_g <- per_cohort_effects %>%
  dplyr::group_by(Gene) %>%
  dplyr::group_modify(~ meta_one_gene(.x, effect_col = "Hedges_g", se_col = "se_Hedges_g")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    meta_k = k,
    meta_Hedges_g = estimate,
    meta_Hedges_g_se = se,
    meta_Hedges_g_z = z,
    meta_Hedges_g_p = p,
    meta_Hedges_g_CI_low = ci_low,
    meta_Hedges_g_CI_high = ci_high,
    meta_Hedges_g_tau2 = tau2,
    meta_Hedges_g_I2 = I2
  )

meta_lfc <- per_cohort_effects %>%
  dplyr::group_by(Gene) %>%
  dplyr::group_modify(~ meta_one_gene(.x, effect_col = "logFC_mean", se_col = "se_logFC")) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    meta_logFC_k = k,
    meta_logFC_mean = estimate,
    meta_logFC_se = se,
    meta_logFC_z = z,
    meta_logFC_p = p,
    meta_logFC_CI_low = ci_low,
    meta_logFC_CI_high = ci_high,
    meta_logFC_tau2 = tau2,
    meta_logFC_I2 = I2
  )

consistency_tbl <- per_cohort_effects %>%
  dplyr::filter(is.finite(Hedges_g)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    n_cohorts_valid = dplyr::n(),
    n_negative = sum(Hedges_g < 0, na.rm = TRUE),
    n_positive = sum(Hedges_g > 0, na.rm = TRUE),
    dominant_direction = dplyr::case_when(
      n_negative > n_positive ~ "Tumor_down",
      n_positive > n_negative ~ "Tumor_up",
      TRUE ~ "Mixed_or_tie"
    ),
    consistency_prop = pmax(n_negative, n_positive) / n_cohorts_valid,
    .groups = "drop"
  )

meta_tn <- meta_g %>%
  dplyr::left_join(meta_lfc, by = "Gene") %>%
  dplyr::left_join(consistency_tbl, by = "Gene") %>%
  dplyr::left_join(expected_dir_table %>% dplyr::select(Gene, expected_dir), by = "Gene") %>%
  dplyr::mutate(
    expected_TN_direction_from_prior = ifelse(expected_dir == "protective_high", "Tumor_down", "Tumor_up"),
    direction_matches_prior = dominant_direction == expected_TN_direction_from_prior,
    MetaTN_base_support = calc_meta_support_score(
      estimate = meta_Hedges_g,
      p = meta_Hedges_g_p,
      consistency_prop = consistency_prop,
      effect_cap = 2,
      p_log_cap = 10
    ),
    MetaTN_prior_gated_support = ifelse(direction_matches_prior, MetaTN_base_support, 0),
    MetaTN_weighted_contribution = 0.15 * MetaTN_prior_gated_support,
    meta_TN_sign = sign(meta_Hedges_g)
  ) %>%
  dplyr::arrange(dplyr::desc(MetaTN_prior_gated_support), meta_Hedges_g_p)

write.csv(meta_tn, file.path(out_dir, "MetaTN_summary_11genes.csv"), row.names = FALSE)

## =========================
## 7. Optional: integrate Meta-TN with GSE39582 external StaBiCut output
## =========================
if (file.exists(external_result_file)) {
  ext <- read.csv(external_result_file, check.names = FALSE, stringsAsFactors = FALSE)

  required_ext <- c("Gene", "HR", "bootstrap_score", "hazard_score", "density_score", "balance_score")
  missing_ext <- setdiff(required_ext, colnames(ext))
  if (length(missing_ext) > 0) {
    warning("External result table missing columns: ", paste(missing_ext, collapse = ", "),
            ". Skipping integrated Meta-TN external scoring.")
  } else {
    ext2 <- ext %>%
      dplyr::left_join(meta_tn %>%
                         dplyr::select(Gene, meta_Hedges_g, meta_Hedges_g_p,
                                       meta_Hedges_g_CI_low, meta_Hedges_g_CI_high,
                                       meta_logFC_mean, meta_logFC_p,
                                       dominant_direction, consistency_prop,
                                       MetaTN_base_support,
                                       MetaTN_prior_gated_support),
                       by = "Gene") %>%
      dplyr::mutate(
        external_HR_sign = sign(log(as.numeric(HR))),
        meta_TN_HR_concordant = is.finite(meta_Hedges_g) & is.finite(as.numeric(HR)) &
          (sign(meta_Hedges_g) * sign(log(as.numeric(HR))) == 1),
        MetaTN_score_for_GSE39582 = ifelse(meta_TN_HR_concordant,
                                           MetaTN_prior_gated_support, 0),
        composite_score_external_with_MetaTN =
          0.40 * as.numeric(bootstrap_score) +
          0.15 * as.numeric(hazard_score) +
          0.15 * as.numeric(MetaTN_score_for_GSE39582) +
          0.20 * as.numeric(density_score) +
          0.10 * as.numeric(balance_score),
        Rank_external_with_MetaTN = rank(-composite_score_external_with_MetaTN,
                                         ties.method = "min", na.last = "keep")
      ) %>%
      dplyr::arrange(Rank_external_with_MetaTN, dplyr::desc(composite_score_external_with_MetaTN))

    write.csv(ext2, file.path(out_dir, "GSE39582_external_results_with_MetaTN.csv"), row.names = FALSE)

    ## Plot: no-TN score vs MetaTN five-layer score.
    ## Use separate plotting data frames because the two panels have different ranks.
    ## Left panel: order by no-TN rank/score.
    ## Right panel: order by Meta-TN rank/score.
    plot_df_left <- ext2 %>%
      dplyr::filter(
        is.finite(composite_score_external_no_TN),
        !is.na(Rank_external_no_TN)
      ) %>%
      dplyr::arrange(Rank_external_no_TN, dplyr::desc(composite_score_external_no_TN)) %>%
      dplyr::mutate(Gene = factor(Gene, levels = rev(Gene)),
                    score = composite_score_external_no_TN)

    plot_df_right <- ext2 %>%
      dplyr::filter(
        is.finite(composite_score_external_with_MetaTN),
        !is.na(Rank_external_with_MetaTN)
      ) %>%
      dplyr::arrange(Rank_external_with_MetaTN, dplyr::desc(composite_score_external_with_MetaTN)) %>%
      dplyr::mutate(Gene = factor(Gene, levels = rev(Gene)),
                    score = composite_score_external_with_MetaTN)

    if (nrow(plot_df_left) > 0 && nrow(plot_df_right) > 0) {
      fill_limits <- range(c(plot_df_left$score, plot_df_right$score), na.rm = TRUE)

      p1 <- ggplot2::ggplot(plot_df_left, ggplot2::aes(x = score, y = Gene, fill = score)) +
        ggplot2::geom_col(width = 0.70) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)),
                           hjust = -0.1, size = 3.4) +
        ggplot2::scale_fill_viridis_c(option = "C", direction = -1, limits = fill_limits) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(x = "GSE39582 external score without TN", y = NULL,
                      fill = "Composite\nscore",
                      title = "GSE39582 no-TN score") +
        ggplot2::coord_cartesian(
          xlim = c(0, max(plot_df_left$score, na.rm = TRUE) * 1.15)
        )

      p2 <- ggplot2::ggplot(plot_df_right, ggplot2::aes(x = score, y = Gene, fill = score)) +
        ggplot2::geom_col(width = 0.70) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)),
                           hjust = -0.1, size = 3.4) +
        ggplot2::scale_fill_viridis_c(option = "C", direction = -1, limits = fill_limits) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(x = "GSE39582 score with external Meta-TN layer", y = NULL,
                      fill = "Composite\nscore",
                      title = "GSE39582 + Meta-TN score") +
        ggplot2::coord_cartesian(
          xlim = c(0, max(plot_df_right$score, na.rm = TRUE) * 1.15)
        )

      ggplot2::ggsave(
        file.path(pdf_dir, "Supplementary_Figure_GSE39582_StaBiCut_external_with_MetaTN.pdf"),
        (p1 | p2) + patchwork::plot_layout(guides = "collect"),
        width = 11,
        height = 5.8
      )

      discovery_rank_file <- file.path(base_dir, "GSE39582_StaBiCut_external_validation",
                                       "TCGA_discovery_StaBiCut_rank_used.csv")
      if (file.exists(discovery_rank_file)) {
        discovery_rank_table <- read.csv(discovery_rank_file, stringsAsFactors = FALSE)
      } else {
        discovery_rank_table <- data.frame(
          Gene = c("CA4", "ZG16", "ABCB1", "CES1", "PYCR1", "CHGB",
                   "SLC37A2", "TNXB", "HMCN2", "MEP1A", "SPOCK2"),
          Rank_discovery_TCGA = seq_len(11),
          stringsAsFactors = FALSE
        )
      }

      rank_comp <- discovery_rank_table %>%
        dplyr::left_join(
          ext2 %>% dplyr::select(Gene, Rank_external_no_TN, Rank_external_with_MetaTN),
          by = "Gene"
        ) %>%
        dplyr::mutate(Gene = factor(Gene, levels = discovery_rank_table$Gene)) %>%
        tidyr::pivot_longer(
          cols = c(Rank_discovery_TCGA, Rank_external_no_TN, Rank_external_with_MetaTN),
          names_to = "ranking_type",
          values_to = "rank"
        ) %>%
        dplyr::filter(!is.na(rank)) %>%
        dplyr::mutate(
          ranking_type = factor(
            ranking_type,
            levels = c("Rank_discovery_TCGA", "Rank_external_no_TN", "Rank_external_with_MetaTN"),
            labels = c("TCGA discovery", "GSE39582 no TN", "GSE39582 + external Meta-TN")
          )
        )

      gene_colors <- c(
        "CA4" = "#280B54FF", "ZG16" = "#0D0887FF",
        "ABCB1" = "#5402A3FF", "CES1" = "#8B0AA5FF",
        "PYCR1" = "#CC4678FF", "CHGB" = "#E97158FF",
        "SLC37A2" = "#FBA139FF", "TNXB" = "#FADA24FF",
        "HMCN2" = "#9FDA3AFF", "MEP1A" = "#2DB27DFF",
        "SPOCK2" = "#40B7ADFF"
      )

      p_rank_shift <- ggplot2::ggplot(rank_comp,
                                      ggplot2::aes(x = ranking_type, y = rank,
                                                    group = Gene, label = Gene, color = Gene)) +
        ggplot2::geom_line(linewidth = 0.65) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::geom_text(size = 3, hjust = -0.05, check_overlap = TRUE, show.legend = FALSE) +
        ggplot2::scale_color_manual(values = gene_colors) +
        ggplot2::scale_y_reverse(breaks = seq_len(11)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(x = NULL, y = "Rank", color = "Gene",
                      title = "Ranking comparison with external Meta-TN biological layer")

      ggplot2::ggsave(
        file.path(pdf_dir, "Supplementary_Figure_GSE39582_rank_shift_with_MetaTN.pdf"),
        p_rank_shift,
        width = 9.2,
        height = 5.8
      )
    }
  }
} else {
  message("GSE39582 external result file not found; Meta-TN summary only was generated.")
}

## =========================
## 8. Workbook export
## =========================
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Dataset_counts")
openxlsx::writeDataTable(wb, "Dataset_counts", sample_counts)
openxlsx::addWorksheet(wb, "Per_cohort_effects")
openxlsx::writeDataTable(wb, "Per_cohort_effects", per_cohort_effects)
openxlsx::addWorksheet(wb, "MetaTN_summary")
openxlsx::writeDataTable(wb, "MetaTN_summary", meta_tn)
openxlsx::addWorksheet(wb, "Probe_mapping")
openxlsx::writeDataTable(wb, "Probe_mapping", probe_maps)
openxlsx::addWorksheet(wb, "Expected_direction")
openxlsx::writeDataTable(wb, "Expected_direction", expected_dir_table)
openxlsx::saveWorkbook(wb, file.path(out_dir, "StaBiCut_v2_Meta_biological_consistency.xlsx"), overwrite = TRUE)

message("\nDone. Outputs saved to: ", out_dir)
