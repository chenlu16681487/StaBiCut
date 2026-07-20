############################################################
# StaBiCut v2 paper analysis
# Script 01: External validation in GSE39582
#
# Purpose:
#   Reproduce the tumor-only external StaBiCut analysis of the 11 candidate
#   genes in GSE39582, including probe annotation, survival-data preparation,
#   direction-constrained cutoff analysis, bootstrap re-support, figures, and
#   supplementary result tables.
#
# Suggested directory layout:
#   paper_analysis/
#     01_external_validation_GSE39582.R
#     data/
#       GSE39582_series_matrix.txt.gz
#       expected_dir_table_11DEGs.csv
#     StaBiCut_v2/
#       modules_core_StaBiCut_v2.R
#       modules_plot_single_StaBiCut_v2.R
#       modules_seed_selection_StaBiCut_v2.R
#       modules_stability_summary_StaBiCut_v2.R
#       Panel_helper_StaBiCut_v2.R
#       run_StaBiCut_v2.R
#
# Output:
#   results/GSE39582_external_validation/
#
# Notes:
#   - All paths are relative to this script by default.
#   - Edit only the configuration block below when using another layout.
#   - The expected-direction and tumor-normal reference values are imported
#     from expected_dir_table_11DEGs.csv; they are not estimated in GSE39582.
############################################################

rm(list = ls())

## =========================
## 0. Paths and configuration
## =========================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

analysis_dir <- get_script_dir()
data_dir     <- file.path(analysis_dir, "data")
stabicut_dir <- file.path(analysis_dir, "StaBiCut_v2")
out_dir      <- file.path(analysis_dir, "results", "GSE39582_external_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Accept a compressed or uncompressed GEO series matrix.
series_file <- file.path(data_dir, "GSE39582_series_matrix.txt.gz")
prior_table_file_config <- file.path(data_dir, "expected_dir_table_11DEGs.csv")

## StaBiCut v2 modules used by the external analysis.
stabicut_core <- file.path(stabicut_dir, "modules_core_StaBiCut_v2.R")
stabicut_plot <- file.path(stabicut_dir, "modules_plot_single_StaBiCut_v2.R")
stabicut_seed <- file.path(stabicut_dir, "modules_seed_selection_StaBiCut_v2.R")
stabicut_sum  <- file.path(stabicut_dir, "modules_stability_summary_StaBiCut_v2.R")
panel_helper  <- file.path(stabicut_dir, "Panel_helper_StaBiCut_v2.R")
stabicut_run  <- file.path(stabicut_dir, "run_StaBiCut_v2.R")

## =========================
## 1. Packages
## =========================
cran_pkgs <- c(
  "survival", "survminer", "ggplot2", "dplyr", "tidyr",
  "patchwork", "openxlsx", "stringr", "tibble", "readr"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## GPL570 annotation packages for the Affymetrix Human Genome U133 Plus 2.0 Array.
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) {
  BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

## =========================
## 2. Load StaBiCut modules
## =========================
files_to_source <- c(
  stabicut_core,
  stabicut_plot,
  stabicut_seed,
  stabicut_sum,
  panel_helper,
  stabicut_run
)

for (f in files_to_source) {
  if (!file.exists(f)) stop("Cannot find StaBiCut module: ", f)
  source(f, encoding = "UTF-8")
}

if (!exists("run_batch_sur_cutpoint_analysis_v2")) {
  stop("StaBiCut main runner was not loaded: run_batch_sur_cutpoint_analysis_v2() not found.")
}

## =========================
## 3. Define genes and prior
## =========================
geneset <- c(
  "SPOCK2", "PYCR1", "CA4", "CES1", "ABCB1", "ZG16",
  "TNXB", "HMCN2", "MEP1A", "SLC37A2", "CHGB"
)

## Discovery-derived direction/TN reference table.
## This file is the official StaBiCut example/prior table exported from the
## TCGA-CRC discovery analysis. GSE39582 is tumor-only, so expected_dir is used
## for the official direction-constrained external cutpoint scan, and TN_log2FC
## is used only for the TCGA-reference-TN sensitivity score below.
resolve_prior_table_file <- function(configured_file, data_dir, out_dir) {
  candidates <- c(
    configured_file,
    file.path(data_dir, "expected_dir_table_11DEGs.csv"),
    file.path(out_dir, "expected_dir_table_11DEGs.csv"),
    file.path(getwd(), "expected_dir_table_11DEGs.csv")
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0) {
    stop(
      "Cannot find expected_dir_table_11DEGs.csv. Tried:\n  ",
      paste(c(
        configured_file,
        file.path(data_dir, "expected_dir_table_11DEGs.csv"),
        file.path(out_dir, "expected_dir_table_11DEGs.csv"),
        file.path(getwd(), "expected_dir_table_11DEGs.csv")
      ), collapse = "\n  ")
    )
  }
  normalizePath(candidates[1], winslash = "/", mustWork = TRUE)
}

prior_table_file <- resolve_prior_table_file(
  configured_file = prior_table_file_config,
  data_dir = data_dir,
  out_dir = out_dir
)
gene_prior_table <- read.csv(prior_table_file, check.names = FALSE, stringsAsFactors = FALSE)

required_prior_cols <- c("Gene", "expected_dir", "TN_log2FC")
missing_prior_cols <- setdiff(required_prior_cols, colnames(gene_prior_table))
if (length(missing_prior_cols) > 0) {
  stop("expected_dir_table_11DEGs.csv is missing required column(s): ",
       paste(missing_prior_cols, collapse = ", "))
}

gene_prior_table <- gene_prior_table[match(geneset, gene_prior_table$Gene), , drop = FALSE]
if (any(is.na(gene_prior_table$Gene))) {
  stop("expected_dir_table_11DEGs.csv does not contain all 11 genes: ",
       paste(geneset[is.na(gene_prior_table$Gene)], collapse = ", "))
}
gene_prior_table$TN_log2FC <- suppressWarnings(as.numeric(gene_prior_table$TN_log2FC))
gene_prior_table$expected_dir <- as.character(gene_prior_table$expected_dir)
message("Loaded discovery prior/TN reference table: ", prior_table_file)

## =========================
## 4. Read GEO series matrix
## =========================
## Resolve common compressed and uncompressed series-matrix filenames.
resolve_series_matrix_file <- function(path) {
  candidates <- unique(c(
    path,
    paste0(path, ".txt"),
    paste0(path, ".gz"),
    paste0(path, ".txt.gz"),
    file.path(dirname(path), "GSE39582_series_matrix"),
    file.path(dirname(path), "GSE39582_series_matrix.txt"),
    file.path(dirname(path), "GSE39582_series_matrix.txt.gz")
  ))

  candidates <- candidates[file.exists(candidates)]
  candidates <- candidates[!file.info(candidates)$isdir]
  if (length(candidates) == 0) {
    stop(
      "Cannot find GSE39582 series matrix. Tried:\n  ",
      paste(unique(c(
        path,
        paste0(path, ".txt"),
        paste0(path, ".gz"),
        paste0(path, ".txt.gz"),
        file.path(dirname(path), "GSE39582_series_matrix"),
        file.path(dirname(path), "GSE39582_series_matrix.txt"),
        file.path(dirname(path), "GSE39582_series_matrix.txt.gz")
      )), collapse = "\n  ")
    )
  }

  ## Prefer the decompressed no-extension file if it exists, then .txt, then .gz.
  preferred <- candidates[order(
    !basename(candidates) == "GSE39582_series_matrix",
    grepl("\\.gz$", candidates, ignore.case = TRUE)
  )]
  normalizePath(preferred[1], winslash = "/", mustWork = TRUE)
}

is_gzip_file <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && all(as.integer(magic) == c(31L, 139L))
}

read_text_lines_auto <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  if (file.info(path)$isdir) stop("Path is a directory, not a file: ", path)
  if (file.access(path, mode = 4) != 0) stop("File is not readable by R: ", path)

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
  file <- resolve_series_matrix_file(file)
  message("Reading GSE39582 series matrix from: ", file)

  x <- read_text_lines_auto(file)

  i1 <- grep("^!series_matrix_table_begin", x)
  i2 <- grep("^!series_matrix_table_end", x)
  if (length(i1) != 1 || length(i2) != 1 || i2 <= i1) {
    stop("Cannot find a valid series matrix table boundary in: ", file)
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

  if (ncol(expr) < 2) stop("Expression table has fewer than two columns after parsing.")
  probe_col <- colnames(expr)[1]
  rownames(expr) <- as.character(expr[[probe_col]])
  expr[[probe_col]] <- NULL
  expr <- as.matrix(expr)
  suppressWarnings(mode(expr) <- "numeric")

  if (anyNA(expr)) {
    warning("NAs were introduced during expression matrix numeric conversion. Check the series matrix table if this is unexpected.")
  }

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

  message("Parsed expression matrix: ", nrow(expr), " probes x ", ncol(expr), " samples")
  message("Parsed clinical annotation: ", nrow(clin), " samples x ", ncol(clin), " columns")

  list(expr_probe = expr, clin_raw = clin, series_file_used = file)
}

series_file <- resolve_series_matrix_file(series_file)
gse <- read_gse_series_matrix(series_file)
expr_probe <- gse$expr_probe
clin_raw <- gse$clin_raw

writeLines(gse$series_file_used, con = file.path(out_dir, "GSE39582_series_matrix_file_used.txt"))
write.csv(clin_raw, file.path(out_dir, "GSE39582_clinical_annotation_raw.csv"), row.names = FALSE)

## =========================
## 5. Probe -> Gene mapping
## =========================
probe_ids <- rownames(expr_probe)

probe_map <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(PROBEID, SYMBOL, .keep_all = TRUE)

write.csv(probe_map, file.path(out_dir, "GSE39582_probe_to_gene_mapping_GPL570.csv"), row.names = FALSE)

## Collapse multiple probes per gene by retaining the probe with the highest
## median expression across all samples. This deterministic rule avoids
## averaging probes with potentially different specificity.
collapse_probe_to_gene <- function(expr_probe, probe_map) {
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
  
  list(expr_gene = expr_gene, selected_probe_map = keep_map)
}

collapsed <- collapse_probe_to_gene(expr_probe, probe_map)
expr_gene <- collapsed$expr_gene
selected_probe_map <- collapsed$selected_probe_map

write.csv(selected_probe_map, file.path(out_dir, "GSE39582_selected_probe_per_gene.csv"), row.names = FALSE)

expr_11 <- expr_gene[intersect(geneset, rownames(expr_gene)), , drop = FALSE]
missing_genes <- setdiff(geneset, rownames(expr_11))

write.csv(
  data.frame(Gene = rownames(expr_11), expr_11, check.names = FALSE),
  file.path(out_dir, "GSE39582_11gene_expression_matrix.csv"),
  row.names = FALSE
)

if (length(missing_genes) > 0) {
  warning("Missing genes after probe mapping: ", paste(missing_genes, collapse = ", "))
}

## =========================
## 6. Auto-detect survival variables
## =========================
clean_clin_names <- function(x) {
  names(x) <- make.names(names(x), unique = TRUE)
  x
}
clin <- clean_clin_names(clin_raw)

## Write a column inventory for troubleshooting. This is useful because
## GEO series-matrix clinical fields are stored as heterogeneous
## !Sample_characteristics_ch1 entries rather than standardized column names.
clin_column_inventory <- data.frame(
  column = colnames(clin),
  non_missing = vapply(clin, function(z) sum(!is.na(z) & nzchar(as.character(z))), integer(1)),
  example_1 = vapply(clin, function(z) {
    z <- unique(as.character(z[!is.na(z) & nzchar(as.character(z))]))
    if (length(z) == 0) NA_character_ else z[1]
  }, character(1)),
  example_2 = vapply(clin, function(z) {
    z <- unique(as.character(z[!is.na(z) & nzchar(as.character(z))]))
    if (length(z) < 2) NA_character_ else z[2]
  }, character(1)),
  stringsAsFactors = FALSE
)
write.csv(clin_column_inventory, file.path(out_dir, "GSE39582_clinical_column_inventory.csv"), row.names = FALSE)

detect_col_by_name <- function(df, patterns) {
  cn <- colnames(df)
  hit <- cn[Reduce(`|`, lapply(patterns, function(p) grepl(p, cn, ignore.case = TRUE)))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

parse_survival_time <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(stringr::str_extract(x, "-?\\d+\\.?\\d*")))
}

parse_event <- function(x) {
  z <- tolower(as.character(x))
  out <- rep(NA_real_, length(z))

  ## explicit text labels first
  out[grepl("dead|death|deceased|event|relapse|recurrence|progress|yes", z)] <- 1
  out[grepl("alive|censored|no event|no relapse|no recurrence|disease free|disease-free|no", z)] <- 0

  ## numeric encodings: keep simple and transparent
  num <- suppressWarnings(as.numeric(stringr::str_extract(z, "-?\\d+\\.?\\d*")))
  out[is.na(out) & !is.na(num)] <- ifelse(num[is.na(out) & !is.na(num)] > 0, 1, 0)
  out
}

score_time_column <- function(x) {
  y <- parse_survival_time(x)
  sum(is.finite(y) & y > 0, na.rm = TRUE)
}

score_event_column <- function(x) {
  y <- parse_event(x)
  sum(!is.na(y), na.rm = TRUE)
}

## GSE39582 uses GEO characteristics columns. Prefer OS-specific fields.
time_name_patterns <- c(
  "os.*time", "overall.*survival.*time", "survival.*time", "survival",
  "follow.*up", "followup", "follow.up", "days.*death", "time.*death",
  "time", "months", "month", "days", "delay"
)
event_name_patterns <- c(
  "os.*event", "overall.*survival.*event", "survival.*event", "vital",
  "death", "dead", "status", "event", "censor", "censored", "relapse", "recurrence"
)

time_candidates <- colnames(clin)[Reduce(`|`, lapply(time_name_patterns, function(p) grepl(p, colnames(clin), ignore.case = TRUE)))]
event_candidates <- colnames(clin)[Reduce(`|`, lapply(event_name_patterns, function(p) grepl(p, colnames(clin), ignore.case = TRUE)))]

## Remove obvious non-survival columns from candidates.
time_candidates <- setdiff(time_candidates, c("sample", "geo_accession", "title", "source_name"))
event_candidates <- setdiff(event_candidates, c("sample", "geo_accession", "title", "source_name"))

if (length(time_candidates) > 0) {
  time_scores <- vapply(time_candidates, function(nm) score_time_column(clin[[nm]]), integer(1))
  time_col <- time_candidates[which.max(time_scores)]
  if (max(time_scores, na.rm = TRUE) == 0) time_col <- NA_character_
} else {
  time_col <- NA_character_
}

if (length(event_candidates) > 0) {
  event_scores <- vapply(event_candidates, function(nm) score_event_column(clin[[nm]]), integer(1))
  event_col <- event_candidates[which.max(event_scores)]
  if (max(event_scores, na.rm = TRUE) == 0) event_col <- NA_character_
} else {
  event_col <- NA_character_
}

message("Detected time column: ", ifelse(is.na(time_col), "<not detected>", time_col))
message("Detected event column: ", ifelse(is.na(event_col), "<not detected>", event_col))

## Manual override if auto-detection is wrong. Uncomment and set after checking
## GSE39582_clinical_column_inventory.csv.
## time_col  <- "your_OS_time_column"
## event_col <- "your_OS_event_column"

if (is.na(time_col) || !time_col %in% colnames(clin)) {
  stop(
    "No usable survival time column was detected. Open this diagnostic file and set time_col manually: ",
    file.path(out_dir, "GSE39582_clinical_column_inventory.csv")
  )
}
if (is.na(event_col) || !event_col %in% colnames(clin)) {
  stop(
    "No usable survival event/status column was detected. Open this diagnostic file and set event_col manually: ",
    file.path(out_dir, "GSE39582_clinical_column_inventory.csv")
  )
}

clin$time_raw <- parse_survival_time(clin[[time_col]])
clin$event <- parse_event(clin[[event_col]])

if (length(clin$time_raw) != nrow(clin) || length(clin$event) != nrow(clin)) {
  stop("Survival parsing returned vectors with unexpected length. Check time_col/event_col.")
}

## Unit detection: if the median time is large, treat as days; otherwise months.
median_time_raw <- stats::median(clin$time_raw[is.finite(clin$time_raw) & clin$time_raw > 0], na.rm = TRUE)
convert_days_to_months <- is.finite(median_time_raw) && median_time_raw > 120
clin$time <- if (convert_days_to_months) clin$time_raw / 30.4375 else clin$time_raw
message("Survival time interpreted as: ", ifelse(convert_days_to_months, "days converted to months", "months"))

surv_parse_summary <- data.frame(
  time_col = time_col,
  event_col = event_col,
  n_samples = nrow(clin),
  n_time_valid = sum(is.finite(clin$time) & clin$time > 0, na.rm = TRUE),
  n_event_valid = sum(!is.na(clin$event), na.rm = TRUE),
  n_events = sum(clin$event == 1, na.rm = TRUE),
  n_censored = sum(clin$event == 0, na.rm = TRUE),
  time_unit_interpreted = ifelse(convert_days_to_months, "days_to_months", "months"),
  stringsAsFactors = FALSE
)
write.csv(surv_parse_summary, file.path(out_dir, "GSE39582_survival_parse_summary.csv"), row.names = FALSE)

clin_surv <- clin %>%
  dplyr::filter(is.finite(time), time > 0, !is.na(event)) %>%
  dplyr::distinct(sample, .keep_all = TRUE)

common_samples <- intersect(colnames(expr_gene), clin_surv$sample)
expr_gene <- expr_gene[, common_samples, drop = FALSE]
clin_surv <- clin_surv[match(common_samples, clin_surv$sample), , drop = FALSE]
rownames(clin_surv) <- clin_surv$sample

message("Valid survival samples: ", nrow(clin_surv))

if (nrow(clin_surv) < 50) {
  stop("Too few valid survival samples after parsing: ", nrow(clin_surv),
       ". Check time_col/event_col using GSE39582_clinical_column_inventory.csv.")
}

write.csv(clin_surv, file.path(out_dir, "GSE39582_clinical_annotation_survival_used.csv"), row.names = FALSE)

## =========================
## 7. External StaBiCut using the official StaBiCut v2 main runner
## =========================
## run_batch_sur_cutpoint_analysis_v2() was originally written for TCGA-like
## expression + clinical objects. To reuse the exact same cutoff-scan,
## bootstrap, Cox and plotting code path for GSE39582, we only adapt the input
## interface here:
##   1) GEO log2 microarray expression is back-transformed to pseudo-linear
##      scale so that the runner's internal log2(x + 1) returns the original
##      log2 expression values.
##   2) GEO sample IDs are temporarily replaced by TCGA-like tumor barcodes
##      with positions 14-15 equal to "01", because the official runner uses
##      substr(colnames(exprset), 14, 15) < 10 to identify tumor samples.
##   3) survival columns are converted to days_to_death,
##      days_to_last_follow_up and vital_status, matching the official runner's
##      expected clinical interface.
##
## No local duplicate of scan_cutpoints_v2(), bootstrap_cutoffs_v2(), or the
## per-gene StaBiCut loop is retained below.

make_tcga_like_ids <- function(n) {
  ## Length before sample-type position is 13; positions 14-15 are "01".
  paste0("GSE39582-", sprintf("%04d", seq_len(n)), "01")
}

geo_to_stabicut_input <- function(expr_gene, clin_surv) {
  common <- intersect(colnames(expr_gene), clin_surv$sample)
  if (length(common) == 0) stop("No overlapping samples between expression and survival clinical table.")

  expr0 <- expr_gene[, common, drop = FALSE]
  clin0 <- clin_surv[match(common, clin_surv$sample), , drop = FALSE]

  pseudo_ids <- make_tcga_like_ids(length(common))

  ## GEO series matrix expression is already log2-like. The official runner
  ## applies log2(exprset + 1), so use 2^x - 1 here to preserve the original scale
  ## after the runner's internal transformation.
  expr_input <- 2^expr0 - 1
  expr_input[expr_input < 0] <- 0
  colnames(expr_input) <- pseudo_ids

  time_days <- as.numeric(clin0$time) * 30.4375
  event <- as.numeric(clin0$event)

  clin_input <- data.frame(
    sample = pseudo_ids,
    original_sample = common,
    days_to_death = ifelse(event == 1, time_days, NA_real_),
    days_to_last_follow_up = ifelse(event == 0, time_days, NA_real_),
    vital_status = ifelse(event == 1, "Dead", "Alive"),
    stringsAsFactors = FALSE
  )
  rownames(clin_input) <- pseudo_ids

  list(expr_input = expr_input, clin_input = clin_input, id_map = clin_input[, c("sample", "original_sample")])
}

stabicut_input <- geo_to_stabicut_input(expr_gene = expr_gene, clin_surv = clin_surv)
write.csv(stabicut_input$id_map, file.path(out_dir, "GSE39582_StaBiCut_synthetic_sample_id_map.csv"), row.names = FALSE)

stabicut_res <- run_batch_sur_cutpoint_analysis_v2(
  exprset = stabicut_input$expr_input,
  geneset = geneset,
  clin = stabicut_input$clin_input,
  gene_prior_table = gene_prior_table,
  force_direction = FALSE,
  n_boot = 1000,
  adjust_method = "BH",
  minprop = 0.25,
  score_threshold = 0,
  save_boot_rds = TRUE,
  boot_dir = file.path(out_dir, "boot_rds"),
  seed = 11,
  save_plots = TRUE,
  plot_dir = file.path(out_dir, "StaBiCut_official_runner_plots")
)

## External GSE39582 has no matched normal samples, so this external validation
## focuses on prognostic cutoff stability rather than tumor-normal consistency.
## The official StaBiCut runner still produces the cutoff scan, bootstrap, Cox,
## KM plotting caches, and hazard-direction consistency using the same code path.
## Here we only add external-specific reporting columns; we do NOT overwrite the
## official composite_score or Rank returned by the runner.

## TCGA discovery ranking from the full StaBiCut v2 framework, including the
## tumor-normal concordance layer. These ranks are carried forward only for
## side-by-side reporting, not recalculated from GSE39582.
discovery_rank_table <- data.frame(
  Gene = c("CA4", "ZG16", "ABCB1", "CES1", "PYCR1", "CHGB",
           "SLC37A2", "TNXB", "HMCN2", "MEP1A", "SPOCK2"),
  Rank_discovery_TCGA = seq_len(11),
  stringsAsFactors = FALSE
)

gene_colors <- c(
  "CA4" = "#280B54FF", "ZG16" = "#0D0887FF",
  "ABCB1" = "#5402A3FF", "CES1" = "#8B0AA5FF",
  "PYCR1" = "#CC4678FF", "CHGB" = "#E97158FF",
  "SLC37A2" = "#FBA139FF", "TNXB" = "#FADA24FF",
  "HMCN2" = "#9FDA3AFF", "MEP1A" = "#2DB27DFF",
  "SPOCK2" = "#40B7ADFF"
)

## Optional TCGA-derived tumor-normal reference layer.
## GSE39582 is tumor-only. Therefore, the primary external score remains
## composite_score_external_no_TN. However, for transparent sensitivity analysis,
## a TCGA-derived tumor-normal shift can be carried forward as a reference TN
## layer and combined with the GSE39582-derived bootstrap/hazard/density/balance
## components. This is NOT a de novo TN estimate in GSE39582; it is explicitly
## labelled as a TCGA-reference-TN sensitivity score.
##
## The TCGA discovery TN reference is imported from expected_dir_table_11DEGs.csv,
## which must contain Gene, expected_dir, and TN_log2FC columns.
load_tcga_tn_reference <- function(gene_prior_table, geneset, source_file) {
  ## Import the TCGA discovery TN_log2FC reference directly from the official
  ## expected_dir_table_11DEGs.csv. No TN value is recomputed from GSE39582.
  ref <- data.frame(
    Gene = geneset,
    expected_dir_ref = gene_prior_table$expected_dir[match(geneset, gene_prior_table$Gene)],
    TCGA_TN_log2FC_ref = gene_prior_table$TN_log2FC[match(geneset, gene_prior_table$Gene)],
    TCGA_tn_score_ref = NA_real_,
    TCGA_TN_source_file = basename(source_file),
    stringsAsFactors = FALSE
  )
  ref
}

calc_tcga_ref_tn_score <- function(tcga_TN_log2FC, external_HR) {
  ## Same sign-agreement logic as calc_direction_concordance():
  ## sign(TN_log2FC) * sign(log(HR)) == 1 is concordant.
  out <- rep(NA_real_, length(tcga_TN_log2FC))
  ok <- is.finite(tcga_TN_log2FC) & is.finite(external_HR) & external_HR > 0
  s1 <- sign(tcga_TN_log2FC[ok])
  s2 <- sign(log(external_HR[ok]))
  concord <- s1 * s2
  out[ok] <- ifelse(concord == 1, 1, ifelse(concord == -1, 0, NA_real_))
  out
}

tcga_tn_reference_table <- load_tcga_tn_reference(
  gene_prior_table = gene_prior_table,
  geneset = geneset,
  source_file = prior_table_file
)

add_external_rankings <- function(results_df,
                                  discovery_rank_table = NULL,
                                  tcga_tn_reference_table = NULL) {
  if (is.null(results_df) || nrow(results_df) == 0) return(results_df)

  ## Keep the official runner output explicit. In GSE39582 this is expected to
  ## be NA when tn_score is NA, because no normal samples are available.
  results_df$composite_score_official_with_TN_formula <- results_df$composite_score
  results_df$Rank_official_with_TN_formula <- results_df$Rank

  results_df$hazard_score <- as.numeric(results_df$hazard_dir_consistency)
  results_df$tn_score <- NA_real_
  results_df$balance_score <- pmin(1, pmax(0, as.numeric(results_df$min_group_prop) * 2))

  w_boot <- 0.40
  w_den  <- 0.20
  w_bal  <- 0.10
  w_haz  <- 0.15
  w_tn   <- 0.15

  denom_no_TN <- w_boot + w_den + w_bal + w_haz

  ## Primary external score: computed only from GSE39582-available layers.
  results_df$composite_score_external_no_TN <- with(
    results_df,
    (w_boot * as.numeric(bootstrap_score) +
       w_haz * as.numeric(hazard_score) +
       w_den  * as.numeric(density_score) +
       w_bal  * as.numeric(balance_score)) / denom_no_TN
  )

  results_df$Rank_external_no_TN <- rank(
    -results_df$composite_score_external_no_TN,
    ties.method = "min",
    na.last = "keep"
  )

  ## Sensitivity score: GSE39582 layers + TCGA-derived reference TN layer.
  if (!is.null(tcga_tn_reference_table)) {
    results_df$TCGA_TN_log2FC_ref <- tcga_tn_reference_table$TCGA_TN_log2FC_ref[
      match(results_df$Gene, tcga_tn_reference_table$Gene)
    ]
    results_df$TCGA_TN_source_file <- tcga_tn_reference_table$TCGA_TN_source_file[
      match(results_df$Gene, tcga_tn_reference_table$Gene)
    ]

    results_df$TCGA_tn_score_ref <- calc_tcga_ref_tn_score(
      tcga_TN_log2FC = results_df$TCGA_TN_log2FC_ref,
      external_HR = as.numeric(results_df$HR)
    )

    ## Full five-layer formula, but TN is explicitly labelled as TCGA-reference.
    results_df$composite_score_external_with_TCGA_ref_TN <- with(
      results_df,
      w_boot * as.numeric(bootstrap_score) +
        w_haz * as.numeric(hazard_score) +
        w_tn  * as.numeric(TCGA_tn_score_ref) +
        w_den * as.numeric(density_score) +
        w_bal * as.numeric(balance_score)
    )

    if (all(is.na(results_df$TCGA_tn_score_ref))) {
      results_df$Rank_external_with_TCGA_ref_TN <- NA_integer_
    } else {
      results_df$Rank_external_with_TCGA_ref_TN <- rank(
        -results_df$composite_score_external_with_TCGA_ref_TN,
        ties.method = "min",
        na.last = "keep"
      )
    }
  }

  if (!is.null(discovery_rank_table)) {
    results_df$Rank_discovery_TCGA <- discovery_rank_table$Rank_discovery_TCGA[
      match(results_df$Gene, discovery_rank_table$Gene)
    ]
  }

  results_df <- results_df[
    order(results_df$Rank_external_no_TN, results_df$P, na.last = TRUE),
    , drop = FALSE
  ]
  results_df
}

stabicut_res$results_df <- add_external_rankings(
  results_df = stabicut_res$results_df,
  discovery_rank_table = discovery_rank_table,
  tcga_tn_reference_table = tcga_tn_reference_table
)

write.csv(stabicut_res$results_df, file.path(out_dir, "GSE39582_StaBiCut_external_results.csv"), row.names = FALSE)
write.csv(discovery_rank_table, file.path(out_dir, "TCGA_discovery_StaBiCut_rank_used.csv"), row.names = FALSE)

## Export the imported TCGA TN reference together with the GSE39582 HR-derived
## TN concordance score used in the sensitivity composite.
tcga_tn_reference_scored <- stabicut_res$results_df %>%
  dplyr::select(Gene, TCGA_TN_log2FC_ref, HR, TCGA_tn_score_ref, TCGA_TN_source_file)
write.csv(tcga_tn_reference_scored, file.path(out_dir, "TCGA_TN_reference_used_for_GSE39582_sensitivity.csv"), row.names = FALSE)

expected_dir_table <- data.frame(
  Gene = geneset,
  expected_dir = gene_prior_table$expected_dir[match(geneset, gene_prior_table$Gene)],
  TN_log2FC = gene_prior_table$TN_log2FC[match(geneset, gene_prior_table$Gene)],
  prior_source_file = basename(prior_table_file),
  stringsAsFactors = FALSE
)
write.csv(expected_dir_table, file.path(out_dir, "GSE39582_expected_direction_and_TCGA_TN_reference_used.csv"), row.names = FALSE)
## Backward-compatible file name for previous downstream scripts.
write.csv(expected_dir_table, file.path(out_dir, "GSE39582_expected_direction_used.csv"), row.names = FALSE)

## =========================
## 8. ZG16 KM
## =========================
zg16_df <- stabicut_res$df_cache[["ZG16"]]
zg16_row <- stabicut_res$results_df %>% dplyr::filter(Gene == "ZG16")

if (!is.null(zg16_df) && nrow(zg16_row) == 1) {
  pdf(file.path(out_dir, "Supplementary_Figure_GSE39582_ZG16_KM.pdf"), width = 5.2, height = 5.2)
  print(
    survminer::ggsurvplot(
      survival::survfit(survival::Surv(time, event) ~ group, data = zg16_df),
      data = zg16_df,
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      legend.title = "ZG16",
      legend.labs = c("Low", "High"),
      xlab = "Months",
      ylab = "Overall survival probability",
      title = paste0(
        "GSE39582 ZG16 KM; HR = ",
        sprintf("%.3f", zg16_row$HR),
        " [", sprintf("%.3f", zg16_row$CI_low), ", ",
        sprintf("%.3f", zg16_row$CI_high), "]"
      )
    )
  )
  dev.off()
}

## =========================
## 9. Univariate and multivariate Cox
## =========================
cox_uni <- lapply(names(stabicut_res$df_cache), function(g) {
  df <- stabicut_res$df_cache[[g]]
  fit <- survival::coxph(survival::Surv(time, event) ~ group, data = df)
  s <- summary(fit)
  data.frame(
    Gene = g,
    HR = s$coefficients[, "exp(coef)"][1],
    CI_low = s$conf.int[, "lower .95"][1],
    CI_high = s$conf.int[, "upper .95"][1],
    P = s$coefficients[, "Pr(>|z|)"][1],
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

write.csv(cox_uni, file.path(out_dir, "GSE39582_univariate_Cox_11genes.csv"), row.names = FALSE)

## Use age and stage as multivariable covariates when valid fields are available.
age_col <- detect_col_by_name(clin_surv, c("^age$", "age"))
stage_col <- detect_col_by_name(clin_surv, c("stage", "tnm"))

cox_multi <- NULL
if (!is.null(zg16_df) && nrow(zg16_df) > 0) {
  multi_df <- zg16_df

  ## Map synthetic StaBiCut IDs back to the original GEO sample IDs for covariates.
  id_map <- stabicut_input$id_map
  multi_df$sample <- id_map$original_sample[match(rownames(multi_df), id_map$sample)]

  if (!is.na(age_col)) {
    multi_df$age <- suppressWarnings(as.numeric(parse_survival_time(clin_surv[[age_col]][match(multi_df$sample, clin_surv$sample)])))
  }
  if (!is.na(stage_col)) {
    multi_df$stage <- as.factor(clin_surv[[stage_col]][match(multi_df$sample, clin_surv$sample)])
  }

  rhs <- "group"
  if ("age" %in% names(multi_df) && sum(is.finite(multi_df$age)) > 20) rhs <- paste(rhs, "+ age")
  if ("stage" %in% names(multi_df) && length(unique(na.omit(multi_df$stage))) > 1) rhs <- paste(rhs, "+ stage")

  fml <- as.formula(paste0("survival::Surv(time, event) ~ ", rhs))
  fit_multi <- tryCatch(survival::coxph(fml, data = multi_df), error = function(e) NULL)

  if (!is.null(fit_multi)) {
    sm <- summary(fit_multi)
    cox_multi <- data.frame(
      Variable = rownames(sm$coefficients),
      HR = sm$coefficients[, "exp(coef)"],
      CI_low = sm$conf.int[, "lower .95"],
      CI_high = sm$conf.int[, "upper .95"],
      P = sm$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
}

if (!is.null(cox_multi)) {
  write.csv(cox_multi, file.path(out_dir, "GSE39582_ZG16_multivariate_Cox.csv"), row.names = FALSE)
}

## =========================
## 10. Supplementary Figure
## =========================
plot_df <- stabicut_res$results_df %>%
  dplyr::arrange(Rank_external_no_TN) %>%
  dplyr::mutate(Gene = factor(Gene, levels = rev(Gene)))

p_rank <- ggplot(plot_df, aes(x = composite_score_external_no_TN, y = Gene, fill = composite_score_external_no_TN)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", composite_score_external_no_TN)), hjust = -0.1, size = 3.5) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  labs(x = "External StaBiCut score without TN layer", y = NULL, title = "GSE39582 external score without TN") +
  coord_cartesian(xlim = c(0, max(plot_df$composite_score_external_no_TN, na.rm = TRUE) * 1.15))

has_tcga_ref_tn <- "composite_score_external_with_TCGA_ref_TN" %in% colnames(plot_df) &&
  any(is.finite(plot_df$composite_score_external_with_TCGA_ref_TN))

if (has_tcga_ref_tn) {
  plot_df_ref <- plot_df %>%
    dplyr::arrange(Rank_external_with_TCGA_ref_TN) %>%
    dplyr::mutate(Gene = factor(as.character(Gene), levels = rev(as.character(Gene))))

  p_rank_ref <- ggplot(plot_df_ref, aes(x = composite_score_external_with_TCGA_ref_TN, y = Gene,
                                        fill = composite_score_external_with_TCGA_ref_TN)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", composite_score_external_with_TCGA_ref_TN)), hjust = -0.1, size = 3.5) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    theme_bw(base_size = 12) +
    guides(fill = "none") +
    labs(x = "External score with TCGA-reference TN layer", y = NULL,
         title = "GSE39582 score with TCGA-reference TN") +
    coord_cartesian(xlim = c(0, max(plot_df_ref$composite_score_external_with_TCGA_ref_TN, na.rm = TRUE) * 1.15))
}

p_hr <- ggplot(plot_df, aes(y = Gene, x = HR, xmin = CI_low, xmax = CI_high, color = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange() +
  scale_x_log10() +
  scale_color_manual(values = gene_colors) +
  theme_bw(base_size = 12) +
  guides(color = "none") +
  labs(x = "Hazard ratio, High vs Low", y = NULL, title = "Univariate Cox")

supp_fig_no_tn <- p_rank | p_hr

ggsave(
  file.path(out_dir, "Supplementary_Figure_GSE39582_11gene_StaBiCut_external_validation.pdf"),
  supp_fig_no_tn,
  width = 11,
  height = 5.8
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_GSE39582_11gene_StaBiCut_Univariate_Cox.pdf"),
  p_hr + labs(title = NULL),
  width = 6,
  height = 5.8
)

if (has_tcga_ref_tn) {
  supp_fig_with_ref_tn <- p_rank | p_rank_ref | p_hr
  ggsave(
    file.path(out_dir, "Supplementary_Figure_GSE39582_11gene_StaBiCut_external_validation_with_TCGA_reference_TN.pdf"),
    supp_fig_with_ref_tn,
    width = 16.5,
    height = 5.8
  )

  rank_comp <- stabicut_res$results_df %>%
    dplyr::select(Gene, Rank_discovery_TCGA, Rank_external_no_TN, Rank_external_with_TCGA_ref_TN,
                  composite_score_external_no_TN, composite_score_external_with_TCGA_ref_TN,
                  TCGA_TN_log2FC_ref, TCGA_tn_score_ref, HR, P) %>%
    tidyr::pivot_longer(
      cols = c(Rank_discovery_TCGA, Rank_external_no_TN, Rank_external_with_TCGA_ref_TN),
      names_to = "ranking_type",
      values_to = "rank"
    ) %>%
    dplyr::filter(!is.na(rank)) %>%
    dplyr::mutate(
      ranking_type = factor(
        ranking_type,
        levels = c("Rank_discovery_TCGA", "Rank_external_no_TN", "Rank_external_with_TCGA_ref_TN"),
        labels = c("TCGA discovery", "GSE39582 no TN", "GSE39582 + TCGA-ref TN")
      )
    )

  p_rank_shift <- ggplot(rank_comp, aes(x = ranking_type, y = rank, group = Gene, label = Gene)) +
    geom_line(alpha = 0.55) +
    geom_point(size = 2) +
    geom_text(size = 3, hjust = -0.05, check_overlap = TRUE) +
    scale_y_reverse(breaks = seq_len(11)) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "Rank", title = "Ranking comparison across discovery and external scoring")

  ggsave(
    file.path(out_dir, "Supplementary_Figure_GSE39582_rank_shift_TCGA_external.pdf"),
    p_rank_shift,
    width = 8.5,
    height = 5.8
  )
}

## ZG16 4-panel generated directly from official StaBiCut cached outputs
if (exists("plot_gene_panel_main_v2") && "ZG16" %in% names(stabicut_res$df_cache) && nrow(zg16_row) == 1) {
  plot_gene_panel_main_v2(
    df = stabicut_res$df_cache[["ZG16"]],
    gene = "ZG16",
    best_cut = zg16_row$Cutoff[1],
    hr = zg16_row$HR[1],
    ci_low = zg16_row$CI_low[1],
    ci_high = zg16_row$CI_high[1],
    cutoff_p = zg16_row$Cutoff_P[1],
    boot = stabicut_res$boot_cache[["ZG16"]],
    scan_df = stabicut_res$scan_cache[["ZG16"]],
    expected_dir = gene_prior_table$expected_dir[match("ZG16", gene_prior_table$Gene)],
    plot_dir = out_dir,
    width = 12,
    height = 10,
    xps_aspect = 0.625,
    km_aspect = 1
  )
}

## Multi-gene sheets from official StaBiCut cached outputs
## The original one-page 9/11-gene sheet is too compressed for reading. Keep ZG16
## as an individual 4-panel figure above, and split the remaining candidate genes
## into two 5-gene supplementary sheets with larger per-row panel sizes.
if (exists("plot_multi_gene_sheet_main_v2") && nrow(stabicut_res$results_df) > 0) {
  plot_results_df <- stabicut_res$results_df

  ## The plotting helper expects composite_score/Rank columns. For GSE39582,
  ## use the external no-TN score/rank only for plotting labels/order; the
  ## original official columns are still preserved in the exported table.
  plot_results_df$composite_score <- plot_results_df$composite_score_external_no_TN
  plot_results_df$Rank <- plot_results_df$Rank_external_no_TN

  plot_multigene_group <- function(gene_order, out_pdf, width = 16, height = 14) {
    gene_order <- intersect(gene_order, names(stabicut_res$df_cache))
    if (length(gene_order) == 0) {
      warning("No genes available for multi-gene sheet: ", out_pdf)
      return(invisible(NULL))
    }

    plot_multi_gene_sheet_main_v2(
      gene_order = gene_order,
      results_df = plot_results_df,
      df_cache = stabicut_res$df_cache,
      boot_cache = stabicut_res$boot_cache,
      scan_cache = stabicut_res$scan_cache,
      expected_dir_table = expected_dir_table,
      out_pdf = file.path(out_dir, out_pdf),
      width = width,
      height = height,
      show_row_tags = TRUE
    )
  }

  ## User-defined layout: ZG16 is not repeated here because it is already saved
  ## separately as a dedicated single-gene 4-panel figure.
  multigene_group1 <- c("CA4", "ABCB1", "CES1", "PYCR1", "CHGB")
  multigene_group2 <- c("SLC37A2", "TNXB", "HMCN2", "MEP1A", "SPOCK2")

  plot_multigene_group(
    gene_order = multigene_group1,
    out_pdf = "Supplementary_Figure_GSE39582_StaBiCut_multigene_sheet_part1_CA4_ABCB1_CES1_PYCR1_CHGB.pdf",
    width = 16,
    height = 14
  )

  plot_multigene_group(
    gene_order = multigene_group2,
    out_pdf = "Supplementary_Figure_GSE39582_StaBiCut_multigene_sheet_part2_SLC37A2_TNXB_HMCN2_MEP1A_SPOCK2.pdf",
    width = 16,
    height = 14
  )

  ## Optional archive version: keep a single combined PDF, but make it much
  ## taller than before so it is readable if opened separately.
  plot_multigene_group(
    gene_order = c(multigene_group1, multigene_group2),
    out_pdf = "Supplementary_Figure_GSE39582_StaBiCut_multigene_sheet_10genes_readable_archive.pdf",
    width = 16,
    height = 27
  )
}

## =========================
## 11. Supplementary Table Excel
## =========================
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "StaBiCut_results")
openxlsx::writeDataTable(wb, "StaBiCut_results", stabicut_res$results_df)

openxlsx::addWorksheet(wb, "Probe_mapping_selected")
openxlsx::writeDataTable(wb, "Probe_mapping_selected", selected_probe_map)

openxlsx::addWorksheet(wb, "11gene_expression")
openxlsx::writeDataTable(
  wb, "11gene_expression",
  data.frame(Gene = rownames(expr_11), expr_11, check.names = FALSE)
)

openxlsx::addWorksheet(wb, "Clinical_used")
openxlsx::writeDataTable(wb, "Clinical_used", clin_surv)

openxlsx::addWorksheet(wb, "StaBiCut_ID_map")
openxlsx::writeDataTable(wb, "StaBiCut_ID_map", stabicut_input$id_map)

openxlsx::addWorksheet(wb, "Expected_direction")
openxlsx::writeDataTable(wb, "Expected_direction", expected_dir_table)

openxlsx::addWorksheet(wb, "TCGA_discovery_rank")
openxlsx::writeDataTable(wb, "TCGA_discovery_rank", discovery_rank_table)

openxlsx::addWorksheet(wb, "Univariate_Cox")
openxlsx::writeDataTable(wb, "Univariate_Cox", cox_uni)

if (!is.null(cox_multi)) {
  openxlsx::addWorksheet(wb, "ZG16_multivariate_Cox")
  openxlsx::writeDataTable(wb, "ZG16_multivariate_Cox", cox_multi)
}

openxlsx::addWorksheet(wb, "Missing_genes")
openxlsx::writeDataTable(wb, "Missing_genes", data.frame(Missing_gene = missing_genes))

openxlsx::saveWorkbook(
  wb,
  file.path(out_dir, "Supplementary_Table_GSE39582_StaBiCut_external_validation.xlsx"),
  overwrite = TRUE
)

saveRDS(stabicut_res, file.path(out_dir, "GSE39582_StaBiCut_external_validation_result_object.rds"))

message("Done. Outputs saved to: ", out_dir)
