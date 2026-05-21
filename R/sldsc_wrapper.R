# Stratified LD Score Regression (S-LDSC) post-processing wrappers around polyfun.
#
# This file provides the post-processing layer for the xqtl-protocol sLDSC pipeline:
# read polyfun outputs per trait, compute Gazal-style standardized tau* and the
# differential per-SNP heritability statistic (EnrichStat), and run DerSimonian-Laird
# random-effects meta-analysis across traits.
#
# Reference panel convention: all LD-derived quantities (baseline LD scores,
# target LD scores, regression weights, allele frequencies) must come from the
# same reference panel. Do not mix files from different panels (e.g. 1000G vs ADSP).
#
# MAF convention: by default we restrict to MAF > 5% per the sLDSC recommendation.
# Pass maf_cutoff = 0 to opt out (not recommended).
#
# Cross-type comparison: tau* (Gazal et al. 2017 standardization) is the
# cross-type comparable statistic. Use tau* to rank or meta-analyze annotations
# that mix binary and continuous types. E (proportion-based enrichment) is
# scale-dependent for continuous annotations and is only comparable within type.

# ---- internal helpers ----

.sldsc_std_cols <- c("CHR", "SNP", "BP", "CM", "A1", "A2", "MAF")

.sldsc_chrom_from_filename <- function(f) {
  bn <- basename(f)
  m  <- regmatches(bn, regexec("\\.([0-9]+)\\.annot\\.gz$", bn))[[1]]
  if (length(m) >= 2) as.integer(m[2]) else NA_integer_
}

.sldsc_detect_annot_cols <- function(file_path) {
  sample <- vroom(file_path, n_max = 5L, show_col_types = FALSE)
  setdiff(names(sample), .sldsc_std_cols)
}


#' @title Read S-LDSC outputs from polyfun for one trait/run
#'
#' @description Reads the regression outputs produced by `polyfun/ldsc.py` for a
#'   single polyfun run (one trait, one annotation set) and returns them as a
#'   tidy list ready for downstream standardization. Hides the underlying file
#'   formats; downstream code consumes only modeling quantities.
#'
#' @param prefix Character. Path prefix to the polyfun outputs for one trait/run.
#'   The function appends `.results`, `.log`, and `.part_delete` to this prefix.
#'   Example: `"/path/to/cwd/CAD_META.filtered.sumstats.gz"`.
#'
#' @return A named list. See `sldsc_postprocessing_pipeline` for components.
#'
#' @examples
#' \dontrun{
#' run <- read_sldsc_trait("/output/CAD_META.filtered.sumstats.gz")
#' run$tau["my_target_annotation"]
#' }
#'
#' @importFrom stats setNames var na.omit
#' @importFrom utils head
#' @importFrom vroom vroom
#' @export
read_sldsc_trait <- function(prefix) {
  results_file <- paste0(prefix, ".results")
  log_file     <- paste0(prefix, ".log")
  delete_file  <- paste0(prefix, ".part_delete")

  for (f in c(results_file, log_file, delete_file)) {
    if (!file.exists(f)) stop("read_sldsc_trait: missing file: ", f)
  }

  results <- vroom(results_file, show_col_types = FALSE)
  cats <- as.character(results$Category)

  log_lines <- readLines(log_file, warn = FALSE)
  h2_line <- grep("Total Observed scale h2:", log_lines, value = TRUE)
  if (length(h2_line) == 0L)
    stop("read_sldsc_trait: could not find 'Total Observed scale h2:' in ", log_file)
  h2g <- suppressWarnings(as.numeric(gsub(".*h2: (-?[0-9.eE+-]+).*", "\\1", h2_line[1])))
  if (is.na(h2g))
    stop("read_sldsc_trait: failed to parse h2g numeric from log line: ", h2_line[1])

  delete_values <- as.matrix(vroom(delete_file, show_col_types = FALSE))
  if (ncol(delete_values) != length(cats)) {
    stop("read_sldsc_trait: .part_delete has ", ncol(delete_values),
         " columns but .results has ", length(cats), " categories.")
  }
  colnames(delete_values) <- cats

  list(
    categories     = cats,
    tau            = setNames(as.numeric(results$Coefficient),                 cats),
    tau_se         = setNames(as.numeric(results[["Coefficient_std_error"]]),  cats),
    enrichment     = setNames(as.numeric(results$Enrichment),                  cats),
    enrichment_se  = setNames(as.numeric(results[["Enrichment_std_error"]]),   cats),
    enrichment_p   = setNames(as.numeric(results[["Enrichment_p"]]),           cats),
    prop_h2        = setNames(as.numeric(results[["Prop._h2"]]),               cats),
    prop_snps      = setNames(as.numeric(results[["Prop._SNPs"]]),             cats),
    h2g            = h2g,
    tau_blocks     = delete_values,
    n_blocks       = nrow(delete_values)
  )
}


#' @title Compute per-annotation standard deviation, MAF-restricted
#'
#' @description Computes the standard deviation of each annotation column in the
#'   target annotation files, restricted to SNPs above a MAF cutoff via PLINK
#'   `.frq` files. Required for internal consistency with polyfun's regression,
#'   which operates on MAF > cutoff SNPs by default.
#'
#' @param target_anno_dir Character. Directory containing target annotation files
#'   (one per chromosome) in polyfun's `.annot.gz` format.
#' @param frqfile_dir Character or NULL. Directory containing PLINK `.frq` files
#'   for the reference panel. Required when `maf_cutoff > 0`; the function
#'   errors if missing.
#' @param plink_name Character. Filename prefix of the `.frq` files
#'   (e.g. `"ADSP_chr"`). Files are expected at `{plink_name}{chr}.frq`.
#' @param maf_cutoff Numeric, default `0.05`.
#' @param annot_cols Character or integer vector, default NULL. Annotation columns
#'   to compute sd for. If NULL, all annotation columns are used.
#'
#' @return Named numeric vector of \eqn{sd_C} values, one per annotation.
#'
#' @importFrom stats setNames var
#' @export
compute_sldsc_annot_sd <- function(target_anno_dir, frqfile_dir = NULL,
                                   plink_name = "ADSP_chr",
                                   maf_cutoff = 0.05, annot_cols = NULL) {
  if (maf_cutoff > 0 && (is.null(frqfile_dir) || !dir.exists(frqfile_dir))) {
    stop("compute_sldsc_annot_sd: maf_cutoff = ", maf_cutoff,
         " requires frqfile_dir, but '", frqfile_dir, "' is not a directory.")
  }
  if (!dir.exists(target_anno_dir)) {
    stop("compute_sldsc_annot_sd: target_anno_dir does not exist: ", target_anno_dir)
  }

  anno_files <- list.files(target_anno_dir, pattern = "\\.annot\\.gz$", full.names = TRUE)
  if (length(anno_files) == 0L)
    stop("compute_sldsc_annot_sd: no .annot.gz files in: ", target_anno_dir)

  detected <- .sldsc_detect_annot_cols(anno_files[1])
  if (is.null(annot_cols)) {
    cols_use <- detected
  } else if (is.numeric(annot_cols)) {
    cols_use <- detected[annot_cols]
  } else {
    cols_use <- annot_cols
  }
  if (length(cols_use) == 0L)
    stop("compute_sldsc_annot_sd: no annotation columns to process.")

  num <- setNames(numeric(length(cols_use)), cols_use)
  den <- 0

  for (anno_file in anno_files) {
    dat <- vroom(anno_file, show_col_types = FALSE)
    if (maf_cutoff > 0) {
      chrom <- .sldsc_chrom_from_filename(anno_file)
      if (is.na(chrom))
        stop("compute_sldsc_annot_sd: could not parse chromosome from: ", anno_file)
      frq_file <- file.path(frqfile_dir, paste0(plink_name, chrom, ".frq"))
      if (!file.exists(frq_file))
        stop("compute_sldsc_annot_sd: .frq file not found: ", frq_file)
      frq <- vroom(frq_file, col_select = c("SNP", "MAF"), show_col_types = FALSE)
      dat <- merge(dat, frq, by = "SNP", all.x = FALSE, all.y = FALSE)
      dat <- dat[!is.na(dat$MAF) & dat$MAF > maf_cutoff, ]
    }
    if (nrow(dat) <= 1L) next
    n_minus_1 <- nrow(dat) - 1L
    for (col in cols_use) {
      vals <- as.numeric(dat[[col]])
      v <- var(vals, na.rm = TRUE)
      if (!is.na(v)) num[col] <- num[col] + n_minus_1 * v
    }
    den <- den + n_minus_1
  }

  if (den <= 0)
    stop("compute_sldsc_annot_sd: zero degrees of freedom after MAF filtering.")
  sqrt(num / den)
}


#' @title Reference-panel SNP count (the M_ref used to standardise tau*)
#'
#' @description `M_ref` is the number of SNPs in the REFERENCE PANEL over which
#'   heritability is partitioned in the sLDSC model
#'   (`h2(C) = sum_{j in M_ref} a_C(j) sum_{C'} tau_{C'} a_{C'}(j)`). It is
#'   panel-defined and is **not** the regression SNP set (HapMap3 ~1M) nor any
#'   HM3-subsetted target output:
#'   \itemize{
#'     \item `maf_cutoff > 0` (Gazal/Finucane convention): count MAF > cutoff
#'       SNPs across all `.frq` files (the same set polyfun's `.l2.M_5_50` sums).
#'     \item `maf_cutoff == 0` (all-M variant): count ALL SNPs across all
#'       `.frq` files (the same set polyfun's `.l2.M` sums).
#'   }
#'   `target_anno_dir` is a fallback only, used when no `.frq` directory is
#'   given; that fallback counts `.l2.ldscore` rows and is WRONG when the target
#'   was HM3-subsetted (it then yields the regression SNP count, not M_ref).
#'
#' @param target_anno_dir Character or NULL. Fallback only - directory of
#'   `.l2.ldscore` files. Used only when `frqfile_dir` is unavailable.
#' @param frqfile_dir Character or NULL. Directory of PLINK `.frq` files; the
#'   preferred (recommended) source of M_ref.
#' @param plink_name Character. Filename prefix of `.frq` files.
#' @param maf_cutoff Numeric, default `0.05`.
#'
#' @return Scalar integer.
#'
#' @export
compute_sldsc_M_ref <- function(target_anno_dir = NULL, frqfile_dir = NULL,
                                plink_name = "ADSP_chr", maf_cutoff = 0.05) {
  ## --- preferred path: count reference-panel SNPs from the .frq files ---
  if (!is.null(frqfile_dir) && dir.exists(frqfile_dir)) {
    pat <- paste0("^", gsub("([.])", "\\\\\\1", plink_name), "[0-9]+\\.frq$")
    frq_files <- list.files(frqfile_dir, pattern = pat, full.names = TRUE)
    if (length(frq_files) == 0L)
      frq_files <- list.files(frqfile_dir, pattern = "\\.frq$", full.names = TRUE)
    if (length(frq_files) > 0L) {
      total <- 0L
      for (f in frq_files) {
        frq <- vroom(f, col_select = "MAF", show_col_types = FALSE)
        total <- total + if (maf_cutoff > 0)
          sum(!is.na(frq$MAF) & frq$MAF > maf_cutoff) else nrow(frq)
      }
      return(as.integer(total))
    }
  }

  ## --- fallback only (no .frq dir): count target .l2.ldscore rows ---
  if (maf_cutoff > 0)
    stop("compute_sldsc_M_ref: maf_cutoff = ", maf_cutoff,
         " requires frqfile_dir (to count MAF>cutoff reference-panel SNPs).")
  if (is.null(target_anno_dir) || !dir.exists(target_anno_dir))
    stop("compute_sldsc_M_ref: need frqfile_dir, or target_anno_dir as fallback.")
  files <- list.files(target_anno_dir,
                      pattern = "\\.l2\\.ldscore\\.(gz|parquet)$",
                      full.names = TRUE)
  if (length(files) == 0L)
    stop("compute_sldsc_M_ref: no .frq files and no .l2.ldscore files found.")
  warning("compute_sldsc_M_ref: no .frq dir given; counting target .l2.ldscore ",
          "rows as M_ref. If the target was HM3-subsetted this UNDERCOUNTS the ",
          "reference panel and shrinks tau*. Pass frqfile_dir instead.")
  total <- 0L
  for (f in files) {
    if (endsWith(f, ".parquet")) {
      if (!requireNamespace("arrow", quietly = TRUE))
        stop("compute_sldsc_M_ref: install 'arrow' to read .parquet files.")
      total <- total + nrow(arrow::read_parquet(f))
    } else {
      total <- total + nrow(vroom(f, show_col_types = FALSE))
    }
  }
  as.integer(total)
}


#' @title Detect whether each annotation is binary or continuous
#'
#' @description Inspects each annotation column and returns whether its values
#'   lie in \{0, 1\} (binary) or take other values (continuous).
#'
#' @param target_anno_dir Character. Directory containing the target `.annot.gz`
#'   files (one per chromosome).
#' @param annot_cols Character or integer vector, default NULL.
#'
#' @return Named logical vector: TRUE for binary, FALSE for continuous.
#'
#' @importFrom stats setNames
#' @export
is_binary_sldsc_annot <- function(target_anno_dir, annot_cols = NULL) {
  anno_files <- list.files(target_anno_dir, pattern = "\\.annot\\.gz$", full.names = TRUE)
  if (length(anno_files) == 0L)
    stop("is_binary_sldsc_annot: no .annot.gz files in: ", target_anno_dir)

  detected <- .sldsc_detect_annot_cols(anno_files[1])
  if (is.null(annot_cols)) {
    cols_use <- detected
  } else if (is.numeric(annot_cols)) {
    cols_use <- detected[annot_cols]
  } else {
    cols_use <- annot_cols
  }

  is_binary <- setNames(rep(TRUE, length(cols_use)), cols_use)

  for (f in anno_files) {
    dat <- vroom(f, col_select = all_of(cols_use), show_col_types = FALSE)
    for (col in cols_use) {
      if (!is_binary[[col]]) next
      vals <- unique(na.omit(as.numeric(dat[[col]])))
      if (any(!(vals %in% c(0, 1)))) is_binary[[col]] <- FALSE
    }
    if (!any(is_binary)) break
  }

  is_binary
}


#' @title Standardize tau and compute EnrichStat for one polyfun run
#'
#' @description Applies the Gazal standardization
#'   \eqn{\tau^*_C = \tau_C \cdot sd_C \cdot M_{ref} / h^2_g} to the point and
#'   to each jackknife block. For `mode = "single"`, additionally computes
#'   EnrichStat and back-solves its standard error from polyfun's reported
#'   `Enrichment_p` using \eqn{|Z| = \Phi^{-1}(1 - p/2)}.
#'
#' @param trait_data List from \code{\link{read_sldsc_trait}}.
#' @param sd_annot Named numeric vector from \code{\link{compute_sldsc_annot_sd}}.
#' @param M_ref Scalar from \code{\link{compute_sldsc_M_ref}}.
#' @param target_categories Character vector or NULL. If NULL, intersects
#'   `trait_data$categories` with `names(sd_annot)`.
#' @param mode Character: `"single"` or `"joint"`.
#'
#' @return A list with `summary` (data frame), `tau_star_blocks` (matrix),
#'   `h2g`, `n_blocks`, `mode`.
#'
#' @importFrom stats qnorm var
#' @export
standardize_sldsc_trait <- function(trait_data, sd_annot, M_ref,
                                    target_categories = NULL,
                                    mode = c("single", "joint")) {
  mode <- match.arg(mode)
  if (is.null(target_categories))
    target_categories <- intersect(trait_data$categories, names(sd_annot))
  if (length(target_categories) == 0L)
    stop("standardize_sldsc_trait: no target categories.")

  target_idx <- match(target_categories, trait_data$categories)
  if (any(is.na(target_idx)))
    stop("standardize_sldsc_trait: missing categories: ",
         paste(target_categories[is.na(target_idx)], collapse = ", "))

  h2g <- trait_data$h2g
  sd_target <- as.numeric(sd_annot[target_categories])
  if (any(is.na(sd_target) | sd_target == 0))
    warning("standardize_sldsc_trait: zero/NA sd for some targets; tau* will be NA/0.")

  tau        <- as.numeric(trait_data$tau[target_categories])
  tau_se     <- as.numeric(trait_data$tau_se[target_categories])
  blocks_target   <- trait_data$tau_blocks[, target_idx, drop = FALSE]

  ts <- standardize_tau_star(tau, blocks_target, sd_target, M_ref, h2g)
  tau_star    <- ts$tau_star
  tau_star_se <- ts$tau_star_se

  tau_star_blocks <- sweep(blocks_target, 2L, sd_target * M_ref / h2g, FUN = "*")

  summary_df <- data.frame(
    target      = target_categories,
    tau         = tau,
    tau_se      = tau_se,
    tau_star    = tau_star,
    tau_star_se = tau_star_se,
    stringsAsFactors = FALSE
  )

  if (mode == "single") {
    enrich    <- as.numeric(trait_data$enrichment[target_categories])
    enrich_se <- as.numeric(trait_data$enrichment_se[target_categories])
    enrich_p  <- as.numeric(trait_data$enrichment_p[target_categories])
    p_h2      <- as.numeric(trait_data$prop_h2[target_categories])
    p_M       <- as.numeric(trait_data$prop_snps[target_categories])

    diff_ratio  <- (p_h2 / p_M) - (1 - p_h2) / (1 - p_M)
    enrichstat  <- (h2g / M_ref) * diff_ratio

    abs_z <- qnorm(1 - enrich_p / 2)
    enrichstat_se <- abs(enrichstat) / abs_z
    enrichstat_se[!is.finite(abs_z) | abs_z <= 0] <- NA_real_

    summary_df$enrichment    <- enrich
    summary_df$enrichment_se <- enrich_se
    summary_df$enrichment_p  <- enrich_p
    summary_df$enrichstat    <- enrichstat
    summary_df$enrichstat_se <- enrichstat_se
  }

  list(
    summary         = summary_df,
    tau_star_blocks = tau_star_blocks,
    h2g             = h2g,
    n_blocks        = nrow(blocks_target),
    mode            = mode
  )
}


#' @title Random-effects meta-analysis of S-LDSC quantities across traits
#'
#' @description DerSimonian-Laird random-effects meta-analysis of one S-LDSC
#'   quantity for one annotation across multiple traits.
#'
#' @details Per-trait \eqn{SE_i} sources:
#'   - `quantity = "tau_star"`: jackknife SE from per-block \eqn{\tau^*}.
#'   - `quantity = "enrichment"`: polyfun-reported `Enrichment_std_error`.
#'   - `quantity = "enrichstat"`: back-solved SE from polyfun's `Enrichment_p`.
#'
#' @param per_trait_estimates Named list of per-trait results (each with a
#'   `summary` data frame).
#' @param category Character. Annotation name to meta-analyze.
#' @param quantity Character: `"tau_star"`, `"enrichment"`, or `"enrichstat"`.
#'
#' @return List with `mean`, `se`, `p`, `n_traits`, `traits_used`, `tau2`.
#'
#' @importFrom stats pnorm
#' @export
meta_sldsc_random <- function(per_trait_estimates, category,
                              quantity = c("tau_star", "enrichment", "enrichstat")) {
  quantity <- match.arg(quantity)
  col_pairs <- list(
    tau_star   = c("tau_star",   "tau_star_se"),
    enrichment = c("enrichment", "enrichment_se"),
    enrichstat = c("enrichstat", "enrichstat_se")
  )
  cols <- col_pairs[[quantity]]
  trait_names <- names(per_trait_estimates)
  if (is.null(trait_names))
    trait_names <- as.character(seq_along(per_trait_estimates))

  means <- numeric(0); ses <- numeric(0); used <- character(0)
  for (i in seq_along(per_trait_estimates)) {
    pt <- per_trait_estimates[[i]]
    if (is.null(pt) || is.null(pt$summary)) next
    df <- pt$summary
    row <- df[df$target == category, , drop = FALSE]
    if (nrow(row) == 0L) next
    if (!all(cols %in% names(row))) next
    m <- as.numeric(row[[cols[1]]])[1]
    s <- as.numeric(row[[cols[2]]])[1]
    if (is.na(m) || is.na(s) || !is.finite(s) || s <= 0) next
    means <- c(means, m); ses <- c(ses, s); used <- c(used, trait_names[i])
  }

  if (length(means) < 2L) {
    return(list(mean = NA_real_, se = NA_real_, p = NA_real_,
                n_traits = length(means), traits_used = used,
                tau2 = NA_real_))
  }
  meta <- meta_random_effects(means, ses)
  z    <- meta$mean / meta$se
  p    <- 2 * pnorm(-abs(z))
  list(
    mean        = meta$mean,
    se          = meta$se,
    p           = as.numeric(p),
    n_traits    = length(means),
    traits_used = used,
    tau2        = meta$tau2
  )
}


# Internal helper: assemble a wide per-trait summary frame with single + joint
# columns side by side.
.sldsc_assemble_trait_summary <- function(single_df, joint_df, target_categories,
                                          is_binary_vec) {
  rows <- if (!is.null(single_df)) single_df$target else
          if (!is.null(joint_df))  joint_df$target  else target_categories
  out <- data.frame(target = rows,
                    is_binary = unname(is_binary_vec[rows]),
                    stringsAsFactors = FALSE)

  add_cols <- function(out, src, suffix) {
    cols_to_add <- c("tau", "tau_se", "tau_star", "tau_star_se",
                     "enrichment", "enrichment_se", "enrichment_p",
                     "enrichstat", "enrichstat_se")
    for (c in cols_to_add) {
      newcol <- paste0(c, "_", suffix)
      if (!is.null(src) && c %in% names(src)) {
        out[[newcol]] <- src[[c]][match(out$target, src$target)]
      } else {
        out[[newcol]] <- NA_real_
      }
    }
    out
  }
  out <- add_cols(out, single_df, "single")
  out <- add_cols(out, joint_df,  "joint")
  out
}


# Internal helper: build a per-trait list view that meta_sldsc_random can read.
# Each list element has a $summary frame with the requested mode's columns
# renamed to the canonical names (tau_star, tau_star_se, enrichment, ...).
.sldsc_view_for_meta <- function(per_trait, suffix) {
  lapply(per_trait, function(pt) {
    if (is.null(pt$summary)) return(NULL)
    df <- pt$summary
    cols_have <- c("tau_star", "tau_star_se", "enrichment", "enrichment_se",
                   "enrichment_p", "enrichstat", "enrichstat_se")
    src_cols <- paste0(cols_have, "_", suffix)
    avail    <- src_cols %in% names(df)
    if (!any(avail)) return(NULL)
    new_df <- data.frame(target = df$target, stringsAsFactors = FALSE)
    for (k in seq_along(cols_have)) {
      if (avail[k]) new_df[[cols_have[k]]] <- df[[src_cols[k]]]
    }
    list(summary = new_df)
  })
}


#' @title End-to-end S-LDSC post-processing across traits, single + joint in one pass
#'
#' @description Top-level orchestration. Reads polyfun outputs (one single-target
#'   run per target plus, when available, one joint run per trait), standardizes
#'   both modes, and runs the default random-effects meta across all traits.
#'
#' @param trait_single_prefixes Named list. For each trait, a character vector
#'   of length \eqn{N} giving the polyfun output prefixes for the \eqn{N}
#'   single-target runs (order must match `target_categories`).
#' @param trait_joint_prefix Named character. For each trait, the polyfun output
#'   prefix for the joint run. Pass `NA` (or `""`) for a trait without a joint run.
#' @param target_anno_dir Character. Directory of target `.annot.gz` files used
#'   for `sd_C` and binary detection (typically the joint-mode dir).
#' @param frqfile_dir Character or NULL.
#' @param plink_name Character. Default `"ADSP_chr"`.
#' @param maf_cutoff Numeric, default `0.05`.
#' @param target_categories Character vector or NULL. Auto-detected from the
#'   first available run if NULL.
#' @param target_labels Character vector or NULL. Optional user-friendly display
#'   names for the target annotations, same length and order as the resolved
#'   `target_categories` (e.g. `c("quantile_eQTL", "eQTL")` to replace the
#'   polyfun `.results` names `c("ANNOT_1_0", "ANNOT_2_0")`). When given, every
#'   `target` column and `tau*`-block column name in the output is renamed;
#'   `params$target_categories` then holds the labels and
#'   `params$target_categories_orig` keeps the original polyfun names. When NULL
#'   (default), nothing is renamed - the original `.results` category names are
#'   used as before.
#'
#' @return List with `per_trait`, `meta` (three frames), `params`.
#'
#' @export
sldsc_postprocessing_pipeline <- function(trait_single_prefixes,
                                          trait_joint_prefix,
                                          target_anno_dir,
                                          frqfile_dir = NULL,
                                          plink_name = "ADSP_chr",
                                          maf_cutoff = 0.05,
                                          target_categories = NULL,
                                          target_labels = NULL) {
  trait_names <- names(trait_single_prefixes)
  if (is.null(trait_names))
    stop("sldsc_postprocessing_pipeline: trait_single_prefixes must be a named list.")

  message("[sldsc] Computing M_ref...")
  M_ref <- compute_sldsc_M_ref(target_anno_dir = target_anno_dir,
                               frqfile_dir = frqfile_dir,
                               plink_name = plink_name,
                               maf_cutoff = maf_cutoff)
  message(sprintf("[sldsc]   M_ref = %d (MAF cutoff %g)", M_ref, maf_cutoff))

  message("[sldsc] Computing per-annotation sd...")
  sd_annot_full <- compute_sldsc_annot_sd(target_anno_dir = target_anno_dir,
                                          frqfile_dir = frqfile_dir,
                                          plink_name = plink_name,
                                          maf_cutoff = maf_cutoff)
  message(sprintf("[sldsc]   sd computed for %d annotation columns",
                  length(sd_annot_full)))

  message("[sldsc] Detecting binary vs continuous annotations...")
  is_binary_full <- is_binary_sldsc_annot(target_anno_dir = target_anno_dir)

  # Polyfun renames target columns to `<col>_0` (file_idx=0 in --ref-ld-chr);
  # mirror that suffix so intersect() with pivot_run$categories matches.
  names(sd_annot_full)  <- paste0(names(sd_annot_full),  "_0")
  names(is_binary_full) <- paste0(names(is_binary_full), "_0")

  # Auto-detect target categories from a representative run.
  if (is.null(target_categories)) {
    pivot_run <- NULL
    if (!is.null(trait_joint_prefix) && length(trait_joint_prefix) > 0) {
      jp <- trait_joint_prefix[[1]]
      if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
        pivot_run <- tryCatch(read_sldsc_trait(jp), error = function(e) NULL)
      }
    }
    if (is.null(pivot_run) &&
        length(trait_single_prefixes) > 0L &&
        length(trait_single_prefixes[[1]]) > 0L) {
      pivot_run <- tryCatch(read_sldsc_trait(trait_single_prefixes[[1]][1]),
                            error = function(e) NULL)
    }
    if (is.null(pivot_run))
      stop("sldsc_postprocessing_pipeline: cannot auto-detect target_categories.")
    target_categories <- intersect(pivot_run$categories, names(sd_annot_full))
    # Fallback when .annot.gz names + "_0" do not match polyfun's .results
    # Category. Happens when the pipeline ran with --snp-list, since
    # ldsc.py --l2 hardcodes the LD score column to "L2" (single annot) or
    # "<annot>L2" (joint), instead of preserving the .annot.gz names.
    # Trust polyfun's invariant that target categories occupy the first
    # length(sd_annot_full) rows of .results, then rename positionally.
    if (length(target_categories) == 0L) {
      n_target   <- length(sd_annot_full)
      n_baseline <- length(pivot_run$categories) - n_target
      old_names  <- names(sd_annot_full)
      target_categories <- pivot_run$categories[seq_len(n_target)]
      names(sd_annot_full)  <- target_categories
      names(is_binary_full) <- target_categories
      baseline_preview <- if (n_baseline > 0L)
        paste(head(pivot_run$categories[-seq_len(n_target)], 3), collapse = ", ")
      else "(none)"
      message(sprintf(paste0(
        "[sldsc] sd_annot/is_binary names did not match polyfun .results categories;\n",
        "        falling back to positional rename (target = first %d rows of .results)\n",
        "        target  (%d): %s -> %s\n",
        "        baseline (%d): %s%s"),
        n_target,
        n_target, paste(old_names, collapse = ", "),
                  paste(target_categories, collapse = ", "),
        n_baseline, baseline_preview,
        if (n_baseline > 3L) ", ..." else ""))
    }
    message(sprintf("[sldsc] Auto-detected %d target categories", length(target_categories)))
  }

  baseline_categories <- character(0)
  if (!is.null(trait_joint_prefix) && length(trait_joint_prefix) > 0L) {
    jp <- trait_joint_prefix[[1]]
    if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
      pivot <- tryCatch(read_sldsc_trait(jp), error = function(e) NULL)
      if (!is.null(pivot))
        baseline_categories <- setdiff(pivot$categories, target_categories)
    }
  }
  if (length(baseline_categories) > 0L) {
    msg_head <- paste(head(baseline_categories, 5), collapse = ", ")
    msg_tail <- if (length(baseline_categories) > 5) ", ..." else ""
    message(sprintf("[sldsc] Detected %d baseline annotations: %s%s",
                    length(baseline_categories), msg_head, msg_tail))
  } else {
    message("[sldsc] No baseline annotations detected (joint-run prefix missing or unreadable).")
  }

  sd_annot <- sd_annot_full[target_categories]
  is_binary <- if (length(is_binary_full) > 0L) is_binary_full[target_categories] else
               setNames(rep(FALSE, length(target_categories)), target_categories)

  message(sprintf("[sldsc] Standardizing %d traits...", length(trait_names)))
  per_trait <- list()

  for (trait in trait_names) {
    # ---- single-mode ----
    single_summaries <- list()
    single_blocks    <- list()
    single_h2gs      <- numeric(0)
    sing_prefs <- trait_single_prefixes[[trait]]
    for (i in seq_along(target_categories)) {
      cat_name <- target_categories[i]
      if (i > length(sing_prefs)) break
      pref <- sing_prefs[i]
      run  <- tryCatch(read_sldsc_trait(pref), error = function(e) {
        warning(sprintf("[sldsc] Failed to read single %s for %s: %s",
                        cat_name, trait, e$message)); NULL
      })
      if (is.null(run)) next
      std <- tryCatch(
        standardize_sldsc_trait(run, sd_annot[cat_name], M_ref,
                                target_categories = cat_name, mode = "single"),
        error = function(e) {
          warning(sprintf("[sldsc] Failed to standardize single %s for %s: %s",
                          cat_name, trait, e$message)); NULL
        })
      if (is.null(std)) next
      single_summaries[[cat_name]] <- std$summary
      single_blocks[[cat_name]]    <- std$tau_star_blocks
      single_h2gs                  <- c(single_h2gs, std$h2g)
    }
    single_df <- if (length(single_summaries) > 0L)
                   do.call(rbind, single_summaries) else NULL
    if (!is.null(single_df)) rownames(single_df) <- NULL
    blocks_single <- if (length(single_blocks) > 0L) do.call(cbind, single_blocks) else NULL

    # ---- joint-mode ----
    joint_df       <- NULL
    blocks_joint   <- NULL
    joint_h2g      <- NA_real_
    n_blocks_trait <- NA_integer_
    if (!is.null(trait_joint_prefix) && trait %in% names(trait_joint_prefix)) {
      jp <- trait_joint_prefix[[trait]]
      if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
        run <- tryCatch(read_sldsc_trait(jp), error = function(e) {
          warning(sprintf("[sldsc] Failed to read joint for %s: %s",
                          trait, e$message)); NULL
        })
        if (!is.null(run)) {
          std <- tryCatch(
            standardize_sldsc_trait(run, sd_annot, M_ref,
                                    target_categories = target_categories,
                                    mode = "joint"),
            error = function(e) {
              warning(sprintf("[sldsc] Failed to standardize joint for %s: %s",
                              trait, e$message)); NULL
            })
          if (!is.null(std)) {
            joint_df       <- std$summary
            blocks_joint   <- std$tau_star_blocks
            joint_h2g      <- std$h2g
            n_blocks_trait <- std$n_blocks
          }
        }
      }
    }

    summary_wide <- .sldsc_assemble_trait_summary(single_df, joint_df,
                                                  target_categories, is_binary)
    per_trait[[trait]] <- list(
      summary                = summary_wide,
      tau_star_blocks_single = blocks_single,
      tau_star_blocks_joint  = blocks_joint,
      h2g                    = if (!is.na(joint_h2g)) joint_h2g
                               else if (length(single_h2gs) > 0L) median(single_h2gs)
                               else NA_real_,
      n_blocks               = n_blocks_trait
    )
  }

  message("[sldsc] Running random-effects meta across traits...")
  pt_view_single <- .sldsc_view_for_meta(per_trait, "single")
  pt_view_joint  <- .sldsc_view_for_meta(per_trait, "joint")

  build_table <- function(quantity, view, label) {
    rows <- list()
    for (cat in target_categories) {
      m <- meta_sldsc_random(view, cat, quantity)
      rows[[cat]] <- data.frame(
        target    = cat,
        is_binary = unname(is_binary[cat]),
        mean      = m$mean,
        se        = m$se,
        p         = m$p,
        n_traits  = m$n_traits,
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    nm_old <- c("mean", "se", "p")
    nm_new <- paste0(label, "_", nm_old)
    names(df)[names(df) %in% nm_old] <- nm_new
    df
  }

  meta_tau_star_single <- build_table("tau_star",   pt_view_single, "single")
  meta_tau_star_joint  <- build_table("tau_star",   pt_view_joint,  "joint")
  meta_E_single        <- build_table("enrichment", pt_view_single, "single")
  meta_ES_single       <- build_table("enrichstat", pt_view_single, "single")

  # Combine tau_star single + joint into one wide frame.
  meta_tau_star <- meta_tau_star_single
  ord <- match(meta_tau_star$target, meta_tau_star_joint$target)
  meta_tau_star$joint_mean <- meta_tau_star_joint$joint_mean[ord]
  meta_tau_star$joint_se   <- meta_tau_star_joint$joint_se[ord]
  meta_tau_star$joint_p    <- meta_tau_star_joint$joint_p[ord]

  # Two-channel enrichment meta: effect/SE from E meta, p from EnrichStat meta.
  meta_enrichment <- meta_E_single
  meta_enrichment$single_p <- meta_ES_single$single_p[match(meta_enrichment$target,
                                                             meta_ES_single$target)]

  # Pure EnrichStat meta (separate frame).
  meta_enrichstat <- meta_ES_single

  res <- list(
    per_trait = per_trait,
    meta = list(
      tau_star   = meta_tau_star,
      enrichment = meta_enrichment,
      enrichstat = meta_enrichstat
    ),
    params = list(
      maf_cutoff          = maf_cutoff,
      M_ref               = M_ref,
      target_categories   = target_categories,
      n_baseline          = length(baseline_categories),
      baseline_categories = baseline_categories,
      trait_names         = trait_names
    )
  )

  # Optional: relabel target categories to user-friendly display names.
  # If target_labels is NULL we keep the original polyfun .results names.
  if (!is.null(target_labels)) {
    target_labels <- as.character(target_labels)
    if (length(target_labels) != length(target_categories))
      stop(sprintf(paste0("sldsc_postprocessing_pipeline: target_labels has length %d but ",
                          "there are %d target categories (%s)."),
                   length(target_labels), length(target_categories),
                   paste(target_categories, collapse = ", ")))
    relab     <- setNames(target_labels, target_categories)
    relab_vec <- function(x) { y <- unname(relab[x]); y[is.na(y)] <- x[is.na(y)]; y }

    for (t in names(res$per_trait)) {
      pt <- res$per_trait[[t]]
      if (!is.null(pt$summary) && "target" %in% names(pt$summary))
        res$per_trait[[t]]$summary$target <- relab_vec(pt$summary$target)
      for (bn in c("tau_star_blocks_single", "tau_star_blocks_joint")) {
        b <- pt[[bn]]
        if (!is.null(b) && !is.null(colnames(b)))
          colnames(res$per_trait[[t]][[bn]]) <- relab_vec(colnames(b))
      }
    }
    for (mn in names(res$meta)) {
      if (!is.null(res$meta[[mn]]) && "target" %in% names(res$meta[[mn]]))
        res$meta[[mn]]$target <- relab_vec(res$meta[[mn]]$target)
    }
    res$params$target_categories_orig <- res$params$target_categories
    res$params$target_categories      <- unname(relab[target_categories])
    message(sprintf("[sldsc] Relabeled target categories: %s",
                    paste(sprintf("%s -> %s", target_categories,
                                  unname(relab[target_categories])), collapse = ", ")))
  }

  res
}
