# Stratified LD Score Regression (S-LDSC) post-processing wrappers around
# polyfun.
#
# This file provides the post-processing layer for the xqtl-protocol sLDSC
# pipeline: read polyfun outputs per trait, compute Gazal-style standardized
# tau* and the differential per-SNP heritability statistic (EnrichStat), and run
# DerSimonian-Laird random-effects meta-analysis across traits.
#
# Reference panel convention: all LD-derived quantities (baseline LD scores,
# target LD scores, regression weights, allele frequencies) must come from the
# same reference panel. Do not mix files from different panels (e.g. 1000G vs
# ADSP).
#
# MAF convention: by default we restrict to MAF > 5% per the sLDSC
# recommendation. Pass maf_cutoff = 0 to opt out (not recommended).
#
# Cross-type comparison: tau* (Gazal et al. 2017 standardization) is the
# cross-type comparable statistic. Use tau* to rank or meta-analyze annotations
# that mix binary and continuous types. E (proportion-based enrichment) is
# scale-dependent for continuous annotations and is only comparable within type.

# ---- internal helpers ----

.sldscStdCols <- c("CHR", "SNP", "BP", "CM", "A1", "A2", "MAF")

.sldscDetectAnnotCols <- function(filePath) {
    sample <- vroom(filePath, n_max = 5L, show_col_types = FALSE)
    setdiff(names(sample), .sldscStdCols)
}


#' @title Read S-LDSC outputs from polyfun for one trait/run
#'
#' @description Reads the regression outputs produced by `polyfun/ldsc.py` for a
#'   single polyfun run (one trait, one annotation set) and returns them as a
#'   tidy list ready for downstream standardization. Hides the underlying file
#'   formats; downstream code consumes only modeling quantities.
#'
#' @param prefix Character. Path prefix to the polyfun outputs for one
#'   trait/run. The function appends `.results`, `.log`, and `.part_delete` to
#'   this prefix. Example: `"/path/to/cwd/CAD_META.filtered.sumstats.gz"`.
#'
#' @return A named list. See `sldscPostprocessingPipeline` for components.
#'
#' @examples
#' prefix <- file.path(
#'   system.file("extdata", "sldsc_trait", package = "pecotmr"),
#'   "sumstats.parquet")
#' run <- readSldscTrait(prefix)
#' names(run)
#'
#' @importFrom stats setNames var na.omit
#' @importFrom utils head
#' @importFrom vroom vroom
#' @importFrom readr read_lines
#' @export
readSldscTrait <- function(prefix) {
    files <- str_c(prefix, c(".results", ".log", ".part_delete"))
    for (f in files) {
        if (!file.exists(f)) {
            msg <- glue("readSldscTrait: missing file: {f}")
            abort(msg)
        }
    }
    results <- vroom(files[1], show_col_types = FALSE)
    cats <- as.character(results$Category)
    h2g <- .readSldscH2g(files[2])
    deleteValues <- .readSldscBlocks(files[3], cats)
    list(
        categories = cats,
        tau = set_names(as.numeric(results$Coefficient), cats),
        tauSe = set_names(as.numeric(results[["Coefficient_std_error"]]), cats),
        enrichment = set_names(as.numeric(results$Enrichment), cats),
        enrichmentSe = set_names(
            as.numeric(results[["Enrichment_std_error"]]),
            cats
        ),
        enrichmentP = set_names(as.numeric(results[["Enrichment_p"]]), cats),
        propH2 = set_names(as.numeric(results[["Prop._h2"]]), cats),
        propSnps = set_names(as.numeric(results[["Prop._SNPs"]]), cats),
        h2g = h2g,
        tauBlocks = deleteValues,
        nBlocks = nrow(deleteValues)
    )
}

# Parse the total observed-scale h2g from an S-LDSC .log file.
# @noRd
.readSldscH2g <- function(logFile) {
    logLines <- read_lines(logFile)
    h2Line <- logLines[str_detect(logLines, "Total Observed scale h2:")]
    if (length(h2Line) == 0L) {
        msg <- glue(
            "readSldscTrait: could not find 'Total Observed scale h2:' in ",
            "{logFile}"
        )
        abort(msg)
    }
    h2g <- suppressWarnings(as.numeric(str_replace_all(
        h2Line[1],
        ".*h2: (-?[0-9.eE+-]+).*",
        "\\1"
    )))
    if (is.na(h2g)) {
        msg <- glue(
            "readSldscTrait: failed to parse h2g numeric from log line: ",
            "{h2Line[1]}"
        )
        abort(msg)
    }
    h2g
}

# Read the jackknife per-block matrix (.part_delete), validating its column
# count against the category set.
# @noRd
.readSldscBlocks <- function(deleteFile, cats) {
    deleteValues <- as.matrix(vroom(deleteFile, show_col_types = FALSE))
    if (ncol(deleteValues) != length(cats)) {
        msg <- glue(
            "readSldscTrait: .part_delete has {ncol(deleteValues)} columns ",
            "but .results has {length(cats)} categories."
        )
        abort(msg)
    }
    colnames(deleteValues) <- cats
    deleteValues
}


#' @title Read target annotation files (.annot.gz) into one table
#'
#' @description Reads the per-chromosome polyfun `.annot.gz` files in a
#'   directory and stacks them into a single \code{data.frame} of \code{CHR},
#'   \code{SNP}, and the annotation columns. This is the I/O step feeding the
#'   \code{annot} slot of \code{\link{SldscData}}; the computation
#'   (\code{\link{computeSldscAnnotSd}}, \code{\link{isBinarySldscAnnot}}) then
#'   runs on the loaded table, not on paths.
#'
#' @param targetAnnoDir Character. Directory of `.annot.gz` files.
#' @param annotCols Character or integer vector, default NULL. Annotation
#'   columns to keep. NULL keeps all non-standard columns (auto-detected).
#' @return A \code{data.frame}: \code{CHR}, \code{SNP}, and annotation columns.
#' @importFrom vroom vroom
#' @importFrom tidyselect all_of
#' @examples
#' sldsc <- system.file("extdata", "sldsc", package = "pecotmr")
#' readSldscAnnot(sldsc)
#' @export
readSldscAnnot <- function(targetAnnoDir, annotCols = NULL) {
    if (!dir.exists(targetAnnoDir)) {
        msg <- glue(
            "readSldscAnnot: targetAnnoDir does not exist: {targetAnnoDir}"
        )
        abort(msg)
    }
    annoFiles <- list.files(
        targetAnnoDir,
        pattern = "\\.annot\\.gz$",
        full.names = TRUE
    )
    if (length(annoFiles) == 0L) {
        msg <- glue("readSldscAnnot: no .annot.gz files in: {targetAnnoDir}")
        abort(msg)
    }

    detected <- .sldscDetectAnnotCols(annoFiles[1])
    colsUse <- if (is.null(annotCols)) {
        detected
    } else if (is.numeric(annotCols)) {
        detected[annotCols]
    } else {
        annotCols
    }
    if (length(colsUse) == 0L) {
        abort("readSldscAnnot: no annotation columns to read.")
    }

    parts <- map(annoFiles, .sldscReadAnnotFile, colsUse = colsUse)
    bind_rows(parts)
}


#' @title Read PLINK allele-frequency files (.frq) into one table
#'
#' @description Reads the per-chromosome PLINK `.frq` files for the reference
#'   panel and stacks them into a single \code{data.frame} of \code{CHR},
#'   \code{SNP}, \code{MAF}. Feeds the \code{frq} slot of
#'   \code{\link{SldscData}}.
#'
#' @param frqfileDir Character. Directory of `.frq` files.
#' @param plinkName Character. Filename prefix (files at
#'   `<plinkName><chr>.frq`). Falls back to all `*.frq` in the directory when
#'   the prefix matches nothing.
#' @return A \code{data.frame}: \code{CHR}, \code{SNP}, \code{MAF}.
#' @importFrom vroom vroom
#' @importFrom tidyselect all_of
#' @examples
#' sldsc <- system.file("extdata", "sldsc", package = "pecotmr")
#' readSldscFrq(sldsc, plinkName = "reference.")
#' @export
readSldscFrq <- function(frqfileDir, plinkName = "ADSP_chr") {
    if (!dir.exists(frqfileDir)) {
        msg <- glue("readSldscFrq: frqfileDir does not exist: {frqfileDir}")
        abort(msg)
    }
    pat <- str_c(
        "^",
        str_replace_all(plinkName, "([.])", "\\\\\\1"),
        "[0-9]+\\.frq$"
    )
    frqFiles <- list.files(frqfileDir, pattern = pat, full.names = TRUE)
    if (length(frqFiles) == 0L) {
        frqFiles <- list.files(
            frqfileDir,
            pattern = "\\.frq$",
            full.names = TRUE
        )
    }
    if (length(frqFiles) == 0L) {
        msg <- glue("readSldscFrq: no .frq files in: {frqfileDir}")
        abort(msg)
    }

    parts <- map(frqFiles, .sldscReadFrqFile)
    bind_rows(parts)
}


#' @title Compute per-annotation standard deviation, MAF-restricted
#'
#' @description Computes the standard deviation of each annotation column in the
#'   target annotation files, restricted to SNPs above a MAF cutoff via PLINK
#'   `.frq` files. Required for internal consistency with polyfun's regression,
#'   which operates on MAF > cutoff SNPs by default.
#'
#' @param sldscData An \code{\link{SldscData}} object (its \code{annot} and
#'   \code{frq} slots supply the annotation values and MAF, respectively).
#' @param mafCutoff Numeric, default `0.05`. Requires frq data when > 0.
#' @param annotCols Character or integer vector, default NULL. Annotation
#'   columns to compute sd for. If NULL, all annotation columns are used.
#'
#' @return Named numeric vector of \eqn{sd_C} values, one per annotation.
#'
#' @importFrom stats setNames var
#' @importFrom methods is
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' computeSldscAnnotSd(sldscData = sd)
#' @importFrom purrr map map_dbl compact reduce
#' @export
computeSldscAnnotSd <- function(sldscData, mafCutoff = 0.05, annotCols = NULL) {
    if (!is(sldscData, "SldscData")) {
        abort("computeSldscAnnotSd: `sldscData` must be an SldscData object.")
    }
    annot <- getAnnotData(sldscData)
    frq <- getFrqData(sldscData)
    if (mafCutoff > 0 && nrow(frq) == 0L) {
        msg <- glue(
            "computeSldscAnnotSd: mafCutoff = {mafCutoff} requires frq ",
            "data (read via readSldscFrq); none present."
        )
        abort(msg)
    }
    colsUse <- .sldscColsUse(sldscData, annotCols)
    # Pool within-chromosome variance (matches polyfun's per-file accumulation).
    contribs <- compact(map(
        unique(annot$CHR),
        .sldscChromVar,
        annot = annot,
        frq = frq,
        mafCutoff = mafCutoff,
        colsUse = colsUse
    ))
    den <- sum(map_dbl(contribs, "den"))
    if (den <= 0) {
        abort(
            "computeSldscAnnotSd: zero degrees of freedom after MAF filtering."
        )
    }
    num <- reduce(
        map(contribs, "num"),
        `+`,
        .init = set_names(numeric(length(colsUse)), colsUse)
    )
    sqrt(num / den)
}

# Resolve the annotation columns to process (all, by index, or by name).
# @noRd
.sldscColsUse <- function(sldscData, annotCols) {
    colsUse <- if (is.null(annotCols)) {
        getAnnotCols(sldscData)
    } else if (is.numeric(annotCols)) {
        getAnnotCols(sldscData)[annotCols]
    } else {
        annotCols
    }
    if (length(colsUse) == 0L) {
        abort("computeSldscAnnotSd: no annotation columns to process.")
    }
    colsUse
}

# Within-chromosome (n-1)-weighted variance contribution per annotation column.
# Returns list(num = named weighted-variance vector, den = n-1), or NULL when
# the chromosome has <= 1 usable variant after MAF filtering.
# @noRd
.sldscChromVar <- function(chrom, annot, frq, mafCutoff, colsUse) {
    dat <- filter(annot, .data$CHR == chrom)
    if (mafCutoff > 0) {
        dat <- inner_join(
            dat,
            select(frq, all_of(c("SNP", "MAF"))),
            by = "SNP"
        )
        dat <- filter(dat, !is.na(.data$MAF) & .data$MAF > mafCutoff)
    }
    if (nrow(dat) <= 1L) {
        return(NULL)
    }
    nMinus1 <- nrow(dat) - 1L
    num <- map_dbl(colsUse, .sldscColVarContrib, dat = dat, nMinus1 = nMinus1)
    list(num = set_names(num, colsUse), den = nMinus1)
}


#' @title Reference-panel SNP count (the M_ref used to standardise tau*)
#'
#' @description `M_ref` is the number of SNPs in the REFERENCE PANEL over which
#'   heritability is partitioned in the sLDSC model
#'   (`h2(C) = sum_(j in M_ref) a_C(j) sum_(C') tau_(C') a_(C')(j)`). It is
#'   panel-defined and is **not** the regression SNP set (HapMap3 ~1M) nor any
#'   HM3-subsetted target output:
#'   \itemize{
#'     \item `mafCutoff > 0` (Gazal/Finucane convention): count MAF > cutoff
#'       SNPs across all `.frq` files (the same set polyfun's `.l2.M_5_50`
#'       sums).
#'     \item `mafCutoff == 0` (all-M variant): count ALL SNPs across all
#'       `.frq` files (the same set polyfun's `.l2.M` sums).
#'   }
#'   When no frq data is present, `mafCutoff == 0` falls back to the number of
#'   annotation rows; `mafCutoff > 0` errors (a MAF-restricted count needs frq).
#'
#' @param sldscData An \code{\link{SldscData}} object (its \code{frq} slot is
#'   the reference-panel SNP set).
#' @param mafCutoff Numeric, default `0.05`.
#'
#' @return Scalar integer.
#'
#' @importFrom methods is
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' computeSldscMRef(sldscData = sd)
#' @export
computeSldscMRef <- function(sldscData, mafCutoff = 0.05) {
    if (!is(sldscData, "SldscData")) {
        abort("computeSldscMRef: `sldscData` must be an SldscData object.")
    }
    frq <- getFrqData(sldscData)
    if (nrow(frq) > 0L) {
        return(as.integer(
            if (mafCutoff > 0) {
                sum(!is.na(frq$MAF) & frq$MAF > mafCutoff)
            } else {
                nrow(frq)
            }
        ))
    }
    if (mafCutoff > 0) {
        msg <- glue(
            "computeSldscMRef: mafCutoff = {mafCutoff} requires frq ",
            "data (read via readSldscFrq); none present."
        )
        abort(msg)
    }
    as.integer(nrow(getAnnotData(sldscData)))
}


#' @title Detect whether each annotation is binary or continuous
#'
#' @description Inspects each annotation column and returns whether its values
#'   lie in \{0, 1\} (binary) or take other values (continuous).
#'
#' @param sldscData An \code{\link{SldscData}} object.
#' @param annotCols Character or integer vector, default NULL.
#'
#' @return Named logical vector: TRUE for binary, FALSE for continuous.
#'
#' @importFrom stats setNames na.omit
#' @importFrom methods is
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' isBinarySldscAnnot(sd)
#' @export
isBinarySldscAnnot <- function(sldscData, annotCols = NULL) {
    if (!is(sldscData, "SldscData")) {
        abort("isBinarySldscAnnot: `sldscData` must be an SldscData object.")
    }
    annot <- getAnnotData(sldscData)
    colsUse <- if (is.null(annotCols)) {
        getAnnotCols(sldscData)
    } else if (is.numeric(annotCols)) {
        getAnnotCols(sldscData)[annotCols]
    } else {
        annotCols
    }

    isBinary <- set_names(rep(TRUE, length(colsUse)), colsUse)
    for (col in colsUse) {
        vals <- unique(na.omit(as.numeric(annot[[col]])))
        if (any(!is_in(vals, c(0, 1)))) isBinary[[col]] <- FALSE
    }
    isBinary
}


#' @title Standardize tau and compute EnrichStat for one polyfun run
#'
#' @description Applies the Gazal standardization \eqn{\tau^*_C = \tau_C \cdot
#'   sd_C \cdot M_{ref} / h^2_g} to the point and to each jackknife block. For
#'   `mode = "single"`, additionally computes EnrichStat and back-solves its
#'   standard error from polyfun's reported `Enrichment_p` using \eqn{|Z| =
#'   \Phi^{-1}(1 - p/2)}.
#'
#' @param sldscData An \code{\link{SldscData}} object (the run is pulled from it
#'   via \code{getTraitRun}).
#' @param trait Character. Trait name (a key of the SldscData traits list).
#' @param mode Character: `"single"` or `"joint"`.
#' @param idx Integer or NULL. For `mode = "single"`, which of the trait's
#'   single-target runs to standardize.
#' @param sdAnnot Named numeric vector from \code{\link{computeSldscAnnotSd}}.
#' @param MRef Scalar from \code{\link{computeSldscMRef}}.
#' @param targetCategories Character vector or NULL. If NULL, intersects the
#'   run's `categories` with `names(sdAnnot)`.
#'
#' @return A list with `summary` (data frame), `tau_star_blocks` (matrix),
#'   `h2g`, `nBlocks`, `mode`.
#'
#' @importFrom stats qnorm var
#' @importFrom methods is
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' sd <- SldscData(annot = annot, frq = frq,
#'   traits = setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY")))
#' sdAnnot <- computeSldscAnnotSd(sd)
#' MRef <- computeSldscMRef(sd)
#' standardizeSldscTrait(sd, "traitX", mode = "single", idx = 1,
#'   sdAnnot = sdAnnot, MRef = MRef, targetCategories = "annot_A_0")
#' @export
standardizeSldscTrait <- function(
    sldscData,
    trait,
    mode = c("single", "joint"),
    idx = NULL,
    sdAnnot,
    MRef,
    targetCategories = NULL
) {
    if (!is(sldscData, "SldscData")) {
        abort("standardizeSldscTrait: `sldscData` must be an SldscData object.")
    }
    mode <- arg_match(mode)
    traitData <- .stdTraitRun(sldscData, trait, mode, idx)
    targetCategories <- .stdTargetCategories(
        traitData,
        sdAnnot,
        targetCategories
    )
    targetIdx <- .stdTargetIdx(traitData, targetCategories)
    h2g <- traitData$h2g
    sdTarget <- .stdSdTarget(sdAnnot, targetCategories)
    tau <- as.numeric(traitData$tau[targetCategories])
    tauSe <- as.numeric(traitData$tauSe[targetCategories])
    blocksTarget <- traitData$tauBlocks[, targetIdx, drop = FALSE]
    ts <- standardizeTauStar(tau, blocksTarget, sdTarget, MRef, h2g)
    summaryDf <- .stdSummaryDf(targetCategories, tau, tauSe, ts)
    if (mode == "single") {
        summaryDf <- .stdEnrichmentCols(
            summaryDf,
            traitData,
            targetCategories,
            h2g,
            MRef
        )
    }
    tauStarBlocks <- sweep(blocksTarget, 2L, sdTarget * MRef / h2g, FUN = "*")
    list(
        summary = summaryDf,
        tau_star_blocks = tauStarBlocks,
        h2g = h2g,
        nBlocks = nrow(blocksTarget),
        mode = mode
    )
}

# Base per-target summary frame (tau + tau* columns).
# @noRd
.stdSummaryDf <- function(targetCategories, tau, tauSe, ts) {
    tibble(
        target = targetCategories,
        tau = unname(tau),
        tauSe = unname(tauSe),
        tauStar = unname(ts$tauStar),
        tauStarSe = unname(ts$tauStarSe)
    )
}

# Fetch the requested trait's single/joint run (error when absent).
# @noRd
.stdTraitRun <- function(sldscData, trait, mode, idx) {
    traitData <- getTraitRun(sldscData, trait, mode, idx)
    if (is.null(traitData)) {
        idxNote <- if (!is.null(idx)) glue(" (idx={idx})") else ""
        msg <- glue(
            "standardizeSldscTrait: no {mode} run for trait ",
            "'{trait}'{idxNote}."
        )
        abort(msg)
    }
    traitData
}

# Resolve the target categories (intersect run categories with sdAnnot when
# unspecified) and require a non-empty set.
# @noRd
.stdTargetCategories <- function(traitData, sdAnnot, targetCategories) {
    if (is.null(targetCategories)) {
        targetCategories <- intersect(traitData$categories, names(sdAnnot))
    }
    if (length(targetCategories) == 0L) {
        abort("standardizeSldscTrait: no target categories.")
    }
    targetCategories
}

# Match target categories to their run positions (error on any missing).
# @noRd
.stdTargetIdx <- function(traitData, targetCategories) {
    targetIdx <- match(targetCategories, traitData$categories)
    if (any(is.na(targetIdx))) {
        msg <- glue(
            "standardizeSldscTrait: missing categories: ",
            "{str_flatten(targetCategories[is.na(targetIdx)], ', ')}"
        )
        abort(msg)
    }
    targetIdx
}

# Per-target annotation SDs (warn on zero/NA, which yield NA/0 tau*).
# @noRd
.stdSdTarget <- function(sdAnnot, targetCategories) {
    sdTarget <- as.numeric(sdAnnot[targetCategories])
    if (any(is.na(sdTarget) | sdTarget == 0)) {
        msg <- glue(
            "standardizeSldscTrait: zero/NA sd for some targets; tau* ",
            "will be NA/0."
        )
        warn(msg)
    }
    sdTarget
}

# Append the single-mode enrichment + enrichment-statistic columns to the
# summary frame (enrichstat SE back-solved from the polyfun enrichment p-value).
# @noRd
.stdEnrichmentCols <- function(
    summaryDf,
    traitData,
    targetCategories,
    h2g,
    MRef
) {
    pH2 <- as.numeric(traitData$propH2[targetCategories])
    pM <- as.numeric(traitData$propSnps[targetCategories])
    enrichstat <- (h2g / MRef) * ((pH2 / pM) - (1 - pH2) / (1 - pM))
    enrichP <- as.numeric(traitData$enrichmentP[targetCategories])
    absZ <- qnorm(1 - enrichP / 2)
    enrichstatSe <- abs(enrichstat) / absZ
    enrichstatSe[!is.finite(absZ) | absZ <= 0] <- NA_real_
    summaryDf$enrichment <- as.numeric(traitData$enrichment[targetCategories])
    summaryDf$enrichmentSe <- as.numeric(
        traitData$enrichmentSe[targetCategories]
    )
    summaryDf$enrichmentP <- enrichP
    summaryDf$enrichstat <- enrichstat
    summaryDf$enrichstatSe <- enrichstatSe
    summaryDf
}


#' @title Random-effects meta-analysis of S-LDSC quantities across traits
#'
#' @description Random-effects meta-analysis (DerSimonian-Laird, via
#'   \code{metafor::rma}) of one S-LDSC quantity for one annotation across
#'   multiple traits.
#'
#' @details Per-trait \eqn{SE_i} sources: - `quantity = "tauStar"`: jackknife SE
#'   from per-block \eqn{\tau^*}. - `quantity = "enrichment"`: polyfun-reported
#'   `Enrichment_std_error`. - `quantity = "enrichstat"`: back-solved SE from
#'   polyfun's `Enrichment_p`.
#'
#' @param perTraitEstimates Named list of per-trait results (each with a
#'   `summary` data frame).
#' @param category Character. Annotation name to meta-analyze.
#' @param quantity Character: `"tauStar"`, `"enrichment"`, or `"enrichstat"`.
#'
#' @return List with `mean`, `se`, `p`, `nTraits`, `traitsUsed`, `tau2`.
#'
#' @importFrom stats pnorm
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' pp <- sldscPostprocessingPipeline(sd)
#' metaSldscRandom(pp$per_trait, category = "annot_A_0",
#'   quantity = "enrichment")
#' @export
metaSldscRandom <- function(
    perTraitEstimates,
    category,
    quantity = c("tauStar", "enrichment", "enrichstat")
) {
    quantity <- arg_match(quantity)
    cols <- .metaColPair(quantity)
    traitNames <- names(perTraitEstimates) %||%
        as.character(seq_along(perTraitEstimates))
    collected <- compact(map(
        seq_along(perTraitEstimates),
        .metaTraitContrib,
        perTraitEstimates = perTraitEstimates,
        category = category,
        cols = cols,
        traitNames = traitNames
    ))
    means <- map_dbl(collected, "mean")
    ses <- map_dbl(collected, "se")
    used <- map_chr(collected, "trait")
    if (length(means) < 2L) {
        return(.metaEmptyResult(length(means), used))
    }
    meta <- .rmaMeta(means, ses)
    list(
        mean = meta$mean,
        se = meta$se,
        p = as.numeric(.zToPvalue(meta$mean / meta$se)),
        nTraits = length(means),
        traitsUsed = used,
        tau2 = meta$tau2
    )
}

# (value, SE) column-name pair in the per-trait summary for a meta quantity.
# @noRd
.metaColPair <- function(quantity) {
    list(
        tauStar = c("tauStar", "tauStarSe"),
        enrichment = c("enrichment", "enrichmentSe"),
        enrichstat = c("enrichstat", "enrichstatSe")
    )[[quantity]]
}

# Extract one trait's (mean, se, trait) for `category`, or NULL when the trait
# lacks a usable finite positive-SE estimate.
# @noRd
.metaExtractTrait <- function(pt, category, cols, traitName) {
    if (is.null(pt) || is.null(pt$summary)) {
        return(NULL)
    }
    row <- pt$summary[pt$summary$target == category, , drop = FALSE]
    if (nrow(row) == 0L || !all(is_in(cols, names(row)))) {
        return(NULL)
    }
    m <- as.numeric(row[[cols[1]]])[1]
    s <- as.numeric(row[[cols[2]]])[1]
    if (is.na(m) || is.na(s) || !is.finite(s) || s <= 0) {
        return(NULL)
    }
    list(mean = m, se = s, trait = traitName)
}

# The all-NA meta result used when fewer than two traits contribute.
# @noRd
.metaEmptyResult <- function(nTraits, used) {
    list(
        mean = NA_real_,
        se = NA_real_,
        p = NA_real_,
        nTraits = nTraits,
        traitsUsed = used,
        tau2 = NA_real_
    )
}


# Append the single/joint enrichment columns (suffix-capitalized) to `out`,
# aligned to out$target; missing sources fill NA.
# @noRd
.sldscAddCols <- function(out, src, suffix) {
    colsToAdd <- c(
        "tau",
        "tauSe",
        "tauStar",
        "tauStarSe",
        "enrichment",
        "enrichmentSe",
        "enrichmentP",
        "enrichstat",
        "enrichstatSe"
    )
    suffixCap <- str_c(str_to_upper(str_sub(suffix, 1, 1)), str_sub(suffix, 2))
    for (c in colsToAdd) {
        newcol <- str_c(c, suffixCap)
        if (!is.null(src) && is_in(c, names(src))) {
            out[[newcol]] <- src[[c]][match(out$target, src$target)]
        } else {
            out[[newcol]] <- NA_real_
        }
    }
    out
}

# Internal helper: assemble a wide per-trait summary frame with single + joint
# columns side by side.
.sldscAssembleTraitSummary <- function(
    singleDf,
    jointDf,
    targetCategories,
    isBinaryVec
) {
    rows <- if (!is.null(singleDf)) {
        singleDf$target
    } else if (!is.null(jointDf)) {
        jointDf$target
    } else {
        targetCategories
    }
    out <- tibble(
        target = rows,
        isBinary = unname(isBinaryVec[rows])
    )

    out <- .sldscAddCols(out, singleDf, "single")
    out <- .sldscAddCols(out, jointDf, "joint")
    out
}


# Internal helper: build a per-trait list view that metaSldscRandom can read.
# Each list element has a $summary frame with the requested mode's columns
# renamed to the canonical names (tauStar, tauStarSe, enrichment, ...).
.sldscViewForMeta <- function(perTrait, suffix) {
    map(perTrait, .sldscTraitMetaView, suffix = suffix)
}

# Meta-analyze `quantity` across all target categories for one view, returning a
# per-category named list of metaSldscRandom results.
# @noRd
.sldscPerCategory <- function(view, quantity, targetCategories) {
    set_names(
        map(
            targetCategories,
            .sldscMetaForCategory,
            view = view,
            quantity = quantity
        ),
        targetCategories
    )
}

#' Random-effects meta-analysis over a subset of sLDSC traits
#'
#' Re-run the random-effects meta-analysis (DerSimonian-Laird, via
#' \code{metafor::rma}) on a chosen subset of the per-trait standardised tables
#' produced by \code{\link{sldscPostprocessingPipeline}} -- no regression is
#' re-run, only the already-standardised per-trait estimates are re-meta'd.
#' Powers a "meta on a subset of traits" workflow: pick traits, pick target
#' annotation categories, and get the per-category tau* / enrichment /
#' enrichstat meta results back.
#'
#' @param postprocessResult The list returned by
#'   \code{\link{sldscPostprocessingPipeline}}. Must carry a \code{$per_trait}
#'   element; \code{$params$target_categories} is used when
#'   \code{targetCategories} is \code{NULL}.
#' @param subsetTraits Character vector of trait ids to meta over; every id must
#'   be present in \code{postprocessResult$per_trait}.
#' @param targetCategories Optional character vector of target annotation names.
#'   Defaults to \code{postprocessResult$params$target_categories}.
#' @return A list with \code{tau_star_single}, \code{tau_star_joint},
#'   \code{enrichment}, and \code{enrichstat}; each is a per-category named list
#'   of \code{\link{metaSldscRandom}} results.
#' @seealso \code{\link{sldscPostprocessingPipeline}},
#'   \code{\link{metaSldscRandom}}
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' pp <- sldscPostprocessingPipeline(sd)
#' sldscSubsetMeta(pp, subsetTraits = "traitX")
#' @export
sldscSubsetMeta <- function(
    postprocessResult,
    subsetTraits,
    targetCategories = NULL
) {
    perTrait <- postprocessResult$per_trait
    if (is.null(perTrait)) {
        abort(
            "sldscSubsetMeta: `postprocessResult` has no `per_trait` element."
        )
    }
    targetCategories <- .sldscSubsetTargets(postprocessResult, targetCategories)
    missingTraits <- setdiff(subsetTraits, names(perTrait))
    if (length(missingTraits) > 0L) {
        msg <- glue(
            "sldscSubsetMeta: trait(s) absent from `per_trait`: ",
            "{str_flatten(missingTraits, ', ')}"
        )
        abort(msg)
    }
    sub <- perTrait[subsetTraits]
    viewSingle <- .sldscViewForMeta(sub, "single")
    viewJoint <- .sldscViewForMeta(sub, "joint")
    list(
        tau_star_single = .sldscPerCategory(
            viewSingle,
            "tauStar",
            targetCategories
        ),
        tau_star_joint = .sldscPerCategory(
            viewJoint,
            "tauStar",
            targetCategories
        ),
        enrichment = .sldscPerCategory(
            viewSingle,
            "enrichment",
            targetCategories
        ),
        enrichstat = .sldscPerCategory(
            viewSingle,
            "enrichstat",
            targetCategories
        )
    )
}

# Resolve the target categories (explicit, else from the pipeline params).
# @noRd
.sldscSubsetTargets <- function(postprocessResult, targetCategories) {
    if (is.null(targetCategories)) {
        targetCategories <- postprocessResult$params$target_categories
    }
    if (is.null(targetCategories) || length(targetCategories) == 0L) {
        msg <- glue(
            "sldscSubsetMeta: no `targetCategories` supplied and none found ",
            "in `postprocessResult$params$target_categories`."
        )
        abort(msg)
    }
    targetCategories
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Read one annotation file's CHR/SNP + requested annotation columns.
# @noRd
.sldscReadAnnotFile <- function(f, colsUse) {
    vroom(
        f,
        col_select = all_of(c("CHR", "SNP", colsUse)),
        show_col_types = FALSE
    )
}

# Read one .frq file's CHR/SNP/MAF columns.
# @noRd
.sldscReadFrqFile <- function(f) {
    vroom(
        f,
        col_select = all_of(c("CHR", "SNP", "MAF")),
        show_col_types = FALSE
    )
}

# (n-1)*Var contribution of annotation column `col` (0 for a constant column).
# @noRd
.sldscColVarContrib <- function(col, dat, nMinus1) {
    v <- var(as.numeric(dat[[col]]), na.rm = TRUE)
    if (is.na(v)) 0 else nMinus1 * v
}

# The (mean, se, trait) contribution of per-trait estimate `i`, or NULL.
# @noRd
.metaTraitContrib <- function(
    i,
    perTraitEstimates,
    category,
    cols,
    traitNames
) {
    .metaExtractTrait(perTraitEstimates[[i]], category, cols, traitNames[i])
}

# One trait's meta view: rename the `suffix`-mode columns to canonical names.
# @noRd
.sldscTraitMetaView <- function(pt, suffix) {
    if (is.null(pt$summary)) {
        return(NULL)
    }
    df <- pt$summary
    colsHave <- c(
        "tauStar",
        "tauStarSe",
        "enrichment",
        "enrichmentSe",
        "enrichmentP",
        "enrichstat",
        "enrichstatSe"
    )
    suffixCap <- str_c(
        str_to_upper(str_sub(suffix, 1, 1)),
        str_sub(suffix, 2)
    )
    srcCols <- str_c(colsHave, suffixCap)
    avail <- is_in(srcCols, names(df))
    if (!any(avail)) {
        return(NULL)
    }
    newDf <- tibble(target = df$target)
    for (k in seq_along(colsHave)) {
        if (avail[k]) newDf[[colsHave[k]]] <- df[[srcCols[k]]]
    }
    list(summary = newDf)
}

# The random-effects meta result for one target category of a view.
# @noRd
.sldscMetaForCategory <- function(category, view, quantity) {
    metaSldscRandom(view, category, quantity)
}
