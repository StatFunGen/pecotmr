#' @title Annotation Handling for Stratified Heritability
#' @description Read and manage genomic annotations for stratified heritability
#'   analysis. Supports BED, BigWig, and LDSC .annot formats.
#' @name pecotmr-h2-annotations
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom stringr str_detect str_to_lower
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @include AllGenerics.R
NULL

# =============================================================================
# Reader method
# =============================================================================

#' @rdname readAnnotations
#' @export
setMethod(
    "readAnnotations",
    signature(path = "character"),
    function(path, snpRanges, annotationMeta = NULL, genome = "hg19") {
        if (is.null(names(path))) {
            msg <- glue(
                "'path' must be a named character vector (names = ",
                "annotation names)"
            )
            abort(msg)
        }
        annotNames <- names(path)
        nAnnots <- length(path)
        annotMat <- .readAnnotMatrix(
            path,
            snpRanges,
            annotNames,
            length(snpRanges),
            nAnnots
        )
        if (is.null(annotationMeta)) {
            annotationMeta <- tibble(
                name = annotNames,
                tier = rep("candidate", nAnnots),
                type = .readAnnotTypes(path)
            )
        }
        AnnotationMatrix(annotMat, snpRanges, annotationMeta, genome)
    }
)

# Auto-detected annotation types (continuous for BigWig, else binary). Unnamed
# so the tibble column carries no element names (tibble, unlike data.frame,
# would otherwise keep the paths' names on the column).
# @noRd
.readAnnotTypes <- function(paths) {
    unname(map_chr(paths, .annotType))
}

# @noRd
.annotType <- function(p) {
    if (.annotDetectFormat(p) == "bigwig") "continuous" else "binary"
}

# @noRd
.readAnnotColumnAt <- function(i, paths, snpRanges, annotNames) {
    .readAnnotColumn(paths[i], snpRanges, annotNames[i])
}

# Build the (SNP x annotation) matrix by reading each annotation column.
# @noRd
.readAnnotMatrix <- function(paths, snpRanges, annotNames, nSnps, nAnnots) {
    matrix(
        unname(list_c(map(
            seq_along(paths),
            .readAnnotColumnAt,
            paths = paths,
            snpRanges = snpRanges,
            annotNames = annotNames
        ))),
        nrow = nSnps,
        ncol = nAnnots,
        dimnames = list(NULL, annotNames)
    )
}

# One annotation column, dispatched by detected format (BigWig / .annot / BED).
# @noRd
.readAnnotColumn <- function(path, snpRanges, annotName) {
    fmt <- .annotDetectFormat(path)
    if (fmt == "bigwig") {
        .readBigwigAtSnps(path, snpRanges)
    } else if (fmt == "ldsc_annot") {
        .readLdscAnnot(path, snpRanges, annotName)
    } else {
        .readBedAnnotation(path, snpRanges)
    }
}

# =============================================================================
# Internal helpers
# =============================================================================

#' @title Detect Annotation File Format
#' @description Detect annotation file format from extension. This is separate
#'   from the genotype format detection (\code{.genotypePathFormat}) because
#'   BED annotation files (genomic intervals for rtracklayer) must be
#'   distinguished from plink BED files.
#' @param path Character, file path.
#' @return Character, one of "bigwig", "ldsc_annot", or "bed".
#' @keywords internal
.annotDetectFormat <- function(path) {
    lpath <- str_to_lower(path)
    if (str_detect(lpath, "\\.annot\\.gz$")) {
        return("ldsc_annot")
    }

    ext <- str_to_lower(file_ext(path))
    switch(
        ext,
        "bw" = ,
        "bigwig" = "bigwig",
        "annot" = "ldsc_annot",
        # Default: treat as BED (genomic interval file for rtracklayer)
        "bed"
    )
}

#' @title Read BigWig Scores at SNP Positions
#' @description Import scores from a BigWig file at specified SNP positions.
#' @param bwPath Character, path to a BigWig file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of scores (length = number of SNPs).
#' @keywords internal
.readBigwigAtSnps <- function(bwPath, snpRanges) {
    bw <- rtracklayer::BigWigFile(bwPath)
    scores <- rtracklayer::import(bw, which = snpRanges, as = "NumericList")
    # Take mean score at each SNP position
    map_dbl(scores, .bwMeanScore)
}

#' @title Read BED Annotation
#' @description Read a BED file and compute binary overlap with SNP positions.
#' @param bedPath Character, path to a BED file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @return Numeric vector of 0/1 values (length = number of SNPs).
#' @keywords internal
.readBedAnnotation <- function(bedPath, snpRanges) {
    regions <- rtracklayer::import(bedPath)
    hits <- findOverlaps(snpRanges, regions)
    result <- replace(rep(0L, length(snpRanges)), queryHits(hits), 1L)
    as.numeric(result)
}

#' @title Read LDSC Annotation File
#' @description Read an S-LDSC .annot[.gz] file and extract a named annotation
#'   column, matched to SNP positions.
#' @param annotPath Character, path to an .annot or .annot.gz file.
#' @param snpRanges A \code{GRanges} object with SNP positions.
#' @param annotName Character, name of the annotation column to extract.
#' @return Numeric vector of annotation values (length = number of SNPs).
#' @keywords internal
.readLdscAnnot <- function(annotPath, snpRanges, annotName) {
    # S-LDSC .annot files are tab-separated with columns: CHR, BP, SNP, CM, ...
    dt <- vroom(annotPath, show_col_types = FALSE)

    if (!is_in(annotName, colnames(dt))) {
        msg <- glue("Annotation column '{annotName}' not found in {annotPath}")
        abort(msg)
    }

    if (!all(is_in(c("CHR", "BP"), colnames(dt)))) {
        abort("LDSC annot file must contain CHR and BP columns")
    }

    # Build GRanges from the annot file positions
    annotGr <- GRanges(
        seqnames = withChrPrefix(dt$CHR),
        ranges = IRanges(start = dt$BP, width = 1L)
    )

    # Match SNPs by genomic position
    hits <- findOverlaps(snpRanges, annotGr)

    # Unmatched SNPs keep the default 0.
    replace(
        rep(0, length(snpRanges)),
        queryHits(hits),
        as.numeric(dt[[annotName]][subjectHits(hits)])
    )
}

# (The AnnotationMatrix() constructor and the baselineAnnotations /
# testedAnnotations tier
# accessors now live in R/AnnotationMatrix.R alongside the class definition.)

# The mean BigWig score at one SNP (0 when the SNP has no overlapping
# intervals).
# @noRd
.bwMeanScore <- function(x) {
    if (length(x) > 0) mean(x) else 0
}

# =============================================================================
# Annotations from fine-mapping results
# =============================================================================

#' @title Annotations From Fine-Mapping Results
#' @description Turn fine-mapping results -- typically molecular QTLs -- into
#'   per-variant annotations over an LD reference's variants, for
#'   \code{\link{estimateH2}} to test for heritability enrichment. Each
#'   (\code{study}, \code{context}, \code{method}) unit of \code{x} pools its
#'   traits (e.g. genes) into one column per requested type:
#'   \describe{
#'     \item{\code{"maxCpp"}}{A variant's largest PIP over the unit's traits
#'       (continuous; Hormozdiari et al. 2018).}
#'     \item{\code{"sumPip"}}{The sum of its PIPs over the traits: the
#'       expected number of the unit's causal effects it carries
#'       (continuous).}
#'     \item{\code{"credibleSet"}}{1 when it sits in a credible set of any
#'       trait at \code{coverage} (and \code{minPurity}), else 0.}
#'     \item{\code{"highPip"}}{1 when its largest PIP is at least
#'       \code{pipThreshold}, else 0.}
#'   }
#'
#'   Variants are matched to the reference allele-aware (id format, ref/alt
#'   swap, strand) by their coordinates and alleles, so a reference naming
#'   its variants by rsID matches too; a fine-mapped variant the reference
#'   lacks contributes nothing, and a reference variant no trait fine-mapped
#'   takes 0. A column that comes out all zero -- a unit none of whose
#'   variants are in the reference, or none of whose PIPs reach
#'   \code{pipThreshold} -- is dropped, with a warning.
#'
#'   The rows follow the reference's variants in its order, as
#'   \code{estimateH2} requires.
#' @param x A \code{\link{QtlFineMappingResult}} or
#'   \code{\link{GwasFineMappingResult}}.
#' @param ldRef The \code{\linkS4class{LdStatistic}} the annotations will be
#'   used with.
#' @param types Which annotation types to build; see above. Default all four.
#' @param coverage Credible-set coverage for \code{"credibleSet"}. Default
#'   \code{0.95}.
#' @param minPurity Minimum credible-set purity for \code{"credibleSet"};
#'   \code{NULL} (default) applies none.
#' @param pipThreshold PIP cutoff for \code{"highPip"}. Default \code{0.5}.
#' @param tier The tier the columns take: \code{"candidate"} (default), so
#'   \code{estimateH2} score-tests each for enrichment over the baseline, or
#'   \code{"baseline"}, to fit them jointly with it.
#' @return An \code{\linkS4class{AnnotationMatrix}} with one row per variant of
#'   \code{ldRef}. Columns are named
#'   \code{"<study>|<context>|<method>|<type>"}, omitting an axis \code{x}
#'   does not have.
#' @examples
#' data(qtlFineMappingExample)
#' fm <- topLoci(qtlFineMappingExample, signalCutoff = 0)
#' ref <- LdScore(
#'     ldBlocks = GenomicRanges::GRanges(
#'         "chr22",
#'         IRanges::IRanges(min(GenomicRanges::start(fm)),
#'             max(GenomicRanges::start(fm)))
#'     ),
#'     snpInfo = data.frame(
#'         SNP = S4Vectors::mcols(fm)$variant_id,
#'         CHR = "chr22",
#'         BP = GenomicRanges::start(fm),
#'         A1 = S4Vectors::mcols(fm)$A1,
#'         A2 = S4Vectors::mcols(fm)$A2
#'     ),
#'     nSamples = 1000L,
#'     inSample = FALSE,
#'     ldScores = matrix(1, length(fm), 1, dimnames = list(NULL, "base_l2")),
#'     ldScoreWeights = rep(1, length(fm))
#' )
#' annot <- fineMappingAnnotations(qtlFineMappingExample, ref, pipThreshold = 0.1)
#' colnames(annot)
#' @importFrom checkmate assertMultiClass assertClass assertNumber
#' @importFrom rlang arg_match arg_match0
#' @importFrom stringr str_extract str_flatten
#' @export
fineMappingAnnotations <- function(
    x,
    ldRef,
    types = c("maxCpp", "sumPip", "credibleSet", "highPip"),
    coverage = 0.95,
    minPurity = NULL,
    pipThreshold = 0.5,
    tier = "candidate"
) {
    assertMultiClass(x, c("QtlFineMappingResult", "GwasFineMappingResult"))
    assertClass(ldRef, "LdStatistic")
    types <- arg_match(types, multiple = TRUE)
    assertNumber(coverage, lower = 0, upper = 1)
    assertNumber(minPurity, lower = 0, upper = 1, null.ok = TRUE)
    assertNumber(pipThreshold, lower = 0, upper = 1)
    tier <- arg_match0(tier, c("candidate", "baseline"))
    # The reference's variants, read off its ranges and alleles.
    refTable <- .variantTable(ldRef)
    tuples <- .enrAnnotationTuples(x)
    columns <- map(
        seq_len(nrow(tuples)),
        .fmaTupleColumns,
        x = x,
        tuples = tuples,
        refTable = refTable,
        types = types,
        coverage = coverage,
        minPurity = minPurity,
        pipThreshold = pipThreshold
    )
    mat <- .fmaDropEmptyColumns(exec(cbind, !!!columns))
    colTypes <- str_extract(colnames(mat), "[^|]+$")
    meta <- data.frame(
        name = colnames(mat),
        tier = rep(tier, ncol(mat)),
        type = if_else(
            is_in(colTypes, c("maxCpp", "sumPip")),
            "continuous",
            "binary"
        )
    )
    AnnotationMatrix(
        mat,
        GenomicRanges::granges(ldRef, use.mcols = FALSE),
        meta,
        genome = NULL
    )
}

# One unit's columns over the reference's variants: its traits' PIPs and
# credible-set members, its variants aligned to the reference in one pass by
# the matching engine (.alignTablesToPanel()).
# @noRd
.fmaTupleColumns <- function(
    k,
    x,
    tuples,
    refTable,
    types,
    coverage,
    minPurity,
    pipThreshold
) {
    rows <- .enrMatchRows(x, .enrAnnotationIdent(tuples, k))
    gr <- .fmrGatherVariants(x, rows)
    md <- mcols(gr, use.names = FALSE)
    at <- .alignTablesToPanel(.variantTable(gr), refTable)
    values <- as.numeric(md$pip %||% rep(NA_real_, length(gr)))
    hit <- !is.na(at) & !is.na(values)
    M <- .variantTableSize(refTable)
    maxPip <- .fmaScatter(M, at[hit], values[hit], "max")
    out <- list(
        maxCpp = maxPip,
        sumPip = .fmaScatter(M, at[hit], values[hit], "sum"),
        credibleSet = if (is_in("credibleSet", types)) {
            .fmaCredibleSetColumn(md, at, M, coverage, minPurity)
        },
        highPip = as.numeric(maxPip >= pipThreshold)
    )[types]
    mat <- exec(cbind, !!!out)
    colnames(mat) <- str_c(.fmaUnitName(tuples, k), "|", types)
    mat
}

# A reference-length vector holding the largest (`f = "max"`) or the sum
# (`"sum"`) of the values at each reference index `at`, 0 where none falls.
# @noRd
.fmaScatter <- function(M, at, values, f) {
    out <- numeric(M)
    if (length(at) == 0L) {
        return(out)
    }
    if (f == "sum") {
        sums <- rowsum(values, at)
        out[as.integer(rownames(sums))] <- sums[, 1L]
        return(out)
    }
    o <- order(at, -values)
    first <- o[!duplicated(at[o])]
    out[at[first]] <- values[first]
    out
}

# 1 for the reference variants (at their indices `at`) in a credible set of
# any of the unit's traits at `coverage`, as credibleSets() reads the
# credible-set columns `md` hold.
# @noRd
.fmaCredibleSetColumn <- function(md, at, M, coverage, minPurity) {
    out <- numeric(M)
    csCol <- str_c("cs_", coverage * 100)
    if (!is_in(csCol, colnames(md))) {
        return(out)
    }
    labels <- as.character(md[[csCol]])
    inCs <- !is.na(labels) &
        str_length(labels) > 0L &
        !str_detect(labels, "_0$")
    purCol <- str_c(csCol, "_purity")
    hasPurity <- is_in(purCol, colnames(md))
    if (!is.null(minPurity) && !hasPurity) {
        warn(glue(
            "credibleSets: no purity column '{purCol}' for coverage ",
            "{coverage}; minPurity filter skipped."
        ))
    }
    if (!is.null(minPurity) && hasPurity) {
        purity <- as.numeric(md[[purCol]])
        inCs <- inCs & !is.na(purity) & purity >= minPurity
    }
    members <- at[inCs]
    replace(out, members[!is.na(members)], 1)
}

# The unit's name: its study, context and method, without an axis the
# collection does not have.
# @noRd
.fmaUnitName <- function(tuples, k) {
    axes <- unlist(tuples[k, c("qtlStudy", "qtlContext", "qtlMethod")])
    str_flatten(axes[!is.na(axes)], "|")
}

# The columns with at least one non-zero entry, warning about the rest.
# @noRd
.fmaDropEmptyColumns <- function(mat) {
    empty <- colSums(mat != 0) == 0
    if (!any(empty)) {
        return(mat)
    }
    warn(glue(
        "fineMappingAnnotations: dropping {sum(empty)} column(s) that are ",
        "zero for every variant of the LD reference: ",
        "{str_flatten(utils::head(colnames(mat)[empty], 5L), ', ')}."
    ))
    mat[, !empty, drop = FALSE]
}
