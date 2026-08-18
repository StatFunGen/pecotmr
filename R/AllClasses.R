# =============================================================================
# AllClasses.R
# -----------------------------------------------------------------------------
# Virtual base classes shared across the package. Concrete subclasses live
# in their own per-class files (QtlSumStats.R, GwasSumStats.R, QtlDataset.R,
# QtlFineMappingResult.R, GwasFineMappingResult.R, etc.).
#
# Per Bioconductor convention this file is loaded first in the Collate
# ordering (the "AllClasses.R" filename sorts to the top of the alphabet),
# and every method-bearing file uses `@include AllClasses.R` so roxygen
# topologically orders the Collate field for us.
# =============================================================================

#' @include AllGenerics.R GenotypeHandle.R
#' @importFrom methods setClass setMethod new is validObject
NULL

# =============================================================================
# SumStatsBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS summary statistics collections.
# Concrete subclasses (QtlSumStats, GwasSumStats) inherit from DFrame and
# share the ldSketch / genome / qcInfo slots. getZ / getN / getMaf / nSnps are
# defined once on SumStatsBase (they only delegate to getSumStats); subsetChr /
# getVarY / getSumStats / getSumstatDf stay on the concrete subclass because
# they rely on the tuple shape (3-tuple QtlSumStats, 1-tuple GwasSumStats).
# =============================================================================

#' @title Summary Statistics Base Class
#' @description Virtual base class for QTL and GWAS summary statistics
#'   collections. Concrete subclasses (\code{QtlSumStats}, \code{GwasSumStats})
#'   inherit from \code{DFrame} and share the \code{ldSketch} / \code{genome} /
#'   \code{qcInfo} slots.
#' @slot ldSketch The \code{GenotypeHandle} the QC pipeline harmonized against,
#'   or \code{NULL}. Optional: LD-free workflows (e.g. mash, which operates
#'   across conditions per variant) carry \code{NULL}; pipelines that need LD
#'   validate its presence when they consume the collection.
#' @slot genome Character, genome build label.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty \code{list()}
#'   on construction; populated by \code{summaryStatsQc()} with a per-step audit
#'   record (filter names, drop counts, liftover target, RAISS settings, etc.).
#'   Fine-mapping and TWAS-weights pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0L} -- the slot serves as both the gating
#'   flag and the audit trail.
#' @export
setClass(
    "SumStatsBase",
    contains = c("VIRTUAL", "DFrame"),
    representation(
        ldSketch = "ANY",
        genome = "character",
        qcInfo = "list"
    )
)

#' @rdname getGenome
#' @examples
#' data(qtlSumStatsExample)
#' getGenome(qtlSumStatsExample)
#' @export
setMethod("getGenome", "SumStatsBase", function(x, ...) x@genome)

#' @rdname getQcInfo
#' @examples
#' data(qtlSumStatsExample)
#' getQcInfo(qtlSumStatsExample)
#' @export
setMethod("getQcInfo", "SumStatsBase", function(x, ...) x@qcInfo)

#' @rdname getQcDiagnostics
#' @examples
#' data(qtlSumStatsExample)
#' getQcDiagnostics(qtlSumStatsExample)
#' @export
setMethod("getQcDiagnostics", "SumStatsBase", function(x, entry = 1L, ...) {
    qc <- x@qcInfo
    if (length(qc) == 0L) {
        return(NULL)
    }
    audits <- qc$entryAudit
    if (is.null(audits)) {
        return(NULL)
    }
    if (is.null(entry)) {
        out <- map(audits, "ldMismatchDiagnostics")
        keep <- !map_lgl(out, is.null)
        if (!any(keep)) {
            return(NULL)
        }
        set_names(out[keep], seq_along(audits)[keep])
    } else {
        if (
            !is.numeric(entry) ||
                length(entry) != 1L ||
                entry < 1L ||
                entry > length(audits)
        ) {
            msg <- glue(
                "`entry` must be a single integer in 1:{length(audits)}."
            )
            abort(msg)
        }
        audits[[as.integer(entry)]]$ldMismatchDiagnostics
    }
})

#' @rdname getLdSketch
#' @examples
#' data(qtlSumStatsExample)
#' getLdSketch(qtlSumStatsExample)
#' @export
setMethod("getLdSketch", "SumStatsBase", function(x, ...) x@ldSketch)

#' @rdname getStudy
#' @examples
#' data(qtlDatasetExample)
#' getStudy(qtlDatasetExample)
#' @export
setMethod("getStudy", "SumStatsBase", function(x) unique(as.character(x$study)))

# getZ / getN / getMaf / nSnps delegate purely to getSumStats (which the
# concrete subclass dispatches with its own tuple shape), so they live once on
# the base and forward `...` through.
#' @rdname getZ
#' @examples
#' data(qtlSumStatsExample)
#' getZ(qtlSumStatsExample)
#' @export
setMethod("getZ", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$Z)

#' @rdname getN
#' @examples
#' data(qtlSumStatsExample)
#' getN(qtlSumStatsExample)
#' @export
setMethod("getN", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$N)

# getP / getBeta / getSe are first-class alongside getZ / getN: they read the
# optional P / BETA / SE mcols and return NULL when the entry does not carry
# them (DataFrame `$` semantics), so a p-value-primary sumstats (e.g. TensorQTL
# cis output) is an equal citizen to a Z-primary GWAS sumstats.
#' @rdname getP
#' @examples
#' data(qtlSumStatsExample)
#' getP(qtlSumStatsExample)
#' @export
setMethod("getP", "SumStatsBase", function(x, ...) mcols(getSumStats(x, ...))$P)

#' @rdname getBeta
#' @examples
#' data(qtlSumStatsExample)
#' getBeta(qtlSumStatsExample)
#' @export
setMethod("getBeta", "SumStatsBase", function(x, ...) {
    mcols(getSumStats(x, ...))$BETA
})

#' @rdname getSe
#' @examples
#' data(qtlSumStatsExample)
#' getSe(qtlSumStatsExample)
#' @export
setMethod("getSe", "SumStatsBase", function(x, ...) {
    mcols(getSumStats(x, ...))$SE
})

#' @rdname getMaf
#' @examples
#' data(qtlDatasetExample)
#' getMaf(qtlDatasetExample)
#' @export
setMethod("getMaf", "SumStatsBase", function(x, ...) {
    mc <- mcols(getSumStats(x, ...))
    if (is_in("MAF", colnames(mc))) mc$MAF else NULL
})

#' @rdname nSnps
#' @examples
#' data(qtlSumStatsExample)
#' nSnps(qtlSumStatsExample)
#' @export
setMethod("nSnps", "SumStatsBase", function(x, ...) length(getSumStats(x, ...)))

# =============================================================================
# FineMappingResultBase
# -----------------------------------------------------------------------------
# Shared parent of the QTL and GWAS fine-mapping result collections.
# Concrete subclasses (QtlFineMappingResult, GwasFineMappingResult) carry
# a DFrame of per-fit rows plus a shared ldSketch slot. Downstream
# pipelines dispatch on FineMappingResultBase for behaviors that apply to
# either flavour, and on the concrete subclass when the tuple shape
# matters.
# =============================================================================

#' @title Fine-Mapping Result Base Class
#' @description Virtual base class for fine-mapping result collections. Concrete
#'   subclasses (\code{QtlFineMappingResult}, \code{GwasFineMappingResult})
#'   carry a \code{DFrame} of per-fit rows and a shared \code{ldSketch} slot.
#'   Downstream pipelines should dispatch on \code{FineMappingResultBase} for
#'   behaviors that apply to either flavour, and on the concrete subclass when
#'   the tuple shape matters.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the fits were computed
#'   against, or \code{NULL} when the fits were derived from individual-level
#'   data (no LD reference). Used downstream for cross-pipeline LD-sketch
#'   identity validation.
#' @export
setClass(
    "FineMappingResultBase",
    contains = c("VIRTUAL", "DFrame"),
    representation(ldSketch = "ANY")
)

#' @rdname getStudy
#' @export
setMethod("getStudy", "FineMappingResultBase", function(x) {
    unique(as.character(x$study))
})

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "FineMappingResultBase", function(x, ...) x@ldSketch)

#' @rdname getMethodNames
#' @examples
#' data(qtlFineMappingExample)
#' getMethodNames(qtlFineMappingExample)
#' @export
setMethod("getMethodNames", "FineMappingResultBase", function(x) {
    unique(as.character(x$method))
})

#' @rdname adjustPips
#' @export
setMethod(
    "adjustPips",
    "FineMappingResultBase",
    function(x, keepVariants, ...) {
        if (nrow(x) == 0L) {
            return(x)
        }
        entries <- x@listData$entry
        for (i in seq_along(entries)) {
            adj <- tryCatch(
                adjustPips(entries[[i]], keepVariants, ...),
                error = function(err) NULL
            )
            if (!is.null(adj)) entries[[i]] <- adj
        }
        x@listData$entry <- entries
        x
    }
)

# Select the single FineMappingEntry addressed by a (tuple / region) key. Each
# concrete subclass implements this with its own row selector; the delegating
# accessors below then live once on the base and route through it.
setGeneric(".fmrSelectEntry", function(x, ...) {
    standardGeneric(".fmrSelectEntry")
})

#' @rdname getCs
#' @export
setMethod(
    "getCs",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        coverage = 0.95,
        minPurity = NULL,
        ...
    ) {
        # Selectors pinning one entry -> that entry's bare credible-set table;
        # no /
        # partial selectors -> aggregate every matching entry's credible sets,
        # prefixed with the row identity (study/context/trait/region_id/method).
        # `minPurity` is an independent CS-quality filter, orthogonal to
        # coverage.
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = getCs,
            coverage = coverage,
            minPurity = minPurity
        )
    }
)

#' @rdname getTopLoci
#' @export
setMethod(
    "getTopLoci",
    "FineMappingResultBase",
    function(
        x,
        type = c("data.frame", "GRanges"),
        signalCutoff = 0.025,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        minPurity = NULL,
        ...
    ) {
        type <- arg_match(type)
        # type = "GRanges" is honored only for a single pinned entry; the
        # aggregate (identity-prefixed) form is data.frame-only.
        if (type == "GRanges") {
            return(.fmrbTopLociGranges(
                x,
                study,
                context,
                trait,
                method,
                region,
                signalCutoff,
                minPurity
            ))
        }
        # data.frame: bare per-variant table when selectors pin one entry, else
        # the matching rows' per-variant tables stacked with row-identity cols.
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = getTopLoci,
            type = "data.frame",
            signalCutoff = signalCutoff,
            minPurity = minPurity
        )
    }
)

# Resolve the single pinned entry and return its GRanges view (aggregating
# across multiple entries requires type = "data.frame").
# @noRd
.fmrbTopLociGranges <- function(
    x,
    study,
    context,
    trait,
    method,
    region,
    signalCutoff,
    minPurity
) {
    sel <- tryCatch(
        .fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ),
        error = function(e) e
    )
    if (inherits(sel, "error")) {
        msg <- glue(
            "getTopLoci: aggregating across multiple entries requires ",
            "type = 'data.frame'."
        )
        abort(msg)
    }
    getTopLoci(
        sel,
        type = "GRanges",
        signalCutoff = signalCutoff,
        minPurity = minPurity
    )
}

#' @rdname getMarginalEffects
#' @export
setMethod(
    "getMarginalEffects",
    "FineMappingResultBase",
    function(
        x,
        maxPval = NULL,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        # Selectors pinning one entry -> that entry's bare marginal table; no /
        # partial selectors -> aggregate every matching entry's marginals,
        # prefixed with the row identity (study/context/trait/region_id/method).
        .fmrAggregateView(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region,
            perEntry = getMarginalEffects,
            maxPval = maxPval
        )
    }
)

#' @rdname getRegion
#' @examples
#' data(qtlFineMappingExample)
#' getRegion(qtlFineMappingExample)
#' @export
setMethod("getRegion", "FineMappingResultBase", function(x, ...) {
    .getRegionColumn(x)
})

#' @rdname getTraitPosition
#' @examples
#' data(qtlDatasetExample)
#' getTraitPosition(qtlDatasetExample)
#' @export
setMethod("getTraitPosition", "FineMappingResultBase", function(x, ...) {
    .getTraitPosColumn(x)
})

#' @rdname getSusieFit
#' @export
setMethod(
    "getSusieFit",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        getSusieFit(.fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)

#' @rdname getVariantIds
#' @export
setMethod(
    "getVariantIds",
    "FineMappingResultBase",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        region = NULL,
        ...
    ) {
        getVariantIds(.fmrSelectEntry(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            region = region
        ))
    }
)
