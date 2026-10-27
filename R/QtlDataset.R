# =============================================================================
# QtlDataset S4 class
# -----------------------------------------------------------------------------
# Single-study individual-level QTL container, built on MultiAssayExperiment:
# one RangedSummarizedExperiment per QTL context, plus a `genotype` experiment
# whose assay reads lazily through a GenotypeHandle. It holds data only: the QC
# and residualization a run applies are the pipeline's settings
# (GenotypeFilterParam, ResidualizationParam), not stored here. The entry point
# for individual-level
# fine-mapping (fineMappingPipeline), TWAS weight learning
# (twasWeightsPipeline), and multi-study composition (MultiStudyQtlDataset).
# =============================================================================

#' @include AllGenerics.R
#' @include MethodParam.R
NULL

#' @title QTL Dataset (individual-level data for one study)
#' @description S4 container for a single QTL study's regional data,
#'   extending \code{MultiAssayExperiment}. Every QTL context is one
#'   \code{RangedSummarizedExperiment}: rows are molecular traits positioned
#'   by \code{rowRanges}, columns are samples, and per-context phenotype
#'   covariates sit in that experiment's \code{colData}. Alongside them a
#'   \code{genotype} experiment carries variants in its \code{rowRanges},
#'   dosages in a lazily-read \code{DelayedArray} assay, and
#'   genotype-derived covariates (e.g., ancestry PCs) in its own
#'   \code{colData}. The \code{sampleMap} records which samples each
#'   experiment observes, so contexts need not share a sample set.
#'
#'   Extending \code{MultiAssayExperiment} means the multi-assay surface
#'   applies directly: \code{experiments()}, \code{colData()},
#'   \code{sampleMap()}, and two-dimensional subsetting
#'   (\code{x[, samples, contexts]}), all of which preserve the class and its
#'   own slots. The genotype assay stays unread until an operation touches it.
#'
#'   One caveat comes with that laziness: \code{longForm()} and
#'   \code{wideFormat()} fail on a delayed assay, with
#'   \code{MultiAssayExperiment} reshaping it to zero rows and reporting
#'   "replacement has 0 rows". Subsetting does \emph{not} sidestep it,
#'   because selecting contexts keeps the genotype experiment. Reshape a
#'   context experiment on its own, or drop the genotype one by coercing
#'   first: \code{longForm(as(x, "MultiAssayExperiment")[, , "brain"])}.
#'   Reshaping genotypes is variants x samples and enormous on real data, so
#'   it is rarely what you want regardless; anything that materialises a tidy
#'   view of that experiment reads its dosages.
#'
#' @slot studyName Character (length 1). Study identifier; used in collection
#'   classes to tag downstream \code{QtlFineMappingResult} / \code{TwasWeights}
#'   entries.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
setClass(
    "QtlDataset",
    contains = "MultiAssayExperiment",
    representation(studyName = "character"),
    validity = function(object) .validateQtlDataset(object)
)

# The reserved experiment name holding genotype dosages; every other
# experiment in the MAE is a QTL context.
# @noRd
.QTL_GENO_EXPERIMENT <- "genotype"

# The per-context phenotype experiments as a plain named list, in
# experiment order. The genotype experiment is not a context.
# @noRd
.qtlPhenotypeList <- function(object) {
    exps <- MultiAssayExperiment::experiments(object)
    nms <- setdiff(names(exps), .QTL_GENO_EXPERIMENT)
    out <- as.list(exps)
    out[nms]
}

# Validity for QtlDataset: the study name, the phenotypes list, and
# cross-context trait-position consistency. Returns TRUE or an error vector.
# @noRd
.validateQtlDataset <- function(object) {
    errors <- c(
        .qtlValidateScalars(object),
        .qtlValidatePhenotypes(object),
        .qtlValidateTraitPositions(object)
    )
    if (length(errors) == 0) TRUE else errors
}

# The study name: one non-empty string.
# @noRd
#' @importFrom checkmate makeAssertCollection assertString
.qtlValidateScalars <- function(object) {
    coll <- makeAssertCollection()
    assertString(
        object@studyName,
        min.chars = 1L,
        .var.name = "studyName",
        add = coll
    )
    coll$getMessages()
}

# What is wrong with the context names, if anything. The checks are ordered:
# a list with no usable names is not also reported as having duplicates.
# @noRd
.qtlContextNameErrors <- function(contextNames) {
    if (
        is.null(contextNames) ||
            any(str_length(contextNames) == 0L, na.rm = TRUE) ||
            any(is.na(contextNames))
    ) {
        return("'phenotypes' must be a named list with non-empty names")
    }
    if (n_distinct(contextNames) < length(contextNames)) {
        return("context names in 'phenotypes' must be unique")
    }
    if (is_in(.QTL_GENO_EXPERIMENT, contextNames)) {
        return(glue(
            "'{.QTL_GENO_EXPERIMENT}' is reserved for the genotype ",
            "experiment and cannot name a context"
        ))
    }
    character(0)
}

# What is wrong with one context's SummarizedExperiment, if anything.
# @noRd
.qtlPhenotypeEntryErrors <- function(ctx, phenotypes) {
    se <- phenotypes[[ctx]]
    if (!methods::is(se, "SummarizedExperiment")) {
        return(glue(
            "phenotypes[[{ctx}]] must be a SummarizedExperiment ",
            "(got {class(se)[[1L]]})"
        ))
    }
    if (is.null(colnames(se))) {
        return(glue(
            "phenotypes[[{ctx}]] has no column names. Which samples ",
            "a context observes is recorded in the sampleMap, so ",
            "every context must name its samples"
        ))
    }
    character(0)
}

# Shape checks on the phenotype list handed to the constructor. Run before
# the MultiAssayExperiment is assembled, because a malformed list makes
# ExperimentList fail first with a message about experiments rather than
# about contexts. Returns an error vector.
# @noRd
.qtlCheckPhenotypeList <- function(phenotypes) {
    emptyError <- if (length(phenotypes) == 0L) {
        "'phenotypes' must not be empty"
    } else {
        character(0)
    }
    # Each check answers with its own message or nothing, so the report is a
    # concatenation rather than a vector appended to in place.
    entryErrors <- map(
        seq_along(phenotypes),
        .qtlPhenotypeEntryErrors,
        phenotypes = phenotypes
    )
    as.character(c(
        emptyError,
        .qtlContextNameErrors(names(phenotypes)),
        .qtlConcatChr(entryErrors)
    ))
}

# Concatenate per-item message vectors, empty-safe.
# @noRd
.qtlConcatChr <- function(pieces) {
    if (length(pieces) == 0L) {
        return(character(0))
    }
    as.character(list_c(pieces))
}

# The MAE must carry the genotype experiment plus at least one context.
# @noRd
.qtlValidatePhenotypes <- function(object) {
    exps <- MultiAssayExperiment::experiments(object)
    missingGeno <- if (is_in(.QTL_GENO_EXPERIMENT, names(exps))) {
        character(0)
    } else {
        glue("experiment '{.QTL_GENO_EXPERIMENT}' is missing")
    }
    noContexts <- if (length(.qtlPhenotypeList(object)) == 0L) {
        "'phenotypes' must not be empty"
    } else {
        character(0)
    }
    as.character(c(missingGeno, noContexts))
}

# TRUE when two GRanges share chrom + start + end.
# @noRd
.qtlSameRange <- function(prev, this) {
    identical(
        as.character(GenomicRanges::seqnames(prev)),
        as.character(GenomicRanges::seqnames(this))
    ) &&
        GenomicRanges::start(prev) == GenomicRanges::start(this) &&
        GenomicRanges::end(prev) == GenomicRanges::end(this)
}

# A trait shared across contexts must have consistent rowRanges everywhere.
# @noRd
.qtlValidateTraitPositions <- function(object) {
    pheno <- .qtlPhenotypeList(object)
    allSe <- length(pheno) > 1L &&
        all(map_lgl(pheno, methods::is, "SummarizedExperiment"))
    if (!allSe) {
        return(character())
    }
    # Every (trait, range) observation in context order, so each one can be
    # compared against the first sighting of its trait without a running map.
    pairs <- .qtlConcat(map(
        seq_along(pheno),
        .qtlTraitRangePairs,
        pheno = pheno
    ))
    tids <- map_chr(pairs, "tid")
    .qtlConcatChr(map(
        seq_along(pairs),
        .qtlTraitRangeError,
        tids = tids,
        ranges = map(pairs, "range")
    ))
}

# Concatenate per-item lists, empty-safe.
# @noRd
.qtlConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# @noRd
.qtlTraitRangePair <- function(i, ids, rr) {
    list(tid = ids[[i]], range = rr[i])
}

# One context's (trait, range) observations, or none when the rowRanges and
# the rownames disagree on length.
# @noRd
.qtlTraitRangePairs <- function(ctx, pheno) {
    se <- pheno[[ctx]]
    rr <- SummarizedExperiment::rowRanges(se)
    ids <- rownames(se)
    if (length(rr) != length(ids)) {
        return(list())
    }
    map(seq_along(ids), .qtlTraitRangePair, ids = ids, rr = rr)
}

# Observation `i` disagrees with the first sighting of the same trait. Earlier
# sightings are the reference, so the first one never reports.
# @noRd
.qtlTraitRangeError <- function(i, tids, ranges) {
    earlier <- which(tids[seq_len(i - 1L)] == tids[[i]])
    if (length(earlier) == 0L) {
        return(character(0))
    }
    if (.qtlSameRange(ranges[[earlier[[1L]]]], ranges[[i]])) {
        return(character(0))
    }
    glue("trait '{tids[[i]]}' has inconsistent rowRanges across contexts")
}

# =============================================================================
# QtlDataset constructor and accessors
# =============================================================================

# Check the constructor's inputs and return the genotype source as given.
# Every phenotype must be a usable SummarizedExperiment, and the genotype
# source must be an unsubset panel (or a bare handle) -- see
# .qtlCheckWholePanel for why the "unsubset" part matters.
# @noRd
.qtlValidateInputs <- function(phenotypes, genotypes) {
    errors <- .qtlCheckPhenotypeList(phenotypes)
    if (length(errors) > 0L) {
        abort(str_flatten(errors, "\n"))
    }
    source <- .openGenotypeSource(genotypes)
    if (is.null(source)) {
        abort(glue(
            "'genotypes' must be a genotype panel from readGenotypes() ",
            "(got {class(genotypes)[[1L]]})"
        ))
    }
    .qtlCheckWholePanel(genotypes)
    .qtlCheckRealGenotypes(genotypes)
    source
}

# A stochastic LD sketch is a panel whose columns are pseudo-samples
# projected from real ones: it carries LD, not anyone's genotypes, so it
# cannot be the genotype experiment a QTL study is fitted on.
# @noRd
.qtlCheckRealGenotypes <- function(genotypes) {
    if (!.ldRefIsProjection(genotypes)) {
        return(invisible(NULL))
    }
    abort(c(
        "'genotypes' is a stochastic LD sketch, not a study's genotypes.",
        i = "Its samples are projections; use it as an `ldReference` instead."
    ))
}

# A panel that has already been narrowed is refused rather than taken as the
# study's genotypes: the dataset is built from the whole panel and narrowed
# afterwards (x[i, j, ]), so what it holds is not decided by which panel
# happened to be passed.
# @noRd
.qtlCheckWholePanel <- function(genotypes) {
    if (!methods::is(genotypes, "RangedSummarizedExperiment")) {
        return(invisible(NULL))
    }
    full <- .panelFullDims(genotypes)
    got <- c(nrow(genotypes), ncol(genotypes))
    if (identical(as.integer(got), as.integer(full))) {
        return(invisible(NULL))
    }
    abort(glue(
        "'genotypes' is a subset panel ({got[[1L]]} x {got[[2L]]} of ",
        "{full[[1L]]} x {full[[2L]]}); pass the whole panel and subset the ",
        "QtlDataset afterwards, with x[i, j, ]."
    ))
}

#' @title Create a QtlDataset Object
#' @description Construct a \code{QtlDataset}: one study's individual-level
#'   QTL data as a \code{MultiAssayExperiment}. The genotype handle becomes a
#'   \code{genotype} experiment whose dosage assay reads lazily through it,
#'   and each named phenotype \code{SummarizedExperiment} becomes one
#'   context experiment. The \code{sampleMap} is derived from the column
#'   names actually present in each, so contexts observing different sample
#'   subsets are recorded rather than assumed away.
#' @param studyName Character (length 1). Study identifier.
#' @param genotypes A genotype panel (see \code{\link{readGenotypes}}).
#' @param phenotypes Named list of \code{SummarizedExperiment} objects, keyed by
#'   context. Each SE must have \code{rowRanges} carrying trait positions and
#'   \code{colData} carrying per-context phenotype covariates. The name
#'   \code{"genotype"} is reserved.
#' @param genotypeCovariates Numeric matrix of genotype-derived covariates
#'   (e.g., ancestry PCs); rows are samples. Becomes the \code{colData} of the
#'   genotype experiment.
#' @details A \code{QtlDataset} holds data only. Which variants and samples
#'   an analysis keeps, and whether residuals are scaled, are the pipeline's
#'   settings (\code{\link{GenotypeFilterParam}} and
#'   \code{\link{ResidualizationParam}} on \code{\link{fineMappingPipeline}},
#'   \code{\link{twasWeightsPipeline}} and \code{\link{colocboostPipeline}}),
#'   applied per extracted block. To work with part of the data, subset the
#'   dataset itself: \code{x[i, j, k]} narrows features, samples and
#'   contexts.
#' @return A \code{QtlDataset} object.
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr")
#' )
#' rng <- GenomicRanges::GRanges(
#'   "chr22", IRanges::IRanges(14600000L, width = 1000L)
#' )
#' names(rng) <- "ENSG1"
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(expression = matrix(
#'     rnorm(ncol(panel)), 1,
#'     dimnames = list("ENSG1", colnames(panel))
#'   )),
#'   rowRanges = rng
#' )
#' QtlDataset(
#'     studyName = "s1", genotypes = panel, phenotypes = list(brain = se)
#' )
#' @export
QtlDataset <- function(
    studyName,
    genotypes,
    phenotypes,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
) {
    handle <- .qtlValidateInputs(phenotypes, genotypes)
    # Every experiment's rowRanges in UCSC style, so the dataset's traits and
    # variants meet every other range pecotmr holds exactly.
    experiments <- map(
        c(
            set_names(
                list(.genotypeExperiment(genotypes, genotypeCovariates)),
                .QTL_GENO_EXPERIMENT
            ),
            as.list(phenotypes)
        ),
        .ucscExperiment
    )
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiments,
        colData = .qtlPrimaryColData(experiments),
        sampleMap = .qtlSampleMap(experiments)
    )
    methods::new(
        "QtlDataset",
        mae,
        studyName = as.character(studyName)
    )
}

#' @describeIn QtlDataset-class Subset by feature, sample and experiment.
#'   Selecting experiments selects among \emph{contexts}: the genotype
#'   experiment is the substrate every context is interpreted against rather
#'   than one of the things being chosen between, so it is always retained.
#' @param x A \code{QtlDataset}.
#' @param i,j,k Feature, sample and experiment subscripts, as for
#'   \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}.
#' @param ... Passed on to row subsetting.
#' @param drop Passed through to \code{MultiAssayExperiment}: when
#'   \code{TRUE}, experiments the subset leaves empty are removed.
#' @return A \code{QtlDataset} narrowed to the requested features, samples
#'   and contexts.
#' @export
setMethod("[", "QtlDataset", function(x, i, j, k, ..., drop = FALSE) {
    if (missing(k)) {
        return(methods::callNextMethod())
    }
    # The experiment axis is applied here rather than by rewriting `k` and
    # deferring: setMethod() moves a method with extra formals into a
    # generated .local(), and callNextMethod() from inside it re-dispatches
    # on the ORIGINAL arguments, so a reassigned `k` never arrives. The
    # remaining axes then go back through this method with no experiment
    # subscript, which does reach the inherited one.
    nms <- names(MultiAssayExperiment::experiments(x))
    out <- MultiAssayExperiment::subsetByAssay(
        x,
        .qtlSelectExperiments(k, nms)
    )
    if (missing(i) && missing(j)) {
        return(out)
    }
    if (missing(i)) {
        return(out[, j, ])
    }
    if (missing(j)) {
        return(out[i, , ])
    }
    out[i, j, ]
})

#' @describeIn QtlDataset-class Reshape the measurements into one long
#'   \code{DataFrame}. Only the QTL contexts are reshaped: the genotype
#'   experiment is the substrate they are interpreted against rather than a
#'   measurement, and being variants x samples it would dwarf them. Pass
#'   \code{genotype = TRUE} to include it, which reads the dosages into
#'   memory -- \code{MultiAssayExperiment} cannot reshape a delayed assay,
#'   and fails with "replacement has 0 rows" if asked to.
#'
#'   \code{wideFormat()} has no such method: it is a plain function rather
#'   than a generic, and it reshapes the \code{ExperimentList} directly, so
#'   there is nothing to dispatch on. Coerce first --
#'   \code{wideFormat(as(x, "MultiAssayExperiment")[, , contexts])}.
#' @param object A \code{QtlDataset}.
#' @param genotype Logical (length 1), default \code{FALSE}. Include the
#'   genotype experiment, reading its dosages into memory to do so.
#' @return A long-format \code{DataFrame}, one row per
#'   (assay, primary, rowname) observation.
#' @importFrom MultiAssayExperiment longForm
#' @export
setMethod("longForm", "QtlDataset", function(object, ..., genotype = FALSE) {
    MultiAssayExperiment::longForm(
        .qtlReshapeSource(object, genotype),
        ...
    )
})

# The object a reshape runs over. Dropping the genotype experiment is both
# the useful default -- a caller asking for the measurements does not mean
# every dosage -- and the working one, since a delayed assay cannot be
# reshaped at all.
# @noRd
.qtlReshapeSource <- function(x, genotype) {
    mae <- methods::as(x, "MultiAssayExperiment")
    if (isTRUE(genotype)) {
        return(.qtlRealizeDosages(mae))
    }
    .qtlMuffleDrop(MultiAssayExperiment::subsetByAssay(mae, contexts(x)))
}

# Silence the warning and the message the drop is guaranteed to raise --
# MultiAssayExperiment announces the same non-news twice. Omitting the
# genotype experiment is this method's documented behaviour, so saying so on
# every call is noise. Anything else still gets through.
# @noRd
.qtlMuffleDrop <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            if (str_detect(conditionMessage(w), "'experiments' dropped")) {
                invokeRestart("muffleWarning")
            }
        },
        message = function(m) {
            if (str_detect(conditionMessage(m), "harmonizing input")) {
                invokeRestart("muffleMessage")
            }
        }
    )
}

# Read the dosages into memory so the reshape can see them.
# @noRd
.qtlRealizeDosages <- function(mae) {
    lazy <- MultiAssayExperiment::experiments(mae)
    bare <- lazy[[.QTL_GENO_EXPERIMENT]]
    se <- SummarizedExperiment::`assay<-`(
        bare,
        "dosage",
        value = as.matrix(SummarizedExperiment::assay(bare, "dosage"))
    )
    exps <- `[[<-`(lazy, .QTL_GENO_EXPERIMENT, value = se)
    MultiAssayExperiment::`experiments<-`(mae, value = exps)
}

# Resolve an experiment subscript to names, always keeping the genotype
# experiment. Dropping it would leave an object that fails its own validity
# and has silently lost the genotype covariates, since those live in that
# experiment's colData.
# @noRd
.qtlSelectExperiments <- function(k, nms) {
    sel <- nms[.qtlExperimentIndex(k, nms)]
    contexts <- setdiff(sel, .QTL_GENO_EXPERIMENT)
    if (length(contexts) == 0L) {
        known <- setdiff(nms, .QTL_GENO_EXPERIMENT)
        abort(glue(
            "subscript selects no QTL context; a QtlDataset must keep at ",
            "least one of: {str_flatten(known, ', ')}"
        ))
    }
    unique(c(.QTL_GENO_EXPERIMENT, contexts))
}

# Normalize a character / logical / numeric experiment subscript to indices.
# @noRd
.qtlExperimentIndex <- function(k, nms) {
    if (is.character(k)) {
        return(match(k, nms))
    }
    if (is.logical(k)) {
        return(which(rep(k, length.out = length(nms))))
    }
    as.integer(k)
}

# Narrow a dataset (or a bare MAE) to a sample set, keeping colData and
# sampleMap in step. Samples the object does not have are ignored rather
# than an error, matching how the retired keepSamples slot was intersected.
# @noRd
.qtlRestrictSamples <- function(x, keepSamples) {
    if (length(keepSamples) == 0L) {
        return(x)
    }
    ids <- rownames(MultiAssayExperiment::colData(x))
    # Index positionally. A character subscript matching nothing is an error
    # in MultiAssayExperiment, while an empty positional one is simply an
    # empty result -- and a keep set disjoint from the panel is a legitimate
    # request for no samples, not a mistake.
    x[, which(is_in(ids, as.character(keepSamples))), ]
}

# Replace the genotype handle by rebuilding the genotype experiment around
# it. The handle lives in exactly one place -- the assay's seed -- so there
# is no second copy to keep in step; genotypeHandle() reads it back.
# @noRd
.qtlWithGenotypeHandle <- function(x, handle) {
    exps <- MultiAssayExperiment::experiments(x)
    gCov <- .qtlColDataMatrix(exps[[.QTL_GENO_EXPERIMENT]])
    withHandle <- `[[<-`(
        exps,
        .QTL_GENO_EXPERIMENT,
        value = .genotypeExperiment(handle, gCov)
    )
    rebuilt <- MultiAssayExperiment::`experiments<-`(x, value = withHandle)
    validObject(rebuilt)
    rebuilt
}

# The primary sample table: every sample any experiment observes, in
# genotype-panel order first so the common case reads naturally.
# @noRd
.qtlPrimaryColData <- function(experiments) {
    ids <- reduce(map(experiments, colnames), union)
    S4Vectors::DataFrame(row.names = as.character(ids))
}

# The sampleMap: one row per (experiment, sample) pair actually present.
# Built from the column names rather than assumed, since contexts need not
# share a sample set.
# @noRd
.qtlSampleMap <- function(experiments) {
    parts <- imap(experiments, .qtlSampleMapPart)
    exec(rbind, !!!unname(parts))
}

# One experiment's slice of the sampleMap.
# @noRd
.qtlSampleMapPart <- function(se, name) {
    ids <- as.character(colnames(se))
    S4Vectors::DataFrame(
        assay = factor(rep(name, length(ids)), levels = name),
        primary = ids,
        colname = ids
    )
}

#' @rdname studyName
#' @export
setMethod("studyName", "QtlDataset", function(x) x@studyName)

#' @rdname contexts
#' @export
setMethod("contexts", "QtlDataset", function(x) {
    names(.qtlPhenotypeList(x))
})

#' @rdname genotypeCovariates
#' @export
setMethod("genotypeCovariates", "QtlDataset", function(x) {
    .qtlSuppliedCovariates(.qtlColDataMatrix(.qtlGenotypeSe(x)))
})

# The genotype experiment's colData has one row per panel sample, because a
# SummarizedExperiment cannot have fewer. A covariate matrix covering only
# some samples therefore leaves the rest NA throughout, and those rows mean
# "not supplied" rather than "supplied as missing". Dropping them keeps the
# residualization design over the samples that actually have covariates,
# which is what the alignment step intersects on.
# @noRd
.qtlSuppliedCovariates <- function(m) {
    if (ncol(m) == 0L || nrow(m) == 0L) {
        return(m)
    }
    m[rowSums(!is.na(m)) > 0L, , drop = FALSE]
}

# The genotype experiment.
# @noRd
.qtlGenotypeSe <- function(x) {
    MultiAssayExperiment::experiments(x)[[.QTL_GENO_EXPERIMENT]]
}

# An experiment's colData as the numeric samples x covariates matrix the
# residualization design expects. A colData with no columns still has to
# carry its sample names, so the design can align on them.
# @noRd
.qtlColDataMatrix <- function(se) {
    cd <- SummarizedExperiment::colData(se)
    `rownames<-`(as.matrix(as.data.frame(cd)), rownames(cd))
}

#' @rdname genotypeHandle
#' @keywords internal
setMethod("genotypeHandle", "QtlDataset", function(x) {
    # Derived, not stored. The handle already lives inside the genotype
    # assay's seed -- that is what lets the dosages read lazily -- so a
    # parallel slot was a second copy that had to be kept in step by hand.
    # Reading it back removes the invariant instead of policing it. NULL for
    # genotypes read through another backend (VCF, GDS) or held in memory.
    .ldRefHandle(.qtlGenotypeSe(x))
})

# The genotype experiment's variants, one row per panel row: SNP (the
# panel's own ids), CHR, BP, A1, A2. Read off its rowRanges, so
# genotypes from any backend answer alike.
# @noRd
.qtlGenotypeVariants <- function(x) {
    .rangesSnpInfo(SummarizedExperiment::rowRanges(.qtlGenotypeSe(x)))
}

# --- Internal: resolve the variant-selection region for the genotype handle.
# Returns a GRanges (one or more ranges). When `traitId` is supplied, expand
# each trait's rowRange by `cisWindow` bp and take the union span (per the
# multi-trait rule: `[min(start) - cisWindow, max(end) + cisWindow]`). When
# `region` is supplied it is taken literally and may contain multiple ranges
# (e.g. for joint multi-region extraction), optionally extended per-range by
# `cisWindow`. Exactly one of (traitId, region) may be supplied; if neither is,
# return NULL meaning "all variants in handle".
.qtlResolveVariantRegion <- function(
    x,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL
) {
    if (!is.null(traitId) && !is.null(region)) {
        abort("Specify either `traitId` or `region`, not both.")
    }
    if (is.null(traitId) && is.null(region)) {
        return(NULL)
    }
    if (!is.null(traitId)) {
        return(.qtlTraitRegion(x, traitId, cisWindow))
    }
    .qtlLiteralRegion(region, cisWindow)
}

# The requested traits' ranges within one context, or NULL when it carries
# none of them.
# @noRd
.qtlTraitRangesInContext <- function(ctx, x, traitId) {
    se <- .qtlContextSe(x, ctx)
    hits <- .qtlPresentIndices(traitId, rownames(se))
    if (length(hits) == 0L) {
        return(NULL)
    }
    SummarizedExperiment::rowRanges(se)[hits]
}

# Positions of `wanted` in `available`, dropping the ones not present.
# @noRd
.qtlPresentIndices <- function(wanted, available) {
    hits <- match(wanted, available)
    hits[!is.na(hits)]
}

# The union span (+/- cisWindow) of the requested traits' rowRanges across all
# contexts. Requires cisWindow and a single shared chromosome.
# @noRd
.qtlTraitRegion <- function(x, traitId, cisWindow) {
    if (is.null(cisWindow) || length(cisWindow) != 1L || cisWindow < 0) {
        msg <- glue(
            "`cisWindow` is required (and must be non-negative) when ",
            "`traitId` is specified."
        )
        abort(msg)
    }
    perTraitRanges <- compact(map(
        contexts(x),
        .qtlTraitRangesInContext,
        x = x,
        traitId = traitId
    ))
    if (length(perTraitRanges) == 0L) {
        abort("None of the requested traitId values were found in any context.")
    }
    allRanges <- exec(c, !!!perTraitRanges)
    chrs <- unique(as.character(GenomicRanges::seqnames(allRanges)))
    if (length(chrs) != 1L) {
        msg <- glue(
            "Multi-trait variant extraction requires all selected traits ",
            "to share a chromosome ",
            "(got: {str_flatten(chrs, ', ')})."
        )
        abort(msg)
    }
    spanStart <- max(1L, min(GenomicRanges::start(allRanges)) - cisWindow)
    spanEnd <- max(GenomicRanges::end(allRanges)) + cisWindow
    GenomicRanges::GRanges(
        seqnames = chrs,
        ranges = IRanges::IRanges(start = spanStart, end = spanEnd)
    )
}

# A literal `region` GRanges, each range optionally extended by cisWindow.
# @noRd
#' @importFrom checkmate assertNumber
.qtlLiteralRegion <- function(region, cisWindow) {
    region <- .assertRegionGRanges(region)
    if (length(region) == 0L) {
        abort("`region` must contain at least one range.")
    }
    if (is.null(cisWindow)) {
        return(region)
    }
    assertNumber(cisWindow, lower = 0)
    GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(region),
        ranges = IRanges::IRanges(
            start = pmax(1L, GenomicRanges::start(region) - cisWindow),
            end = GenomicRanges::end(region) + cisWindow
        )
    )
}

# Per-trait genomic position: each trait's OWN rowRanges span across contexts
# (union min-start to max-end), WITHOUT the cisWindow expansion. One range per
# traitId; a trait absent from every context is an error, as a lookup by an
# unknown name is for rowRanges().
# @noRd
.qtlTraitPos <- function(x, traitIds) {
    # Extract chrom/start/end as plain vectors per trait (union span across
    # contexts), then build ONE fresh GRanges at the end. Combining per-context
    # GRanges with do.call(c, .) can trip S4 seqinfo reconciliation in some
    # GenomeInfoDb builds, so we avoid it entirely.
    spans <- map(traitIds, .qtlTraitSpan, x = x, contexts = contexts(x))
    absent <- traitIds[map_lgl(spans, is.null)]
    if (length(absent) > 0L) {
        abort(glue(
            "traitPosition: trait(s) in no context of this QtlDataset: ",
            "{str_flatten(absent, ', ')}."
        ))
    }
    starts <- map_int(spans, "start")
    # `set_names()` is vector-only, so name the GRanges through `names<-`
    # applied as a function -- still a copy, no binding rewritten.
    `names<-`(
        GenomicRanges::GRanges(
            map_chr(spans, "chr"),
            IRanges::IRanges(
                start = starts,
                end = pmax(map_int(spans, "end"), starts)
            )
        ),
        traitIds
    )
}

# @noRd
.qtlFirstSeqname <- function(rr) {
    as.character(GenomicRanges::seqnames(rr))[[1L]]
}

# One trait's union span across every context that carries it. The chromosome
# is the last matching context's, which is what the running assignment left
# behind; the span is the widest across all of them.
# @noRd
.qtlTraitSpan <- function(tid, x, contexts) {
    ranges <- compact(map(
        contexts,
        .qtlTraitRangesInContext,
        x = x,
        traitId = tid
    ))
    if (length(ranges) == 0L) {
        return(NULL)
    }
    chrs <- map_chr(ranges, .qtlFirstSeqname)
    list(
        chr = chrs[[length(chrs)]],
        start = as.integer(min(map_dbl(ranges, .qtlRangeMinStart))),
        end = as.integer(max(map_dbl(ranges, .qtlRangeMaxEnd)))
    )
}

# @noRd
.qtlRangeMinStart <- function(rr) {
    min(GenomicRanges::start(rr))
}

# @noRd
.qtlRangeMaxEnd <- function(rr) {
    max(GenomicRanges::end(rr))
}

#' @rdname traitPosition
#' @export
setMethod("traitPosition", "QtlDataset", function(x, trait = NULL) {
    tids <- if (is.null(trait)) {
        unique(list_c(map(.qtlPhenotypeList(x), rownames)))
    } else {
        as.character(trait)
    }
    .qtlTraitPos(x, tids)
})

# Internal: map a GRanges region (one or more ranges) into 1-based snpIdx into
# the genotype experiment's rows. Indices are unioned across ranges in range
# order (first occurrence wins), so overlapping ranges contribute each
# variant once.
.qtlVariantIndices <- function(x, region = NULL) {
    variants <- .qtlGenotypeVariants(x)
    if (is.null(region)) {
        return(seq_len(nrow(variants)))
    }
    siChr <- variants$CHR
    bp <- variants$BP
    rChr <- as.character(GenomicRanges::seqnames(region))
    rStart <- GenomicRanges::start(region)
    rEnd <- GenomicRanges::end(region)
    if (length(region) == 0L) {
        return(integer(0))
    }
    unique(list_c(map(
        seq_along(region),
        .qtlRegionHitIndices,
        siChr = siChr,
        bp = bp,
        rChr = rChr,
        rStart = rStart,
        rEnd = rEnd
    )))
}

# Rows of the SNP table falling inside region `i`.
# @noRd
.qtlRegionHitIndices <- function(i, siChr, bp, rChr, rStart, rEnd) {
    which(siChr == rChr[i] & bp >= rStart[i] & bp <= rEnd[i])
}

# A GenotypeFilterParam() bundle with every field resolved: an unset field
# (absent from the bundle) takes its default -- 0 for each cutoff, no
# restriction for each keep-list, TRUE for keepIndel. This is the form a
# pipeline threads to every extraction it makes; a dataset records no filter
# of its own, so the accessors extract with the all-default one.
# @noRd
.qtlResolveFilter <- function(genotypeFilterParam = GenotypeFilterParam()) {
    .assertMethodParam(genotypeFilterParam, "GenotypeFilterParam")
    list(
        mafCutoff = genotypeFilterParam$mafCutoff %||% 0,
        macCutoff = genotypeFilterParam$macCutoff %||% 0,
        xvarCutoff = genotypeFilterParam$xvarCutoff %||% 0,
        imissCutoff = genotypeFilterParam$imissCutoff %||% 0,
        keepSamples = genotypeFilterParam$keepSamples %||% character(0),
        keepVariants = genotypeFilterParam$keepVariants %||% character(0),
        keepIndel = genotypeFilterParam$keepIndel %||% TRUE
    )
}

# Drop samples whose missingness across the block exceeds the filter's
# imissCutoff. A cutoff of 0 (or an empty block) keeps every sample.
# @noRd
.qtlDropMissingSamples <- function(dosage, filter) {
    cutoff <- filter$imissCutoff
    if (cutoff <= 0 || nrow(dosage) == 0L || ncol(dosage) == 0L) {
        return(dosage)
    }
    dosage[rowMeans(is.na(dosage)) <= cutoff, , drop = FALSE]
}

# Internal: extract the dosage block (samples x variants) for the variants
# in `window` (every variant when NULL), narrow to `samples`, and apply the
# resolved genotype `filter` (.qtlResolveFilter()): per-sample imiss, then
# per-variant max(mafCutoff, macCutoff / (2 * n)) and xvarCutoff. The
# accessors extract the whole dataset with the all-default filter; a
# pipeline passes the window, samples and filter its run settled on, so
# genotypes, allele frequencies and residualized genotypes share one
# selection.
#
# Returns a list:
#   geno       : numeric matrix (kept samples x kept variants)
#   variantIds : character vector of kept variant IDs (= colnames(geno))
#   sampleIds  : character vector of kept sample IDs (= rownames(geno))
#   maf        : numeric vector of per-variant MAF for kept variants
#   af         : numeric vector of per-variant effect-allele (A1) frequency
#                for kept variants. Directional (NOT folded to the minor
#                allele): the frequency of the dosage-counted allele, which
#                is A1 by the same convention the marginal betas use. `maf`
#                is `pmin(af, 1 - af)`.
.qtlExtractBlock <- function(
    x,
    window = NULL,
    samples = NULL,
    filter = .qtlResolveFilter()
) {
    # A pipeline threads its resolved filter to every extraction; a NULL one
    # means a call site was missed, which must not pass for "no filter".
    if (!is.list(filter) || is.null(filter$mafCutoff)) {
        abort(
            "No resolved genotype filter reached this extraction.",
            .internal = TRUE
        )
    }
    inRegion <- .qtlVariantIndices(x, window)
    if (length(inRegion) == 0L) {
        return(.qtlEmptyBlockAllSamples(x))
    }
    # Apply keepVariants + indel restrictions before materialization so we do
    # not extract dosage we would immediately drop.
    snpIdx <- .qtlNarrowSnpIdx(x, inRegion, filter)
    if (length(snpIdx) == 0L) {
        return(.qtlEmptyBlock())
    }
    se <- .qtlGenotypeSe(x)
    # Labelled with the panel's own variant ids, as the rows a block reports.
    dosage <- `colnames<-`(
        .ldRefDosage(se, snpIdx, meanImpute = FALSE),
        .ldRefVariantIds(se)[snpIdx]
    )
    keep <- .qtlResolveSamples(dosage, x, samples)
    if (length(keep) == 0L) {
        return(.qtlEmptyBlockNoSamples(dosage))
    }
    kept <- dosage[keep, , drop = FALSE]
    filtered <- .qtlVariantFilters(
        .qtlDropMissingSamples(kept, filter),
        filter
    )
    imputed <- .qtlMeanImpute(filtered$dosage)
    list(
        geno = imputed,
        variantIds = colnames(imputed),
        sampleIds = rownames(imputed),
        maf = filtered$maf,
        af = filtered$af
    )
}

# Empty block preserving the full panel sample set (no variants selected).
# @noRd
.qtlEmptyBlockAllSamples <- function(x) {
    samples <- as.character(colnames(.qtlGenotypeSe(x)))
    list(
        geno = matrix(
            numeric(0),
            nrow = length(samples),
            ncol = 0L,
            dimnames = list(samples, character(0))
        ),
        variantIds = character(0),
        sampleIds = samples,
        maf = numeric(0),
        af = numeric(0)
    )
}

# Fully empty block (no variants, no samples).
# @noRd
.qtlEmptyBlock <- function() {
    list(
        geno = matrix(
            numeric(0),
            nrow = 0L,
            ncol = 0L,
            dimnames = list(character(0), character(0))
        ),
        variantIds = character(0),
        sampleIds = character(0),
        maf = numeric(0),
        af = numeric(0)
    )
}

# Empty block preserving the selected variants (no samples survived).
# @noRd
.qtlEmptyBlockNoSamples <- function(dosage) {
    list(
        geno = dosage[integer(0), , drop = FALSE],
        variantIds = colnames(dosage),
        sampleIds = character(0),
        maf = rep(NA_real_, ncol(dosage)),
        af = rep(NA_real_, ncol(dosage))
    )
}

# Narrow the selected variant indices by keepVariants (matched by chrom/pos/
# allele) and, unless kept, by dropping indels. Done before materialization.
# @noRd
.qtlNarrowSnpIdx <- function(x, snpIdx, filter) {
    si <- .qtlGenotypeVariants(x)
    kept <- if (length(filter$keepVariants) == 0L) {
        snpIdx
    } else {
        snpAll <- si$SNP[snpIdx]
        km <- matchVariants(snpAll, as.character(filter$keepVariants))
        snpIdx[replace(logical(length(snpAll)), km$idxA, TRUE)]
    }
    if (length(kept) == 0L || isTRUE(filter$keepIndel)) {
        return(kept)
    }
    # which() (not the mask) so an NA mask drops the variant rather than
    # injecting an NA index.
    snpMask <- str_length(as.character(si$A1[kept])) == 1L &
        str_length(as.character(si$A2[kept])) == 1L
    kept[which(snpMask)]
}

# Resolve the sample set: panel samples intersected with the dataset's
# primary colData -- which is what narrows a subset dataset -- and with the
# per-call `samples` arg.
# @noRd
.qtlResolveSamples <- function(dosage, x, samples) {
    inDataset <- intersect(
        rownames(dosage),
        rownames(MultiAssayExperiment::colData(x))
    )
    keep <- if (is.null(samples)) {
        inDataset
    } else {
        intersect(inDataset, as.character(samples))
    }
    keep
}

# Per-variant MAF / MAC / X-variance filters against the post-narrowing sample
# count. MAF is computed from the un-imputed dosage (A1 = the effect allele
# gives directional `af`; `maf` folds to the minor allele). Returns
# list(dosage, maf, af).
# @noRd
.qtlVariantFilters <- function(dosage, filter) {
    # No zero-variant guard: the only caller returns early when snpIdx is
    # empty, and extraction preserves the column count (.restoreRequestedOrder
    # aborts on a mismatch). The body is 0-column safe regardless -- colSums
    # gives numeric(0) and the mask logical(0), yielding the same empty result.
    nSamp <- nrow(dosage)
    nObs <- colSums(!is.na(dosage))
    sumD <- colSums(dosage, na.rm = TRUE)
    p <- if_else(nObs > 0L, sumD / (2 * nObs), NA_real_)
    afVec <- p
    mafVec <- pmin(p, 1 - p)
    effectiveMaf <- max(
        filter$mafCutoff,
        if (nSamp > 0L) filter$macCutoff / (2 * nSamp) else 0
    )
    byMaf <- !is.na(mafVec) & mafVec >= effectiveMaf
    keepVarMask <- if (filter$xvarCutoff <= 0 || nSamp <= 1L) {
        byMaf
    } else {
        mu <- if_else(nObs > 0L, sumD / nObs, 0)
        # mu is finite everywhere, so the centered NAs are exactly the
        # dosage NAs; they contribute nothing to the variance.
        centered <- replace(
            sweep(dosage, 2L, mu, FUN = "-"),
            is.na(dosage),
            0
        )
        varVec <- colSums(centered * centered) / (nSamp - 1L)
        byMaf & varVec >= filter$xvarCutoff
    }
    list(
        dosage = dosage[, keepVarMask, drop = FALSE],
        maf = mafVec[keepVarMask],
        af = afVec[keepVarMask]
    )
}

# Mean-impute remaining missing dosage cells (per column) so downstream linear
# algebra is well-defined; MAF was already computed pre-imputation.
# @noRd
.qtlMeanImpute <- function(dosage) {
    if (!anyNA(dosage)) {
        return(dosage)
    }
    naMask <- is.na(dosage)
    # An all-NA column means NaN either way, matching mean(numeric(0)).
    means <- colMeans(dosage, na.rm = TRUE)
    # One fill over the whole matrix instead of a copy per column: dosage is
    # variants x samples, so rebuilding it column by column is the expensive
    # way to say this.
    replace(dosage, naMask, means[col(dosage)[naMask]])
}

#' @rdname genotypes
#' @export
setMethod("genotypes", "QtlDataset", function(x) {
    # The whole dataset, unfiltered: narrow it first with x[i, j, ].
    .qtlExtractBlock(x)$geno
})

#' @rdname maf
#' @export
setMethod("maf", "QtlDataset", function(x) {
    block <- .qtlExtractBlock(x)
    set_names(block$maf, block$variantIds)
})

#' @rdname af
#' @export
setMethod("af", "QtlDataset", function(x) {
    block <- .qtlExtractBlock(x)
    set_names(block$af, block$variantIds)
})

#' @rdname molecularTraits
#' @export
setMethod(
    "molecularTraits",
    "QtlDataset",
    function(x, context, trait = NULL) {
        .qtlMolecularTraits(x, context, trait)
    }
)

# molecularTraits() for the pipelines, which may also keep only the traits
# positioned inside an analysis `window` (a GRanges) and apply the run's
# resolved phenotype QC (.qtlResolvePhenotypeFilter()). The accessor passes
# neither, so it returns the data as stored. Always a list keyed by context.
# @noRd
.qtlMolecularTraits <- function(
    x,
    contexts,
    traitId = NULL,
    window = NULL,
    phenotypeFilter = .qtlKeepAllPhenotypes()
) {
    .qtlValidateContexts(x, contexts)
    .qtlPhenotypeList(x)[contexts] |>
        .qtlFilterPhenotypes(
            contexts,
            traitId,
            window,
            phenotypeFilter$naAction,
            phenotypeFilter$outlierAction,
            phenotypeFilter$outlierPvalThreshold,
            phenotypeFilter$outlierArgs
        )
}

# One context's phenotype SummarizedExperiment as stored, for the internal
# callers that hold a single context name and read its traits or ranges.
# @noRd
.qtlContextSe <- function(x, context) {
    .qtlValidateContexts(x, context)
    .qtlPhenotypeList(x)[[context]]
}

# `contexts` is required and must all be known phenotype contexts.
# @noRd
.qtlValidateContexts <- function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        msg <- glue(
            "`context` is required for molecularTraits(QtlDataset). Pass a ",
            "character vector of one or more context names; use ",
            "contexts(x) to list the available contexts."
        )
        abort(msg)
    }
    available <- contexts(x)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        msg <- glue(
            "Unknown context(s): {str_flatten(bad, ', ')}. ",
            "Available: {str_flatten(available, ', ')}"
        )
        abort(msg)
    }
}

# Restrict each context's SE to `traitId`, warning about traits absent in a
# context.
# @noRd
.qtlFilterTraits <- function(out, contexts, traitId) {
    filtered <- map(
        seq_along(out),
        .qtlFilterTraitSe,
        out = out,
        traitId = traitId
    )
    set_names(filtered, contexts)
}

# Apply the optional trait / window / NA / outlier filters to the per-context
# phenotype SEs.
# @noRd
.qtlFilterPhenotypes <- function(
    out,
    contexts,
    traitId,
    window,
    naAction,
    outlierAction,
    outlierPvalThreshold,
    outlierArgs = list()
) {
    byTrait <- if (is.null(traitId)) {
        out
    } else {
        .qtlFilterTraits(out, contexts, traitId)
    }
    inRegion <- if (is.null(window)) {
        byTrait
    } else {
        set_names(
            map(
                byTrait,
                .qtlSeInRegion,
                region = .assertRegionGRanges(window, "`window`")
            ),
            contexts
        )
    }
    naHandled <- if (naAction == "keep") {
        inRegion
    } else {
        set_names(
            map(inRegion, .qtlApplyPhenoNaAction, naAction = naAction),
            contexts
        )
    }
    if (outlierAction == "keep") {
        return(naHandled)
    }
    set_names(
        map(
            naHandled,
            .qtlApplyPhenoOutliers,
            action = outlierAction,
            pvalThreshold = outlierPvalThreshold,
            outlierArgs = outlierArgs
        ),
        contexts
    )
}

# Internal: apply naAction to a SummarizedExperiment slice. SE assay rows
# are traits and columns are samples.
#   "drop"   -> drop samples (cols) where any selected trait is NA
#   "impute" -> mean-impute each trait (row) independently over its
#               non-NA sample values
# Operates jointly over the rows currently in `se` -- the caller is
# expected to have already subset the SE to the user's requested
# (traitId, region) subset.
.qtlApplyPhenoNaAction <- function(se, naAction) {
    assayName <- SummarizedExperiment::assayNames(se)[[1L]]
    Y <- SummarizedExperiment::assay(se, assayName)
    if (length(Y) == 0L) {
        return(se)
    }
    if (naAction == "drop") {
        return(se[, colSums(is.na(Y)) == 0L, drop = FALSE])
    }
    if (naAction != "impute" || !anyNA(Y)) {
        return(se)
    }
    SummarizedExperiment::`assay<-`(
        se,
        assayName,
        value = .qtlRowMeanImpute(Y)
    )
}

# Fill each row's missing values with that row's mean, or 0 when the row has
# nothing observed. One fill over the whole matrix rather than a copy per row.
# @noRd
.qtlRowMeanImpute <- function(Y) {
    naMask <- is.na(Y)
    # rowMeans of an all-NA row is NaN, which is the `else 0` case.
    means <- rowMeans(Y, na.rm = TRUE)
    filled <- if_else(is.nan(means), 0, means)
    replace(Y, naMask, filled[row(Y)[naMask]])
}

# Multivariate-outlier keep mask via Mahalanobis distance against a
# (preferably robust) centre / covariance estimate. Returns a logical
# vector of length nrow(Y); TRUE = keep, FALSE = drop.
#
# When the `robustbase` package is installed, the centre and covariance
# come from `robustbase::covMcd` (minimum-covariance-determinant) so the
# detector itself is resistant to the outliers it's trying to find.
# Without robustbase we fall back to `colMeans` / `cov` with a one-shot
# message; the test then still works but its estimates are pulled by
# the very outliers it should be flagging.
#
#' @title Options for the Robust Outlier Covariance Estimator
#' @description Build a checked record of extra arguments for
#'   \code{robustbase::covMcd()}, the minimum-covariance-determinant
#'   estimator behind \code{outlierAction = "drop"} on
#'   \code{\link{molecularTraits}} and
#'   \code{\link{residualizedPhenotypes}}.
#' @param ... Arguments for \code{robustbase::covMcd()} -- in practice
#'   \code{alpha} (the subset fraction the determinant is minimised over,
#'   which sets the breakdown point), \code{nsamp}, \code{use.correction}
#'   and \code{nmini} / \code{kmini}. \code{x} is the trait matrix pecotmr
#'   assembles and is refused, as is \code{seed}: pecotmr owns RNG
#'   reproducibility, and \code{covMcd()}'s \code{seed} restores the ambient
#'   \code{.Random.seed} on exit, which would silently undo the caller's
#'   stream.
#' @return A \code{MethodOptions} record for the \code{outlierArgs} argument.
#' @seealso \code{\link{molecularTraits}},
#'   \code{\link{residualizedPhenotypes}}
#' @examples
#' CovMcdOptions(alpha = 0.75)
#' @export
CovMcdOptions <- function(...) {
    extra <- list(...)
    .configRefuseOwned(
        extra,
        c(
            x = "the trait matrix pecotmr assembles",
            seed = "pecotmr's own RNG handling"
        ),
        "CovMcdOptions"
    )
    .newMethodOptions(
        "robustbase::covMcd",
        defaults = list(),
        extra = extra,
        label = "CovMcdOptions",
        engine = "covMcd"
    )
}

# The centre and covariance the Mahalanobis distance is measured against:
# robustbase's MCD when it is installed and converges, else the non-robust
# colMeans/cov pair. Both fallbacks are announced -- a silently non-robust
# estimate would make the outlier rule mean something different from what
# the caller asked for. `singularReported` records that robustbase has
# already warned its MCD covariance is singular, so the caller does not say
# it a second time.
# @noRd
#' @importFrom rlang try_fetch
.qtlOutlierCenterCov <- function(Y, outlierArgs = list()) {
    # Engine-call exception: robustbase::covMcd is a covariance estimator,
    # interchangeable with the colMeans/stats::cov fallback above -- a numeric
    # primitive for outlier detection while a QtlDataset is built.
    fallback <- list(
        center = colMeans(Y),
        cov = stats::cov(Y),
        singularReported = FALSE
    )
    if (!requireNamespace("robustbase", quietly = TRUE)) {
        msg <- glue(
            "outlier detection: install 'robustbase' for an MCD-based ",
            "estimator; falling back to non-robust colMeans/cov."
        )
        inform(msg)
        return(fallback)
    }
    mcd <- try_fetch(
        exec(robustbase::covMcd, Y, !!!as.list(outlierArgs)),
        error = function(cnd) {
            msg <- glue(
                "outlier detection: robustbase's MCD estimator failed; ",
                "falling back to non-robust colMeans/cov."
            )
            inform(msg, parent = cnd)
            NULL
        }
    )
    if (is.null(mcd)) {
        return(fallback)
    }
    list(
        center = mcd$center,
        cov = mcd$cov,
        singularReported = !is.null(mcd$singularity)
    )
}

# Significance: per-sample chi-squared(p) p-value with Bonferroni
# correction over the sample count. A sample is flagged when its
# corrected p-value falls below `pvalThreshold`. With single-trait Y
# (ncol == 1) this reduces to the standard z-test on (y - center)/sd.
#
# Returns all-TRUE (no-op) when there are too few samples to support
# a covariance estimate (n < p + 2).
#' @importFrom rlang try_fetch
.qtlOutlierKeepMask <- function(Y, pvalThreshold, outlierArgs = list()) {
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    if (n == 0L || p == 0L) {
        return(rep(TRUE, n))
    }
    if (n < p + 2L) {
        msg <- glue(
            "outlier detection skipped: {n} samples < {p} traits + 2 ",
            "needed for a covariance estimate."
        )
        warn(msg)
        return(rep(TRUE, n))
    }
    est <- .qtlOutlierCenterCov(Y, outlierArgs)
    ctr <- est$center
    covMat <- est$cov
    invCov <- try_fetch(
        solve(covMat),
        error = function(cnd) {
            # When robustbase has already warned that its covariance is
            # singular, that warning is the one signal the caller gets.
            if (!est$singularReported) {
                msg <- glue(
                    "outlier detection: the trait covariance is singular; ",
                    "using a Moore-Penrose pseudo-inverse instead."
                )
                inform(msg, parent = cnd)
            }
            MASS::ginv(covMat)
        }
    )
    Yc <- sweep(Y, 2L, ctr)
    d2 <- rowSums((Yc %*% invCov) * Yc)
    raw <- stats::pchisq(d2, df = p, lower.tail = FALSE)
    raw >= (pvalThreshold / n)
}

# Wrapper: apply the keep-mask to a SummarizedExperiment slice. SE
# columns are samples; transpose the assay (traits x samples) before
# calling .qtlOutlierKeepMask which expects samples x traits.
.qtlApplyPhenoOutliers <- function(
    se,
    action,
    pvalThreshold,
    outlierArgs = list()
) {
    if (action == "keep") {
        return(se)
    }
    assayName <- SummarizedExperiment::assayNames(se)[[1L]]
    Y <- t(SummarizedExperiment::assay(se, assayName))
    keep <- .qtlOutlierKeepMask(Y, pvalThreshold, outlierArgs)
    if (all(keep)) {
        return(se)
    }
    se[, keep, drop = FALSE]
}

#' @rdname phenotypeCovariates
#' @export
setMethod("phenotypeCovariates", "QtlDataset", function(x, context) {
    if (missing(context) || is.null(context) || length(context) == 0L) {
        abort("`context` is required.")
    }
    available <- contexts(x)
    bad <- setdiff(context, available)
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    set_names(map(context, .qtlContextColData, x = x), context)
})

# Whether to scale residual columns to unit variance: what the bundle sets,
# else TRUE -- the run's setting, since a dataset records none. A NULL bundle
# answers NULL for the field, so it resolves to TRUE as well.
#
# Unset (NULL) rather than FALSE means "not asked", the same three-state
# convention `residualizePhenotype` / `residualizeGenotype` use -- a flat
# formal defaulting to TRUE could not tell "FALSE" from "not given".
# @noRd
.qtlScaleResiduals <- function(residualizationParam) {
    residualizationParam$scaleResiduals %||% TRUE
}

# Internal: residualize a numeric matrix Y (n x k) against a covariate
# matrix C (n x p) via pivoted QR decomposition. Adds an intercept column
# to C. When C is rank-deficient (e.g., union of all contexts' phenotype
# covariates includes collinear / duplicate columns), the pivoted QR drops
# the redundant columns automatically. Optionally rescales each residual
# column to unit standard deviation; constant-valued columns are left
# unchanged.
.qtlResidualizeQr <- function(Y, C, scaleResiduals = TRUE) {
    # A trait filter that keeps none of a context's traits leaves nothing to
    # residualize; lm.fit() cannot take a zero-column response.
    if (ncol(Y) == 0L) {
        return(Y)
    }
    X <- if (is.null(C) || ncol(C) == 0L) {
        matrix(
            1,
            nrow = nrow(Y),
            ncol = 1L,
            dimnames = list(rownames(Y), "intercept")
        )
    } else {
        cbind(intercept = 1, C)
    }
    # `qr.resid` does not support LAPACK pivoted QR, so use `lm.fit`. It
    # handles rank-deficient designs gracefully via base-R's pivoted QR
    # internally -- same effect the LAPACK path was meant to deliver.
    res <- `dimnames<-`(
        as.matrix(stats::lm.fit(x = X, y = Y)$residuals),
        list(rownames(Y), colnames(Y))
    )
    if (isTRUE(scaleResiduals)) {
        # `sds == 0` exact-zero test is unreliable for residuals coming out of
        # lm.fit on a constant Y: roundoff gives sd ~ 1e-16 instead of 0, and
        # dividing the (also-tiny) residuals by it amplifies floating-point
        # noise to unit-scale. Treat anything below sqrt(.Machine$double.eps)
        # as effectively zero (column is constant) and skip rescaling.
        rawSds <- apply(res, 2L, stats::sd, na.rm = TRUE)
        nearZero <- !is.finite(rawSds) | rawSds < sqrt(.Machine$double.eps)
        # One multiplier per column: 1/sd where it is meaningful, 0 where the
        # column is constant, which zeroes that column outright.
        return(sweep(res, 2L, if_else(nearZero, 0, 1 / rawSds), FUN = "*"))
    }
    res
}

# Resolve the phenotype-covariate selection for one context: NULL requested ->
# all available covariates; otherwise validate the requested names are present.
# @noRd
.qtlResolveOne <- function(ctx, requested, x) {
    se <- .qtlContextSe(x, ctx)
    avail <- colnames(SummarizedExperiment::colData(se))
    if (is.null(requested)) {
        return(avail)
    }
    keep <- intersect(requested, avail)
    if (length(keep) != length(requested)) {
        missingNames <- setdiff(requested, avail)
        msg <- glue(
            "residualizationParam$phenotypeCovariates: context '{ctx}' has no ",
            "covariate(s) named: {str_flatten(missingNames, ', ')}"
        )
        abort(msg)
    }
    keep
}

# Internal: validate and resolve a covariate-selection field of
# ResidualizationParam() against a
# set of contexts and the covariates actually present in those contexts'
# colData. Accepts either NULL (use all), a character vector (apply to all
# listed contexts), or a named list keyed by context. Returns a named list
# keyed by context giving the actual character vector of covariate names
# to use for that context (or character(0) if none). Errors when:
#   - a named-list key is not in `contexts`
#   - `contexts` contains entries missing from a supplied named-list
#     (per the rule: named-list keys must equal `contexts`)
#   - an explicitly requested name matches no actual covariate
.qtlResolvePhenoSelection <- function(x, contexts, toResidualize) {
    if (is.null(toResidualize)) {
        return(set_names(
            map(contexts, .qtlResolveOne, requested = NULL, x = x),
            contexts
        ))
    }
    if (is.list(toResidualize)) {
        return(.qtlPhenoSelectionList(x, contexts, toResidualize))
    }
    if (is.character(toResidualize)) {
        return(set_names(
            map(contexts, .qtlResolveOne, requested = toResidualize, x = x),
            contexts
        ))
    }
    msg <- glue(
        "`phenotypeCovariates` must be NULL, a character vector, ",
        "or a named list keyed by context."
    )
    abort(msg)
}

# List-form `phenotypeCovariates`: must be named with EXACTLY the
# `contexts` set. Resolves each context's selection.
# @noRd
.qtlPhenoSelectionList <- function(x, contexts, toResidualize) {
    toResNames <- names(toResidualize)
    if (
        is.null(toResNames) || any(str_length(toResNames) == 0L, na.rm = TRUE)
    ) {
        msg <- glue(
            "`phenotypeCovariates`: when supplied as a list, it ",
            "must be named with context names."
        )
        abort(msg)
    }
    badKeys <- setdiff(names(toResidualize), contexts)
    if (length(badKeys) > 0L) {
        msg <- glue(
            "`phenotypeCovariates`: list key(s) not in ",
            "`context`: {str_flatten(badKeys, ', ')}"
        )
        abort(msg)
    }
    missingKeys <- setdiff(contexts, names(toResidualize))
    if (length(missingKeys) > 0L) {
        msg <- glue(
            "`phenotypeCovariates`: list does not cover all ",
            "`context`. Per-context lists must have exactly the same context ",
            "set as `context`. Missing keys: ",
            "{str_flatten(missingKeys, ', ')}"
        )
        abort(msg)
    }
    set_names(
        map(contexts, .qtlResolveContext, toResidualize = toResidualize, x = x),
        contexts
    )
}

# Internal: validate the genotype-covariate selection vector. Returns
# character(0) when nothing selected, the resolved set otherwise.
.qtlResolveGenoSelection <- function(x, toResidualize) {
    avail <- colnames(genotypeCovariates(x)) %||% character(0)
    if (is.null(toResidualize)) {
        return(avail)
    }
    keep <- intersect(toResidualize, avail)
    if (length(keep) != length(toResidualize)) {
        missingNames <- setdiff(toResidualize, avail)
        msg <- glue(
            "`genotypeCovariates`: no covariate(s) named: ",
            "{str_flatten(missingNames, ', ')}"
        )
        abort(msg)
    }
    keep
}

# Internal: build the covariate matrix used for residualization, given a
# set of contexts, the resolved per-context phenotype selections, and the
# resolved genotype-covariate selection. Honors the inclusion flags. For
# `length(contexts) == 1` (per-context mode), the per-context phenotype
# covariates and the genotype covariates are taken with no cross-context
# alignment. For `length(contexts) >= 2` (joint mode), per-context
# phenotype covariates from all listed contexts are concatenated
# (prefixed with "{context}." to keep same-named columns distinct) and
# the sample set is intersected across all contributing matrices.
# Returns a single matrix (with rownames = sample IDs) or NULL.
.qtlBuildResidualizationDesign <- function(
    x,
    contexts,
    phenoSelection,
    genoSelection,
    includePheno,
    includeGeno
) {
    perContext <- if (includePheno) {
        .qtlPhenoCovBlocks(x, contexts, phenoSelection)
    } else {
        list()
    }
    gCov <- if (includeGeno && length(genoSelection) > 0L) {
        genotypeCovariates(x)[, genoSelection, drop = FALSE]
    } else {
        matrix(numeric(0), nrow = 0, ncol = 0)
    }
    haveAny <- length(perContext) > 0L || (!is.null(gCov) && ncol(gCov) > 0L)
    if (!haveAny) {
        return(NULL)
    }
    design <- .qtlAlignCovariates(perContext, gCov)
    # NULL here means covariates were asked for but no sample carries them
    # all. Returning it would residualize against nothing and hand back the
    # raw values, so the caller hears about it instead.
    if (is.null(design)) {
        abort(glue(
            "No samples in common among the covariate blocks requested ",
            "for contexts: {str_flatten(contexts, ', ')}"
        ))
    }
    design
}

# Per-context phenotype covariate matrices (colData columns from
# phenoSelection), column names prefixed with the context.
# @noRd
.qtlPhenoCovBlocks <- function(x, contexts, phenoSelection) {
    compact(set_names(
        map(
            contexts,
            .qtlPhenoCovBlock,
            x = x,
            phenoSelection = phenoSelection
        ),
        contexts
    ))
}

# One context's selected covariate columns, context-qualified, or NULL when
# nothing is selected there.
#
# Sample names survive the coercion because the constructor requires every
# context to name its samples, which SummarizedExperiment carries into
# rownames(colData(se)).
# @noRd
.qtlPhenoCovBlock <- function(ctx, x, phenoSelection) {
    keep <- phenoSelection[[ctx]]
    if (length(keep) == 0L) {
        return(NULL)
    }
    se <- .qtlContextSe(x, ctx)
    cd <- as.matrix(as.data.frame(SummarizedExperiment::colData(se)))
    block <- cd[, keep, drop = FALSE]
    `colnames<-`(block, str_c(ctx, ".", colnames(block)))
}

# Intersect the covariate blocks to their common samples and column-bind them
# into one design matrix; NULL when no samples are shared.
# @noRd
.qtlAlignCovariates <- function(perContext, gCov) {
    sampleSets <- c(
        map(perContext, .qtlBlockSamples),
        if (!is.null(gCov) && ncol(gCov) > 0L) {
            list(.qtlBlockSamples(gCov))
        }
    )
    common <- if (length(sampleSets) == 0L) {
        character(0)
    } else {
        reduce(sampleSets, intersect)
    }
    if (length(common) == 0L) {
        return(NULL)
    }
    genoBlock <- if (!is.null(gCov) && ncol(gCov) > 0L) {
        list(gCov[common, , drop = FALSE])
    } else {
        list()
    }
    blocks <- c(
        map(perContext, .qtlRestrictToSamples, common = common),
        genoBlock
    )
    exec(cbind, !!!blocks)
}

# @noRd
.qtlRestrictToSamples <- function(mat, common) {
    mat[common, , drop = FALSE]
}

# The samples a covariate block covers, or NULL when it does not say. A block
# with no rows covers none of them, which base R reports as NULL rownames
# rather than an empty character vector; reading that as "unconstrained"
# would let a covariate matrix sharing no samples with the panel fall through
# to an out-of-bounds subscript instead of resolving to no design.
# @noRd
.qtlBlockSamples <- function(mat) {
    if (nrow(mat) == 0L) {
        return(character(0))
    }
    rownames(mat)
}

# Internal: resolve missing values in the residualization covariate matrix
# `C` (samples x covariates) before it reaches `stats::lm.fit`, which does
# not tolerate NA in the design and would otherwise error. Two strategies:
#   - "impute" (default): replace each covariate column's NA cells with that
#     column's observed (non-NA) mean. A wholly-missing column has an
#     undefined mean and is filled with 0, so it contributes nothing to the
#     fit rather than poisoning every row.
#   - "drop": complete-case -- remove any sample (row) carrying an NA in any
#     covariate. The caller's downstream sample intersection then narrows the
#     response matrix (G or Y) to the retained samples.
# Returns the cleaned matrix (fewer rows possible under "drop"), or `C`
# unchanged when it is NULL or already NA-free.
.qtlHandleCovariateNa <- function(C, action = c("impute", "drop")) {
    action <- arg_match(action)
    if (is.null(C) || !anyNA(C)) {
        return(C)
    }
    if (action == "drop") {
        keep <- rowSums(is.na(C)) == 0L
        return(C[keep, , drop = FALSE])
    }
    naMask <- is.na(C)
    # colMeans of an all-NA column is NaN, which is the non-finite `else 0`.
    means <- colMeans(C, na.rm = TRUE)
    filled <- if_else(is.finite(means), means, 0)
    replace(C, naMask, filled[col(C)[naMask]])
}

# Intersect a genotype matrix G and covariate design C to their common samples
# (no-op when C is NULL); errors if they share none.
# @noRd
.qtlAlignGC <- function(G, C, contexts) {
    if (is.null(C)) {
        return(list(G = G, C = C))
    }
    common <- intersect(rownames(G), rownames(C))
    if (length(common) == 0L) {
        msg <- glue(
            "No samples in common between the genotype matrix and the ",
            "covariate matrix for contexts: ",
            "{str_flatten(contexts, ', ')}"
        )
        abort(msg)
    }
    list(G = G[common, , drop = FALSE], C = C[common, , drop = FALSE])
}

# residualizedGenotypes() for the pipelines: the genotypes of `window` and
# `samples` under the run's resolved genotype `filter`, residualized on the
# covariates of `contexts`. The accessor passes the whole dataset and the
# all-default filter.
# @noRd
.qtlResidualizedGenotypes <- function(
    x,
    contexts,
    window = NULL,
    samples = NULL,
    filter = .qtlResolveFilter(),
    residualizationParam = ResidualizationParam()
) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        msg <- glue(
            "`context` is required for ",
            "residualizedGenotypes(QtlDataset). ",
            "Use contexts(x) to list the available contexts. ",
            "Pass a single context for per-context mode or multiple ",
            "contexts for joint mode (sample intersection)."
        )
        abort(msg)
    }
    .assertMethodParam(
        residualizationParam,
        "ResidualizationParam"
    )
    bad <- setdiff(contexts, contexts(x))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, \', \')}")
        abort(msg)
    }
    G <- .qtlExtractBlock(x, window, samples, filter)$geno
    if (ncol(G) == 0L) {
        return(G)
    }
    C <- .qtlResidCovariateMatrix(x, contexts, residualizationParam)
    aligned <- .qtlAlignGC(G, C, contexts)
    .qtlResidualizeQr(
        aligned$G,
        aligned$C,
        scaleResiduals = .qtlScaleResiduals(residualizationParam)
    )
}

#' @rdname residualizedGenotypes
#' @export
setMethod(
    "residualizedGenotypes",
    "QtlDataset",
    function(x, context, residualizationParam = ResidualizationParam()) {
        .qtlResidualizedGenotypes(
            x,
            context,
            residualizationParam = residualizationParam
        )
    }
)

# NA-handled raw phenotypes, a list keyed by context.
# @noRd
.qtlResidPhenoY <- function(x, contexts, traitId, window, naAction) {
    # Missing values only: outliers are judged on the residualized values.
    .qtlMolecularTraits(
        x,
        contexts = contexts,
        traitId = traitId,
        window = window,
        phenotypeFilter = list_modify(
            .qtlKeepAllPhenotypes(),
            naAction = naAction
        )
    )
}

# Phenotype and covariate matrices restricted to the samples they share.
# A NULL covariate matrix leaves the phenotypes whole.
# @noRd
.qtlAlignPhenoCovariates <- function(allY, C, ctx) {
    if (is.null(C)) {
        return(list(Y = allY, C = NULL))
    }
    common <- intersect(rownames(allY), rownames(C))
    if (length(common) == 0L) {
        abort(glue(
            "context '{ctx}': no samples shared between phenotype data ",
            "and the resolved covariate matrix."
        ))
    }
    list(
        Y = allY[common, , drop = FALSE],
        C = C[common, , drop = FALSE]
    )
}

# Residualize one context's phenotypes against the covariate design
# (intersected to common samples) and drop residual-scale outliers.
# @noRd
.qtlResidualizeContextPheno <- function(
    se,
    C,
    ctx,
    outlierAction,
    outlierPvalThreshold,
    scaleResiduals,
    outlierArgs = list()
) {
    allY <- t(SummarizedExperiment::assay(se)) # samples x traits
    aligned <- .qtlAlignPhenoCovariates(allY, C, ctx)
    allRes <- .qtlResidualizeQr(
        aligned$Y,
        aligned$C,
        scaleResiduals = scaleResiduals
    )
    if (outlierAction == "keep") {
        return(allRes)
    }
    # Outlier detection on the residualized scale.
    keep <- .qtlOutlierKeepMask(allRes, outlierPvalThreshold, outlierArgs)
    if (all(keep)) {
        return(allRes)
    }
    allRes[keep, , drop = FALSE]
}

# residualizedPhenotypes() worker: resolve covariate inclusion + selections,
# NA-handle Y, build the covariate design, and per-context residualize +
# outlier-filter, under the run's resolved phenotype QC. The pipelines may
# also keep only the traits inside an analysis `window`; the accessor passes
# no window and no QC. Always a list keyed by context.
# @noRd
.qtlResidualizedPhenotypes <- function(
    x,
    contexts,
    traitId = NULL,
    window = NULL,
    phenotypeFilter = .qtlKeepAllPhenotypes(),
    residualizationParam = ResidualizationParam()
) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        abort("`context` is required for residualizedPhenotypes().")
    }
    .assertMethodParam(
        residualizationParam,
        "ResidualizationParam"
    )
    bad <- setdiff(contexts, contexts(x))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, \', \')}")
        abort(msg)
    }
    C <- .qtlResidCovariateMatrix(x, contexts, residualizationParam)
    Yraw <- .qtlResidPhenoY(
        x,
        contexts,
        traitId,
        window,
        phenotypeFilter$naAction
    )
    .qtlResidPerContext(
        contexts,
        Yraw = Yraw,
        C = C,
        x = x,
        outlierAction = phenotypeFilter$outlierAction,
        outlierPvalThreshold = phenotypeFilter$outlierPvalThreshold,
        outlierArgs = phenotypeFilter$outlierArgs,
        scaleResiduals = .qtlScaleResiduals(residualizationParam)
    )
}

#' @rdname residualizedPhenotypes
#' @export
setMethod(
    "residualizedPhenotypes",
    "QtlDataset",
    function(
        x,
        context,
        trait = NULL,
        residualizationParam = ResidualizationParam()
    ) {
        .qtlResidualizedPhenotypes(
            x,
            context,
            trait,
            residualizationParam = residualizationParam
        )
    }
)

# What the phenotypes are regressed against: resolve the convenience vs
# precise inclusion flags, pick the phenotype- and genotype-covariate sets
# they select, build the design, and apply the NA policy to it. One step,
# because none of the four halves means anything without the others.
# @noRd
.qtlResidCovariateMatrix <- function(x, contexts, residualizationParam) {
    # arg_match() reports on its argument by name, so it must be handed a
    # symbol: resolve the bundle field first, then validate it.
    covariateNaAction <- residualizationParam$covariateNaAction %||% "impute"
    covariateNaAction <- arg_match(covariateNaAction, c("impute", "drop"))
    phenoSel <- .qtlResolvePhenoSelection(
        x,
        contexts,
        residualizationParam$phenotypeCovariates
    )
    design <- .qtlBuildResidualizationDesign(
        x,
        contexts = contexts,
        phenoSelection = phenoSel,
        genoSelection = .qtlResolveGenoSelection(
            x,
            residualizationParam$genotypeCovariates
        ),
        # Unset means "residualize on this side": the accessors' own
        # historical default, now expressed as the absence of a value
        # rather than as a formal defaulting to TRUE.
        includePheno = residualizationParam$residualizePhenotype %||% TRUE,
        includeGeno = residualizationParam$residualizeGenotype %||% TRUE
    )
    .qtlHandleCovariateNa(design, covariateNaAction)
}

# Residualize every requested context, unwrapping a single-context request
# the way the rest of the QtlDataset accessors do.
# @noRd
.qtlResidPerContext <- function(
    contexts,
    Yraw,
    C,
    x,
    outlierAction,
    outlierPvalThreshold,
    outlierArgs,
    scaleResiduals
) {
    set_names(
        map(
            contexts,
            .qtlResidualizeContext,
            Yraw = Yraw,
            C = C,
            x = x,
            outlierAction = outlierAction,
            outlierPvalThreshold = outlierPvalThreshold,
            outlierArgs = outlierArgs,
            scaleResiduals = scaleResiduals
        ),
        contexts
    )
}

#' @rdname show-methods
#' @export
setMethod("show", "QtlDataset", function(object) {
    pheno <- .qtlPhenotypeList(object)
    nCtx <- length(pheno)
    ctxNames <- names(pheno)
    totalTraits <- length(unique(unname(list_c(map(pheno, rownames)))))
    cat(glue(
        "QtlDataset for study '{object@studyName}'\n",
        .trim = FALSE
    ))
    cat(glue(
        "  {nCtx} context(s): {str_flatten(ctxNames, ', ')}\n",
        .trim = FALSE
    ))
    cat(glue("  {totalTraits} unique traits across contexts\n", .trim = FALSE))
    cat(glue(
        "  Genotypes: {.ldRefLabel(.qtlGenotypeSe(object))}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Genotype covariates: ",
        "{ncol(genotypeCovariates(object))} cols\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Samples: {nrow(MultiAssayExperiment::colData(object))}\n",
        .trim = FALSE
    ))
})

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Restrict context `i`'s SE to the requested traits (warns about absent ones).
# @noRd
.qtlFilterTraitSe <- function(i, out, traitId) {
    se <- out[[i]]
    ctx <- names(out)[[i]]
    present <- intersect(traitId, rownames(se))
    missing <- setdiff(traitId, rownames(se))
    if (length(missing) > 0L) {
        msg <- glue(
            "context '{ctx}' is missing trait(s): ",
            "{str_flatten(missing, ', ')}"
        )
        warn(msg)
    }
    se[present, , drop = FALSE]
}

# One context's SE restricted to features overlapping `region`. A context with
# no trait on the window's chromosomes keeps none, decided before
# overlapsAny(), which would warn that the two share no seqlevel.
# @noRd
.qtlSeInRegion <- function(se, region) {
    rr <- SummarizedExperiment::rowRanges(se)
    if (
        !any(is_in(as.character(seqnames(rr)), as.character(seqnames(region))))
    ) {
        return(se[integer(0), , drop = FALSE])
    }
    se[IRanges::overlapsAny(rr, region), , drop = FALSE]
}

# One context's covariate matrix (colData of its phenotype SE).
# @noRd
.qtlContextColData <- function(ctx, x) {
    se <- .qtlContextSe(x, ctx)
    cd <- SummarizedExperiment::colData(se)
    as.matrix(as.data.frame(cd))
}

# Resolve the per-context residualization spec for context `ctx`.
# @noRd
.qtlResolveContext <- function(ctx, toResidualize, x) {
    .qtlResolveOne(ctx, toResidualize[[ctx]], x)
}

# Residualize + filter one context's raw phenotype matrix against covariates C.
# @noRd
.qtlResidualizeContext <- function(
    ctx,
    Yraw,
    C,
    x,
    outlierAction,
    outlierPvalThreshold,
    outlierArgs = list(),
    scaleResiduals = TRUE
) {
    .qtlResidualizeContextPheno(
        Yraw[[ctx]],
        C,
        ctx,
        outlierAction,
        outlierPvalThreshold,
        scaleResiduals,
        outlierArgs
    )
}

# Which covariates to residualize on: a character vector of names applied to
# every context, a list keyed by context name, or NULL to take each context's
# own. .qtlResolvePhenoSelection() branches on exactly these three.
setClassUnion("CovariateSelection", c("character", "list", "NULL"))

#' @rdname ResidualizationParam
#' @aliases ResidualizationParam-class
#' @exportClass ResidualizationParam
setClass(
    "ResidualizationParam",
    contains = "MethodParam",
    slots = c(
        phenotypeCovariates = "CovariateSelection",
        genotypeCovariates = "CovariateSelection",
        residualizePhenotype = "logical_OR_NULL",
        residualizeGenotype = "logical_OR_NULL",
        scaleResiduals = "logical_OR_NULL",
        covariateNaAction = "character_OR_NULL"
    )
)

#' @title Covariate Residualization Settings
#' @description What is regressed out of the phenotype and genotype before
#'   fitting. Shared by \code{\link{fineMappingPipeline}} and
#'   \code{\link{twasWeightsPipeline}}, which hold the same four settings.
#'
#'   \code{fineMappingPipeline}'s \code{usePCA} / \code{nPCs} are
#'   \strong{not} here. They do not residualize anything: they PCA-reduce a
#'   multi-trait context's phenotype matrix and fine-map each top principal
#'   component \emph{as a trait}. That is an analysis mode, a sibling of the
#'   univariate and multivariate dispatch paths, not a covariate setting.
#' @section Summary-statistics inputs:
#'   A \code{QtlSumStats} / \code{GwasSumStats} input carries no genotypes
#'   or covariates, so nothing here applies and the whole bundle is
#'   \strong{ignored} --- not refused.
#'
#'   That is the opposite of \code{\link{CrossValidationParam}}, which
#'   \emph{is} refused on those inputs
#'   (\code{fineMappingPipeline(QtlSumStats, crossValidation = ...)} errors).
#'   The difference is the default: CV is off unless asked for
#'   (\code{folds = 0}), so a non-default value on a summary-statistics run
#'   is an explicit request for something impossible. Residualization is on
#'   by default (\code{residualizePhenotype} and \code{residualizeGenotype}
#'   are both \code{TRUE}), so refusing a non-default would reject the
#'   \emph{default} bundle and force every sumstats caller to unset it.
#'
#'   One bundle therefore travels to either input kind unchanged, the same
#'   way \code{CrossValidationParam}'s \code{weightMethods} /
#'   \code{maxVariants} are carried but ignored by
#'   \code{fineMappingPipeline}.
#' @param phenotypeCovariates Covariates to residualize the phenotype on, or
#'   \code{NULL} (default) for the dataset's own.
#' @param genotypeCovariates Covariates to residualize the genotype on, or
#'   \code{NULL} (default) for the dataset's own.
#' @param residualizePhenotype Logical. Residualize on the phenotype
#'   covariates. \code{NULL} (default) means unset, which the accessors
#'   read as \code{TRUE}.
#' @param residualizeGenotype Logical. Residualize on the genotype
#'   covariates. \code{NULL} (default) means unset, read as \code{TRUE}.
#' @param scaleResiduals Logical. Scale each residual column to unit
#'   variance, which is what puts a fitted effect on the standardized
#'   genotype scale. \code{NULL} (default) means unset, read as
#'   \code{TRUE}. It is the only scale knob the pipelines need: it reaches
#'   \code{\link{residualizedGenotypes}}, the one point every
#'   individual-level fit draws its genotypes through.
#' @param covariateNaAction How missing covariate values are handled when the
#'   design is assembled: \code{"impute"} or \code{"drop"}. \code{NULL}
#'   (default) leaves the accessor's own default in place.
#' @return A \code{ResidualizationParam} object, a \code{\link{MethodParam}}.
#' @examples
#' ResidualizationParam(residualizeGenotype = FALSE)
#' @export
ResidualizationParam <- function(
    phenotypeCovariates = NULL,
    genotypeCovariates = NULL,
    residualizePhenotype = NULL,
    residualizeGenotype = NULL,
    scaleResiduals = NULL,
    covariateNaAction = NULL
) {
    # NULL means "unset", which is why the accessors no longer need a
    # `missing()` companion for each flag: FALSE and "not given" are both
    # falsy, so a flat formal defaulting to TRUE could not tell them apart,
    # and each needed a second formal recording whether it was supplied.
    new(
        "ResidualizationParam",
        phenotypeCovariates = phenotypeCovariates,
        genotypeCovariates = genotypeCovariates,
        residualizePhenotype = residualizePhenotype,
        residualizeGenotype = residualizeGenotype,
        scaleResiduals = scaleResiduals,
        covariateNaAction = covariateNaAction
    )
}

# Slots a QtlDataset held before its run settings moved to the pipelines'
# Params (GenotypeFilterParam, ResidualizationParam). An object saved then
# still carries them as attributes, which updateObject() drops.
# @noRd
.QTL_RETIRED_SLOTS <- c(
    "scaleResiduals",
    "mafCutoff",
    "macCutoff",
    "xvarCutoff",
    "imissCutoff",
    "keepVariants",
    "keepIndel"
)

# A QtlDataset saved before its study's slot was named `studyName` holds it
# as `study`; one saved with run settings has them dropped
# (.QTL_RETIRED_SLOTS); and its genotype panel is brought up to date
# (.qtlUpdateGenotypes()).
#' @rdname updateObject-methods
#' @export
setMethod("updateObject", "QtlDataset", function(object, ..., verbose = FALSE) {
    if (isTRUE(verbose)) {
        message("updateObject(object = 'QtlDataset')")
    }
    trimmed <- reduce(
        .QTL_RETIRED_SLOTS,
        .qtlDropAttribute,
        .init = .renameSavedSlots(object, c(study = "studyName"))
    )
    # Updated here, not as callNextMethod()'s argument: the inherited update
    # first reads its argument inside try(), so an error there would be
    # swallowed and the update run a second time.
    updated <- .qtlUpdateGenotypes(trimmed)
    methods::callNextMethod(updated, ..., verbose = verbose)
})

# @noRd
.qtlDropAttribute <- function(object, nm) {
    attr(object, nm) <- NULL
    object
}

# A dataset whose genotype panel was saved before panels took today's shape
# has it read again (.updatePanel()). Unlike an LD reference it cannot be
# kept as saved when that fails -- the dataset is nothing without its
# genotypes -- so the error stands.
# @noRd
.qtlUpdateGenotypes <- function(x) {
    panel <- .qtlGenotypeSe(x)
    updated <- .updatePanel(panel)
    if (identical(updated, panel)) {
        return(x)
    }
    .qtlWithGenotypeHandle(x, updated)
}

#' @rdname GenotypeFilterParam
#' @aliases GenotypeFilterParam-class
#' @exportClass GenotypeFilterParam
setClass(
    "GenotypeFilterParam",
    contains = "MethodParam",
    slots = c(
        mafCutoff = "numeric_OR_NULL",
        macCutoff = "numeric_OR_NULL",
        xvarCutoff = "numeric_OR_NULL",
        imissCutoff = "numeric_OR_NULL",
        keepSamples = "character_OR_NULL",
        keepVariants = "character_OR_NULL",
        keepIndel = "logical_OR_NULL"
    )
)

#' @title Genotype Filtering Options
#' @description Which variants and samples survive when a genotype matrix is
#'   assembled, as one checked bundle: the genotype QC of a run of
#'   \code{\link{fineMappingPipeline}}, \code{\link{twasWeightsPipeline}}
#'   or \code{\link{colocboostPipeline}} on a \code{\link{QtlDataset}}.
#'   A \code{QtlDataset} records no filter of its own; \code{keepSamples}
#'   narrows it once at the start of the run, and the other fields apply to
#'   each extracted block, against that block's samples.
#' @section Unset fields:
#'   Every field defaults to \code{NULL}, meaning \strong{not set}: an unset
#'   field is absent from the result rather than carried as \code{NULL}, and
#'   takes its default when the run resolves the bundle --- \code{0} for
#'   each cutoff, no restriction for each \code{keep}, \code{TRUE} for
#'   \code{keepIndel}. So \code{GenotypeFilterParam()} filters nothing.
#' @param mafCutoff Minor-allele-frequency floor.
#' @param macCutoff Minor-allele-count floor; the stricter of this and
#'   \code{mafCutoff} applies.
#' @param xvarCutoff Genotype-variance floor.
#' @param imissCutoff Per-sample missingness ceiling: a sample missing more
#'   than this fraction of an extracted block's genotypes is dropped from
#'   it.
#' @param keepSamples Sample ids to restrict to; \code{character(0)} keeps
#'   all samples.
#' @param keepVariants Variant ids to restrict to; \code{character(0)} keeps
#'   all variants.
#' @param keepIndel Logical. Retain insertions and deletions.
#' @return A \code{GenotypeFilterParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{PanelFilterParam}}, \code{\link{SumstatsFilterParam}}
#' @examples
#' GenotypeFilterParam(mafCutoff = 0.01, keepIndel = FALSE)
#' @export
GenotypeFilterParam <- function(
    mafCutoff = NULL,
    macCutoff = NULL,
    xvarCutoff = NULL,
    imissCutoff = NULL,
    keepSamples = NULL,
    keepVariants = NULL,
    keepIndel = NULL
) {
    new(
        "GenotypeFilterParam",
        mafCutoff = mafCutoff,
        macCutoff = macCutoff,
        xvarCutoff = xvarCutoff,
        imissCutoff = imissCutoff,
        keepSamples = keepSamples,
        keepVariants = keepVariants,
        keepIndel = keepIndel
    )
}

# Each cutoff a single finite non-negative number and keepIndel a flag,
# checked where the filter is built, so a bad value is reported against the
# argument that set it rather than surfacing inside an extraction. An unset
# field has nothing to check.
setValidity("GenotypeFilterParam", function(object) {
    coll <- makeAssertCollection()
    walk(
        c("mafCutoff", "macCutoff", "xvarCutoff", "imissCutoff"),
        .gfpCheckCutoff,
        object = object,
        coll = coll
    )
    if (!is.null(object@keepIndel)) {
        assertFlag(object@keepIndel, .var.name = "keepIndel", add = coll)
    }
    problems <- coll$getMessages()
    if (length(problems) == 0L) TRUE else problems
})

# @noRd
.gfpCheckCutoff <- function(nm, object, coll) {
    value <- methods::slot(object, nm)
    if (is.null(value)) {
        return(invisible(NULL))
    }
    assertNumber(value, lower = 0, finite = TRUE, .var.name = nm, add = coll)
}

#' @rdname PhenotypeFilterParam
#' @aliases PhenotypeFilterParam-class
#' @exportClass PhenotypeFilterParam
setClass(
    "PhenotypeFilterParam",
    contains = "MethodParam",
    slots = c(
        naAction = "character_OR_NULL",
        outlierAction = "character_OR_NULL",
        outlierPvalThreshold = "numeric_OR_NULL",
        outlierArgs = "MethodOptions_OR_NULL"
    )
)

#' @title Phenotype Quality-Control Options
#' @description How a run of \code{\link{fineMappingPipeline}},
#'   \code{\link{twasWeightsPipeline}} or \code{\link{colocboostPipeline}} on
#'   a \code{\link{QtlDataset}} treats the phenotypes it reads: missing
#'   values, and multivariate outlier samples. The phenotype counterpart of
#'   \code{\link{GenotypeFilterParam}}. A \code{QtlDataset}'s own accessors
#'   (\code{\link{molecularTraits}}, \code{\link{residualizedPhenotypes}})
#'   apply none of it and return the data as stored.
#' @section Unset fields:
#'   Every field defaults to \code{NULL}, meaning \strong{not set}, and takes
#'   its default when the run resolves the bundle. \code{naAction} resolves
#'   to the pipeline's own policy: fine-mapping and TWAS fits cannot take a
#'   missing phenotype, so they drop the affected samples (and refuse
#'   \code{"keep"}), while colocboost receives the phenotypes as they are.
#'   Outliers are kept unless \code{outlierAction = "drop"}.
#' @param naAction How missing phenotype values are handled: \code{"drop"}
#'   the samples with any missing value among the selected traits,
#'   \code{"impute"} each trait's mean, or \code{"keep"} them as they are.
#' @param outlierAction \code{"keep"} or \code{"drop"} samples flagged as
#'   multivariate outliers by a robust (MCD) Mahalanobis distance, computed
#'   on the residualized phenotypes.
#' @param outlierPvalThreshold Two-sided p-value below which a sample is an
#'   outlier. Default \code{1e-3} when unset.
#' @param outlierArgs A \code{\link{CovMcdOptions}} record of extra arguments
#'   for \code{robustbase::covMcd()}, the outlier detector.
#' @return A \code{PhenotypeFilterParam} object, a \code{\link{MethodParam}}.
#' @seealso \code{\link{GenotypeFilterParam}}
#' @examples
#' PhenotypeFilterParam(naAction = "impute", outlierAction = "drop")
#' @importFrom checkmate assertChoice
#' @export
PhenotypeFilterParam <- function(
    naAction = NULL,
    outlierAction = NULL,
    outlierPvalThreshold = NULL,
    outlierArgs = NULL
) {
    new(
        "PhenotypeFilterParam",
        naAction = naAction,
        outlierAction = outlierAction,
        outlierPvalThreshold = outlierPvalThreshold,
        outlierArgs = outlierArgs
    )
}

# Each set field one of its allowed values, checked where the bundle is
# built, so a bad value is reported against the argument that set it. An
# unset field has nothing to check.
setValidity("PhenotypeFilterParam", function(object) {
    coll <- makeAssertCollection()
    if (!is.null(object@naAction)) {
        assertChoice(
            object@naAction,
            c("drop", "impute", "keep"),
            .var.name = "naAction",
            add = coll
        )
    }
    if (!is.null(object@outlierAction)) {
        assertChoice(
            object@outlierAction,
            c("keep", "drop"),
            .var.name = "outlierAction",
            add = coll
        )
    }
    if (!is.null(object@outlierPvalThreshold)) {
        assertNumber(
            object@outlierPvalThreshold,
            lower = 0,
            upper = 1,
            .var.name = "outlierPvalThreshold",
            add = coll
        )
    }
    # Any MethodOptions fits the slot; only a covMcd record means anything to
    # the outlier detector, which would otherwise fail on foreign arguments.
    if (!is.null(object@outlierArgs)) {
        engine <- metadata(object@outlierArgs)$engine
        if (!identical(engine, "covMcd")) {
            coll$push(glue(
                "`outlierArgs` must be built with CovMcdOptions() (got a ",
                "record for {engine %||% 'no engine'})."
            ))
        }
    }
    problems <- coll$getMessages()
    if (length(problems) == 0L) TRUE else problems
})

# The run's phenotype QC with every unset field at its default, as a plain
# list the extraction reads. `naDefault` is the pipeline's NA policy and
# `naAllowed` what its fits can take; a choice outside it is refused here,
# naming the pipeline, rather than failing inside a fit.
# @noRd
.qtlResolvePhenotypeFilter <- function(
    phenotypeFilterParam,
    naDefault,
    naAllowed,
    caller
) {
    .assertMethodParam(phenotypeFilterParam, "PhenotypeFilterParam")
    naAction <- phenotypeFilterParam$naAction %||% naDefault
    if (!is_in(naAction, naAllowed)) {
        abort(glue(
            "{caller}: naAction = \"{naAction}\" is not available here; its ",
            "fits cannot take a missing phenotype. Use ",
            "{str_flatten(str_c('\"', naAllowed, '\"'), ' or ')}."
        ))
    }
    outlierArgs <- phenotypeFilterParam$outlierArgs %||% CovMcdOptions()
    .assertMethodOptions(outlierArgs, "CovMcdOptions", "outlierArgs")
    list(
        naAction = naAction,
        outlierAction = phenotypeFilterParam$outlierAction %||% "keep",
        outlierPvalThreshold = phenotypeFilterParam$outlierPvalThreshold %||%
            1e-3,
        outlierArgs = outlierArgs
    )
}

# What an accessor applies: nothing. Missing values and outliers are kept, so
# molecularTraits() and residualizedPhenotypes() return the data as stored.
# @noRd
.qtlKeepAllPhenotypes <- function() {
    list(
        naAction = "keep",
        outlierAction = "keep",
        outlierPvalThreshold = 1e-3,
        outlierArgs = CovMcdOptions()
    )
}

# =============================================================================
# Variant / sample filtering bundles
# -----------------------------------------------------------------------------
# Which rows and columns of the DATA survive, as opposed to the stage bundles
# above, which describe how a stage behaves. Three filters run in this
# package, on three different objects, and they are deliberately three
# constructors rather than one:
#
#   GenotypeFilterParam  -- a study's own genotype matrix (QtlDataset and the
#                          pipelines that build one)
#   PanelFilterParam     -- an LD reference panel's variants
#   SumstatsFilterParam  -- rows of a summary-statistics table
#
# Their fields overlap without meaning the same thing, and their defaults
# genuinely differ: `imissCutoff` resolves to 0 on the genotype path and 1 on
# the panel path. One union bundle would have to pick one of those, and would
# accept `removeIndels` where nothing reads it -- the failure mode these
# constructors exist to prevent.
#
# GenotypeFilterParam defaults every field to NULL ("not set"), resolved to
# its default by .qtlResolveFilter() when a run threads it to its
# extractions; the other two keep ordinary defaults.
#
# SumstatsCleaningParam sits with them because it runs on the same pass over
# the same table, but it is coercion and normalisation (coerceNumeric,
# normalizeChr, clampSmallP) that also drops rows -- not a filter.
# =============================================================================
