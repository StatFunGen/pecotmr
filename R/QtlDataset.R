# =============================================================================
# QtlDataset S4 class
# -----------------------------------------------------------------------------
# Single-study individual-level QTL container, built on MultiAssayExperiment:
# one RangedSummarizedExperiment per QTL context, plus a `genotype` experiment
# whose assay reads lazily through a GenotypeHandle. Adds the handle itself and
# constructor-level QC knobs (mafCutoff, macCutoff, xvarCutoff, imissCutoff,
# keepVariants, keepIndel), which the getGenotypes / getResidualizedGenotypes
# accessors apply at extraction time. The entry point for individual-level
# fine-mapping (fineMappingPipeline), TWAS weight learning
# (twasWeightsPipeline), and multi-study composition (MultiStudyQtlDataset).
# =============================================================================

#' @include AllGenerics.R
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
#' @slot study Character (length 1). Study identifier; used in collection
#'   classes to tag downstream \code{FineMappingResult} / \code{TwasWeights}
#'   entries.
#'   The \code{genotype} experiment's assay reads through this handle; the
#'   extraction accessors read it directly, so that QC can be applied per
#'   block.
#' @slot scaleResiduals Logical (length 1). Whether residualization accessors
#'   scale residuals to unit variance.
#' @slot mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time inside
#'   \code{getGenotypes()} / \code{getResidualizedGenotypes()}. Default 0 (no
#'   filter).
#' @slot macCutoff Numeric (length 1). Minor allele count threshold; converted
#'   to a MAF threshold using \code{max(mafCutoff, macCutoff / (2 * n))} where
#'   \code{n} is the post-narrowing sample count of the extracted block. Default
#'   0 (no filter).
#' @slot xvarCutoff Numeric (length 1). Per-variant genotype variance threshold;
#'   variants with column variance below this are dropped at extraction time.
#'   Default 0 (no filter).
#' @slot imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are dropped at
#'   extraction time. Default 0 (no filter).
#' @slot keepVariants Character vector of variant identifiers to retain prior to
#'   per-block QC. Length 0 means no restriction.
#' @slot keepIndel Logical (length 1). When \code{FALSE}, indel variants
#'   (alleles that are not single nucleotides) are dropped at extraction time.
#'   Default \code{TRUE}.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
setClass(
    "QtlDataset",
    contains = "MultiAssayExperiment",
    representation(
        study = "character",
        scaleResiduals = "logical",
        mafCutoff = "numeric",
        macCutoff = "numeric",
        xvarCutoff = "numeric",
        imissCutoff = "numeric",
        keepVariants = "character",
        keepIndel = "logical"
    ),
    prototype = prototype(keepIndel = TRUE),
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

# Validity for QtlDataset: scalar/cutoff slots, the phenotypes list, and
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

# study / scaleResiduals / keepIndel scalars + the four non-negative cutoffs.
# @noRd
.qtlValidateScalars <- function(object) {
    errors <- character()
    if (length(object@study) != 1L || str_length(object@study) == 0L) {
        errors <- c(
            errors,
            "'study' must be a single non-empty character string"
        )
    }
    if (length(object@scaleResiduals) != 1L) {
        errors <- c(errors, "'scaleResiduals' must be a single logical value")
    }
    if (length(object@keepIndel) != 1L || is.na(object@keepIndel)) {
        errors <- c(errors, "'keepIndel' must be a single logical value")
    }
    for (nm in c("mafCutoff", "macCutoff", "xvarCutoff", "imissCutoff")) {
        v <- methods::slot(object, nm)
        if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 0) {
            errors <- c(
                errors,
                glue("'{nm}' must be a single finite non-negative numeric")
            )
        }
    }
    errors
}

# Shape checks on the phenotype list handed to the constructor. Run before
# the MultiAssayExperiment is assembled, because a malformed list makes
# ExperimentList fail first with a message about experiments rather than
# about contexts. Returns an error vector.
# @noRd
.qtlCheckPhenotypeList <- function(phenotypes) {
    errors <- character()
    if (length(phenotypes) == 0L) {
        errors <- c(errors, "'phenotypes' must not be empty")
    }
    contextNames <- names(phenotypes)
    if (
        is.null(contextNames) ||
            any(str_length(contextNames) == 0L, na.rm = TRUE) ||
            any(is.na(contextNames))
    ) {
        errors <- c(
            errors,
            "'phenotypes' must be a named list with non-empty names"
        )
    } else if (n_distinct(contextNames) < length(contextNames)) {
        errors <- c(errors, "context names in 'phenotypes' must be unique")
    } else if (is_in(.QTL_GENO_EXPERIMENT, contextNames)) {
        errors <- c(
            errors,
            glue(
                "'{.QTL_GENO_EXPERIMENT}' is reserved for the genotype ",
                "experiment and cannot name a context"
            )
        )
    }
    for (ctx in seq_along(phenotypes)) {
        se <- phenotypes[[ctx]]
        if (!methods::is(se, "SummarizedExperiment")) {
            errors <- c(
                errors,
                glue(
                    "phenotypes[[{ctx}]] must be a SummarizedExperiment ",
                    "(got {class(se)[[1L]]})"
                )
            )
        } else if (is.null(colnames(se))) {
            errors <- c(
                errors,
                glue(
                    "phenotypes[[{ctx}]] has no column names. Which samples ",
                    "a context observes is recorded in the sampleMap, so ",
                    "every context must name its samples"
                )
            )
        }
    }
    errors
}

# The MAE must carry the genotype experiment plus at least one context.
# @noRd
.qtlValidatePhenotypes <- function(object) {
    exps <- MultiAssayExperiment::experiments(object)
    errors <- character()
    if (!is_in(.QTL_GENO_EXPERIMENT, names(exps))) {
        errors <- c(
            errors,
            glue("experiment '{.QTL_GENO_EXPERIMENT}' is missing")
        )
    }
    if (length(.qtlPhenotypeList(object)) == 0L) {
        errors <- c(errors, "'phenotypes' must not be empty")
    }
    errors
}

# TRUE when two GRanges share canonical chrom + start + end.
# @noRd
.qtlSameRange <- function(prev, this) {
    isTRUE(all.equal(
        canonChrom(GenomicRanges::seqnames(prev)),
        canonChrom(GenomicRanges::seqnames(this))
    )) &&
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
    errors <- character()
    traitToRange <- list()
    for (ctx in seq_along(pheno)) {
        se <- pheno[[ctx]]
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- rownames(se)
        if (length(rr) != length(ids)) {
            next
        }
        for (i in seq_along(ids)) {
            tid <- ids[[i]]
            prev <- traitToRange[[tid]]
            if (is.null(prev)) {
                traitToRange[[tid]] <- rr[i]
            } else if (!.qtlSameRange(prev, rr[i])) {
                errors <- c(
                    errors,
                    glue(
                        "trait '{tid}' has inconsistent rowRanges across ",
                        "contexts"
                    )
                )
            }
        }
    }
    errors
}

# =============================================================================
# QtlDataset constructor and accessors
# =============================================================================

# Check the constructor's inputs and return the handle behind `genotypes`.
# Every phenotype must be a usable SummarizedExperiment, and the genotype
# source must be an unsubset panel (or a bare handle) -- see
# .qtlCheckWholePanel for why the "unsubset" part matters.
# @noRd
.qtlValidateInputs <- function(phenotypes, genotypes) {
    errors <- .qtlCheckPhenotypeList(phenotypes)
    if (length(errors) > 0L) {
        abort(str_flatten(errors, "\n"))
    }
    handle <- .openGenotypeHandle(genotypes)
    if (is.null(handle)) {
        abort(glue(
            "'genotypes' must be a genotype panel from readGenotypes() ",
            "(got {class(genotypes)[[1L]]})"
        ))
    }
    .qtlCheckWholePanel(genotypes, handle)
    handle
}

# A QtlDataset indexes variants positionally into the handle's whole snpInfo
# (see .qtlVariantIndices), so a panel that has already been narrowed would
# leave extraction reading the wrong rows. Refuse it rather than silently
# widening back to the file: narrowing belongs to keepVariants / keepSamples,
# which the extraction path honours.
# @noRd
.qtlCheckWholePanel <- function(genotypes, handle) {
    if (!methods::is(genotypes, "RangedSummarizedExperiment")) {
        return(invisible(NULL))
    }
    full <- c(nrow(getSnpInfo(handle)), getNSamples(handle))
    got <- c(nrow(genotypes), ncol(genotypes))
    if (identical(as.integer(got), as.integer(full))) {
        return(invisible(NULL))
    }
    abort(glue(
        "'genotypes' is a subset panel ({got[[1L]]} x {got[[2L]]} of ",
        "{full[[1L]]} x {full[[2L]]}); pass the whole panel and narrow with ",
        "'keepVariants' / 'keepSamples' instead."
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
#' @param study Character (length 1). Study identifier.
#' @param genotypes A genotype panel (see \code{\link{readGenotypes}}).
#' @param phenotypes Named list of \code{SummarizedExperiment} objects, keyed by
#'   context. Each SE must have \code{rowRanges} carrying trait positions and
#'   \code{colData} carrying per-context phenotype covariates. The name
#'   \code{"genotype"} is reserved.
#' @param genotypeCovariates Numeric matrix of genotype-derived covariates
#'   (e.g., ancestry PCs); rows are samples. Becomes the \code{colData} of the
#'   genotype experiment.
#' @param scaleResiduals Logical (length 1). Default \code{TRUE}.
#' @param mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time inside
#'   \code{getGenotypes()} / \code{getResidualizedGenotypes()}. Default 0 (no
#'   filter).
#' @param macCutoff Numeric (length 1). Minor allele count threshold; converted
#'   to a MAF threshold using \code{max(mafCutoff, macCutoff / (2 * n))} where
#'   \code{n} is the post-narrowing sample count of the extracted block. Default
#'   0 (no filter).
#' @param xvarCutoff Numeric (length 1). Per-variant genotype variance
#'   threshold; variants with column variance below this are dropped at
#'   extraction time. Default 0 (no filter).
#' @param imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are dropped at
#'   extraction time. Default 0 (no filter).
#' @param keepSamples Character vector of sample identifiers to retain. The
#'   dataset is subset to them, narrowing \code{colData} and
#'   \code{sampleMap} together, so the sample set has one home rather than
#'   two. Length 0 means no restriction.
#' @param keepVariants Character vector of variant identifiers to retain prior
#'   to per-block QC. Length 0 means no restriction.
#' @param keepIndel Logical (length 1). When \code{FALSE}, variants whose
#'   alleles are not single nucleotides (indels) are dropped at extraction.
#'   Default \code{TRUE} (keep all variants).
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
#' QtlDataset(study = "s1", genotypes = panel, phenotypes = list(brain = se))
#' @export
QtlDataset <- function(
    study,
    genotypes,
    phenotypes,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0),
    scaleResiduals = TRUE,
    mafCutoff = 0,
    macCutoff = 0,
    xvarCutoff = 0,
    imissCutoff = 0,
    keepSamples = character(0),
    keepVariants = character(0),
    keepIndel = TRUE
) {
    handle <- .qtlValidateInputs(phenotypes, genotypes)
    experiments <- c(
        set_names(
            list(.genotypeExperiment(genotypes, genotypeCovariates)),
            .QTL_GENO_EXPERIMENT
        ),
        as.list(phenotypes)
    )
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiments,
        colData = .qtlPrimaryColData(experiments),
        sampleMap = .qtlSampleMap(experiments)
    )
    obj <- methods::new(
        "QtlDataset",
        .qtlRestrictSamples(mae, keepSamples),
        study = as.character(study),
        scaleResiduals = isTRUE(scaleResiduals),
        mafCutoff = as.numeric(mafCutoff),
        macCutoff = as.numeric(macCutoff),
        xvarCutoff = as.numeric(xvarCutoff),
        imissCutoff = as.numeric(imissCutoff),
        keepVariants = as.character(keepVariants),
        keepIndel = isTRUE(keepIndel)
    )
    validObject(obj)
    obj
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
    .qtlMuffleDrop(MultiAssayExperiment::subsetByAssay(mae, getContexts(x)))
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
    exps <- MultiAssayExperiment::experiments(mae)
    se <- exps[[.QTL_GENO_EXPERIMENT]]
    SummarizedExperiment::assay(se, "dosage") <- as.matrix(
        SummarizedExperiment::assay(se, "dosage")
    )
    exps[[.QTL_GENO_EXPERIMENT]] <- se
    MultiAssayExperiment::experiments(mae) <- exps
    mae
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
# is no second copy to keep in step; getGenotypeHandle() reads it back.
# @noRd
.qtlWithGenotypeHandle <- function(x, handle) {
    exps <- MultiAssayExperiment::experiments(x)
    gCov <- .qtlColDataMatrix(exps[[.QTL_GENO_EXPERIMENT]])
    exps[[.QTL_GENO_EXPERIMENT]] <- .genotypeExperiment(handle, gCov)
    MultiAssayExperiment::experiments(x) <- exps
    validObject(x)
    x
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

#' @rdname getStudy
#' @export
setMethod("getStudy", "QtlDataset", function(x) x@study)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlDataset", function(x) {
    names(.qtlPhenotypeList(x))
})

#' @rdname getGenotypeCovariates
#' @export
setMethod("getGenotypeCovariates", "QtlDataset", function(x) {
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
    out <- as.matrix(as.data.frame(cd))
    rownames(out) <- rownames(cd)
    out
}

#' @rdname getScaleResiduals
#' @export
setMethod("getScaleResiduals", "QtlDataset", function(x) x@scaleResiduals)

#' @rdname getGenotypeHandle
#' @keywords internal
setMethod("getGenotypeHandle", "QtlDataset", function(x) {
    # Derived, not stored. The handle already lives inside the genotype
    # assay's seed -- that is what lets the dosages read lazily -- so a
    # parallel slot was a second copy that had to be kept in step by hand.
    # Reading it back removes the invariant instead of policing it.
    .ldSketchHandle(
        MultiAssayExperiment::experiments(x)[[.QTL_GENO_EXPERIMENT]]
    )
})

#' @rdname qtlDatasetFilters
#' @export
setMethod("getMafCutoff", "QtlDataset", function(x, ...) x@mafCutoff)

#' @rdname qtlDatasetFilters
#' @export
setMethod("getMacCutoff", "QtlDataset", function(x, ...) x@macCutoff)

#' @rdname qtlDatasetFilters
#' @export
setMethod("getXvarCutoff", "QtlDataset", function(x, ...) x@xvarCutoff)

#' @rdname qtlDatasetFilters
#' @export
setMethod("getImissCutoff", "QtlDataset", function(x, ...) x@imissCutoff)

#' @rdname qtlDatasetFilters
#' @export
setMethod("getKeepVariants", "QtlDataset", function(x, ...) x@keepVariants)

#' @rdname qtlDatasetFilters
#' @export
setMethod("getKeepIndel", "QtlDataset", function(x, ...) x@keepIndel)

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

# The union span (+/- cisWindow) of a trait's rowRanges across all contexts.
# Requires cisWindow and a single shared chromosome.
# @noRd
.qtlTraitRegion <- function(x, traitId, cisWindow) {
    if (is.null(cisWindow) || length(cisWindow) != 1L || cisWindow < 0) {
        msg <- glue(
            "`cisWindow` is required (and must be non-negative) when ",
            "`traitId` is specified."
        )
        abort(msg)
    }
    perTraitRanges <- list()
    for (ctx in getContexts(x)) {
        se <- getPhenotypes(x, ctx)
        hits <- match(traitId, rownames(se))
        hits <- hits[!is.na(hits)]
        if (length(hits) > 0) {
            rr <- SummarizedExperiment::rowRanges(se)
            perTraitRanges[[length(perTraitRanges) + 1L]] <- rr[hits]
        }
    }
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
.qtlLiteralRegion <- function(region, cisWindow) {
    if (!methods::is(region, "GRanges")) {
        abort("`region` must be a GRanges object.")
    }
    if (length(region) == 0L) {
        abort("`region` must contain at least one range.")
    }
    if (is.null(cisWindow)) {
        return(region)
    }
    if (length(cisWindow) != 1L || cisWindow < 0) {
        abort("`cisWindow` must be a single non-negative value.")
    }
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
# traitId (a 0-width chrUn sentinel for a trait absent from every context), for
# threading trait-position provenance onto QtlFineMappingResult / TwasWeights.
# @noRd
.qtlTraitPos <- function(x, traitIds) {
    # Extract chrom/start/end as plain vectors per trait (union span across
    # contexts), then build ONE fresh GRanges at the end. Combining per-context
    # GRanges with do.call(c, .) can trip S4 seqinfo reconciliation in some
    # GenomeInfoDb builds, so we avoid it entirely.
    n <- length(traitIds)
    chrs <- rep("chrUn", n)
    starts <- rep(1L, n)
    ends <- rep(1L, n)
    for (i in seq_len(n)) {
        tid <- traitIds[[i]]
        st <- Inf
        en <- -Inf
        ch <- NA_character_
        for (ctx in getContexts(x)) {
            se <- getPhenotypes(x, ctx)
            h <- match(tid, rownames(se))
            h <- h[!is.na(h)]
            if (length(h) == 0L) {
                next
            }
            rr <- SummarizedExperiment::rowRanges(se)[h]
            ch <- as.character(GenomicRanges::seqnames(rr))[1L]
            st <- min(st, GenomicRanges::start(rr))
            en <- max(en, GenomicRanges::end(rr))
        }
        if (!is.na(ch)) {
            chrs[i] <- ch
            starts[i] <- as.integer(st)
            ends[i] <- as.integer(en)
        }
    }
    gr <- GenomicRanges::GRanges(
        chrs,
        IRanges::IRanges(start = starts, end = pmax(ends, starts))
    )
    names(gr) <- traitIds
    gr
}

#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "QtlDataset", function(x, traitId = NULL, ...) {
    tids <- if (is.null(traitId)) {
        unique(unlist(map(.qtlPhenotypeList(x), rownames)))
    } else {
        as.character(traitId)
    }
    .qtlTraitPos(x, tids)
})

# Internal: map a GRanges region (one or more ranges) into 1-based snpIdx into
# handle@snpInfo. Indices are unioned across ranges in range order (first
# occurrence wins), so overlapping ranges contribute each variant once.
.qtlVariantIndices <- function(x, region = NULL) {
    handle <- getGenotypeHandle(x)
    if (is.null(region)) {
        return(seq_len(nrow(getSnpInfo(handle))))
    }
    snpInfo <- getSnpInfo(handle)
    siChr <- canonChrom(snpInfo$CHR)
    bp <- as.integer(snpInfo$BP)
    rChr <- canonChrom(GenomicRanges::seqnames(region))
    rStart <- GenomicRanges::start(region)
    rEnd <- GenomicRanges::end(region)
    idx <- integer(0)
    for (i in seq_along(region)) {
        idx <- c(idx, which(siChr == rChr[i] & bp >= rStart[i] & bp <= rEnd[i]))
    }
    unique(idx)
}

# Internal: keepIndel slot read, tolerant of QtlDataset objects serialized
# before the slot existed (treat a missing slot as TRUE = keep indels).
.qtlKeepIndel <- function(x) {
    isTRUE(tryCatch(getKeepIndel(x), error = function(e) TRUE))
}

# Internal: return a copy of a QtlDataset with the supplied filter cutoffs /
# keep-lists REPLACING the stored slot values (NULL = leave the stored value
# untouched). This lets a pipeline accept per-call filter overrides as ordinary
# arguments instead of forcing callers to mutate @slots directly (which bypasses
# the class's validity checks). Applied against a validated copy.
.qtlApplyFilterOverrides <- function(
    data,
    mafCutoff = NULL,
    macCutoff = NULL,
    xvarCutoff = NULL,
    imissCutoff = NULL,
    keepIndel = NULL,
    keepSamples = NULL,
    keepVariants = NULL
) {
    if (!is.null(mafCutoff)) {
        data@mafCutoff <- as.numeric(mafCutoff)
    }
    if (!is.null(macCutoff)) {
        data@macCutoff <- as.numeric(macCutoff)
    }
    if (!is.null(xvarCutoff)) {
        data@xvarCutoff <- as.numeric(xvarCutoff)
    }
    if (!is.null(imissCutoff)) {
        data@imissCutoff <- as.numeric(imissCutoff)
    }
    if (!is.null(keepIndel)) {
        data@keepIndel <- as.logical(keepIndel)
    }
    if (!is.null(keepSamples)) {
        data <- .qtlRestrictSamples(data, keepSamples)
    }
    if (!is.null(keepVariants)) {
        data@keepVariants <- as.character(keepVariants)
    }
    methods::validObject(data)
    data
}

# Internal: extract the panel dosage block (samples x variants) for the
# requested region, narrow to the requested sample set, and apply lazy QC
# (per-sample imiss filter, then per-variant max(mafCutoff,
# macCutoff / (2 * n)) and xvarCutoff filters). Used by getGenotypes,
# getResidualizedGenotypes (via getGenotypes), and getMaf so all three
# share a single variant/sample selection result.
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
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    samples = NULL
) {
    gr <- .qtlResolveVariantRegion(
        x,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow
    )
    snpIdx <- .qtlVariantIndices(x, gr)
    if (length(snpIdx) == 0L) {
        return(.qtlEmptyBlockAllSamples(x))
    }
    # Apply keepVariants + indel restrictions before materialization so we do
    # not extract dosage we would immediately drop.
    snpIdx <- .qtlNarrowSnpIdx(x, snpIdx)
    if (length(snpIdx) == 0L) {
        return(.qtlEmptyBlock())
    }
    handle <- getGenotypeHandle(x)
    dosage <- .dosageMatrix(handle, snpIdx, meanImpute = FALSE)
    keep <- .qtlResolveSamples(dosage, x, samples)
    if (length(keep) == 0L) {
        return(.qtlEmptyBlockNoSamples(dosage))
    }
    dosage <- dosage[keep, , drop = FALSE]
    # Per-sample missingness filter.
    if (getImissCutoff(x) > 0 && nrow(dosage) > 0L && ncol(dosage) > 0L) {
        dosage <- dosage[
            rowMeans(is.na(dosage)) <= getImissCutoff(x),
            ,
            drop = FALSE
        ]
    }
    filtered <- .qtlVariantFilters(dosage, x)
    dosage <- .qtlMeanImpute(filtered$dosage)
    list(
        geno = dosage,
        variantIds = colnames(dosage),
        sampleIds = rownames(dosage),
        maf = filtered$maf,
        af = filtered$af
    )
}

# Empty block preserving the full panel sample set (no variants selected).
# @noRd
.qtlEmptyBlockAllSamples <- function(x) {
    handle <- getGenotypeHandle(x)
    list(
        geno = matrix(
            numeric(0),
            nrow = getNSamples(handle),
            ncol = 0L,
            dimnames = list(getSampleIds(handle), character(0))
        ),
        variantIds = character(0),
        sampleIds = getSampleIds(handle),
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
.qtlNarrowSnpIdx <- function(x, snpIdx) {
    handle <- getGenotypeHandle(x)
    if (length(getKeepVariants(x)) > 0L) {
        snpAll <- as.character(getSnpInfo(handle)$SNP[snpIdx])
        km <- matchVariants(snpAll, as.character(getKeepVariants(x)))
        keepMask <- logical(length(snpAll))
        keepMask[km$idxA] <- TRUE
        snpIdx <- snpIdx[keepMask]
    }
    if (length(snpIdx) > 0L && !.qtlKeepIndel(x)) {
        si <- getSnpInfo(handle)
        # which() (not the mask) so an NA mask drops the variant rather than
        # injecting an NA index.
        snpMask <- str_length(as.character(si$A1[snpIdx])) == 1L &
            str_length(as.character(si$A2[snpIdx])) == 1L
        snpIdx <- snpIdx[which(snpMask)]
    }
    snpIdx
}

# Resolve the sample set: panel samples intersected with the dataset's
# primary colData -- which is what narrows a subset dataset -- and with the
# per-call `samples` arg.
# @noRd
.qtlResolveSamples <- function(dosage, x, samples) {
    keep <- intersect(
        rownames(dosage),
        rownames(MultiAssayExperiment::colData(x))
    )
    if (!is.null(samples)) {
        keep <- intersect(keep, as.character(samples))
    }
    keep
}

# Per-variant MAF / MAC / X-variance filters against the post-narrowing sample
# count. MAF is computed from the un-imputed dosage (A1 = the effect allele
# gives directional `af`; `maf` folds to the minor allele). Returns
# list(dosage, maf, af).
# @noRd
.qtlVariantFilters <- function(dosage, x) {
    nSamp <- nrow(dosage)
    if (ncol(dosage) == 0L) {
        # nocov start (unreachable: zero-variant paths return early above)
        return(list(dosage = dosage, maf = numeric(0), af = numeric(0)))
        # nocov end
    }
    nObs <- colSums(!is.na(dosage))
    sumD <- colSums(dosage, na.rm = TRUE)
    p <- if_else(nObs > 0L, sumD / (2 * nObs), NA_real_)
    afVec <- p
    mafVec <- pmin(p, 1 - p)
    effectiveMaf <- max(
        getMafCutoff(x),
        if (nSamp > 0L) getMacCutoff(x) / (2 * nSamp) else 0
    )
    keepVarMask <- !is.na(mafVec) & mafVec >= effectiveMaf
    if (getXvarCutoff(x) > 0 && nSamp > 1L) {
        mu <- if_else(nObs > 0L, sumD / nObs, 0)
        centered <- sweep(dosage, 2L, mu, FUN = "-")
        centered[is.na(centered)] <- 0
        varVec <- colSums(centered * centered) / (nSamp - 1L)
        keepVarMask <- keepVarMask & varVec >= getXvarCutoff(x)
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
    for (j in seq_len(ncol(dosage))) {
        col <- dosage[, j]
        na <- is.na(col)
        if (any(na)) {
            col[na] <- mean(col[!na])
            dosage[, j] <- col
        }
    }
    dosage
}

#' @rdname getGenotypes
#' @export
setMethod(
    "getGenotypes",
    "QtlDataset",
    function(
        x,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        ...
    ) {
        .qtlExtractBlock(
            x,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )$geno
    }
)

#' @rdname getMaf
#' @export
setMethod(
    "getMaf",
    "QtlDataset",
    function(x, region = NULL, cisWindow = NULL, samples = NULL, ...) {
        block <- .qtlExtractBlock(
            x,
            traitId = NULL,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )
        out <- block$maf
        names(out) <- block$variantIds
        out
    }
)

#' @rdname getAf
#' @export
setMethod(
    "getAf",
    "QtlDataset",
    function(
        x,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        ...
    ) {
        block <- .qtlExtractBlock(
            x,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )
        out <- block$af
        names(out) <- block$variantIds
        out
    }
)

#' @rdname getPhenotypes
#' @export
setMethod(
    "getPhenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        naAction = c("keep", "drop", "impute"),
        outlierAction = c("keep", "drop"),
        outlierPvalThreshold = 1e-3,
        ...
    ) {
        naAction <- arg_match(naAction)
        outlierAction <- arg_match(outlierAction)
        .qtlValidateContexts(x, contexts)
        out <- .qtlPhenotypeList(x)[contexts]
        out <- .qtlFilterPhenotypes(
            out,
            contexts,
            traitId,
            region,
            naAction,
            outlierAction,
            outlierPvalThreshold
        )
        if (length(contexts) == 1L) out[[1L]] else out
    }
)

# `contexts` is required and must all be known phenotype contexts.
# @noRd
.qtlValidateContexts <- function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        msg <- glue(
            "`contexts` is required for getPhenotypes(QtlDataset). Pass a ",
            "character vector of one or more context names; use ",
            "getContexts(x) to list the available contexts."
        )
        abort(msg)
    }
    available <- getContexts(x)
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

# Apply the optional trait / region / NA / outlier filters to the per-context
# phenotype SEs.
# @noRd
.qtlFilterPhenotypes <- function(
    out,
    contexts,
    traitId,
    region,
    naAction,
    outlierAction,
    outlierPvalThreshold
) {
    if (!is.null(traitId)) {
        out <- .qtlFilterTraits(out, contexts, traitId)
    }
    if (!is.null(region)) {
        out <- set_names(
            map(out, .qtlSeInRegion, region = region),
            contexts
        )
    }
    if (naAction != "keep") {
        out <- set_names(
            map(out, .qtlApplyPhenoNaAction, naAction = naAction),
            contexts
        )
    }
    if (outlierAction != "keep") {
        out <- set_names(
            map(
                out,
                .qtlApplyPhenoOutliers,
                action = outlierAction,
                pvalThreshold = outlierPvalThreshold
            ),
            contexts
        )
    }
    out
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
        keepSamp <- colSums(is.na(Y)) == 0L
        se <- se[, keepSamp, drop = FALSE]
    } else if (naAction == "impute") {
        if (anyNA(Y)) {
            for (j in seq_len(nrow(Y))) {
                row <- Y[j, ]
                na <- is.na(row)
                if (any(na)) {
                    obs <- row[!na]
                    row[na] <- if (length(obs) > 0L) mean(obs) else 0
                    Y[j, ] <- row
                }
            }
            SummarizedExperiment::assay(se, assayName) <- Y
        }
    }
    se
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
# Significance: per-sample chi-squared(p) p-value with Bonferroni
# correction over the sample count. A sample is flagged when its
# corrected p-value falls below `pvalThreshold`. With single-trait Y
# (ncol == 1) this reduces to the standard z-test on (y - center)/sd.
#
# Returns all-TRUE (no-op) when there are too few samples to support
# a covariance estimate (n < p + 2).
.qtlOutlierKeepMask <- function(Y, pvalThreshold) {
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
    if (requireNamespace("robustbase", quietly = TRUE)) {
        mcd <- tryCatch(robustbase::covMcd(Y), error = function(e) NULL)
        if (!is.null(mcd)) {
            ctr <- mcd$center
            covMat <- mcd$cov
        } else {
            ctr <- colMeans(Y)
            covMat <- stats::cov(Y)
        }
    } else {
        msg <- glue(
            "outlier detection: install 'robustbase' for an MCD-based ",
            "estimator; falling back to non-robust colMeans/cov."
        )
        inform(msg)
        ctr <- colMeans(Y)
        covMat <- stats::cov(Y)
    }
    invCov <- tryCatch(solve(covMat), error = function(e) MASS::ginv(covMat))
    Yc <- sweep(Y, 2L, ctr)
    d2 <- rowSums((Yc %*% invCov) * Yc)
    raw <- stats::pchisq(d2, df = p, lower.tail = FALSE)
    raw >= (pvalThreshold / n)
}

# Wrapper: apply the keep-mask to a SummarizedExperiment slice. SE
# columns are samples; transpose the assay (traits x samples) before
# calling .qtlOutlierKeepMask which expects samples x traits.
.qtlApplyPhenoOutliers <- function(se, action, pvalThreshold) {
    if (action == "keep") {
        return(se)
    }
    assayName <- SummarizedExperiment::assayNames(se)[[1L]]
    Y <- t(SummarizedExperiment::assay(se, assayName))
    keep <- .qtlOutlierKeepMask(Y, pvalThreshold)
    if (all(keep)) {
        return(se)
    }
    se[, keep, drop = FALSE]
}

#' @rdname getPhenotypeCovariates
#' @export
setMethod("getPhenotypeCovariates", "QtlDataset", function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        abort("`contexts` is required.")
    }
    available <- getContexts(x)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    out <- map(contexts, .qtlContextColData, x = x)
    names(out) <- contexts
    out
})

# Internal: residualize a numeric matrix Y (n x k) against a covariate
# matrix C (n x p) via pivoted QR decomposition. Adds an intercept column
# to C. When C is rank-deficient (e.g., union of all contexts' phenotype
# covariates includes collinear / duplicate columns), the pivoted QR drops
# the redundant columns automatically. Optionally rescales each residual
# column to unit standard deviation; constant-valued columns are left
# unchanged.
.qtlResidualizeQr <- function(Y, C, scaleResiduals = TRUE) {
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
    res <- stats::lm.fit(x = X, y = Y)$residuals
    res <- as.matrix(res)
    rownames(res) <- rownames(Y)
    colnames(res) <- colnames(Y)
    if (isTRUE(scaleResiduals)) {
        sds <- apply(res, 2L, stats::sd, na.rm = TRUE)
        # `sds == 0` exact-zero test is unreliable for residuals coming out of
        # lm.fit on a constant Y: roundoff gives sd ~ 1e-16 instead of 0, and
        # dividing the (also-tiny) residuals by it amplifies floating-point
        # noise to unit-scale. Treat anything below sqrt(.Machine$double.eps)
        # as effectively zero (column is constant) and skip rescaling.
        nearZero <- !is.finite(sds) | sds < sqrt(.Machine$double.eps)
        sds[nearZero] <- 1
        res[, nearZero] <- 0
        res <- sweep(res, 2L, sds, FUN = "/")
    }
    res
}

# Resolve the phenotype-covariate selection for one context: NULL requested ->
# all available covariates; otherwise validate the requested names are present.
# @noRd
.qtlResolveOne <- function(ctx, requested, x) {
    se <- getPhenotypes(x, ctx)
    avail <- colnames(SummarizedExperiment::colData(se))
    if (is.null(requested)) {
        return(avail)
    }
    keep <- intersect(requested, avail)
    if (length(keep) != length(requested)) {
        missingNames <- setdiff(requested, avail)
        msg <- glue(
            "phenotypeCovariatesToResidualize: context '{ctx}' has no ",
            "covariate(s) named: {str_flatten(missingNames, ', ')}"
        )
        abort(msg)
    }
    keep
}

# Internal: validate and resolve the `*ToResidualize` argument against a
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
        "phenotypeCovariatesToResidualize must be NULL, a character vector, ",
        "or a named list keyed by context."
    )
    abort(msg)
}

# List-form phenotypeCovariatesToResidualize: must be named with EXACTLY the
# `contexts` set. Resolves each context's selection.
# @noRd
.qtlPhenoSelectionList <- function(x, contexts, toResidualize) {
    toResNames <- names(toResidualize)
    if (
        is.null(toResNames) || any(str_length(toResNames) == 0L, na.rm = TRUE)
    ) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: when supplied as a list, it ",
            "must be named with context names."
        )
        abort(msg)
    }
    badKeys <- setdiff(names(toResidualize), contexts)
    if (length(badKeys) > 0L) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: list key(s) not in ",
            "`contexts`: {str_flatten(badKeys, ', ')}"
        )
        abort(msg)
    }
    missingKeys <- setdiff(contexts, names(toResidualize))
    if (length(missingKeys) > 0L) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: list does not cover all ",
            "`contexts`. Per-context lists must have exactly the same context ",
            "set as `contexts`. Missing keys: ",
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
    avail <- colnames(getGenotypeCovariates(x))
    if (is.null(avail)) {
        avail <- character(0)
    }
    if (is.null(toResidualize)) {
        return(avail)
    }
    keep <- intersect(toResidualize, avail)
    if (length(keep) != length(toResidualize)) {
        missingNames <- setdiff(toResidualize, avail)
        msg <- glue(
            "genotypeCovariatesToResidualize: no covariate(s) named: ",
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
        getGenotypeCovariates(x)[, genoSelection, drop = FALSE]
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
# Sample names survive the coercion because the constructor requires every
# context to name its samples, which SummarizedExperiment carries into
# rownames(colData(se)).
# @noRd
.qtlPhenoCovBlocks <- function(x, contexts, phenoSelection) {
    perContext <- list()
    for (ctx in contexts) {
        keep <- phenoSelection[[ctx]]
        if (length(keep) == 0L) {
            next
        }
        se <- getPhenotypes(x, ctx)
        cd <- as.matrix(as.data.frame(SummarizedExperiment::colData(se)))
        cdMat <- cd[, keep, drop = FALSE]
        colnames(cdMat) <- str_c(ctx, ".", colnames(cdMat))
        perContext[[ctx]] <- cdMat
    }
    perContext
}

# Intersect the covariate blocks to their common samples and column-bind them
# into one design matrix; NULL when no samples are shared.
# @noRd
.qtlAlignCovariates <- function(perContext, gCov) {
    sampleSets <- map(perContext, .qtlBlockSamples)
    if (!is.null(gCov) && ncol(gCov) > 0L) {
        sampleSets <- c(sampleSets, list(.qtlBlockSamples(gCov)))
    }
    common <- if (length(sampleSets) == 0L) {
        character(0)
    } else {
        reduce(sampleSets, intersect)
    }
    if (length(common) == 0L) {
        return(NULL)
    }
    blocks <- list()
    for (mat in perContext) {
        blocks[[length(blocks) + 1L]] <- mat[common, , drop = FALSE]
    }
    if (!is.null(gCov) && ncol(gCov) > 0L) {
        blocks[[length(blocks) + 1L]] <- gCov[common, , drop = FALSE]
    }
    exec(cbind, !!!blocks)
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
    for (j in seq_len(ncol(C))) {
        col <- C[, j]
        na <- is.na(col)
        if (any(na)) {
            mu <- mean(col[!na])
            col[na] <- if (is.finite(mu)) mu else 0
            C[, j] <- col
        }
    }
    C
}

# Internal: resolve a (convenience, precise) flag pair to a single boolean.
# `missing*` arguments are passed as the result of `missing()` evaluated in
# the calling method to detect whether the user explicitly set the value.
# Rules:
#   - both missing: returns TRUE (the documented default)
#   - only convenience set: returns convenience
#   - only precise set: returns precise
#   - both set: must agree, else error
.qtlResolveResidualizationFlag <- function(
    conveniencePassed,
    convenienceMissing,
    precisePassed,
    preciseMissing,
    convenienceName,
    preciseName
) {
    if (preciseMissing && convenienceMissing) {
        return(TRUE)
    }
    if (preciseMissing) {
        return(isTRUE(conveniencePassed))
    }
    if (convenienceMissing) {
        return(isTRUE(precisePassed))
    }
    if (isTRUE(conveniencePassed) != isTRUE(precisePassed)) {
        msg <- glue(
            "Conflicting values: `{convenienceName}` = {conveniencePassed} ",
            "and `{preciseName}` = {precisePassed}. ",
            "Set only one, or pass consistent values."
        )
        abort(msg)
    }
    isTRUE(precisePassed)
}

#' @rdname getResidualizedGenotypes
#' @export
setMethod(
    "getResidualizedGenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        residualizePhenotypeCovariatesFromGenotypes = NULL,
        residualizeGenotypeCovariatesFromGenotypes = NULL,
        covariateNaAction = c("impute", "drop"),
        ...
    ) {
        if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
            msg <- glue(
                "`contexts` is required for ",
                "getResidualizedGenotypes(QtlDataset). ",
                "Use getContexts(x) to list the available contexts. ",
                "Pass a single context for per-context mode or multiple ",
                "contexts for joint mode (sample intersection)."
            )
            abort(msg)
        }
        covariateNaAction <- arg_match(covariateNaAction)
        convPhenoMissing <- missing(residualizePhenotypeCovariates)
        convGenoMissing <- missing(residualizeGenotypeCovariates)
        precPhenoMissing <- missing(
            residualizePhenotypeCovariatesFromGenotypes
        ) ||
            is.null(residualizePhenotypeCovariatesFromGenotypes)
        precGenoMissing <- missing(
            residualizeGenotypeCovariatesFromGenotypes
        ) ||
            is.null(residualizeGenotypeCovariatesFromGenotypes)
        p <- as.list(environment())
        p$dots <- list(...)
        .qtlResidualizedGenotypesImpl(p)
    }
)

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

# getResidualizedGenotypes worker: resolve the covariate inclusion flags +
# selections, extract genotypes, build + NA-handle the covariate design, and
# QR-residualize. `p` holds the setMethod args + the precomputed `missing()`
# booleans.
# @noRd
# The two residualization flags, each reconciled against its convenience and
# precision spellings. Paired here because the reconciliation rule is the same
# for both and only the argument names differ.
# @noRd
.qtlResidualizationFlags <- function(p) {
    list(
        pheno = .qtlResolveResidualizationFlag(
            p$residualizePhenotypeCovariates,
            p$convPhenoMissing,
            p$residualizePhenotypeCovariatesFromGenotypes,
            p$precPhenoMissing,
            "residualizePhenotypeCovariates",
            "residualizePhenotypeCovariatesFromGenotypes"
        ),
        geno = .qtlResolveResidualizationFlag(
            p$residualizeGenotypeCovariates,
            p$convGenoMissing,
            p$residualizeGenotypeCovariatesFromGenotypes,
            p$precGenoMissing,
            "residualizeGenotypeCovariates",
            "residualizeGenotypeCovariatesFromGenotypes"
        )
    )
}

.qtlResidualizedGenotypesImpl <- function(p) {
    bad <- setdiff(p$contexts, getContexts(p$x))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    include <- .qtlResidualizationFlags(p)
    includePheno <- include$pheno
    includeGeno <- include$geno
    phenoSel <- .qtlResolvePhenoSelection(
        p$x,
        p$contexts,
        p$phenotypeCovariatesToResidualize
    )
    genoSel <- .qtlResolveGenoSelection(p$x, p$genotypeCovariatesToResidualize)
    G <- getGenotypes(
        p$x,
        traitId = p$traitId,
        region = p$region,
        cisWindow = p$cisWindow,
        samples = p$samples
    )
    if (ncol(G) == 0L) {
        return(G)
    }
    C <- .qtlBuildResidualizationDesign(
        p$x,
        contexts = p$contexts,
        phenoSelection = phenoSel,
        genoSelection = genoSel,
        includePheno = includePheno,
        includeGeno = includeGeno
    )
    C <- .qtlHandleCovariateNa(C, p$covariateNaAction)
    aligned <- .qtlAlignGC(G, C, p$contexts)
    .qtlResidualizeQr(
        aligned$G,
        aligned$C,
        scaleResiduals = getScaleResiduals(p$x)
    )
}

#' @rdname getResidualizedPhenotypes
#' @export
setMethod(
    "getResidualizedPhenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        residualizePhenotypeCovariatesFromPhenotypes = NULL,
        residualizeGenotypeCovariatesFromPhenotypes = NULL,
        naAction = c("keep", "drop", "impute"),
        covariateNaAction = c("impute", "drop"),
        outlierAction = c("keep", "drop"),
        outlierPvalThreshold = 1e-3,
        ...
    ) {
        if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
            abort("`contexts` is required for getResidualizedPhenotypes().")
        }
        naAction <- arg_match(naAction)
        covariateNaAction <- arg_match(covariateNaAction)
        outlierAction <- arg_match(outlierAction)
        convPhenoMissing <- missing(residualizePhenotypeCovariates)
        convGenoMissing <- missing(residualizeGenotypeCovariates)
        precPhenoMissing <- missing(
            residualizePhenotypeCovariatesFromPhenotypes
        ) ||
            is.null(residualizePhenotypeCovariatesFromPhenotypes)
        precGenoMissing <- missing(
            residualizeGenotypeCovariatesFromPhenotypes
        ) ||
            is.null(residualizeGenotypeCovariatesFromPhenotypes)
        p <- as.list(environment())
        p$dots <- list(...)
        .qtlResidualizedPhenotypesImpl(p)
    }
)

# NA-handled raw phenotypes per context (re-wrapped to a list so single- and
# multi-context callers see the same shape).
# @noRd
.qtlResidPhenoY <- function(x, contexts, traitId, region, naAction) {
    Yraw <- getPhenotypes(
        x,
        contexts = contexts,
        traitId = traitId,
        region = region,
        naAction = naAction
    )
    if (length(contexts) == 1L) {
        Yraw <- set_names(list(Yraw), contexts)
    }
    Yraw
}

# Residualize one context's phenotypes against the covariate design (intersected
# to common samples) and drop residual-scale outliers.
# @noRd
.qtlResidualizeContextPheno <- function(
    se,
    C,
    ctx,
    outlierAction,
    outlierPvalThreshold,
    scaleResiduals
) {
    Y <- t(SummarizedExperiment::assay(se)) # samples x traits
    Cctx <- NULL
    if (!is.null(C)) {
        common <- intersect(rownames(Y), rownames(C))
        if (length(common) == 0L) {
            msg <- glue(
                "context '{ctx}': no samples shared between phenotype data ",
                "and the resolved covariate matrix."
            )
            abort(msg)
        }
        Y <- Y[common, , drop = FALSE]
        Cctx <- C[common, , drop = FALSE]
    }
    Yres <- .qtlResidualizeQr(Y, Cctx, scaleResiduals = scaleResiduals)
    # Outlier detection on the residualized scale.
    if (outlierAction != "keep") {
        keep <- .qtlOutlierKeepMask(Yres, outlierPvalThreshold)
        if (!all(keep)) {
            Yres <- Yres[keep, , drop = FALSE]
        }
    }
    Yres
}

# Resolve the phenotype/genotype covariate inclusion flags (convenience vs
# precise `*FromPhenotypes`) for getResidualizedPhenotypes.
# @noRd
.qtlResidPhenoFlags <- function(p) {
    list(
        includePheno = .qtlResolveResidualizationFlag(
            p$residualizePhenotypeCovariates,
            p$convPhenoMissing,
            p$residualizePhenotypeCovariatesFromPhenotypes,
            p$precPhenoMissing,
            "residualizePhenotypeCovariates",
            "residualizePhenotypeCovariatesFromPhenotypes"
        ),
        includeGeno = .qtlResolveResidualizationFlag(
            p$residualizeGenotypeCovariates,
            p$convGenoMissing,
            p$residualizeGenotypeCovariatesFromPhenotypes,
            p$precGenoMissing,
            "residualizeGenotypeCovariates",
            "residualizeGenotypeCovariatesFromPhenotypes"
        )
    )
}

# getResidualizedPhenotypes worker: resolve covariate inclusion + selections,
# NA-handle Y, build the covariate design, and per-context residualize +
# outlier-filter. `p` holds the setMethod args + precomputed missing() flags.
# @noRd
.qtlResidualizedPhenotypesImpl <- function(p) {
    bad <- setdiff(p$contexts, getContexts(p$x))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    flags <- .qtlResidPhenoFlags(p)
    includePheno <- flags$includePheno
    includeGeno <- flags$includeGeno
    phenoSel <- .qtlResolvePhenoSelection(
        p$x,
        p$contexts,
        p$phenotypeCovariatesToResidualize
    )
    genoSel <- .qtlResolveGenoSelection(p$x, p$genotypeCovariatesToResidualize)
    Yraw <- .qtlResidPhenoY(
        p$x,
        p$contexts,
        p$traitId,
        p$region,
        p$naAction
    )
    C <- .qtlBuildResidualizationDesign(
        p$x,
        contexts = p$contexts,
        phenoSelection = phenoSel,
        genoSelection = genoSel,
        includePheno = includePheno,
        includeGeno = includeGeno
    )
    C <- .qtlHandleCovariateNa(C, p$covariateNaAction)
    out <- set_names(
        map(p$contexts, .qtlResidualizeContext, Yraw = Yraw, C = C, p = p),
        p$contexts
    )
    if (length(p$contexts) == 1L) out[[1L]] else out
}


#' @rdname show-methods
#' @export
setMethod("show", "QtlDataset", function(object) {
    pheno <- .qtlPhenotypeList(object)
    nCtx <- length(pheno)
    ctxNames <- names(pheno)
    totalTraits <- length(unique(unlist(
        map(pheno, rownames),
        use.names = FALSE
    )))
    cat(glue("QtlDataset for study '{object@study}'\n", .trim = FALSE))
    cat(glue(
        "  {nCtx} context(s): {str_flatten(ctxNames, ', ')}\n",
        .trim = FALSE
    ))
    cat(glue("  {totalTraits} unique traits across contexts\n", .trim = FALSE))
    gh <- getGenotypeHandle(object)
    cat(glue(
        "  Genotypes: {getFormat(gh)} @ {getPath(gh)}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Genotype covariates: ",
        "{ncol(getGenotypeCovariates(object))} cols\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Samples: {nrow(MultiAssayExperiment::colData(object))}\n",
        .trim = FALSE
    ))
    cat(glue("  Scale residuals: {object@scaleResiduals}\n", .trim = FALSE))
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

# One context's SE restricted to features overlapping `region`.
# @noRd
.qtlSeInRegion <- function(se, region) {
    rr <- SummarizedExperiment::rowRanges(se)
    se[IRanges::overlapsAny(rr, region), , drop = FALSE]
}

# One context's covariate matrix (colData of its phenotype SE).
# @noRd
.qtlContextColData <- function(ctx, x) {
    se <- getPhenotypes(x, ctx)
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
.qtlResidualizeContext <- function(ctx, Yraw, C, p) {
    .qtlResidualizeContextPheno(
        Yraw[[ctx]],
        C,
        ctx,
        p$outlierAction,
        p$outlierPvalThreshold,
        getScaleResiduals(p$x)
    )
}
