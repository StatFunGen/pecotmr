# =============================================================================
# GenotypeMixture S4 class
# -----------------------------------------------------------------------------
# A mixture of genotype panels used as one LD reference: a
# MultiAssayExperiment with one experiment per cohort, that cohort's panel.
# Its LD is the weighted average of the panels' LD, sum_k w_k R_k. The panels
# stay as they are -- lazy, on disk -- and are combined only when a step
# reads LD for some variants, so each keeps answering for its own samples:
# allele frequencies, sample size, missingness.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title Mixture of Genotype Panels
#' @description A \code{GenotypeMixture} uses several genotype panels (see
#'   \code{\link{readGenotypes}}) as one LD reference whose LD is the
#'   weighted average of theirs, \eqn{\sum_k w_k R_k}. That is the LD of
#'   summary statistics meta-analysed across cohorts, with each cohort's
#'   weight its share of the GWAS sample size.
#'
#'   Unlike an \code{\link{LdMatrix}} mixture, which is combined once when it
#'   is built, the panels stay lazy and on disk: LD is read from each and
#'   combined only for the variants a step asks about, and every panel keeps
#'   answering for its own samples. A step that needs allele frequencies,
#'   missingness or minor-allele counts judges each panel on its own and
#'   drops a variant that fails in any of them.
#'
#'   It \strong{is} a \code{MultiAssayExperiment} with one experiment per
#'   cohort: that cohort's genotype panel, cut to the variants every panel
#'   carries, in the first panel's order. The mixture's variants, as
#'   \code{\link{variantIds}} gives them, are the first panel's. Each panel
#'   keeps its own samples and allele coding, and its \code{rowData} column
#'   \code{mixtureSign} is \code{1} where it codes a variant's alleles as the
#'   first panel does and \code{-1} where it codes them the other way round.
#'   The \code{sampleMap} names each sample \code{<cohort>_<sampleId>}, since
#'   cohorts can reuse sample ids, and the weights are
#'   \code{metadata(x)$mixtureWeights} (see \code{\link{mixtureWeights}}).
#'
#'   Row subsetting, \code{x[i, ]} by position, variant id or \code{GRanges},
#'   cuts every panel to the same variants, and column subsetting narrows
#'   the panels' samples. A panel cannot be dropped, since the weights
#'   describe all of them: build a new \code{GenotypeMixture} over the panels
#'   wanted instead.
#' @seealso \code{\link{GenotypeMixture}} to build one, \code{\link{LdMatrix}}.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
setClass("GenotypeMixture", contains = "MultiAssayExperiment")

setValidity("GenotypeMixture", function(object) {
    problems <- .gmCheckComponents(object)
    if (length(problems) == 0L) {
        problems <- c(.gmCheckSign(object), .gmCheckWeights(object))
    }
    if (length(problems) == 0L) TRUE else problems
})

# @noRd
.gmCheckComponents <- function(object) {
    comps <- .gmComponents(object)
    if (length(comps) < 2L || !all(map_lgl(comps, .gmIsPanel))) {
        return("the experiments must be at least two genotype panels.")
    }
    if (length(unique(map_int(comps, nrow))) != 1L) {
        return("every panel must hold exactly the mixture's variants.")
    }
    character()
}

# A panel as a mixture holds it: the RSE readGenotypes() returns.
# @noRd
.gmIsPanel <- function(x) {
    methods::is(x, "RangedSummarizedExperiment")
}

# @noRd
.gmCheckSign <- function(object) {
    signs <- map(.gmComponents(object), .gmComponentSign)
    if (all(map_lgl(signs, .gmIsSign)) && all(signs[[1L]] == 1)) {
        return(character())
    }
    str_c(
        "every panel's rowData must carry `mixtureSign`, +1 or -1 per ",
        "variant, and the first panel's must be all +1."
    )
}

# @noRd
.gmIsSign <- function(s) {
    !is.null(s) && all(is_in(s, c(-1, 1)))
}

# @noRd
.gmCheckWeights <- function(object) {
    w <- mixtureWeights(object)
    named <- identical(names(w), names(.gmComponents(object)))
    c(
        .ldmCheckWeights(w),
        if (length(w) == 0L || !named) {
            "`mixtureWeights` must hold one weight per panel, named after it."
        }
    )
}

#' @title Build a GenotypeMixture
#' @description Combine genotype panels into one LD reference whose LD is
#'   the weighted average of theirs. See \code{\link{GenotypeMixture-class}}.
#'
#'   The panels are matched by position and alleles and put into the first
#'   panel's allele coding, so they need not share variant ids or allele
#'   orientation. The mixture covers the variants every panel carries; a
#'   variant some panel lacks is left out, with a message.
#' @param panels A named list of at least two genotype panels (see
#'   \code{\link{readGenotypes}}). The names are the cohorts: they name the
#'   mixture's experiments and prefix its samples.
#' @param mixtureWeights One non-negative weight per panel, summing to 1;
#'   typically each cohort's share of the GWAS sample size. Matched to the
#'   panels by name when named, by position otherwise.
#' @return A \code{GenotypeMixture}.
#' @examples
#' panel <- readGenotypes(plink2Prefix = file.path(system.file("extdata",
#'   "ld_reference", "chr22", package = "pecotmr"),
#'   "protocol_example.LD.chr22"))
#' half <- seq_len(ncol(panel) / 2)
#' mix <- GenotypeMixture(list(EUR = panel[, half], AFR = panel[, -half]),
#'   mixtureWeights = c(EUR = 0.7, AFR = 0.3))
#' mix
#' @rdname GenotypeMixture-constructor
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
GenotypeMixture <- function(panels, mixtureWeights) {
    .gmAssertPanels(panels)
    panels <- map(map(panels, .asLdRef), .ucscExperiment)
    if (!all(map_lgl(panels, .gmIsPanel))) {
        abort("Every element of `panels` must be a genotype panel.")
    }
    w <- .gmNamedWeights(mixtureWeights, names(panels))
    ref <- panels[[1L]]
    matched <- map(panels, .ldmMatchTo, ref = ref)
    common <- sort(reduce(map(matched, "refIdx"), intersect))
    if (length(common) == 0L) {
        abort("The genotype panels share no variant to combine over.")
    }
    .ldmReportDropped(
        length(.ldRefRanges(ref)) - length(common),
        length(common),
        what = "GenotypeMixture"
    )
    aligned <- map(matched, .ldmAlignedOn, common = common)
    experiments <- map2(panels, aligned, .gmAlignedComponent)
    mae <- MultiAssayExperiment(
        experiments = experiments,
        colData = .gmColData(experiments),
        sampleMap = .gmSampleMap(experiments),
        metadata = list(mixtureWeights = w)
    )
    methods::new("GenotypeMixture", mae)
}

# @noRd
.gmAssertPanels <- function(panels) {
    if (!is.list(panels) || length(panels) < 2L) {
        abort("`panels` must be a list of at least two genotype panels.")
    }
    nms <- names(panels)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
        abort(c(
            str_c(
                "`panels` must be a named list: each name is a cohort, which ",
                "names its experiment and prefixes its samples."
            ),
            i = "For example, list(EUR = eurPanel, AFR = afrPanel)."
        ))
    }
    invisible(NULL)
}

# The weights, one per panel and named after it: matched by name when they
# are named, by position otherwise.
# @noRd
.gmNamedWeights <- function(mixtureWeights, cohorts) {
    given <- names(mixtureWeights)
    w <- .ldmAssertWeights(mixtureWeights, length(cohorts))
    if (is.null(given)) {
        return(set_names(w, cohorts))
    }
    if (anyDuplicated(given) || !setequal(given, cohorts)) {
        abort(glue(
            "Named `mixtureWeights` must name each panel once: ",
            "{str_flatten(cohorts, ', ')}."
        ))
    }
    set_names(w, given)[cohorts]
}

# A panel cut to the shared variants in the first panel's order, recording
# how it codes each one's alleles relative to the first panel.
# @noRd
.gmAlignedComponent <- function(panel, aligned) {
    .gmWithSign(.ldRefSubset(panel, aligned$idx), aligned$sign)
}

# @noRd
.gmWithSign <- function(panel, sign) {
    rd <- SummarizedExperiment::rowData(panel)
    rd$mixtureSign <- as.numeric(sign)
    SummarizedExperiment::`rowData<-`(panel, value = rd)
}

# The primary name of each of a cohort's samples. Cohorts can reuse sample
# ids, and every stochastic sketch does, so the cohort is part of the name.
# @noRd
.gmPrimaryIds <- function(cohort, ids) {
    str_c(cohort, "_", ids)
}

# One row per sample of every panel, recording its cohort.
# @noRd
.gmColData <- function(experiments) {
    exec(rbind, !!!unname(imap(experiments, .gmColDataPart)))
}

# @noRd
.gmColDataPart <- function(panel, cohort) {
    ids <- as.character(colnames(panel))
    S4Vectors::DataFrame(
        cohort = rep(cohort, length(ids)),
        row.names = .gmPrimaryIds(cohort, ids)
    )
}

# The sampleMap: each panel's samples under their cohort-prefixed names.
# @noRd
.gmSampleMap <- function(experiments) {
    exec(rbind, !!!unname(imap(experiments, .gmSampleMapPart)))
}

# @noRd
.gmSampleMapPart <- function(panel, cohort) {
    ids <- as.character(colnames(panel))
    S4Vectors::DataFrame(
        assay = factor(rep(cohort, length(ids)), levels = cohort),
        primary = .gmPrimaryIds(cohort, ids),
        colname = ids
    )
}

# ---- subsetting -------------------------------------------------------------

# Row subsetting cuts every panel to the same variants. Every subscript is
# turned into positions on the first panel, whose ids and alleles the
# mixture's are, and those positions are applied to all of them: the panels'
# own ids can differ in allele order, so matching each on its own could cut
# them to different variants.
# @noRd
.gmSubsetByRow <- function(x, y, ...) {
    .gmRowsAt(x, .gmRowIndex(x, y))
}

# Registered for List and list as well as ANY: a GRanges is a List, and for
# either MultiAssayExperiment's own methods would otherwise win the dispatch
# and cut each panel on its own -- a GRanges by overlap with that panel, a
# list by one index per panel, which can leave the panels holding different
# variants.
#' @importFrom MultiAssayExperiment subsetByRow
#' @importClassesFrom S4Vectors List
setMethod("subsetByRow", c("GenotypeMixture", "ANY"), .gmSubsetByRow)
setMethod("subsetByRow", c("GenotypeMixture", "List"), .gmSubsetByRow)
setMethod("subsetByRow", c("GenotypeMixture", "list"), .gmSubsetByRow)

# The weights describe every panel, so none can be dropped.
#' @importFrom MultiAssayExperiment subsetByAssay
setMethod("subsetByAssay", "GenotypeMixture", function(x, y) {
    abort(c(
        str_c(
            "A GenotypeMixture's panels cannot be subset: its weights ",
            "describe all of them."
        ),
        i = "Build a new GenotypeMixture() over the panels wanted."
    ))
})

# The positions of the variants a row subscript names: a position, a logical
# mask, a variant id (the first panel's row names) or a GRanges.
# @noRd
#' @importFrom IRanges overlapsAny
.gmRowIndex <- function(x, y) {
    ranges <- .gmRanges(x)
    if (methods::is(y, "GenomicRanges")) {
        return(which(overlapsAny(ranges, y)))
    }
    if (is.character(y)) {
        at <- match(y, names(ranges))
        return(at[!is.na(at)])
    }
    if (is.logical(y) || is.numeric(y)) {
        return(seq_along(ranges)[y])
    }
    abort(str_c(
        "A GenotypeMixture's rows are subset by position, logical mask, ",
        "variant id or GRanges."
    ))
}

# The mixture with every panel cut to the variants at positions `i`, in
# that order.
# @noRd
#' @importFrom MultiAssayExperiment ExperimentList
.gmRowsAt <- function(x, i) {
    comps <- map(.gmComponents(x), .gmComponentRowsAt, i = i)
    MultiAssayExperiment::`experiments<-`(x, value = ExperimentList(comps))
}

# One panel's rows `i`. With no rows left it is held in memory with its
# samples, which the sampleMap names, so an empty mixture does not carry its
# panels' seeds (see .emptyPanel()).
# @noRd
.gmComponentRowsAt <- function(panel, i) {
    sign <- .gmComponentSign(panel)[i]
    rows <- if (length(i) == 0L) {
        .emptyPanel(panel, keepSamples = TRUE)
    } else {
        .ldRefSubset(panel, i)
    }
    .gmWithSign(rows, sign)
}

# ---- accessors --------------------------------------------------------------

#' @rdname mixtureWeights
#' @export
setMethod("mixtureWeights", "GenotypeMixture", function(x) {
    S4Vectors::metadata(x)$mixtureWeights
})

# The canonical chrom:pos:A2:A1 id every class renders; the panels' own labels
# stay in the SNP column of the mixture's ranges.
#' @rdname variantIds
#' @export
setMethod("variantIds", "GenotypeMixture", function(x, ...) {
    .grVariantIds(.gmRanges(x))
})

#' @rdname GenotypeMixture-class
#' @param object A \code{GenotypeMixture}.
#' @return \code{show} is called for its side effect and returns
#'   \code{invisible(NULL)}.
#' @export
setMethod("show", "GenotypeMixture", function(object) {
    w <- mixtureWeights(object)
    comps <- .gmComponents(object)
    cat(glue(
        "GenotypeMixture: {.gmNVariants(object)} variants from ",
        "{length(comps)} genotype panels\n",
        .trim = FALSE
    ))
    cat(glue(
        "  panels {str_flatten(names(comps), ', ')}; weights ",
        "{str_flatten(signif(w, 3), ', ')}; panel sizes ",
        "{str_flatten(map_int(comps, ncol), ', ')}\n",
        .trim = FALSE
    ))
})

# ---- the LD reference answers, combined over the panels ---------------------

# The panels, named after their cohorts.
# @noRd
#' @importFrom MultiAssayExperiment experiments
.gmComponents <- function(x) {
    as.list(experiments(x))
}

# @noRd
.gmComponentSign <- function(panel) {
    SummarizedExperiment::rowData(panel)$mixtureSign
}

# A +1 / -1 matrix with one row per variant and one column per panel: each
# panel's allele coding relative to the first's.
# @noRd
.gmSign <- function(x) {
    exec(cbind, !!!unname(map(.gmComponents(x), .gmComponentSign)))
}

# The mixture's variants: the first panel's ranges, with the SNP / A1 / A2
# columns that identify them.
# @noRd
.gmRanges <- function(x) {
    rr <- SummarizedExperiment::rowRanges(.gmComponents(x)[[1L]])
    S4Vectors::`mcols<-`(
        rr,
        value = S4Vectors::mcols(rr)[, c("SNP", "A1", "A2")]
    )
}

# @noRd
.gmNVariants <- function(x) {
    nrow(.gmComponents(x)[[1L]])
}

# 1 / sum(w^2 / n_k) over the panels' sizes as `nOf` measures them: the size
# of one panel whose LD would be as noisy as the weighted average.
# @noRd
.gmPooledN <- function(x, nOf) {
    n <- map_int(.gmComponents(x), nOf)
    as.integer(round(1 / sum(mixtureWeights(x)^2 / n)))
}

# The mixture's LD among its variants at `idx`: each panel's LD, re-signed
# into the first panel's allele coding, averaged with the weights. One
# panel's block is held at a time.
# @noRd
.gmLdAt <- function(x, idx) {
    reduce(
        seq_along(.gmComponents(x)),
        .gmAddLd,
        .init = NULL,
        x = x,
        idx = idx
    )
}

# @noRd
.gmAddLd <- function(acc, k, x, idx) {
    s <- .gmSign(x)[idx, k]
    term <- mixtureWeights(x)[[k]] *
        unname(.ldRefLdAt(.gmComponents(x)[[k]], idx)) *
        outer(s, s)
    if (is.null(acc)) term else acc + term
}

# Per-variant statistics a panel filter judges, from each panel on its own
# terms: a variant's MAF is its smallest in any panel and its missingness its
# largest, so a cutoff drops a variant that fails in any of them.
# @noRd
.gmVariantStats <- function(x, idx, needMissing, step) {
    each <- map(
        .gmComponents(x),
        .ldRefVariantStats,
        idx = idx,
        needMissing = needMissing,
        step = step
    )
    list(
        maf = exec(pmin, !!!unname(map(each, "maf"))),
        missRate = if (needMissing) exec(pmax, !!!unname(map(each, "missRate")))
    )
}

# The mixture's A1 frequency, sum_k w_k p_k with each panel's frequency put
# into the first panel's coding, or NULL unless every panel records one.
# @noRd
.gmAf <- function(x) {
    afs <- map(.gmComponents(x), .ldRefAf)
    if (any(map_lgl(afs, is.null))) {
        return(NULL)
    }
    signs <- .gmSign(x)
    aligned <- map(seq_along(afs), .gmAlignedAf, afs = afs, signs = signs)
    reduce(map2(unname(mixtureWeights(x)), aligned, `*`), `+`)
}

# @noRd
.gmAlignedAf <- function(k, afs, signs) {
    if_else(signs[, k] > 0, afs[[k]], 1 - afs[[k]])
}

# Per-variant genotype variance averaged over the panels with the weights --
# the average within-panel variance, which is what the averaged LD describes.
# The variance of the pooled frequency, 2p(1 - p), would also count the
# spread of frequencies between panels, which the mixture LD leaves out.
# @noRd
.gmVariance <- function(x, idx) {
    each <- map(.gmComponents(x), .ldRefVariance, idx = idx)
    reduce(map2(unname(mixtureWeights(x)), each, `*`), `+`)
}

# A factor of the mixture's LD over `idx`: each panel's standardized
# dosages, re-signed and scaled by sqrt(w_k / (m_k - 1)), stacked by rows,
# then scaled so the crossproduct is (N - 1) R for N rows -- the same
# convention as a single panel's scaled dosages. Its crossproduct is the
# weighted average of the panels' LD exactly, .gmLdAt()'s.
# @noRd
.gmFactor <- function(x, idx) {
    blocks <- map(seq_along(.gmComponents(x)), .gmFactorBlock, x = x, idx = idx)
    stacked <- exec(
        rbind,
        !!!blocks,
        !!!map2(blocks, unname(mixtureWeights(x)), .gmMonomorphicRows)
    )
    stacked * sqrt(nrow(stacked) - 1)
}

# A variant monomorphic in a panel has an all-zero column in that panel's
# block, which would leave the mixture's LD short of a unit diagonal. A
# panel's LD treats such a variant as independent of every other (see
# .ldCleanCorrelation()), so in the average it contributes w_k to the
# variant's own variance and nothing else: two rows per such variant, holding
# +/- sqrt(w_k / 2) in its column, so the column keeps a zero mean and an
# engine that centres the factor's rows itself still sees the same LD.
# @noRd
.gmMonomorphicRows <- function(block, w) {
    dead <- which(colSums(block^2) == 0)
    .ldPaddingRows(dead, ncol(block), sqrt(w / 2))
}

# @noRd
.gmFactorBlock <- function(k, x, idx) {
    scaled <- .ldRefStandardized(.gmComponents(x)[[k]], idx)
    signed <- sweep(scaled, 2L, .gmSign(x)[idx, k], `*`)
    signed * sqrt(mixtureWeights(x)[[k]] / (nrow(scaled) - 1))
}

# TRUE when two mixtures average the same panels -- the same samples, in the
# same order -- with the same weights. Their variant sets may differ, as two
# QC'd subsets of one mixture's do.
# @noRd
.gmSameMixture <- function(x, y) {
    cx <- .gmComponents(x)
    cy <- .gmComponents(y)
    length(cx) == length(cy) &&
        isTRUE(all.equal(mixtureWeights(x), mixtureWeights(y))) &&
        identical(
            map(cx, .ldRefSampleIds),
            map(cy, .ldRefSampleIds)
        )
}
