# =============================================================================
# LdMatrix S4 class
# -----------------------------------------------------------------------------
# An LD correlation matrix the caller supplies, over the variants it
# describes. This is one of the two ways LD enters pecotmr. The other is a
# genotype panel (readGenotypes()), which is read lazily from disk; pecotmr
# never reads a pre-computed matrix from disk, so a matrix is built in memory
# by the caller and handed over whole.
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title LD Correlation Matrix
#' @description An \code{LdMatrix} carries one LD correlation matrix together
#'   with the variants it describes, and can stand wherever pecotmr takes an
#'   LD reference (the \code{ldReference} of a summary-statistics collection,
#'   \code{\link{raiss}}, the heritability builders). It \strong{is} a
#'   \code{GRanges} over those variants, in the order of the matrix's rows,
#'   with metadata columns \code{SNP} (the variant id), \code{A1} (the allele
#'   whose dosage the correlations are signed by), \code{A2} and, optionally,
#'   \code{AF} (the frequency of \code{A1}).
#'
#'   The other kind of LD reference is a genotype panel (see
#'   \code{\link{readGenotypes}}), which pecotmr reads from disk on demand. A
#'   matrix carries no genotypes, so the steps that need them are not
#'   available with one: the panel missingness filter, and RAiSS's
#'   singular-value route (RAiSS falls back to the matrix itself). Steps that
#'   need the panel's sample size or allele frequencies -- MAF / MAC panel
#'   filters, DENTIST, ldsc's unbiased \eqn{r^2}, HDL -- use \code{nSamples}
#'   and \code{AF} when they are supplied and refuse otherwise.
#'
#'   Validity checks the matrix's shape and labels, that its diagonal is 1,
#'   that its entries lie in \eqn{[-1, 1]}, and its symmetry on an evenly
#'   spaced subset of entries. It does not check positive semi-definiteness,
#'   which costs an eigendecomposition; run \code{\link{checkLd}} for that.
#'
#' @slot correlation The correlation matrix, its dimnames the variant ids.
#' @slot nSamples Integer, the number of independent samples the LD is
#'   worth, or \code{NA} when unknown. For one panel this is its sample size.
#'   For a mixture it is \eqn{1 / \sum_k w_k^2 / n_k}, the size of a single
#'   panel whose LD would be as noisy.
#' @slot mixtureWeights The weights a mixture's components were combined
#'   with; empty for a single panel's LD.
#' @seealso \code{\link{LdMatrix}} to build one, \code{\link{checkLd}}.
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom S4Vectors extractROWS normalizeSingleBracketSubscript
#' @export
setClass(
    "LdMatrix",
    contains = "GRanges",
    slots = c(
        correlation = "matrix",
        nSamples = "integer",
        mixtureWeights = "numeric"
    ),
    prototype = prototype(
        correlation = matrix(numeric(0), 0L, 0L),
        nSamples = NA_integer_,
        mixtureWeights = numeric(0)
    )
)

setValidity("LdMatrix", function(object) {
    problems <- c(
        .ldmCheckColumns(object),
        .ldmCheckShape(object),
        .ldmCheckValues(object@correlation),
        .ldmCheckAf(object),
        .ldmCheckN(object@nSamples),
        .ldmCheckWeights(object@mixtureWeights)
    )
    if (length(problems) == 0L) TRUE else problems
})

# @noRd
.ldmCheckColumns <- function(object) {
    missingCols <- setdiff(
        c("SNP", "A1", "A2"),
        colnames(S4Vectors::mcols(object))
    )
    if (length(missingCols) == 0L) {
        return(character())
    }
    str_c(
        "the variants carry no ",
        str_flatten(missingCols, ", "),
        " column(s)."
    )
}

# The matrix must be square over exactly the object's variants, labelled by
# their ids: every consumer looks LD up by id, so a matrix whose labels drift
# from the variant table would hand back some other variant's LD.
# @noRd
.ldmCheckShape <- function(object) {
    R <- object@correlation
    p <- length(object)
    if (!is.numeric(R) || nrow(R) != p || ncol(R) != p) {
        return(glue(
            "`correlation` is {nrow(R)} x {ncol(R)} but there are {p} ",
            "variant(s); it must be numeric and square over them."
        ))
    }
    # An empty matrix has no labels to check: subsetting to zero variants
    # leaves NULL dimnames, which is as labelled as nothing gets.
    ids <- as.character(S4Vectors::mcols(object)$SNP)
    if (p == 0L) {
        return(character())
    }
    if (
        !identical(rownames(R), ids) ||
            !identical(colnames(R), ids)
    ) {
        return(
            "`correlation`'s row and column names must be the variant ids."
        )
    }
    character()
}

# A correlation matrix has a unit diagonal, entries in [-1, 1] and is
# symmetric. Symmetry is checked on an evenly spaced block of rows and
# columns: comparing a large matrix with its transpose costs another copy of
# it, and the strided block catches a wrong input just as well.
# @noRd
.ldmCheckValues <- function(R, tol = 1e-6) {
    p <- nrow(R)
    if (!is.numeric(R) || p == 0L || nrow(R) != ncol(R)) {
        return(character())
    }
    if (anyNA(R)) {
        return("`correlation` holds missing values.")
    }
    span <- range(R)
    problems <- c(
        if (span[[1L]] < -1 - tol || span[[2L]] > 1 + tol) {
            "`correlation` holds entries outside [-1, 1]."
        },
        if (max(abs(diag(R) - 1)) > tol) {
            "`correlation`'s diagonal is not 1."
        }
    )
    probe <- unique(round(seq(1, p, length.out = min(p, 200L))))
    block <- R[probe, probe, drop = FALSE]
    if (max(abs(block - t(block))) > tol) {
        problems <- c(problems, "`correlation` is not symmetric.")
    }
    problems
}

# @noRd
.ldmCheckAf <- function(object) {
    af <- S4Vectors::mcols(object)$AF
    if (is.null(af)) {
        return(character())
    }
    if (!is.numeric(af) || any(af < 0 | af > 1, na.rm = TRUE)) {
        return("`AF` must be numeric allele frequencies in [0, 1].")
    }
    character()
}

# @noRd
.ldmCheckN <- function(n) {
    if (length(n) == 1L && (is.na(n) || n >= 1L)) {
        return(character())
    }
    "`nSamples` must be one positive count, or NA when unknown."
}

# @noRd
.ldmCheckWeights <- function(w) {
    if (length(w) == 0L) {
        return(character())
    }
    if (anyNA(w) || any(w < 0) || abs(sum(w) - 1) > 1e-8) {
        return("`mixtureWeights` must be non-negative and sum to 1.")
    }
    character()
}

# Documented on its own page, LdMatrix-constructor.Rd: the default LdMatrix.Rd
# would collide with the ldMatrix() generic's ldMatrix.Rd on a
# case-insensitive filesystem (macOS), and one would overwrite the other.
#' @title Build an LdMatrix
#' @description Wrap an LD correlation matrix, or combine several into a
#'   mixture, for use as an LD reference. See \code{\link{LdMatrix-class}}.
#'
#'   \strong{One panel's LD.} Pass the matrix and the variants it covers.
#'   The variants may be given as a \code{GRanges} (metadata columns
#'   \code{A1}, \code{A2}, and optionally \code{SNP} and \code{AF}), as a
#'   data frame with columns \code{chrom}, \code{pos}, \code{A1}, \code{A2}
#'   (and optionally \code{SNP}, \code{AF}), or omitted when the matrix's
#'   dimnames are \code{chrom:pos:A2:A1} ids, from which they are parsed. The
#'   correlations must be signed by the dosage of \code{A1}.
#'
#'   \strong{A mixture.} Pass a list of \code{LdMatrix} objects and their
#'   \code{mixtureWeights}. The result is their weighted average
#'   \eqn{\sum_k w_k R_k}, computed once here, over the variants every
#'   component covers; a variant some component lacks is left out, with a
#'   message. Components are matched by position and alleles and put into
#'   the first component's allele coding, so they need not share variant ids
#'   or allele orientation. The mixture's \code{AF} is \eqn{\sum_k w_k AF_k}
#'   when every component records \code{AF}, and its \code{nSamples} is
#'   \eqn{1 / \sum_k w_k^2 / n_k} when every component records a size.
#'
#' @param correlation A numeric correlation matrix, or a list of
#'   \code{LdMatrix} objects to combine into a mixture.
#' @param variants The variants the matrix covers, in its row order: a
#'   \code{GRanges} or a data frame (see Description), or \code{NULL} to
#'   parse them from the matrix's dimnames. Not used for a mixture.
#' @param nSamples The size of the panel the LD was estimated from, or
#'   \code{NA} when unknown. Not used for a mixture.
#' @param mixtureWeights For a mixture, one non-negative weight per
#'   component, summing to 1.
#' @return An \code{LdMatrix}.
#' @rdname LdMatrix-constructor
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' ld <- LdMatrix(cor(X), nSamples = nrow(X))
#' ld
#' ldMatrix(ld)[1:3, 1:3]
#'
#' # A mixture of two panels' LD, weighted 3:1.
#' half <- seq_len(nrow(X) / 2)
#' ldA <- LdMatrix(cor(X[half, ]), nSamples = length(half))
#' ldB <- LdMatrix(cor(X[-half, ]), nSamples = nrow(X) - length(half))
#' LdMatrix(list(ldA, ldB), mixtureWeights = c(0.75, 0.25))
#' @importFrom checkmate assertNumber
#' @export
LdMatrix <- function(
    correlation,
    variants = NULL,
    nSamples = NA_integer_,
    mixtureWeights = NULL
) {
    if (is.list(correlation) && !is.data.frame(correlation)) {
        return(.ldmMixture(correlation, mixtureWeights))
    }
    if (!is.null(mixtureWeights)) {
        abort(c(
            "`mixtureWeights` is for combining LD matrices.",
            i = "Pass a list of LdMatrix objects as `correlation`."
        ))
    }
    assertNumber(nSamples, na.ok = TRUE, lower = 1)
    R <- .ldmAsMatrix(correlation)
    gr <- .ldmVariants(variants, R)
    .ldmNew(gr, .ldmLabelled(R, gr), nSamples, numeric(0))
}

# @noRd
.ldmNew <- function(gr, R, nSamples, mixtureWeights) {
    methods::new(
        "LdMatrix",
        gr,
        correlation = R,
        nSamples = as.integer(round(nSamples)),
        mixtureWeights = as.numeric(mixtureWeights)
    )
}

# A plain numeric matrix: a Matrix-package matrix is made dense, as every
# consumer works on one.
# @noRd
#' @importFrom checkmate assertMatrix
.ldmAsMatrix <- function(correlation) {
    if (methods::is(correlation, "Matrix")) {
        correlation <- as.matrix(correlation)
    }
    assertMatrix(
        correlation,
        mode = "numeric",
        any.missing = FALSE,
        .var.name = "correlation"
    )
    correlation
}

# The matrix labelled by the variant ids. A matrix that already carries
# labels must carry exactly these, in this order: a mismatch means the matrix
# and the variant table describe different variants, or the same ones in a
# different order.
# @noRd
.ldmLabelled <- function(R, gr) {
    ids <- as.character(S4Vectors::mcols(gr)$SNP)
    if (nrow(R) != length(ids) || ncol(R) != length(ids)) {
        abort(glue(
            "`correlation` is {nrow(R)} x {ncol(R)} but `variants` lists ",
            "{length(ids)} variant(s)."
        ))
    }
    labelled <- !is.null(rownames(R)) || !is.null(colnames(R))
    if (
        labelled &&
            (!identical(rownames(R), ids) || !identical(colnames(R), ids))
    ) {
        abort(c(
            "`correlation`'s dimnames do not match the variant ids.",
            i = str_c(
                "Order the matrix as the variants are ordered, or drop its ",
                "dimnames to label it by `variants`."
            )
        ))
    }
    `dimnames<-`(R, list(ids, ids))
}

# The variant table, normalised to the GRanges an LdMatrix is: chr-prefixed
# seqnames (as a genotype panel's), names and `SNP` the ids, `A1` / `A2`, and
# `AF` when supplied.
# @noRd
.ldmVariants <- function(variants, R) {
    tbl <- if (is.null(variants)) {
        .ldmVariantsFromIds(rownames(R))
    } else if (methods::is(variants, "GRanges")) {
        .ldmVariantsFromRanges(variants, rownames(R))
    } else if (is.data.frame(variants)) {
        .ldmVariantsFromTable(variants, rownames(R))
    } else {
        abort(glue(
            "`variants` must be a GRanges, a data frame or NULL ",
            "(got {class(variants)[[1L]]})."
        ))
    }
    .ldmRangesOf(tbl)
}

# @noRd
.ldmVariantsFromIds <- function(ids) {
    if (is.null(ids)) {
        abort(c(
            "`variants` is required when the matrix has no dimnames.",
            i = "Pass the variants, or label the matrix chrom:pos:A2:A1."
        ))
    }
    parsed <- parseVariantId(ids)
    if (anyNA(parsed$pos) || anyNA(parsed$A1)) {
        abort(c(
            "The matrix's dimnames are not chrom:pos:A2:A1 variant ids.",
            i = "Pass `variants` to say which variants it covers."
        ))
    }
    tibble(
        SNP = ids,
        chrom = parsed$chrom,
        pos = parsed$pos,
        A1 = parsed$A1,
        A2 = parsed$A2
    )
}

# @noRd
.ldmVariantsFromRanges <- function(gr, ids) {
    mc <- as.data.frame(S4Vectors::mcols(gr))
    .ldmVariantsFromTable(
        mutate(
            as_tibble(mc),
            chrom = as.character(GenomicRanges::seqnames(gr)),
            pos = GenomicRanges::start(gr),
            SNP = .ldmIdsOf(mc, names(gr), ids)
        ),
        ids
    )
}

# @noRd
.ldmVariantsFromTable <- function(tbl, ids) {
    missingCols <- setdiff(c("chrom", "pos", "A1", "A2"), colnames(tbl))
    if (length(missingCols) > 0L) {
        abort(glue(
            "`variants` lacks column(s) {str_flatten(missingCols, ', ')}."
        ))
    }
    snp <- if (is_in("SNP", colnames(tbl))) {
        as.character(tbl$SNP)
    } else {
        .ldmIdsOf(tbl, NULL, ids)
    }
    out <- tibble(
        SNP = snp %||%
            formatVariantId(
                as.character(tbl$chrom),
                as.integer(tbl$pos),
                as.character(tbl$A2),
                as.character(tbl$A1)
            ),
        chrom = as.character(tbl$chrom),
        pos = as.integer(tbl$pos),
        A1 = as.character(tbl$A1),
        A2 = as.character(tbl$A2)
    )
    if (is_in("AF", colnames(tbl))) {
        mutate(out, AF = as.numeric(tbl$AF))
    } else {
        out
    }
}

# The ids a variant table names, in order of preference: its own SNP or
# variant_id column, then the ranges' names, then the matrix's dimnames.
# NULL when none of those exist.
# @noRd
.ldmIdsOf <- function(tbl, rangeNames, ids) {
    if (is_in("SNP", colnames(tbl))) {
        return(as.character(tbl$SNP))
    }
    if (is_in("variant_id", colnames(tbl))) {
        return(as.character(tbl$variant_id))
    }
    if (!is.null(rangeNames)) {
        return(as.character(rangeNames))
    }
    ids
}

# @noRd
.ldmRangesOf <- function(tbl) {
    if (anyDuplicated(tbl$SNP) > 0L) {
        abort("The variant ids must be unique.")
    }
    gr <- GenomicRanges::GRanges(
        seqnames = withChrPrefix(as.character(tbl$chrom)),
        ranges = IRanges::IRanges(as.integer(tbl$pos), width = 1L)
    )
    S4Vectors::`mcols<-`(
        `names<-`(gr, tbl$SNP),
        value = S4Vectors::DataFrame(
            select(tbl, -all_of(c("chrom", "pos"))),
            row.names = NULL
        )
    )
}

# ---- mixtures ---------------------------------------------------------------

# The weighted average of several LD matrices, over the variants they all
# cover, in the first component's allele coding.
# @noRd
.ldmMixture <- function(components, mixtureWeights) {
    .ldmAssertComponents(components)
    w <- .ldmAssertWeights(mixtureWeights, length(components))
    ref <- components[[1L]]
    matched <- map(components, .ldmMatchTo, ref = ref)
    common <- sort(reduce(map(matched, "refIdx"), intersect))
    if (length(common) == 0L) {
        abort("The LD matrices share no variant to combine over.")
    }
    .ldmReportDropped(length(ref) - length(common), length(common))
    aligned <- map(matched, .ldmAlignedOn, common = common)
    R <- reduce(
        seq_along(components),
        .ldmAddComponent,
        .init = NULL,
        components = components,
        aligned = aligned,
        w = w
    )
    ids <- rownames(ref@correlation)[common]
    gr <- .ldmMixtureVariants(ref[common], components, aligned, w)
    methods::new(
        "LdMatrix",
        gr,
        correlation = `dimnames<-`(R, list(ids, ids)),
        nSamples = .ldmMixtureN(map_int(components, nSamples), w),
        mixtureWeights = w
    )
}

# @noRd
.ldmAssertComponents <- function(components) {
    isLd <- map_lgl(components, methods::is, class2 = "LdMatrix")
    if (length(components) < 2L || !all(isLd)) {
        abort(c(
            "A mixture needs a list of at least two LdMatrix objects.",
            i = str_c(
                "Build each component with LdMatrix(), then combine them ",
                "with LdMatrix(list(...), mixtureWeights = ...)."
            )
        ))
    }
}

# @noRd
.ldmAssertWeights <- function(w, k) {
    if (
        is.null(w) ||
            !is.numeric(w) ||
            length(w) != k ||
            anyNA(w) ||
            any(w < 0) ||
            abs(sum(w) - 1) > 1e-8
    ) {
        abort(glue(
            "`mixtureWeights` must be {k} non-negative weights summing to 1, ",
            "one per component."
        ))
    }
    as.numeric(w)
}

# Where each of `ref`'s variants sits in `x`, and with which allele sign.
# @noRd
.ldmMatchTo <- function(x, ref) {
    m <- matchVariants(
        .ldRefMatchIds(ref),
        .ldRefMatchIds(x),
        removeStrandAmbiguous = FALSE
    )
    list(refIdx = m$idxA, idx = m$idxB, sign = m$sign)
}

# A component's positions and signs for the shared variants, in `ref` order.
# @noRd
.ldmAlignedOn <- function(matched, common) {
    at <- match(common, matched$refIdx)
    list(idx = matched$idx[at], sign = matched$sign[at])
}

# A component's matrix over the shared variants, re-signed into the first
# component's allele coding: a swapped variant's dosage counts the other
# allele, which negates each of its correlations.
# @noRd
.ldmAlignedR <- function(x, aligned) {
    R <- unname(x@correlation[aligned$idx, aligned$idx, drop = FALSE])
    if (all(aligned$sign > 0)) {
        return(R)
    }
    R * outer(aligned$sign, aligned$sign)
}

# Add the k-th component's weighted, aligned matrix to the running sum, so only
# one aligned copy exists at a time rather than every component's at once.
# @noRd
.ldmAddComponent <- function(acc, k, components, aligned, w) {
    term <- w[[k]] * .ldmAlignedR(components[[k]], aligned[[k]])
    if (is.null(acc)) term else acc + term
}

# The mixture's variant table: the first component's, with AF the weighted
# average of the components' (each put into the first component's coding)
# when every component records one.
# @noRd
.ldmMixtureVariants <- function(gr, components, aligned, w) {
    afs <- map2(components, aligned, .ldmAlignedAf)
    mc <- S4Vectors::mcols(gr)
    mc$AF <- if (any(map_lgl(afs, is.null))) {
        NULL
    } else {
        reduce(map2(w, afs, `*`), `+`)
    }
    S4Vectors::`mcols<-`(gr, value = mc)
}

# @noRd
.ldmAlignedAf <- function(x, aligned) {
    af <- S4Vectors::mcols(x)$AF
    if (is.null(af)) {
        return(NULL)
    }
    af <- as.numeric(af)[aligned$idx]
    if_else(aligned$sign > 0, af, 1 - af)
}

# 1 / sum(w^2 / n): the size of one panel whose LD is as noisy as the mixture.
# @noRd
.ldmMixtureN <- function(n, w) {
    if (anyNA(n)) {
        return(NA_integer_)
    }
    as.integer(round(1 / sum(w^2 / n)))
}

# @noRd
.ldmReportDropped <- function(nDropped, nKept, what = "LdMatrix") {
    if (nDropped == 0L) {
        return(invisible(NULL))
    }
    inform(glue(
        "{what}: {nDropped} variant(s) of the first component are not in ",
        "every component and are left out; the mixture covers {nKept}."
    ))
}

# ---- subsetting -------------------------------------------------------------

# Subsetting the variants subsets the matrix with them, so the two never fall
# out of step.
setMethod("extractROWS", "LdMatrix", function(x, i) {
    i <- normalizeSingleBracketSubscript(i, x)
    ans <- methods::callNextMethod(x, i)
    ans@correlation <- x@correlation[i, i, drop = FALSE]
    ans
})

# ---- accessors --------------------------------------------------------------

#' @rdname ldMatrix
#' @export
setMethod("ldMatrix", "LdMatrix", function(x) x@correlation)

#' @rdname ldMatrix
#' @importFrom S4Arrays as.matrix
#' @export
setMethod("as.matrix", "LdMatrix", function(x, ...) ldMatrix(x))

#' @rdname ldMatrix
#' @name coerce,LdMatrix,matrix-method
#' @aliases coerce,LdMatrix,matrix-method
#' @usage NULL
#' @exportMethod coerce
setAs("LdMatrix", "matrix", function(from) ldMatrix(from))

#' @rdname nSamples
#' @export
setMethod("nSamples", "LdMatrix", function(x, ...) x@nSamples)

#' @rdname mixtureWeights
#' @export
setMethod("mixtureWeights", "LdMatrix", function(x) {
    if (length(x@mixtureWeights) == 0L) NULL else x@mixtureWeights
})

# The canonical chrom:pos:A2:A1 id every class renders; the labels the matrix
# was given stay in mcols(x)$SNP.
#' @rdname variantIds
#' @export
setMethod("variantIds", "LdMatrix", function(x, ...) {
    .grVariantIds(methods::as(x, "GRanges"))
})

#' @rdname af
#' @export
setMethod("af", "LdMatrix", function(x, ...) {
    af <- S4Vectors::mcols(x)$AF
    ids <- as.character(S4Vectors::mcols(x)$SNP)
    set_names(if (is.null(af)) rep(NA_real_, length(ids)) else af, ids)
})

#' @rdname LdMatrix-class
#' @param object An \code{LdMatrix}.
#' @return \code{show} is called for its side effect and returns
#'   \code{invisible(NULL)}.
#' @export
setMethod("show", "LdMatrix", function(object) {
    n <- if (is.na(object@nSamples)) "unknown" else object@nSamples
    cat(glue(
        "LdMatrix: {length(object)} variants, panel size {n}\n",
        .trim = FALSE
    ))
    w <- object@mixtureWeights
    if (length(w) > 0L) {
        cat(glue(
            "  mixture of {length(w)} panels, weights ",
            "{str_flatten(signif(w, 3), ', ')}\n",
            .trim = FALSE
        ))
    }
    if (!is_in("AF", colnames(S4Vectors::mcols(object)))) {
        cat("  no allele frequencies\n")
    }
})
