#' @include MethodParam.R
NULL

#' @rdname QtlEnrichmentParam
#' @aliases QtlEnrichmentParam-class
#' @exportClass QtlEnrichmentParam
setClass(
    "QtlEnrichmentParam",
    contains = "MethodParam",
    slots = c(
        numGwas = "numeric_OR_NULL",
        piQtl = "numeric_OR_NULL",
        lambda = "numeric",
        impN = "numeric",
        doubleShrinkage = "logical",
        besselCorrection = "logical"
    )
)

#' @title Arguments For The QTL Enrichment Estimator
#' @description The algorithmic settings \code{\link{qtlEnrichment}} takes,
#'   as one checked bundle. Everything else that function accepts --- the two
#'   data arguments, \code{BPPARAM}, \code{seed}, \code{verbose}, and the
#'   \code{alignNames} switch --- is the pipeline's own business and is not
#'   settable here.
#' @param numGwas The total number of variants the GWAS tested (fastenloc's
#'   \code{-total_variants}), over which \code{piGwas} and \code{piQtl} are
#'   estimated. \code{NULL} (default) uses the number of distinct variants in
#'   the outcome PIP vector and the annotation fits, with a warning; that is
#'   biased unless the outcome PIP vector really is genome-wide.
#' @param piQtl Per-variant prior of being a causal QTL variant. \code{NULL}
#'   (default) estimates it as fastenloc does: the annotation fits' summed
#'   PIPs over the variant total (see \code{numGwas}).
#' @param lambda Shrinkage parameter, shrinking the enrichment estimate
#'   towards zero as in ridge regression. \code{0} applies no shrinkage.
#'   Default \code{1.0}.
#' @param impN Rounds of multiple imputation to draw QTL from. Default
#'   \code{25}.
#' @param doubleShrinkage Logical. Apply the double-shrinkage correction to
#'   the enrichment estimate. Default \code{FALSE}.
#' @param besselCorrection Logical. Apply Bessel's correction when estimating
#'   the sampling variance. Default \code{TRUE}.
#' @return A \code{QtlEnrichmentParam} object, a \code{\link{MethodParam}}.
#' @examples
#' QtlEnrichmentParam(lambda = 0, impN = 50)
#' @export
QtlEnrichmentParam <- function(
    numGwas = NULL,
    piQtl = NULL,
    lambda = 1.0,
    impN = 25,
    doubleShrinkage = FALSE,
    besselCorrection = TRUE
) {
    new(
        "QtlEnrichmentParam",
        numGwas = numGwas,
        piQtl = piQtl,
        lambda = lambda,
        impN = impN,
        doubleShrinkage = doubleShrinkage,
        besselCorrection = besselCorrection
    )
}

#' @title QTL Enrichment Pipeline (Genome-Wide)
#' @description Genome-wide pipeline that computes per-pair enrichment
#'   estimates by passing an outcome PIP vector and a set of annotation
#'   credible-set posteriors to \code{\link{qtlEnrichment}}. The returned table
#'   feeds \code{\link{colocPipeline}} via its \code{enrichment} argument.
#'
#'   Either side may be a \code{\link{QtlFineMappingResult}} or a
#'   \code{\link{GwasFineMappingResult}}, so QTL-in-GWAS, QTL-in-QTL and
#'   GWAS-in-GWAS enrichment all run through one code path. The argument names
#'   keep the QTL / GWAS wording of the common case; what they mean generally
#'   is that \code{gwasFineMappingResult} is the \strong{outcome} whose PIPs
#'   are scanned and \code{qtlFineMappingResult} is the \strong{annotation}
#'   whose region fits are tested for enrichment within them.
#'
#' \strong{Not gene-parallelisable}: the enrichment estimator runs over the full
#' genome of outcome PIPs and the full collection of annotation fits at once.
#'
#' @section Inputs:
#' \itemize{
#'   \item \code{gwasFineMappingResult}: the genome-wide outcome collection.
#'     Each row's stored fit (\code{\link{susieFit}}) must carry a \code{pip}
#'     vector. One PIP vector is built per outcome fit, across all its LD
#'     blocks -- keyed by (\code{study}, \code{method}) for a GWAS collection
#'     (one trait per study) and by (\code{study}, \code{context},
#'     \code{trait}, \code{method}) for a QTL one, since a variant's PIP
#'     differs between molecular traits and between fine-mapping methods, and
#'     pooling them would collide.
#'   \item \code{qtlFineMappingResult}: the annotation collection. Each row's
#'     stored fit must carry \code{alpha}, \code{pip}, and
#'     prior-variance fields (\code{V}). Its region fits are pooled per
#'     (\code{study}, \code{context}, \code{method}) -- per (\code{study},
#'     \code{method}) for a GWAS collection, which has no context axis -- so
#'     each method's fits make their own estimate rather than counting every
#'     region once per method.
#' }
#' Both sides' fits are read whole, so a row narrowed since its fit warns
#' (condition class \code{pecotmrFitCoverage}).
#'
#' @section LD reference compatibility check: Every row of a GWAS outcome
#'   collection must have an \code{ldReference} (it should be RSS-derived).
#'   For each QTL row and GWAS row whose regions overlap, where both carry
#'   one, the two must be the same LD reference; rows that never overlap may
#'   sit on different references.
#'   The two must be the same kind of reference. Two genotype panels must
#'   share their samples and allele orientation; two \code{LdMatrix} objects
#'   must agree on the LD they share; two \code{GenotypeMixture} objects must
#'   average the same panels with the same weights. They need NOT carry the
#'   same variants: QC-ing the two sides separately trims each reference to
#'   its own surviving variants, and a partial overlap only warns. A row with
#'   no reference on either side (an individual-level fit) skips the check.
#'
#' @param gwasFineMappingResult The outcome side; see above.
#' @param qtlFineMappingResult The annotation side; see above.
#' @param qtlEnrichmentParam The estimator's algorithmic settings, built with
#'   \code{\link{QtlEnrichmentParam}}: \code{numGwas}, \code{piQtl},
#'   \code{lambda}, \code{impN}, \code{doubleShrinkage} and
#'   \code{besselCorrection}. A bare list is refused, since it cannot be
#'   checked.
#' @param BPPARAM A \code{BiocParallelParam} handed to
#'   \code{\link{qtlEnrichment}}, whose multiple-imputation threads number
#'   its workers. Default \code{SerialParam()}.
#' @param seed Integer or \code{NULL}. Base random seed forwarded to
#'   \code{\link{qtlEnrichment}} for reproducible multiple imputation; every
#'   (outcome, annotation) pair is run with this same seed. \code{NULL}
#'   (default) draws a nondeterministic seed per pair.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param duplicatePips How to resolve a variant to which the outcome's LD
#'   blocks give different PIPs, as adjacent blocks sharing a boundary variant
#'   do when each block's fit gives it its own. Copies that agree are always
#'   merged. \code{"max"} (default) keeps the largest PIP. \code{"blockEnd"}
#'   keeps the copy from the block whose window was given (see
#'   \code{\link{genomicRegion}}) and ends at the variant --- the upstream of
#'   two blocks sharing it --- and errors when no single block qualifies.
#'   \code{"error"} refuses. fastenloc does none of these: it keeps the last
#'   copy's PIP but counts every copy towards the GWAS prior.
#' @return A tibble with one row per (outcome fit, annotation unit) pair.
#'   The identity columns are \code{gwasStudy}, \code{gwasContext},
#'   \code{gwasTrait}, \code{gwasMethod}, \code{qtlStudy}, \code{qtlContext},
#'   \code{qtlMethod}; the axes a side does not have are \code{NA}
#'   (\code{gwasContext} / \code{gwasTrait} for a GWAS outcome,
#'   \code{qtlContext} for a GWAS annotation). Suitable as the
#'   \code{enrichment} argument to \code{\link{colocPipeline}}, which joins on
#'   those columns.
#'
#'   The estimates are \code{enrichmentLogOdds}, the enrichment parameter
#'   \eqn{a_1} on the log-odds scale, with its standard error
#'   \code{enrichmentSe}; \code{enrichment} is the same quantity as a
#'   multiplicative factor, \eqn{e^{a_1} - 1}, reported for interpretation
#'   (so \eqn{a_1 = 0} means no enrichment).
#'   \code{enrichmentLogOddsNoShrinkage} and
#'   \code{enrichmentSeNoShrinkage} are the same estimate before shrinkage,
#'   \code{intercept} / \code{interceptSe} are \eqn{a_0}, and \code{colocP1},
#'   \code{colocP2}, \code{colocP12} are the enrichment-informed priors the
#'   estimator derives from \eqn{(a_0, a_1)} exactly as fastenloc's
#'   \code{set_enrich_params} does; those three are what
#'   \code{\link{colocPipeline}} scores with when given this table.
#'   \code{effectiveMiRounds} is how many multiple-imputation rounds survived
#'   outlier filtering.
#' @examples
#' data(gwasFineMappingExample)
#' data(qtlFineMappingExample)
#' qtlEnrichmentPipeline(
#'   gwasFineMappingResult = gwasFineMappingExample,
#'   qtlFineMappingResult = qtlFineMappingExample
#' )
#' @export
qtlEnrichmentPipeline <- function(
    gwasFineMappingResult,
    qtlFineMappingResult,
    qtlEnrichmentParam = QtlEnrichmentParam(),
    BPPARAM = SerialParam(),
    seed = NULL,
    verbose = TRUE,
    duplicatePips = c("max", "error", "blockEnd")
) {
    .assertMethodParam(qtlEnrichmentParam, "QtlEnrichmentParam")
    .assertBpparam(BPPARAM)
    duplicatePips <- arg_match(duplicatePips)
    .enrValidateInputs(gwasFineMappingResult, qtlFineMappingResult)
    prep <- .enrPrepare(
        gwasFineMappingResult,
        qtlFineMappingResult,
        duplicatePips
    )
    gwasTuples <- prep$gwasTuples
    qtlTuples <- prep$qtlTuples
    results <- list_flatten(map(
        seq_len(nrow(gwasTuples)),
        .enrScoreOutcomeTuple,
        prep = prep,
        qtlEnrichmentParam = qtlEnrichmentParam,
        BPPARAM = BPPARAM,
        seed = seed,
        verbose = verbose,
        gwasFineMappingResult = gwasFineMappingResult,
        gwasTuples = gwasTuples,
        qtlFineMappingResult = qtlFineMappingResult,
        qtlTuples = qtlTuples
    ))
    .enrAssemble(results)
}

# Validate the input classes + LD-sketch presence / identity.
# @noRd
#' @importFrom checkmate assertMultiClass
.enrValidateInputs <- function(gwasFineMappingResult, qtlFineMappingResult) {
    # assertMultiClass rather than assertClass on the virtual parent: naming
    # both concrete subclasses is more use to a caller than
    # "FineMappingResultBase", and these two are its only subclasses.
    fmrClasses <- c("GwasFineMappingResult", "QtlFineMappingResult")
    assertMultiClass(gwasFineMappingResult, fmrClasses)
    assertMultiClass(qtlFineMappingResult, fmrClasses)
    ldFree <- is.na(.ldRefKeys(gwasFineMappingResult))
    if (
        any(ldFree) &&
            methods::is(gwasFineMappingResult, "GwasFineMappingResult")
    ) {
        msg <- glue(
            "qtlEnrichmentPipeline: every GwasFineMappingResult row must have ",
            "an ldReference (it should be RSS-derived); {sum(ldFree)} of ",
            "{length(ldFree)} have none."
        )
        abort(msg)
    }
    # Per row pair, over the QTL and GWAS rows whose regions overlap -- the
    # ones the enrichment combines. Lenient rather than qtl-required: a QTL
    # outcome side may be an individual-level fit carrying no panel, and the
    # requirement above already covers the RSS-derived GWAS case.
    .requireOverlapLdRefs(
        qtlFineMappingResult,
        gwasFineMappingResult,
        pipelineName = "qtlEnrichmentPipeline",
        label = "qtlFineMappingResult",
        nullPolicy = "lenient"
    )
    invisible(NULL)
}

# Hoist the outcome-independent work out of the double loop: each outcome
# tuple's PIPs with its variants, the union of those variants, each annotation
# tuple's regions with theirs, and each tuple's one-time alignment to the union
# (errors captured as values). Variants travel as variant tables
# (.variantTable()) -- chromosome, position, alleles -- and are matched by the
# matching engine, never through rendered id strings.
# @noRd
.enrPrepare <- function(
    gwasFineMappingResult,
    qtlFineMappingResult,
    duplicatePips = "max"
) {
    gwasTuples <- .enrOutcomeTuples(gwasFineMappingResult)
    qtlTuples <- .enrAnnotationTuples(qtlFineMappingResult)
    if (nrow(gwasTuples) == 0L || nrow(qtlTuples) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no (outcome, annotation) pairs to ",
            "compute (one of the inputs has zero rows)."
        )
        abort(msg)
    }
    # Each row's fit is read whole, so a row narrowed since its fit warns.
    .fitCoverageCheckedAll(
        list(gwasFineMappingResult, qtlFineMappingResult),
        "qtlEnrichmentPipeline()"
    )
    outcomes <- map(
        seq_len(nrow(gwasTuples)),
        .enrGwasOutcomeForRow,
        gwasTuples = gwasTuples,
        fmr = gwasFineMappingResult,
        duplicatePips = duplicatePips
    )
    union <- .enrUnionPanel(outcomes)
    annotations <- map(
        seq_len(nrow(qtlTuples)),
        .enrQtlTupleForRow,
        qtlTuples = qtlTuples,
        fmr = qtlFineMappingResult
    )
    list(
        gwasTuples = gwasTuples,
        qtlTuples = qtlTuples,
        outcomes = outcomes,
        union = union,
        annotations = annotations,
        aligned = map(annotations, .enrAlignTupleSafe, panel = union$table)
    )
}

# The outcome side's per-fit keys: one PIP vector is built per key. A GWAS
# collection carries one trait per study, so its QTL-only axes read NA and the
# key collapses to (study, method). The method is part of it: two methods'
# fits of one trait are two PIP vectors, and pooling them would put two
# different PIPs on every variant.
# @noRd
.enrOutcomeTuples <- function(fmr) {
    distinct(tibble(
        gwasStudy = .fmrIdentityColumn(fmr, "study"),
        gwasContext = .fmrIdentityColumn(fmr, "context"),
        gwasTrait = .fmrIdentityColumn(fmr, "trait"),
        gwasMethod = .fmrIdentityColumn(fmr, "method")
    ))
}

# The annotation side's keys: (study, context, method), the unit whose region
# fits are pooled into one enrichment estimate. The joint key matters --
# context alone would silently merge two studies sharing a context label, and
# without the method each region would count once per method that fit it --
# and a GWAS annotation pools its blocks under one NA-context key.
# @noRd
.enrAnnotationTuples <- function(fmr) {
    distinct(tibble(
        qtlStudy = .fmrIdentityColumn(fmr, "study"),
        qtlContext = .fmrIdentityColumn(fmr, "context"),
        qtlMethod = .fmrIdentityColumn(fmr, "method")
    ))
}

# One tuple-table row as an identity list addressed by the collection's OWN
# column names; the table prefixes them so the result can name both sides.
# @noRd
.enrOutcomeIdent <- function(gwasTuples, k) {
    list(
        study = gwasTuples$gwasStudy[[k]],
        context = gwasTuples$gwasContext[[k]],
        trait = gwasTuples$gwasTrait[[k]],
        method = gwasTuples$gwasMethod[[k]]
    )
}

# @noRd
.enrAnnotationIdent <- function(qtlTuples, k) {
    list(
        study = qtlTuples$qtlStudy[[k]],
        context = qtlTuples$qtlContext[[k]],
        method = qtlTuples$qtlMethod[[k]]
    )
}

# Row indices of `fmr` matching an identity tuple. An axis the collection does
# not have is NA on both sides and matches, rather than excluding every row.
# @noRd
.enrMatchRows <- function(fmr, ident) {
    hits <- map(names(ident), .enrColumnMatches, fmr = fmr, ident = ident)
    which(reduce(hits, `&`, .init = rep(TRUE, nrow(fmr))))
}

# @noRd
.enrColumnMatches <- function(column, fmr, ident) {
    values <- .fmrIdentityColumn(fmr, column)
    wanted <- ident[[column]]
    if (is.na(wanted)) {
        return(is.na(values))
    }
    !is.na(values) & values == wanted
}

# The outcome PIPs and variants for the k-th outcome tuple.
# @noRd
.enrGwasOutcomeForRow <- function(k, gwasTuples, fmr, duplicatePips = "max") {
    .enrGwasOutcome(fmr, .enrOutcomeIdent(gwasTuples, k), duplicatePips)
}

# Annotation regions and variants for the k-th (study, context) tuple.
# @noRd
.enrQtlTupleForRow <- function(k, qtlTuples, fmr) {
    .enrQtlTuple(fmr, .enrAnnotationIdent(qtlTuples, k))
}

# One tuple's variants aligned to the union of the outcome variants: each
# variant's index in the union, or NA. Any error is captured as a value
# (re-raised + skipped per (gwas, tuple) below, never aborting).
# @noRd
#' @importFrom rlang try_fetch
.enrAlignTupleSafe <- function(annotation, panel) {
    try_fetch(
        .alignTablesToPanel(annotation$table, panel),
        error = function(cnd) cnd
    )
}

# Score one outcome trait against every annotation tuple -> enrichment records
# (empty when the outcome has no usable PIPs).
# @noRd
.enrScoreOutcomeTuple <- function(
    gi,
    prep,
    qtlEnrichmentParam,
    BPPARAM,
    seed,
    verbose,
    gwasFineMappingResult = gwasFineMappingResult,
    gwasTuples = gwasTuples,
    qtlFineMappingResult = qtlFineMappingResult,
    qtlTuples = qtlTuples
) {
    outcome <- prep$outcomes[[gi]]
    if (length(outcome$pip) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable PIPs for ",
            "{.enrOutcomeLabel(gwasFineMappingResult, gwasTuples, gi)}; ",
            "skipping."
        )
        warn(msg)
        return(list())
    }
    # Where each union variant sits in this outcome's PIP vector.
    outcomeAt <- rep(NA_integer_, .variantTableSize(prep$union$table))
    outcomeAt[prep$union$of[[gi]]] <- seq_along(outcome$pip)
    compact(map(
        seq_len(nrow(qtlTuples)),
        .enrScoreTuple,
        gi = gi,
        outcome = outcome,
        outcomeAt = outcomeAt,
        prep = prep,
        qtlEnrichmentParam = qtlEnrichmentParam,
        BPPARAM = BPPARAM,
        seed = seed,
        verbose = verbose,
        gwasFineMappingResult = gwasFineMappingResult,
        gwasTuples = gwasTuples,
        qtlFineMappingResult = qtlFineMappingResult,
        qtlTuples = qtlTuples
    ))
}

# Score one (outcome trait, annotation tuple) pair -> an enrichment record, or
# NULL when the tuple has no regions or the enrichment fails.
# @noRd
.enrScoreTuple <- function(
    k,
    gi,
    outcome,
    outcomeAt,
    prep,
    qtlEnrichmentParam,
    BPPARAM,
    seed,
    verbose,
    gwasFineMappingResult = gwasFineMappingResult,
    gwasTuples = gwasTuples,
    qtlFineMappingResult = qtlFineMappingResult,
    qtlTuples = qtlTuples
) {
    if (length(prep$annotations[[k]]$regions) == 0L) {
        msg <- glue(
            "qtlEnrichmentPipeline: no usable regions for ",
            "{.enrAnnotationLabel(qtlFineMappingResult, qtlTuples, k)}; ",
            "skipping."
        )
        warn(msg)
        return(NULL)
    }
    enr <- .enrRunEnrichment(
        gi,
        outcome,
        outcomeAt,
        k,
        prep = prep,
        qtlEnrichmentParam = qtlEnrichmentParam,
        BPPARAM = BPPARAM,
        seed = seed,
        verbose = verbose,
        gwasFineMappingResult = gwasFineMappingResult,
        gwasTuples = gwasTuples,
        qtlFineMappingResult = qtlFineMappingResult,
        qtlTuples = qtlTuples
    )
    if (is.null(enr)) {
        return(NULL)
    }
    c(
        .enrFlattenEnrichment(enr),
        as.list(gwasTuples[gi, , drop = FALSE]),
        as.list(qtlTuples[k, , drop = FALSE])
    )
}

# Human-readable identities for the warnings above, naming each side by its own
# flavour and only the axes it has.
# @noRd
.enrOutcomeLabel <- function(gwasFineMappingResult, gwasTuples, gi) {
    .fmrTupleLabel(
        .fmrSideName(gwasFineMappingResult),
        .enrOutcomeIdent(gwasTuples, gi)
    )
}

# @noRd
.enrAnnotationLabel <- function(qtlFineMappingResult, qtlTuples, k) {
    .fmrTupleLabel(
        .fmrSideName(qtlFineMappingResult),
        .enrAnnotationIdent(qtlTuples, k)
    )
}

# Run the enrichment for a pair -- the tuple's regions coded against this
# outcome through the precomputed alignment -- warning + NULL on failure.
# @noRd
.enrRunEnrichment <- function(
    gi,
    outcome,
    outcomeAt,
    k,
    prep,
    qtlEnrichmentParam,
    BPPARAM,
    seed,
    verbose,
    gwasFineMappingResult,
    gwasTuples,
    qtlFineMappingResult,
    qtlTuples
) {
    aligned <- prep$aligned[[k]]
    annotation <- prep$annotations[[k]]
    try_fetch(
        {
            if (inherits(aligned, "condition")) {
                cnd_signal(aligned)
            }
            nullFits <- .enrValidateFits(annotation$regions)
            coded <- .enrCodeAligned(
                annotation,
                outcomeAt[aligned],
                aligned,
                prep$union$table
            )
            exec(
                .qtlEnrichmentCoded,
                gwasPip = outcome$pip,
                susieQtlRegions = annotation$regions,
                coded = coded,
                nullFits = nullFits,
                gwasVariants = outcome$table,
                BPPARAM = BPPARAM,
                seed = seed,
                verbose = verbose,
                !!!qtlEnrichmentParam
            )
        },
        error = function(cnd) {
            eMsg <- conditionMessage(cnd)
            msg <- glue(
                "qtlEnrichmentPipeline: qtlEnrichment failed for ",
                "{.enrOutcomeLabel(gwasFineMappingResult, gwasTuples, gi)} x ",
                "{.enrAnnotationLabel(qtlFineMappingResult, qtlTuples, k)}: ",
                "{eMsg}"
            )
            warn(msg)
            NULL
        }
    )
}

# Row-bind the enrichment records into the id-first result table.
# @noRd
.enrAssemble <- function(results) {
    if (length(results) == 0L) {
        return(.enrEmptyResult())
    }
    out <- bind_rows(results)
    select(out, all_of(.enrIdCols()), everything())
}

# The identity columns of a result row: the outcome trait's tuple, then the
# annotation unit's.
# @noRd
.enrIdCols <- function() {
    c(
        "gwasStudy",
        "gwasContext",
        "gwasTrait",
        "gwasMethod",
        "qtlStudy",
        "qtlContext",
        "qtlMethod"
    )
}

# The empty enrichment result table, with the same columns a populated one has.
# @noRd
.enrEmptyResult <- function() {
    idCols <- set_names(
        rep(list(character(0)), length(.enrIdCols())),
        .enrIdCols()
    )
    valueNames <- names(.enrNaEnrichment())
    valueCols <- set_names(
        rep(list(numeric(0)), length(valueNames)),
        valueNames
    )
    as_tibble(c(idCols, valueCols))
}

# =============================================================================
# Internal helpers
# =============================================================================

# Align one tuple's QTL regions to the GWAS naming convention: relabel each
# pip name that matches a panel variant to the panel's id, keeping the rest.
# The tuple's names are matched in one pass, not one per region, so the panel
# is read once per tuple however many genes it has. Used by qtlEnrichment(),
# whose interface is named vectors; the pipeline aligns variant tables
# (.alignTablesToPanel()).
# @noRd
.enrAlignRegions <- function(
    regions,
    panel,
    panelTable = .lazyVariantTable(panel)
) {
    ids <- map(regions, .enrRegionIds)
    flat <- as.character(unlist(ids, use.names = FALSE))
    idx <- .alignIdsToPanel(flat, panel, panelTable)
    matched <- !is.na(idx)
    aligned <- replace(flat, matched, panel[idx[matched]])
    regionOf <- factor(
        rep(seq_along(regions), lengths(ids)),
        levels = seq_along(regions)
    )
    map2(regions, unname(split(aligned, regionOf)), .enrRenamePip)
}

# A region with its pip renamed (an unnamed pip has nothing to rename).
# @noRd
.enrRenamePip <- function(x, ids) {
    if (is.null(names(x$pip))) {
        return(x)
    }
    list_assign(x, pip = set_names(x$pip, ids))
}

# Rows `idx`'s PIPs with their variants, read for all of them at once:
# list(pip, table, row), `row` being each PIP's row of `gwasFmr`. A row that
# carries no fit, no PIPs, or variants that do not line up with them gives
# none. The variants are the names the fit gives its PIPs, else the row's own
# when it holds all of the fit (.fmrRowsFitVariants()).
# @noRd
.enrPipTablesAt <- function(gwasFmr, idx) {
    fits <- .fmrFitsAt(gwasFmr, idx)
    usable <- which(map_lgl(fits, .enrHasPip))
    vars <- .fmrRowsFitVariants(
        gwasFmr,
        idx[usable],
        fits[usable],
        unplaced = function(i) NULL
    )
    pips <- map(fits[usable], function(fit) as.numeric(fit$pip))
    n <- vars$n
    lined <- !is.na(n) & n == lengths(pips)
    held <- rep(lined[!is.na(n)], n[!is.na(n)])
    list(
        pip = unlist(pips[lined], use.names = FALSE),
        table = .variantTableRows(vars$table, which(held)),
        row = rep(idx[usable][lined], n[lined])
    )
}

# @noRd
.enrHasPip <- function(fit) {
    is.list(fit) && !is.null(fit$pip)
}

# One row's PIP vector keyed by variant id (see .enrPipTablesAt()).
# @noRd
.enrRowPipVector <- function(i, gwasFmr) {
    piece <- .enrPipTablesAt(gwasFmr, i)
    if (length(piece$pip) == 0L) {
        return(NULL)
    }
    set_names(piece$pip, .variantTableIds(piece$table))
}

# One outcome fit's PIPs and variants: every row of the collection carrying
# that identity -- the fit's LD blocks -- joined. A variant several rows carry
# is resolved by `duplicatePips` -- which is why the identity is the full fit
# tuple: pooling two molecular traits of one study, or two methods' fits of
# one trait, would collide on every variant they share. list(pip, table).
# @noRd
.enrGwasOutcome <- function(gwasFmr, ident, duplicatePips = "max") {
    idx <- .enrMatchRows(gwasFmr, ident)
    pieces <- .enrPipTablesAt(gwasFmr, idx)
    if (length(pieces$pip) == 0L) {
        return(list(
            pip = numeric(0),
            table = .variantTableRows(pieces$table, 0L)
        ))
    }
    .enrCollapseDuplicateVariants(pieces, duplicatePips, gwasFmr)
}

# The outcome PIP vector of one fit keyed by variant id (see
# .enrGwasOutcome()).
# @noRd
.enrBuildGwasPipVector <- function(gwasFmr, ident, duplicatePips = "max") {
    outcome <- .enrGwasOutcome(gwasFmr, ident, duplicatePips)
    if (length(outcome$pip) == 0L) {
        return(numeric(0))
    }
    set_names(outcome$pip, .variantTableIds(outcome$table))
}

# Collapse variants several GWAS blocks carry, keeping each where it first
# appears. Copies whose PIPs agree (to 12 digits) merge; for a variant whose
# copies conflict, `duplicatePips` decides (see qtlEnrichmentPipeline()).
# A variant is its chromosome, position and alleles as given.
#' @importFrom vctrs vec_group_id vec_unique_loc
#' @noRd
.enrCollapseDuplicateVariants <- function(pieces, duplicatePips, gwasFmr) {
    group <- vec_group_id(.variantKeyFrame(pieces$table))
    first <- which(!duplicated(group))
    if (length(first) == length(group)) {
        return(list(pip = pieces$pip, table = pieces$table))
    }
    out <- list(
        pip = pieces$pip[first],
        table = .variantTableRows(pieces$table, first)
    )
    copied <- which(is_in(group, group[duplicated(group)]))
    distinctPips <- copied[vec_unique_loc(data.frame(
        group = group[copied],
        pip = round(pieces$pip[copied], 12)
    ))]
    nPips <- tabulate(group[distinctPips], nbins = attr(group, "n"))
    conflicting <- which(nPips > 1L)
    if (length(conflicting) == 0L) {
        return(out)
    }
    copies <- copied[is_in(group[copied], conflicting)]
    resolved <- switch(
        duplicatePips,
        error = .enrAbortConflictingPips(
            .variantTableIds(.variantTableRows(pieces$table, copies[[1L]]))
        ),
        max = .enrMaxPips(pieces, copies, group),
        blockEnd = .enrBlockEndPips(pieces, copies, group, gwasFmr)
    )
    # `out` holds each group once, in the order of first appearance, which
    # is the order of the group ids.
    out$pip[as.integer(names(resolved))] <- resolved
    out
}

# @noRd
.enrAbortConflictingPips <- function(id) {
    msg <- glue(
        "qtlEnrichmentPipeline: variant '{id}' appears with conflicting PIPs ",
        "across GWAS blocks; the GWAS fine-mapping must produce a consistent ",
        "PIP per variant, or pass duplicatePips = \"max\" or \"blockEnd\"."
    )
    abort(msg)
}

# The largest PIP of each conflicting variant, named by its group.
# @noRd
.enrMaxPips <- function(pieces, copies, group) {
    vapply(split(pieces$pip[copies], group[copies]), max, numeric(1))
}

# The PIP of each conflicting variant's copy from the one block whose given
# window ends at the variant, named by its group; a variant with no such
# block, or several, aborts.
#' @importFrom GenomicRanges end seqnames
#' @noRd
.enrBlockEndPips <- function(pieces, copies, group, gwasFmr) {
    rows <- pieces$row[copies]
    windows <- genomicRegion(gwasFmr)[rows]
    given <- as.logical(.tupleColumn(gwasFmr, "regionGiven"))[rows]
    atEnd <- !is.na(given) &
        given &
        as.character(seqnames(windows)) == pieces$table$chrom[copies] &
        end(windows) == pieces$table$pos[copies]
    nAtEnd <- vapply(split(atEnd, group[copies]), sum, integer(1))
    undecided <- names(nAtEnd)[nAtEnd != 1L]
    if (length(undecided) > 0L) {
        which1 <- copies[group[copies] == as.integer(undecided[[1L]])][[1L]]
        id <- .variantTableIds(.variantTableRows(pieces$table, which1))
        msg <- glue(
            "qtlEnrichmentPipeline: variant '{id}' appears with ",
            "conflicting PIPs across GWAS blocks, and duplicatePips = ",
            "\"blockEnd\" cannot choose between them: ",
            "{nAtEnd[[undecided[[1L]]]]} of its blocks were given a window ",
            "ending at it, where one must be."
        )
        abort(msg)
    }
    set_names(pieces$pip[copies][atEnd], group[copies][atEnd])
}

# The union of the outcome tuples' variants, and each tuple's variants' index
# in it: list(table, of).
#' @importFrom vctrs vec_match vec_unique_loc
#' @noRd
.enrUnionPanel <- function(outcomes) {
    tables <- map(outcomes, "table")
    if (length(tables) == 1L) {
        n <- .variantTableSize(tables[[1L]])
        return(list(table = tables[[1L]], of = list(seq_len(n))))
    }
    all <- .bindVariantTables(tables)
    frame <- .variantKeyFrame(all)
    union <- .variantTableRows(all, vec_unique_loc(frame))
    unionFrame <- .variantKeyFrame(union)
    list(
        table = union,
        of = map(tables, function(t) vec_match(.variantKeyFrame(t), unionFrame))
    )
}

# One annotation tuple's region fits for qtlEnrichment -- alpha, pip and the
# prior variances -- and their variants as one variant table, region after
# region: list(regions, table, n). Filters on BOTH study and context so
# entries from different studies that happen to share a context label are not
# pooled into one enrichment estimate. Regions are named by their rows' own
# names ("study|context|trait|method"), so a fit the estimator refuses is
# reported by name rather than position.
# @noRd
.enrQtlTuple <- function(qtlFmr, ident) {
    idx <- .enrMatchRows(qtlFmr, ident)
    fits <- .fmrFitsAt(qtlFmr, idx)
    usable <- which(map_lgl(fits, .enrHasRegion))
    vars <- .fmrRowsFitVariants(
        qtlFmr,
        idx[usable],
        fits[usable],
        unplaced = function(i) .partsFitIds(.fmrRowParts(qtlFmr, i))
    )
    regions <- map(fits[usable], .enrRegionOf)
    names(regions) <- names(qtlFmr)[idx[usable]]
    list(regions = regions, table = vars$table, n = vars$n)
}

# Rows `idx`'s region payloads for qtlEnrichment with their PIPs named by
# variant id; NULL for a row that lacks a fit, an alpha, PIPs, or a prior
# variance (see .enrQtlTuple()).
# @noRd
.enrRegionsAt <- function(qtlFmr, idx) {
    fits <- .fmrFitsAt(qtlFmr, idx)
    usable <- which(map_lgl(fits, .enrHasRegion))
    out <- vector("list", length(idx))
    if (length(usable) == 0L) {
        return(out)
    }
    vars <- .fmrRowsFitVariants(
        qtlFmr,
        idx[usable],
        fits[usable],
        unplaced = function(i) .partsFitIds(.fmrRowParts(qtlFmr, i))
    )
    ids <- split(
        .variantTableIds(vars$table),
        factor(rep(seq_along(usable), vars$n), levels = seq_along(usable))
    )
    out[usable] <- map2(
        map(fits[usable], .enrRegionOf),
        unname(ids),
        function(region, ids) {
            list_assign(region, pip = set_names(region$pip, ids))
        }
    )
    out
}

# One row's region payload for qtlEnrichment (see .enrRegionsAt()).
# @noRd
.enrRowRegion <- function(i, qtlFmr) {
    .enrRegionsAt(qtlFmr, i)[[1L]]
}

# The per-(study, context) list of region fits with named PIPs, as
# qtlEnrichment() takes them (see .enrQtlTuple()).
# @noRd
.enrBuildQtlRegionsList <- function(qtlFmr, ident) {
    idx <- .enrMatchRows(qtlFmr, ident)
    if (length(idx) == 0L) {
        return(list())
    }
    regions <- .enrRegionsAt(qtlFmr, idx)
    names(regions) <- names(qtlFmr)[idx]
    compact(regions)
}

# Whether `fit` carries what a region payload needs.
# @noRd
.enrHasRegion <- function(fit) {
    is.list(fit) &&
        !is.null(fit$alpha) &&
        !is.null(fit$pip) &&
        !is.null(fit$V %||% fit$prior_variance)
}

# @noRd
.enrRegionOf <- function(fit) {
    list(
        alpha = fit$alpha,
        pip = fit$pip,
        prior_variance = fit$V %||% fit$prior_variance
    )
}

# An annotation tuple coded against one outcome: each region variable's
# 0-based index in the outcome's PIP vector (`outcomeIdx`, through the union
# alignment `unionIdx`), or a negative id shared by every region carrying the
# same absent variant -- the union's variant where the union holds it, else
# the tuple's own. Coded as .enrCodeRegions() codes names.
#' @importFrom vctrs vec_group_id
#' @noRd
.enrCodeAligned <- function(annotation, outcomeIdx, unionIdx, unionTable) {
    absent <- which(is.na(outcomeIdx))
    codes <- outcomeIdx - 1L
    nQtlOnly <- 0L
    if (length(absent) > 0L) {
        own <- .variantTableRows(annotation$table, absent)
        inUnion <- !is.na(unionIdx[absent])
        key <- own
        if (any(inUnion)) {
            fromUnion <- .variantTableRows(
                unionTable,
                unionIdx[absent][inUnion]
            )
            for (nm in names(key)) {
                key[[nm]][inUnion] <- fromUnion[[nm]]
            }
        }
        group <- vec_group_id(.variantKeyFrame(key))
        codes[absent] <- -as.vector(group)
        nQtlOnly <- attr(group, "n")
    }
    regionOf <- factor(
        rep(seq_along(annotation$regions), annotation$n),
        levels = seq_along(annotation$regions)
    )
    list(
        fits = map2(
            annotation$regions,
            unname(split(codes, regionOf)),
            .enrCodedFit
        ),
        unmatched = NULL,
        nQtlOnly = nQtlOnly
    )
}

# Pull one enrichment field from qtlEnrichment's list output as a scalar numeric
# (NA when absent).
# @noRd
.enrPickScalar <- function(field, enr) {
    v <- enr[[field]]
    if (is.null(v)) {
        NA_real_
    } else {
        as.numeric(v[[1L]])
    }
}

# Project qtlEnrichment()'s output onto the columns this pipeline publishes.
#
# The field names are the estimator's own, written in src/qtl_enrichment.h and
# shared verbatim with the original fastenloc's enloc.enrich.out, so they are
# matched literally rather than guessed at. The shrinkage estimates are the
# ones reported, matching the original fastenloc, which likewise switches to
# the shrunk a1 before deriving the coloc priors.
#
# `enrichment` is expm1 of the log-odds a1 rather than a1 itself: the
# relative change in an annotated variant's odds of being a QTL, 0 for an
# unenriched annotation and bounded below by -1. colocPipeline reports it
# beside each pair; its priors are the coloc p1 / p2 / p12 the estimator
# derives (colocP1, colocP2, colocP12), not a rescaling of this column.
# @noRd
.enrFlattenEnrichment <- function(enr) {
    if (!is.list(enr) || !is_in("Enrichment (w/ shrinkage)", names(enr))) {
        msg <- glue(
            "qtlEnrichmentPipeline: the enrichment estimator returned no ",
            "'Enrichment (w/ shrinkage)' field, so every estimate is ",
            "reported as NA. Expected the field names written by ",
            "qtl_enrichment.h."
        )
        warn(msg)
        return(.enrNaEnrichment())
    }
    logOdds <- .enrPickScalar("Enrichment (w/ shrinkage)", enr)
    list(
        enrichment = expm1(logOdds),
        enrichmentSe = .enrPickScalar("sd (w/ shrinkage)", enr),
        enrichmentLogOdds = logOdds,
        enrichmentLogOddsNoShrinkage = .enrPickScalar(
            "Enrichment (no shrinkage)",
            enr
        ),
        enrichmentSeNoShrinkage = .enrPickScalar("sd (no shrinkage)", enr),
        intercept = .enrPickScalar("Intercept", enr),
        interceptSe = .enrPickScalar("sd (intercept)", enr),
        colocP1 = .enrPickScalar("Alternative (coloc) p1", enr),
        colocP2 = .enrPickScalar("Alternative (coloc) p2", enr),
        colocP12 = .enrPickScalar("Alternative (coloc) p12", enr),
        effectiveMiRounds = .enrPickScalar("Effective MI rounds", enr)
    )
}

# The value columns, all unmeasured. Also the single source of the value-column
# schema, so the empty result cannot drift from the populated one.
# @noRd
.enrNaEnrichment <- function() {
    list(
        enrichment = NA_real_,
        enrichmentSe = NA_real_,
        enrichmentLogOdds = NA_real_,
        enrichmentLogOddsNoShrinkage = NA_real_,
        enrichmentSeNoShrinkage = NA_real_,
        intercept = NA_real_,
        interceptSe = NA_real_,
        colocP1 = NA_real_,
        colocP2 = NA_real_,
        colocP12 = NA_real_,
        effectiveMiRounds = NA_real_
    )
}


# =============================================================================
# qtlEnrichment: low-level enrichment estimation
# -----------------------------------------------------------------------------
# Per-(GWAS, QTL-region-list) enrichment estimator. Called per-(gwasStudy,
# qtlContext) pair by qtlEnrichmentPipeline above. Uses the fastenloc-style
# C++ kernel (qtlEnrichmentCpp) under the hood.
# =============================================================================
#' @title Implementation of enrichment analysis described in
#'   https://doi.org/10.1371/journal.pgen.1006646
#'
#' @description Largely follows from fastenloc
#'   https://github.com/xqwen/fastenloc but uses `susieR` fitted objects as
#'   input to estimate prior for use with `coloc` package (coloc v5, aka
#'   SuSiE-coloc). The main differences are 1) now enrichment is based on all
#'   QTL variants whether or not they are inside signal clusters; 2) Causal QTL
#'   are sampled from SuSiE single effects, not signal clusters; 3) Allow a
#'   variant to be QTL for not only multiple conditions (eg cell types) but also
#'   multiple regions (eg genes). Other minor improvements include 1) Make GSL
#'   RNG thread-safe; 2) Release memory from QTL binary annotation samples
#'   immediately after they are used.
#' @details Uses output of \code{\link[susieR]{susie}} from the \code{susieR}
#'   package.
#'
#'   As in fastenloc, the variants the estimate is taken over are every
#'   variant of \code{gwasPip} plus every QTL variant \code{gwasPip} does not
#'   carry. The latter have GWAS PIP 0, so each one sampled in an imputation
#'   round counts as a QTL variant that is not causal for the GWAS trait. The
#'   variant total is \code{numGwas}, raised to the number of those variants
#'   when smaller; \code{piGwas} is \code{sum(gwasPip)} over it and, unless
#'   supplied, \code{piQtl} is the fits' summed PIPs over it. The intercept is
#'   derived from \code{piGwas}, \code{piQtl} and the unshrunk enrichment,
#'   and the coloc priors from that intercept and the shrunk enrichment, both
#'   as fastenloc does. One difference remains: fastenloc's EM keeps the
#'   previous iteration's intercept at convergence where this one updates it,
#'   which changes only \code{sd (intercept)}.
#'
#' @param gwasPip A named vector of GWAS PIPs, one per variant, ideally
#'   genome-wide.
#' @param susieQtlRegions A named list of SuSiE fits, one per QTL unit
#'   analyzed (e.g. a gene), each with \code{alpha}, a named \code{pip} and
#'   \code{prior_variance}.
#' @param numGwas The total number of variants the GWAS tested (fastenloc's
#'   \code{-total_variants}). It matters whenever \code{gwasPip} does not
#'   carry every variant tested, e.g. when it covers only fine-mapped blocks.
#'   \code{NULL} (default) uses the number of distinct variants in
#'   \code{gwasPip} and the QTL fits, with a warning; a value below that number
#'   is raised to it.
#' @param piQtl Per-variant prior of being a causal QTL variant. \code{NULL}
#'   (default) estimates it as fastenloc does, from the fits' summed PIPs over
#'   the variant total.
#' @param lambda Similar to the shrinkage parameter used in ridge regression. It
#'   takes any non-negative value and shrinks the enrichment estimate towards 0.
#'   When it is set to 0, no shrinkage will be applied. A large value indicates
#'   strong shrinkage. The default value is set to 1.0.
#' @param impN Rounds of multiple imputation to draw QTL from, default is 25.
#' @param BPPARAM A \code{BiocParallelParam} whose workers number the threads
#'   the multiple imputation runs on when they run on this machine
#'   (\code{SerialParam}, \code{MulticoreParam}); any other backend gives one
#'   thread. Default \code{SerialParam()}.
#' @param alignNames Logical; when TRUE (default) a QTL pip name not spelled
#'   as any \code{gwasPip} name is matched, allele-aware (id format, ref/alt
#'   swap, strand), to the GWAS variants at its position, and takes the name
#'   of the one it matches. Set FALSE when the caller has already aligned them
#'   (e.g. \code{qtlEnrichmentPipeline} aligns each QTL tuple once against
#'   the union GWAS panel rather than re-aligning per GWAS study); the names
#'   are then matched to \code{gwasPip} exactly as given.
#' @param doubleShrinkage Logical. Apply the double-shrinkage correction to the
#'   enrichment estimate. Default \code{FALSE}.
#' @param besselCorrection Logical. Apply Bessel's correction when estimating
#'   the sampling variance. Default \code{TRUE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param seed Integer or \code{NULL}. Base random seed for the multiple-
#'   imputation sampler; each imputation round derives its own seed from it, so
#'   a fixed \code{seed} gives reproducible results. \code{NULL} (default) draws
#'   a nondeterministic seed.
#' @return A named list of enrichment parameter estimates, carrying the
#'   fields the C++ estimator writes -- \code{Intercept},
#'   \code{sd (intercept)}, \code{Enrichment (no shrinkage)},
#'   \code{Enrichment (w/ shrinkage)}, \code{sd (no shrinkage)},
#'   \code{sd (w/ shrinkage)}, \code{Alternative (coloc) p1} / \code{p2} /
#'   \code{p12} and \code{Effective MI rounds} -- plus diagnostics:
#'   \code{Missing proportion per MI round}, the share of each round's sampled
#'   QTL variants that \code{gwasPip} does not carry (they enter with GWAS
#'   PIP 0);
#'   \code{EM converged per MI round}; \code{EM iterations per MI round},
#'   how many EM iterations each round took; \code{Null QTL fits}, the number of
#'   fits with no effect of positive prior variance, which contribute no
#'   sampled variants; and \code{unused_xqtl_variants}, the QTL variants of
#'   each region that no GWAS variant matched. The estimate names are the
#'   original fastenloc's; \code{\link{qtlEnrichmentPipeline}} is what renames
#'   them to a tidy table.
#'
#' @examples
#' set.seed(1)
#' # GWAS PIPs for 2000 variants of a GWAS that tested one million: mostly
#' # near zero, with two strong signals.
#' ids <- sprintf("chr1:%d:A:G", 1000 * seq_len(2000))
#' gwasPip <- stats::setNames(stats::rbeta(2000, 0.05, 100), ids)
#' gwasPip[c(300, 1200)] <- c(0.95, 0.6)
#' # A SuSiE-like fit per gene: two single effects over the gene's variants,
#' # each concentrated on a few of them.
#' simulateFit <- function(idx, hits) {
#'   alpha <- matrix(1e-4, length(hits), length(idx))
#'   alpha[cbind(seq_along(hits), hits)] <- 1
#'   alpha <- alpha / rowSums(alpha)
#'   list(
#'     alpha = alpha,
#'     pip = stats::setNames(1 - apply(1 - alpha, 2, prod), ids[idx]),
#'     prior_variance = c(0.4, 0.1)
#'   )
#' }
#' susieFits <- list(
#'   gene1 = simulateFit(251:350, c(50, 10)),
#'   gene2 = simulateFit(1151:1250, c(50, 80)),
#'   gene3 = simulateFit(1601:1700, c(20, 70))
#' )
#' en <- qtlEnrichment(
#'   gwasPip, susieFits,
#'   numGwas = 1e6, impN = 10, seed = 1,
#'   BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
#' )
#' en[["Enrichment (w/ shrinkage)"]]
#'
#' @seealso \code{\link[susieR]{susie}}
#' @useDynLib pecotmr, .registration = TRUE
#' @export
#'
qtlEnrichment <- function(
    gwasPip,
    susieQtlRegions,
    numGwas = NULL,
    piQtl = NULL,
    lambda = 1.0,
    impN = 25,
    doubleShrinkage = FALSE,
    besselCorrection = TRUE,
    BPPARAM = SerialParam(),
    verbose = TRUE,
    alignNames = TRUE,
    seed = NULL
) {
    .assertBpparam(BPPARAM)
    .enrValidateNames(gwasPip, susieQtlRegions)
    nullFits <- .enrValidateFits(susieQtlRegions)
    # Relabel each region's pip names to the GWAS convention, then code every
    # variable against the GWAS vector, so the kernel works on indices.
    aligned <- .enrAlignPipNames(susieQtlRegions, gwasPip, alignNames)
    .qtlEnrichmentCoded(
        gwasPip = gwasPip,
        susieQtlRegions = susieQtlRegions,
        coded = .enrCodeRegions(aligned, names(gwasPip)),
        nullFits = nullFits,
        numGwas = numGwas,
        piQtl = piQtl,
        lambda = lambda,
        impN = impN,
        doubleShrinkage = doubleShrinkage,
        besselCorrection = besselCorrection,
        BPPARAM = BPPARAM,
        verbose = verbose,
        seed = seed
    )
}

# The estimate from regions already coded against the GWAS PIP vector
# (`coded`: the kernel-ready fits, the number of distinct QTL-only variants
# and each region's unmatched ids). Shared by qtlEnrichment(), which codes
# by name, and qtlEnrichmentPipeline(), which codes by variant
# (.enrCodeAligned()); `gwasVariants`, the PIPs' variants, is the
# pipeline's and only informs.
# @noRd
.qtlEnrichmentCoded <- function(
    gwasPip,
    susieQtlRegions,
    coded,
    nullFits,
    numGwas = NULL,
    piQtl = NULL,
    lambda = 1.0,
    impN = 25,
    doubleShrinkage = FALSE,
    besselCorrection = TRUE,
    BPPARAM = SerialParam(),
    verbose = TRUE,
    seed = NULL,
    gwasVariants = NULL
) {
    # As in fastenloc, the variant universe is every GWAS variant plus every
    # QTL variant the GWAS vector does not carry, which has GWAS PIP 0.
    totalVariants <- .enrTotalVariants(
        numGwas,
        length(gwasPip) + coded$nQtlOnly
    )
    piGwas <- .enrEstimatePiGwas(gwasPip, totalVariants, numGwas, verbose)
    piQtl <- .enrEstimatePiQtl(susieQtlRegions, piQtl, totalVariants, verbose)
    .enrValidatePi(piGwas, piQtl)
    # cpp11 requires exact integer types for int parameters.
    en <- qtlEnrichmentCpp(
        rGwasPip = as.numeric(gwasPip),
        rQtlSusieFit = coded$fits,
        piGwas = piGwas,
        piQtl = piQtl,
        totalVariants = totalVariants,
        ImpN = as.integer(impN),
        shrinkageLambda = lambda,
        doubleShrinkage = doubleShrinkage,
        besselCorrection = besselCorrection,
        numThreads = .bpLocalThreads(BPPARAM, "qtlEnrichment"),
        seed = if (is.null(seed)) NULL else as.integer(seed)
    )
    .enrCheckRounds(en)
    .enrReportDiagnostics(en, nullFits, length(susieQtlRegions), verbose)
    list_assign(
        en,
        `Invalid MI rounds` = zap(),
        `Null QTL fits` = nullFits,
        unused_xqtl_variants = coded$unmatched
    )
}

# The variant total both priors are taken over: fastenloc's -total_variants,
# which may not fall below the number of distinct variants the data carry (the
# GWAS PIP vector plus the QTL variants it lacks). Without numGwas that number
# is used, with a warning, since a GWAS PIP vector covering only fine-mapped
# blocks makes it far smaller than the variants the GWAS tested.
#' @importFrom checkmate assertNumber
#' @noRd
.enrTotalVariants <- function(numGwas, universe) {
    if (is.null(numGwas)) {
        warn(glue(
            "numGwas is not provided. Using the {universe} distinct variants ",
            "of gwasPip and the QTL fits as the variant total, which biases ",
            "piGwas and piQtl upward unless gwasPip is genome-wide. fastenloc ",
            "requires this total as -total_variants."
        ))
        return(universe)
    }
    assertNumber(numGwas, lower = 1, finite = TRUE)
    if (numGwas < universe) {
        warn(glue(
            "numGwas ({numGwas}) is smaller than the {universe} distinct ",
            "variants of gwasPip and the QTL fits; using {universe}, as ",
            "fastenloc does."
        ))
        return(universe)
    }
    numGwas
}

# piGwas = sum(gwasPip) / the variant total (reported when the total was
# estimated).
# @noRd
.enrEstimatePiGwas <- function(gwasPip, totalVariants, numGwas, verbose) {
    piGwas <- sum(gwasPip) / totalVariants
    if (is.null(numGwas) && verbose) {
        piGwasR <- round(piGwas, 5)
        inform(glue("Estimated piGwas: {piGwasR}\n", .trim = FALSE))
    }
    piGwas
}

# piQtl = the QTL fits' summed PIPs / the variant total, as in fastenloc
# (estimated, with a warning, when piQtl is absent).
# @noRd
.enrEstimatePiQtl <- function(susieQtlRegions, piQtl, totalVariants, verbose) {
    if (!is.null(piQtl)) {
        return(piQtl)
    }
    msg <- glue(
        "piQtl is not provided. Estimating piQtl as the QTL fits' summed PIPs ",
        "over the variant total. Note that this estimate may be biased if ",
        "either 1) the variant total is not the number of variants tested ",
        "(see numGwas), or 2) the single effects only include variables ",
        "inside of credible sets or signal clusters."
    )
    warn(msg)
    piQtl <- sum(map_dbl(susieQtlRegions, .enrPipSum)) / totalVariants
    if (verbose) {
        piQtlR <- round(piQtl, 5)
        inform(glue("Estimated piQtl: {piQtlR}\n", .trim = FALSE))
    }
    piQtl
}

# @noRd
.enrPipSum <- function(x) {
    sum(x$pip)
}

# Each prior probability must lie in (0, 1), and so must their sum: the EM's
# starting cells and its intercept are only defined then.
# @noRd
.enrValidatePi <- function(piGwas, piQtl) {
    if (!(piGwas > 0)) {
        msg <- glue(
            "Cannot perform enrichment analysis. No association signal found ",
            "in GWAS data."
        )
        abort(msg)
    }
    if (!(piQtl > 0)) {
        msg <- glue(
            "Cannot perform enrichment analysis. No QTL associated with the ",
            "molecular phenotype."
        )
        abort(msg)
    }
    if (!(piGwas + piQtl < 1)) {
        piGwasR <- signif(piGwas, 3)
        piQtlR <- signif(piQtl, 3)
        msg <- glue(
            "qtlEnrichment: piGwas ({piGwasR}) and piQtl ({piQtlR}) must sum ",
            "to less than 1. Supply numGwas, the total number of variants ",
            "the GWAS tested."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Both gwasPip and every region's pip must carry variant names, and gwasPip
# one PIP per variant: the variant total counts its names.
# @noRd
.enrValidateNames <- function(gwasPip, susieQtlRegions) {
    if (is.null(names(gwasPip))) {
        msg <- glue(
            "Variant names are missing in gwasPip. Please provide named ",
            "gwasPip data."
        )
        abort(msg)
    }
    dup <- anyDuplicated(names(gwasPip))
    if (dup > 0L) {
        msg <- glue(
            "gwasPip names variant '{names(gwasPip)[[dup]]}' more than once; ",
            "give each variant one PIP."
        )
        abort(msg)
    }
    if (!all(map_lgl(susieQtlRegions, .enrHasPipNames))) {
        msg <- glue(
            "Variant names are missing in susieQtlRegions$pip. Please provide ",
            "susieQtlRegions with named pip data."
        )
        abort(msg)
    }
    invisible(NULL)
}

# TRUE when a region's pip vector is named.
# @noRd
.enrHasPipNames <- function(x) {
    !is.null(names(x$pip))
}

# Every region fit must be one the kernel can sample from: alpha a finite,
# non-negative matrix with a column per pip variant and a row per prior
# variance, and each effect of positive prior variance a distribution (its row
# sums to 1). A malformed fit used to read out of bounds in the kernel (a
# null_weight fit has one more alpha column than variants) or slip NaN past its
# row-sum check, so it is refused here, by name. A fit with no effect of
# positive prior variance is a null fit: it samples nothing, but a tuple made
# only of null fits has nothing to estimate from. Returns the null-fit count.
# @noRd
.enrValidateFits <- function(susieQtlRegions) {
    problems <- compact(imap(susieQtlRegions, .enrFitProblem))
    if (length(problems) > 0L) {
        abort(c(
            "qtlEnrichment: malformed QTL fit(s):",
            set_names(utils::head(unlist(problems), 5L), "x")
        ))
    }
    nullFit <- map_lgl(susieQtlRegions, .enrIsNullFit)
    if (length(nullFit) > 0L && all(nullFit)) {
        abort(glue(
            "qtlEnrichment: none of the {length(nullFit)} QTL fits has an ",
            "effect with prior variance > 0, so there is nothing to sample."
        ))
    }
    sum(nullFit)
}

# The first problem with one region fit, or NULL.
# @noRd
.enrFitProblem <- function(x, label) {
    alpha <- x$alpha
    priorVariance <- x$prior_variance
    if (!is.matrix(alpha) || !is.numeric(alpha)) {
        return(glue("{label}: alpha is not a numeric matrix"))
    }
    if (ncol(alpha) != length(x$pip)) {
        return(glue(
            "{label}: alpha has {ncol(alpha)} columns for {length(x$pip)} ",
            "pip variants"
        ))
    }
    if (nrow(alpha) != length(priorVariance)) {
        return(glue(
            "{label}: alpha has {nrow(alpha)} rows for ",
            "{length(priorVariance)} prior variances"
        ))
    }
    if (!all(is.finite(priorVariance))) {
        return(glue("{label}: prior variances must be finite"))
    }
    if (!all(is.finite(alpha)) || any(alpha < 0)) {
        return(glue("{label}: alpha must be finite and non-negative"))
    }
    rowDev <- abs(rowSums(alpha[priorVariance > 0, , drop = FALSE]) - 1)
    if (any(rowDev > 1e-6)) {
        return(glue(
            "{label}: an alpha row of positive prior variance does not sum ",
            "to 1"
        ))
    }
    NULL
}

# TRUE when no effect of the fit has positive prior variance.
# @noRd
.enrIsNullFit <- function(x) {
    !any(x$prior_variance > 0)
}

# Relabel each region's pip names to the GWAS convention, unless the caller
# already did (alignNames = FALSE).
# @noRd
.enrAlignPipNames <- function(susieQtlRegions, gwasPip, alignNames) {
    if (!alignNames) {
        return(susieQtlRegions)
    }
    .enrAlignRegions(susieQtlRegions, names(gwasPip))
}

# Code every region's variables against the GWAS names in one pass: the
# 0-based index of the same variant in the GWAS PIP vector, or a negative id
# shared by every region carrying the same absent variant. One match() over the
# concatenated names hashes the GWAS panel once, where a lookup per region
# would rehash it for every gene. Also returns each region's absent names and
# the number of distinct absent variants.
# @noRd
.enrCodeRegions <- function(regions, gwasNames) {
    ids <- map(regions, .enrRegionIds)
    flat <- as.character(unlist(ids, use.names = FALSE))
    idx <- match(flat, gwasNames)
    absent <- is.na(idx)
    qtlOnly <- unique(flat[absent])
    codes <- idx - 1L
    codes[absent] <- -match(flat[absent], qtlOnly)
    regionOf <- factor(
        rep(seq_along(regions), lengths(ids)),
        levels = seq_along(regions)
    )
    unmatched <- map(
        unname(split(flat[absent], regionOf[absent])),
        .enrNullIfEmpty
    )
    list(
        fits = map2(regions, unname(split(codes, regionOf)), .enrCodedFit),
        unmatched = set_names(unmatched, names(regions)),
        nQtlOnly = length(qtlOnly)
    )
}

# @noRd
.enrRegionIds <- function(x) {
    names(x$pip)
}

# @noRd
.enrNullIfEmpty <- function(x) {
    if (length(x) == 0L) NULL else x
}

# One region in the form the kernel reads: alpha (double), prior variances
# and the variables' integer codes.
# @noRd
.enrCodedFit <- function(x, codes) {
    alpha <- x$alpha
    if (!is.double(alpha)) {
        storage.mode(alpha) <- "double"
    }
    list(
        alpha = alpha,
        prior_variance = as.numeric(x$prior_variance),
        codes = as.integer(codes)
    )
}

# A round whose (not QTL, not causal) EM cell came out non-positive had a
# variant total too small for its data; its estimates are NaN.
# @noRd
.enrCheckRounds <- function(en) {
    invalid <- en[["Invalid MI rounds"]]
    if (invalid > 0L) {
        nRounds <- length(en[["EM converged per MI round"]])
        msg <- glue(
            "qtlEnrichment: the variant total is too small for the data in ",
            "{invalid} of {nRounds} MI round(s), leaving the EM no variants ",
            "that are neither QTL nor causal. Supply numGwas, the total ",
            "number of variants the GWAS tested."
        )
        abort(msg)
    }
    invisible(NULL)
}

# The kernel no longer prints from its worker threads; its per-round
# diagnostics come back and are reported here. A round whose EM ran out of
# iterations is always worth a warning; the rest is progress, shown only when
# verbose.
# @noRd
.enrReportDiagnostics <- function(en, nullFits, nFits, verbose) {
    converged <- en[["EM converged per MI round"]]
    if (!all(converged)) {
        warn(glue(
            "qtlEnrichment: the EM did not converge within 1000 iterations ",
            "in {sum(!converged)} of {length(converged)} MI round(s)."
        ))
    }
    if (!verbose) {
        return(invisible(NULL))
    }
    if (nullFits > 0L) {
        inform(glue(
            "{nullFits} of {nFits} QTL fits have no effect with prior ",
            "variance > 0 and contribute no sampled variants."
        ))
    }
    missing <- en[["Missing proportion per MI round"]]
    inform(glue(
        "Proportion of sampled xQTL variants missing from the GWAS PIP ",
        "vector: {signif(mean(missing), 3)} (mean over {length(missing)} MI ",
        "rounds); they enter the estimate with GWAS PIP 0."
    ))
    iterations <- en[["EM iterations per MI round"]]
    if (length(iterations) > 0L) {
        inform(glue(
            "EM iterations per MI round: median {stats::median(iterations)} ",
            "(range {min(iterations)}-{max(iterations)})."
        ))
    }
    removed <- length(missing) - en[["Effective MI rounds"]]
    if (removed > 0) {
        inform(glue("Outlier filtering removed {removed} MI round(s)."))
    }
    invisible(NULL)
}
