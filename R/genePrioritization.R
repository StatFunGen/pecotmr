#' Prioritize genes by their features (PoPS)
#'
#' @description Scores each gene by how well its features predict its TWAS
#'   or MR association, as the Polygenic Priority Score (PoPS; Weeks et al.
#'   2023) does with MAGMA's gene z-scores. The genes' statistics from a
#'   \code{\link{causalInferencePipeline}} result are regressed on their
#'   features by generalized least squares, whitened by the correlation the
#'   LD between their variants induces (\code{\link{geneCorrelation}}), and
#'   each gene is scored from its own features by the fit that left its
#'   chromosome out.
#'
#' @details For each GWAS study and each chromosome, the genes with a
#'   statistic on the other chromosomes (outside the HLA region, by default)
#'   are the training genes. Their statistics and features are whitened by
#'   the genes' correlation, with PoPS's ridge on its diagonal; the whitened
#'   intercept is projected out of both, and the gene covariates out of the
#'   statistics. \code{fitMethod} is then fitted to them, with no intercept of
#'   its own, and each gene of the held-out chromosome is scored as its
#'   features times the fitted effects. Every gene of \code{features} gets a
#'   score, whether or not it has a statistic.
#'
#'   Unlike PoPS, no feature is screened out before the fit: a screen on
#'   every gene would choose features by the held-out genes' own statistics.
#'   The fit methods are those whose engine can fit without an intercept:
#'   \code{"ridge"} (default, PoPS's), \code{"lasso"} and \code{"enet"}
#'   (glmnet), \code{"l0learn"}, \code{"abess"}, \code{"mrash"},
#'   \code{"sparsevb"}, \code{"susie"}, and the DPR methods \code{"dprVb"},
#'   \code{"dprGibbs"} and \code{"dprAdaptiveGibbs"}, which take the whitened
#'   intercept as their covariate. Sparse methods suit features with few
#'   effects; ridge suits effects spread over many.
#'
#'   With \code{covariates = TRUE}, the gene covariates projected out are, for
#'   TWAS, the number of variants with a nonzero weight, the weights'
#'   variance under the LD (\eqn{w' R w}), the logarithms of both, and the
#'   weights' cross-validated \eqn{r^2}; for MR, the number of instruments and
#'   its logarithm. A covariate missing for any gene (a CV \eqn{r^2} the
#'   weights do not record, say) is left out, with a message. They play the
#'   part of MAGMA's gene size and density in PoPS.
#'
#'   The features are read as a dense matrix of the training genes for each
#'   held-out chromosome, and whitening keeps them dense.
#'
#' @inheritParams geneCorrelation
#' @param features A \code{RangedSummarizedExperiment} of genes by features:
#'   its rows the genes, named as \code{x}'s \code{trait}, with their
#'   positions (the chromosome each is held out with, and the HLA region);
#'   its first assay the features, binary or continuous, dense or sparse.
#' @param fitMethod The method fitted to the whitened data; see Details.
#' @param methodArgs Optional options for \code{fitMethod}, built with its
#'   constructor: \code{\link{GlmnetOptions}} for \code{"ridge"},
#'   \code{"lasso"} and \code{"enet"}, \code{\link{L0learnOptions}},
#'   \code{\link{AbessOptions}}, \code{\link{MrashOptions}},
#'   \code{\link{SparsevbOptions}}, \code{\link{SusieOptions}} or
#'   \code{\link{DprOptions}}. The intercept is always off; an option that
#'   turns it on is refused.
#' @param covariates Logical: project the gene covariates (see Details) out
#'   of the statistics. Default \code{TRUE}.
#' @param excludeHla Logical: leave the genes of the HLA region (a TSS on
#'   chromosome 6 between 20 and 40 Mb, as PoPS) out of every fit. They are
#'   still scored. Default \code{TRUE}.
#' @param credibleSetParam For \code{fitMethod = "susie"}, built with
#'   \code{\link{CredibleSetParam}}: its \code{maxNumSingleEffects} is
#'   SuSiE's \code{L} and its \code{maxNumSingleEffectsGreedy} its
#'   \code{L_greedy}. Its other settings shape credible sets, which the fit
#'   does not use. Naming it for another method is an error.
#' @param seed Optional integer seed for the fits' random numbers (their
#'   cross-validation folds, or a sampler's draws).
#' @return A list with one element per GWAS study, named by study, each a
#'   list of:
#'   \describe{
#'     \item{\code{scores}}{A \code{GRanges}, one range per gene of
#'       \code{features} (its \code{rowRanges}, named by gene), with columns
#'       \code{y}, the gene's statistic on \code{scale} (\code{NA} without
#'       one); \code{score}, its prediction from the fit that held out its
#'       chromosome; and \code{chromosome}, that chromosome.}
#'     \item{\code{coefficients}}{A features by held-out chromosomes matrix
#'       of the fitted effects.}
#'   }
#' @seealso \code{\link{geneCorrelation}},
#'   \code{\link{causalInferencePipeline}}
#' @examples
#' set.seed(1)
#' # 40 genes, 20 on each of two chromosomes, each with 5 variants of its own.
#' chr <- rep(c("chr1", "chr2"), each = 100)
#' pos <- rep(100L * seq_len(100), 2)
#' ids <- paste(chr, pos, "C", "T", sep = ":")
#' gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'     SNP = ids, A1 = "T", A2 = "C", Z = rnorm(200, sd = 2)
#' )
#' ld <- LdMatrix(
#'     diag(200),
#'     variants = data.frame(
#'         chrom = chr, pos = pos, SNP = ids, A1 = "T", A2 = "C"
#'     )
#' )
#' gwas <- GwasSumStats(
#'     "gwas", list(gr), "hg19",
#'     ldReference = ld, qcOptions = list(qc = "done")
#' )
#' genes <- sprintf("gene%d", 1:40)
#' first <- 5L * (0:39) + 1L
#' tw <- TwasWeights(
#'     rep("qtl", 40), rep("ctx", 40), genes, rep("ridge", 40),
#'     lapply(first, function(i) TwasWeightsRow(ids[i + 0:4], rnorm(5)))
#' )
#' res <- causalInferencePipeline(
#'     gwas,
#'     twasWeights = tw,
#'     weightSelectionParam = WeightSelectionParam(methods = "all")
#' )
#' features <- SummarizedExperiment::SummarizedExperiment(
#'     list(features = matrix(rbinom(40 * 8, 1, 0.5), 40,
#'         dimnames = list(genes, sprintf("set%d", 1:8)))),
#'     rowRanges = GenomicRanges::GRanges(
#'         chr[first],
#'         IRanges::IRanges(pos[first], width = 500),
#'         names = genes
#'     )
#' )
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'     out <- genePrioritization(
#'         res, features, gwas,
#'         twasWeights = tw, covariates = FALSE, seed = 1
#'     )
#'     out$gwas$scores
#' }
#' @export
genePrioritization <- function(
    x,
    features,
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    analysis = c("twas", "mr"),
    method = NULL,
    context = NULL,
    scale = c("oneSided", "signed", "absolute", "squared"),
    fitMethod = "ridge",
    methodArgs = NULL,
    covariates = TRUE,
    excludeHla = TRUE,
    mrParam = MrParam(),
    alleleFlip = TRUE,
    credibleSetParam = CredibleSetParam(),
    seed = NULL
) {
    analysis <- arg_match(analysis)
    scale <- arg_match(scale)
    fitMethod <- arg_match(fitMethod, names(.gpFitMethods))
    cfg <- .gpConfig(
        features,
        scale,
        fitMethod,
        methodArgs,
        credibleSetParam,
        !missing(credibleSetParam),
        covariates,
        excludeHla,
        analysis,
        seed
    )
    studies <- .gcCompute(
        x,
        gwasSumStats,
        twasWeights,
        fineMappingResult,
        analysis,
        method,
        context,
        mrParam,
        alleleFlip
    )
    .gpScoreStudies(studies, cfg, twasWeights %||% fineMappingResult)
}

# The run's settings, checked: the genes of `features` and their features,
# the fit method and its options, and the rest of genePrioritization()'s
# arguments the fits read.
# @noRd
.gpConfig <- function(
    features,
    scale,
    fitMethod,
    methodArgs,
    credibleSetParam,
    supplied,
    covariates,
    excludeHla,
    analysis,
    seed
) {
    .gpAssertArgs(
        covariates,
        excludeHla,
        seed,
        fitMethod,
        credibleSetParam,
        supplied
    )
    list(
        genes = .gpFeatureGenes(features, excludeHla),
        features = SummarizedExperiment::assay(features),
        fitMethod = fitMethod,
        options = .gpOptions(fitMethod, methodArgs, credibleSetParam),
        scale = scale,
        covariates = covariates,
        analysis = analysis,
        seed = seed
    )
}

# Each study's scores (.gpStudy()), under `cfg`'s seed when it has one. The
# CV r^2 covariate is read off the TWAS weights `source`.
# @noRd
.gpScoreStudies <- function(studies, cfg, source) {
    rsq <- if (cfg$covariates && cfg$analysis == "twas") .gpCvRsq(source)
    if (!is.null(cfg$seed)) {
        withr::local_seed(cfg$seed)
    }
    map(studies, .gpStudy, cfg = cfg, rsq = rsq)
}

# The fit methods: the Options constructor each takes, and the option that
# turns its intercept off (none for DPR, which takes the whitened intercept
# as its covariate instead).
# @noRd
.gpFitMethods <- list(
    ridge = list(ctor = "GlmnetOptions", off = list(intercept = FALSE)),
    lasso = list(ctor = "GlmnetOptions", off = list(intercept = FALSE)),
    enet = list(ctor = "GlmnetOptions", off = list(intercept = FALSE)),
    l0learn = list(ctor = "L0learnOptions", off = list(intercept = FALSE)),
    abess = list(ctor = "AbessOptions", off = list(fit.intercept = FALSE)),
    mrash = list(ctor = "MrashOptions", off = list(intercept = FALSE)),
    sparsevb = list(ctor = "SparsevbOptions", off = list(intercept = FALSE)),
    susie = list(ctor = "SusieOptions", off = list(intercept = FALSE)),
    dprVb = list(ctor = "DprOptions", off = list()),
    dprGibbs = list(ctor = "DprOptions", off = list()),
    dprAdaptiveGibbs = list(ctor = "DprOptions", off = list())
)

# The flags and seed, and `credibleSetParam`, which only SuSiE reads: naming
# it (`supplied`) for another method is refused.
# @noRd
#' @importFrom checkmate assertClass assertFlag assertInt
.gpAssertArgs <- function(
    covariates,
    excludeHla,
    seed,
    fitMethod,
    credibleSetParam,
    supplied
) {
    assertFlag(covariates)
    assertFlag(excludeHla)
    assertInt(seed, null.ok = TRUE)
    assertClass(credibleSetParam, "CredibleSetParam")
    if (supplied && fitMethod != "susie") {
        abort(c(
            glue(
                "genePrioritization: `credibleSetParam` sets SuSiE's L, and ",
                "the {fitMethod} fit has none."
            ),
            i = "Leave it unset, or use `fitMethod = \"susie\"`."
        ))
    }
    invisible(NULL)
}

# The caller's options for `fitMethod`, rebuilt by its constructor (which
# refuses another engine's settings) with the intercept off; an option that
# turns it on, or DPR's covariates, which the whitening sets, is refused.
# SuSiE's come back as its argument list, with credibleSetParam's L and
# L_greedy seeded as fine mapping seeds them (.fmSeedTokenDefaults()).
# @noRd
.gpOptions <- function(fitMethod, methodArgs, credibleSetParam) {
    spec <- .gpFitMethods[[fitMethod]]
    .assertMethodOptions(methodArgs %||% list(), spec$ctor, "methodArgs")
    constructor <- if (.isMethodOptions(methodArgs)) {
        S4Vectors::metadata(methodArgs)$constructor
    }
    if (!is.null(constructor) && !identical(constructor, spec$ctor)) {
        abort(glue(
            "genePrioritization: `methodArgs` for {fitMethod} must be built ",
            "with {spec$ctor}(), not {constructor}()."
        ))
    }
    given <- as.list(methodArgs %||% list())
    owned <- c(names(spec$off), if (length(spec$off) == 0L) "w")
    clash <- intersect(names(given), owned)
    turnsOn <- clash[!map_lgl(given[clash], isFALSE)]
    if (length(turnsOn) > 0L) {
        abort(glue(
            "genePrioritization: the {fitMethod} fit has no intercept of its ",
            "own (the whitened intercept is projected out), so ",
            "`methodArgs` cannot set {str_flatten(turnsOn, ', ')}."
        ))
    }
    built <- exec(spec$ctor, !!!list_modify(given, !!!spec$off))
    if (fitMethod != "susie") {
        return(built)
    }
    .fmSeedTokenDefaults(
        "susie",
        list(susie = as.list(built)),
        credibleSetParam$maxNumSingleEffects %||% 10L,
        credibleSetParam$maxNumSingleEffectsGreedy
    )
}

# The genes of `features`: their ids, chromosomes and whether each lies in
# the HLA region, which `excludeHla` keeps out of the fits.
# @noRd
#' @importFrom checkmate assertClass
.gpFeatureGenes <- function(features, excludeHla) {
    assertClass(features, "RangedSummarizedExperiment")
    ids <- rownames(features)
    if (is.null(ids) || anyDuplicated(ids) > 0L) {
        abort(str_c(
            "genePrioritization: `features` must have its genes as row ",
            "names, each gene once."
        ))
    }
    if (length(SummarizedExperiment::assays(features)) == 0L) {
        abort("genePrioritization: `features` has no assay of features.")
    }
    ranges <- .ucscRanges(SummarizedExperiment::rowRanges(features))
    chromosome <- as.character(GenomicRanges::seqnames(ranges))
    tss <- GenomicRanges::start(GenomicRanges::resize(ranges, 1L))
    hla <- chromosome == "chr6" & tss >= 20e6 & tss <= 40e6
    list(
        ids = ids,
        ranges = unname(ranges),
        chromosome = chromosome,
        train = !(excludeHla & hla)
    )
}

# Each weight set's cross-validated r^2, by its key (.cellKey()); NA where
# its weights record none.
# @noRd
.gpCvRsq <- function(source) {
    cells <- .weightCells(source, WeightSelectionParam())
    set_names(cells$rsq, .cellKey(cells))
}

# One GWAS study's scores (.gcStudyPairs()'s `study`).
# @noRd
.gpStudy <- function(study, cfg, rsq) {
    response <- .gpResponse(study$genes, cfg$scale, cfg$genes$ids)
    correlation <- .gcScaledMatrix(study, cfg$scale)
    nuisance <- if (cfg$covariates) {
        .gpCovariates(study, cfg$analysis, rsq)
    }
    folds <- unique(cfg$genes$chromosome)
    fits <- map(
        folds,
        .gpFold,
        response = response,
        correlation = correlation,
        nuisance = nuisance,
        cfg = cfg
    )
    coefficients <- `dimnames<-`(
        exec(cbind, !!!map(fits, "coef")),
        list(colnames(cfg$features), folds)
    )
    score <- set_names(
        list_c(map(fits, "score")),
        list_c(map(fits, "genes"))
    )
    list(
        scores = .gpScores(cfg$genes, response, score),
        coefficients = coefficients
    )
}

# Each gene's statistic on `scale`, for the genes `features` holds; the
# others are dropped, with a message.
# @noRd
.gpResponse <- function(studyGenes, scale, featureIds) {
    ids <- map_chr(studyGenes, "gene")
    z <- map_dbl(studyGenes, "statistic")
    value <- switch(
        scale,
        signed = z,
        absolute = abs(z),
        squared = z^2,
        oneSided = .zOneSided(z)
    )
    known <- is_in(ids, featureIds)
    if (!all(known)) {
        inform(glue(
            "genePrioritization: {sum(!known)} gene(s) with a statistic are ",
            "not in `features` and are left out."
        ))
    }
    set_names(value[known], ids[known])
}

# The gene covariates (see genePrioritization()'s Details), a row per gene;
# a column missing for any gene is left out, with a message.
# @noRd
.gpCovariates <- function(study, analysis, rsq) {
    ids <- map_chr(study$genes, "gene")
    table <- if (analysis == "twas") {
        nonzero <- map_int(study$genes, .gpNonzeroWeights)
        cbind(
            nVariants = nonzero,
            logNVariants = log(nonzero),
            weightVariance = study$variance,
            logWeightVariance = log(study$variance),
            cvRsq = unname(rsq[map_chr(study$genes, "cellKey")])
        )
    } else {
        instruments <- map_dbl(study$genes, "nIV")
        cbind(
            nInstruments = instruments,
            logNInstruments = log(instruments)
        )
    }
    complete <- apply(is.finite(table), 2L, all)
    if (!all(complete)) {
        inform(glue(
            "genePrioritization: covariate(s) ",
            "{str_flatten(colnames(table)[!complete], ', ')} are missing ",
            "for some genes and are left out."
        ))
    }
    if (!any(complete)) {
        return(NULL)
    }
    `rownames<-`(table[, complete, drop = FALSE], ids)
}

# @noRd
.gpNonzeroWeights <- function(gene) {
    sum(gene$coef != 0)
}

# The fit that holds out chromosome `fold`: trained on the genes with a
# statistic elsewhere (and outside HLA, when asked), and the scores of the
# genes on `fold`.
# @noRd
.gpFold <- function(fold, response, correlation, nuisance, cfg) {
    genes <- cfg$genes
    features <- cfg$features
    held <- genes$ids[genes$chromosome == fold]
    trainIds <- genes$ids[genes$chromosome != fold & genes$train]
    train <- response[is_in(names(response), trainIds)]
    if (length(train) < 2L) {
        abort(glue(
            "genePrioritization: holding out {fold} leaves {length(train)} ",
            "gene(s) with a statistic to fit on."
        ))
    }
    whitened <- .glsWhiten(
        train,
        as.matrix(features[names(train), , drop = FALSE]),
        correlation,
        covariates = nuisance
    )
    coef <- .gpFit(cfg$fitMethod, whitened, cfg$options)
    list(
        coef = coef,
        genes = held,
        score = as.numeric(features[held, , drop = FALSE] %*% coef)
    )
}

# `fitMethod`'s effects of the whitened features on the whitened statistics.
# @noRd
.gpFit <- function(fitMethod, whitened, options) {
    X <- whitened$X
    y <- whitened$y
    weights <- switch(
        fitMethod,
        ridge = ridgeWeights(X, y, methodArgs = options),
        lasso = lassoWeights(X, y, backend = "glmnet", methodArgs = options),
        enet = enetWeights(X, y, backend = "glmnet", methodArgs = options),
        l0learn = l0learnWeights(
            X,
            y,
            backend = "l0learn",
            methodArgs = options
        ),
        abess = l0learnWeights(X, y, backend = "abess", methodArgs = options),
        mrash = mrashWeights(X, y, methodArgs = options),
        sparsevb = sparsevbWeights(X, y, methodArgs = options),
        susie = .gpSusie(X, y, options),
        .gpDpr(fitMethod, X, y, whitened$intercept, options)
    )
    as.numeric(weights)
}

# SuSiE's posterior mean effects, fitted as fine-mapping fits it.
# @noRd
.gpSusie <- function(X, y, options) {
    fit <- .fmFitSusieIndiv(X, y, "susie", userArgs = as.list(options))
    susieWeights(X = X, susieFit = fit)
}

# A DPR method's effects, the whitened intercept its covariate.
# @noRd
.gpDpr <- function(fitMethod, X, y, intercept, options) {
    withIntercept <- exec(
        "DprOptions",
        !!!list_modify(as.list(options), w = intercept)
    )
    fn <- switch(
        fitMethod,
        dprVb = dprVbWeights,
        dprGibbs = dprGibbsWeights,
        dprAdaptiveGibbs = dprAdaptiveGibbsWeights
    )
    fn(X, y, methodArgs = withIntercept)
}

# The scores as GRanges: every gene of `features`, its statistic (NA without
# one), its score and the chromosome it was held out with.
# @noRd
.gpScores <- function(genes, response, score) {
    scored <- S4Vectors::`mcols<-`(
        genes$ranges,
        value = S4Vectors::DataFrame(
            y = unname(response[genes$ids]),
            score = unname(score[genes$ids]),
            chromosome = genes$chromosome
        )
    )
    `names<-`(scored, genes$ids)
}

# =============================================================================
# Gene correlation, and whitening by it
# =============================================================================

#' Correlation of genes' TWAS or MR statistics
#'
#' @description For the genes of a \code{\link{causalInferencePipeline}}
#'   result, the correlation between their TWAS z-scores, or their MR
#'   estimates, that the LD between their variants induces, against each GWAS
#'   study. It is the error correlation a regression of gene-level statistics
#'   on gene features whitens by (generalized least squares), as PoPS does
#'   with MAGMA's gene correlations.
#'
#' @details A gene's TWAS z-score is linear in the GWAS z-scores,
#'   \eqn{z_g = w_g' z / \sqrt{w_g' R w_g}}, for its weights \eqn{w_g} on the
#'   GWAS variants, matched and allele-aligned as the pipeline matches them,
#'   and the GWAS LD \eqn{R}. Two genes' z-scores therefore correlate as
#'   \deqn{\rho_{gh} = \frac{w_g' R w_h}{\sqrt{w_g' R w_g \; w_h' R w_h}}}{
#'   rho_gh = w_g' R w_h / sqrt(w_g' R w_g w_h' R w_h)}
#'   exactly, whatever the GWAS signal. An MR estimate pools Wald ratios by
#'   inverse-variance weights. With those weights held fixed it is linear in
#'   the GWAS z-scores too, \eqn{\hat\beta_g = a_g' z}, and correlates as
#'   above with \eqn{a_g} in place of \eqn{w_g}. That is a first-order
#'   approximation: the weights themselves depend on the GWAS z-scores.
#'
#'   Two genes correlate only when their variants fall in a common row of the
#'   GWAS study (an LD block, or a chromosome), the rows the pipeline read
#'   each gene's statistic from; genes that share no row get 0. The LD is the
#'   study's \code{ldReference}, the one each statistic was computed with.
#'   Each gene's statistic is recomputed from the weights and the GWAS and
#'   compared with \code{x}'s, so that inputs other than the ones \code{x}
#'   came from are refused rather than giving a correlation for other
#'   statistics.
#'
#'   \code{scale} gives the correlation of the statistic a regression would
#'   take as its response, for jointly normal statistics of correlation
#'   \eqn{\rho}:
#'   \describe{
#'     \item{\code{"oneSided"}}{The one-sided z-score of the two-sided
#'       p-value, \eqn{\Phi^{-1}(1 - p)}, the scale of MAGMA's gene statistic
#'       that PoPS regresses: \eqn{\sum_k b_k^2 |\rho|^k} over the Hermite
#'       coefficients \eqn{b_k} of that transformation, computed numerically.}
#'     \item{\code{"signed"}}{The z-scores themselves: \eqn{\rho}.}
#'     \item{\code{"absolute"}}{\eqn{|z|}: \eqn{(\sqrt{1 - \rho^2} +
#'       \rho \arcsin \rho - 1) / (\pi / 2 - 1)}.}
#'     \item{\code{"squared"}}{\eqn{z^2}: \eqn{\rho^2}.}
#'   }
#'   All but \code{"signed"} depend on \eqn{|\rho|} only. Each keeps the
#'   matrix positive semi-definite, as a power series in \eqn{\rho} with
#'   non-negative coefficients.
#'
#' @param x A \code{\link{causalInferencePipeline}} result.
#' @param gwasSumStats The \code{\link{GwasSumStats}} \code{x} was computed
#'   from.
#' @param twasWeights,fineMappingResult The weights \code{x} was computed
#'   from: for TWAS rows \code{twasWeights}, else the
#'   \code{fineMappingResult}, as the pipeline chooses; for MR rows the
#'   \code{fineMappingResult}.
#' @param analysis Which rows of \code{x}: \code{"twas"} (default), their
#'   TWAS z-scores, or \code{"mr"}, their Wald-ratio MR estimates.
#' @param method,context Optional: the one method and the one context whose
#'   rows are used. After them each gene (\code{trait}) must have one row per
#'   GWAS study; without them \code{x} must have one already.
#' @param scale The statistic the correlation is of: \code{"oneSided"}
#'   (default), \code{"signed"}, \code{"absolute"} or \code{"squared"}; see
#'   Details.
#' @param mrParam,alleleFlip As passed to \code{causalInferencePipeline} for
#'   \code{x}.
#' @return A list with one element per GWAS study, named by study: a sparse
#'   symmetric \code{Matrix} of the genes' correlations, its rows and columns
#'   named by \code{trait} in genomic order. A gene with no statistic in
#'   \code{x} has no row.
#' @seealso \code{\link{causalInferencePipeline}}, \code{\link{twasZ}}
#' @examples
#' ids <- sprintf("chr1:%d:C:T", 100L * seq_len(8))
#' ld <- LdMatrix(
#'     0.5^abs(outer(1:8, 1:8, "-")),
#'     variants = data.frame(
#'         chrom = "chr1", pos = 100L * seq_len(8), SNP = ids,
#'         A1 = "T", A2 = "C"
#'     )
#' )
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L * 1:8, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'     SNP = ids, A1 = "T", A2 = "C", Z = c(3, 2, 1, 0, -1, 0.5, 1, 2)
#' )
#' gwas <- GwasSumStats(
#'     "gwas", list(gr), "hg19",
#'     ldReference = ld, qcOptions = list(qc = "done")
#' )
#' tw <- TwasWeights(
#'     rep("qtl", 2), rep("ctx", 2), c("geneA", "geneB"), rep("lasso", 2),
#'     list(
#'         TwasWeightsRow(ids[1:5], c(0.4, 0.2, 0, 0.1, 0.3)),
#'         TwasWeightsRow(ids[3:8], c(0.2, 0.3, 0.1, 0, 0.2, 0.1))
#'     )
#' )
#' res <- causalInferencePipeline(
#'     gwas,
#'     twasWeights = tw,
#'     weightSelectionParam = WeightSelectionParam(methods = "all")
#' )
#' geneCorrelation(res, gwas, twasWeights = tw)
#' @export
geneCorrelation <- function(
    x,
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    analysis = c("twas", "mr"),
    method = NULL,
    context = NULL,
    scale = c("oneSided", "signed", "absolute", "squared"),
    mrParam = MrParam(),
    alleleFlip = TRUE
) {
    analysis <- arg_match(analysis)
    scale <- arg_match(scale)
    studies <- .gcCompute(
        x,
        gwasSumStats,
        twasWeights,
        fineMappingResult,
        analysis,
        method,
        context,
        mrParam,
        alleleFlip
    )
    map(studies, .gcScaledMatrix, scale = scale)
}

# Each GWAS study's genes (.gcGene()) and the signed correlations between
# them (.gcStudyPairs()), for the rows of `x` the arguments choose.
# @noRd
#' @importFrom checkmate assertClass assertFlag assertString
.gcCompute <- function(
    x,
    gwasSumStats,
    twasWeights,
    fineMappingResult,
    analysis,
    method,
    context,
    mrParam,
    alleleFlip
) {
    assertClass(x, "GRanges")
    assertClass(gwasSumStats, "GwasSumStats")
    assertClass(twasWeights, "TwasWeights", null.ok = TRUE)
    assertClass(fineMappingResult, "QtlFineMappingResult", null.ok = TRUE)
    assertString(method, null.ok = TRUE)
    assertString(context, null.ok = TRUE)
    .assertMethodParam(mrParam, "MrParam")
    assertFlag(alleleFlip)
    .inputsBuild(
        list(x, gwasSumStats, twasWeights, fineMappingResult),
        "geneCorrelation"
    )
    source <- .gcWeightSource(analysis, twasWeights, fineMappingResult)
    rows <- .gcRows(x, analysis, method, context)
    keys <- .weightCellIds(source)
    genes <- map(
        seq_len(nrow(rows)),
        .gcGene,
        rows = rows,
        cellIndex = match(rows$cellKey, .cellKey(keys)),
        keys = keys,
        source = source,
        analysis = analysis,
        gwasSumStats = gwasSumStats,
        gwasIndex = .rowSpans(gwasSumStats),
        mrParam = mrParam,
        alleleFlip = alleleFlip
    )
    studies <- factor(rows$gwasStudy, levels = unique(rows$gwasStudy))
    imap(
        split(genes, studies),
        .gcStudyPairs,
        gwasSumStats = gwasSumStats,
        analysis = analysis
    )
}

# The weights `x`'s rows of `analysis` were computed from: for TWAS the
# pipeline's choice, `twasWeights` else `fineMappingResult`; for MR the
# fine-mapping result.
# @noRd
.gcWeightSource <- function(analysis, twasWeights, fineMappingResult) {
    if (analysis == "mr") {
        if (is.null(fineMappingResult)) {
            abort(str_c(
                "geneCorrelation: MR rows need the `fineMappingResult` `x` ",
                "was computed from."
            ))
        }
        return(fineMappingResult)
    }
    source <- twasWeights %||% fineMappingResult
    if (is.null(source)) {
        abort(str_c(
            "geneCorrelation: TWAS rows need the weights `x` was scored ",
            "with: pass `twasWeights`, or the `fineMappingResult` they came ",
            "from."
        ))
    }
    source
}

# The rows of `x` the correlation is over, in genomic order: those of
# `analysis` (and of `method` and `context`, when given) that hold a
# statistic, one per gene and GWAS study. Each carries the value its
# recomputed statistic must match (`target`: the TWAS z, or the MR estimate)
# and its weight set's key (.cellKey()).
# @noRd
.gcRows <- function(x, analysis, method, context) {
    mc <- as_tibble(as.data.frame(S4Vectors::mcols(x)))
    .gcRequireColumns(mc)
    target <- if (analysis == "twas") mc$twasZ else mc$waldRatio
    statistic <- if (analysis == "twas") {
        mc$twasZ
    } else {
        mc$waldRatio / mc$waldRatioSe
    }
    keep <- which(
        mc$analysis == analysis &
            is.finite(target) &
            is.finite(statistic) &
            .gcMatches(mc$method, method) &
            .gcMatches(mc$context, context)
    )
    if (length(keep) == 0L) {
        abort(glue(
            "geneCorrelation: `x` holds no {analysis} row with a statistic",
            "{.gcFilterLabel(method, context)}."
        ))
    }
    ordered <- keep[order(x[keep])]
    chosen <- mc[ordered, ]
    cells <- tibble(
        study = chosen$qtlStudy,
        context = chosen$context,
        trait = chosen$trait,
        method = chosen$method,
        !!!select(chosen, any_of(.JOINT_COLUMNS))
    )
    rows <- tibble::add_column(
        chosen,
        target = target[ordered],
        statistic = statistic[ordered],
        cellKey = .cellKey(cells)
    )
    .gcAssertOneRowPerGene(rows, analysis)
    rows
}

# @noRd
.gcRequireColumns <- function(mc) {
    need <- c(
        "qtlStudy",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "analysis",
        "twasZ",
        "waldRatio",
        "waldRatioSe",
        "nIV"
    )
    absent <- setdiff(need, names(mc))
    if (length(absent) > 0L) {
        abort(glue(
            "geneCorrelation: `x` is not a causalInferencePipeline() ",
            "result; it lacks column(s) {str_flatten(absent, ', ')}."
        ))
    }
    invisible(NULL)
}

# Whether each of `values` is `wanted`; all are when nothing is wanted.
# @noRd
.gcMatches <- function(values, wanted) {
    if (is.null(wanted)) {
        return(rep(TRUE, length(values)))
    }
    is_in(values, wanted)
}

# @noRd
.gcFilterLabel <- function(method, context) {
    parts <- c(
        if (!is.null(method)) glue("method '{method}'"),
        if (!is.null(context)) glue("context '{context}'")
    )
    if (length(parts) == 0L) "" else str_c(" for ", str_flatten(parts, ", "))
}

# @noRd
.gcAssertOneRowPerGene <- function(rows, analysis) {
    key <- str_c(rows$trait, rows$gwasStudy, sep = "\r")
    twice <- unique(rows$trait[duplicated(key)])
    if (length(twice) == 0L) {
        return(invisible(NULL))
    }
    abort(c(
        glue(
            "geneCorrelation: `x` holds more than one {analysis} row against ",
            "one GWAS study for gene(s) {str_flatten(head(twice, 5L), ', ')}",
            "{if (length(twice) > 5L) ', ...' else ''}."
        ),
        i = str_c(
            "Each gene needs one statistic: pass `method` and `context`, or ",
            "subset `x` (to one QTL study, say)."
        )
    ))
}

# Row `i` of `rows` as a gene: its coefficients on the GWAS variants (the
# TWAS weights, or the MR estimate's linear weights), the GWAS z-scores there,
# the rows of its GWAS study the pipeline read, and the statistic `x` holds.
# @noRd
.gcGene <- function(
    i,
    rows,
    cellIndex,
    keys,
    source,
    analysis,
    gwasSumStats,
    gwasIndex,
    mrParam,
    alleleFlip
) {
    row <- rows[i, ]
    label <- glue(
        "(study={row$qtlStudy}, context={row$context}, ",
        "trait={row$trait}, method={row$method}) x GWAS '{row$gwasStudy}'"
    )
    idx <- cellIndex[[i]]
    input <- if (is.na(idx)) {
        NULL
    } else if (analysis == "twas") {
        .cipTwasTuple(idx, keys, source)
    } else {
        .cipMrTuple(idx, keys, source)
    }
    gene <- .gcGeneCoefficients(
        input,
        row$gwasStudy,
        analysis,
        gwasSumStats,
        gwasIndex,
        mrParam,
        alleleFlip
    )
    if (is.null(gene)) {
        .gcAbortMismatch(label)
    }
    list_assign(
        gene,
        gene = row$trait,
        label = label,
        target = row$target,
        statistic = row$statistic,
        nIV = row$nIV,
        cellKey = row$cellKey
    )
}

# A weight set's coefficients against GWAS study `gStudy`, over the study's
# statistics in the window of its variants, as the pipeline reads them; NULL
# when it has none there.
# @noRd
.gcGeneCoefficients <- function(
    input,
    gStudy,
    analysis,
    gwasSumStats,
    gwasIndex,
    mrParam,
    alleleFlip
) {
    if (is.null(input)) {
        return(NULL)
    }
    window <- range(input$variants)
    gRows <- .cipStudyRows(gwasIndex, gwasSumStats, window)[[gStudy]]
    if (is.null(gRows)) {
        return(NULL)
    }
    gdf <- .cipGwasDf(gwasSumStats, gRows, window, gStudy)
    coef <- if (analysis == "twas") {
        .gcTwasCoefficients(input$weightsInfo, gdf, alleleFlip)
    } else {
        .gcMrCoefficients(input$fmrEntry, gdf, mrParam, alleleFlip)
    }
    if (is.null(coef)) NULL else list_assign(coef, rows = gRows)
}

# A TWAS weight set's weights on the GWAS variants (.cipAlignWeights()); NULL
# with fewer than two, which the pipeline scores no z for.
# @noRd
.gcTwasCoefficients <- function(weightsInfo, gdf, alleleFlip) {
    aligned <- .cipAlignWeights(
        weightsInfo$weights,
        weightsInfo$variantIds,
        gdf,
        alleleFlip
    )
    if (length(aligned$variantIds) < 2L) {
        return(NULL)
    }
    list(
        variantIds = aligned$variantIds,
        coef = aligned$weights,
        z = aligned$z
    )
}

# An MR estimate as a linear form in the GWAS z-scores, with its
# inverse-variance weights held fixed: the per-variant IVW pool of
# .cipComputeMr(), or the per-credible-set one of .cipComputeMrCsAware().
# NULL when it has no instrument.
# @noRd
.gcMrCoefficients <- function(fmrEntry, gdf, mrParam, alleleFlip) {
    if (mrParam$method == "csAware") {
        return(.gcCsAwareCoefficients(
            fmrEntry,
            gdf,
            mrParam$cpipCutoff,
            alleleFlip
        ))
    }
    .gcIvwCoefficients(
        .cipMrIvwInputs(fmrEntry, gdf, mrParam$pipCutoff, alleleFlip)
    )
}

# The per-variant IVW estimate sum(w * ratio) / sum(w) with ratio = z * seY /
# betaX, as coefficients on z; the instruments .ivwPool() drops get none.
# @noRd
.gcIvwCoefficients <- function(iv) {
    if (is.null(iv)) {
        return(NULL)
    }
    w <- 1 / iv$rSe^2
    ok <- is.finite(w) & w > 0
    if (!any(ok)) {
        return(NULL)
    }
    .gcByVariant(
        iv$variantIds[ok],
        w[ok] * iv$seY[ok] / iv$betaX[ok] / sum(w[ok]),
        iv$z[ok]
    )
}

# The credible-set-aware estimate: each set's PIP-weighted composite ratio,
# pooled across sets by inverse variance, as coefficients on z.
# @noRd
.gcCsAwareCoefficients <- function(fmrEntry, gdf, cpipCutoff, alleleFlip) {
    inst <- .cipCsAwareInstruments(fmrEntry, gdf, alleleFlip)
    if (is.null(inst)) {
        return(NULL)
    }
    comp <- .cipCsComposites(inst, cpipCutoff)
    w <- 1 / comp$sbhat^2
    ok <- is.finite(w) & w > 0
    if (!any(ok)) {
        return(NULL)
    }
    pieces <- list_rbind(map2(
        comp$cs[ok],
        w[ok] / sum(w[ok]),
        .gcCsCoefficients,
        inst = inst
    ))
    .gcByVariant(pieces$variantId, pieces$coef, pieces$z)
}

# Credible set `cs`'s share `share` of the pooled estimate, spread over its
# variants as their composite ratio weights them.
# @noRd
.gcCsCoefficients <- function(cs, share, inst) {
    inCs <- which(inst$cs == cs)
    pNorm <- inst$pip[inCs] / sum(inst$pip[inCs])
    tibble(
        variantId = inst$variantIds[inCs],
        coef = share * pNorm * inst$sbhatY[inCs] / inst$bhatX[inCs],
        z = inst$z[inCs]
    )
}

# Coefficients summed over a variant that enters more than once.
# @noRd
.gcByVariant <- function(variantIds, coef, z) {
    ids <- unique(variantIds)
    list(
        variantIds = ids,
        coef = as.numeric(rowsum(coef, variantIds, reorder = FALSE)),
        z = z[match(ids, variantIds)]
    )
}

# @noRd
.gcAbortMismatch <- function(labels) {
    abort(c(
        glue(
            "geneCorrelation: the statistic of {str_flatten(head(labels, 3L), ",
            "'; ')}{if (length(labels) > 3L) '; ...' else ''} does not match ",
            "the one these inputs give."
        ),
        i = str_c(
            "Pass the GWAS, weights, `mrParam` and `alleleFlip` `x` was ",
            "computed with."
        )
    ))
}

# One GWAS study's genes `genes`, each one's variance (w'Rw, or a'Ra), and
# the signed correlation of each pair that shares a row, i < j.
# @noRd
.gcStudyPairs <- function(genes, gStudy, gwasSumStats, analysis) {
    blockRows <- sort(unique(list_c(map(genes, "rows"))))
    # A pair whose genes share several rows is computed once per row, the same
    # each time; the first is kept.
    computed <- list_rbind(map(
        blockRows,
        .gcRowCovariances,
        genes = genes,
        gwasSumStats = gwasSumStats,
        gStudy = gStudy
    ))
    pairs <- computed[!duplicated(computed[c("i", "j")]), ]
    diagonal <- filter(pairs, .data$i == .data$j)
    variance <- diagonal$cov[match(seq_along(genes), diagonal$i)]
    .gcCheckStatistics(genes, variance, analysis)
    off <- filter(pairs, .data$i < .data$j)
    list(
        genes = genes,
        variance = variance,
        i = off$i,
        j = off$j,
        rho = off$cov / sqrt(variance[off$i] * variance[off$j])
    )
}

# A study's correlation matrix on `scale` (.gcStudyPairs()).
# @noRd
.gcScaledMatrix <- function(study, scale) {
    n <- length(study$genes)
    ids <- map_chr(study$genes, "gene")
    Matrix::sparseMatrix(
        i = c(seq_len(n), study$i),
        j = c(seq_len(n), study$j),
        x = c(rep(1, n), .gcScaleCorrelation(study$rho, scale)),
        dims = c(n, n),
        dimnames = list(ids, ids),
        symmetric = TRUE
    )
}

# The covariances of the genes whose variants fall in GWAS row `r`, every pair
# with i <= j, over the LD reference of the rows they read.
# @noRd
.gcRowCovariances <- function(r, genes, gwasSumStats, gStudy) {
    members <- which(map_lgl(genes, .gcInRow, r = r))
    ref <- .cipRowsLdRef(
        gwasSumStats,
        sort(unique(list_c(map(genes[members], "rows")))),
        gStudy
    )
    ids <- unique(list_c(map(genes[members], "variantIds")))
    W <- exec(cbind, !!!map(genes[members], .gcColumn, ids = ids))
    cov <- .ldCrossForms(.gcLdSource(ref, ids), W)
    upper <- which(upper.tri(cov, diag = TRUE), arr.ind = TRUE)
    tibble(
        i = members[upper[, 1L]],
        j = members[upper[, 2L]],
        cov = cov[upper]
    )
}

# @noRd
.gcInRow <- function(gene, r) {
    is_in(r, gene$rows)
}

# A gene's coefficients as a column over the variants `ids`.
# @noRd
.gcColumn <- function(gene, ids) {
    replace(numeric(length(ids)), match(gene$variantIds, ids), gene$coef)
}

# The LD over `ids` as the pipeline reads it for a TWAS z (.cipTwasZ()): a
# genotype reference's factor, or an LdMatrix's packed block.
# @noRd
.gcLdSource <- function(ldReference, ids) {
    if (.ldIsGenotypeRef(ldReference)) {
        return(list(
            ldFactor = .ldFactorFromRef(ldReference, ids, "geneCorrelation")
        ))
    }
    .ldPackedFromRef(ldReference, ids, "geneCorrelation")
}

# Each gene's statistic from its coefficients -- the TWAS z, w'z /
# sqrt(w'Rw), or the MR estimate, a'z -- against the one `x` holds.
# @noRd
.gcCheckStatistics <- function(genes, variance, analysis) {
    value <- map_dbl(genes, .gcLinearValue)
    statistic <- if (analysis == "twas") value / sqrt(variance) else value
    target <- map_dbl(genes, "target")
    bad <- abs(statistic - target) > 1e-6 * pmax(1, abs(target))
    if (any(bad)) {
        .gcAbortMismatch(map_chr(genes[bad], "label"))
    }
    invisible(NULL)
}

# @noRd
.gcLinearValue <- function(gene) {
    sum(gene$coef * gene$z)
}

# The correlation of the statistics on `scale` from that of the z-scores.
# @noRd
.gcScaleCorrelation <- function(rho, scale) {
    r <- pmax(pmin(rho, 1), -1)
    switch(
        scale,
        signed = r,
        squared = r^2,
        absolute = (sqrt(1 - r^2) + r * asin(r) - 1) / (pi / 2 - 1),
        oneSided = .oneSidedZCorrelation(r)
    )
}

# ---------- The one-sided z-score's correlation ----------
#
# The one-sided z of a z-score's two-sided p-value, g(z) = Phi^-1(1 - p), is
# standard normal under the null. For z-scores of correlation rho, Mehler's
# expansion gives cor(g(z1), g(z2)) = sum_k b_k^2 rho^k, over the coefficients
# b_k = E[g(Z) h_k(Z)] of g on the normalized Hermite polynomials h_k; g is
# even, so only even k enter, and the b_k^2 sum to Var g(Z) = 1. The b_k^2
# decay slowly (as a power of k, from g's log singularity at 0), so the
# series is summed to a fixed order and the rest is a power law fitted to
# its upper half, scaled so the whole sums to 1.

# g(z), read off the chi-square on the log scale, from whichever tail keeps
# its precision.
# @noRd
.zOneSided <- function(z) {
    q <- z^2
    ifelse(
        q < 1,
        stats::qnorm(stats::pchisq(q, 1, log.p = TRUE), log.p = TRUE),
        stats::qnorm(
            stats::pchisq(q, 1, lower.tail = FALSE, log.p = TRUE),
            lower.tail = FALSE,
            log.p = TRUE
        )
    )
}

# cor(g(z1), g(z2)) for z-scores of correlation `rho`.
# @noRd
.oneSidedZCorrelation <- function(rho) {
    series <- .oneSidedHermiteSeries()
    r <- pmin(abs(rho), 1)
    body <- as.numeric(outer(r, series$k, "^") %*% series$b2) +
        .oneSidedTail(r, series)
    ifelse(r >= 1, 1, ifelse(r <= 0, 0, body))
}

# The squared coefficients b_k^2 of the even orders up to `maxOrder`, and the
# exponent alpha of the power law b_k^2 ~ k^-alpha fitted to the upper half
# of them.
# @noRd
.oneSidedHermiteSeries <- function(maxOrder = 200L) {
    k <- seq(2L, maxOrder, by = 2L)
    # The coefficients' integrands vanish past here: h_k is negligible beyond
    # 2 sqrt(k) standard deviations, and the normal density well before.
    upper <- 2 * sqrt(maxOrder) + 12
    b2 <- map_dbl(k, .oneSidedHermiteCoef, upper = upper)^2
    upperHalf <- k > maxOrder / 2
    fit <- stats::lm.fit(cbind(1, log(k[upperHalf])), log(b2[upperHalf]))
    list(
        k = k,
        b2 = b2,
        alpha = -fit$coefficients[[2L]],
        maxOrder = maxOrder
    )
}

# b_k = E[g(Z) h_k(Z)], twice the integral over the positive half (k even).
# @noRd
.oneSidedHermiteCoef <- function(k, upper) {
    stats::integrate(
        .oneSidedHermiteIntegrand,
        0,
        upper,
        k = k,
        subdivisions = 2000L,
        rel.tol = 1e-10
    )$value
}

# @noRd
.oneSidedHermiteIntegrand <- function(x, k) {
    2 * stats::dnorm(x) * .zOneSided(x) * .hermiteNormalized(x, k)
}

# h_k(x) = He_k(x) / sqrt(k!), by its three-term recurrence, which stays on
# the scale of the values where He_k overflows.
# @noRd
.hermiteNormalized <- function(x, k) {
    reduce(
        seq_len(k - 1L),
        .hermiteStep,
        x = x,
        .init = list(previous = rep(1, length(x)), current = x)
    )$current
}

# @noRd
.hermiteStep <- function(state, j, x) {
    list(
        previous = state$current,
        current = (x * state$current - sqrt(j) * state$previous) / sqrt(j + 1)
    )
}

# The series past its last order: sum_k c k^-alpha r^k over the even k
# beyond it, as an integral (an upper incomplete gamma), with c set so that
# it is the mass the computed orders leave at r = 1. With u = -K log r and
# a = 1 - alpha in (-1, 0), Gamma(a, u) = (Gamma(a + 1, u) - u^a e^-u) / a.
# @noRd
.oneSidedTail <- function(r, series) {
    u <- -log(r) * series$maxOrder
    a <- 1 - series$alpha
    upperGamma <- (gamma(a + 1) *
        stats::pgamma(u, a + 1, lower.tail = FALSE) -
        u^a * exp(-u)) /
        a
    (1 - sum(series$b2)) *
        (series$alpha - 1) *
        u^(series$alpha - 1) *
        upperGamma
}

# ---------- Whitening by the gene correlation ----------

# Gene-level response `y` and features `X` whitened by the genes'
# correlation, as PoPS fits its regression: with the correlation's Cholesky
# factor L (computed per independent block), L^-1 y and L^-1 X, with the
# whitened intercept L^-1 1 projected out of both, and the whitened gene
# covariates out of y only. `regularize` adds PoPS's ridge to the diagonal
# first (.glsRidge()). `y` is named by gene; `X` and `covariates` have a row
# for each of its genes, and `correlation` a row and column. The whitened
# intercept is returned too, for an engine that takes covariates.
# @noRd
.glsWhiten <- function(
    y,
    X,
    correlation,
    covariates = NULL,
    regularize = TRUE
) {
    genes <- .glsGenes(y, X, correlation, covariates)
    sigma <- correlation[genes, genes]
    blocks <- .glsBlocks(sigma)
    ridge <- if (regularize) .glsRidge(sigma, blocks, y[genes]) else 0
    factors <- map(blocks, .glsBlockFactor, sigma = sigma, ridge = ridge)
    intercept <- .glsApply(matrix(1, length(genes)), blocks, factors)
    yW <- .glsApply(matrix(y[genes]), blocks, factors)
    nuisance <- if (is.null(covariates)) {
        intercept
    } else {
        cbind(
            intercept,
            .glsApply(covariates[genes, , drop = FALSE], blocks, factors)
        )
    }
    XW <- .glsApply(X[genes, , drop = FALSE], blocks, factors)
    list(
        y = set_names(as.numeric(qr.resid(qr(nuisance), yW)), genes),
        intercept = intercept,
        X = `dimnames<-`(
            qr.resid(qr(intercept), XW),
            list(genes, colnames(X))
        )
    )
}

# The genes of `y`, which every other input must cover.
# @noRd
.glsGenes <- function(y, X, correlation, covariates) {
    genes <- names(y)
    if (is.null(genes) || anyDuplicated(genes) > 0L) {
        abort("`y` must be named by gene, each gene once.")
    }
    .glsRequireGenes(genes, rownames(correlation), "the gene correlation")
    .glsRequireGenes(genes, rownames(X), "`X`")
    if (!is.null(covariates)) {
        .glsRequireGenes(genes, rownames(covariates), "`covariates`")
    }
    genes
}

# @noRd
.glsRequireGenes <- function(genes, have, what) {
    absent <- setdiff(genes, have)
    if (length(absent) > 0L) {
        abort(glue(
            "{length(absent)} gene(s) of `y` are not in {what}: ",
            "{str_flatten(head(absent, 5L), ', ')}",
            "{if (length(absent) > 5L) ', ...' else ''}."
        ))
    }
    invisible(NULL)
}

# The correlation's independent blocks: the connected components of its
# nonzero pattern, as index vectors.
# @noRd
.glsBlocks <- function(sigma) {
    entries <- Matrix::summary(methods::as(sigma, "TsparseMatrix"))
    linked <- entries[entries$x != 0 & entries$i != entries$j, ]
    graph <- igraph::make_graph(
        as.vector(rbind(linked$i, linked$j)),
        n = nrow(sigma),
        directed = FALSE
    )
    unname(split(seq_len(nrow(sigma)), igraph::components(graph)$membership))
}

# PoPS's ridge on the correlation's diagonal (regularize_error_cov()): enough
# to lift its smallest eigenvalue to zero, 0.05 more, and 0.9 of the response's
# variance beyond the unit variance the correlation describes, which a null
# correlation leaves out.
# @noRd
.glsRidge <- function(sigma, blocks, y) {
    smallest <- min(0, map_dbl(blocks, .glsMinEigen, sigma = sigma))
    abs(smallest) + 0.05 + 0.9 * max(0, stats::var(y) - 1)
}

# @noRd
.glsMinEigen <- function(idx, sigma) {
    block <- as.matrix(sigma[idx, idx, drop = FALSE])
    min(eigen(block, symmetric = TRUE, only.values = TRUE)$values)
}

# The upper Cholesky factor of one block of the (ridged) correlation.
# @noRd
.glsBlockFactor <- function(idx, sigma, ridge) {
    chol(as.matrix(sigma[idx, idx, drop = FALSE]) + diag(ridge, length(idx)))
}

# L^-1 M, block by block, its rows in the genes' order.
# @noRd
.glsApply <- function(M, blocks, factors) {
    pieces <- map2(blocks, factors, .glsBlockSolve, M = M)
    exec(rbind, !!!pieces)[order(list_c(blocks)), , drop = FALSE]
}

# @noRd
.glsBlockSolve <- function(idx, factor, M) {
    backsolve(
        factor,
        as.matrix(M[idx, , drop = FALSE]),
        transpose = TRUE
    )
}
