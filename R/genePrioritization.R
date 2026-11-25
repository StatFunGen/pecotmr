#' Prioritize genes by their features (PoPS)
#'
#' @description Scores each gene by how well its features predict its TWAS,
#'   MR or variance-component association, as the Polygenic Priority Score
#'   (PoPS; Weeks et al. 2023) does with MAGMA's gene z-scores. The genes'
#'   statistics from a \code{\link{causalInferencePipeline}} or
#'   \code{\link{geneVarianceComponentTest}} result are regressed on their
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
#'   its logarithm; for variance components, the number of variants tested,
#'   the effective number of eigenvalues of the statistic's null distribution
#'   (\eqn{(\sum_k \lambda_k)^2 / \sum_k \lambda_k^2}), the logarithms of
#'   both, and for squared TWAS weights their cross-validated \eqn{r^2}. A
#'   covariate missing for any gene (a CV \eqn{r^2} the weights do not
#'   record, say) is left out, with a message. They play the part of MAGMA's
#'   gene size and density in PoPS.
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
#'   \code{\link{causalInferencePipeline}},
#'   \code{\link{geneVarianceComponentTest}},
#'   \code{\link{geneSetEnrichment}}
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
    analysis = c("twas", "mr", "varianceComponent"),
    method = NULL,
    context = NULL,
    scale = NULL,
    fitMethod = "ridge",
    methodArgs = NULL,
    covariates = TRUE,
    excludeHla = TRUE,
    mrParam = MrParam(),
    alleleFlip = TRUE,
    credibleSetParam = CredibleSetParam(),
    seed = NULL
) {
    cfg <- .gpConfig(
        features,
        analysis,
        scale,
        fitMethod,
        methodArgs,
        credibleSetParam,
        !missing(credibleSetParam),
        covariates,
        excludeHla,
        seed
    )
    studies <- .gcCompute(
        x,
        gwasSumStats,
        twasWeights,
        fineMappingResult,
        cfg$analysis,
        method,
        context,
        mrParam,
        alleleFlip
    )
    .gpScoreStudies(studies, cfg, twasWeights %||% fineMappingResult)
}

# The run's settings, checked: the analysis, its scale (its default when
# NULL) and the fit method chosen; the genes of `features` and
# their features; the fit method's options; and the rest of
# genePrioritization()'s arguments the fits read.
# @noRd
.gpConfig <- function(
    features,
    analysis,
    scale,
    fitMethod,
    methodArgs,
    credibleSetParam,
    supplied,
    covariates,
    excludeHla,
    seed
) {
    chosen <- arg_match(analysis, c("twas", "mr", "varianceComponent"))
    engine <- arg_match(fitMethod, names(.gpFitMethods))
    .gpAssertArgs(
        covariates,
        excludeHla,
        seed,
        engine,
        credibleSetParam,
        supplied
    )
    list(
        genes = .gpFeatureGenes(features, excludeHla),
        features = SummarizedExperiment::assay(features),
        fitMethod = engine,
        options = .gpOptions(engine, methodArgs, credibleSetParam),
        scale = .gcResolveScale(scale, chosen),
        covariates = covariates,
        analysis = chosen,
        seed = seed
    )
}

# Each study's scores (.gpStudy()), under `cfg`'s seed when it has one. The
# CV r^2 covariate is read off the TWAS weights `source`.
# @noRd
.gpScoreStudies <- function(studies, cfg, source) {
    rsq <- .gpCvRsqIfRead(studies, cfg$analysis, cfg$covariates, source)
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
.gpFeatureGenes <- function(
    features,
    excludeHla,
    caller = "genePrioritization"
) {
    assertClass(features, "RangedSummarizedExperiment")
    ids <- rownames(features)
    if (is.null(ids) || anyDuplicated(ids) > 0L) {
        abort(glue(
            "{caller}: `features` must have its genes as row names, each ",
            "gene once."
        ))
    }
    if (length(SummarizedExperiment::assays(features)) == 0L) {
        abort(glue("{caller}: `features` has no assay of features."))
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

# The CV r^2 of `source`'s weight sets (.gpCvRsq()) when the covariates read
# it: those of TWAS z-scores and of squared TWAS weights.
# @noRd
.gpCvRsqIfRead <- function(studies, analysis, covariates, source) {
    weights <- map_chr(
        list_flatten(map(studies, "genes")),
        .gpWeightKind,
        analysis = analysis
    )
    if (covariates && any(is_in(weights, c("twas", "squaredWeight")))) {
        .gpCvRsq(source)
    }
}

# A gene's kind of weight: its analysis, or its variant weights.
# @noRd
.gpWeightKind <- function(gene, analysis) {
    gene$variantWeights %||% analysis
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
.gpResponse <- function(
    studyGenes,
    scale,
    featureIds,
    caller = "genePrioritization"
) {
    ids <- map_chr(studyGenes, "gene")
    z <- map_dbl(studyGenes, "statistic")
    value <- switch(
        scale,
        signed = z,
        absolute = abs(z),
        squared = z^2,
        oneSided = .zOneSided(z),
        standardized = z,
        probit = map_dbl(studyGenes, "vcZ")
    )
    known <- is_in(ids, featureIds)
    if (!all(known)) {
        inform(glue(
            "{caller}: {sum(!known)} gene(s) with a statistic are not in ",
            "`features` and are left out."
        ))
    }
    set_names(value[known], ids[known])
}

# The gene covariates (see genePrioritization()'s Details), a row per gene;
# a column missing for any gene is left out, with a message.
# @noRd
.gpCovariates <- function(
    study,
    analysis,
    rsq,
    caller = "genePrioritization"
) {
    table <- switch(
        analysis,
        twas = .gpTwasCovariates(study, rsq),
        mr = .gpMrCovariates(study),
        varianceComponent = .gpVcCovariates(study, rsq)
    )
    complete <- apply(is.finite(table), 2L, all)
    if (!all(complete)) {
        inform(glue(
            "{caller}: covariate(s) ",
            "{str_flatten(colnames(table)[!complete], ', ')} are missing ",
            "for some genes and are left out."
        ))
    }
    if (!any(complete)) {
        return(NULL)
    }
    `rownames<-`(
        table[, complete, drop = FALSE],
        map_chr(study$genes, "gene")
    )
}

# @noRd
.gpTwasCovariates <- function(study, rsq) {
    nonzero <- map_int(study$genes, .gpNonzeroWeights)
    cbind(
        nVariants = nonzero,
        logNVariants = log(nonzero),
        weightVariance = study$variance,
        logWeightVariance = log(study$variance),
        cvRsq = unname(rsq[map_chr(study$genes, "cellKey")])
    )
}

# @noRd
.gpMrCovariates <- function(study) {
    instruments <- map_dbl(study$genes, "nIV")
    cbind(nInstruments = instruments, logNInstruments = log(instruments))
}

# The variants tested and the effective number of eigenvalues,
# (sum lambda)^2 / sum lambda^2, with sum lambda^2 the study's variance
# (.gcStudyPairs()); and squared TWAS weights' CV r^2.
# @noRd
.gpVcCovariates <- function(study, rsq) {
    tested <- map_int(study$genes, .gpNonzeroWeights)
    effective <- map_dbl(study$genes, "vcMean")^2 / study$variance
    squared <- map_chr(study$genes, "variantWeights") == "squaredWeight"
    cbind(
        nVariants = tested,
        logNVariants = log(tested),
        effectiveDf = effective,
        logEffectiveDf = log(effective),
        cvRsq = if (all(squared)) {
            unname(rsq[map_chr(study$genes, "cellKey")])
        }
    )
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
# Gene-set enrichment
# =============================================================================

#' Gene-set enrichment of genes' TWAS, MR or variance-component statistics
#'
#' @description Tests whether the genes of each set have larger statistics
#'   than other genes, adjusted for the gene covariates and for the
#'   correlation the LD between their variants induces: MAGMA's competitive
#'   gene-set test, a generalized least squares regression of the genes'
#'   statistics on the set's indicator, for the TWAS z-scores or MR
#'   estimates of a \code{\link{causalInferencePipeline}} result or the
#'   variance components of a \code{\link{geneVarianceComponentTest}}
#'   result.
#'
#' @details For each GWAS study, the genes of \code{features} with a
#'   statistic are whitened by their correlation
#'   (\code{\link{geneCorrelation}}), with PoPS's ridge on its diagonal as
#'   \code{\link{genePrioritization}} whitens them. The whitened intercept,
#'   and with \code{covariates = TRUE} the gene covariates
#'   \code{\link{genePrioritization}} describes, are projected out of both
#'   the statistics and each set's column, and the set's effect is their
#'   least-squares slope. Its standard error takes the residual variance from
#'   the fit, so association in excess of the correlation's unit variance (a
#'   polygenic trait's) widens it rather than inflating the test. The
#'   p-value is one-sided, the upper tail of the \eqn{t} statistic on the
#'   residual degrees of freedom: a set is enriched when its genes'
#'   statistics are larger. Each set is tested on its own; a column may be
#'   membership (0/1) or a continuous gene annotation. Genes of the HLA
#'   region are kept.
#'
#' @inheritParams genePrioritization
#' @param features A \code{RangedSummarizedExperiment} of genes by gene
#'   sets: its rows the genes, named as \code{x}'s \code{trait}; its first
#'   assay the sets, with column names, binary or continuous, dense or
#'   sparse.
#' @return A list with one element per GWAS study, named by study: a tibble
#'   with a row per set and columns \code{set}; \code{nGenes}, the set's
#'   genes with a statistic (nonzero entries); \code{beta}, \code{se},
#'   \code{t} and \code{df}; and \code{pval}, one-sided. A set with no
#'   variation left once the covariates are projected out (every gene, or
#'   none) has \code{NA} for its effect.
#' @seealso \code{\link{geneCorrelation}}, \code{\link{genePrioritization}},
#'   \code{\link{geneVarianceComponentTest}}
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
#' gwas <- GwasSumStats(
#'     "gwas", list(gr), "hg19",
#'     ldReference = LdMatrix(
#'         diag(200),
#'         variants = data.frame(
#'             chrom = chr, pos = pos, SNP = ids, A1 = "T", A2 = "C"
#'         )
#'     ),
#'     qcOptions = list(qc = "done")
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
#' sets <- SummarizedExperiment::SummarizedExperiment(
#'     list(sets = matrix(rbinom(40 * 3, 1, 0.3), 40,
#'         dimnames = list(genes, sprintf("set%d", 1:3)))),
#'     rowRanges = GenomicRanges::GRanges(
#'         chr[first],
#'         IRanges::IRanges(pos[first], width = 500),
#'         names = genes
#'     )
#' )
#' geneSetEnrichment(res, sets, gwas, twasWeights = tw, covariates = FALSE)
#' @importFrom checkmate assertFlag
#' @export
geneSetEnrichment <- function(
    x,
    features,
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    analysis = c("twas", "mr", "varianceComponent"),
    method = NULL,
    context = NULL,
    scale = NULL,
    covariates = TRUE,
    mrParam = MrParam(),
    alleleFlip = TRUE
) {
    analysis <- arg_match(analysis)
    assertFlag(covariates)
    cfg <- .gseConfig(
        features,
        .gcResolveScale(scale, analysis),
        covariates,
        analysis
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
    rsq <- .gpCvRsqIfRead(
        studies,
        analysis,
        covariates,
        twasWeights %||% fineMappingResult
    )
    map(studies, .gseStudy, cfg = cfg, rsq = rsq)
}

# The run's settings: the genes of `features` and their sets, checked, and
# the scale, covariates and analysis the tests read.
# @noRd
.gseConfig <- function(features, scale, covariates, analysis) {
    genes <- .gpFeatureGenes(features, FALSE, "geneSetEnrichment")
    sets <- SummarizedExperiment::assay(features)
    if (is.null(colnames(sets)) || anyDuplicated(colnames(sets)) > 0L) {
        abort(str_c(
            "geneSetEnrichment: `features` must name its sets (its assay's ",
            "columns), each set once."
        ))
    }
    list(
        ids = genes$ids,
        sets = sets,
        scale = scale,
        covariates = covariates,
        analysis = analysis
    )
}

# One GWAS study's tests, a set at a time in chunks of the sets.
# @noRd
.gseStudy <- function(study, cfg, rsq) {
    response <- .gpResponse(
        study$genes,
        cfg$scale,
        cfg$ids,
        "geneSetEnrichment"
    )
    nuisance <- if (cfg$covariates) {
        .gpCovariates(study, cfg$analysis, rsq, "geneSetEnrichment")
    }
    fit <- .gseNullFit(
        response,
        .gcScaledMatrix(study, cfg$scale),
        nuisance
    )
    starts <- seq.int(1L, ncol(cfg$sets), by = .gseChunkSize())
    list_rbind(map(starts, .gseChunk, sets = cfg$sets, fit = fit))
}

# How many sets are whitened together.
# @noRd
.gseChunkSize <- function() 256L

# The model without a set: the statistics `y` whitened by the ridged
# correlation (.glsWhiten()'s), and the whitened intercept and covariates
# (the nuisance) projected out of them; the whitening and the nuisance's QR
# decomposition, for the sets'; and the residual degrees of freedom once a
# set is added.
# @noRd
.gseNullFit <- function(y, correlation, covariates) {
    genes <- names(y)
    .glsRequireGenes(genes, rownames(correlation), "the gene correlation")
    sigma <- correlation[genes, genes]
    blocks <- .glsBlocks(sigma)
    ridge <- .glsRidge(sigma, blocks, y)
    factors <- map(blocks, .glsBlockFactor, sigma = sigma, ridge = ridge)
    nuisance <- .glsApply(
        cbind(rep(1, length(genes)), covariates[genes, , drop = FALSE]),
        blocks,
        factors
    )
    basis <- qr(nuisance)
    df <- length(genes) - basis$rank - 1L
    if (df < 1L) {
        abort(glue(
            "geneSetEnrichment: {length(genes)} gene(s) with a statistic ",
            "leave no residual degrees of freedom for a set's effect."
        ))
    }
    list(
        genes = genes,
        blocks = blocks,
        factors = factors,
        basis = basis,
        y = qr.resid(basis, .glsApply(matrix(y), blocks, factors)),
        df = df
    )
}

# The tests of the sets from column `start` on, a chunk of them: each set's
# whitened column with the nuisance projected out, regressed on alone.
# @noRd
#' @importFrom stats pt
.gseChunk <- function(start, sets, fit) {
    cols <- seq.int(start, min(start + .gseChunkSize() - 1L, ncol(sets)))
    raw <- as.matrix(sets[fit$genes, cols, drop = FALSE])
    whitened <- .glsApply(raw, fit$blocks, fit$factors)
    X <- qr.resid(fit$basis, whitened)
    sxx <- colSums(X^2)
    sxy <- as.numeric(crossprod(X, fit$y))
    # A set the nuisance spans (every gene, or none) has nothing left.
    identifiable <- sxx > sqrt(.Machine$double.eps) * colSums(whitened^2)
    beta <- if_else(identifiable, sxy / sxx, NA_real_)
    se <- sqrt((sum(fit$y^2) - sxy * beta) / fit$df / sxx)
    tibble(
        set = colnames(sets)[cols],
        nGenes = as.integer(colSums(raw != 0)),
        beta = beta,
        se = se,
        t = beta / se,
        df = fit$df,
        pval = pt(beta / se, fit$df, lower.tail = FALSE)
    )
}

# =============================================================================
# Variance-component test of each gene's QTL-weighted variants
# =============================================================================

#' Variance-component test of each gene's QTL-weighted variants
#'
#' @description For each QTL weight set and each GWAS study, tests whether
#'   the GWAS associations at the gene's cis variants exceed their spread
#'   under the null, each variant weighted by how strongly the QTL analysis
#'   ties it to the gene. It is the score test of a local variance component
#'   whose per-variant effect variance is proportional to that weight.
#'   MAGMA's SNP-wise mean test is the case of equal weights over a window of
#'   positions; this one takes its weights from the QTL fits instead.
#'
#' @details With the GWAS z-scores \eqn{z} at the gene's variants, their
#'   weights \eqn{a} and the GWAS LD \eqn{R}, the statistic is
#'   \eqn{Q = \sum_j a_j z_j^2 = z' A z}. Under the null,
#'   \eqn{z \sim N(0, R)}, \eqn{Q} is distributed as
#'   \eqn{\sum_k \lambda_k \chi^2_1} over the eigenvalues \eqn{\lambda_k} of
#'   \eqn{A^{1/2} R A^{1/2}}. Its p-value is computed exactly by
#'   Farebrother's algorithm (\pkg{CompQuadForm}) and, below \eqn{10^{-12}}
#'   where that loses its precision, by Kuonen's saddlepoint approximation
#'   (\pkg{survey}). Two gene-level statistics are reported for a regression
#'   on gene features: the probit \eqn{\Phi^{-1}(1 - p)}, MAGMA's scale, and
#'   the standardized \eqn{(Q - E Q) / SD(Q)}, with \eqn{E Q = \sum_k
#'   \lambda_k} and \eqn{Var(Q) = 2 \sum_k \lambda_k^2}. Under the
#'   variance-component model the standardized statistic's mean is linear in
#'   the gene's weighted local heritability, and two genes' standardized
#'   statistics correlate exactly as \code{\link{geneCorrelation}} gives.
#'
#'   \code{variantWeights} chooses \eqn{a}:
#'   \describe{
#'     \item{\code{"pip"}}{Each variant's posterior inclusion probability
#'       from the fine-mapping row, over every variant the row holds (the
#'       gene's cis window).}
#'     \item{\code{"flooredPip"}}{\code{pipFloor} plus the PIP, which gives
#'       every variant of the window some weight, as MAGMA's window does.}
#'     \item{\code{"squaredWeight"}}{The squared TWAS weights \eqn{w_j^2},
#'       the sum-of-squares counterpart of the TWAS z-score
#'       \eqn{w' z / \sqrt{w' R w}}: from \code{twasWeights}, else the
#'       \code{fineMappingResult}'s posterior effects, every weight set
#'       (method) either holds.}
#'   }
#'   The weights are matched to the GWAS variants as
#'   \code{\link{causalInferencePipeline}} matches TWAS weights, and only
#'   variants with a positive weight enter. Each weight set is tested once
#'   per GWAS study, over the study's rows its variants overlap, with the LD
#'   reference those rows share.
#'
#'   Each test takes the eigenvalues of a \eqn{p \times p} matrix, \eqn{p}
#'   the gene's variants with a weight; with a genotype reference of \eqn{m}
#'   samples, of an \eqn{m \times m} one when \eqn{m < p}.
#'
#' @param gwasSumStats A QC'd \code{\link{GwasSumStats}}.
#' @param twasWeights,fineMappingResult The QTL fits: \code{"pip"} and
#'   \code{"flooredPip"} read the \code{fineMappingResult} (and refuse
#'   \code{twasWeights}, which they do not read); \code{"squaredWeight"} reads
#'   \code{twasWeights}, else the \code{fineMappingResult}.
#' @param variantWeights The variant weights: \code{"pip"} (default),
#'   \code{"squaredWeight"} or \code{"flooredPip"}; see Details.
#' @param pipFloor For \code{"flooredPip"} only, the weight added to every
#'   variant's PIP: a positive number.
#' @param alleleFlip Logical, default \code{TRUE}: match variants to the GWAS
#'   with ref/alt swaps recognized, as \code{\link{causalInferencePipeline}}
#'   does.
#' @return A \code{GRanges}, one range per weight set and GWAS study, at the
#'   span of the variants tested, with columns \code{qtlStudy},
#'   \code{context}, \code{trait}, \code{method} (and the joint-fit columns
#'   when a weight set came from a joint fit), \code{gwasStudy},
#'   \code{analysis} (\code{"varianceComponent"}), \code{variantWeights},
#'   \code{pipFloor} (\code{NA} unless floored), \code{nVariants} (the
#'   variants tested), \code{vcStatistic} (\eqn{Q}), \code{vcMean} and
#'   \code{vcSd} (its null mean and standard deviation),
#'   \code{vcStandardized}, \code{vcPval} and \code{vcZ} (the probit of the
#'   p-value). Its \code{seqinfo()} names the inputs' genome build.
#' @seealso \code{\link{geneCorrelation}}, \code{\link{genePrioritization}},
#'   \code{\link{geneSetEnrichment}}
#' @examples
#' ids <- sprintf("chr1:%d:C:T", 100L * seq_len(8))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L * 1:8, width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'     SNP = ids, A1 = "T", A2 = "C", Z = c(3, 2, 1, 0, -1, 0.5, 1, 2)
#' )
#' gwas <- GwasSumStats(
#'     "gwas", list(gr), "hg19",
#'     ldReference = LdMatrix(
#'         0.5^abs(outer(1:8, 1:8, "-")),
#'         variants = data.frame(
#'             chrom = "chr1", pos = 100L * seq_len(8), SNP = ids,
#'             A1 = "T", A2 = "C"
#'         )
#'     ),
#'     qcOptions = list(qc = "done")
#' )
#' tw <- TwasWeights(
#'     rep("qtl", 2), rep("ctx", 2), c("geneA", "geneB"), rep("lasso", 2),
#'     list(
#'         TwasWeightsRow(ids[1:5], c(0.4, 0.2, 0, 0.1, 0.3)),
#'         TwasWeightsRow(ids[3:8], c(0.2, 0.3, 0.1, 0, 0.2, 0.1))
#'     )
#' )
#' if (requireNamespace("CompQuadForm", quietly = TRUE)) {
#'     geneVarianceComponentTest(
#'         gwas,
#'         twasWeights = tw,
#'         variantWeights = "squaredWeight"
#'     )
#' }
#' @importFrom checkmate assertFlag
#' @export
geneVarianceComponentTest <- function(
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    variantWeights = c("pip", "squaredWeight", "flooredPip"),
    pipFloor = NULL,
    alleleFlip = TRUE
) {
    variantWeights <- arg_match(variantWeights)
    assertFlag(alleleFlip)
    .vcAssertFloor(variantWeights, pipFloor)
    build <- .inputsBuild(
        list(gwasSumStats, twasWeights, fineMappingResult),
        "geneVarianceComponentTest"
    )
    source <- .vcWeightSource(
        variantWeights,
        twasWeights,
        fineMappingResult,
        gwasSumStats
    )
    keys <- .weightCellIds(source)
    entries <- compact(map(
        seq_len(nrow(keys)),
        .vcEntry,
        keys = keys,
        source = source,
        variantWeights = variantWeights,
        pipFloor = pipFloor
    ))
    .checkGwasBlockReach(
        set_names(map(entries, "variants"), map_chr(entries, .cipTupleLabel)),
        gwasSumStats,
        "geneVarianceComponentTest",
        "QTL weight sets"
    )
    gwasIndex <- .rowSpans(gwasSumStats)
    records <- list_flatten(map(
        entries,
        .vcScoreEntry,
        gwasSumStats = gwasSumStats,
        gwasIndex = gwasIndex,
        alleleFlip = alleleFlip
    ))
    .withGenome(.vcResult(records), build)
}

# `pipFloor` is a positive number with "flooredPip", and absent otherwise.
# @noRd
#' @importFrom checkmate assertNumber
.vcAssertFloor <- function(variantWeights, pipFloor) {
    if (variantWeights != "flooredPip") {
        if (!is.null(pipFloor)) {
            abort(glue(
                "geneVarianceComponentTest: `pipFloor` applies to ",
                "`variantWeights = \"flooredPip\"` only, not ",
                "\"{variantWeights}\"."
            ))
        }
        return(invisible(NULL))
    }
    if (is.null(pipFloor)) {
        abort(str_c(
            "geneVarianceComponentTest: `variantWeights = \"flooredPip\"` ",
            "needs `pipFloor`, the weight every variant's PIP is floored by."
        ))
    }
    assertNumber(pipFloor, lower = 0, finite = TRUE)
    if (pipFloor == 0) {
        abort("geneVarianceComponentTest: `pipFloor` must be positive.")
    }
    invisible(NULL)
}

# The fits `variantWeights` reads, checked against the GWAS they are tested
# with: the fine-mapping result for PIPs; the TWAS weights, else the
# fine-mapping result, for squared weights.
# @noRd
#' @importFrom checkmate assertClass
.vcWeightSource <- function(
    variantWeights,
    twasWeights,
    fineMappingResult,
    gwasSumStats
) {
    assertClass(gwasSumStats, "GwasSumStats")
    assertClass(twasWeights, "TwasWeights", null.ok = TRUE)
    assertClass(fineMappingResult, "QtlFineMappingResult", null.ok = TRUE)
    if (length(qcOptions(gwasSumStats)) == 0L) {
        abort(str_c(
            "geneVarianceComponentTest: gwasSumStats has no QC record ",
            "(qcOptions() is empty). Call summaryStatsQc() first."
        ))
    }
    source <- .vcSourceFor(variantWeights, twasWeights, fineMappingResult)
    .requireOverlapLdRefs(
        source,
        gwasSumStats,
        pipelineName = "geneVarianceComponentTest",
        label = class(source)[[1L]]
    )
    source
}

# @noRd
.vcSourceFor <- function(variantWeights, twasWeights, fineMappingResult) {
    if (variantWeights == "squaredWeight") {
        source <- twasWeights %||% fineMappingResult
        if (is.null(source)) {
            abort(str_c(
                "geneVarianceComponentTest: squared weights need ",
                "`twasWeights`, or a `fineMappingResult`."
            ))
        }
        return(source)
    }
    if (is.null(fineMappingResult)) {
        abort(glue(
            "geneVarianceComponentTest: \"{variantWeights}\" weights are ",
            "PIPs, read from a `fineMappingResult`."
        ))
    }
    if (!is.null(twasWeights)) {
        abort(glue(
            "geneVarianceComponentTest: \"{variantWeights}\" weights are ",
            "read from the `fineMappingResult` alone; `twasWeights` is not ",
            "read, so leave it unset."
        ))
    }
    fineMappingResult
}

# Weight set `i` of `keys` (.weightCellIds()): its identity (.cipCellTuple()),
# its variants and their weights, and how those were chosen. NULL when it has
# none to test.
# @noRd
.vcEntry <- function(i, keys, source, variantWeights, pipFloor) {
    weighted <- if (variantWeights == "squaredWeight") {
        .vcSquaredWeights(i, keys, source)
    } else {
        .vcPips(i, keys, source, pipFloor %||% 0)
    }
    if (is.null(weighted)) {
        return(NULL)
    }
    list_assign(
        weighted,
        variantWeights = variantWeights,
        pipFloor = pipFloor %||% NA_real_
    )
}

# @noRd
.vcSquaredWeights <- function(i, keys, source) {
    input <- .cipTwasTuple(i, keys, source)
    if (is.null(input)) {
        return(NULL)
    }
    list(
        tuple = input$tuple,
        variants = input$variants,
        variantIds = input$weightsInfo$variantIds,
        weights = input$weightsInfo$weights^2
    )
}

# @noRd
.vcPips <- function(i, keys, source, floor) {
    input <- .cipMrTuple(i, keys, source)
    pip <- .fmrRowPip(.fmrRowParts(input$fmrEntry, 1L))
    if (length(pip) == 0L) {
        warn(glue(
            "geneVarianceComponentTest: {.cipTupleLabel(input)} holds no ",
            "PIPs; skipping."
        ))
        return(NULL)
    }
    list(
        tuple = input$tuple,
        variants = input$variants,
        variantIds = names(pip),
        weights = unname(pip) + floor
    )
}

# One weight set's tests, one per GWAS study its variants overlap.
# @noRd
.vcScoreEntry <- function(entry, gwasSumStats, gwasIndex, alleleFlip) {
    window <- range(entry$variants)
    compact(imap(
        .cipStudyRows(gwasIndex, gwasSumStats, window),
        .vcScorePair,
        entry = entry,
        gwasSumStats = gwasSumStats,
        window = window,
        alleleFlip = alleleFlip
    ))
}

# One weight set's test against GWAS study `gStudy`, over its rows `rows`;
# NULL when none of its weighted variants is in the study.
# @noRd
.vcScorePair <- function(
    rows,
    gStudy,
    entry,
    gwasSumStats,
    window,
    alleleFlip
) {
    gdf <- .cipGwasDf(gwasSumStats, rows, window, gStudy)
    aligned <- .vcAlign(entry$weights, entry$variantIds, gdf, alleleFlip)
    if (length(aligned$variantIds) == 0L) {
        return(NULL)
    }
    ref <- .rowsSharedLdRef(
        gwasSumStats,
        rows,
        glue(
            "geneVarianceComponentTest: the rows of GWAS study '{gStudy}' a ",
            "QTL weight set overlaps"
        )
    )
    stats <- .vcStatistics(
        aligned$coef,
        aligned$z,
        .gcLdSource(ref, aligned$variantIds, "geneVarianceComponentTest")
    )
    .vcRecord(
        entry,
        gStudy,
        stats,
        gdf$chrom[[aligned$idxB[[1L]]]],
        range(gdf$pos[aligned$idxB])
    )
}

# Weights on the GWAS variants they match, those with a positive weight:
# the weights (unsigned, so an allele swap leaves them), the GWAS z-scores
# and ids there, and their GWAS rows (`idxB`).
# @noRd
.vcAlign <- function(weights, variantIds, gdf, alleleFlip) {
    m <- matchVariants(variantIds, gdf$variant_id, allowFlip = alleleFlip)
    weight <- weights[m$idxA]
    keep <- weight > 0
    list(
        variantIds = gdf$variant_id[m$idxB][keep],
        coef = weight[keep],
        z = gdf$z[m$idxB][keep],
        idxB = m$idxB[keep]
    )
}

# Q = sum(a z^2), the eigenvalues of its null distribution
# (.vcEigenvalues()), its null mean and SD, p-value and the two gene-level
# statistics.
# @noRd
.vcStatistics <- function(a, z, ld) {
    lambda <- .vcEigenvalues(ld, a)
    statistic <- sum(a * z^2)
    expected <- sum(lambda)
    spread <- sqrt(2 * sum(lambda^2))
    pval <- .quadFormUpperTail(statistic, lambda)
    list(
        nVariants = length(a),
        vcStatistic = statistic,
        vcMean = expected,
        vcSd = spread,
        vcStandardized = (statistic - expected) / spread,
        vcPval = pval,
        vcZ = stats::qnorm(pval, lower.tail = FALSE)
    )
}

# The nonzero eigenvalues of A^1/2 R A^1/2, A = diag(a); from a genotype
# factor F, those of the smaller of (F A^1/2)'(F A^1/2) and its transpose.
# Eigenvalues within rounding of zero are dropped.
# @noRd
.vcEigenvalues <- function(ld, a) {
    root <- sqrt(a)
    values <- if (!is.null(ld$ldFactor)) {
        .vcFactorEigenvalues(ld$ldFactor, root)
    } else {
        eigen(
            outer(root, root) * .ldSourceMatrix(ld),
            symmetric = TRUE,
            only.values = TRUE
        )$values
    }
    values[values > max(values) * length(values) * .Machine$double.eps]
}

# @noRd
.vcFactorEigenvalues <- function(ldFactor, root) {
    scaled <- sweep(ldFactor, 2L, root, "*") / sqrt(nrow(ldFactor) - 1)
    gram <- if (nrow(scaled) < ncol(scaled)) {
        tcrossprod(scaled)
    } else {
        crossprod(scaled)
    }
    eigen(gram, symmetric = TRUE, only.values = TRUE)$values
}

# One result record: the weight set's identity against `gStudy`, its
# statistics, and the span `at` (two positions) on `chrom` it sits at.
# @noRd
.vcRecord <- function(entry, gStudy, stats, chrom, at) {
    tuple <- entry$tuple
    c(
        list(
            qtlStudy = tuple$qStudy,
            context = tuple$qContext,
            trait = tuple$qTrait,
            method = tuple$qMethod,
            jointStudies = tuple$qJoint$jointStudies,
            jointContexts = tuple$qJoint$jointContexts,
            jointTraits = tuple$qJoint$jointTraits,
            gwasStudy = gStudy,
            analysis = "varianceComponent",
            variantWeights = entry$variantWeights,
            pipFloor = entry$pipFloor
        ),
        stats,
        list(chrom = chrom, startPos = at[[1L]], endPos = at[[2L]])
    )
}

# The records as the result GRanges; with none, an empty one with the same
# columns, and a message.
# @noRd
.vcResult <- function(records) {
    if (length(records) == 0L) {
        inform(str_c(
            "geneVarianceComponentTest: no weight set has a weighted variant ",
            "in the GWAS; nothing was tested."
        ))
    }
    df <- if (length(records) == 0L) .vcNoRecords() else bind_rows(records)
    gr <- GenomicRanges::GRanges(
        seqnames = df$chrom,
        ranges = IRanges::IRanges(
            start = as.integer(df$startPos),
            end = as.integer(df$endPos)
        )
    )
    columns <- select(
        df,
        all_of(c(
            "qtlStudy",
            "context",
            "trait",
            "method",
            keep(.JOINT_COLUMNS, .cipColumnHasValues, df = df),
            "gwasStudy",
            "analysis",
            "variantWeights",
            "pipFloor",
            .VC_STATISTICS
        ))
    )
    S4Vectors::`mcols<-`(gr, value = S4Vectors::DataFrame(columns))
}

# The statistic columns of a result, in order.
# @noRd
.VC_STATISTICS <- c(
    "nVariants",
    "vcStatistic",
    "vcMean",
    "vcSd",
    "vcStandardized",
    "vcPval",
    "vcZ"
)

# A result's columns with no rows.
# @noRd
.vcNoRecords <- function() {
    tibble(
        qtlStudy = character(),
        context = character(),
        trait = character(),
        method = character(),
        gwasStudy = character(),
        analysis = character(),
        variantWeights = character(),
        pipFloor = double(),
        nVariants = integer(),
        vcStatistic = double(),
        vcMean = double(),
        vcSd = double(),
        vcStandardized = double(),
        vcPval = double(),
        vcZ = double(),
        chrom = character(),
        startPos = integer(),
        endPos = integer()
    )
}

# =============================================================================
# Gene correlation, and whitening by it
# =============================================================================

#' Correlation of genes' TWAS, MR or variance-component statistics
#'
#' @description For the genes of a \code{\link{causalInferencePipeline}}
#'   result, the correlation between their TWAS z-scores, or their MR
#'   estimates, that the LD between their variants induces, against each GWAS
#'   study; or, for a \code{\link{geneVarianceComponentTest}} result, between
#'   their variance-component statistics. It is the error correlation a
#'   regression of gene-level statistics on gene features whitens by
#'   (generalized least squares), as PoPS does with MAGMA's gene
#'   correlations.
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
#'   A gene's variance-component statistic is a quadratic form,
#'   \eqn{Q_g = z' A_g z} with \eqn{A_g} the diagonal of its variant weights,
#'   and two genes' statistics correlate as
#'   \deqn{\rho_{gh} = \frac{a_g' (R \circ R) a_h}{\sqrt{a_g' (R \circ R)
#'   a_g \; a_h' (R \circ R) a_h}}}{rho_gh = a_g' (R o R) a_h /
#'   sqrt(a_g' (R o R) a_g a_h' (R o R) a_h)}
#'   exactly, \eqn{R \circ R} the elementwise square of the LD (MAGMA's gene
#'   correlation, for weights from QTL fits).
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
#'   non-negative coefficients. These are the scales of TWAS and MR
#'   statistics; \code{"oneSided"} is their default. A variance-component
#'   statistic takes one of two others, \code{"standardized"} its default:
#'   \describe{
#'     \item{\code{"standardized"}}{\eqn{(Q - E Q) / SD(Q)}, whose
#'       correlation is \eqn{\rho} above, exactly.}
#'     \item{\code{"probit"}}{\eqn{\Phi^{-1}(1 - p)} of its p-value, MAGMA's
#'       gene z-score, with \eqn{\rho} above as its correlation, as MAGMA
#'       takes it. That approximation is close for a gene whose statistic
#'       spreads over many eigenvalues and overstates the correlation of one
#'       dominated by a single one.}
#'   }
#'
#' @param x A \code{\link{causalInferencePipeline}} result, or for
#'   \code{analysis = "varianceComponent"} a
#'   \code{\link{geneVarianceComponentTest}} result.
#' @param gwasSumStats The \code{\link{GwasSumStats}} \code{x} was computed
#'   from.
#' @param twasWeights,fineMappingResult The weights \code{x} was computed
#'   from: for TWAS rows \code{twasWeights}, else the
#'   \code{fineMappingResult}, as the pipeline chooses; for MR rows the
#'   \code{fineMappingResult}; for variance-component rows the fits
#'   \code{\link{geneVarianceComponentTest}} read.
#' @param analysis Which rows of \code{x}: \code{"twas"} (default), their
#'   TWAS z-scores, \code{"mr"}, their Wald-ratio MR estimates, or
#'   \code{"varianceComponent"}, the statistics of a
#'   \code{\link{geneVarianceComponentTest}} result.
#' @param method,context Optional: the one method and the one context whose
#'   rows are used. After them each gene (\code{trait}) must have one row per
#'   GWAS study; without them \code{x} must have one already.
#' @param scale The statistic the correlation is of: for TWAS and MR rows
#'   \code{"oneSided"}, \code{"signed"}, \code{"absolute"} or
#'   \code{"squared"}; for variance-component rows \code{"standardized"} or
#'   \code{"probit"}. \code{NULL} (default) takes the analysis's default,
#'   \code{"oneSided"} or \code{"standardized"}. See Details.
#' @param mrParam,alleleFlip As passed to \code{causalInferencePipeline} for
#'   \code{x}.
#' @return A list with one element per GWAS study, named by study: a sparse
#'   symmetric \code{Matrix} of the genes' correlations, its rows and columns
#'   named by \code{trait} in genomic order. A gene with no statistic in
#'   \code{x} has no row.
#' @seealso \code{\link{causalInferencePipeline}},
#'   \code{\link{geneVarianceComponentTest}}, \code{\link{twasZ}}
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
    analysis = c("twas", "mr", "varianceComponent"),
    method = NULL,
    context = NULL,
    scale = NULL,
    mrParam = MrParam(),
    alleleFlip = TRUE
) {
    analysis <- arg_match(analysis)
    resolved <- .gcResolveScale(scale, analysis)
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
    map(studies, .gcScaledMatrix, scale = resolved)
}

# The scale of `analysis`'s statistics: the one given, which must be one of
# its own, or its default when NULL.
# @noRd
.gcResolveScale <- function(scale, analysis) {
    quadratic <- analysis == "varianceComponent"
    if (is.null(scale)) {
        return(if (quadratic) "standardized" else "oneSided")
    }
    own <- if (quadratic) {
        c("standardized", "probit")
    } else {
        c("oneSided", "signed", "absolute", "squared")
    }
    if (length(scale) != 1L || !is_in(scale, own)) {
        abort(glue(
            "`scale` for {analysis} statistics must be one of ",
            "{str_flatten(str_c('\"', own, '\"'), ', ')}."
        ))
    }
    scale
}

# Each GWAS study's genes (.gcGene()) and the signed correlations between
# them (.gcStudyPairs()), for the rows of `x` the arguments choose. The TWAS
# and MR weights are known before the rows are read, a variance component's
# from its rows' kind of weight.
# @noRd
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
    .gcAssertInputs(
        x,
        gwasSumStats,
        twasWeights,
        fineMappingResult,
        method,
        context,
        mrParam,
        alleleFlip
    )
    early <- if (analysis != "varianceComponent") {
        .gcWeightSource(analysis, twasWeights, fineMappingResult, NULL)
    }
    rows <- .gcRows(x, analysis, method, context)
    source <- early %||%
        .gcWeightSource(analysis, twasWeights, fineMappingResult, rows)
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

# The inputs' classes and flags, and their genome builds, which must agree.
# @noRd
#' @importFrom checkmate assertClass assertFlag assertString
.gcAssertInputs <- function(
    x,
    gwasSumStats,
    twasWeights,
    fineMappingResult,
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
    invisible(NULL)
}

# The weights `x`'s rows of `analysis` were computed from: for TWAS the
# pipeline's choice, `twasWeights` else `fineMappingResult`; for MR the
# fine-mapping result; for variance components what their `rows`' variant
# weights read (.gcVcWeightSource()).
# @noRd
.gcWeightSource <- function(analysis, twasWeights, fineMappingResult, rows) {
    if (analysis == "varianceComponent") {
        return(.gcVcWeightSource(rows, twasWeights, fineMappingResult))
    }
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

# The fits variance-component `rows` read, by their one kind of variant
# weight (as geneVarianceComponentTest() reads them): squared weights from
# `twasWeights` else `fineMappingResult`, PIPs from `fineMappingResult`.
# @noRd
.gcVcWeightSource <- function(rows, twasWeights, fineMappingResult) {
    kinds <- unique(rows$variantWeights)
    if (length(kinds) > 1L) {
        abort(str_c(
            "geneCorrelation: `x` holds variance-component rows of more than ",
            "one kind of variant weight; subset it to one."
        ))
    }
    source <- if (kinds == "squaredWeight") {
        twasWeights %||% fineMappingResult
    } else {
        fineMappingResult
    }
    if (is.null(source)) {
        needed <- if (kinds == "squaredWeight") {
            "`twasWeights` (or the `fineMappingResult`)"
        } else {
            "`fineMappingResult`"
        }
        abort(glue(
            "geneCorrelation: \"{kinds}\" rows need the {needed} `x` was ",
            "computed from."
        ))
    }
    source
}

# The rows of `x` the correlation is over, in genomic order: those of
# `analysis` (and of `method` and `context`, when given) that hold a
# statistic, one per gene and GWAS study. Each carries the value its
# recomputed statistic must match (`target`: the TWAS z, the MR estimate, or
# the variance component's Q), the statistic a regression takes (`statistic`:
# the z, the estimate's z, or the standardized Q) and its weight set's key
# (.cellKey()).
# @noRd
.gcRows <- function(x, analysis, method, context) {
    mc <- as_tibble(as.data.frame(S4Vectors::mcols(x)))
    .gcRequireColumns(mc, analysis)
    target <- switch(
        analysis,
        twas = mc$twasZ,
        mr = mc$waldRatio,
        varianceComponent = mc$vcStatistic
    )
    statistic <- switch(
        analysis,
        twas = mc$twasZ,
        mr = mc$waldRatio / mc$waldRatioSe,
        varianceComponent = mc$vcStandardized
    )
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

# The columns `x` needs for `analysis`: a causalInferencePipeline() result's
# for TWAS and MR, a geneVarianceComponentTest() result's for variance
# components.
# @noRd
.gcRequireColumns <- function(mc, analysis) {
    quadratic <- analysis == "varianceComponent"
    need <- c(
        "qtlStudy",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "analysis",
        if (quadratic) {
            c("variantWeights", "pipFloor", .VC_STATISTICS)
        } else {
            c("twasZ", "waldRatio", "waldRatioSe", "nIV")
        }
    )
    absent <- setdiff(need, names(mc))
    if (length(absent) > 0L) {
        producer <- if (quadratic) {
            "geneVarianceComponentTest()"
        } else {
            "causalInferencePipeline()"
        }
        abort(glue(
            "geneCorrelation: `x` is not a {producer} result; it lacks ",
            "column(s) {str_flatten(absent, ', ')}."
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
    input <- if (!is.na(idx)) .gcInput(idx, keys, source, analysis, row)
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
        cellKey = row$cellKey,
        !!!.gcRowExtras(row, analysis)
    )
}

# Weight set `idx`'s inputs, as the result `x` came from reads them.
# @noRd
.gcInput <- function(idx, keys, source, analysis, row) {
    switch(
        analysis,
        twas = .cipTwasTuple(idx, keys, source),
        mr = .cipMrTuple(idx, keys, source),
        varianceComponent = .vcEntry(
            idx,
            keys,
            source,
            row$variantWeights,
            if (!is.na(row$pipFloor)) row$pipFloor
        )
    )
}

# What else a gene carries from its row: the instruments of an MR estimate
# (its covariate), or a variance component's probit z, null mean and kind of
# weight.
# @noRd
.gcRowExtras <- function(row, analysis) {
    if (analysis != "varianceComponent") {
        return(list(nIV = row$nIV))
    }
    list(
        vcZ = row$vcZ,
        vcMean = row$vcMean,
        variantWeights = row$variantWeights
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
    coef <- switch(
        analysis,
        twas = .gcTwasCoefficients(input$weightsInfo, gdf, alleleFlip),
        mr = .gcMrCoefficients(input$fmrEntry, gdf, mrParam, alleleFlip),
        varianceComponent = .gcVcCoefficients(input, gdf, alleleFlip)
    )
    if (is.null(coef)) NULL else list_assign(coef, rows = gRows)
}

# A variance component's variant weights on the GWAS variants
# (.vcAlign()); NULL when none is there.
# @noRd
.gcVcCoefficients <- function(entry, gdf, alleleFlip) {
    aligned <- .vcAlign(entry$weights, entry$variantIds, gdf, alleleFlip)
    if (length(aligned$variantIds) == 0L) {
        return(NULL)
    }
    aligned[c("variantIds", "coef", "z")]
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

# One GWAS study's genes `genes`, each one's variance (w'Rw, or a'Ra; for a
# variance component a'(R o R)a, half its Q's), and the signed correlation of
# each pair that shares a row, i < j.
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
        gStudy = gStudy,
        quadratic = analysis == "varianceComponent"
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
# with i <= j, over the LD reference of the rows they read: W'RW for linear
# statistics, W'(R o R)W for `quadratic` ones (half their covariances).
# @noRd
.gcRowCovariances <- function(r, genes, gwasSumStats, gStudy, quadratic) {
    members <- which(map_lgl(genes, .gcInRow, r = r))
    ref <- .cipRowsLdRef(
        gwasSumStats,
        sort(unique(list_c(map(genes[members], "rows")))),
        gStudy
    )
    ids <- unique(list_c(map(genes[members], "variantIds")))
    W <- exec(cbind, !!!map(genes[members], .gcColumn, ids = ids))
    ld <- .gcLdSource(ref, ids)
    cov <- if (quadratic) .ldSquaredCrossForms(ld, W) else .ldCrossForms(ld, W)
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
.gcLdSource <- function(ldReference, ids, label = "geneCorrelation") {
    if (.ldIsGenotypeRef(ldReference)) {
        return(list(ldFactor = .ldFactorFromRef(ldReference, ids, label)))
    }
    .ldPackedFromRef(ldReference, ids, label)
}

# Each gene's statistic from its coefficients -- the TWAS z, w'z /
# sqrt(w'Rw), the MR estimate, a'z, or the variance component's Q,
# sum(a z^2) -- against the one `x` holds.
# @noRd
.gcCheckStatistics <- function(genes, variance, analysis) {
    statistic <- switch(
        analysis,
        twas = map_dbl(genes, .gcLinearValue) / sqrt(variance),
        mr = map_dbl(genes, .gcLinearValue),
        varianceComponent = map_dbl(genes, .gcQuadraticValue)
    )
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

# @noRd
.gcQuadraticValue <- function(gene) {
    sum(gene$coef * gene$z^2)
}

# The correlation of the statistics on `scale` from that of the z-scores, or
# for a variance component's scales, of its Q's (see geneCorrelation()).
# @noRd
.gcScaleCorrelation <- function(rho, scale) {
    r <- pmax(pmin(rho, 1), -1)
    switch(
        scale,
        signed = r,
        squared = r^2,
        absolute = (sqrt(1 - r^2) + r * asin(r) - 1) / (pi / 2 - 1),
        oneSided = .oneSidedZCorrelation(r),
        standardized = r,
        probit = r
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
