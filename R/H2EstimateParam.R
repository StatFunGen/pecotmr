# =============================================================================
# H2EstimateParam: the settings of one heritability model
# -----------------------------------------------------------------------------
# Every heritability estimator pecotmr implements is a setting of one model of
# the marginal statistics: within a window, their covariance is the LD-weighted
# genetic covariance plus a noise term. LAVA, HDL-L, FastHer, GREML, LDER,
# gLDSC, S-LDSC and HDL differ only in which statistic they fit, by which
# criterion, how the noise and the LD reference's own sampling noise enter,
# and the procedures that guard each fit. H2EstimateParam holds those choices
# as typed fields, none of whose values names a method; the documentation
# lists the settings that reproduce each published one.
#
# The nested records hold one part of the procedures each, as MashPriorParam
# nests its component records: H2TruncationParam (which eigen-directions a fit
# keeps), and one record per criterion -- H2LikelihoodParam,
# H2RegressionParam and H2GlsParam. GREML's arguments are qgg's, so they
# travel as a MethodOptions record, GremlOptions(), not as a Param.
# =============================================================================

#' @include AllGenerics.R MethodParam.R MethodOptions.R
NULL

setClass(
    "H2TruncationParam",
    contains = "MethodParam",
    slots = c(
        fraction = "numeric",
        of = "character",
        maxPerSample = "numeric_OR_NULL",
        minimum = "numeric",
        keepTies = "logical"
    ),
    validity = function(object) .h2TruncationProblems(object)
)

setClass(
    "H2LikelihoodParam",
    contains = "MethodParam",
    slots = c(
        h2Bounds = "numeric",
        interceptBounds = "numeric",
        crossInterceptBounds = "numeric",
        jointFit = "logical",
        varianceFloor = "numeric",
        regressionStart = "logical",
        h2Starts = "numeric",
        interceptStarts = "numeric",
        mustBeatNull = "logical",
        halfBoundaryP = "logical",
        nullSpace = "logical",
        profileScale = "logical",
        maxit = "numeric"
    ),
    validity = function(object) .h2LikelihoodProblems(object)
)

setClass(
    "H2RegressionParam",
    contains = "MethodParam",
    slots = c(
        weights = "character",
        fitWeightIntercept = "logical",
        initialWeights = "character",
        ldRedundancy = "logical",
        iterations = "numeric_OR_NULL",
        twoStage = "logical",
        chi2Threshold = "numeric_OR_NULL",
        chi2Filter = "numeric_OR_NULL",
        jackknifeBlocks = "numeric_OR_NULL",
        eigenvalueTaper = "character",
        interceptFloor = "numeric_OR_NULL"
    ),
    validity = function(object) .h2RegressionProblems(object)
)

setClass(
    "H2GlsParam",
    contains = "MethodParam",
    slots = c(
        ridge = "numeric",
        omegaFraction = "numeric",
        plugIn = "character"
    ),
    validity = function(object) .h2GlsProblems(object)
)

# A statistic other than the eigen one keeps no truncation record, and a
# criterion's procedures are one of the three records, GREML's options, or
# none (the moments).
setClassUnion("H2TruncationParam_OR_NULL", c("H2TruncationParam", "NULL"))
setClassUnion(
    "H2Procedures",
    c(
        "H2LikelihoodParam",
        "H2RegressionParam",
        "H2GlsParam",
        "MethodOptions",
        "NULL"
    )
)

#' @rdname H2EstimateParam
#' @aliases H2EstimateParam-class
#' @exportClass H2EstimateParam
setClass(
    "H2EstimateParam",
    contains = "MethodParam",
    slots = c(
        statistic = "character",
        criterion = "character",
        noise = "character",
        crossNoise = "character",
        effects = "character",
        alpha = "numeric",
        referenceNoise = "character",
        truncation = "H2TruncationParam_OR_NULL",
        sampleSize = "character",
        inference = "character",
        procedures = "H2Procedures"
    ),
    validity = function(object) .h2EstimateParamProblems(object)
)

# The levels of each enumerated field, in the order the constructors list
# them: the first is the default.
# @noRd
.H2_LEVELS <- list(
    statistic = c("eigen", "variant", "phenotype"),
    criterion = c("likelihood", "moments", "regression", "gls", "reml"),
    noise = c("free", "tied", "floored", "supplied"),
    crossNoise = c("window", "supplied", "tied", "zero"),
    effects = c("random", "fixed"),
    referenceNoise = c(
        "eigenTerm",
        "overlapFactor",
        "none",
        "unbiasedR2",
        "weightsOnly"
    ),
    sampleSize = c("median", "mean", "harmonic"),
    inference = c("lrt", "beta", "information", "jackknife"),
    of = c("variants", "eigenvalues"),
    weights = c("inverseVariance", "fixed"),
    initialWeights = c("flat", "moments"),
    eigenvalueTaper = c("none", "linear", "logistic"),
    plugIn = "firstPass"
)

#' @title Heritability Model Settings
#' @description The settings of one heritability model, as
#'   \code{\link{estimateLocalH2}}, \code{\link{estimateStratifiedH2}} and
#'   \code{\link{estimateH2}} fit it, and \code{\link{estimateLocalRg}} its
#'   bivariate form. Within
#'   a window, the marginal statistics' covariance is the LD-weighted genetic
#'   covariance plus a noise term; every estimator pecotmr implements is a
#'   setting of that model. The fields say which statistic is fitted, by
#'   which criterion, how the noise and the LD reference's own sampling noise
#'   enter, how the eigen-directions are truncated, which per-variant sample
#'   size is used and how the uncertainty is reported. The procedures that
#'   guard one criterion's fit are a nested record.
#'
#'   No setting is named after a method. These reproduce the published ones:
#'   \describe{
#'     \item{HDL-L (Li et al. 2025)}{\code{H2EstimateParam()}, the defaults:
#'       eigen, likelihood, free noise, the eigenvalue term, the median sample
#'       size, a likelihood-ratio test, truncated at 99\% of the variant
#'       count with no cap, with \code{H2LikelihoodParam()}'s defaults; the
#'       shared noise of a genetic correlation fitted in the window.}
#'     \item{LAVA (Werme et al. 2022)}{eigen, moments, tied noise, the
#'       overlap factor, the mean sample size and a Beta test, truncated at
#'       99\% of the positive eigenvalues, at most \eqn{0.75 n} directions
#'       and at least 2.}
#'     \item{FastHer (Svishcheva 2025)}{eigen, likelihood, tied noise, no
#'       reference term, the harmonic-mean sample size and the information,
#'       truncated at 99.999\% of the positive eigenvalues keeping ties, with
#'       the null space and the profiled scale on and \eqn{h^2 \le 1 -
#'       10^{-9}}.}
#'     \item{GREML}{the phenotype, REML (\code{qgg::greml}, with
#'       \code{\link{GremlOptions}}), tied noise, no reference term and the
#'       information.}
#'     \item{S-LDSC (Finucane et al. 2015)}{the variants, regression, free
#'       noise, the unbiased \eqn{r^2}, the median sample size and the
#'       jackknife, with the moments' first weights and the LD-redundancy
#'       correction; unstratified, two stages with the intercept fitted below
#'       a chi-square of 30 and the weights updated twice; stratified, fixed
#'       weights and the chi-squares above \eqn{\max(0.001 N, 80)}
#'       dropped.}
#'     \item{g-LDSC (Xiong et al. 2024)}{the variants, GLS, free noise, no
#'       reference term, the median sample size and the jackknife.}
#'     \item{LDER (Song et al. 2022)}{eigen, regression, noise floored at an
#'       intercept of 1, the weights-only reference, every direction, the
#'       median sample size and the jackknife: two stages at LDER's
#'       threshold, and weights that start flat, hold the intercept at 1, are
#'       updated twice and take the linear taper for a reference from the
#'       study's own samples, the logistic one otherwise.}
#'     \item{HDL (Ning et al. 2020), and sHDL stratified}{eigen, likelihood,
#'       free noise, the eigenvalue term, every direction, the median sample
#'       size and the jackknife, from the eigen regression's solution and
#'       from \eqn{h^2 = 0.1} with an intercept of 1, within \eqn{[0, 10]},
#'       a variance floor of \eqn{e^{-10}} and no test against the null.
#'       Like every likelihood here it is fitted by Fisher scoring, which
#'       reaches the maximum HDL's own L-BFGS-B stops short of by up to
#'       \eqn{10^{-4}} of the estimate.}
#'   }
#'   Three local settings follow the genome-wide estimators' criteria in one
#'   window: eigen, regression (inverse-variance weights to convergence),
#'   free noise and the information (LDER's model); the variants, GLS (the
#'   noise-weighting plugged in from a first pass in the window), free noise
#'   and the information (gLDSC's); and the variants, regression
#'   (inverse-variance weights from the moments' first ones, with the
#'   LD-redundancy correction, a jackknife over contiguous blocks of the
#'   window), free noise and the unbiased \eqn{r^2} (S-LDSC's, as MESC fits
#'   it per gene).
#' @param statistic What is fitted: \code{"eigen"} (the default), the
#'   marginal correlations rotated onto the window's LD eigenvectors;
#'   \code{"variant"}, each variant's squared statistic against its LD score;
#'   \code{"phenotype"}, the phenotype itself (a \code{QtlDataset} only).
#' @param criterion How it is fitted: \code{"likelihood"} (the default), the
#'   Gaussian likelihood; \code{"moments"}, closed form;
#'   \code{"regression"}, weighted least squares on the statistics' expected
#'   values; \code{"gls"}, generalized least squares with a plugged-in
#'   covariance of the squares; \code{"reml"}, restricted maximum likelihood
#'   on the phenotype.
#' @param noise Each trait's noise term: \code{"free"} (the default), fitted
#'   as an intercept; \code{"tied"}, the residual \eqn{1 - h^2};
#'   \code{"floored"}, fitted and then held at no less than
#'   \code{interceptFloor} (\code{H2RegressionParam()}); \code{"supplied"},
#'   a value given to the estimator.
#' @param crossNoise For a genetic correlation
#'   (\code{\link{estimateLocalRg}}), the noise two conditions' statistics
#'   share through shared samples: \code{"window"} (the default), fitted in
#'   the window (not under the moments); \code{"supplied"}, the correlation
#'   it is given; \code{"tied"}, the conditions' phenotypic
#'   correlation on their shared donors less their genetic covariance, which
#'   needs individual-level data (a \code{QtlDataset}), and none between
#'   conditions that share no donors; or \code{"zero"}, none. On two LD
#'   references the shared noise's LD is the average of the two. A
#'   heritability estimate does not read it.
#' @param effects How the moments read a genetic correlation's effects:
#'   \code{"random"} (the default), the per-SD effects drawn, as the other
#'   criteria model them, the covariance on two references scaled to the kept
#'   components' share of the per-SD genetic covariance, and its test and
#'   interval drawn with the components' genetic variance; or \code{"fixed"},
#'   LAVA's reading (its preset's), the window's genetic components fixed and
#'   only the noise drawn, so that the covariance is that of the components
#'   the LD holds (on two references, \eqn{\beta_1' R_1^{1/2} R_2^{1/2}
#'   \beta_2} over the two references' kept components), its test holding
#'   each condition's components at their observed directions and
#'   heritabilities. On one reference the estimates agree and only the
#'   inference differs; between two cohorts the fixed reading's test ran at
#'   0.069 at a nominal 0.05, the random one's at 0.035. The other criteria
#'   fit the random-effects model whatever it says, and a heritability does
#'   not read it.
#' @param alpha The frequency exponent of the effects' variance: each
#'   variant's per-SD effect has variance proportional to \eqn{[2p(1 -
#'   p)]^{1 + \alpha}}, \eqn{p} its allele frequency in the LD reference
#'   (an \code{LdMatrix}'s \code{AF}, the \code{MAF} an \code{LdEigen} or
#'   \code{LdScore} records, a panel's or a \code{QtlDataset}'s dosages),
#'   scaled to mean 1 over the window, or genome-wide over the whole
#'   reference, so that the heritability keeps its meaning. \code{-1} (the
#'   default) weighs every variant alike, the per-SD model; \code{0} is the
#'   per-allele model, every allele's effect alike; LDAK's \code{-0.25} lies
#'   between. Across two references a genetic correlation then reads per-SD
#'   effects (\code{-1}) or per-allele effects (\code{0}). A reference
#'   without frequencies is refused for any other value, and a variant it
#'   holds monomorphic weighs nothing. The moments' fixed-effects estimates
#'   do not depend on it.
#' @param referenceNoise The LD reference's own sampling noise:
#'   \code{"eigenTerm"} (the default), HDL's \eqn{-h^2 \lambda / N_{ref}} in
#'   each direction's variance; \code{"overlapFactor"}, LAVA's factor for the
#'   reference's distance from the study's LD, read from its size and the
#'   share of the study's samples in it; \code{"none"}; \code{"unbiasedR2"},
#'   LD scores from the unbiased \eqn{r^2}, \eqn{r^2 - (1 - r^2) / (N_{ref} -
#'   2)}; \code{"weightsOnly"}, a reference that changes only the
#'   regression's weights (its eigenvalue taper).
#' @param truncation For the eigen statistic, which eigen-directions are
#'   kept: an \code{H2TruncationParam()} record. Any other statistic keeps
#'   none.
#' @param sampleSize How the per-variant sample sizes are summarized:
#'   \code{"median"} (the default), \code{"mean"} or \code{"harmonic"}.
#' @param inference The uncertainty reported: \code{"lrt"} (the default),
#'   the likelihood-ratio test and, for a heritability, the interval that
#'   inverts it, the intercept refitted along it; \code{"beta"}, the moments'
#'   Beta test of no heritability; \code{"information"}, a standard error
#'   from the information, a regression's from a sandwich whose middle is its
#'   rows' covariance under the model at the fit (across every pair of
#'   variants, which LD couples, for the variant statistic);
#'   \code{"jackknife"}, a standard error from leaving out blocks.
#' @param procedures The criterion's procedures: an
#'   \code{H2LikelihoodParam()}, \code{H2RegressionParam()} or
#'   \code{H2GlsParam()} record, or \code{\link{GremlOptions}} for REML.
#'   \code{NULL} (the default) is that criterion's default record; the
#'   moments take none.
#' @details Settings that do not belong together are refused:
#'   \itemize{
#'     \item the phenotype goes with REML, and REML with the phenotype;
#'     \item the moments need the eigen statistic and tied noise, and cannot
#'       fit the shared noise in the window;
#'     \item the likelihood needs the eigen statistic, and GLS the variants;
#'     \item the overlap factor and the Beta test go with the moments; the
#'       eigenvalue term with the eigen statistic; the unbiased \eqn{r^2}
#'       with the variants; the weights-only reference with the eigen
#'       regression;
#'     \item the likelihood-ratio test goes with the likelihood, and the
#'       jackknife with a regression, GLS or likelihood;
#'     \item floored noise needs a regression with an \code{interceptFloor},
#'       and two stages free or floored noise;
#'     \item the LD-redundancy correction needs the variants, and the
#'       eigenvalue taper the eigen statistic;
#'     \item the null space and the profiled scale need tied noise, and the
#'       null space no reference term;
#'     \item the procedures must be the criterion's own record.
#'   }
#' @return An \code{H2EstimateParam}, a \code{\link{MethodParam}}.
#' @references Werme J, van der Sluis S, Posthuma D, de Leeuw CA (2022). An
#'   integrated framework for local genetic correlation analysis. Nature
#'   Genetics 54, 274-282.
#'
#'   Li Y, Pawitan Y, Shen X (2025). An enhanced framework for local genetic
#'   correlation analysis. Nature Genetics. doi:10.1038/s41588-025-02123-3.
#'
#'   Svishcheva GR (2025). FastHer: a fast and accurate estimator of local
#'   heritability from GWAS summary statistics. bioRxiv.
#'   doi:10.64898/2025.12.24.696369.
#' @seealso \code{\link{H2TruncationParam}} and the other nested records;
#'   \code{\link{estimateLocalH2}}.
#' @examples
#' # HDL-L's settings, the defaults.
#' H2EstimateParam()
#' # LAVA's.
#' H2EstimateParam(
#'     criterion = "moments",
#'     noise = "tied",
#'     crossNoise = "supplied",
#'     effects = "fixed",
#'     referenceNoise = "overlapFactor",
#'     truncation = H2TruncationParam(
#'         of = "eigenvalues",
#'         maxPerSample = 0.75,
#'         minimum = 2
#'     ),
#'     sampleSize = "mean",
#'     inference = "beta"
#' )
#' @export
H2EstimateParam <- function(
    statistic = c("eigen", "variant", "phenotype"),
    criterion = c("likelihood", "moments", "regression", "gls", "reml"),
    noise = c("free", "tied", "floored", "supplied"),
    crossNoise = c("window", "supplied", "tied", "zero"),
    effects = c("random", "fixed"),
    alpha = -1,
    referenceNoise = c(
        "eigenTerm",
        "overlapFactor",
        "none",
        "unbiasedR2",
        "weightsOnly"
    ),
    truncation = H2TruncationParam(),
    sampleSize = c("median", "mean", "harmonic"),
    inference = c("lrt", "beta", "information", "jackknife"),
    procedures = NULL
) {
    stat <- arg_match(statistic)
    crit <- arg_match(criterion)
    new(
        "H2EstimateParam",
        statistic = stat,
        criterion = crit,
        noise = arg_match(noise),
        crossNoise = arg_match(crossNoise),
        effects = arg_match(effects),
        alpha = alpha,
        referenceNoise = arg_match(referenceNoise),
        truncation = .h2ResolveTruncation(
            truncation,
            stat,
            missing(truncation)
        ),
        sampleSize = arg_match(sampleSize),
        inference = arg_match(inference),
        procedures = procedures %||% .h2DefaultProcedures(crit)
    )
}

# The eigen statistic's truncation as given; any other statistic keeps none,
# so its default record is dropped and one given outright is refused.
# @noRd
.h2ResolveTruncation <- function(truncation, statistic, defaulted) {
    if (identical(statistic, "eigen")) {
        return(truncation %||% H2TruncationParam())
    }
    if (defaulted || is.null(truncation)) {
        return(NULL)
    }
    abort(glue(
        "H2EstimateParam: `truncation` truncates eigen-directions, and the ",
        "\"{statistic}\" statistic has none; leave it unset."
    ))
}

# Each criterion's default procedures: the setting the documentation pairs
# with it. The moments take none.
# @noRd
.h2DefaultProcedures <- function(criterion) {
    switch(
        criterion,
        moments = NULL,
        regression = H2RegressionParam(),
        gls = H2GlsParam(),
        likelihood = H2LikelihoodParam(),
        reml = GremlOptions()
    )
}

#' @title Procedures Of A Heritability Model's Fit
#' @description The records nested in an \code{\link{H2EstimateParam}}:
#'   which eigen-directions a fit keeps, and the procedures that guard each
#'   criterion's fit. Each constructor's defaults are those of the published
#'   setting the documentation pairs it with: HDL-L's truncation and
#'   likelihood, the local eigen regression's and the local GLS's.
#' @param fraction The share of the eigenvalues' sum (\code{of =
#'   "eigenvalues"}, over the positive ones) or of the variant count
#'   (\code{of = "variants"}, the trace) the kept directions first reach.
#'   \code{1} keeps every direction, those with eigenvalue 0 included.
#'   Default 0.99.
#' @param of \code{"variants"} (the default) or \code{"eigenvalues"}; see
#'   \code{fraction}.
#' @param maxPerSample At most this share of the sample size, as a whole
#'   number of directions and never fewer than \code{minimum}; \code{NULL}
#'   (the default) for no cap.
#' @param minimum The fewest directions an estimate is made from; fewer give
#'   none. Default 1.
#' @param keepTies Also keep the directions tied with the last one kept.
#'   Default \code{FALSE}.
#' @param h2Bounds The heritability's range. Default \code{c(0, 1)}.
#' @param interceptBounds A fitted intercept's range. Default \code{c(0,
#'   20)}.
#' @param crossInterceptBounds For a genetic correlation, the range of the
#'   intercept of the noise two conditions' statistics share, on the scale
#'   where it is their noise correlation times the root of the product of
#'   their own intercepts. Default \code{c(-20, 20)}.
#' @param jointFit For a genetic correlation, fit every pair of a group of
#'   three or more conditions at once (\code{\link{estimateLocalRg}}): one
#'   likelihood of all of them, whose genetic correlation matrix is built
#'   from partial correlations so that it is one. Default \code{FALSE}, each
#'   pair on its own, as HDL-L fits them.
#' @param varianceFloor The smallest modelled variance of a direction.
#'   Default \code{exp(-18)}.
#' @param regressionStart Also start from the eigen regression's solution:
#'   the squared rotated statistics regressed on their variance's columns,
#'   then weighted at that fit. Default \code{TRUE}.
#' @param h2Starts,interceptStarts Further starting values; every pair is
#'   tried. Defaults \code{c(0, 0.5)} and \code{c(1, 0.5, 1.5, 0)}.
#' @param mustBeatNull Keep only a fit whose likelihood is no lower than the
#'   intercept-only fit's (for a genetic correlation, no covariance's): one
#'   at the null's boundary is the null. Default \code{TRUE}.
#' @param halfBoundaryP Refer the likelihood ratio to the equal mixture of a
#'   chi-square on one degree of freedom and a point mass at 0, as no
#'   heritability is on the boundary of its range. Default \code{TRUE}.
#' @param nullSpace With tied noise and LD from the study's own samples, add
#'   the directions the genotypes do not span. Default \code{FALSE}.
#' @param profileScale With tied noise, profile out a scale on the
#'   phenotypic variance. Default \code{FALSE}.
#' @param maxit The most Fisher-scoring steps a fit takes from one start,
#'   each solving the expected information against the score and halved
#'   until the likelihood rises. Default \code{1000}.
#' @param weights \code{"inverseVariance"} (the default): one over the
#'   variance of each squared statistic at the current fit;
#'   \code{"fixed"}: set once from a first guess and never updated.
#' @param fitWeightIntercept Whether the inverse-variance weights read the
#'   fitted intercept (\code{TRUE}, the default) or hold it at 1, no
#'   inflation, as LDER's do.
#' @param initialWeights The weights of the first fit: \code{"flat"} (the
#'   default), all equal; or \code{"moments"}, the inverse variance at a
#'   method-of-moments heritability and no inflation, as S-LDSC's.
#' @param ldRedundancy Also divide each variant's weight by its LD score
#'   among the regression's variants, so that variants in LD are not counted
#'   several times. Default \code{FALSE}.
#' @param iterations How many times the weights are updated before the last
#'   solve; \code{NULL} (the default) updates them to convergence. Fixed
#'   weights take none.
#' @param twoStage Fit the intercept on the statistics below a threshold
#'   first, then the heritability with it held. Default \code{FALSE}.
#' @param chi2Threshold The two stages' threshold on the chi-square;
#'   \code{NULL} (the default) is LDER's, \eqn{3.5 \sqrt{m / N} + 2.5}, or
#'   the top \eqn{5 / m} when that first fit's inflation reaches 1.05.
#' @param chi2Filter Drop statistics whose chi-square exceeds \eqn{\max(a N,
#'   b)}, given as \code{c(a, b)}; \code{NULL} (the default) drops none.
#' @param jackknifeBlocks The number of contiguous blocks a window's
#'   jackknife leaves out in turn; \code{NULL} (the default) is the LD
#'   reference's own blocks.
#' @param eigenvalueTaper A further weight on each eigen-direction:
#'   \code{"none"} (the default), \code{"linear"}, \eqn{\min(1, \lambda)},
#'   or \code{"logistic"}, \eqn{1 / (1 + e^{-5 (\lambda - 1)})}.
#' @param interceptFloor With floored noise, the least intercept held, on
#'   the scale where 1 is no inflation; \code{NULL} (the default) for none.
#' @param ridge The ridge on the GLS normal equations, as a share of each
#'   equation's diagonal entry. Default \code{1e-8}.
#' @param omegaFraction The share of the plugged-in covariance's eigenvalue
#'   mass kept in its inverse. Default 0.99.
#' @param plugIn Where the covariance's heritability comes from:
#'   \code{"firstPass"}, a method-of-moments fit of the same statistics.
#' @return A record of the class of the same name, a
#'   \code{\link{MethodParam}}.
#' @examples
#' # LAVA's truncation: 99% of the positive eigenvalues' sum, at most 0.75 n
#' # directions and at least 2.
#' H2TruncationParam(of = "eigenvalues", maxPerSample = 0.75, minimum = 2)
#' # FastHer's likelihood.
#' H2LikelihoodParam(
#'     h2Bounds = c(0, 1 - 1e-9),
#'     regressionStart = FALSE,
#'     mustBeatNull = FALSE,
#'     halfBoundaryP = FALSE,
#'     nullSpace = TRUE,
#'     profileScale = TRUE
#' )
#' # S-LDSC's weights in one window, with a jackknife over 10 blocks.
#' H2RegressionParam(
#'     initialWeights = "moments",
#'     ldRedundancy = TRUE,
#'     iterations = 2,
#'     jackknifeBlocks = 10
#' )
#' H2GlsParam()
#' @name H2TruncationParam
#' @export
H2TruncationParam <- function(
    fraction = 0.99,
    of = c("variants", "eigenvalues"),
    maxPerSample = NULL,
    minimum = 1,
    keepTies = FALSE
) {
    new(
        "H2TruncationParam",
        fraction = fraction,
        of = arg_match(of),
        maxPerSample = maxPerSample,
        minimum = minimum,
        keepTies = keepTies
    )
}

#' @rdname H2TruncationParam
#' @export
H2LikelihoodParam <- function(
    h2Bounds = c(0, 1),
    interceptBounds = c(0, 20),
    crossInterceptBounds = c(-20, 20),
    jointFit = FALSE,
    varianceFloor = exp(-18),
    regressionStart = TRUE,
    h2Starts = c(0, 0.5),
    interceptStarts = c(1, 0.5, 1.5, 0),
    mustBeatNull = TRUE,
    halfBoundaryP = TRUE,
    nullSpace = FALSE,
    profileScale = FALSE,
    maxit = 1000
) {
    new(
        "H2LikelihoodParam",
        h2Bounds = h2Bounds,
        interceptBounds = interceptBounds,
        crossInterceptBounds = crossInterceptBounds,
        jointFit = jointFit,
        varianceFloor = varianceFloor,
        regressionStart = regressionStart,
        h2Starts = h2Starts,
        interceptStarts = interceptStarts,
        mustBeatNull = mustBeatNull,
        halfBoundaryP = halfBoundaryP,
        nullSpace = nullSpace,
        profileScale = profileScale,
        maxit = maxit
    )
}

#' @rdname H2TruncationParam
#' @export
H2RegressionParam <- function(
    weights = c("inverseVariance", "fixed"),
    fitWeightIntercept = TRUE,
    initialWeights = c("flat", "moments"),
    ldRedundancy = FALSE,
    iterations = NULL,
    twoStage = FALSE,
    chi2Threshold = NULL,
    chi2Filter = NULL,
    jackknifeBlocks = NULL,
    eigenvalueTaper = c("none", "linear", "logistic"),
    interceptFloor = NULL
) {
    new(
        "H2RegressionParam",
        weights = arg_match(weights),
        fitWeightIntercept = fitWeightIntercept,
        initialWeights = arg_match(initialWeights),
        ldRedundancy = ldRedundancy,
        iterations = iterations,
        twoStage = twoStage,
        chi2Threshold = chi2Threshold,
        chi2Filter = chi2Filter,
        jackknifeBlocks = jackknifeBlocks,
        eigenvalueTaper = arg_match(eigenvalueTaper),
        interceptFloor = interceptFloor
    )
}

#' @rdname H2TruncationParam
#' @export
H2GlsParam <- function(
    ridge = 1e-8,
    omegaFraction = 0.99,
    plugIn = "firstPass"
) {
    new(
        "H2GlsParam",
        ridge = ridge,
        omegaFraction = omegaFraction,
        plugIn = plugIn
    )
}

#' @title Arguments For GREML's REML Fit
#' @description Arguments for \code{qgg::greml()}, which fits REML when an
#'   \code{\link{H2EstimateParam}}'s criterion is \code{"reml"}. pecotmr
#'   supplies the phenotype and the relationship matrix itself.
#' @param ... Any \code{qgg::greml()} argument, under qgg's own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' GremlOptions(maxit = 200)
#' @export
GremlOptions <- function(...) {
    extra <- list(...)
    fromData <- "supplied from the data by the pipeline"
    .configRefuseOwned(
        extra,
        c(
            y = fromData,
            GRM = fromData,
            GRMlist = fromData,
            fm = fromData,
            data = fromData
        ),
        "GremlOptions"
    )
    .newMethodOptions(
        "qgg::greml",
        defaults = list(),
        extra = extra,
        label = "GremlOptions",
        engine = "greml"
    )
}

# ---- validity ----------------------------------------------------------------

# A problem unless `value` is one non-missing number in [lower, upper].
# @noRd
.h2NumberProblem <- function(value, field, lower = -Inf, upper = Inf) {
    if (length(value) == 1L && !is.na(value)) {
        if (value >= lower && value <= upper) {
            return(character(0))
        }
    }
    as.character(glue(
        "{field} must be one number in [{lower}, {upper}]."
    ))
}

# A problem unless `value` is one whole number of at least `lower`; NULL
# passes where `nullable`.
# @noRd
.h2CountProblem <- function(value, field, lower, nullable = FALSE) {
    if (nullable && is.null(value)) {
        return(character(0))
    }
    whole <- length(value) == 1L &&
        !is.na(value) &&
        value == round(value) &&
        value >= lower
    if (whole) {
        return(character(0))
    }
    as.character(glue("{field} must be one whole number of at least {lower}."))
}

# A problem unless `value` is one non-missing logical.
# @noRd
.h2FlagProblem <- function(value, field) {
    if (length(value) == 1L && !is.na(value)) {
        return(character(0))
    }
    as.character(glue("{field} must be TRUE or FALSE."))
}

# A problem unless `value` is one of `field`'s levels.
# @noRd
.h2LevelProblem <- function(value, field) {
    levels <- .H2_LEVELS[[field]]
    if (length(value) == 1L && is_in(value, levels)) {
        return(character(0))
    }
    as.character(glue(
        "{field} must be one of {str_flatten(levels, ', ')}."
    ))
}

# A problem unless `value` is an increasing pair of numbers.
# @noRd
.h2BoundsProblem <- function(value, field) {
    if (length(value) == 2L && !anyNA(value) && value[[1L]] < value[[2L]]) {
        return(character(0))
    }
    as.character(glue("{field} must be two increasing numbers."))
}

# @noRd
.h2TruncationProblems <- function(object) {
    problems <- c(
        if (
            length(object@fraction) != 1L ||
                is.na(object@fraction) ||
                object@fraction <= 0 ||
                object@fraction > 1
        ) {
            "fraction must be one number in (0, 1]."
        },
        .h2LevelProblem(object@of, "of"),
        if (!is.null(object@maxPerSample)) {
            .h2NumberProblem(object@maxPerSample, "maxPerSample", 0, 1)
        },
        .h2CountProblem(object@minimum, "minimum", 1),
        .h2FlagProblem(object@keepTies, "keepTies")
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.h2LikelihoodProblems <- function(object) {
    flags <- c(
        "jointFit",
        "regressionStart",
        "mustBeatNull",
        "halfBoundaryP",
        "nullSpace",
        "profileScale"
    )
    problems <- c(
        .h2BoundsProblem(object@h2Bounds, "h2Bounds"),
        .h2BoundsProblem(object@interceptBounds, "interceptBounds"),
        .h2BoundsProblem(object@crossInterceptBounds, "crossInterceptBounds"),
        .h2NumberProblem(object@varianceFloor, "varianceFloor", 0),
        list_c(map(flags, .h2SlotFlagProblem, object = object)),
        if (anyNA(object@h2Starts) || anyNA(object@interceptStarts)) {
            "h2Starts and interceptStarts must be numbers."
        },
        .h2CountProblem(object@maxit, "maxit", 1)
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.h2SlotFlagProblem <- function(field, object) {
    .h2FlagProblem(slot(object, field), field)
}

# @noRd
.h2RegressionProblems <- function(object) {
    problems <- c(
        .h2LevelProblem(object@weights, "weights"),
        .h2FlagProblem(object@fitWeightIntercept, "fitWeightIntercept"),
        .h2LevelProblem(object@initialWeights, "initialWeights"),
        .h2FlagProblem(object@ldRedundancy, "ldRedundancy"),
        .h2CountProblem(object@iterations, "iterations", 0, nullable = TRUE),
        if (identical(object@weights, "fixed") && !is.null(object@iterations)) {
            "fixed weights are never updated, so iterations must be NULL."
        },
        .h2FlagProblem(object@twoStage, "twoStage"),
        if (!is.null(object@chi2Threshold)) {
            .h2NumberProblem(object@chi2Threshold, "chi2Threshold", 0)
        },
        if (
            !is.null(object@chi2Filter) &&
                (length(object@chi2Filter) != 2L ||
                    anyNA(object@chi2Filter) ||
                    any(object@chi2Filter < 0))
        ) {
            "chi2Filter must be two non-negative numbers, c(a, b)."
        },
        .h2CountProblem(
            object@jackknifeBlocks,
            "jackknifeBlocks",
            2,
            nullable = TRUE
        ),
        .h2LevelProblem(object@eigenvalueTaper, "eigenvalueTaper"),
        if (!is.null(object@interceptFloor)) {
            .h2NumberProblem(object@interceptFloor, "interceptFloor")
        }
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.h2GlsProblems <- function(object) {
    problems <- c(
        .h2NumberProblem(object@ridge, "ridge", 0),
        if (
            length(object@omegaFraction) != 1L ||
                is.na(object@omegaFraction) ||
                object@omegaFraction <= 0 ||
                object@omegaFraction > 1
        ) {
            "omegaFraction must be one number in (0, 1]."
        },
        .h2LevelProblem(object@plugIn, "plugIn")
    )
    if (length(problems) == 0L) TRUE else problems
}

# Every field's level, then the settings that do not belong together.
# @noRd
.h2EstimateParamProblems <- function(object) {
    fields <- c(
        "statistic",
        "criterion",
        "noise",
        "crossNoise",
        "effects",
        "referenceNoise",
        "sampleSize",
        "inference"
    )
    levels <- list_c(map(fields, .h2SlotLevelProblem, object = object))
    if (length(levels) > 0L) {
        return(levels)
    }
    problems <- c(
        .h2AlphaProblem(object@alpha),
        .h2PairingProblems(object),
        .h2ReferencePairingProblems(object),
        .h2ProceduresProblems(object)
    )
    if (length(problems) == 0L) TRUE else problems
}

# A problem unless the frequency exponent is one finite number.
# @noRd
.h2AlphaProblem <- function(alpha) {
    if (length(alpha) == 1L && is.finite(alpha)) {
        return(character(0))
    }
    "alpha must be one finite number."
}

# @noRd
.h2SlotLevelProblem <- function(field, object) {
    .h2LevelProblem(slot(object, field), field)
}

# The statistic, criterion, noise and inference that go together.
# @noRd
.h2PairingProblems <- function(object) {
    stat <- object@statistic
    crit <- object@criterion
    c(
        if (identical(stat, "phenotype") != identical(crit, "reml")) {
            "the phenotype statistic goes with the reml criterion, and only it."
        },
        if (identical(crit, "moments") && !identical(stat, "eigen")) {
            "the moments criterion needs the eigen statistic."
        },
        if (identical(crit, "moments") && !identical(object@noise, "tied")) {
            "the moments criterion needs tied noise."
        },
        if (
            identical(crit, "moments") &&
                identical(object@crossNoise, "window")
        ) {
            str_c(
                "the moments criterion has no eigenvalue structure to fit ",
                "the shared noise in the window: use another crossNoise."
            )
        },
        if (identical(crit, "likelihood") && !identical(stat, "eigen")) {
            "the likelihood criterion needs the eigen statistic."
        },
        if (identical(crit, "gls") && !identical(stat, "variant")) {
            "the gls criterion needs the variant statistic."
        },
        if (identical(crit, "gls") && !identical(object@noise, "free")) {
            "the gls criterion fits a free noise term only."
        },
        .h2InferenceProblems(object)
    )
}

# @noRd
.h2InferenceProblems <- function(object) {
    crit <- object@criterion
    switch(
        object@inference,
        beta = if (!identical(crit, "moments")) {
            "the beta test goes with the moments criterion."
        },
        lrt = if (!identical(crit, "likelihood")) {
            "the likelihood-ratio test goes with the likelihood criterion."
        },
        information = if (identical(crit, "moments")) {
            "the moments criterion reports the beta test, not the information."
        },
        jackknife = if (!is_in(crit, c("regression", "gls", "likelihood"))) {
            "the jackknife goes with a regression, gls or likelihood fit."
        }
    )
}

# The reference-noise settings and the noise floor, with what reads them.
# @noRd
.h2ReferencePairingProblems <- function(object) {
    ref <- object@referenceNoise
    c(
        if (
            identical(ref, "overlapFactor") &&
                !identical(object@criterion, "moments")
        ) {
            "the overlapFactor reference correction goes with the moments."
        },
        if (
            identical(ref, "eigenTerm") &&
                !identical(object@statistic, "eigen")
        ) {
            "the eigenTerm reference correction needs the eigen statistic."
        },
        if (
            identical(ref, "unbiasedR2") &&
                !identical(object@statistic, "variant")
        ) {
            "the unbiasedR2 reference correction needs the variant statistic."
        },
        if (
            identical(ref, "weightsOnly") &&
                !(identical(object@statistic, "eigen") &&
                    identical(object@criterion, "regression"))
        ) {
            "the weightsOnly reference goes with the eigen regression."
        },
        if (
            identical(object@noise, "floored") &&
                (!is(object@procedures, "H2RegressionParam") ||
                    is.null(object@procedures@interceptFloor))
        ) {
            "floored noise needs a regression with an interceptFloor."
        }
    )
}

# The procedures record must be the criterion's own; the null space and the
# profiled scale need tied noise, and the null space no reference term; a
# regression's two stages fit a free or floored noise term first, its
# LD-redundancy correction needs variants and its eigenvalue taper eigen
# directions.
# @noRd
.h2ProceduresProblems <- function(object) {
    want <- c(
        moments = "NULL",
        regression = "H2RegressionParam",
        gls = "H2GlsParam",
        likelihood = "H2LikelihoodParam",
        reml = "MethodOptions"
    )[[object@criterion]]
    proc <- object@procedures
    if (!is(proc, want)) {
        return(as.character(glue(
            "the {object@criterion} criterion takes ",
            "{if (want == 'NULL') 'no procedures' else want} as its ",
            "procedures."
        )))
    }
    c(
        if (
            identical(want, "MethodOptions") &&
                !identical(metadata(proc)$engine, "greml")
        ) {
            "the reml criterion takes GremlOptions() as its procedures."
        },
        if (is(proc, "H2LikelihoodParam")) .h2LikelihoodPairing(object),
        if (is(proc, "H2RegressionParam")) .h2RegressionPairing(object)
    )
}

# @noRd
.h2RegressionPairing <- function(object) {
    proc <- object@procedures
    c(
        if (proc@twoStage && !is_in(object@noise, c("free", "floored"))) {
            "two stages fit a free or floored noise term first."
        },
        if (proc@ldRedundancy && !identical(object@statistic, "variant")) {
            str_c(
                "the LD-redundancy correction weighs variants: it needs the ",
                "variant statistic."
            )
        },
        if (
            !identical(proc@eigenvalueTaper, "none") &&
                !identical(object@statistic, "eigen")
        ) {
            str_c(
                "the eigenvalue taper weighs eigen directions: it needs the ",
                "eigen statistic."
            )
        }
    )
}

# @noRd
.h2LikelihoodPairing <- function(object) {
    proc <- object@procedures
    c(
        if (
            (proc@nullSpace || proc@profileScale) &&
                !identical(object@noise, "tied")
        ) {
            "the null space and the profiled scale need tied noise."
        },
        if (proc@nullSpace && !identical(object@referenceNoise, "none")) {
            str_c(
                "the null space needs LD from the study's own samples, ",
                "which carries no reference noise: use referenceNoise = ",
                "\"none\"."
            )
        }
    )
}

# ---- the settings of the published estimators --------------------------------

# The settings that reproduce each per-trait estimator the pipelines still
# name (`h2` until their arguments take an H2EstimateParam): HDL-L's are the
# defaults, H2EstimateParam()'s.
# @noRd
.h2LocalSettings <- function(method) {
    switch(
        method,
        lava = H2EstimateParam(
            criterion = "moments",
            noise = "tied",
            crossNoise = "supplied",
            effects = "fixed",
            referenceNoise = "overlapFactor",
            truncation = H2TruncationParam(
                of = "eigenvalues",
                maxPerSample = 0.75,
                minimum = 2
            ),
            sampleSize = "mean",
            inference = "beta"
        ),
        hdlL = H2EstimateParam(
            criterion = "likelihood",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "eigenTerm",
            truncation = H2TruncationParam(
                of = "variants",
                maxPerSample = NULL,
                minimum = 1
            ),
            sampleSize = "median",
            inference = "lrt"
        ),
        fastHer = .h2FastHerSettings(),
        greml = H2EstimateParam(
            statistic = "phenotype",
            criterion = "reml",
            noise = "tied",
            crossNoise = "window",
            referenceNoise = "none",
            sampleSize = "mean",
            inference = "information"
        )
    )
}

# FastHer's settings (Svishcheva 2025).
# @noRd
.h2FastHerSettings <- function() {
    H2EstimateParam(
        criterion = "likelihood",
        noise = "tied",
        crossNoise = "window",
        referenceNoise = "none",
        truncation = H2TruncationParam(
            fraction = 0.99999,
            of = "eigenvalues",
            maxPerSample = NULL,
            minimum = 1,
            keepTies = TRUE
        ),
        sampleSize = "harmonic",
        inference = "information",
        procedures = H2LikelihoodParam(
            h2Bounds = c(0, 1 - 1e-9),
            varianceFloor = 0,
            regressionStart = FALSE,
            mustBeatNull = FALSE,
            halfBoundaryP = FALSE,
            nullSpace = TRUE,
            profileScale = TRUE
        )
    )
}

# The settings of a genome-wide estimate by `method`, stratified or not: the
# ones estimateH2() fits by that method, and the published ones the
# documentation lists for estimateStratifiedH2().
# @noRd
.h2GenomeWideSettings <- function(method, stratified, inSample) {
    switch(
        method,
        sldsc = .h2SldscSettings(stratified),
        gldsc = H2EstimateParam(
            statistic = "variant",
            criterion = "gls",
            noise = "free",
            crossNoise = "window",
            referenceNoise = "none",
            sampleSize = "median",
            inference = "jackknife"
        ),
        lder = .h2LderSettings(inSample),
        hdl = .h2HdlSettings()
    )
}

# S-LDSC's, as ldsc runs it: the first weights at its aggregate (moments)
# estimate; unstratified, two steps with the intercept fitted below a
# chi-square of 30, the weights updated twice; stratified, the weights set
# once and the chi-squares above max(0.001 N, 80) dropped.
# @noRd
.h2SldscSettings <- function(stratified) {
    procedures <- if (stratified) {
        H2RegressionParam(
            weights = "fixed",
            initialWeights = "moments",
            ldRedundancy = TRUE,
            chi2Filter = c(0.001, 80)
        )
    } else {
        H2RegressionParam(
            initialWeights = "moments",
            ldRedundancy = TRUE,
            iterations = 2,
            twoStage = TRUE,
            chi2Threshold = 30
        )
    }
    H2EstimateParam(
        statistic = "variant",
        criterion = "regression",
        noise = "free",
        crossNoise = "window",
        referenceNoise = "unbiasedR2",
        sampleSize = "median",
        inference = "jackknife",
        procedures = procedures
    )
}

# LDER's (Song et al. 2022): every direction, two stages with its own
# threshold, the intercept floored at no inflation, and weights that start
# flat, hold the intercept at 1 and are updated twice with the taper the
# reference's kind sets.
# @noRd
.h2LderSettings <- function(inSample) {
    H2EstimateParam(
        criterion = "regression",
        noise = "floored",
        crossNoise = "window",
        referenceNoise = "weightsOnly",
        truncation = H2TruncationParam(
            fraction = 1,
            of = "eigenvalues",
            maxPerSample = NULL,
            minimum = 1
        ),
        sampleSize = "median",
        inference = "jackknife",
        procedures = H2RegressionParam(
            fitWeightIntercept = FALSE,
            initialWeights = "flat",
            iterations = 2,
            twoStage = TRUE,
            eigenvalueTaper = if (inSample) "linear" else "logistic",
            interceptFloor = 1
        )
    )
}

# HDL's (Ning et al. 2020): every direction, from the eigen regression's
# solution and from h2 0.1 with the intercept 1, within [0, 10], stratified
# or not.
# @noRd
.h2HdlSettings <- function() {
    H2EstimateParam(
        criterion = "likelihood",
        noise = "free",
        crossNoise = "window",
        referenceNoise = "eigenTerm",
        truncation = H2TruncationParam(
            fraction = 1,
            of = "variants",
            maxPerSample = NULL,
            minimum = 1
        ),
        sampleSize = "median",
        inference = "jackknife",
        procedures = H2LikelihoodParam(
            interceptBounds = c(0, 10),
            varianceFloor = exp(-10),
            regressionStart = TRUE,
            h2Starts = 0.1,
            interceptStarts = 1,
            mustBeatNull = FALSE,
            halfBoundaryP = FALSE,
            maxit = 1000
        )
    )
}
