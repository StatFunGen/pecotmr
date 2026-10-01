# =============================================================================
# enloc scoring
# -----------------------------------------------------------------------------
# The colocalisation probabilities fastenloc computes, for one eQTL signal
# cluster, under the DAP-1 approximation.
#
# fastenloc enumerates a p x p grid over (GWAS hit m, causal eQTL n) and sums
# it (controller.cc, "case 4"). That is quadratic in the cluster size, but the
# per-cell prior takes only two values -- one on the diagonal, one off it --
# so the off-diagonal total factorises and the whole grid collapses to a
# handful of vector operations:
#
#   diagSum = priD * sum(bf * pip)
#   offSum  = priO * (sum(bf) * sum(pip) - sum(bf * pip))
#
# The two agree exactly; this file computes the second.
#
# The four configurations fastenloc enumerates are coloc's hypotheses, so both
# parameterisations come out of one pass: RCP is PP.H4 (both traits causal at
# the same variant) and LCP is PP.H3 + PP.H4 (both causal, same variant or
# not). What makes this enloc rather than coloc is the prior, which depends on
# whether the GWAS hit IS the eQTL -- pi1_e = expit(a0 + a1) on the diagonal
# against pi1_ne = expit(a0) off it. coloc applies one p12 to every shared
# configuration and so cannot express that dependence.
# =============================================================================

#' @include AllGenerics.R
NULL

# The GWAS Bayes factors fastenloc recovers from fine-mapped PIPs, contrasting
# against "every variant in the cluster has no effect" (controller.cc):
#   log10bf = log10(pip) - log10(P_gwas) + log10(1 - P_gwas) - log10(1 - cpip)
# A zero PIP is floored, as upstream does, so the log is finite.
# @noRd
.enlocGwasLog10Bf <- function(gwasPip, pGwas, clusterPip) {
    floored <- replace(gwasPip, gwasPip == 0, 1e-10)
    log10(floored) - log10(pGwas) + (log10(1 - pGwas) - log10(1 - clusterPip))
}

# The two priors a configuration can carry. `pi1e` applies when the GWAS hit
# is itself the causal eQTL, `pi1ne` when it is not; the enrichment is
# precisely the gap between them.
# @noRd
.enlocConfigPriors <- function(pi1e, pi1ne, p) {
    list(
        diagonal = pi1e * (1 - pi1ne)^(p - 1),
        offDiagonal = pi1ne * (1 - pi1e) * (1 - pi1ne)^(p - 2)
    )
}

#' @title enloc Colocalisation Probabilities For One Signal Cluster
#' @description The quantities fastenloc reports for a single eQTL signal
#'   cluster, given the GWAS Bayes factors and the eQTL PIPs over the same
#'   variants.
#'
#'   Both parameterisations are returned because they come from one
#'   computation: \code{rcp} is \code{ppH4} and \code{lcp} is
#'   \code{ppH3 + ppH4}.
#' @param gwasLog10Bf Numeric vector of per-variant GWAS log10 Bayes factors.
#' @param qtlPip Numeric vector of per-variant eQTL PIPs, same length and
#'   order.
#' @param pi1e Prior probability of GWAS causality for a variant that is a
#'   causal eQTL, \code{expit(a0 + a1)}.
#' @param pi1ne The same for a variant that is not, \code{expit(a0)}.
#' @param p12 Prior probability of colocalisation, used only for the
#'   small-probability floor fastenloc applies to each SNP-level term.
#' @param applyFloor Logical. Apply that floor. fastenloc skips it on the
#'   summary-statistics path, where not all variants are present.
#' @return A list with \code{ppH0} to \code{ppH4}, \code{rcp}, \code{lcp} and
#'   \code{scp} (the per-variant colocalisation probabilities whose sum is
#'   \code{rcp}).
#' @examples
#' enlocClusterProbabilities(
#'   gwasLog10Bf = c(2, 0.1, 0.1),
#'   qtlPip = c(0.7, 0.2, 0.05),
#'   pi1e = 0.02,
#'   pi1ne = 0.001
#' )
#' @export
enlocClusterProbabilities <- function(
    gwasLog10Bf,
    qtlPip,
    pi1e,
    pi1ne,
    p12 = 0,
    applyFloor = FALSE
) {
    p <- length(gwasLog10Bf)
    if (p != length(qtlPip)) {
        abort(glue(
            "enlocClusterProbabilities: `gwasLog10Bf` and `qtlPip` must ",
            "describe the same variants ({p} vs {length(qtlPip)})."
        ))
    }
    bf <- 10^gwasLog10Bf
    # The cluster's total eQTL posterior, capped as upstream does so that
    # log(1 - clusterPip) stays finite.
    clusterPip <- min(sum(qtlPip), 1 - 1e-8)
    pri <- .enlocConfigPriors(pi1e, pi1ne, p)
    # case 4, split by whether the GWAS hit is the causal eQTL
    perVariant <- pri$diagonal * bf * qtlPip
    if (applyFloor) {
        perVariant <- pmax(perVariant, p12)
    }
    diagSum <- sum(perVariant)
    offSum <- pri$offDiagonal *
        (sum(bf) * sum(qtlPip) - sum(bf * qtlPip))
    # cases 1-3: neither, GWAS only, eQTL only
    nullBoth <- (1 - pi1ne)^p * (1 - clusterPip)
    gwasOnly <- pi1ne * (1 - pi1ne)^(p - 1) * sum(bf) * (1 - clusterPip)
    qtlOnly <- (1 - pi1e) * (1 - pi1ne)^(p - 1) * sum(qtlPip)
    nc <- nullBoth + gwasOnly + qtlOnly + diagSum + offSum
    list(
        ppH0 = nullBoth / nc,
        ppH1 = gwasOnly / nc,
        ppH2 = qtlOnly / nc,
        ppH3 = offSum / nc,
        ppH4 = diagSum / nc,
        rcp = diagSum / nc,
        lcp = (offSum + diagSum) / nc,
        scp = perVariant / nc
    )
}

# enloc's priors from the enrichment parameters, as
# fastenloc's set_enrich_params(a0, a1) computes them. `a1` is capped first,
# matching the --cap_a1 option upstream.
# @noRd
.enlocPriorsFromEnrichment <- function(a0, a1, pQtl, capA1 = Inf) {
    a1 <- min(a1, capA1)
    pi1e <- stats::plogis(a0 + a1)
    pi1ne <- stats::plogis(a0)
    list(
        pi1e = pi1e,
        pi1ne = pi1ne,
        p1 = (1 - pQtl) * pi1ne,
        p2 = pQtl * (1 - pi1e),
        p12 = pQtl * pi1e
    )
}

# The inverse: the enrichment parameters implied by a p1/p2/p12 triple, as
# fastenloc's set_enrich_params(p1, p2, p12) computes them. Used to accept an
# enrichment table that carries the priors rather than a0/a1.
# @noRd
.enlocEnrichmentFromPriors <- function(p1, p2, p12) {
    rest <- 1 - p1 - p2 - p12
    list(
        a0 = log(p1 / rest),
        a1 = log(p12 * rest / (p1 * p2))
    )
}

# Per-variant posteriors within one signal cluster. A SuSiE effect's log Bayes
# factors become posterior inclusion probabilities by softmax under the
# uniform prior over variants -- which is what SuSiE's own `alpha` is -- so a
# single effect row is exactly enloc's notion of a signal cluster.
# @noRd
.enlocClusterPip <- function(log10Bf) {
    shifted <- log10Bf - max(log10Bf)
    w <- 10^shifted
    w / sum(w)
}

# Score one (QTL effect, GWAS effect) pair the way fastenloc scores a signal
# cluster. Everything the scorer needs is implied by the priors: the two
# marginals are P_gwas = p1 + p12 and P_eqtl = p2 + p12, and (a0, a1) invert
# out of the triple, so an enrichment table carrying p1/p2/p12 is sufficient.
# @noRd
.enlocScoreEffectPair <- function(qtlLog10Bf, gwasLog10Bf, priors) {
    pGwas <- priors$p1 + priors$p12
    pQtl <- priors$p2 + priors$p12
    en <- .enlocEnrichmentFromPriors(priors$p1, priors$p2, priors$p12)
    pr <- .enlocPriorsFromEnrichment(en$a0, en$a1, pQtl)
    qtlPip <- .enlocClusterPip(qtlLog10Bf)
    gwasPip <- .enlocClusterPip(gwasLog10Bf)
    # A SuSiE effect's alpha is conditional on that effect being present, so
    # it sums to exactly 1 -- the cluster is certain to hold a causal variant.
    # fastenloc's own cap on locus_epip is what keeps log(1 - cpip) finite
    # there, and it has to be applied before the Bayes factors are recovered,
    # not only inside the scorer.
    clusterPip <- min(sum(qtlPip), 1 - 1e-8)
    enlocClusterProbabilities(
        gwasLog10Bf = .enlocGwasLog10Bf(gwasPip, pGwas, clusterPip),
        qtlPip = qtlPip,
        pi1e = pr$pi1e,
        pi1ne = pr$pi1ne,
        p12 = priors$p12,
        applyFloor = FALSE
    )
}
