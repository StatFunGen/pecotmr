# In-memory S-LDSC fixtures shared by test_SldscData.R, test_sldscWrapper.R
# and test_sldscPostprocessingPipeline.R.

# =============================================================================
# In-memory builders for SldscData / compute / pipeline tests (no file I/O)
# =============================================================================

# Build one readSldscTrait-shaped run list (the in-memory shape the pipeline
# and standardizeSldscTrait consume).
.sldscMkRun <- function(
    cats,
    h2g = 0.3,
    nBlocks = 10L,
    tau = 1e-7,
    enrichment = 2.0,
    enrichmentP = 0.01
) {
    n <- length(cats)
    list(
        categories = cats,
        tau = setNames(rep(tau, n), cats),
        tauSe = setNames(rep(abs(tau) * 0.3, n), cats),
        enrichment = setNames(rep(enrichment, n), cats),
        enrichmentSe = setNames(rep(enrichment * 0.2, n), cats),
        enrichmentP = setNames(rep(enrichmentP, n), cats),
        propH2 = setNames(rep(0.2, n), cats),
        propSnps = setNames(rep(0.1, n), cats),
        h2g = h2g,
        # Per-block tau has to vary: a constant matrix gives a jackknife
        # tauStarSe of exactly 0, which drops every trait from the
        # DerSimonian-Laird meta (nTraits == 0 and all-NA meta rows).
        annotationJackknifeCoefs = matrix(
            rep(tau * seq(0.9, 1.1, length.out = nBlocks), n),
            nBlocks,
            n,
            dimnames = list(NULL, cats)
        ),
        nBlocks = nBlocks
    )
}

# Build a small, valid in-memory SldscData: 2 traits, 2 target annotations
# (annot_A binary, annot_B continuous), single + optional joint runs.
.sldscMkData <- function(
    withJoint = TRUE,
    withFrq = TRUE,
    traitNames = c("traitX", "traitY")
) {
    annot <- data.frame(
        CHR = c(1, 1, 1, 2, 2, 2),
        SNP = paste0("rs", 1:6),
        annot_A = c(1, 0, 1, 0, 1, 0),
        annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0),
        stringsAsFactors = FALSE
    )
    frq <- if (withFrq) {
        data.frame(
            CHR = c(1, 1, 1, 2, 2, 2),
            SNP = paste0("rs", 1:6),
            MAF = rep(0.2, 6),
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }
    mkTrait <- function() {
        tr <- list(
            single = list(
                .sldscMkRun(c("annot_A_0", "baselineLD_0")),
                .sldscMkRun(c("annot_B_0", "baselineLD_0"))
            )
        )
        tr$joint <- if (withJoint) {
            .sldscMkRun(c("annot_A_0", "annot_B_0", "baselineLD_0"))
        } else {
            NULL
        }
        tr
    }
    traits <- setNames(map(traitNames, function(.) mkTrait()), traitNames)
    SldscData(annotData = annot, frqData = frq, traits = traits)
}
