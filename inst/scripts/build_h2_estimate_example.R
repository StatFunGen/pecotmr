# =============================================================================
# Rebuild h2EstimateExample: a stratified estimateH2() result.
# -----------------------------------------------------------------------------
# g-LDSC on the shipped ldScoreExample (inst/scripts/build_h2_examples.R),
# for a GWAS simulated on that reference's LD in which one baseline
# annotation, over the first half of the variants, carries three times the
# per-variant heritability of the rest. A second
# annotation, over every third variant, is a candidate scored against the
# fit. It needs only the package's own data:
#   Rscript inst/scripts/build_h2_estimate_example.R
# =============================================================================

devtools::load_all(".", quiet = TRUE)
data(ldScoreExample, package = "pecotmr")

M <- length(ldScoreExample)
nGwas <- 10000
inA <- seq_len(M) <= M / 2
perSnpVar <- ifelse(inA, 0.225 / sum(inA), 0.075 / sum(!inA))

# z ~ N(0, N R diag(v) R + R), block by block.
set.seed(1)
z <- numeric(M)
for (block in ldMatrixList(ldScoreExample)) {
    idx <- block$snpIdx
    p <- length(idx)
    sigma <- nGwas * (block$R %*% diag(perSnpVar[idx], p) %*% block$R) + block$R
    e <- eigen((sigma + t(sigma)) / 2, symmetric = TRUE)
    z[idx] <- as.vector(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(p)))
}

genome <- unique(unname(GenomeInfoDb::genome(ldScoreExample)))
gr <- as(ldScoreExample, "GRanges")
S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    SNP = names(ldScoreExample),
    A1 = S4Vectors::mcols(ldScoreExample)$A1,
    A2 = S4Vectors::mcols(ldScoreExample)$A2,
    Z = z,
    N = as.integer(nGwas)
)
gwas <- GwasSumStats(
    studyName = "trait1",
    entry = list(gr),
    genome = genome,
    ldReference = readGenotypes(
        system.file("extdata", "toy_ref.bed", package = "pecotmr")
    )
)
annot <- AnnotationMatrix(
    cbind(
        annotA = as.integer(inA),
        annotB = as.integer(seq_len(M) %% 3L == 0L)
    ),
    as(ldScoreExample, "GRanges"),
    data.frame(
        name = c("annotA", "annotB"),
        tier = c("baseline", "candidate"),
        type = "binary"
    ),
    genome = genome
)

h2EstimateExample <- estimateH2(
    gwas,
    ldScoreExample,
    method = "gldsc",
    annotations = annot
)
print(h2EstimateExample)
stopifnot(
    length(h2EstimateExample) == 1L,
    nrow(stratifiedHeritabilityEnrichment(h2EstimateExample)) == 1L,
    is.matrix(annotationJackknifeCoefs(h2EstimateExample)[[1L]]),
    !is.null(scoreStats(h2EstimateExample)[[1L]])
)

usethis::use_data(h2EstimateExample, overwrite = TRUE, compress = "xz")
