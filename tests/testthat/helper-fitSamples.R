# Real PLINK2 genotypes under a dataset whose contexts observe only 30 of the
# panel's 100 samples, and one of whose traits (g2) is missing a value for five
# of those: what a fit's samples are differs from the panel's and between
# traits. Each context draws its own phenotypes: identical ones would make a
# joint fit's residual covariance singular, and mr.mash's result on it depends
# on the platform's rounding.
.fitSamplesDataset <- function(
    panel = readGenotypes(
        plink2Prefix = file.path(test_path("test_data"), "test_variants")
    )
) {
    samples <- colnames(panel)[1:30]
    traits <- c("g1", "g2", "g3")
    set.seed(1)
    brain <- .fitSamplesExpression(traits, samples)
    age <- rnorm(length(samples))
    liver <- .fitSamplesExpression(traits, samples)
    at <- GenomicRanges::start(SummarizedExperiment::rowRanges(panel))[[100L]]
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = brain),
        rowRanges = `names<-`(
            GenomicRanges::GRanges(
                "chr21",
                IRanges::IRanges(rep(at, length(traits)), width = 100L)
            ),
            traits
        ),
        colData = S4Vectors::DataFrame(age = age, row.names = samples)
    )
    seLiver <- se
    SummarizedExperiment::assay(seLiver, "expression") <- liver
    QtlDataset(
        studyName = "s",
        genotypes = panel,
        phenotypes = list(brain = se, liver = seLiver)
    )
}

# One context's phenotypes, g2 missing for the first five samples.
.fitSamplesExpression <- function(traits, samples) {
    y <- matrix(
        rnorm(length(traits) * length(samples)),
        nrow = length(traits),
        dimnames = list(traits, samples)
    )
    y["g2", 1:5] <- NA
    y
}

# The same genotypes, their file naming every variant by an rsID.
.rsIdPanel <- function() {
    stem <- file.path(test_path("test_data"), "test_variants")
    dir <- tempfile("rsid_")
    dir.create(dir)
    out <- file.path(dir, "panel")
    file.copy(str_c(stem, ".pgen"), str_c(out, ".pgen"))
    file.copy(str_c(stem, ".psam"), str_c(out, ".psam"))
    pvar <- utils::read.delim(
        str_c(stem, ".pvar"),
        check.names = FALSE,
        colClasses = "character"
    )
    utils::write.table(
        mutate(pvar, ID = str_c("rs", row_number())),
        str_c(out, ".pvar"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    readGenotypes(plink2Prefix = out)
}

# Two individual-level studies of the same genotypes and phenotypes.
.fitSamplesTwoStudies <- function() {
    first <- .fitSamplesDataset()
    second <- QtlDataset(
        studyName = "s2",
        genotypes = pecotmr:::.qtlGenotypeSe(first),
        phenotypes = pecotmr:::.qtlPhenotypeList(first)
    )
    MultiStudyQtlDataset(qtlDatasets = list(s = first, s2 = second))
}

# The window spanning the fixture's variants.
.fitSamplesWindow <- function() {
    GenomicRanges::GRanges("chr21:17513228-17592874")
}

# Runs `expr` counting the dosage reads and the residualized X blocks it builds;
# `samples` holds what each read asked for (NULL: every sample of the panel).
.countBlocks <- function(expr) {
    reads <- list()
    builds <- 0L
    realRead <- pecotmr:::.handleDosage
    realBuild <- pecotmr:::.qtlResidualizedGenotypes
    testthat::local_mocked_bindings(
        .handleDosage = function(handle, snpIdx, sampleIdx = NULL) {
            reads <<- c(reads, list(sampleIdx))
            realRead(handle, snpIdx, sampleIdx)
        },
        .qtlResidualizedGenotypes = function(...) {
            builds <<- builds + 1L
            realBuild(...)
        },
        .package = "pecotmr"
    )
    result <- suppressMessages(expr)
    list(
        result = result,
        reads = length(reads),
        samples = reads,
        builds = builds
    )
}

# X as a fit on `samples` should have it: on those samples, each column scaled
# to unit variance on them, and the variants that pass the filter on them.
.expectFitOn <- function(X, samples, expected) {
    testthat::expect_setequal(rownames(X), samples)
    testthat::expect_equal(
        unname(apply(X, 2L, stats::sd)),
        rep(1, ncol(X)),
        tolerance = 1e-8
    )
    testthat::expect_equal(colnames(X), colnames(expected))
}
