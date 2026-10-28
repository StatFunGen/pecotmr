# Tests for the LdMatrix class (R/LdMatrix.R).

# Eight chr22 variants from the bundled eQTL region, with their sample LD.
.ldmFixture <- function() {
    data(eqtlRegionExample, envir = environment())
    X <- eqtlRegionExample$X[, 1:8]
    list(X = X, R = cor(X), ids = colnames(X))
}

# The variant table of `ids` (chrom:pos:A2:A1) as a data frame.
.ldmTable <- function(ids, af = NULL) {
    p <- parseVariantId(ids)
    tbl <- data.frame(
        chrom = p$chrom,
        pos = p$pos,
        A1 = p$A1,
        A2 = p$A2,
        SNP = ids
    )
    if (!is.null(af)) {
        tbl$AF <- af
    }
    tbl
}

test_that("LdMatrix parses its variants from chrom:pos:A2:A1 dimnames", {
    f <- .ldmFixture()
    ld <- LdMatrix(f$R, nSamples = nrow(f$X))
    expect_s4_class(ld, "LdMatrix")
    expect_s4_class(ld, "GRanges")
    expect_equal(length(ld), 8L)
    expect_identical(variantIds(ld), f$ids)
    expect_identical(names(ld), f$ids)
    expect_identical(as.character(GenomicRanges::seqnames(ld))[[1L]], "chr22")
    p <- parseVariantId(f$ids)
    expect_identical(S4Vectors::mcols(ld)$A1, p$A1)
    expect_identical(S4Vectors::mcols(ld)$A2, p$A2)
    expect_identical(ldMatrix(ld), f$R)
    expect_identical(nSamples(ld), nrow(f$X))
    expect_null(mixtureWeights(ld))
})

test_that("LdMatrix takes its variants from a table or a GRanges", {
    f <- .ldmFixture()
    af <- seq(0.05, 0.4, length.out = 8L)
    fromTable <- LdMatrix(unname(f$R), variants = .ldmTable(f$ids, af))
    expect_identical(rownames(ldMatrix(fromTable)), f$ids)
    expect_equal(unname(af(fromTable)), af)
    expect_identical(names(af(fromTable)), f$ids)
    expect_true(is.na(nSamples(fromTable)))

    gr <- GenomicRanges::GRanges(
        "22",
        IRanges::IRanges(parseVariantId(f$ids)$pos, width = 1L),
        A1 = parseVariantId(f$ids)$A1,
        A2 = parseVariantId(f$ids)$A2
    )
    fromRanges <- LdMatrix(f$R, variants = gr, nSamples = 100)
    # Unnamed ranges with no SNP column take the matrix's labels.
    expect_identical(variantIds(fromRanges), f$ids)
    expect_identical(
        as.character(GenomicRanges::seqnames(fromRanges))[[1L]],
        "chr22"
    )
})

test_that("as.matrix() and as(x, 'matrix') return the correlation matrix", {
    f <- .ldmFixture()
    ld <- LdMatrix(cor(f$X), nSamples = nrow(f$X))
    expect_identical(as.matrix(ld), ldMatrix(ld))
    expect_identical(as(ld, "matrix"), ldMatrix(ld))
})

test_that("an unlabelled matrix without variants is refused", {
    f <- .ldmFixture()
    expect_error(LdMatrix(unname(f$R)), "`variants` is required")
    bad <- `dimnames<-`(f$R, list(letters[1:8], letters[1:8]))
    expect_error(LdMatrix(bad), "not chrom:pos:A2:A1")
})

test_that("variant ids come from a variant_id column or the ranges' names", {
    f <- .ldmFixture()
    p <- parseVariantId(f$ids)
    tbl <- data.frame(
        chrom = p$chrom,
        pos = p$pos,
        A1 = p$A1,
        A2 = p$A2,
        variant_id = f$ids
    )
    expect_identical(variantIds(LdMatrix(unname(f$R), variants = tbl)), f$ids)
    gr <- GenomicRanges::GRanges(
        "22",
        IRanges::IRanges(p$pos, width = 1L),
        A1 = p$A1,
        A2 = p$A2
    )
    named <- LdMatrix(unname(f$R), variants = `names<-`(gr, f$ids))
    expect_identical(variantIds(named), f$ids)
    expect_identical(rownames(ldMatrix(named)), f$ids)
})

test_that("variants that cannot describe the matrix are refused", {
    f <- .ldmFixture()
    p <- parseVariantId(f$ids)
    expect_error(
        LdMatrix(f$R, variants = "chr22"),
        "must be a GRanges, a data frame or NULL \\(got character\\)"
    )
    expect_error(
        LdMatrix(
            unname(f$R),
            variants = data.frame(chrom = p$chrom, pos = p$pos)
        ),
        "lacks column\\(s\\) A1, A2"
    )
    expect_error(
        LdMatrix(
            unname(f$R),
            variants = .ldmTable(replace(f$ids, 2L, f$ids[[1L]]))
        ),
        "variant ids must be unique"
    )
})

test_that("a matrix labelled differently from its variants is refused", {
    f <- .ldmFixture()
    tbl <- .ldmTable(rev(f$ids))
    expect_error(
        LdMatrix(f$R, variants = tbl),
        "dimnames do not match the variant ids"
    )
    expect_error(
        LdMatrix(f$R[1:3, 1:3], variants = .ldmTable(f$ids)),
        "lists 8 variant"
    )
})

test_that("validity rejects a matrix that is not a correlation matrix", {
    f <- .ldmFixture()
    offDiag <- f$R
    offDiag[1, 2] <- offDiag[2, 1] <- 1.5
    expect_error(LdMatrix(offDiag), "outside \\[-1, 1\\]")
    badDiag <- `diag<-`(f$R, 0.9)
    expect_error(LdMatrix(badDiag), "diagonal is not 1")
    asym <- f$R
    asym[1, 2] <- 0.9
    expect_error(LdMatrix(asym), "not symmetric")
    withNa <- f$R
    withNa[1, 2] <- withNa[2, 1] <- NA
    expect_error(LdMatrix(withNa), "missing")
    expect_error(LdMatrix(f$R, nSamples = 0), "nSamples")
    expect_error(
        LdMatrix(unname(f$R), variants = .ldmTable(f$ids, rep(2, 8))),
        "allele frequencies in \\[0, 1\\]"
    )
})

test_that("validity catches a part changed after construction", {
    # The constructor never builds these, but `@<-` and `mcols<-` assign a
    # part on its own, and validObject() is what catches the result.
    f <- .ldmFixture()
    ld <- LdMatrix(f$R, nSamples = 30)
    noA1 <- ld
    S4Vectors::mcols(noA1)$A1 <- NULL
    expect_error(validObject(noA1), "carry no A1 column")
    small <- ld
    small@correlation <- f$R[1:2, 1:2]
    expect_error(validObject(small), "is 2 x 2 but there are 8 variant")
    relabelled <- ld
    relabelled@correlation <- `dimnames<-`(f$R, list(rev(f$ids), rev(f$ids)))
    expect_error(validObject(relabelled), "names must be the variant ids")
    withNa <- ld
    withNa@correlation[1, 2] <- NA
    expect_error(validObject(withNa), "holds missing values")
    noSize <- ld
    noSize@nSamples <- 0L
    expect_error(validObject(noSize), "one positive count")
    mix <- LdMatrix(list(ld, ld), mixtureWeights = c(0.5, 0.5))
    mix@mixtureWeights <- c(0.5, 0.6)
    expect_error(validObject(mix), "non-negative and sum to 1")
})

test_that("subsetting an LdMatrix subsets its matrix with it", {
    f <- .ldmFixture()
    ld <- LdMatrix(f$R, nSamples = 50)
    sub <- ld[c(4, 2)]
    expect_s4_class(sub, "LdMatrix")
    expect_identical(variantIds(sub), f$ids[c(4, 2)])
    expect_identical(ldMatrix(sub), f$R[c(4, 2), c(4, 2)])
    expect_identical(nSamples(sub), 50L)
    empty <- ld[integer(0)]
    expect_equal(length(empty), 0L)
    expect_equal(dim(ldMatrix(empty)), c(0L, 0L))
    expect_true(validObject(empty))
})

test_that("a mixture is the weighted average of its components", {
    f <- .ldmFixture()
    half <- seq_len(nrow(f$X) / 2)
    rA <- cor(f$X[half, ])
    rB <- cor(f$X[-half, ])
    ldA <- LdMatrix(
        unname(rA),
        variants = .ldmTable(f$ids, rep(0.2, 8)),
        nSamples = length(half)
    )
    ldB <- LdMatrix(
        unname(rB),
        variants = .ldmTable(f$ids, rep(0.4, 8)),
        nSamples = nrow(f$X) - length(half)
    )
    mix <- LdMatrix(list(ldA, ldB), mixtureWeights = c(0.75, 0.25))
    expect_equal(unname(ldMatrix(mix)), unname(0.75 * rA + 0.25 * rB))
    expect_identical(mixtureWeights(mix), c(0.75, 0.25))
    expect_equal(unname(af(mix)), rep(0.75 * 0.2 + 0.25 * 0.4, 8))
    nA <- length(half)
    nB <- nrow(f$X) - nA
    expect_identical(
        nSamples(mix),
        as.integer(round(1 / (0.75^2 / nA + 0.25^2 / nB)))
    )
})

test_that("mixture components are put into the first one's allele coding", {
    f <- .ldmFixture()
    ldA <- LdMatrix(f$R, variants = .ldmTable(f$ids, rep(0.3, 8)))
    # Component B codes variant 2 the other way round: its id, its alleles,
    # its frequency and the sign of every correlation it takes part in flip.
    p <- parseVariantId(f$ids[[2L]])
    swapped <- str_c("chr22:", p$pos, ":", p$A1, ":", p$A2)
    idsB <- replace(f$ids, 2L, swapped)
    sgn <- replace(rep(1, 8), 2L, -1)
    rB <- `dimnames<-`(f$R * outer(sgn, sgn), list(idsB, idsB))
    ldB <- LdMatrix(
        rB,
        variants = .ldmTable(idsB, replace(rep(0.3, 8), 2L, 0.7))
    )
    mix <- LdMatrix(list(ldA, ldB), mixtureWeights = c(0.5, 0.5))
    expect_equal(ldMatrix(mix), f$R)
    expect_equal(unname(af(mix)), rep(0.3, 8))
    expect_identical(variantIds(mix), f$ids)
})

test_that("a mixture keeps only the variants every component covers", {
    f <- .ldmFixture()
    ldA <- LdMatrix(f$R)
    ldB <- LdMatrix(f$R[2:8, 2:8])
    expect_message(
        mix <- LdMatrix(list(ldA, ldB), mixtureWeights = c(0.5, 0.5)),
        "1 variant\\(s\\) of the first component are not in every component"
    )
    expect_identical(variantIds(mix), f$ids[2:8])
    # A component without a size or frequencies leaves the mixture without.
    expect_true(is.na(nSamples(mix)))
    expect_null(S4Vectors::mcols(mix)$AF)
})

test_that("a mixture needs LdMatrix components and valid weights", {
    f <- .ldmFixture()
    ld <- LdMatrix(f$R)
    expect_error(
        LdMatrix(list(ld), mixtureWeights = 1),
        "at least two LdMatrix"
    )
    expect_error(
        LdMatrix(list(ld, f$R), mixtureWeights = c(0.5, 0.5)),
        "at least two LdMatrix"
    )
    expect_error(
        LdMatrix(list(ld, ld), mixtureWeights = c(0.5, 0.6)),
        "summing to 1"
    )
    expect_error(LdMatrix(list(ld, ld)), "summing to 1")
    expect_error(
        LdMatrix(f$R, mixtureWeights = c(0.5, 0.5)),
        "for combining LD matrices"
    )
    disjoint <- LdMatrix(f$R[1:2, 1:2])
    expect_error(
        LdMatrix(
            list(disjoint, LdMatrix(f$R[3:4, 3:4])),
            mixtureWeights = c(0.5, 0.5)
        ),
        "share no variant"
    )
})

test_that("a Matrix-package matrix is made dense", {
    skip_if_not_installed("Matrix")
    f <- .ldmFixture()
    ld <- LdMatrix(Matrix::Matrix(f$R), nSamples = 10)
    expect_true(is.matrix(ldMatrix(ld)))
    expect_equal(ldMatrix(ld), f$R)
})

test_that("show reports the size, a mixture and missing frequencies", {
    f <- .ldmFixture()
    out <- capture.output(show(LdMatrix(f$R, nSamples = 30)))
    expect_true(any(grepl("8 variants, panel size 30", out)))
    expect_true(any(grepl("no allele frequencies", out)))
    unknown <- capture.output(show(LdMatrix(f$R)))
    expect_true(any(grepl("panel size unknown", unknown)))
    mix <- LdMatrix(
        list(LdMatrix(f$R), LdMatrix(f$R)),
        mixtureWeights = c(0.5, 0.5)
    )
    expect_true(any(grepl("mixture of 2 panels", capture.output(show(mix)))))
})
