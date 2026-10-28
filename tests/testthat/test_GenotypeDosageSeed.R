context("GenotypeDosageSeed")

# The seed's reads, against dosages known from the records written for each
# test (helper-vcf.R): every way one VCF record can overlap another, both
# FORMAT/DS shapes, both SNPRelate storage orders, and SeqArray's 2-bit
# planes.

.vcfOverlapDs <- function() {
    `[<-`(.vcfOverlapGt(), 1L, , value = c(0.1, 0.9, 1.8))
}

samples3 <- c("s1", "s2", "s3")

.seedOf <- function(panel) {
    DelayedArray::seed(SummarizedExperiment::assay(panel, "dosage"))
}

# .vcfOverlapRecords() with one DS value per record: its first ALT allele's.
.oneValueDs <- function() {
    recs <- .vcfOverlapRecords()
    recs[samples3] <- lapply(
        recs[samples3],
        sub,
        pattern = ",.*$",
        replacement = ""
    )
    recs
}

test_that("a VCF seed picks each record out of what its range returns", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    panel <- suppressWarnings(readGenotypes(vcf))
    seed <- .seedOf(panel)
    expect_s4_class(seed, "GenotypeDosageSeed")
    expect_equal(pecotmr:::.gdSeedField(seed), "vcfGt")
    expect_equal(dim(seed), c(7L, 3L))
    expect_equal(dimnames(seed)[[2L]], samples3)
    expect_equal(.dose(panel), .vcfOverlapGt())
    # One record at a time: its range also returns every record overlapping
    # it -- the SNP inside the deletion, the other half of the split site --
    # and none of them may stand in for it.
    for (i in seq_len(7L)) {
        expect_equal(
            unname(S4Arrays::extract_array(seed, list(i, NULL))),
            .vcfOverlapGt()[i, , drop = FALSE]
        )
    }
})

test_that("a VCF seed reads FORMAT/DS, one value per ALT allele", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(
        .vcfOverlapRecords(),
        samples3,
        format = c("GT", "DS"),
        ds = "A"
    )
    panel <- suppressWarnings(readGenotypes(vcf))
    expect_equal(pecotmr:::.gdSeedField(.seedOf(panel)), "vcfDs")
    expect_equal(.dose(panel), .vcfOverlapDs())
})

test_that("a one-value FORMAT/DS is its record's first ALT allele's", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(
        .oneValueDs(),
        samples3,
        format = c("GT", "DS"),
        ds = "1"
    )
    panel <- suppressWarnings(readGenotypes(vcf))
    # The two-ALT record's second ALT allele has no value of its own.
    expect_equal(.dose(panel), `[<-`(.vcfOverlapDs(), 7L, , value = NA))
})

test_that("a seed lays a block out as asked, repeats included", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    seed <- .seedOf(suppressWarnings(readGenotypes(vcf)))
    i <- c(5L, 2L, 2L, 7L)
    j <- c(3L, 1L, 3L)
    expect_equal(
        unname(S4Arrays::extract_array(seed, list(i, j))),
        .vcfOverlapGt()[i, j]
    )
    empty <- S4Arrays::extract_array(seed, list(integer(0), NULL))
    expect_equal(dim(empty), c(0L, 3L))
    expect_equal(
        dim(S4Arrays::extract_array(seed, list(NULL, integer(0)))),
        c(7L, 0L)
    )
})

test_that("a missing GT call reads as NA", {
    skip_if_not_installed("VariantAnnotation")
    records <- data.frame(
        CHROM = "21",
        POS = c(1000L, 2000L),
        ID = c("chr1:100:A:G", "chr1:200:A:G"),
        REF = c("A", "C"),
        ALT = c("G", "T"),
        S1 = c("0/1", "0|0"),
        S2 = c("./.", "1|1"),
        S3 = c("1/1", "./."),
        stringsAsFactors = FALSE
    )
    panel <- readGenotypes(.writeTestVcf(records, c("S1", "S2", "S3")))
    expect_equal(.dose(panel), matrix(c(1, NA, 2, 0, 2, NA), 2L, byrow = TRUE))
})

test_that("a VCF that no longer holds the panel's records is refused", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    seed <- .seedOf(suppressWarnings(readGenotypes(vcf)))
    # A record ranked past those starting at its position is not in the
    # file, as after the file was rewritten.
    records <- pecotmr:::.gdSeedRecords(seed)
    S4Vectors::mcols(records)$posRank[[1L]] <- 2L
    stale <- methods::initialize(seed, records = records)
    expect_error(
        S4Arrays::extract_array(stale, list(1L, NULL)),
        "no longer holds the records"
    )
})

test_that("a seed cut to some rows carries only their records", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    seed <- .seedOf(suppressWarnings(readGenotypes(vcf)))
    keep <- c(2L, 6L, 7L)
    cut <- pecotmr:::.gdSeedRows(seed, keep)
    expect_equal(dim(cut), c(3L, 3L))
    # Rows 6 and 7 are the two ALT alleles of one record.
    expect_equal(length(pecotmr:::.gdSeedRecords(cut)), 2L)
    expect_equal(dimnames(cut)[[1L]], dimnames(seed)[[1L]][keep])
    expect_equal(
        unname(S4Arrays::extract_array(cut, list(NULL, NULL))),
        .vcfOverlapGt()[keep, ]
    )
})

test_that("a saved VCF or GDS panel reads the same once loaded again", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    gds <- test_path("test_data", "test_variants.gds")
    for (panel in list(
        suppressWarnings(readGenotypes(vcf)),
        readGenotypes(gds)
    )) {
        f <- tempfile(fileext = ".rds")
        saveRDS(panel, f)
        expect_equal(.dose(readRDS(f)), .dose(panel))
    }
})

test_that("a SNPRelate GDS seed reads either storage order, 3 as missing", {
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("GDSArray")
    codes <- matrix(c(0L, 1L, 2L, 3L, 1L, 0L), nrow = 2L) # SNPs x samples
    for (snpFirst in c(TRUE, FALSE)) {
        f <- tempfile(fileext = ".gds")
        SNPRelate::snpgdsCreateGeno(
            f,
            genmat = if (snpFirst) codes else t(codes),
            sample.id = c("a", "b", "c"),
            snp.id = c("v1", "v2"),
            snp.chromosome = c(1L, 1L),
            snp.position = c(10L, 20L),
            snp.allele = c("A/G", "C/T"),
            snpfirstdim = snpFirst
        )
        panel <- readGenotypes(f)
        expect_equal(
            pecotmr:::.gdSeedField(.seedOf(panel)),
            if (snpFirst) "snpGdsSnpOrder" else "snpGdsSampleOrder"
        )
        expect_equal(
            .dose(panel),
            matrix(c(0, 2, 1, 1, NA, 0), 2L, byrow = TRUE)
        )
        # The codes count each SNP's first allele, which is therefore A1.
        expect_equal(SummarizedExperiment::rowRanges(panel)$A1, c("A", "C"))
    }
})

test_that("a SeqArray GDS seed decodes its 2-bit planes", {
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    panel <- suppressWarnings(readGenotypes(.writeTestSeqArray(vcf)))
    expect_equal(pecotmr:::.gdSeedField(.seedOf(panel)), "seqArrayGt")
    expect_equal(.dose(panel), .vcfOverlapGt())
})

test_that("a record with too many alleles for one plane is read across two", {
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    # Four alleles need a second 2-bit plane, and a missing allele sets every
    # bit of both. SeqArray stores a haploid call in a diploid file with a
    # missing second allele, so it is missing too, as SeqArray reads it.
    records <- data.frame(
        CHROM = "22",
        POS = c(100L, 500L, 600L),
        ID = c("rs1", "rs9", "."),
        REF = "A",
        ALT = c("G", "C,G,T", "C"),
        s1 = c("0/.", "./.", "0/1"),
        s2 = c("1", "0/.", "1/1"),
        s3 = c("./.", "1|3", "0/0"),
        stringsAsFactors = FALSE
    )
    gds <- .writeTestSeqArray(.writeTestVcf(records, samples3))
    panel <- suppressWarnings(readGenotypes(gds))
    expect_equal(
        .dose(panel),
        matrix(
            c(NA, NA, NA, NA, NA, 1, NA, NA, 0, NA, NA, 1, 1, 2, 0),
            ncol = 3L,
            byrow = TRUE
        )
    )
})

test_that("a SeqArray GDS seed reads FORMAT/DS per ALT allele", {
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(
        .vcfOverlapRecords(),
        samples3,
        format = c("GT", "DS"),
        ds = "A"
    )
    panel <- suppressWarnings(readGenotypes(.writeTestSeqArray(vcf, ds = TRUE)))
    expect_equal(pecotmr:::.gdSeedField(.seedOf(panel)), "seqArrayDs")
    # SeqArray stores DS as 32-bit floats.
    expect_equal(.dose(panel), .vcfOverlapDs(), tolerance = 1e-6)
})

test_that("a SeqArray DS row with no value of its own is missing", {
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(
        .oneValueDs(),
        samples3,
        format = c("GT", "DS"),
        ds = "1"
    )
    panel <- suppressWarnings(readGenotypes(.writeTestSeqArray(vcf, ds = TRUE)))
    want <- `[<-`(.vcfOverlapDs(), 7L, , value = NA)
    expect_equal(.dose(panel), want, tolerance = 1e-6)
    # A block of only such rows is answered without reading the file.
    seed <- .seedOf(panel)
    expect_true(all(is.na(S4Arrays::extract_array(seed, list(7L, NULL)))))
})

test_that("a GT call counts one ALT allele, and is missing if any allele is", {
    # Only "./." used to read as missing, so ".|." read as 2 and "0|." and
    # "." as 1; and any non-reference allele counted, so "0/2" read as one
    # copy of the first ALT allele.
    gt <- c("./.", ".|.", "0|.", ".", NA, "0/0", "0/1", "1|1", "1", "0/2")
    count <- function(allele) {
        pecotmr:::vcfGtDosageCpp(
            matrix(gt, ncol = 1L),
            rep(allele, length(gt))
        )[, 1L]
    }
    expect_equal(count(1L), c(NA, NA, NA, NA, NA, 0, 1, 2, 1, 0))
    expect_equal(count(2L), c(NA, NA, NA, NA, NA, 0, 0, 0, 0, 1))
    # Any ploidy, multi-digit allele indices, and anything that is not a GT
    # call read as missing.
    odd <- matrix(c("1/1/0", "10|12", "", "1/a", "1/"), ncol = 1L)
    expect_equal(
        pecotmr:::vcfGtDosageCpp(odd, c(1L, 12L, 1L, 1L, 1L))[, 1L],
        c(2, 1, NA, NA, NA)
    )
    # Each row counts its own allele, across every sample.
    block <- matrix(c("0/1", "1/2", "2/2", "0/0"), 2L)
    expect_equal(
        pecotmr:::vcfGtDosageCpp(block, c(1L, 2L)),
        matrix(c(1, 1, 0, 0), 2L)
    )
})


# ---- PLINK1 ----------------------------------------------------------------

test_that("a PLINK1 seed reads any request through PlinkMatrix", {
    # Each distinct row and sample is read once, in file order, then laid
    # out as asked, repeats and all.
    skip_if_not_installed("PlinkMatrix")
    stem <- test_path("test_data", "test_variants")
    seed <- pecotmr:::.panelSeed(readGenotypes(plink1Prefix = stem))
    expect_identical(pecotmr:::.gdSeedField(seed), "plinkBed")
    all <- S4Arrays::extract_array(seed, list(NULL, NULL))
    got <- S4Arrays::extract_array(seed, list(c(5L, 2L, 5L), c(9L, 1L)))
    expect_equal(got, all[c(5L, 2L, 5L), c(9L, 1L)])
    expect_equal(
        dim(S4Arrays::extract_array(seed, list(integer(0), NULL))),
        c(0L, ncol(all))
    )
})
