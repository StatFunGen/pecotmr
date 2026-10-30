context("variantId")

# ===========================================================================
# parseVariantId
# ===========================================================================

test_that("parseVariantId: canonical colon format with chr prefix", {
    res <- parseVariantId(c("chr1:100:A:G", "chr2:200:T:C"))
    expect_s3_class(res, "data.frame")
    expect_equal(res$chrom, c("chr1", "chr2"))
    expect_equal(res$pos, c(100L, 200L))
    expect_equal(res$A2, c("A", "T"))
    expect_equal(res$A1, c("G", "C"))
    expect_equal(attr(res, "convention")$alleleSep, ":")
})

test_that("parseVariantId: canonical colon format without chr prefix", {
    res <- parseVariantId(c("1:100:A:G", "2:200:T:C"))
    expect_equal(res$chrom, c("chr1", "chr2"))
    expect_equal(res$pos, c(100L, 200L))
    expect_equal(res$A2, c("A", "T"))
    expect_equal(res$A1, c("G", "C"))
    expect_equal(attr(res, "convention")$hasChr, FALSE)
})

test_that("parseVariantId: underscore separator format", {
    res <- parseVariantId(c("chr1_100_A_G", "1_200_T_C"))
    expect_equal(res$chrom, c("chr1", "chr1"))
    expect_equal(res$pos, c(100L, 200L))
    expect_equal(res$A2, c("A", "T"))
    expect_equal(res$A1, c("G", "C"))
})

test_that("parseVariantId: mixed colon-underscore format", {
    res <- parseVariantId("chr1:100_A_G")
    expect_equal(res$chrom, "chr1")
    expect_equal(res$pos, 100L)
    expect_equal(res$A2, "A")
    expect_equal(res$A1, "G")
    expect_equal(attr(res, "convention")$alleleSep, "_")
})

test_that("parseVariantId: strips build suffix", {
    res <- parseVariantId(c("chr1:100:A:G:b38", "chr1:200:T:C_b37"))
    # The suffix is dropped, not carried as a fifth column.
    expect_equal(colnames(res), c("chrom", "pos", "A2", "A1"))
    expect_equal(res$pos, c(100L, 200L))
    expect_equal(res$A2, c("A", "T"))
    expect_equal(res$A1, c("G", "C"))
    expect_equal(attr(res, "convention")$hasBuild, TRUE)
})

test_that("parseVariantId: convention attribute records hasChr", {
    with_chr <- parseVariantId("chr1:100:A:G")
    no_chr <- parseVariantId("1:100:A:G")
    expect_equal(attr(with_chr, "convention")$hasChr, TRUE)
    expect_equal(attr(no_chr, "convention")$hasChr, FALSE)
})

test_that("parseVariantId: a chrom/pos/A2/A1 frame is read in UCSC style", {
    df <- data.frame(
        chrom = c("chr1", "2"),
        pos = c(100, 200),
        A2 = c("A", "C"),
        A1 = c("G", "T"),
        stringsAsFactors = FALSE
    )
    res <- parseVariantId(df)
    # A bare chromosome gains the prefix; a double position becomes integer.
    expect_equal(res$chrom, c("chr1", "chr2"))
    expect_equal(res$pos, c(100L, 200L))
    expect_equal(attr(res, "convention")$hasChr, TRUE)
})

test_that("parseVariantId: data.frame input with >=4 columns assigns positional names", {
    df <- data.frame(
        c1 = "1",
        c2 = "100",
        c3 = "A",
        c4 = "G",
        extra = "x",
        stringsAsFactors = FALSE
    )
    res <- parseVariantId(df)
    # The first four columns are renamed; the rest pass through.
    expect_equal(names(res), c("chrom", "pos", "A2", "A1", "extra"))
    expect_equal(res$chrom, "chr1")
    expect_equal(res$pos, 100L)
    expect_equal(res$A2, "A")
    expect_equal(res$A1, "G")
})

# ===========================================================================
# .repairVariantIds
# ===========================================================================

test_that(".repairVariantIds: rebuilds a tag-allele id from the allele columns", {
    ids <- c("chr21:13988152:INS:T", "chr21:13988153:DEL:T")
    expect_equal(
        pecotmr:::.repairVariantIds(ids, A2 = c("A", "GT"), A1 = c("AT", "G")),
        c("chr21:13988152:A:AT", "chr21:13988153:GT:G")
    )
})

test_that(".repairVariantIds: leaves clean, lower-case and rsID labels alone", {
    ids <- c(
        "chr1:100:A:G", # plain SNP
        "chr1:788757:T:TAATGG", # real insertion, spelled out
        "chr1:200:a:g", # lower case is DNA, not a tag
        "chr1:300_A_G", # underscore convention
        "rs12345", # does not parse to the allele form
        NA
    )
    out <- pecotmr:::.repairVariantIds(
        ids,
        A2 = c("A", "T", "A", "A", "C", "A"),
        A1 = c("G", "TAATGG", "G", "G", "T", "G")
    )
    expect_equal(out, as.character(ids))
})

test_that(".repairVariantIds: preserves the ids' own naming convention", {
    ids <- c("1:100_A_G", "1:200_INS_G")
    expect_equal(
        pecotmr:::.repairVariantIds(ids, A2 = c("A", "A"), A1 = c("G", "AG")),
        c("1:100_A_G", "1:200_A_AG")
    )
})

test_that(".repairVariantIds: keeps the id when the alleles cannot replace it", {
    ids <- c("chr1:100:INS:T", "chr1:200:INS:T", "chr1:300:INS:T")
    expect_equal(
        pecotmr:::.repairVariantIds(
            ids,
            A2 = c(NA, "", "A"),
            A1 = c("AT", "AT", "AT")
        ),
        c("chr1:100:INS:T", "chr1:200:INS:T", "chr1:300:A:AT")
    )
    # A mis-shaped allele vector is a no-op rather than a recycling error.
    expect_equal(pecotmr:::.repairVariantIds(ids, "A", "AT"), ids)
    expect_equal(
        pecotmr:::.repairVariantIds(character(0), character(0), character(0)),
        character(0)
    )
})

test_that(".repairVariantIds: non-DNA alleles on both sides are still kept", {
    # Nothing better can be said about these, but the point of the repair is
    # that a variant is never DROPPED for how its id spells an allele.
    ids <- c("chr1:100:A:*", "chr1:200:N:G", "chr1:300:A:<DEL>")
    out <- pecotmr:::.repairVariantIds(
        ids,
        A2 = c("A", "N", "A"),
        A1 = c("*", "G", "<DEL>")
    )
    expect_equal(out, ids)
})

# ===========================================================================
# classifyVariantType
# ===========================================================================

test_that("classifyVariantType: identifies SNPs, indels and MNPs", {
    ids <- c(
        "chr1:100:A:G",
        "chr1:200:ATG:A",
        "chr1:300:A:ATG",
        "chr1:400:AT:GC"
    )
    result <- classifyVariantType(ids)
    expect_equal(result, c("SNP", "deletion", "insertion", "MNP"))
})

test_that("classifyVariantType: accepts data.frame input with A2/A1 columns", {
    df <- data.frame(
        A2 = c("A", "AT", "ATG"),
        A1 = c("G", "GC", "A"),
        stringsAsFactors = FALSE
    )
    expect_equal(classifyVariantType(df), c("SNP", "MNP", "deletion"))
})

test_that("classifyVariantType: errors when given an unsupported input", {
    expect_error(
        classifyVariantType(list(a = 1)),
        "character vector of variant IDs or a data.frame"
    )
    expect_error(classifyVariantType(data.frame(foo = 1)), "A2 and A1 columns")
})

# ===========================================================================
# formatVariantId
# ===========================================================================

test_that("formatVariantId uses convention parameter automatically", {
    expect_equal(pecotmr:::formatVariantId(1, 100, "A", "G"), "chr1:100:A:G")

    conv_mixed <- list(hasChr = TRUE, alleleSep = "_")
    expect_equal(
        pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_mixed),
        "chr1:100_A_G"
    )

    conv_nochr <- list(hasChr = FALSE, alleleSep = "_")
    expect_equal(
        pecotmr:::formatVariantId(1, 100, "A", "G", convention = conv_nochr),
        "1:100_A_G"
    )

    expect_equal(
        pecotmr:::formatVariantId(
            1,
            100,
            "A",
            "G",
            chrPrefix = FALSE,
            convention = conv_mixed
        ),
        "chr1:100_A_G"
    )
})

test_that("formatVariantId constructs canonical IDs", {
    expect_equal(
        pecotmr:::formatVariantId(
            c(1, 2),
            c(100, 200),
            c("A", "C"),
            c("G", "T")
        ),
        c("chr1:100:A:G", "chr2:200:C:T")
    )
    expect_equal(
        pecotmr:::formatVariantId(1, 100, "A", "G", chrPrefix = FALSE),
        "1:100:A:G"
    )
    expect_equal(
        pecotmr:::formatVariantId("chr1", 100, "A", "G"),
        "chr1:100:A:G"
    )
})

# =============================================================================
# detectVariantConvention
# =============================================================================

test_that("detectVariantConvention detects chr prefix and allele separators", {
    conv <- pecotmr:::detectVariantConvention(c("chr1:100:A:G", "chr2:200:C:T"))
    expect_true(conv$hasChr)
    expect_equal(conv$alleleSep, ":")
    expect_false(conv$hasBuild)

    conv2 <- pecotmr:::detectVariantConvention(c("1_100_A_G", "2_200_C_T"))
    expect_false(conv2$hasChr)
    expect_equal(conv2$alleleSep, "_")

    conv3 <- pecotmr:::detectVariantConvention(c("chr1:100:A:G:b38"))
    expect_true(conv3$hasBuild)

    conv4 <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
    expect_true(conv4$hasChr)
    expect_equal(conv4$alleleSep, "_")

    conv5 <- pecotmr:::detectVariantConvention(c("1:100_A_G"))
    expect_false(conv5$hasChr)
    expect_equal(conv5$alleleSep, "_")
})

test_that("detectVariantConvention returns defaults for all-NA input", {
    result <- detectVariantConvention(c(NA, NA, NA))
    expect_false(result$hasChr)
    expect_equal(result$alleleSep, ":")
    expect_false(result$hasBuild)
    expect_true(is.na(result$example))
})

# =============================================================================
# normalizeVariantId
# =============================================================================

test_that("normalizeVariantId normalizes various formats", {
    expect_equal(normalizeVariantId("1_100_A_G"), "chr1:100:A:G")
    expect_equal(normalizeVariantId("chr1:100:A:G"), "chr1:100:A:G")
    expect_equal(normalizeVariantId("1:100:A:G"), "chr1:100:A:G")
    expect_equal(
        normalizeVariantId("chr1:100:A:G", chrPrefix = FALSE),
        "1:100:A:G"
    )
    expect_equal(normalizeVariantId("chr1:100:A:G:b38"), "chr1:100:A:G")
    expect_equal(normalizeVariantId("chr1:100_A_G"), "chr1:100:A:G")
    conv <- pecotmr:::detectVariantConvention(c("chr1:100_A_G"))
    expect_equal(
        normalizeVariantId("1:200:C:T", convention = conv),
        "chr1:200_C_T"
    )
})

test_that("normalizeVariantId leaves unparseable ids (rsIDs) unchanged", {
    expect_equal(
        normalizeVariantId(c("rs123", "1:100:A:G", "chr2:200_C_T")),
        c("rs123", "chr1:100:A:G", "chr2:200:C:T")
    ) # rsID passes through
    expect_equal(normalizeVariantId("rs1"), "rs1") # all-unparseable input
})

# ===========================================================================
# canonChrom / chromOrder / non-autosomal chromosomes (string chrom contract)
# ===========================================================================

test_that("canonChrom normalizes prefixes, case, and synonyms to strings", {
    expect_equal(
        pecotmr:::canonChrom(c(
            "chr1",
            "CHR2",
            "chrX",
            "x",
            "23",
            "24",
            "M",
            "chrMT",
            "MT"
        )),
        c("1", "2", "X", "X", "X", "Y", "MT", "MT", "MT")
    )
    # never coerces to integer -- non-autosomes survive
    expect_type(pecotmr:::canonChrom("chrX"), "character")
    # A zero-padded number (UK Biobank's BGEN files write 01) is the same
    # chromosome; a lone 0 and non-numeric labels keep their zeros.
    expect_equal(
        pecotmr:::canonChrom(c("01", "chr01", "001", "10", "0", "023")),
        c("1", "1", "1", "10", "0", "X")
    )
})

test_that("parseVariantId keeps non-autosomal chromosomes as strings", {
    res <- parseVariantId(c(
        "chrX:100:A:G",
        "chr23:200:T:C",
        "chrM:5:G:A",
        "chrY:9:C:T"
    ))
    expect_equal(res$chrom, c("chrX", "chrX", "chrM", "chrY"))
    expect_equal(res$pos, c(100L, 200L, 5L, 9L))
})

test_that("formatVariantId emits non-autosomal ids without 'chrNA'", {
    expect_equal(pecotmr:::formatVariantId("X", 100, "A", "G"), "chrX:100:A:G")
    expect_equal(pecotmr:::formatVariantId("23", 100, "A", "G"), "chrX:100:A:G")
    expect_false(grepl("NA", pecotmr:::formatVariantId("chrMT", 5, "G", "A")))
})

test_that("normalizeVariantId round-trips non-autosomal chromosomes", {
    expect_equal(normalizeVariantId("23:100:A:G"), "chrX:100:A:G")
    expect_equal(normalizeVariantId("chrX:100:A:G"), "chrX:100:A:G")
    expect_equal(normalizeVariantId("chrM:5:G:A"), "chrM:5:G:A")
    expect_equal(normalizeVariantId("MT:5:G:A"), "chrM:5:G:A")
})

test_that("chromOrder sorts autosomes then X/Y/MT", {
    ord <- pecotmr:::chromOrder(c("MT", "X", "2", "1", "22", "Y"))
    expect_equal(as.character(sort(ord)), c("1", "2", "22", "X", "Y", "MT"))
})

# ===========================================================================
# matchVariants -- allele-aware matcher (match on chrom/pos/ref/alt)
# ===========================================================================

test_that("matchVariants matches exact pairs and reports sign +1", {
    a <- c("chr1:100:A:G", "chr2:200:C:T")
    b <- c("chr2:200:C:T", "chr1:100:A:G")
    m <- pecotmr:::matchVariants(a, b)
    ord <- order(m$idxA)
    expect_equal(m$idxA[ord], c(1L, 2L))
    expect_equal(m$idxB[ord], c(2L, 1L)) # a[1]->b[2], a[2]->b[1]
    expect_equal(m$sign[ord], c(1, 1))
})

test_that("matchVariants is robust to chr-prefix and separator differences", {
    m <- pecotmr:::matchVariants("chr1:100:A:G", "1_100_A_G")
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 1L)
    expect_equal(m$sign, 1)
})

test_that("matchVariants reports sign -1 for an allele swap", {
    m <- pecotmr:::matchVariants("chr1:100:A:G", "chr1:100:G:A")
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 1L)
    expect_equal(m$sign, -1)
})

test_that("matchVariants judges each strand-ambiguous (A/T) swap on its own", {
    # An A/T vs T/A pair is a swap or a strand flip and nothing says which:
    # it is dropped by default, whatever else the batch holds ...
    m <- pecotmr:::matchVariants("chr1:100:A:T", "chr1:100:T:A")
    expect_length(m$idxA, 0L)
    a <- c("chr1:100:A:C", "chr2:200:A:T")
    b <- c("chr1:100:T:G", "chr2:200:T:A")
    expect_equal(pecotmr:::matchVariants(a, b)$idxA, 1L)
    # ... and kept as a swap (sign -1) when ambiguity checking is disabled.
    m <- pecotmr:::matchVariants(
        "chr1:100:A:T",
        "chr1:100:T:A",
        removeStrandAmbiguous = FALSE
    )
    expect_equal(m$idxA, 1L)
    expect_equal(m$sign, -1)
    expect_equal(
        sort(pecotmr:::matchVariants(a, b, removeStrandAmbiguous = FALSE)$idxA),
        c(1L, 2L)
    )
    # An exact A/T pair is no ambiguity at all.
    expect_equal(
        pecotmr:::matchVariants("chr1:100:A:T", "chr1:100:A:T")$sign,
        1
    )
})

test_that("matchVariants matches within maxDistance only when alleles agree", {
    a <- c("chr1:100:A:G", "chr1:500:C:T")
    b <- c("chr1:103:A:G", "chr1:520:C:T", "chr1:101:A:C")
    # The same position by default.
    expect_length(pecotmr:::matchVariants(a, b)$idxA, 0L)
    # Within 5 bp: the A:G copy, not the A:C variant beside it.
    m <- pecotmr:::matchVariants(a, b, maxDistance = 5L)
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 1L)
    m <- pecotmr:::matchVariants(a, b, maxDistance = 20L)
    expect_equal(m$idxA, c(1L, 2L))
    expect_equal(m$idxB, c(1L, 2L))
    # A swapped copy in the window matches with sign -1; a window never
    # crosses chromosomes.
    expect_equal(
        pecotmr:::matchVariants(
            "chr1:100:A:G",
            "chr1:102:G:A",
            maxDistance = 2L
        )$sign,
        -1
    )
    expect_length(
        pecotmr:::matchVariants(
            "chr1:100:A:G",
            "chr2:100:A:G",
            maxDistance = 10L
        )$idxA,
        0L
    )
})

test_that("matchVariants takes GRanges of variants", {
    ids <- c("chr1:100:A:G", "chr1:200:C:T")
    # Ensembl-style seqnames; the first variant stored with its alleles
    # swapped.
    gr <- GenomicRanges::GRanges(
        c("1", "1"),
        IRanges::IRanges(c(100L, 200L), width = 1L),
        A2 = c("G", "C"),
        A1 = c("A", "T")
    )
    m <- pecotmr:::matchVariants(gr, ids)
    expect_equal(m$idxA, 1:2)
    expect_equal(m$sign, c(-1, 1))
    expect_equal(pecotmr:::matchVariants(gr, ids, allowFlip = FALSE)$idxA, 2L)
    expect_error(
        pecotmr:::matchVariants(GenomicRanges::granges(gr), ids),
        "A1 and A2"
    )
})

test_that("matchVariants matches non-autosomal variants (chrX, chr23 synonym)", {
    m <- pecotmr:::matchVariants(
        c("chrX:100:A:G", "chr23:200:C:T"),
        c("X:100:A:G", "X:200:C:T")
    )
    ord <- order(m$idxA)
    expect_equal(m$idxA[ord], c(1L, 2L))
    expect_equal(m$idxB[ord], c(1L, 2L))
})

test_that("matchVariants falls back to exact string match for rsIDs", {
    m <- pecotmr:::matchVariants(c("rs1", "rs2", "rs3"), c("rs3", "rs1"))
    ord <- order(m$idxA)
    expect_equal(m$idxA[ord], c(1L, 3L))
    expect_equal(m$idxB[ord], c(2L, 1L))
    expect_equal(m$sign[ord], c(1, 1))
})

test_that("matchVariants disambiguates a multi-allelic site by allele", {
    m <- pecotmr:::matchVariants(
        "chr1:100:A:G",
        c("chr1:100:A:T", "chr1:100:A:G")
    )
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 2L) # the A:G ref, not A:T
    expect_equal(m$sign, 1)
})

test_that("matchVariants returns empty on no overlap", {
    m <- pecotmr:::matchVariants("chr1:100:A:G", "chr2:200:C:T")
    expect_length(m$idxA, 0)
    expect_length(m$idxB, 0)
    expect_length(m$sign, 0)
})

test_that("matchVariants accepts data.frame (chrom/pos/A2/A1) input", {
    a <- data.frame(
        chrom = "1",
        pos = 100L,
        A2 = "A",
        A1 = "G",
        stringsAsFactors = FALSE
    )
    m <- pecotmr:::matchVariants(a, "chr1:100:A:G")
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 1L)
    expect_equal(m$sign, 1)
})

test_that("matchVariants allowFlip = FALSE matches exact alleles only (no swap)", {
    # chr-prefix / separator still reconcile ...
    m1 <- pecotmr:::matchVariants(
        "chr1:100:A:G",
        "1_100_A_G",
        allowFlip = FALSE
    )
    expect_equal(m1$idxA, 1L)
    expect_equal(m1$sign, 1)
    # ... but a ref/alt swap is NOT matched when flipping is disabled.
    m2 <- pecotmr:::matchVariants(
        "chr1:100:A:G",
        "chr1:100:G:A",
        allowFlip = FALSE
    )
    expect_length(m2$idxA, 0)
})

# ===========================================================================
# harmonizeAlleles -- direct tests for branches no caller exercises (every
# caller uses flipStrand = FALSE, removeUnmatched = TRUE, a valid colToFlip,
# and matchMinProp = 0, so these paths are otherwise unreached).
# ===========================================================================

.vid_df <- function(chrom, pos, A2, A1, ...) {
    data.frame(
        chrom = as.character(chrom),
        pos = as.integer(pos),
        A2 = A2,
        A1 = A1,
        ...,
        stringsAsFactors = FALSE
    )
}

test_that("harmonizeAlleles warns and returns empty when nothing overlaps", {
    # Asserted, not suppressed: the warning is half of what this test names.
    expect_warning(
        res <- pecotmr:::harmonizeAlleles(
            .vid_df("1", 100, "A", "G"),
            .vid_df("1", 999, "A", "G")
        ),
        "No matching variants found"
    )
    expect_equal(nrow(res$harmonizedData), 0L)
    expect_equal(attr(res, "qcCounts")$considered, 0L)
})

test_that("harmonizeAlleles removeUnmatched = FALSE retains unmatched target variants", {
    tgt <- .vid_df(c("1", "1"), c(100, 200), c("A", "A"), c("G", "G"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        .vid_df("1", 100, "A", "G"),
        removeUnmatched = FALSE,
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 2L)
})

test_that("harmonizeAlleles flipStrand = TRUE runs the strand-flip branch", {
    # target T/C is the unambiguous strand flip of ref A/G.
    res <- pecotmr:::harmonizeAlleles(
        .vid_df("1", 100, "T", "C"),
        .vid_df("1", 100, "A", "G"),
        flipStrand = TRUE,
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 1L)
})

test_that("harmonizeAlleles complements colToComplement on an allele swap", {
    # target G/A is the ref/alt swap of ref A/G, so an effect-allele freq
    # complements to 1 - af while the signed colToFlip column is negated.
    tgt <- .vid_df("1", 100, "G", "A", af = 0.3, z = 2)
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        .vid_df("1", 100, "A", "G"),
        colToFlip = "z",
        colToComplement = "af",
        matchMinProp = 0
    )
    expect_equal(res$harmonizedData$af, 0.7)
    expect_equal(res$harmonizedData$z, -2)
})

test_that("harmonizeAlleles errors when colToFlip / colToComplement are absent", {
    tgt <- .vid_df("1", 100, "A", "G")
    expect_error(
        pecotmr:::harmonizeAlleles(
            tgt,
            tgt,
            colToFlip = "nope",
            matchMinProp = 0
        ),
        "not found in targetData"
    )
    expect_error(
        pecotmr:::harmonizeAlleles(
            tgt,
            tgt,
            colToComplement = "nope",
            matchMinProp = 0
        ),
        "not found in targetData"
    )
})

test_that("harmonizeAlleles errors when too few variants clear matchMinProp", {
    ref <- .vid_df("1", c(100, 200, 300, 400, 500), "A", "G")
    expect_error(
        pecotmr:::harmonizeAlleles(
            .vid_df("1", 100, "A", "G"),
            ref,
            matchMinProp = 0.5
        ),
        "Not enough variants"
    )
})

test_that("harmonizeAlleles errors on duplicate variant IDs when removeDups = FALSE", {
    tgt <- .vid_df(c("1", "1"), c(100, 100), c("A", "A"), c("G", "G"))
    expect_error(
        pecotmr:::harmonizeAlleles(
            tgt,
            .vid_df("1", 100, "A", "G"),
            matchMinProp = 0
        ),
        "Duplicated variant IDs remain"
    )
})

test_that("harmonizeAlleles removeDups = TRUE drops duplicates with a warning", {
    tgt <- .vid_df(c("1", "1"), c(100, 100), c("A", "A"), c("G", "G"))
    expect_warning(
        res <- pecotmr:::harmonizeAlleles(
            tgt,
            .vid_df("1", 100, "A", "G"),
            removeDups = TRUE,
            matchMinProp = 0
        ),
        "Removed 1 duplicate"
    )
    expect_equal(nrow(res$harmonizedData), 1L)
})

test_that("harmonizeAlleles removeIndels = TRUE drops indels", {
    tgt <- .vid_df(c("1", "1"), c(100, 200), c("A", "AT"), c("G", "A"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        tgt,
        removeIndels = TRUE,
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 1L)
})

# ===========================================================================
# harmonizeAlleles: one row per target variant
#
# The join is on (chrom, pos), so a reference carrying both orientations of a
# variant used to return that target TWICE with opposite signs -- a phantom
# sign-flipped copy that every matchVariants-based LD lookup downstream then
# refused to match. The resolution now happens here, where the surviving
# variant set is decided.
# ===========================================================================

test_that("harmonizeAlleles drops a target the reference cannot decide", {
    tgt <- .vid_df("1", 100, "A", "G", Z = 2.2)
    ref <- .vid_df(c("1", "1"), c(100, 100), c("A", "G"), c("G", "A"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 0L)
})

test_that("harmonizeAlleles sees a strand-recorded second orientation too", {
    # A:G and C:T are the same ambiguity as A:G and G:A -- only a literal
    # string key would miss it.
    tgt <- .vid_df("1", 100, "A", "G", Z = 2.2)
    ref <- .vid_df(c("1", "1"), c(100, 100), c("A", "C"), c("G", "T"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 0L)
})

test_that("harmonizeAlleles keeps one row when the matches agree", {
    # A second reference row at the same position with reconcilable alleles
    # and the SAME sign is not ambiguous, so the target survives -- once.
    tgt <- .vid_df("1", 100, "A", "G", Z = 2.2)
    ref <- .vid_df(c("1", "1"), c(100, 100), c("A", "T"), c("G", "C"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 1L)
    expect_equal(res$harmonizedData$Z, 2.2)
})

test_that("harmonizeAlleles leaves a plain multi-allelic position alone", {
    # A:G and A:T at one position are two different variants, not two records
    # of one, so the target matches exactly one of them.
    tgt <- .vid_df("1", 100, "A", "G", Z = 2.2)
    ref <- .vid_df(c("1", "1"), c(100, 100), c("A", "A"), c("G", "T"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 1L)
    expect_equal(res$harmonizedData$Z, 2.2)
})

test_that("harmonizeAlleles counts what it kept, not what it joined", {
    # The counts feed the summaryStatsQc log line; they are taken after the
    # resolution, so a dropped target is not also reported as a correction.
    tgt <- .vid_df("1", 100, "A", "G", Z = 2.2)
    ref <- .vid_df(c("1", "1"), c(100, 100), c("A", "G"), c("G", "A"))
    res <- pecotmr:::harmonizeAlleles(
        tgt,
        ref,
        colToFlip = "Z",
        matchMinProp = 0
    )
    counts <- attr(res, "qcCounts")
    expect_equal(counts$kept, 0L)
    expect_equal(counts$signFlip, 0L)
})

test_that("harmonizeAlleles surfaces sign/strand/dropped in qcCounts", {
    # 4 shared positions: 100 exact, 200 sign-flip, 300 strand-flip (A/G
    # unambiguous), 400 allele mismatch (dropped).
    target <- data.frame(
        chrom = c(1, 1, 1, 1),
        pos = c(100, 200, 300, 400),
        A2 = c("A", "A", "A", "A"),
        A1 = c("G", "G", "G", "G"),
        z = c(1, 2, 3, 4),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = c(1, 1, 1, 1),
        pos = c(100, 200, 300, 400),
        A2 = c("A", "G", "T", "C"),
        A1 = c("G", "A", "C", "A"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(
        target,
        ref,
        colToFlip = "z",
        matchMinProp = 0
    )
    # Default return shape unchanged.
    expect_named(res, c("harmonizedData", "qcSummary"))
    expect_equal(nrow(res$harmonizedData), 3L)
    cnt <- attr(res, "qcCounts")
    expect_false(is.null(cnt))
    expect_equal(cnt$considered, 4L)
    expect_equal(cnt$signFlip, 1L)
    expect_equal(cnt$strandFlip, 1L)
    expect_equal(cnt$kept, 3L)
    expect_equal(cnt$dropped, 1L)
})

test_that("harmonizeAlleles: accepts a bare variant-id vector as target", {
    res <- pecotmr:::harmonizeAlleles(
        c("chr1:100:A:G", "chr1:200:C:T"),
        c("chr1:100:A:G", "chr1:200:C:T"),
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 2L)
})

test_that("harmonizeAlleles: strips a conflicting variant_id from target", {
    target <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "C"),
        A1 = c("G", "T"),
        variant_id = c("old1", "old2"),
        z = c(1, 2),
        stringsAsFactors = FALSE
    )
    ref <- data.frame(
        chrom = c(1, 1),
        pos = c(100, 200),
        A2 = c("A", "C"),
        A1 = c("G", "T"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(
        target,
        ref,
        colToFlip = "z",
        matchMinProp = 0
    )
    expect_equal(nrow(res$harmonizedData), 2L)
    # The input variant_id was stripped; the returned id is rebuilt from QC'd alleles.
    expect_false(any(c("old1", "old2") %in% res$harmonizedData$variant_id))
})

test_that("harmonizeAlleles renames an empty target column to unnamed_N", {
    # 5 columns whose first four are positional (NOT literally named
    # chrom/pos/A2/A1), so harmonizeAlleles routes through parseVariantId --
    # which preserves the extra column -- rather than the select()-based path.
    # The 5th column has an empty name, so the joined matchResult carries an
    # empty-named column and .sanitizeNames() rewrites it to "unnamed_1".
    target <- data.frame(
        CHR = c(1, 1, 1),
        POS = c(100, 200, 300),
        ref = c("A", "C", "G"),
        alt = c("G", "T", "A"),
        extra = c(0.1, 0.2, 0.3),
        stringsAsFactors = FALSE
    )
    names(target)[5] <- ""
    ref <- data.frame(
        chrom = c(1, 1, 1),
        pos = c(100, 200, 300),
        A2 = c("A", "C", "G"),
        A1 = c("G", "T", "A"),
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::harmonizeAlleles(target, ref, matchMinProp = 0)
    expect_equal(nrow(res$harmonizedData), 3L)
    nm <- colnames(res$harmonizedData)
    expect_false(any(nm == "" | is.na(nm))) # the empty name was repaired
    expect_true(any(grepl("^unnamed_", nm))) # ... to an "unnamed_*" label
})

# ===========================================================================
# withChrPrefix
# ===========================================================================

test_that("withChrPrefix adds a lowercase chr prefix case-insensitively", {
    expect_equal(
        pecotmr:::withChrPrefix(c("1", "chr2", "CHR3", "X")),
        c("chr1", "chr2", "chr3", "chrX")
    )
})

test_that("withChrPrefix spells chromosomes as UCSC does", {
    expect_equal(
        pecotmr:::withChrPrefix(
            c("01", "23", "24", "x", "M", "MT", "chrM", "chrMT", NA)
        ),
        c("chr1", "chrX", "chrY", "chrX", "chrM", "chrM", "chrM", "chrM", NA)
    )
    # A contig outside the standard set keeps its spelling behind the prefix.
    expect_equal(
        pecotmr:::withChrPrefix(c("chrUn_gl000192", "GL000192.1", "XY")),
        c("chrUn_gl000192", "chrGL000192.1", "chrXY")
    )
})

test_that("withChrPrefix agrees with GenomeInfoDb's UCSC style", {
    ncbi <- c(as.character(1:22), "X", "Y", "MT")
    expect_identical(
        pecotmr:::withChrPrefix(ncbi),
        unname(GenomeInfoDb::mapSeqlevels(ncbi, "UCSC"))
    )
    gr <- GenomicRanges::GRanges(
        pecotmr:::withChrPrefix(ncbi),
        IRanges::IRanges(1L, width = 1L)
    )
    expect_identical(GenomeInfoDb::seqlevelsStyle(gr), "UCSC")
})

test_that(".ucscRanges renames to UCSC and keeps the rest", {
    # Anything but ranges passes through for its caller to refuse.
    expect_null(pecotmr:::.ucscRanges(NULL))
    expect_identical(pecotmr:::.ucscRanges("chr1"), "chr1")
    ucsc <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5L, width = 1L))
    expect_identical(pecotmr:::.ucscRanges(ucsc), ucsc)
    gr <- GenomicRanges::GRanges(
        c("1", "MT"),
        IRanges::IRanges(c(100L, 5L), width = 1L, names = c("a", "b")),
        strand = c("+", "-"),
        SNP = c("rs1", "rs2")
    )
    GenomeInfoDb::genome(gr) <- "hg19"
    out <- pecotmr:::.ucscRanges(gr)
    expect_identical(
        as.character(GenomicRanges::seqnames(out)),
        c("chr1", "chrM")
    )
    expect_identical(names(out), c("a", "b"))
    expect_identical(as.character(GenomicRanges::strand(out)), c("+", "-"))
    expect_identical(out$SNP, c("rs1", "rs2"))
    expect_identical(unname(GenomeInfoDb::genome(out)), c("hg19", "hg19"))
})

test_that(".ucscRanges merges two spellings of one chromosome", {
    gr <- GenomicRanges::GRanges(
        c("1", "chr1", "2"),
        IRanges::IRanges(c(10L, 20L, 30L), width = 1L),
        strand = c("+", "*", "-"),
        SNP = c("a", "b", "c")
    )
    GenomeInfoDb::genome(gr) <- "hg38"
    out <- pecotmr:::.ucscRanges(gr)
    expect_identical(GenomeInfoDb::seqlevels(out), c("chr1", "chr2"))
    expect_identical(
        as.character(GenomicRanges::seqnames(out)),
        c("chr1", "chr1", "chr2")
    )
    expect_identical(as.character(GenomicRanges::strand(out)), c("+", "*", "-"))
    expect_identical(out$SNP, c("a", "b", "c"))
    expect_identical(unname(GenomeInfoDb::genome(out)), c("hg38", "hg38"))
})

test_that(".ucscEntries renames every element and keeps a NULL one", {
    out <- pecotmr:::.ucscEntries(list(
        a = GenomicRanges::GRanges("22", IRanges::IRanges(1L, width = 1L)),
        b = NULL
    ))
    expect_identical(names(out), c("a", "b"))
    expect_identical(as.character(GenomicRanges::seqnames(out$a)), "chr22")
    expect_null(out$b)
})

test_that(".ucscExperiment renames an experiment's rowRanges only", {
    m <- matrix(1, 1L, 1L, dimnames = list("v", "s"))
    plain <- SummarizedExperiment::SummarizedExperiment(list(x = m))
    expect_identical(pecotmr:::.ucscExperiment(plain), plain)
    rse <- SummarizedExperiment::SummarizedExperiment(
        list(x = m),
        rowRanges = GenomicRanges::GRanges(
            "X",
            IRanges::IRanges(9L, width = 1L)
        )
    )
    out <- pecotmr:::.ucscExperiment(rse)
    expect_identical(GenomeInfoDb::seqlevels(out), "chrX")
    expect_identical(SummarizedExperiment::assay(out), m)
    expect_identical(pecotmr:::.ucscExperiment(out), out)
})

# ===========================================================================
# matchVariants: data.frame input
# ===========================================================================

test_that("matchVariants allowFlip = FALSE accepts data.frame input", {
    a <- .vid_df("1", 100, "A", "G")
    b <- .vid_df("1", 100, "A", "G")
    m <- pecotmr:::matchVariants(a, b, allowFlip = FALSE)
    expect_equal(m$idxA, 1L)
    expect_equal(m$idxB, 1L)
    expect_equal(m$sign, 1)
})

test_that("matchVariants returns empty for unparseable data.frame input", {
    # NA chrom / pos -> unparseable; with a data.frame there is no id string to
    # fall back to, so the matcher returns empty.
    a <- .vid_df(NA_character_, NA_integer_, "A", "G")
    m <- pecotmr:::matchVariants(a, .vid_df("1", 100, "A", "G"))
    expect_length(m$idxA, 0)
})

# ===========================================================================
# matchVariants: a reference carrying both orientations of one variant
# ===========================================================================

test_that("matchVariants drops a flip-pair target, whatever the row order", {
    # Two distinct indels at one position that are each other's allele flip is
    # rare but biologically real. A target then matches one exactly (+1) and
    # the other as a swap (-1), and nothing can say which it is.
    # The bug this replaces: `!duplicated()` kept whichever match came first,
    # so the SIGN -- and with it the effect direction -- was decided by the
    # panel's row order.
    panel <- c("chr1:100:A:AT", "chr1:100:AT:A")
    forward <- matchVariants("chr1:100:A:AT", panel)
    reversed <- matchVariants("chr1:100:A:AT", rev(panel))
    expect_equal(forward$sign, reversed$sign)
    expect_length(forward$idxA, 0L)
})

test_that("a single orientation still matches, flipped or not", {
    expect_equal(
        matchVariants("chr1:100:A:AT", "chr1:100:A:AT")$sign,
        1
    )
    # An unambiguous flip is a match, not an error: only disagreement is.
    expect_equal(
        matchVariants("chr1:100:A:AT", "chr1:100:AT:A")$sign,
        -1
    )
})

test_that("duplicate reference entries that agree are still matched", {
    # Same sign, so there is nothing to be ambiguous about -- dropping here
    # would lose variants for no reason.
    m <- matchVariants("chr1:300:A:G", c("chr1:300:A:G", "chr1:300:A:G"))
    expect_length(m$idxA, 1L)
    expect_equal(m$sign, 1)
})

test_that("only the ambiguous target is dropped from a mixed set", {
    sumstats <- c("chr1:100:A:AT", "chr1:400:A:G", "chr1:500:C:T")
    panel <- c(
        "chr1:100:A:AT",
        "chr1:100:AT:A",
        "chr1:400:A:G",
        "chr1:500:T:C"
    )
    m <- matchVariants(sumstats, panel)
    expect_equal(sumstats[m$idxA], c("chr1:400:A:G", "chr1:500:C:T"))
    expect_equal(m$sign, c(1, -1))
})

test_that("the rule applies to SNP flip pairs, not just indels", {
    # The ambiguity is in the reference carrying both orientations, which is
    # not specific to indels.
    panel <- c("chr1:200:A:G", "chr1:200:G:A")
    expect_length(matchVariants("chr1:200:A:G", panel)$idxA, 0L)
})

# ===========================================================================
# Remaining helpers and argument guards
# ===========================================================================

test_that("a convention without an allele separator keeps the default", {
    f <- pecotmr:::formatVariantId
    expect_equal(f("1", 100L, "A", "G"), "chr1:100:A:G")
    # hasChr honoured, alleleSep absent -> ":" retained.
    expect_equal(
        f("1", 100L, "A", "G", convention = list(hasChr = TRUE)),
        "chr1:100:A:G"
    )
    # Both honoured.
    expect_equal(
        f(
            "1",
            100L,
            "A",
            "G",
            convention = list(hasChr = FALSE, alleleSep = "_")
        ),
        "1:100_A_G"
    )
})

test_that("columns with no names at all are given placeholder names", {
    df <- as.data.frame(matrix(1:4, 2L, 2L))
    colnames(df) <- NULL
    out <- pecotmr:::.sanitizeNames(df)
    # make.unique de-duplicates the repeated placeholder.
    expect_equal(colnames(out), c("unnamed", "unnamed_1"))
})

test_that(".harmonizeQcCounts attributes drops with and without INDEL", {
    mk <- function(withIndel) {
        d <- data.frame(
            keep = c(TRUE, FALSE, TRUE),
            sign_flip = c(FALSE, FALSE, TRUE),
            strand_flip = c(FALSE, TRUE, FALSE),
            strand_unambiguous = c(TRUE, FALSE, TRUE),
            stringsAsFactors = FALSE
        )
        if (withIndel) {
            d$INDEL <- c(FALSE, TRUE, FALSE)
        }
        d
    }
    bare <- pecotmr:::.harmonizeQcCounts(mk(FALSE))
    expect_equal(bare$considered, 3L)
    expect_equal(bare$kept, 2L)
    # With no INDEL column nothing can be attributed to indels, so the same
    # dropped variant is counted as strand-ambiguous instead.
    expect_equal(bare$droppedIndel, 0L)
    expect_equal(bare$droppedAmbiguous, 1L)
    expect_equal(bare$droppedOther, 0L)
    withIndel <- pecotmr:::.harmonizeQcCounts(mk(TRUE))
    expect_equal(withIndel$droppedIndel, 1L)
    expect_equal(withIndel$droppedAmbiguous, 0L)
    expect_equal(withIndel$droppedOther, 0L)
})

test_that(".harmonizeRestoreUnmatched is a no-op when nothing is unmatched", {
    matchResult <- data.frame(
        A1.ref = c("G", "T"),
        A2.ref = c("A", "C"),
        variants_id_qced = c("chr1:100:A:G", "chr1:200:C:T"),
        keep = c(TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    result <- data.frame(
        variants_id_original = c("chr1:100:A:G", "chr1:200:C:T"),
        stringsAsFactors = FALSE
    )
    targetData <- data.frame(
        chrom = c("chr1", "chr1"),
        pos = c(100L, 200L),
        A1 = c("G", "T"),
        A2 = c("A", "C"),
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.harmonizeRestoreUnmatched(result, matchResult, targetData)
    # Every target variant is already in the result, so nothing is appended.
    expect_identical(out$result, result)
    expect_setequal(names(out$qcSummary), c("A1", "A2", "variant_id"))
    # Control: an unmatched target variant IS appended.
    withExtra <- rbind(
        targetData,
        data.frame(
            chrom = "chr1",
            pos = 300L,
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        )
    )
    appended <- pecotmr:::.harmonizeRestoreUnmatched(
        result,
        matchResult,
        withExtra
    )
    expect_equal(nrow(appended$result), 3L)
})

test_that("parseVariantId: rejects input that is neither character nor df", {
    expect_error(parseVariantId(1:3), "One of the following must apply")
})

test_that("normalizeVariantId: argument guards fire", {
    ids <- "chr1:100:A:G"
    expect_error(normalizeVariantId(1:3), "Must be of type 'character'")
    expect_error(
        normalizeVariantId(ids, chrPrefix = NA),
        "chrPrefix.*May not be NA"
    )
    expect_error(
        normalizeVariantId(ids, chrPrefix = c(TRUE, FALSE)),
        "chrPrefix.*Must have length 1"
    )
    expect_error(
        normalizeVariantId(ids, convention = "nope"),
        "convention.*Must be of type 'list'"
    )
})
