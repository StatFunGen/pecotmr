test_that(".annotDetectFormat() detects bigwig, bed, ldsc_annot", {
    detect <- pecotmr:::.annotDetectFormat

    expect_equal(detect("path/to/file.bw"), "bigwig")
    expect_equal(detect("path/to/file.bigwig"), "bigwig")
    expect_equal(detect("path/to/file.annot.gz"), "ldsc_annot")
    expect_equal(detect("path/to/file.annot"), "ldsc_annot")
    expect_equal(detect("path/to/file.bed"), "bed")
    # Unknown extensions default to bed
    expect_equal(detect("path/to/file.txt"), "bed")
})

# =============================================================================
# .readBedAnnotation (h2Annotations.R)
# =============================================================================

test_that(".readBedAnnotation returns correct binary overlap vector", {
    skip_if_not_installed("rtracklayer")

    bed_content <- paste(
        "chr1\t100\t500\tregion1",
        "chr1\t800\t900\tregion2",
        sep = "\n"
    )
    bed_file <- tempfile(fileext = ".bed")
    on.exit(unlink(bed_file), add = TRUE)
    writeLines(bed_content, bed_file)

    # SNPs at positions 50, 200, 300, 600, 850
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 5),
        ranges = IRanges::IRanges(start = c(50, 200, 300, 600, 850), width = 1)
    )

    result <- pecotmr:::.readBedAnnotation(bed_file, snp_gr)

    expect_equal(length(result), 5)
    # Position 50 is outside [100,500) and [800,900) -> 0
    expect_equal(result[1], 0)
    # Positions 200 and 300 are inside [100,500) -> 1
    expect_equal(result[2], 1)
    expect_equal(result[3], 1)
    # Position 600 is outside both regions -> 0
    expect_equal(result[4], 0)
    # Position 850 is inside [800,900) -> 1
    expect_equal(result[5], 1)
})


test_that(".readBedAnnotation returns all zeros when no overlaps", {
    skip_if_not_installed("rtracklayer")

    bed_content <- "chr1\t100\t200\tregion1"
    bed_file <- tempfile(fileext = ".bed")
    on.exit(unlink(bed_file), add = TRUE)
    writeLines(bed_content, bed_file)

    # All SNPs outside the region
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3),
        ranges = IRanges::IRanges(start = c(50, 300, 500), width = 1)
    )

    result <- pecotmr:::.readBedAnnotation(bed_file, snp_gr)
    expect_equal(result, c(0, 0, 0))
})

# =============================================================================
# .readLdscAnnot (h2Annotations.R)
# =============================================================================

test_that(".readLdscAnnot reads annotation column and matches SNPs", {
    annot_df <- data.frame(
        CHR = c(1, 1, 1, 1),
        BP = c(100, 200, 300, 400),
        SNP = paste0("rs", 1:4),
        CM = 0,
        my_annot = c(1, 0, 1, 0),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # SNP positions matching BP 100, 200, 300, 400
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(100, 200, 300, 400), width = 1)
    )

    result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "my_annot")

    expect_equal(length(result), 4)
    expect_equal(result, c(1, 0, 1, 0))
})


test_that(".readLdscAnnot returns 0 for unmatched SNP positions", {
    annot_df <- data.frame(
        CHR = c(1, 1),
        BP = c(100, 200),
        SNP = c("rs1", "rs2"),
        CM = 0,
        score = c(5, 10),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # SNP at position 100 matches, position 999 does not
    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(100, 999), width = 1)
    )

    result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "score")
    expect_equal(result, c(5, 0))
})


test_that(".readLdscAnnot errors when annotation column not found", {
    annot_df <- data.frame(
        CHR = c(1),
        BP = c(100),
        SNP = c("rs1"),
        CM = 0,
        my_annot = c(1),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    snp_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100, width = 1)
    )

    expect_error(
        pecotmr:::.readLdscAnnot(annot_file, snp_gr, "nonexistent_col"),
        "Annotation column.*nonexistent_col.*not found"
    )
})


test_that(".readLdscAnnot handles chr-prefixed CHR column", {
    annot_df <- data.frame(
        CHR = c("chr1", "chr1"),
        BP = c(100, 200),
        SNP = c("rs1", "rs2"),
        CM = 0,
        val = c(3, 7),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(100, 200), width = 1)
    )

    result <- pecotmr:::.readLdscAnnot(annot_file, snp_gr, "val")
    expect_equal(result, c(3, 7))
})

# =============================================================================
# readAnnotations method (h2Annotations.R)
# =============================================================================

test_that("readAnnotations errors when path is not named", {
    snp_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100, width = 1)
    )

    expect_error(
        readAnnotations(c("/some/file.bed"), snp_gr),
        "path.*must be a named"
    )
})


test_that("readAnnotations with BED file creates AnnotationMatrix", {
    skip_if_not_installed("rtracklayer")

    bed_content <- "chr1\t100\t500\tregion1"
    bed_file <- tempfile(fileext = ".bed")
    on.exit(unlink(bed_file), add = TRUE)
    writeLines(bed_content, bed_file)

    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3),
        ranges = IRanges::IRanges(start = c(50, 200, 600), width = 1)
    )

    paths <- c(enhancer = bed_file)
    result <- readAnnotations(paths, snp_gr, genome = "hg38")

    expect_s4_class(result, "AnnotationMatrix")
    expect_equal(nrow(assay(result, "annotations")), 3)
    expect_equal(ncol(assay(result, "annotations")), 1)
    expect_equal(colnames(assay(result, "annotations")), "enhancer")
    expect_equal(unname(GenomeInfoDb::genome(result)), "hg38")
    # Auto-detected meta should be binary, candidate tier
    expect_equal(SummarizedExperiment::colData(result)$type, "binary")
    expect_equal(SummarizedExperiment::colData(result)$tier, "candidate")
    # Overlap check: only position 200 is inside [100,500)
    expect_equal(as.numeric(assay(result, "annotations")[, 1]), c(0, 1, 0))
})


test_that("readAnnotations with LDSC annot file creates AnnotationMatrix", {
    annot_df <- data.frame(
        CHR = c(1, 1, 1),
        BP = c(100, 200, 300),
        SNP = paste0("rs", 1:3),
        CM = 0,
        my_score = c(1, 0, 1),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3),
        ranges = IRanges::IRanges(start = c(100, 200, 300), width = 1)
    )

    paths <- c(my_score = annot_file)
    result <- readAnnotations(paths, snp_gr)

    expect_s4_class(result, "AnnotationMatrix")
    expect_equal(as.numeric(assay(result, "annotations")[, 1]), c(1, 0, 1))
    expect_equal(colnames(assay(result, "annotations")), "my_score")
})


test_that("readAnnotations uses the provided annotationMeta", {
    annot_df <- data.frame(
        CHR = c(1, 1),
        BP = c(100, 200),
        SNP = c("rs1", "rs2"),
        CM = 0,
        custom = c(1, 0),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(100, 200), width = 1)
    )

    custom_meta <- data.frame(
        name = "custom",
        tier = "baseline",
        type = "continuous",
        stringsAsFactors = FALSE
    )

    paths <- c(custom = annot_file)
    result <- readAnnotations(paths, snp_gr, annotationMeta = custom_meta)

    expect_s4_class(result, "AnnotationMatrix")
    # Should use the provided meta, not auto-detected
    expect_equal(SummarizedExperiment::colData(result)$tier, "baseline")
    expect_equal(SummarizedExperiment::colData(result)$type, "continuous")
})


test_that("readAnnotations with multiple files creates multi-column matrix", {
    # Two LDSC annot files
    annot_df1 <- data.frame(
        CHR = c(1, 1),
        BP = c(100, 200),
        SNP = c("rs1", "rs2"),
        CM = 0,
        ann1 = c(1, 0),
        stringsAsFactors = FALSE
    )
    annot_df2 <- data.frame(
        CHR = c(1, 1),
        BP = c(100, 200),
        SNP = c("rs1", "rs2"),
        CM = 0,
        ann2 = c(0, 1),
        stringsAsFactors = FALSE
    )

    file1 <- tempfile(fileext = ".annot")
    file2 <- tempfile(fileext = ".annot")
    on.exit(unlink(c(file1, file2)), add = TRUE)
    write.table(annot_df1, file1, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(annot_df2, file2, sep = "\t", row.names = FALSE, quote = FALSE)

    snp_gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(100, 200), width = 1)
    )

    paths <- c(ann1 = file1, ann2 = file2)
    result <- readAnnotations(paths, snp_gr)

    expect_s4_class(result, "AnnotationMatrix")
    expect_equal(ncol(assay(result, "annotations")), 2)
    expect_equal(colnames(assay(result, "annotations")), c("ann1", "ann2"))
    expect_equal(as.numeric(assay(result, "annotations")[, "ann1"]), c(1, 0))
    expect_equal(as.numeric(assay(result, "annotations")[, "ann2"]), c(0, 1))
})

# =============================================================================
# BigWig reader, LDSC CHR/BP guard
# =============================================================================

test_that(".readBigwigAtSnps returns the mean BigWig score at each SNP", {
    skip_if_not_installed("rtracklayer")
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(1, 201, 401), end = c(200, 400, 600)),
        score = c(1.5, 2.5, 3.5)
    )
    GenomeInfoDb::seqlengths(gr) <- c(chr1 = 1000L)
    bw <- tempfile(fileext = ".bw")
    on.exit(unlink(bw), add = TRUE)
    rtracklayer::export(gr, bw, format = "bigWig")

    snp_gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(100, 300, 500, 800), width = 1)
    )
    result <- pecotmr:::.readBigwigAtSnps(bw, snp_gr)
    # Positions 100/300/500 fall in the three scored intervals; 800 has no
    # coverage and defaults to 0.
    expect_equal(result, c(1.5, 2.5, 3.5, 0))
})

test_that("readAnnotations reads a BigWig file as a continuous annotation", {
    skip_if_not_installed("rtracklayer")
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(1, 201, 401), end = c(200, 400, 600)),
        score = c(1.5, 2.5, 3.5)
    )
    GenomeInfoDb::seqlengths(gr) <- c(chr1 = 1000L)
    bw <- tempfile(fileext = ".bw")
    on.exit(unlink(bw), add = TRUE)
    rtracklayer::export(gr, bw, format = "bigWig")

    snp_gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(100, 300, 800), width = 1)
    )
    result <- readAnnotations(c(cons = bw), snp_gr)
    expect_s4_class(result, "AnnotationMatrix")
    # A .bw file is auto-detected as a continuous annotation.
    expect_equal(SummarizedExperiment::colData(result)$type, "continuous")
    expect_equal(as.numeric(assay(result, "annotations")[, 1]), c(1.5, 2.5, 0))
})

test_that(".readLdscAnnot errors when CHR/BP columns are absent", {
    # The annotation column is present but the positional CHR/BP columns are not.
    annot_df <- data.frame(
        SNP = c("rs1", "rs2"),
        CM = 0,
        my_annot = c(1, 0),
        stringsAsFactors = FALSE
    )
    annot_file <- tempfile(fileext = ".annot")
    on.exit(unlink(annot_file), add = TRUE)
    write.table(
        annot_df,
        annot_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    snp_gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = 100, width = 1)
    )
    expect_error(
        pecotmr:::.readLdscAnnot(annot_file, snp_gr, "my_annot"),
        "must contain CHR and BP columns"
    )
})

# =============================================================================
# fineMappingAnnotations (h2Annotations.R)
# =============================================================================

# An LD-score reference over `n` variants at chr1:100, 200, ..., named by
# rsID, in `nBlocks` blocks with identity LD. `flip` lists variants it
# stores with the alleles the other way round.
.fma_ref <- function(n = 6L, nBlocks = 2L, flip = integer(0)) {
    a1 <- rep("G", n)
    a2 <- rep("A", n)
    a1[flip] <- "A"
    a2[flip] <- "G"
    per <- n / nBlocks
    LdScore(
        ldBlocks = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(
                start = (seq_len(nBlocks) - 1L) * per * 100L + 1L,
                end = seq_len(nBlocks) * per * 100L
            )
        ),
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(n)),
            CHR = "chr1",
            BP = seq_len(n) * 100L,
            A1 = a1,
            A2 = a2
        ),
        nSamples = 1000L,
        inSample = FALSE,
        genome = "hg19",
        ldScores = matrix(1, n, 1, dimnames = list(NULL, "base_l2")),
        ldScoreWeights = rep(1, n),
        ldMatrixList = lapply(seq_len(nBlocks), function(b) {
            list(
                R = diag(per),
                snpIdx = as.integer((b - 1L) * per + seq_len(per))
            )
        })
    )
}

# One gene's row over variants `idx` (ids chr1:<pos>:A:G) with PIPs `pip`,
# 95% credible-set members `cs` and their purity.
.fma_row <- function(idx, pip, cs, purity) {
    ids <- sprintf("chr1:%d:A:G", idx * 100L)
    FineMappingRow(
        variantIds = ids,
        susieFit = list(pip = stats::setNames(pip, ids)),
        topLoci = data.frame(
            variant_id = ids,
            pip = pip,
            cs_95 = ifelse(seq_along(idx) %in% cs, "susie_1", "susie_0"),
            cs_95_purity = ifelse(seq_along(idx) %in% cs, purity, 0)
        )
    )
}

# Two genes of one context sharing variants 3 and 4.
.fma_qtl <- function() {
    QtlFineMappingResult(
        studyName = c("Q", "Q"),
        context = c("ctx", "ctx"),
        trait = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(
            .fma_row(1:4, c(0.7, 0.2, 0.05, 0.05), cs = 1:2, purity = 0.9),
            .fma_row(3:6, c(0.6, 0.1, 0.2, 0.1), cs = c(1L, 3L), purity = 0.3)
        )
    )
}

.fma_col <- function(annot, type) {
    unname(SummarizedExperiment::assay(annot, "annotations")[,
        paste0("Q|ctx|susie|", type)
    ])
}

test_that("fineMappingAnnotations builds max and summed PIP over the traits", {
    annot <- fineMappingAnnotations(.fma_qtl(), .fma_ref())
    # Variant 3 carries both genes: 0.05 and 0.6.
    expect_equal(.fma_col(annot, "maxCpp"), c(0.7, 0.2, 0.6, 0.1, 0.2, 0.1))
    expect_equal(.fma_col(annot, "sumPip"), c(0.7, 0.2, 0.65, 0.15, 0.2, 0.1))
    expect_equal(.fma_col(annot, "highPip"), c(1, 0, 1, 0, 0, 0))
    expect_equal(.fma_col(annot, "credibleSet"), c(1, 1, 1, 0, 1, 0))
})

test_that("fineMappingAnnotations rows follow the reference", {
    ref <- .fma_ref()
    annot <- fineMappingAnnotations(.fma_qtl(), ref)
    expect_identical(nrow(annot), length(ref))
    expect_identical(
        GenomicRanges::start(SummarizedExperiment::rowRanges(annot)),
        GenomicRanges::start(ref)
    )
    meta <- as.data.frame(SummarizedExperiment::colData(annot))
    expect_identical(
        meta$name,
        paste0("Q|ctx|susie|", c("maxCpp", "sumPip", "credibleSet", "highPip"))
    )
    expect_identical(unique(meta$tier), "candidate")
    expect_identical(
        meta$type,
        c("continuous", "continuous", "binary", "binary")
    )
    baseline <- fineMappingAnnotations(.fma_qtl(), ref, tier = "baseline")
    expect_identical(
        unique(SummarizedExperiment::colData(baseline)$tier),
        "baseline"
    )
})

test_that("fineMappingAnnotations matches a reference variant stored flipped", {
    annot <- fineMappingAnnotations(.fma_qtl(), .fma_ref(flip = 1L))
    expect_equal(.fma_col(annot, "maxCpp")[[1L]], 0.7)
})

test_that("fineMappingAnnotations takes the coverage, purity and threshold", {
    qtl <- .fma_qtl()
    ref <- .fma_ref()
    pure <- fineMappingAnnotations(
        qtl,
        ref,
        types = "credibleSet",
        minPurity = 0.5
    )
    expect_equal(.fma_col(pure, "credibleSet"), c(1, 1, 0, 0, 0, 0))
    low <- fineMappingAnnotations(
        qtl,
        ref,
        types = "highPip",
        pipThreshold = 0.2
    )
    expect_equal(.fma_col(low, "highPip"), c(1, 1, 1, 0, 1, 0))
    # No credible sets were stored at 70% coverage.
    expect_warning(
        none <- fineMappingAnnotations(
            qtl,
            ref,
            types = c("maxCpp", "credibleSet"),
            coverage = 0.7
        ),
        "dropping 1 column"
    )
    expect_identical(colnames(none), "Q|ctx|susie|maxCpp")
})

test_that("fineMappingAnnotations refuses bad arguments", {
    expect_error(
        fineMappingAnnotations(.fma_qtl(), .fma_ref(), types = "meanPip"),
        "must be one of"
    )
    expect_error(
        fineMappingAnnotations(.fma_qtl(), .fma_ref(), tier = "tested"),
        "must be one of"
    )
    expect_error(fineMappingAnnotations(.fma_qtl(), list()), "LdStatistic")
})

test_that("estimateH2 takes fineMappingAnnotations as candidates", {
    n <- 40L
    ref <- .fma_ref(n = n, nBlocks = 4L)
    set.seed(5)
    df <- data.frame(
        SNP = paste0("rs", seq_len(n)),
        CHR = "chr1",
        BP = seq_len(n) * 100L,
        A1 = "G",
        A2 = "A",
        Z = stats::rnorm(n, sd = 1.5),
        N = 50000
    )
    ss <- GwasSumStats(
        studyName = "t",
        entry = list(.dfToGwasGr(df)),
        genome = "hg19",
        ldReference = makeTestGwasGenotypeHandle()
    )
    qtl <- QtlFineMappingResult(
        "Q",
        "ctx",
        "g1",
        "susie",
        list(.fma_row(
            seq(2L, n, by = 3L),
            rep(0.3, length(seq(2L, n, by = 3L))),
            cs = 1:3,
            purity = 0.9
        ))
    )
    annot <- fineMappingAnnotations(
        qtl,
        ref,
        types = c("maxCpp", "credibleSet")
    )
    for (method in c("sldsc", "gldsc")) {
        z <- scoreStats(estimateH2(
            ss,
            ref,
            method = method,
            annotations = annot
        ))$z
        expect_length(z, 2L)
        expect_true(all(is.finite(z)), info = method)
    }
})
