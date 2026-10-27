# Fine-mapping and GWAS results that the writer, reader and loader tests
# publish and read back.

# One block's fit: three variants from `offset`, its effects' V and the
# effects that formed credible sets.
.rw_blockFit <- function(offset, V, csIndex) {
    FineMappingRow(
        variantIds = paste0("chr1:", offset + 100 * (1:3), ":A:G"),
        susieFit = list(V = V, sets = list(cs_index = csIndex)),
        topLoci = .sc_makeTopLoci(3, offset = offset)
    )
}

.rw_blockedGwasFmr <- function(
    dataType = NULL,
    blockId = c("chr1_1_5000", "chr1_5001_20000")
) {
    GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(
            .rw_blockFit(0L, c(0.1, 0.2), 1L),
            .rw_blockFit(10000L, c(0.3, 0.4), 2L)
        ),
        blockId = blockId,
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(1L, 5001L), c(5000L, 20000L))
        ),
        dataType = dataType
    )
}

# The same two blocks of one fit, the first given its window and the second
# not, so its window is its variants' span.
.rw_mixedGwasFmr <- function() {
    combineFineMappingResults(
        GwasFineMappingResult(
            studyName = "g1",
            method = "susie",
            entry = list(.rw_blockFit(0L, c(0.1, 0.2), 1L)),
            blockId = "chr1_1_5000",
            region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 5000L))
        ),
        GwasFineMappingResult(
            studyName = "g1",
            method = "susie",
            entry = list(.rw_blockFit(10000L, c(0.3, 0.4), 2L)),
            blockId = "b2"
        )
    )
}

.rw_sampleLines <- function(path) {
    grep("^##SAMPLE", readLines(path), value = TRUE)
}

# @noRd
.ezss <- function(supplied) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100L, 200L), width = 1L)
    )
    cols <- S4Vectors::DataFrame(
        SNP = c("chr1:100:A:G", "chr1:200:C:T"),
        A1 = c("G", "T"),
        A2 = c("A", "C"),
        Z = c(2.5, -1.2),
        N = c(1000L, 1000L),
        BETA = c(0.25, -0.12),
        SE = c(0.1, 0.1)
    )
    if (!is.null(supplied)) {
        cols$z_supplied <- rep(supplied, 2L)
    }
    GwasSumStats(
        studyName = "s",
        entry = list(S4Vectors::`mcols<-`(gr, value = cols)),
        genome = "hg38"
    )
}
