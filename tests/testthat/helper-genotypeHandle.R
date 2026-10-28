# The internal handle reader.
#
# readGenotypes() is the exported entry point and returns a genotype panel: a
# RangedSummarizedExperiment whose dosage assay reads lazily through a handle.
# The handle is the seed layer beneath it and is not public.
#
# Tests that exercise the reader machinery itself -- format detection, snpInfo,
# sample ids, sharded routing -- assert on the handle, so they read one
# directly rather than unwrapping a panel on every line.
# @noRd
readGenotypeHandle <- function(path, format = NULL, ...) {
    pecotmr:::.readGenotypeHandle(path, format = format, ...)
}

.gh_pgenHandle <- function() {
    skip_if_not_installed("pgenlibr")
    readGenotypeHandle(test_path("test_data", "test_variants.pgen"))
}

.sk_panel <- function() {
    readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink1"
    )
}

.ssQ_makeHandle <- function(snp_n = 8L, n_samples = 60L) {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("chr1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}

# A panel whose SNP ids follow the formatVariantId convention (chr:pos:A2:A1),
# as real genotype LD sketches do, so a re-keyed entry SNP resolves against it.
.ssQ_makeHandleVid <- function(snp_n = 8L, n_samples = 60L) {
    h <- .ssQ_makeHandle(snp_n, n_samples)
    h@snpInfo$SNP <- paste0("chr1:", snpInfo(h)$BP, ":G:A")
    h
}

.dose <- function(panel) {
    unname(as.matrix(SummarizedExperiment::assay(panel, "dosage")))
}
