context("file_utils")
library(tidyverse)

test_that("readPvar reads real pvar file", {
    skip_if_not_installed("pgenlibr")
    pvar_path <- file.path(test_path("test_data"), "test_variants.pvar")
    res <- pecotmr:::readPvar(pvar_path)
    expect_equal(colnames(res), c("chrom", "id", "pos", "A2", "A1"))
    expect_equal(nrow(res), 349L)
    expect_true(all(res$chrom == "21"))
})

test_that("Test loadGenotypeRegion", {
    res <- loadGenotypeRegion(
        "test_data/protocol_example.genotype"
    )
    sample_ids <- read_delim(
        "test_data/protocol_example.genotype.fam",
        delim = "\t",
        col_names = F
    ) %>%
        pull(X1)
    expect_equal(nrow(res), length(sample_ids))
    expect_equal(rownames(res), sample_ids)
})

test_that("Test loadGenotypeRegion no indels", {
    res <- loadGenotypeRegion(
        "test_data/protocol_example.genotype",
        keepIndel = F
    )
    bim_file <- read_delim(
        "test_data/protocol_example.genotype.bim",
        delim = "\t",
        col_names = F
    )
    sample_ids <- read_delim(
        "test_data/protocol_example.genotype.fam",
        delim = "\t",
        col_names = F
    ) %>%
        pull(X1)
    expect_equal(nrow(res), length(sample_ids))
    expect_equal(rownames(res), sample_ids)
    indels <- with(
        bim_file,
        grepl("[^ATCG]", X5) |
            grepl("[^ATCG]", X6) |
            nchar(X5) > 1 |
            nchar(X6) > 1
    )
    expect_equal(
        nrow(bim_file[!indels, ]),
        ncol(res)
    )
})

test_that("Test loadGenotypeRegion with region", {
    res <- loadGenotypeRegion(
        "test_data/protocol_example.genotype",
        which = GenomicRanges::GRanges("chr22:20689453-20845958")
    )
    sample_ids <- read_delim(
        "test_data/protocol_example.genotype.fam",
        delim = "\t",
        col_names = F
    ) %>%
        pull(X1)
    snp_ids <- read_delim(
        "test_data/protocol_example.genotype.bim",
        delim = "\t",
        col_names = F
    ) %>%
        pull(X2)
    expect_equal(nrow(res), length(sample_ids))
    expect_equal(rownames(res), sample_ids)
    expect_equal(ncol(res), 8)
    expect_equal(colnames(res), snp_ids[1:8])
})

test_that("Test loadGenotypeRegion with region and no indels", {
    res <- loadGenotypeRegion(
        "test_data/protocol_example.genotype",
        which = GenomicRanges::GRanges("chr22:20689453-20845958"),
        keepIndel = F
    )
    bim_file <- read_delim(
        "test_data/protocol_example.genotype.bim",
        delim = "\t",
        col_names = F
    )[1:8, ]
    sample_ids <- read_delim(
        "test_data/protocol_example.genotype.fam",
        delim = "\t",
        col_names = F
    ) %>%
        pull(X1)
    expect_equal(nrow(res), length(sample_ids))
    expect_equal(rownames(res), sample_ids)
    indels <- with(
        bim_file,
        grepl("[^ATCG]", X5) |
            grepl("[^ATCG]", X6) |
            nchar(X5) > 1 |
            nchar(X6) > 1
    )
    expect_equal(
        nrow(bim_file[!indels, ]),
        ncol(res)
    )
    expect_equal(colnames(res), bim_file[!indels, ]$X2)
})

# --- findStochasticMeta tests ---

test_that("findStochasticMeta finds the sidecar beside a file or prefix", {
    dir <- withr::local_tempdir()
    file.create(file.path(dir, "sketch.afreq"))
    files <- c(
        "sketch",
        "sketch.bgen",
        "sketch.gds",
        "sketch.vcf.gz",
        "sketch.vcf.bgz",
        "sketch.VCF.GZ"
    )
    for (file in files) {
        expect_equal(
            pecotmr:::findStochasticMeta(file.path(dir, file)),
            file.path(dir, "sketch.afreq"),
            info = file
        )
    }
})

test_that("findStochasticMeta returns NULL when no sidecar exists", {
    td <- test_path("test_data")
    result <- pecotmr:::findStochasticMeta(file.path(
        td,
        "protocol_example.genotype"
    ))
    expect_null(result)
})

# --- readStochasticMeta tests ---

test_that("readStochasticMeta reads generic format", {
    td <- test_path("test_data")
    path <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
    result <- pecotmr:::readStochasticMeta(path)
    expect_true(is.data.frame(result))
    expect_equal(colnames(result), c("id", "u_min", "u_max"))
    expect_equal(nrow(result), 8L)
    expect_true(is.numeric(result$u_min))
    expect_true(is.numeric(result$u_max))
})

test_that("readStochasticMeta reads afreq format", {
    td <- test_path("test_data")
    path <- file.path(td, "test_harmonize_regions.afreq")
    result <- pecotmr:::readStochasticMeta(path)
    expect_true(is.data.frame(result))
    expect_equal(colnames(result), c("id", "u_min", "u_max"))
    expect_equal(nrow(result), 8L)
    expect_true(all(grepl("^chr21_", result$id)))
    # The generic sidecar beside it holds the same values.
    generic <- pecotmr:::readStochasticMeta(
        file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
    )
    expect_equal(result$id, generic$id)
    expect_equal(result$u_min, generic$u_min)
    expect_equal(result$u_max, generic$u_max)
})

test_that("readStochasticMeta reads afreq.zst format", {
    td <- test_path("test_data")
    path <- file.path(td, "test_harmonize_regions.afreq.zst")
    result <- pecotmr:::readStochasticMeta(path)
    expect_true(is.data.frame(result))
    expect_equal(colnames(result), c("id", "u_min", "u_max"))
    expect_equal(nrow(result), 8L)
    # Should produce identical results to the plain afreq
    plain <- pecotmr:::readStochasticMeta(file.path(
        td,
        "test_harmonize_regions.afreq"
    ))
    expect_equal(result, plain)
})

test_that("findStochasticMeta prefers afreq over afreq.zst", {
    td <- test_path("test_data")
    # Both .afreq and .afreq.zst exist; findStochasticMeta should return .afreq first
    result <- pecotmr:::findStochasticMeta(file.path(
        td,
        "test_harmonize_regions"
    ))
    expect_true(grepl("\\.afreq$", result))
})

test_that("readStochasticMeta reads a forced afreq file under any name", {
    # The override reads the file it is given; it once appended .afreq to
    # the name and read <path>.afreq instead, so nothing came back.
    src <- test_path("test_data", "test_harmonize_regions.afreq")
    renamed <- file.path(withr::local_tempdir(), "sidecar.txt")
    file.copy(src, renamed)
    expect_equal(
        pecotmr:::readStochasticMeta(renamed, format = "afreq"),
        pecotmr:::readStochasticMeta(src)
    )
    expect_named(
        pecotmr:::readStochasticMeta(renamed, format = "afreq"),
        c("id", "u_min", "u_max")
    )
})

test_that("readStochasticMeta respects format override", {
    # A generic id/u_min/u_max table saved under an .afreq name: left to its
    # extension it goes to the PLINK2 .afreq reader, which cannot parse it.
    generic <- test_path(
        "test_data",
        "test_harmonize_regions.stochastic_meta.tsv"
    )
    path <- file.path(withr::local_tempdir(), "sidecar.afreq")
    file.copy(generic, path)
    expect_error(pecotmr:::readStochasticMeta(path), "#CHROM")
    # Forcing "generic" reads it as the table it is.
    expect_equal(
        pecotmr:::readStochasticMeta(path, format = "generic"),
        pecotmr:::readStochasticMeta(generic)
    )
})

test_that("readStochasticMeta returns NULL for afreq without U_MIN/U_MAX", {
    td <- test_path("test_data")
    # test_variants.afreq has no U_MIN/U_MAX columns
    path <- file.path(td, "test_variants.afreq")
    result <- pecotmr:::readStochasticMeta(path)
    expect_null(result)
})

test_that("readStochasticMeta returns NULL for nonexistent file", {
    result <- pecotmr:::readStochasticMeta("/nonexistent/file.tsv")
    expect_null(result)
})

# --- loadGenotypeRegion stochastic inversion test ---

test_that("loadGenotypeRegion applies stochastic inversion with explicit sidecar", {
    td <- test_path("test_data")
    metaPath <- file.path(td, "test_harmonize_regions.stochastic_meta.tsv")
    smeta <- pecotmr:::readStochasticMeta(metaPath)

    # Load with explicit sidecar - inversion transforms the integer dosages
    res <- loadGenotypeRegion(
        file.path(td, "test_harmonize_regions"),
        returnVariantInfo = TRUE,
        stochasticMetaPath = metaPath
    )

    expect_equal(ncol(res$X), 8L)
    # u_min/u_max should be attached to variant_info
    expect_true("u_min" %in% colnames(res$variant_info))
    expect_true("u_max" %in% colnames(res$variant_info))
    expect_equal(res$variant_info$u_min, smeta$u_min)
    expect_equal(res$variant_info$u_max, smeta$u_max)

    # Verify inversion math: for a dosage value d with u_min/u_max,
    # inverted = d * (u_max - u_min) / 2 + u_min
    # Check the first variant's first sample manually
    raw <- loadGenotypeRegion(
        file.path(td, "protocol_example.genotype"),
        which = GenomicRanges::GRanges("chr22:20689453-20845958")
    )
    # protocol_example has no sidecar, so raw values are unchanged (integer dosages)
    expect_true(all(raw == round(raw), na.rm = TRUE))

    # The inverted matrix should NOT be all integers (u_min != 0 or u_max != 2)
    expect_false(all(res$X == round(res$X), na.rm = TRUE))
})

# ===========================================================================
# readBim vroom-based tests
# ===========================================================================

test_that("readBim returns correct columns and types", {
    bim_path <- tempfile(fileext = ".bim")
    cat("22\trs100\t0\t50000\tA\tG\n", file = bim_path)
    cat("22\trs200\t0\t60000\tT\tC\n", file = bim_path, append = TRUE)
    cat("22\trs300\t0\t70000\tC\tA\n", file = bim_path, append = TRUE)

    bed_path <- sub("\\.bim$", ".bed", bim_path)
    file.copy(bim_path, bim_path)
    res <- readBim(bed_path)
    expect_equal(nrow(res), 3)
    expect_equal(colnames(res), c("chrom", "id", "gpos", "pos", "a1", "a0"))
    expect_equal(res$id, c("rs100", "rs200", "rs300"))
    expect_equal(res$pos, c(50000, 60000, 70000))
    file.remove(bim_path)
})

# ===========================================================================
# .colocFilterCsByConcentration
# ===========================================================================

# ---- invertMinmaxScaling ----

test_that("invertMinmaxScaling exactly recovers original U", {
    set.seed(42)
    n <- 500
    k <- 4
    # Simulate original U with arbitrary values
    U_original <- matrix(rnorm(n * k, mean = 0.5, sd = 0.3), n, k)

    # Apply the same min-max scaling as rss_ld_sketch
    u_min <- apply(U_original, 2, min)
    u_max <- apply(U_original, 2, max)
    denom <- u_max - u_min
    U_scaled <- sweep(sweep(U_original, 2, u_min, "-"), 2, denom, "/") * 2

    # Verify scaled is in [0, 2]
    expect_true(all(U_scaled >= 0 & U_scaled <= 2))

    # Invert
    U_recovered <- invertMinmaxScaling(U_scaled, u_min, u_max)

    # Must be exactly the original (up to floating point)
    expect_equal(U_recovered, U_original, tolerance = 1e-12)
})

test_that("invertMinmaxScaling handles monomorphic variant", {
    # Column 1 is monomorphic (u_min == u_max). The sketch writer divides it
    # by 1 instead of 0, which leaves it all 0, and inverting has to give its
    # constant back rather than NaN.
    U <- cbind(rep(0.5, 3), c(-1, 0, 3))
    u_min <- c(0.5, -1)
    u_max <- c(0.5, 3)
    X <- cbind(rep(0, 3), c(0, 0.5, 2)) # 2 * (U - u_min) / c(1, 4)
    expect_equal(invertMinmaxScaling(X, u_min, u_max), U)
})

test_that("invertMinmaxScaling errors on mismatched lengths", {
    X <- matrix(1:6, ncol = 2)
    expect_error(
        invertMinmaxScaling(X, c(0, 0, 0), c(1, 1, 1)),
        "Length of u_min"
    )
})

# ===========================================================================
# refVariantInfo
# ===========================================================================

# refVariantInfo() through an LD-metadata file naming `source` for all of
# chr21, written into test_data so its relative path resolves, and removed
# once read.
.gioRefInfoViaMeta <- function(source, which) {
    meta <- withr::local_tempfile(
        tmpdir = test_path("test_data"),
        fileext = ".tsv"
    )
    writeLines(c("chrom\tstart\tend\tpath", str_c("21\t0\t0\t", source)), meta)
    refVariantInfo(meta, which = which)
}

test_that("refVariantInfo reads a genotype file or prefix given directly", {
    # It used to read every source as LD metadata, so a direct path failed
    # whatever its format.
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    stem <- str_remove(
        system.file("extdata", "toy_ref.bed", package = "pecotmr"),
        "\\.bed$"
    )
    bim <- read.table(str_c(stem, ".bim"))
    plink1 <- refVariantInfo(stem)
    expect_equal(nrow(plink1), nrow(bim))
    expect_equal(plink1$pos, bim[[4L]])
    gds <- refVariantInfo(test_path("test_data", "test_variants.gds"))
    expect_equal(nrow(gds), 349L)
})

test_that("refVariantInfo reads PLINK2, VCF and GDS sources via LD metadata", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    # A window over the first 8 of the 349 variants. PLINK2 reads the allele
    # frequency from its .afreq sidecar, VCF and GDS from the dosages.
    window <- GenomicRanges::GRanges("chr21:17513228-17515000")
    sources <- c(
        plink2 = "test_variants",
        vcf = "test_variants.vcf.gz",
        gds = "test_variants.gds"
    )
    infos <- map(sources, .gioRefInfoViaMeta, which = window)
    for (format in names(infos)) {
        info <- infos[[format]]
        expect_true(is.data.frame(info), info = format)
        expect_true(
            all(
                c("chrom", "id", "pos", "A2", "A1", "allele_freq") %in%
                    names(info)
            ),
            info = format
        )
        expect_true(
            all(info$allele_freq > 0 & info$allele_freq < 1),
            info = format
        )
        expect_equal(nrow(info), 8L, info = format)
        expect_true(
            all(info$pos >= 17513228 & info$pos <= 17515000),
            info = format
        )
        expect_equal(info$pos, infos$plink2$pos, info = format)
    }
})

# ===========================================================================
# readAfreq
# ===========================================================================

test_that("readAfreq returns correct structure from .afreq file", {
    td <- test_path("test_data")
    af <- readAfreq(file.path(td, "test_variants"))
    expect_true(is.data.frame(af))
    expect_equal(nrow(af), 349L)
    expect_true(all(
        c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)
    ))
})

test_that("readAfreq returns correct types", {
    td <- test_path("test_data")
    af <- readAfreq(file.path(td, "test_variants"))
    expect_type(af$alt_freq, "double")
    expect_true(all(af$alt_freq >= 0 & af$alt_freq <= 1))
    expect_true(all(af$obs_ct > 0))
})

test_that("readAfreq returns NULL when no afreq file exists", {
    af <- readAfreq(file.path(tempdir(), "nonexistent_prefix"))
    expect_null(af)
})

test_that("readAfreq reads .afreq.zst in-process, preferring it to .afreq", {
    td <- test_path("test_data")
    # Empty PATH: no zstd (or any) CLI is reachable. archive decompresses the
    # .zst in-process via libarchive, so the read must still succeed.
    withr::local_envvar(PATH = "")
    # test_harmonize_regions has both .afreq and .afreq.zst; readAfreq
    # prefers the .zst.
    af <- readAfreq(file.path(td, "test_harmonize_regions"))
    expect_true(is.data.frame(af))
    expect_true(all(
        c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct") %in% colnames(af)
    ))
    # This afreq has U_MIN/U_MAX columns
    expect_true(all(c("u_min", "u_max") %in% colnames(af)))
    expect_equal(nrow(af), 8L)
})

test_that("readAfreq reads plain .afreq with U_MIN/U_MAX", {
    # A copy with no .zst beside it, so readAfreq falls through to plain
    # .afreq without the shared fixture being touched.
    dir <- withr::local_tempdir()
    file.copy(
        test_path("test_data", "test_harmonize_regions.afreq"),
        file.path(dir, "g.afreq")
    )
    af <- readAfreq(file.path(dir, "g"))
    expect_true(is.data.frame(af))
    expect_true(all(c("u_min", "u_max") %in% colnames(af)))
    expect_equal(nrow(af), 8L)
})

test_that("readAfreq IDs match pvar IDs", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    af <- readAfreq(file.path(td, "test_variants"))
    pvar <- readPvar(file.path(td, "test_variants.pvar"))
    expect_equal(af$id, pvar$id)
})

# ===========================================================================
# matchVariantsToKeep
# ===========================================================================

test_that("matchVariantsToKeep filters to specified variants", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    vi <- pecotmr:::.snpInfoToVariantInfo(snpInfo(handle))

    # Write a keep file as tab-delimited with chrom/pos columns
    keep_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(keep_file), add = TRUE)
    keep_df <- vi[c(1, 5, 10), c("chrom", "pos", "A2", "A1")]
    vroom::vroom_write(keep_df, keep_file, delim = "\t")

    mask <- pecotmr:::matchVariantsToKeep(vi, keep_file)
    expect_type(mask, "logical")
    expect_equal(sum(mask), 3L)
    expect_true(mask[1])
    expect_true(mask[5])
    expect_true(mask[10])
})

test_that("matchVariantsToKeep returns all FALSE for non-matching variants", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    vi <- pecotmr:::.snpInfoToVariantInfo(snpInfo(handle))

    keep_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(keep_file), add = TRUE)
    keep_df <- data.frame(
        chrom = c(1L, 2L),
        pos = c(999L, 888L),
        A2 = c("A", "C"),
        A1 = c("T", "G")
    )
    vroom::vroom_write(keep_df, keep_file, delim = "\t")

    mask <- pecotmr:::matchVariantsToKeep(vi, keep_file)
    expect_true(all(!mask))
})

# ===========================================================================
# matchVariantsToKeep (additional coverage)
# ===========================================================================

test_that("matchVariantsToKeep works with single-column variant ID file", {
    # The variant info is pecotmr's own (UCSC); the keep file may spell its
    # chromosomes either way.
    vi <- data.frame(
        chrom = c("chr1", "chr1", "chr1"),
        pos = c(100L, 200L, 300L),
        A2 = c("A", "C", "G"),
        A1 = c("G", "T", "A"),
        stringsAsFactors = FALSE
    )
    keep_file <- tempfile(fileext = ".txt")
    on.exit(unlink(keep_file), add = TRUE)
    writeLines(c("1:100:A:G", "1:300:G:A"), keep_file)

    mask <- pecotmr:::matchVariantsToKeep(vi, keep_file)
    expect_type(mask, "logical")
    expect_equal(sum(mask), 2L)
    expect_true(mask[1])
    expect_true(mask[3])
})

test_that("matchVariantsToKeep uses position-only matching when no alleles", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    vi <- pecotmr:::.snpInfoToVariantInfo(snpInfo(handle))

    keep_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(keep_file), add = TRUE)
    # Write keep file with chrom/pos only (no alleles)
    keep_df <- vi[c(1, 5), c("chrom", "pos")]
    vroom::vroom_write(keep_df, keep_file, delim = "\t")

    mask <- pecotmr:::matchVariantsToKeep(vi, keep_file)
    expect_type(mask, "logical")
    expect_equal(sum(mask), 2L)
    expect_true(mask[1])
    expect_true(mask[5])
})

# ===========================================================================
# readGenotypes + extractBlockGenotypes: plink2 tests
# ===========================================================================

test_that("readGenotypeHandle reads a plink2 fileset with all its variants", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    expect_s4_class(handle, "GenotypeHandle")
    expect_equal(genotypeFormat(handle), "plink2")
    expect_equal(nSamples(handle), 100L)
    expect_equal(nrow(snpInfo(handle)), 349L)
    rse <- extractBlockGenotypes(handle, seq_len(nrow(snpInfo(handle))))
    expect_s4_class(rse, "SummarizedExperiment")
    dosage <- SummarizedExperiment::assay(rse, "dosage")
    expect_equal(nrow(dosage), 349L)
    expect_equal(ncol(dosage), 100L)
})

test_that("loadGenotypeRegion filters by keep_variants_path for plink2", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    vi <- pecotmr:::.snpInfoToVariantInfo(snpInfo(handle))

    keep_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(keep_file), add = TRUE)
    keep_df <- vi[c(1, 3, 7), c("chrom", "pos", "A2", "A1")]
    vroom::vroom_write(keep_df, keep_file, delim = "\t")

    result <- loadGenotypeRegion(
        file.path(td, "test_variants"),
        keepVariantsPath = keep_file
    )
    expect_equal(ncol(result), 3L)
})

test_that("loadGenotypeRegion attaches afreq info for plink2", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    result <- loadGenotypeRegion(
        file.path(td, "test_variants"),
        returnVariantInfo = TRUE
    )
    vi <- result$variant_info
    expect_true("alt_freq" %in% colnames(vi))
    expect_true("obs_ct" %in% colnames(vi))
    expect_true(all(vi$alt_freq >= 0 & vi$alt_freq <= 1))
})

test_that("readGenotypes plink2 sample names match psam IIDs", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    expect_true(all(grepl("^(HG|NA)\\d+", sampleIds(handle))))
    expect_equal(length(unique(sampleIds(handle))), 100L)
})

# Tests for genotype loading via readGenotypes + extractBlockGenotypes,
# and the loadGenotypeRegion dispatcher.

# Fixtures: 100 samples x 349 variants on chr21:17513228-17592874
test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")
vcf_path <- file.path(test_data_dir, "test_variants.vcf.gz")
gds_path <- file.path(test_data_dir, "test_variants.gds")
# A window over the first 8 variants.
region_sub <- GenomicRanges::GRanges("chr21:17513228-17515000")

# Expected dimensions
n_samples <- 100L
n_variants <- 349L

# Shared helper: validate the output structure from loadGenotypeRegion
# (with returnVariantInfo=TRUE)
checkGenotypeResult <- function(
    result,
    expected_nrow = n_samples,
    expected_ncol = n_variants,
    label = ""
) {
    expect_true(is.list(result), label = paste(label, "is list"))
    expect_named(result, c("X", "variant_info"), ignore.order = TRUE)
    expect_true(is.matrix(result$X))
    expect_true(is.numeric(result$X))
    expect_equal(nrow(result$X), expected_nrow)
    expect_equal(ncol(result$X), expected_ncol)
    expect_true(is.data.frame(result$variant_info))
    expect_true(all(
        c("chrom", "id", "pos", "A2", "A1") %in% names(result$variant_info)
    ))
    expect_equal(nrow(result$variant_info), expected_ncol)
    # Column names of X match variant IDs
    expect_equal(colnames(result$X), result$variant_info$id)
    # Dosage values should be non-negative integers (0, 1, 2 for biallelic;
    # multiallelic VCF sites can have higher values)
    vals <- result$X[!is.na(result$X)]
    expect_true(all(vals >= 0), label = paste(label, "dosage non-negative"))
    expect_true(
        all(vals == round(vals)),
        label = paste(label, "dosage integer-valued")
    )
}

# --- readGenotypes: PLINK1, VCF and GDS panels -------------------------------

test_that("readGenotypes reads PLINK1, VCF and GDS files as panels", {
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    sources <- list(
        plink1 = list(plink1Prefix = plink_prefix),
        vcf = list(vcf_path),
        gds = list(gds_path)
    )
    for (format in names(sources)) {
        panel <- exec(readGenotypes, !!!sources[[format]])
        expect_s4_class(panel, "RangedSummarizedExperiment")
        expect_equal(dim(panel), c(n_variants, n_samples), info = format)
        expect_equal(
            S4Vectors::metadata(panel)$genotypeSource$format,
            format,
            info = format
        )
        # No handle: the panel reads through its own GenotypeDosageSeed, and
        # records where it came from itself.
        expect_s4_class(pecotmr:::.panelSeed(panel), "GenotypeDosageSeed")
        expect_null(pecotmr:::.ldRefHandle(panel), info = format)
    }
})

# --- loadGenotypeRegion -------------------------------------------------------

test_that("loadGenotypeRegion reads every format by its path, alike", {
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    # A PLINK2 prefix, a PLINK1 .bed, a VCF and a GDS file of the same
    # genotypes: each is told by its path and read to the same variants.
    paths <- c(
        plink2 = plink_prefix,
        plink1 = paste0(plink_prefix, ".bed"),
        vcf = vcf_path,
        gds = gds_path
    )
    loaded <- map(paths, loadGenotypeRegion, returnVariantInfo = TRUE)
    for (format in names(paths)) {
        checkGenotypeResult(loaded[[format]], label = format)
        expect_equal(
            loaded[[format]]$variant_info$pos,
            loaded$plink2$variant_info$pos,
            info = format
        )
    }
})

test_that("loadGenotypeRegion filters each format to a region", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    for (path in c(plink_prefix, vcf_path, gds_path)) {
        result <- loadGenotypeRegion(
            path,
            which = region_sub,
            returnVariantInfo = TRUE
        )
        checkGenotypeResult(result, expected_ncol = 8L, label = path)
        expect_true(
            all(
                result$variant_info$pos >= 17513228 &
                    result$variant_info$pos <= 17515000
            ),
            info = path
        )
    }
})

test_that("loadGenotypeRegion drops every format's indels", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    # The first 25 variants, 3 of them indels: a VCF is read a record at a
    # time, so the filter is run over a window rather than the whole file.
    window <- GenomicRanges::GRanges("chr21:17513228-17518122")
    for (path in c(plink_prefix, vcf_path, gds_path)) {
        result <- loadGenotypeRegion(
            path,
            which = window,
            keepIndel = FALSE,
            returnVariantInfo = TRUE
        )
        expect_equal(ncol(result$X), 22L, info = path)
        expect_true(all(nchar(result$variant_info$A1) == 1), info = path)
        expect_true(all(nchar(result$variant_info$A2) == 1), info = path)
    }
})

test_that("loadGenotypeRegion errors on an empty region", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("SNPRelate")
    skip_if_not_installed("gdsfmt")
    # Off the panel's variants, or on a chromosome it does not hold.
    for (path in c(plink_prefix, gds_path)) {
        for (window in c("chr21:1-2", "chr1:1-2")) {
            expect_error(
                loadGenotypeRegion(
                    path,
                    which = GenomicRanges::GRanges(window)
                ),
                "No SNPs found",
                class = "NoSnpsError",
                info = str_c(path, " ", window)
            )
        }
    }
})

test_that("loadGenotypeRegion returns matrix when returnVariantInfo=FALSE", {
    skip_if_not_installed("pgenlibr")
    result <- loadGenotypeRegion(plink_prefix)
    expect_true(is.matrix(result))
    expect_equal(nrow(result), n_samples)
    expect_equal(ncol(result), n_variants)
})

test_that("loadGenotypeRegion errors on unrecognized format", {
    expect_error(
        loadGenotypeRegion("/nonexistent/file.xyz"),
        "loadGenotypeRegion: .* is neither a genotype file nor a PLINK prefix"
    )
})

# --- Cross-format consistency -------------------------------------------------

test_that("PLINK1 and PLINK2 readGenotypes return consistent alleles", {
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("pgenlibr")
    # The fixture's two PLINK copies hold the same variants and calls, read
    # by different readers: PlinkMatrix and pgenlibr.
    p1 <- readGenotypes(plink1Prefix = plink_prefix)
    p2 <- readGenotypes(plink2Prefix = plink_prefix)
    rr1 <- S4Vectors::mcols(SummarizedExperiment::rowRanges(p1))
    rr2 <- S4Vectors::mcols(SummarizedExperiment::rowRanges(p2))
    expect_equal(rr1$A1, rr2$A1)
    expect_equal(rr1$A2, rr2$A2)
    expect_equal(
        as.matrix(SummarizedExperiment::assay(p1, "dosage")),
        as.matrix(SummarizedExperiment::assay(p2, "dosage"))
    )
})

# ===========================================================================
# Tests migrated from test_dataStructures.R (extractBlockGenotypes)
# ===========================================================================

test_that("extractBlockGenotypes returns SummarizedExperiment", {
    skip_if_not_installed("pgenlibr")
    stem <- test_path("test_data", "test_variants")
    handle <- readGenotypeHandle(stem, format = "plink2")
    n_snps <- nrow(snpInfo(handle))
    skip_if(n_snps == 0, "No SNPs in handle")

    rse <- extractBlockGenotypes(handle, seq_len(min(5L, n_snps)))
    expect_s4_class(rse, "SummarizedExperiment")
    expect_true("dosage" %in% SummarizedExperiment::assayNames(rse))
    dosage <- SummarizedExperiment::assay(rse, "dosage")
    # Bioc convention: variants x samples
    expect_equal(nrow(dosage), min(5L, n_snps))
    expect_equal(ncol(dosage), nSamples(handle))
    # rowRanges should have variant info
    rr <- SummarizedExperiment::rowRanges(rse)
    expect_true("A1" %in% names(S4Vectors::mcols(rr)))
    expect_true("A2" %in% names(S4Vectors::mcols(rr)))
})

test_that("extractBlockGenotypes mean-imputes a missing call on request", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # The BGEN fixture flags sample 2 missing at variant 3.
    handle <- readGenotypeHandle(.bgw_fixture()$path)
    dosageOf <- function(meanImpute) {
        se <- extractBlockGenotypes(handle, 1:5, meanImpute = meanImpute)
        t(as.matrix(SummarizedExperiment::assay(se, "dosage")))
    }
    raw <- dosageOf(FALSE)
    expect_true(is.na(raw[2L, 3L]))
    imputed <- dosageOf(TRUE)
    expect_equal(imputed[2L, 3L], mean(raw[-2L, 3L]))
    expect_equal(imputed[-2L, ], raw[-2L, ])
})

# ===========================================================================
# One-file-per-chromosome (sharded) handle: extraction routing.
# extractBlockGenotypes() routes each region to its per-chromosome payload and
# assembles cross-chromosome requests. Fixtures: chr21 test_variants.* + chr22
# test_variants_chr22.* (same 100 samples; CHR relabelled 22, SNP ids "_c22"),
# available in all four backends.
# ===========================================================================
.shardDose <- function(h, idx) {
    unname(as.matrix(SummarizedExperiment::assay(
        extractBlockGenotypes(h, idx),
        "dosage"
    )))
}

# Per-backend single-file constructor + per-chromosome genoMeta payloads.
# Each payload names its file: the fixtures ship a .bed and a .pgen at one
# stem, and a bare stem with both is read as PLINK2. A PLINK1 set is not a
# sharded handle but stacked panels, tested with the VCF and GDS sets.
.shardBackends <- list(
    plink2 = list(
        pkg = "pgenlibr",
        ref = function(sfx) {
            GenotypeHandle(
                plink2Prefix = file.path(
                    test_data_dir,
                    paste0("test_variants", sfx)
                )
            )
        },
        p21 = file.path(test_data_dir, "test_variants.pgen"),
        p22 = file.path(test_data_dir, "test_variants_chr22.pgen")
    )
)

# Register the per-backend routing tests for one backend spec.
.shardRoutingTests <- function(spec, label) {
    test_that(sprintf("[%s] sharded handle routes per chromosome", label), {
        skip_if_not_installed(spec$pkg)
        ref21 <- spec$ref("")
        ref22 <- spec$ref("_chr22")
        shard <- GenotypeHandle(genoMeta = c("21" = spec$p21, "22" = spec$p22))
        expect_s4_class(shard, "GenotypeHandle")
        expect_equal(genotypeFormat(shard), genotypeFormat(ref21))
        expect_equal(sort(names(chromPaths(shard))), c("chr21", "chr22"))
        n21 <- nrow(snpInfo(ref21))
        n22 <- nrow(snpInfo(ref22))
        expect_equal(nrow(snpInfo(shard)), n21 + n22)
        # A chr21 request routes to the chr21 payload (== single-file chr21).
        expect_equal(.shardDose(shard, 1:5), .shardDose(ref21, 1:5))
        # A chr22 request routes to the chr22 payload, and the global->local index
        # conversion lines up with the single-file chr22 handle. (PLINK2 is
        # positional, so this is the key check for that backend.)
        expect_equal(
            .shardDose(shard, (n21 + 1):(n21 + 5)),
            .shardDose(ref22, 1:5)
        )
    })

    test_that(
        sprintf("[%s] cross-shard request assembles both chromosomes", label),
        {
            skip_if_not_installed(spec$pkg)
            ref21 <- spec$ref("")
            ref22 <- spec$ref("_chr22")
            shard <- GenotypeHandle(
                genoMeta = c("21" = spec$p21, "22" = spec$p22)
            )
            n21 <- nrow(snpInfo(ref21))
            em <- extractBlockGenotypes(shard, c(1L, 2L, n21 + 1L, n21 + 2L))
            dm <- unname(as.matrix(SummarizedExperiment::assay(em, "dosage")))
            expect_equal(nrow(dm), 4L)
            expect_equal(
                as.character(GenomeInfoDb::seqnames(
                    SummarizedExperiment::rowRanges(em)
                )),
                c("chr21", "chr21", "chr22", "chr22")
            )
            # Requested order preserved; each half matches its single-file source.
            expect_equal(dm[1:2, ], .shardDose(ref21, 1:2))
            expect_equal(dm[3:4, ], .shardDose(ref22, 1:2))
        }
    )
}

for (.bk in names(.shardBackends)) {
    .shardRoutingTests(.shardBackends[[.bk]], .bk)
}

test_that("single-shard sharded handle equals the single-file handle", {
    skip_if_not_installed("pgenlibr")
    ref <- GenotypeHandle(
        plink2Prefix = file.path(test_data_dir, "test_variants")
    )
    sh <- GenotypeHandle(
        genoMeta = c("21" = file.path(test_data_dir, "test_variants.pgen"))
    )
    expect_equal(length(chromPaths(sh)), 1L)
    expect_equal(nrow(snpInfo(sh)), nrow(snpInfo(ref)))
    expect_equal(.shardDose(sh, 1:10), .shardDose(ref, 1:10))
})

test_that("readGenotypes(genoMeta=) forwards `format` to the shard reader", {
    skip_if_not_installed("pgenlibr")
    # Regression: the path = "missing" method bound `format` to its own formal
    # and then called GenotypeHandle(...) without it, so an explicit format was
    # silently dropped and every shard fell back to extension probing.
    seen <- new.env(parent = emptyenv())
    seen$fmt <- NA_character_
    local_mocked_bindings(
        .resolveGenotypeShard = function(p, format = NULL) {
            seen$fmt <- if (is.null(format)) "<NULL>" else format
            pecotmr:::.makePlink2Handle(p)
        },
        .package = "pecotmr"
    )
    readGenotypes(
        genoMeta = c("21" = file.path(test_data_dir, "test_variants.pgen")),
        format = "plink2"
    )
    expect_identical(seen$fmt, "plink2")
})

test_that("genoMeta meta-file form matches the named-vector form", {
    skip_if_not_installed("pgenlibr")
    td_abs <- normalizePath(test_data_dir)
    metafile <- tempfile(fileext = ".tsv")
    writeLines(
        c(
            "#chr\tpath",
            paste0("21\t", file.path(td_abs, "test_variants.pgen")),
            paste0("22\t", file.path(td_abs, "test_variants_chr22.pgen"))
        ),
        metafile
    )
    hFile <- GenotypeHandle(genoMeta = metafile)
    hVec <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(td_abs, "test_variants.pgen"),
            "22" = file.path(td_abs, "test_variants_chr22.pgen")
        )
    )
    expect_equal(nrow(snpInfo(hFile)), nrow(snpInfo(hVec)))
    expect_equal(
        sort(names(chromPaths(hFile))),
        sort(names(chromPaths(hVec)))
    )
    expect_equal(.shardDose(hFile, 1:5), .shardDose(hVec, 1:5))
})

# ===========================================================================
# Coverage top-ups for R/genotypeIo.R: dispatch/handle error stops, empty- and
# missing-chromosome block paths, format-specific extractor edge cases, the
# GDS LD path, format detection / stem branches, resolvePlink2Paths
# validation, readStochasticMeta validation, the plink1 refVariantInfo path,
# a multi-range refVariantInfo region, and the non-integer dosages
# loadGenotypeRegion keeps without a sidecar.
# ===========================================================================

# --- readGenotypes dispatch + handle-constructor file-not-found stops --------

test_that("readGenotypes errors on unsupported format", {
    expect_error(
        readGenotypeHandle("/whatever/path", format = "bogusFormat"),
        "Unsupported genotype format: bogusFormat"
    )
})

test_that("a PLINK1 panel errors when its trio is absent", {
    skip_if_not_installed("PlinkMatrix")
    expect_error(
        readGenotypes(
            plink1Prefix = file.path(tempdir(), "missingPlink1Prefix")
        ),
        "Plink file.*File does not exist"
    )
})

# --- .extractBlockSharded: empty block + missing-chromosome stop -------------

test_that("extractBlockGenotypes on a sharded handle handles an empty block", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    shard <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(td, "test_variants.pgen"),
            "22" = file.path(td, "test_variants_chr22.pgen")
        )
    )
    se <- extractBlockGenotypes(shard, integer(0))
    expect_s4_class(se, "SummarizedExperiment")
    dosage <- SummarizedExperiment::assay(se, "dosage")
    expect_equal(nrow(dosage), 0L)
    expect_equal(ncol(dosage), nSamples(shard))
    expect_equal(colnames(dosage), sampleIds(shard))
})

test_that("extractBlockGenotypes answers an empty index without a reader", {
    # A reader need not accept an empty request -- snpStats segfaulted on
    # one when it read PLINK1 -- so none is ever made.
    gh <- .gh_pgenHandle()
    se <- extractBlockGenotypes(gh, snpIdx = integer(0))
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(nrow(se), 0L)
    expect_equal(ncol(se), nSamples(gh))
})

test_that("sharded extraction errors for a chromosome with no payload", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    shard <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(td, "test_variants.pgen"),
            "22" = file.path(td, "test_variants_chr22.pgen")
        )
    )
    # Drop the chr22 payload but keep its variants in @snpInfo, so routing a
    # chr22 request finds no per-chromosome file.
    shard@chromPaths <- chromPaths(shard)["21"]
    idx22 <- which(pecotmr:::canonChrom(snpInfo(shard)$CHR) == "22")[1]
    expect_error(
        extractBlockGenotypes(shard, idx22),
        "no per-chromosome file for chromosome"
    )
})

# --- .extractBlockPlink2: stale-pointer reopen path --------------------------

test_that("plink2 extraction reopens a stale (deserialized) pgen pointer", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    handle <- readGenotypeHandle(
        file.path(td, "test_variants"),
        format = "plink2"
    )
    ref <- SummarizedExperiment::assay(
        extractBlockGenotypes(handle, 1:5),
        "dosage"
    )
    # Round-tripping through saveRDS/readRDS staleness the cached @pgenPtr; the
    # first ReadList errors and the tryCatch reopen path takes over.
    tf <- tempfile(fileext = ".rds")
    on.exit(unlink(tf), add = TRUE)
    saveRDS(handle, tf)
    reloaded <- readRDS(tf)
    got <- SummarizedExperiment::assay(
        extractBlockGenotypes(reloaded, 1:5),
        "dosage"
    )
    expect_equal(got, ref)
})

# --- computeLd on a panel: single-variant and general path ------------------

# --- genotype format detection ----------------------------------------------

test_that("a PLINK prefix is told by the PLINK file beside it", {
    td <- test_path("test_data")
    # protocol_example.genotype has .bed/.bim/.fam but no .pgen.
    expect_equal(
        pecotmr:::.genotypePathFormat(
            file.path(td, "protocol_example.genotype")
        ),
        "plink1"
    )
    # test_variants has both: PLINK2, which can hold dosages, wins.
    expect_equal(
        pecotmr:::.genotypePathFormat(file.path(td, "test_variants")),
        "plink2"
    )
    # A dotted prefix's last part is not taken for an extension.
    dotted <- file.path(withr::local_tempdir(), "ADSP.R4.EUR.chr21")
    file.create(str_c(dotted, c(".pgen", ".pvar", ".psam")))
    expect_equal(pecotmr:::.genotypePathFormat(dotted), "plink2")
    # A GDS file is named by its own path, never by a stem beside it.
    stem <- tempfile("gioGdsStem_")
    file.create(paste0(stem, ".gds"))
    on.exit(unlink(paste0(stem, ".gds")), add = TRUE)
    expect_null(pecotmr:::.genotypePathFormat(stem))
})

test_that("a path naming no genotype format is refused by its caller", {
    for (path in file.path(tempdir(), c("thing.qzx", "noSuchStemXyz"))) {
        expect_error(
            pecotmr:::.requireGenotypeFormat(path, "caller"),
            "caller: .* is neither a genotype file nor a PLINK prefix"
        )
    }
    expect_error(
        readGenotypes(file.path(tempdir(), "thing.qzx")),
        "readGenotypes: .* is neither a genotype file nor a PLINK prefix"
    )
})

# --- .genotypeStem strips the extension naming a format -----------------------

test_that(".genotypeStem strips the extension naming a format", {
    expect_equal(
        pecotmr:::.genotypeStem(c(
            "a/x.vcf.gz",
            "x.VCF.BGZ",
            "x.bcf",
            "x.gds",
            "x.bgen",
            "x.bed",
            "x.bim",
            "x.pgen",
            "x.pvar.zst",
            "x.Psam",
            "plink/stem",
            "dotted.prefix.chr21"
        )),
        c(
            "a/x",
            "x",
            "x",
            "x",
            "x",
            "x",
            "x",
            "x",
            "x",
            "x",
            "plink/stem",
            "dotted.prefix.chr21"
        )
    )
})

# --- resolvePlink2Paths: .pvar.zst fallback + missing pvar/psam stops --------

.gioMakePlink2Stub <- function(which = c("pgen", "pvar", "pvarZst", "psam")) {
    dir <- tempfile("gioP2_")
    dir.create(dir)
    prefix <- file.path(dir, "g")
    if ("pgen" %in% which) {
        file.create(paste0(prefix, ".pgen"))
    }
    if ("pvar" %in% which) {
        file.create(paste0(prefix, ".pvar"))
    }
    if ("pvarZst" %in% which) {
        file.create(paste0(prefix, ".pvar.zst"))
    }
    if ("psam" %in% which) {
        file.create(paste0(prefix, ".psam"))
    }
    prefix
}

test_that("resolvePlink2Paths falls back to .pvar.zst", {
    prefix <- .gioMakePlink2Stub(c("pgen", "pvarZst", "psam"))
    paths <- pecotmr:::resolvePlink2Paths(prefix)
    expect_equal(paths$pvar, paste0(prefix, ".pvar.zst"))
    expect_equal(paths$pgen, paste0(prefix, ".pgen"))
    expect_equal(paths$psam, paste0(prefix, ".psam"))
})

test_that("resolvePlink2Paths errors when no .pvar/.pvar.zst exists", {
    prefix <- .gioMakePlink2Stub(c("pgen", "psam"))
    expect_error(
        pecotmr:::resolvePlink2Paths(prefix),
        "PLINK2 .pvar\\[.zst\\] file not found"
    )
})

test_that("resolvePlink2Paths errors when .psam is missing", {
    prefix <- .gioMakePlink2Stub(c("pgen", "pvar"))
    expect_error(
        pecotmr:::resolvePlink2Paths(prefix),
        "PLINK2 .psam file.*File does not exist"
    )
})

# --- readStochasticMeta: generic format missing required columns -------------

test_that("readStochasticMeta errors when generic file lacks required columns", {
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(c("foo\tbar", "1\t2"), tmp)
    expect_error(
        pecotmr:::readStochasticMeta(tmp),
        "must contain columns"
    )
})

# --- refVariantInfo: plink1 .bed path + multi-row region filter -----------

test_that("refVariantInfo reads variant info from a plink1 source", {
    td <- test_path("test_data")
    meta_file <- file.path(td, "ld_meta_refinfo_plink1_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("22", "0", "0", "protocol_example.genotype", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    info <- refVariantInfo(meta_file) # no region -> all variants
    expect_true(is.data.frame(info))
    expect_true(all(c("chrom", "id", "pos", "A2", "A1") %in% names(info)))
    expect_gt(nrow(info), 0L)
})

test_that("refVariantInfo applies a multi-range region filter", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    meta_file <- file.path(td, "ld_meta_refinfo_multirow_tmp.tsv")
    on.exit(unlink(meta_file), add = TRUE)
    writeLines(paste("chrom", "start", "end", "path", sep = "\t"), meta_file)
    cat(
        paste("21", "0", "0", "test_variants", sep = "\t"),
        "\n",
        file = meta_file,
        append = TRUE
    )
    # A two-range window (both chr21) drives the multi-range filter; together
    # the ranges cover the full block.
    twoRanges <- GenomicRanges::GRanges(
        "21",
        IRanges::IRanges(c(17513228L, 17550001L), c(17550000L, 17592874L))
    )
    info <- refVariantInfo(meta_file, which = twoRanges)
    expect_true(is.data.frame(info))
    expect_equal(nrow(info), 349L)
    expect_true(all(info$pos >= 17513228 & info$pos <= 17592874))
})

# --- loadGenotypeRegion: non-integer dosages without a sidecar --------------

test_that("loadGenotypeRegion keeps imputed dosages without a sidecar", {
    skip_if_not_installed("pgenlibr")
    td <- test_path("test_data")
    # dummy_data is a plink2 fixture with no .afreq / .stochastic_meta sidecar.
    # A dosage need not be a whole number, so fractional dosages without a
    # sidecar are an imputed file's, not a sketch's: returned untouched, and
    # without a warning about restoring a scale they never had.
    testthat::local_mocked_bindings(
        .extractBlockPlink2 = function(handle, snpIdx) {
            matrix(0.5, nrow = nSamples(handle), ncol = length(snpIdx))
        },
        .package = "pecotmr"
    )
    expect_no_warning(X <- loadGenotypeRegion(file.path(td, "dummy_data")))
    expect_true(all(X == 0.5))
})

# ---- fileIdx + snpInfo subsetting (memory-safe LD-sketch trimming) -----------

test_that(".withFileIdx + readers attach a sequential fileIdx to snpInfo", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    expect_true("fileIdx" %in% names(si))
    expect_identical(si$fileIdx, seq_len(nrow(si)))
})

test_that(".subsetGenotypeHandle keeps PLINK2 reads correct for kept variants", {
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    si <- snpInfo(h)
    want <- c(2L, 4L, 5L, 9L)
    dosFull <- pecotmr:::.dosageMatrix(h, want, meanImpute = TRUE)
    # Keep a reordered superset so the wanted variants move to NEW row positions;
    # fileIdx must rescue the (now-shifted) positional PLINK2 read.
    keep <- c(9L, 5L, 4L, 2L, 7L, 1L)
    hSub <- pecotmr:::.subsetGenotypeHandle(h, keep)
    expect_equal(nrow(snpInfo(hSub)), length(keep))
    expect_identical(snpInfo(hSub)$fileIdx, keep) # original file order retained
    wantSub <- match(want, keep)
    dosSub <- pecotmr:::.dosageMatrix(hSub, wantSub, meanImpute = TRUE)
    expect_equal(unname(dosFull), unname(dosSub))
})

test_that(".subsetGenotypeHandle is NULL-safe and a no-op when nothing dropped", {
    expect_null(pecotmr:::.subsetGenotypeHandle(NULL, TRUE))
    h <- readGenotypeHandle(
        test_path("test_data/test_variants"),
        format = "plink2"
    )
    expect_identical(
        pecotmr:::.subsetGenotypeHandle(h, seq_len(nrow(snpInfo(h)))),
        h
    )
})

# ---------------------------------------------------------------------------
# extractBlockGenotypes: an unsorted request must not mislabel the block.
#
# snpStats (select.snps=) and SNPRelate (snp.id=) return variants in FILE
# order whatever order was asked for, while the block is labelled from
# snpInfo[snpIdx] -- the REQUESTED order. Both are the right length, so the
# mismatch was silent: names described one variant and dosages another.
# ---------------------------------------------------------------------------

# @noRd
.gio_orderCheck <- function(fmt, path) {
    handle <- readGenotypeHandle(path, format = fmt)
    perm <- c(4L, 1L, 6L, 2L, 5L, 3L)
    si <- snpInfo(handle)
    seq6 <- extractBlockGenotypes(handle, 1:6, meanImpute = FALSE)
    seP <- extractBlockGenotypes(handle, perm, meanImpute = FALSE)
    a6 <- as.matrix(SummarizedExperiment::assay(seq6, "dosage"))
    aP <- as.matrix(SummarizedExperiment::assay(seP, "dosage"))
    # Names follow the request ...
    expect_identical(rownames(seP), si$SNP[perm])
    # ... and so do the dosages behind them.
    expect_equal(unname(aP), unname(a6[perm, , drop = FALSE]))
}

test_that("extractBlockGenotypes honours request order (plink2)", {
    skip_if_not_installed("pgenlibr")
    .gio_orderCheck("plink2", test_path("test_data", "test_variants"))
})

test_that("a PLINK1, GDS or delayed PLINK2 read honours request order", {
    # snpStats returned file order when it read PLINK1; the panel's seed
    # reads each block in file order and lays it out as asked. DelayedArray
    # hands the index straight to extract_array(), so a handle's delayed
    # array inherits whatever the reader does.
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    skip_if_not_installed("pgenlibr")
    stem <- test_path("test_data", "test_variants")
    arrays <- list(
        plink1 = SummarizedExperiment::assay(
            readGenotypes(plink1Prefix = stem),
            "dosage"
        ),
        gds = SummarizedExperiment::assay(
            readGenotypes(str_c(stem, ".gds")),
            "dosage"
        ),
        plink2 = genotypeDelayedArray(
            readGenotypeHandle(stem, format = "plink2")
        )
    )
    perm <- c(4L, 1L, 6L, 2L, 5L, 3L)
    for (format in names(arrays)) {
        a <- arrays[[format]]
        expect_equal(
            unname(as.matrix(a[perm, 1:4])),
            unname(as.matrix(a[1:6, 1:4])[perm, , drop = FALSE]),
            info = format
        )
    }
})

test_that(".restoreRequestedOrder refuses a short block", {
    # A backend returning fewer variants than asked for cannot be labelled
    # from snpInfo[snpIdx]; that must be loud, not silently truncated.
    expect_error(
        pecotmr:::.restoreRequestedOrder(matrix(0, 3, 2), 1:3),
        "cannot be labelled"
    )
})

test_that("ordering uses FILE position, not the snpInfo row number", {
    # A row-subset handle's snpInfo rows are no longer file positions, so
    # sorting the request by row number would ask the backend for the wrong
    # ascending order. fileIdx carries the real mapping.
    skip_if_not_installed("pgenlibr")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    keep <- c(9L, 2L, 7L, 4L)
    sub <- pecotmr:::.subsetGenotypeHandle(handle, keep)
    expect_identical(snpInfo(sub)$fileIdx, keep)
    got <- as.matrix(SummarizedExperiment::assay(
        extractBlockGenotypes(sub, seq_along(keep), meanImpute = FALSE),
        "dosage"
    ))
    want <- as.matrix(SummarizedExperiment::assay(
        extractBlockGenotypes(handle, keep, meanImpute = FALSE),
        "dosage"
    ))
    expect_equal(unname(got), unname(want))
})

# ---------------------------------------------------------------------------
# readGenotypes(): the exported entry point returns a genotype PANEL.
#
# The handle is the seed layer beneath it and is no longer public. These pin
# the panel contract; the tests above exercise the reader machinery through
# readGenotypeHandle() and assert on the handle itself.
# ---------------------------------------------------------------------------

test_that("readGenotypes returns a panel, not a handle", {
    skip_if_not_installed("pgenlibr")
    panel <- readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    expect_s4_class(panel, "RangedSummarizedExperiment")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    # variants x samples, the Bioconductor orientation.
    expect_equal(nrow(panel), nrow(snpInfo(handle)))
    expect_equal(
        ncol(panel),
        length(sampleIds(
            pecotmr:::.ghSeedHandle(
                DelayedArray::seed(SummarizedExperiment::assay(panel, "dosage"))
            )
        ))
    )
})

test_that("the panel's assay is delayed and its seed is the handle", {
    skip_if_not_installed("pgenlibr")
    panel <- readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    a <- SummarizedExperiment::assay(panel, "dosage")
    expect_s4_class(a, "DelayedMatrix")
    handle <- readGenotypeHandle(
        test_path("test_data", "test_variants"),
        format = "plink2"
    )
    # Equal, not identical: each read opens its own pgen pointer.
    expect_equal(pecotmr:::.ldRefHandle(panel), handle)
})

test_that("a panel answers panel questions without the handle", {
    skip_if_not_installed("PlinkMatrix")
    panel <- readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink1"
    )
    # The RSE spellings that replace the retired handle accessors.
    expect_type(rownames(panel), "character")
    expect_type(colnames(panel), "character")
    gr <- SummarizedExperiment::rowRanges(panel)
    expect_s4_class(gr, "GRanges")
    cn <- colnames(S4Vectors::mcols(gr))
    expect_true(all(is_in(c("SNP", "A1", "A2"), cn)))
})

test_that("readGenotypes names a panel by keyword when path is missing", {
    # GenotypeHandle is internal now, so its source vocabulary (prefixes,
    # explicit triplets, genoMeta) has to stay reachable from the public
    # entry point or it is gone from the API entirely.
    skip_if_not_installed("PlinkMatrix")
    stem <- test_path("test_data", "test_variants")
    byPath <- readGenotypes(stem, format = "plink1")
    byPrefix <- readGenotypes(plink1Prefix = stem)
    expect_equal(dim(byPrefix), dim(byPath))
    expect_equal(rownames(byPrefix), rownames(byPath))
    expect_equal(colnames(byPrefix), colnames(byPath))
})

test_that("readGenotypes reaches the explicit-triplet source too", {
    skip_if_not_installed("PlinkMatrix")
    stem <- test_path("test_data", "test_variants")
    panel <- readGenotypes(
        bed = str_c(stem, ".bed"),
        bim = str_c(stem, ".bim"),
        fam = str_c(stem, ".fam")
    )
    expect_s4_class(panel, "RangedSummarizedExperiment")
    expect_equal(dim(panel), dim(readGenotypes(stem, format = "plink1")))
})

test_that("readGenotypes reaches the explicit PLINK2 triplet source", {
    skip_if_not_installed("pgenlibr")
    stem <- test_path("test_data", "test_variants")
    panel <- readGenotypes(
        pgen = str_c(stem, ".pgen"),
        pvar = str_c(stem, ".pvar"),
        psam = str_c(stem, ".psam")
    )
    byPrefix <- readGenotypes(plink2Prefix = stem)
    expect_s4_class(panel, "RangedSummarizedExperiment")
    expect_equal(rownames(panel), rownames(byPrefix))
    expect_equal(colnames(panel), colnames(byPrefix))
    expect_equal(
        as.matrix(SummarizedExperiment::assay(panel, "dosage")[1:5, ]),
        as.matrix(SummarizedExperiment::assay(byPrefix, "dosage")[1:5, ])
    )
})

# ===========================================================================
# Format detection and handle/covariate helpers
#
# The extension->format switch had only two of its fourteen arms exercised,
# and the NULL-safe / legacy-handle / covariate-alignment branches none.
# ===========================================================================

test_that("every file extension maps to its genotype format", {
    f <- function(path) pecotmr:::.genotypeFileFormat(path) %||% "none"
    paths <- c(
        vcf = "x.vcf",
        vcf = "x.vcf.gz",
        vcf = "x.VCF.BGZ",
        vcf = "x.bcf",
        gds = "x.gds",
        bgen = "x.bgen",
        plink1 = "x.bed",
        plink1 = "x.bim",
        plink1 = "x.fam",
        plink2 = "x.pgen",
        plink2 = "x.pvar",
        plink2 = "x.pvar.zst",
        plink2 = "x.psam",
        # An unrecognised extension is NULL, not an error: the caller decides.
        none = "x.zzz",
        none = "x.annot.gz"
    )
    expect_equal(purrr::map_chr(unname(paths), f), names(paths))
})

test_that("the unsupported-format signaller names its context when given one", {
    expect_error(
        pecotmr:::.abortUnsupportedFormat("weird", "myContext"),
        "Unsupported format in myContext: weird"
    )
    expect_error(
        pecotmr:::.abortUnsupportedFormat("weird"),
        "Unsupported genotype format: weird"
    )
})

test_that("the handle helpers are NULL-safe", {
    expect_null(pecotmr:::.emptyGenotypeHandle(NULL))
    expect_null(pecotmr:::.subsetGenotypeHandle(NULL, 1L))
})

test_that("a legacy handle without fileIdx is not subset", {
    # Without fileIdx there is no way to keep on-disk positions correct after
    # dropping rows, so the handle is returned whole rather than silently
    # mis-indexed.
    h <- new(
        "GenotypeHandle",
        path = "/tmp/x",
        format = "plink2",
        snpInfo = data.frame(
            SNP = c("a", "b", "c"),
            CHR = "chr1",
            BP = 1:3,
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 2L,
        sampleIds = c("s1", "s2"),
        pgenPtr = NULL
    )
    sub <- pecotmr:::.subsetGenotypeHandle(h, 1L)
    expect_equal(nrow(pecotmr:::snpInfo(sub)), 3L)
})

test_that("covariate colData aligns to the panel's samples", {
    f <- pecotmr:::.genotypeColData
    ids <- c("s1", "s2")
    # No covariates at all: an empty DataFrame that still names the samples.
    empty <- f(matrix(numeric(0), 0L, 0L), ids)
    expect_equal(rownames(empty), ids)
    expect_equal(ncol(empty), 0L)
    # Unnamed rows of the right length are taken positionally.
    cd <- f(matrix(1:2, 2L, 1L, dimnames = list(NULL, "pc1")), ids)
    expect_equal(rownames(cd), ids)
    expect_equal(as.numeric(cd$pc1), c(1, 2))
    # Unnamed rows of the WRONG length are an error, not a silent recycle.
    expect_error(
        f(matrix(1, 3L, 1L), ids),
        "has 3 rows but the panel has 2 samples"
    )
})

test_that("a panel with no covariates gets an empty colData", {
    prefix <- sub(
        "\\.bed$",
        "",
        system.file("extdata", "toy_ref.bed", package = "pecotmr")
    )
    panel <- readGenotypes(plink1Prefix = prefix)
    out <- pecotmr:::.genotypeExperimentCovariates(panel, NULL)
    cd <- SummarizedExperiment::colData(out)
    expect_equal(ncol(cd), 0L)
    # The samples are still named, so the assay keeps a column for each.
    expect_equal(rownames(cd), colnames(panel))
})

test_that("block extraction rejects a format it has no reader for", {
    # The constructor rejects an unknown format outright, so this guard can
    # only be reached by a direct slot set -- which is exactly what it exists
    # to catch.
    h <- new(
        "GenotypeHandle",
        path = "/tmp/x",
        format = "plink2",
        snpInfo = data.frame(
            SNP = "a",
            CHR = "chr1",
            BP = 1L,
            A1 = "A",
            A2 = "G",
            fileIdx = 1L,
            stringsAsFactors = FALSE
        ),
        nSamples = 1L,
        sampleIds = "s1",
        pgenPtr = NULL
    )
    expect_error(
        new(
            "GenotypeHandle",
            path = "/tmp/x",
            format = "zzz",
            snpInfo = snpInfo(h),
            nSamples = 1L,
            sampleIds = "s1",
            pgenPtr = NULL
        ),
        "'format' must be one of"
    )
    h@format <- "zzz"
    expect_error(
        pecotmr:::.extractBlockByFormat(h, 1L),
        "Unsupported format in extractBlockGenotypes: zzz"
    )
})

test_that("stochastic inversion is a no-op when no id matches the metadata", {
    # The metadata exists and parses, but names other variants entirely, so
    # there is nothing to invert and the result passes through untouched.
    meta <- withr::local_tempfile(fileext = ".stochastic_meta.tsv")
    readr::write_tsv(
        data.frame(
            id = c("other1", "other2"),
            u_min = c(0, 0),
            u_max = c(1, 1)
        ),
        meta
    )
    X <- matrix(
        1:4,
        2L,
        2L,
        dimnames = list(c("s1", "s2"), c("v1", "v2"))
    )
    res <- list(X = X)
    out <- pecotmr:::.loadGenoInvertStochastic(
        res,
        genotype = "/tmp/does_not_matter",
        stochasticMetaPath = meta,
        stochasticMetaFormat = "generic"
    )
    expect_identical(out$X, X)
})

test_that("readGenotypes names its keyword sources, not `...`", {
    # The sources used to ride `...` into GenotypeHandle(), so a misspelled
    # one was simply unset and the call failed later with "no source given".
    expect_error(readGenotypes(plink1Prefx = "x"), "unused argument")
    # Every GenotypeHandle source is reachable by name: passing one gets past
    # argument matching and fails (if at all) on the source itself, never on
    # "unused argument".
    sources <- setdiff(
        names(formals(pecotmr:::GenotypeHandle)),
        c("path", "format", "vcfArgs")
    )
    for (nm in sources) {
        err <- tryCatch(
            {
                exec(readGenotypes, !!!set_names(list("nope"), nm))
                NULL
            },
            error = function(e) conditionMessage(e)
        )
        expect_false(
            isTRUE(str_detect(err %||% "", "unused argument")),
            label = nm
        )
    }
})

# ===========================================================================
# The BGEN reader (and src/bgen.cpp). Synthetic files come from the writer in
# helper-bgen.R; complex.bgen and haplotypes.bgen are bgen's own examples
# (test_data/bgen, Boost Software License 1.0). A BGEN file is read through
# its .bgi index, an SQLite database, so every test reading one needs DBI and
# RSQLite.
# ===========================================================================

.bgenPanelDosage <- function(panel) {
    t(as.matrix(SummarizedExperiment::assay(panel, "dosage")))
}

test_that("a BGEN file reads as the dosages its probabilities give", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # Every compression and a bit depth that does not fill bytes: the
    # dosage is the expected count of the second allele, NA where a sample
    # is flagged missing.
    for (compression in 0:2) {
        for (bits in c(3L, 8L, 16L)) {
            fx <- .bgw_fixture(bits = bits, compression = compression)
            panel <- readGenotypes(fx$path)
            label <- str_c("compression ", compression, ", bits ", bits)
            expect_equal(dim(panel), c(5L, 6L), label = label)
            expect_equal(colnames(panel), fx$ids, label = label)
            expect_equal(
                unname(.bgenPanelDosage(panel)),
                fx$dosage,
                tolerance = 1e-12,
                label = label
            )
        }
    }
})

test_that("a BGEN file's variants read as pecotmr's convention", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # A2 the first allele, A1 the second; the record's rsid as its id when
    # it has one, else chr:pos:A2:A1; a zero-padded chromosome is canonical.
    fx <- .bgw_fixture()
    handle <- pecotmr:::.readGenotypeHandle(fx$path)
    si <- snpInfo(handle)
    expect_equal(genotypeFormat(handle), "bgen")
    expect_equal(si$A2, rep("A", 5L))
    expect_equal(si$A1, rep("G", 5L))
    expect_equal(si$CHR, rep("chr1", 5L))
    expect_equal(si$SNP, c("chr1:100:A:G", sprintf("rs%d", 1:4)))
    expect_equal(si$BP, 100L * (1:5))
})

test_that("a BGEN file without sample ids reads them from its .sample", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    fx <- .bgw_fixture(storeIds = FALSE)
    panel <- readGenotypes(fx$path)
    expect_equal(colnames(panel), fx$ids)
    expect_equal(unname(.bgenPanelDosage(panel)), fx$dosage, tolerance = 1e-12)
    file.remove(str_c(str_remove(fx$path, "\\.bgen$"), ".sample"))
    expect_error(readGenotypes(fx$path), "stores no sample ids")
})

test_that("a BGEN 1.1 (layout 1) file reads as its probabilities give", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # Diploid and biallelic, three 16-bit probabilities per sample out of
    # 32768, zlib-compressed or not; a sample stored as all zeros is missing.
    for (compression in 0:1) {
        fx <- .bgw_fixture(compression = compression, layout = 1L)
        panel <- readGenotypes(fx$path)
        label <- str_c("layout 1, compression ", compression)
        expect_equal(colnames(panel), fx$ids, label = label)
        expect_equal(
            unname(.bgenPanelDosage(panel)),
            fx$dosage,
            tolerance = 1e-12,
            label = label
        )
    }
})

test_that("a BGEN file needs its .bgi index and a layout the spec defines", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    fx <- .bgw_fixture()
    file.remove(str_c(fx$path, ".bgi"))
    expect_error(readGenotypes(fx$path), "is read through its .bgi index")
    # The header is refused before any record is read.
    header <- function(layout, compression) {
        path <- withr::local_tempfile(
            fileext = ".bgen",
            .local_envir = parent.frame()
        )
        .bgw_write(
            path,
            fx$variants[0L, ],
            list(),
            matrix(FALSE, 6L, 0L),
            fx$ids,
            compression = compression,
            layout = layout
        )
        path
    }
    expect_error(readGenotypes(header(3L, 1L)), "BGEN layout 3")
    expect_error(
        readGenotypes(header(1L, 2L)),
        "compressed with zstd, which the spec does not allow"
    )
})

test_that("a corrupt BGEN header or sample block is refused", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # The fixture's bytes: a 4-byte offset, the 20-byte header (its length,
    # the variant and sample counts, the magic number, the flags), then the
    # sample block (its length, the sample count again, the ids).
    fx <- .bgw_fixture()
    bytes <- readBin(fx$path, "raw", file.size(fx$path))
    dir <- withr::local_tempdir()
    rewritten <- function(name, bytes) {
        path <- file.path(dir, str_c(name, ".bgen"))
        writeBin(bytes, path)
        file.copy(str_c(fx$path, ".bgi"), str_c(path, ".bgi"))
        path
    }
    expect_error(
        readGenotypes(rewritten("short", bytes[1:10])),
        "too short to be a BGEN file"
    )
    expect_error(
        readGenotypes(rewritten(
            "magic",
            replace(bytes, 17:20, charToRaw("bgex"))
        )),
        "no 'bgen' magic number"
    )
    # The first flag byte's low two bits are the compression; 3 is undefined.
    expect_error(
        readGenotypes(rewritten("flags", replace(bytes, 21L, as.raw(3L + 8L)))),
        "unknown BGEN compression \\(3\\)"
    )
    expect_error(
        readGenotypes(rewritten("count", replace(bytes, 29L, as.raw(7L)))),
        "sample block's sample count disagrees with its header"
    )
})

test_that("a BGEN file's sample ids and records must match its header", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # A .sample file listing fewer samples than the header declares.
    fx <- .bgw_fixture(storeIds = FALSE)
    writeLines(
        c(
            "ID_1 ID_2 missing",
            "0 0 0",
            str_c(fx$ids[-1], fx$ids[-1], "0", sep = " ")
        ),
        str_c(str_remove(fx$path, "\\.bgen$"), ".sample")
    )
    expect_error(
        readGenotypes(fx$path),
        "declares 6 samples but 5 sample ids were found"
    )
    # A layout-1 record holding five samples' probabilities under a header
    # that declares six: refused when its dosages are read.
    path <- file.path(withr::local_tempdir(), "short.bgen")
    .bgw_write(
        path,
        fx$variants[1L, ],
        list(matrix(c(10000L, 20000L), 5L, 2L, byrow = TRUE)),
        matrix(FALSE, 5L, 1L),
        fx$ids,
        layout = 1L
    )
    panel <- readGenotypes(path)
    expect_error(
        as.matrix(SummarizedExperiment::assay(panel, "dosage")),
        "layout-1 BGEN record holds the wrong amount of data"
    )
})

test_that("bgen's complex example reads as its VCF gives it", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # Haploid and diploid samples, phased and unphased records, multi-
    # allelic records split one variant per non-first allele, and records
    # of ploidy above 2.
    path <- test_path("test_data", "bgen", "complex.bgen")
    expect_warning(
        handle <- pecotmr:::.readGenotypeHandle(path),
        "carry more than two alleles"
    )
    si <- snpInfo(handle)
    m4 <- which(si$BP == 4L)
    expect_equal(si$A1[m4], c("G", "T"))
    dosageOf <- function(rows) {
        pecotmr:::.extractBlockBgen(handle, rows)
    }
    # V1, V2: unphased, haploid sample 1; V3: phased; M4: G and T.
    want <- cbind(
        V1 = c(0, 0, 0, 1),
        V2 = c(0, 0, 0, 1),
        V3 = c(0, 2, 0, 1),
        M4G = c(0, 0, 1, 0),
        M4T = c(0, 0, 0, 1)
    )
    rows <- c(which(si$BP %in% 1:3), m4)
    expect_equal(unname(dosageOf(rows)), unname(want))
    # M5: phased, ploidies 1, 3, 3 and 2 (0, 0|0|0, 0|0|1, 0|1); M10:
    # unphased tetraploid (0/0/0/0 to 0/1/1/1). A dosage runs from 0 to the
    # sample's ploidy.
    expect_equal(
        unname(dosageOf(which(si$BP %in% c(5L, 10L)))),
        cbind(c(0, 0, 1, 1), c(0, 1, 2, 3))
    )
})

test_that("bgen's haplotype example reads as its .haps gives it", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    path <- test_path("test_data", "bgen", "haplotypes.bgen")
    panel <- readGenotypes(path)
    # Columns: chromosome, id, rsid, position, the two alleles, then two
    # haplotypes per sample.
    haps <- read.table(test_path("test_data", "bgen", "haplotypes.haps"))
    alleles <- as.matrix(haps[, -(1:6)])
    want <- alleles[, c(TRUE, FALSE)] + alleles[, c(FALSE, TRUE)]
    expect_equal(
        unname(as.matrix(SummarizedExperiment::assay(panel, "dosage"))),
        unname(want)
    )
})

test_that("a BGEN file reads as bigsnpr reads it", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # An independent reader of the UK Biobank layout (8-bit zlib, biallelic),
    # which rounds its dosages to two decimals.
    skip_if_not_installed("bigsnpr")
    fx <- .bgw_fixture(n = 40L, nVar = 12L, bits = 8L, compression = 1L)
    ours <- .bgenPanelDosage(readGenotypes(fx$path))
    ids <- str_c(
        fx$variants$chrom,
        fx$variants$pos,
        fx$variants$a1,
        fx$variants$a2,
        sep = "_"
    )
    rds <- bigsnpr::snp_readBGEN(
        fx$path,
        backingfile = withr::local_tempfile(),
        list_snp_id = list(ids),
        read_as = "dosage",
        ncores = 1
    )
    theirs <- bigsnpr::snp_attach(rds)$genotypes[]
    expect_equal(is.na(theirs), unname(is.na(ours)))
    present <- !is.na(ours)
    expect_lt(max(abs(unname(ours)[present] - theirs[present])), 0.005 + 1e-9)
})

test_that("a BGEN file serves every genotype route", {
    skip_if_not_installed("DBI")
    skip_if_not_installed("RSQLite")
    # The format is one of the GenotypeHandle's, so every route that reads
    # genotypes reads it: a file, a per-chromosome manifest, an LD-meta
    # table, the region loader, LD, and a mixture of cohorts.
    fx <- .bgw_fixture(n = 30L, nVar = 8L, bits = 16L)
    want <- fx$dosage
    panel <- readGenotypes(fx$path)
    sharded <- readGenotypes(genoMeta = c(`1` = fx$path))
    expect_equal(unname(.bgenPanelDosage(sharded)), want, tolerance = 1e-12)
    meta <- withr::local_tempfile(fileext = ".tsv")
    writeLines(
        c("#chr\tstart\tend\tpath", str_c("chr1\t1\t10000\t", fx$path)),
        meta
    )
    viaMeta <- readGenotypes(
        ldMeta = meta,
        block = GenomicRanges::GRanges("chr1:1-10000")
    )
    expect_equal(unname(.bgenPanelDosage(viaMeta)), want, tolerance = 1e-12)
    loaded <- loadGenotypeRegion(
        fx$path,
        which = GenomicRanges::GRanges("chr1:1-10000")
    )
    expect_equal(unname(loaded), want, tolerance = 1e-12)
    expect_equal(
        refVariantInfo(fx$path)$pos,
        fx$variants$pos
    )
    # LD is the samples' correlation, the missing call imputed with the
    # variant's mean as everywhere else.
    imputed <- want
    imputed[is.na(imputed)] <- colMeans(want, na.rm = TRUE)[
        col(want)[is.na(want)]
    ]
    expect_equal(unname(computeLd(panel)), cor(imputed), tolerance = 1e-10)
    mix <- GenotypeMixture(
        list(a = panel[, 1:15], b = panel[, 16:30]),
        mixtureWeights = c(0.5, 0.5)
    )
    expect_equal(pecotmr:::.ldRefNVariants(mix), 8L)
})

test_that("the BGEN dosage kernel checks what it is given", {
    expect_error(
        pecotmr:::bgenDosageCpp(as.raw(1:5), 1L, 2L),
        "truncated"
    )
    # A record claiming more genotypes than its data could hold is refused
    # before they are enumerated: ploidy 63 with 65535 alleles, one byte.
    corrupt <- as.raw(c(1, 0, 0, 0, 255, 255, 63, 63, 63, 0, 8, 0))
    expect_error(pecotmr:::bgenDosageCpp(corrupt, 1L, 1L), "truncated")
})

test_that("the BGEN dosage kernel reads any ploidy and allele count", {
    # One unphased triploid sample with three alleles: ten genotypes in colex
    # order -- (3,0,0) (2,1,0) (1,2,0) (0,3,0) (2,0,1) (1,1,1) (0,2,1)
    # (1,0,2) (0,1,2) (0,0,3) -- the last one's probability implied.
    triploid <- function(genotype) {
        stored <- integer(9L)
        if (genotype <= 9L) {
            stored[[genotype]] <- 255L
        }
        as.raw(c(1, 0, 0, 0, 3, 0, 3, 3, 3, 0, 8, stored))
    }
    # All the mass on (0,2,1): two copies of allele 2, one of allele 3.
    expect_equal(
        pecotmr:::bgenDosageCpp(triploid(7L), 1L, c(2L, 3L)),
        matrix(c(2, 1), 1L)
    )
    # On the implied last genotype, (0,0,3).
    expect_equal(
        pecotmr:::bgenDosageCpp(triploid(10L), 1L, c(1L, 3L)),
        matrix(c(0, 3), 1L)
    )
})

# ===========================================================================
# The VCF and GDS panel builders: what a VCF or GDS panel's variants are
# called, what is refused before anything is read, and how a per-chromosome
# set stacks. The seed's reads are tested in test_GenotypeDosageSeed.R.
# ===========================================================================

samples3 <- c("s1", "s2", "s3")

# A GDS file built node by node with gdsfmt, of root format `fileFormat`
# (NULL records none), by `build(gds)`.
.writeTestGds <- function(fileFormat, build) {
    f <- tempfile(fileext = ".gds")
    gds <- gdsfmt::createfn.gds(f)
    if (!is.null(fileFormat)) {
        gdsfmt::put.attr.gdsn(gds$root, "FileFormat", fileFormat)
    }
    build(gds)
    gdsfmt::closefn.gds(gds)
    f
}

# A minimal SeqArray file of one biallelic record and two samples, all REF
# calls: its genotype calls stored as `storage`, with `nExtra` calls stored
# beyond its ploidy, and no annotation/id.
.writeFakeSeqArray <- function(storage = "bit2", nExtra = 0L) {
    .writeTestGds("SEQ_ARRAY", function(gds) {
        gdsfmt::add.gdsn(gds, "sample.id", c("a", "b"))
        gdsfmt::add.gdsn(gds, "variant.id", 1L)
        gdsfmt::add.gdsn(gds, "position", 100L)
        gdsfmt::add.gdsn(gds, "chromosome", "1")
        gdsfmt::add.gdsn(gds, "allele", "A,G")
        geno <- gdsfmt::addfolder.gdsn(gds, "genotype")
        gdsfmt::add.gdsn(
            geno,
            "data",
            array(0L, c(2L, 2L, 1L)),
            storage = storage
        )
        gdsfmt::add.gdsn(geno, "@data", 1L, storage = "uint8", visible = FALSE)
        gdsfmt::add.gdsn(
            geno,
            "extra.index",
            matrix(1L, 3L, nExtra),
            storage = "int32"
        )
    })
}

# ---- VCF -------------------------------------------------------------------

test_that("a VCF panel names its variants and records where it came from", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    expect_warning(
        panel <- readGenotypes(vcf),
        "1 VCF record\\(s\\) of .* carry more than one ALT allele"
    )
    rr <- SummarizedExperiment::rowRanges(panel)
    # A record keeps its own ID when it has one and a single ALT allele;
    # every other row is named by its alleles.
    expect_equal(
        rr$SNP,
        c(
            "chr22:100:A:G",
            "rs2",
            "rs3",
            "rs4",
            "rs4",
            "chr22:400:A:C",
            "chr22:400:A:G"
        )
    )
    expect_equal(rr$A1, c("G", "A", "T", "T", "G", "C", "G"))
    expect_equal(rr$A2, c("A", "ACGTA", "G", "C", "C", "A", "A"))
    expect_equal(as.character(GenomicRanges::seqnames(rr)), rep("chr22", 7L))
    expect_equal(colnames(panel), samples3)
    expect_equal(
        S4Vectors::metadata(panel)$genotypeSource,
        list(
            format = "vcf",
            path = normalizePath(vcf),
            chromPaths = character(0)
        )
    )
})

test_that("a VCF record with no ALT allele holds no variant", {
    skip_if_not_installed("VariantAnnotation")
    records <- data.frame(
        CHROM = "22",
        POS = c(100L, 500L),
        ID = c("rs1", "rs5"),
        REF = "A",
        ALT = c("G", "."),
        s1 = c("0/1", "0/0"),
        s2 = c("1/1", "0/0"),
        s3 = c("0/0", "0/0"),
        stringsAsFactors = FALSE
    )
    panel <- readGenotypes(.writeTestVcf(records, samples3))
    expect_equal(rownames(panel), "rs1")
    expect_equal(.dose(panel), matrix(c(1, 2, 0), 1L))
})

test_that("vcfArgs reach the read of the records", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    # Read without row names, no record has an ID of its own.
    panel <- suppressWarnings(
        readGenotypes(vcf, vcfArgs = list(row.names = FALSE))
    )
    expect_equal(
        SummarizedExperiment::rowRanges(panel)$SNP[2:3],
        c("chr22:200:ACGTA:A", "chr22:202:G:T")
    )
})

test_that("a VCF that cannot be read lazily is refused before any read", {
    skip_if_not_installed("VariantAnnotation")
    expect_error(
        readGenotypes("/no/such/file.vcf.gz"),
        "VCF file.*File does not exist"
    )
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), samples3)
    unindexed <- file.path(tempfile("vcf_"), "copy.vcf.gz")
    dir.create(dirname(unindexed))
    file.copy(vcf, unindexed)
    expect_error(readGenotypes(unindexed), "tabix index")
    bcf <- tempfile(fileext = ".bcf")
    file.copy(vcf, bcf)
    expect_error(readGenotypes(bcf), "BCF file")
})

test_that("a VCF's dosage field must be one pecotmr can read", {
    skip_if_not_installed("VariantAnnotation")
    records <- .vcfOverlapRecords()
    perGenotype <- .writeTestVcf(
        records,
        samples3,
        format = c("GT", "DS"),
        ds = "G"
    )
    expect_error(readGenotypes(perGenotype), "FORMAT/DS with Number=G")
    noGenotypes <- .writeTestVcf(
        .vcfGtOnly(records),
        samples3,
        format = "AD",
        gt = FALSE
    )
    expect_error(
        suppressWarnings(readGenotypes(noGenotypes)),
        "neither FORMAT/DS nor FORMAT/GT"
    )
})

test_that("a sketch's .afreq beside a VCF is recorded on its panel", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords())[2:3, ], samples3)
    writeLines(
        c(
            "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\tU_MIN\tU_MAX",
            "22\trs2\tACGTA\tA\t0.25\t40\t-1\t1",
            "22\trs3\tG\tT\t0.5\t38\t-1\t1"
        ),
        sub("\\.vcf\\.gz$", ".afreq", vcf)
    )
    panel <- readGenotypes(vcf)
    expect_equal(SummarizedExperiment::rowRanges(panel)$AF, c(0.25, 0.5))
    expect_equal(S4Vectors::metadata(panel)$ldProjection$nSamples, 20L)
})

# ---- GDS -------------------------------------------------------------------

test_that("a GDS file that holds no genotypes pecotmr reads is refused", {
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    expect_error(
        readGenotypes("/no/such/file.gds"),
        "GDS file.*File does not exist"
    )
    bare <- .writeTestGds(NULL, function(gds) {
        gdsfmt::add.gdsn(gds, "x", 1L)
    })
    expect_error(readGenotypes(bare), "GDS file of format 'unknown'")
})

test_that("a SNPRelate GDS whose genotypes disagree with its SNPs is refused", {
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    f <- .writeTestGds("SNP_ARRAY", function(gds) {
        gdsfmt::add.gdsn(gds, "sample.id", c("a", "b", "c"))
        gdsfmt::add.gdsn(gds, "snp.id", c("v1", "v2"))
        gdsfmt::add.gdsn(gds, "snp.position", c(10L, 20L))
        gdsfmt::add.gdsn(gds, "snp.chromosome", c(1L, 1L))
        gdsfmt::add.gdsn(gds, "snp.allele", c("A/G", "C/T"))
        node <- gdsfmt::add.gdsn(
            gds,
            "genotype",
            matrix(0L, 3L, 3L),
            storage = "bit2"
        )
        gdsfmt::put.attr.gdsn(node, "sample.order")
    })
    expect_error(
        readGenotypes(f),
        "genotype matrix of .* is 3 x 3, but the file lists 2 SNPs"
    )
})

test_that("a SeqArray GDS without IDs names its variants by their alleles", {
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    panel <- readGenotypes(.writeFakeSeqArray())
    expect_equal(rownames(panel), "chr1:100:A:G")
    expect_equal(.dose(panel), matrix(0, 1L, 2L))
})

test_that("SeqArray calls stored otherwise than as read are refused", {
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    expect_error(
        readGenotypes(.writeFakeSeqArray(storage = "int8")),
        "stores its genotype calls as Int8"
    )
    expect_error(
        readGenotypes(.writeFakeSeqArray(nExtra = 1L)),
        "1 genotype call\\(s\\) beyond its ploidy"
    )
})

test_that("a SeqArray FORMAT/DS of another Number is refused", {
    skip_if_not_installed("SeqArray")
    skip_if_not_installed("GDSArray")
    vcf <- .writeTestVcf(
        .vcfOverlapRecords(),
        samples3,
        format = c("GT", "DS"),
        ds = "A"
    )
    gds <- .writeTestSeqArray(vcf, ds = TRUE)
    opened <- gdsfmt::openfn.gds(gds, readonly = FALSE)
    gdsfmt::put.attr.gdsn(
        gdsfmt::index.gdsn(opened, "annotation/format/DS"),
        "Number",
        "G"
    )
    gdsfmt::closefn.gds(opened)
    expect_error(readGenotypes(gds), "stores FORMAT/DS with Number=G")
})

# ---- PLINK1 ----------------------------------------------------------------

test_that("a PLINK1 panel reads the dosages snpStats reads", {
    # toy_ref has missing calls. snpStats counts the .bim's sixth allele, so
    # its count of A1, the fifth, is two minus that.
    skip_if_not_installed("PlinkMatrix")
    skip_if_not_installed("snpStats")
    stem <- str_remove(
        system.file("extdata", "toy_ref.bed", package = "pecotmr"),
        "\\.bed$"
    )
    sm <- snpStats::read.plink(
        str_c(stem, ".bed"),
        str_c(stem, ".bim"),
        str_c(stem, ".fam")
    )
    want <- t(2 - methods::as(sm$genotypes, "numeric"))
    panel <- readGenotypes(plink1Prefix = stem)
    expect_true(anyNA(want))
    expect_equal(.dose(panel), unname(want))
    expect_identical(rownames(panel), normalizeVariantId(sm$map$snp.name))
    expect_identical(colnames(panel), as.character(sm$fam$member))
    expect_identical(
        S4Vectors::metadata(panel)$genotypeSource$path,
        file.path(normalizePath(dirname(stem)), basename(stem))
    )
})

# ---- Per-chromosome sets ---------------------------------------------------

test_that("a per-chromosome VCF, GDS or PLINK1 set stacks into one panel", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("GDSArray")
    skip_if_not_installed("PlinkMatrix")
    for (ext in c(".vcf.gz", ".gds", ".bed")) {
        p21 <- file.path(test_data_dir, paste0("test_variants", ext))
        p22 <- file.path(test_data_dir, paste0("test_variants_chr22", ext))
        one21 <- readGenotypes(p21)
        one22 <- readGenotypes(p22)
        stacked <- readGenotypes(genoMeta = c("21" = p21, "22" = p22))
        n21 <- nrow(one21)
        expect_equal(dim(stacked), c(n21 + nrow(one22), ncol(one21)))
        expect_equal(rownames(stacked), c(rownames(one21), rownames(one22)))
        # A VCF is read a record at a time, so the dosages are compared on
        # each shard's first, second and last rows, either side of the join.
        rows21 <- c(1L, 2L, n21)
        rows22 <- c(1L, 2L, nrow(one22))
        dose22 <- .dose(one22[rows22, ])
        expect_equal(
            .dose(stacked[c(rows21, n21 + rows22), ]),
            rbind(.dose(one21[rows21, ]), dose22)
        )
        source <- S4Vectors::metadata(stacked)$genotypeSource
        # A PLINK1 file is recorded by its stem, as a PLINK prefix is.
        recorded <- str_remove(normalizePath(c(p21, p22)), "\\.bed$")
        expect_equal(
            source$chromPaths,
            c(chr21 = recorded[[1L]], chr22 = recorded[[2L]])
        )
        # A stacked panel has no single seed to read through.
        expect_null(pecotmr:::.panelSeed(stacked))
        # Only the files of the chromosomes asked for are read.
        only22 <- readGenotypes(
            genoMeta = c("21" = p21, "22" = p22),
            chroms = "22"
        )
        expect_equal(rownames(only22), rownames(one22))
        expect_equal(.dose(only22[rows22, ]), dose22)
    }
})

test_that("a set's shards hold the same samples, one chromosome each", {
    skip_if_not_installed("VariantAnnotation")
    records <- .vcfGtOnly(.vcfOverlapRecords())
    chr22 <- .writeTestVcf(records, samples3)
    other <- .writeTestVcf(
        `[<-`(records, "CHROM", value = "21")[c(names(records)[1:5], "s1")],
        "s1"
    )
    expect_error(
        suppressWarnings(readGenotypes(
            genoMeta = c("22" = chr22, "21" = other)
        )),
        "identical sample IDs"
    )
    again <- .writeTestVcf(records, samples3)
    expect_error(
        suppressWarnings(readGenotypes(genoMeta = c(a = chr22, b = again))),
        "appears in more than one per-chromosome file"
    )
})

# ===========================================================================
# Panels saved by an earlier version (.updatePanel)
# ===========================================================================

# The saved panels are built by helper-savedObjects.R.

test_that(".updatePanel reads a saved PLINK1 panel again through its backend", {
    stem <- test_path("test_data", "test_variants")
    fresh <- readGenotypes(stem, format = "plink1")
    saved <- .savedHandlePanel(fresh, stem, "plink1")
    SummarizedExperiment::colData(saved) <- S4Vectors::DataFrame(
        age = seq_len(ncol(fresh)),
        row.names = colnames(fresh)
    )
    # Saved as it was held: a few variants and samples, out of file order.
    saved <- saved[c(4L, 2L, 7L), c(3L, 1L)]
    # As saved, its dosages go through a reader handles no longer have.
    expect_error(
        as.matrix(SummarizedExperiment::assay(saved)),
        "Unsupported format.*plink1"
    )
    up <- .updatePanel(saved)
    expect_s4_class(.panelSeed(up), "GenotypeDosageSeed")
    expect_identical(dimnames(up), dimnames(saved))
    expect_equal(
        as.matrix(SummarizedExperiment::assay(up)),
        as.matrix(SummarizedExperiment::assay(fresh))[c(4L, 2L, 7L), c(3L, 1L)]
    )
    expect_identical(
        SummarizedExperiment::colData(up),
        SummarizedExperiment::colData(saved)
    )
    expect_identical(.ldRefSource(up)$format, "plink1")
    # Brought up to date, it is current, and a second update leaves it be.
    expect_identical(.updatePanel(up), up)
})

test_that(".updatePanel gives a saved stochastic sketch its real n and AF", {
    skip_if_not_installed("pgenlibr")
    ref <- "pecotmr://extdata/ld_reference/chr22/protocol_example.LD.chr22"
    fresh <- readGenotypes(.resolveGenotypeResourcePath(ref), format = "plink2")
    saved <- .savedHandlePanel(fresh, ref, "plink2")
    # Saved before sketches recorded them, it passes for real genotypes.
    expect_false(.ldRefIsProjection(saved))
    up <- .updatePanel(saved)
    expect_true(.ldRefIsProjection(up))
    expect_identical(.ldRefRealNSamples(up), .ldRefRealNSamples(fresh))
    expect_identical(.ldRefAf(up), .ldRefAf(fresh))
    expect_identical(dimnames(up), dimnames(saved))
    expect_identical(.updatePanel(up), up)
})

test_that(".updatePanel reads a saved per-chromosome panel again", {
    files <- c(
        chr21 = test_path("test_data", "test_variants.bed"),
        chr22 = test_path("test_data", "test_variants_chr22.bed")
    )
    fresh <- readGenotypes(genoMeta = files)
    saved <- .savedHandlePanel(fresh, "per-chromosome", "plink1", files)
    # One variant from each chromosome.
    keep <- c(2L, nrow(fresh) - 1L)
    up <- .updatePanel(saved[keep, ])
    expect_identical(dimnames(up), dimnames(fresh[keep, ]))
    expect_equal(
        as.matrix(SummarizedExperiment::assay(up)),
        as.matrix(SummarizedExperiment::assay(fresh)[keep, ])
    )
    expect_identical(.updatePanel(up), up)
})

test_that(".updatePanel takes a saved bare handle as the panel over it", {
    stem <- test_path("test_data", "test_variants")
    fresh <- readGenotypes(stem, format = "plink1")
    handle <- .ldRefHandle(.savedHandlePanel(fresh, stem, "plink1"))
    up <- .updatePanel(handle)
    expect_s4_class(.panelSeed(up), "GenotypeDosageSeed")
    expect_identical(dimnames(up), dimnames(fresh))
})

test_that(".updatePanel leaves a current panel, or no panel, as it is", {
    expect_null(.updatePanel(NULL))
    stem <- test_path("test_data", "test_variants")
    plink1 <- readGenotypes(stem, format = "plink1")
    expect_identical(.updatePanel(plink1), plink1)
    skip_if_not_installed("pgenlibr")
    plink2 <- readGenotypes(stem, format = "plink2")
    expect_identical(.updatePanel(plink2), plink2)
})

test_that(".updatePanel refuses a saved panel it cannot read again", {
    stem <- test_path("test_data", "test_variants")
    fresh <- readGenotypes(stem, format = "plink1")
    gone <- .savedHandlePanel(fresh, file.path(tempdir(), "nowhere"), "plink1")
    cnd <- tryCatch(.updatePanel(gone), error = identity)
    expect_s3_class(cnd, "pecotmrUnreadablePanel")
    expect_match(conditionMessage(cnd), "nowhere.*cannot be read")
    # Its files are where they were, but no longer hold one of its variants.
    moved <- .savedHandlePanel(fresh, stem, "plink1")
    rr <- SummarizedExperiment::rowRanges(moved)
    rr[1L] <- GenomicRanges::shift(rr[1L], 1L)
    SummarizedExperiment::rowRanges(moved) <- rr
    cnd <- tryCatch(.updatePanel(moved), error = identity)
    expect_s3_class(cnd, "pecotmrUnreadablePanel")
    expect_match(conditionMessage(cnd), "no longer hold every variant")
})
