context("GenotypeHandle constructor")

# ===========================================================================
# Fixtures: 100 samples x 349 variants on chr21:17513228-17592874.
# Available in all formats (plink1, plink2, vcf.gz, gds).
# ===========================================================================

test_data_dir <- test_path("test_data")
plink_prefix <- file.path(test_data_dir, "test_variants")
gds_path <- file.path(test_data_dir, "test_variants.gds")
vcf_path <- file.path(test_data_dir, "test_variants.vcf.gz")

# ===========================================================================
# Source-counting & input-validation branches (no file reads)
# ===========================================================================

test_that("GenotypeHandle: zero sources errors", {
    expect_error(GenotypeHandle(), "Exactly one of `path`")
})

test_that("GenotypeHandle: two simultaneous sources errors", {
    expect_error(
        GenotypeHandle(path = gds_path, plink2Prefix = plink_prefix),
        "Exactly one of `path`"
    )
})

test_that("readGenotypes: partial bed/bim/fam triplet errors", {
    # GenotypeHandle() no longer takes a PLINK1 triplet: a BED fileset is
    # read as a panel, through readGenotypes(), which checks the triplet.
    expect_error(
        readGenotypes(bed = "x.bed", bim = "x.bim"),
        "all three must be provided"
    )
})

test_that("GenotypeHandle: partial pgen/pvar/psam triplet errors", {
    expect_error(
        GenotypeHandle(pgen = "x.pgen", pvar = "x.pvar"),
        "all three must be provided"
    )
})

test_that("GenotypeHandle: ldMeta without which errors", {
    expect_error(
        GenotypeHandle(ldMeta = "x.tsv"),
        "`ldMeta` requires `block`"
    )
})

test_that("GenotypeHandle: which without ldMeta errors", {
    expect_error(
        GenotypeHandle(
            path = gds_path,
            block = GenomicRanges::GRanges("chr1:1-100")
        ),
        "`block` is only meaningful when `ldMeta`"
    )
})

# ===========================================================================
# Format-specific constructor branches (real file reads)
# ===========================================================================

test_that("GenotypeHandle refuses a VCF, GDS or PLINK1 path, read as a panel", {
    # A VCF, GDS or PLINK1 file reads through its own DelayedArray seed, not
    # a handle; asking for a handle points at readGenotypes() rather than
    # failing somewhere inside a reader.
    expect_error(GenotypeHandle(path = gds_path), "read as a genotype panel")
    expect_error(GenotypeHandle(path = vcf_path), "read as a genotype panel")
    expect_error(
        GenotypeHandle(path = paste0(plink_prefix, ".bed")),
        "read as a genotype panel"
    )
    expect_error(
        GenotypeHandle(path = gds_path, format = "gds"),
        "readGenotypes"
    )
})

test_that("GenotypeHandle: plink2Prefix builds a plink2 handle", {
    skip_if_not_installed("pgenlibr")
    h <- GenotypeHandle(plink2Prefix = plink_prefix)
    expect_s4_class(h, "GenotypeHandle")
    expect_equal(genotypeFormat(h), "plink2")
    expect_equal(nSamples(h), 100L)
})

# ===========================================================================
# Explicit triplet constructors
# ===========================================================================

test_that(".plink1TripletStem: the stem the three files share", {
    expect_equal(
        pecotmr:::.plink1TripletStem(
            bed = paste0(plink_prefix, ".bed"),
            bim = paste0(plink_prefix, ".bim"),
            fam = paste0(plink_prefix, ".fam")
        ),
        plink_prefix
    )
})

test_that(".plink1TripletStem: errors when stems disagree", {
    expect_error(
        pecotmr:::.plink1TripletStem(
            bed = "a/x.bed",
            bim = "b/y.bim",
            fam = "c/z.fam"
        ),
        "must share a common path stem"
    )
})

test_that(".plink1TripletStem: errors on non-character input", {
    expect_error(
        pecotmr:::.plink1TripletStem(
            bed = 1L,
            bim = "x.bim",
            fam = "x.fam"
        ),
        "Must be of type 'string'"
    )
})

test_that("GenotypeHandle: pgen/pvar/psam triplet with matching stems builds plink2 handle", {
    skip_if_not_installed("pgenlibr")
    h <- GenotypeHandle(
        pgen = paste0(plink_prefix, ".pgen"),
        pvar = paste0(plink_prefix, ".pvar"),
        psam = paste0(plink_prefix, ".psam")
    )
    expect_s4_class(h, "GenotypeHandle")
    expect_equal(genotypeFormat(h), "plink2")
})

test_that(".genotypeHandleFromPlink2Triplet: accepts .pvar.zst by stripping the .zst", {
    skip_if_not_installed("pgenlibr")
    # The stems agree once the .zst is stripped, so the call gets past the
    # stem check and fails only in the reader, on the files not existing.
    expect_error(
        pecotmr:::.genotypeHandleFromPlink2Triplet(
            pgen = "/tmp/x.pgen",
            pvar = "/tmp/x.pvar.zst",
            psam = "/tmp/x.psam"
        ),
        "PLINK2 .pgen file.*File does not exist"
    )
})

test_that(".genotypeHandleFromPlink2Triplet: errors when stems disagree", {
    expect_error(
        pecotmr:::.genotypeHandleFromPlink2Triplet(
            pgen = "a/x.pgen",
            pvar = "b/y.pvar",
            psam = "c/z.psam"
        ),
        "must share a common path stem"
    )
})

test_that(".genotypeHandleFromPlink2Triplet: errors on non-character input", {
    expect_error(
        pecotmr:::.genotypeHandleFromPlink2Triplet(
            pgen = 1L,
            pvar = "x.pvar",
            psam = "x.psam"
        ),
        "must be a single file path"
    )
})

# ===========================================================================
# .genotypeHandleFromLdMeta: dispatch by file extension
# ===========================================================================

# An LD-meta file whose one row names `path` for chr21:17000000-18000000.
.gh_makeLdMeta <- function(path) {
    meta <- data.frame(
        chrom = "chr21",
        start = 17000000L,
        end = 18000000L,
        path = path,
        stringsAsFactors = FALSE
    )
    f <- tempfile(fileext = ".tsv")
    write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
    f
}

test_that("an LD-meta row naming a GDS, VCF or PLINK1 file reads as a panel", {
    skip_if_not_installed("gdsfmt")
    skip_if_not_installed("GDSArray")
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("PlinkMatrix")
    region <- GenomicRanges::GRanges("chr21:17513228-17592874")
    files <- c(
        gds = gds_path,
        vcf = vcf_path,
        plink1 = paste0(plink_prefix, ".bed")
    )
    for (format in names(files)) {
        f <- .gh_makeLdMeta(files[[format]])
        on.exit(unlink(f), add = TRUE)
        expect_error(
            GenotypeHandle(ldMeta = f, block = region),
            "read as a genotype panel",
            info = format
        )
        panel <- readGenotypes(ldMeta = f, block = region)
        expect_s4_class(panel, "RangedSummarizedExperiment")
        expect_equal(dim(panel), c(349L, 100L), info = format)
        expect_equal(
            S4Vectors::metadata(panel)$genotypeSource$format,
            format,
            info = format
        )
    }
})

test_that("an LD-meta row's file is read by the reader its format names", {
    # A PLINK file or prefix an LD-meta row names is read through a handle,
    # as before; the panel route wraps it rather than reading it another way.
    skip_if_not_installed("pgenlibr")
    panel <- pecotmr:::.ghLdPathToPanel(paste0(plink_prefix, ".pgen"))
    expect_equal(S4Vectors::metadata(panel)$genotypeSource$format, "plink2")
    prefixed <- pecotmr:::.ghLdPathToPanel(plink_prefix)
    expect_equal(dim(prefixed), dim(panel))
    expect_error(
        pecotmr:::.ghLdPathToPanel(tempfile(fileext = ".txt")),
        "ldMeta: .* is neither a genotype file nor a PLINK prefix"
    )
})

test_that("GenotypeHandle ldMeta: dispatches to plink2 reader for .pgen path", {
    skip_if_not_installed("pgenlibr")
    f <- .gh_makeLdMeta(paste0(plink_prefix, ".pgen"))
    on.exit(unlink(f), add = TRUE)
    h <- GenotypeHandle(
        ldMeta = f,
        block = GenomicRanges::GRanges("chr21:17513228-17592874")
    )
    expect_s4_class(h, "GenotypeHandle")
    expect_equal(genotypeFormat(h), "plink2")
})

test_that("GenotypeHandle ldMeta: .cor.xz payload is rejected (out of scope)", {
    # findValidFilePath needs the target file to exist (otherwise it falls
    # back to returning the meta-TSV path itself, which would dead-end the
    # downstream format dispatch elsewhere). Materialize an empty .cor.xz so
    # we exercise the proper "cor.xz payload" rejection branch.
    corPath <- tempfile(fileext = ".cor.xz")
    file.create(corPath)
    meta <- .gh_makeLdMeta(corPath)
    on.exit(unlink(c(meta, corPath)), add = TRUE)
    expect_error(
        GenotypeHandle(
            ldMeta = meta,
            block = GenomicRanges::GRanges("chr21:17513228-17592874")
        ),
        "points at a pre-computed correlation matrix"
    )
})

test_that("GenotypeHandle ldMeta: region with no covering row errors", {
    f <- .gh_makeLdMeta(gds_path)
    on.exit(unlink(f), add = TRUE)
    # The upstream helper emits "No data for chromosome ..." when nothing on
    # the region's chromosome exists in the meta TSV; either phrasing
    # indicates the same "region uncovered" failure mode.
    expect_error(
        GenotypeHandle(
            ldMeta = f,
            block = GenomicRanges::GRanges("chr22:1-1000")
        ),
        "No data for chromosome|no LD-meta row covers"
    )
})

test_that("GenotypeHandle ldMeta: the window to look up is a GRanges", {
    f <- .gh_makeLdMeta(gds_path)
    on.exit(unlink(f), add = TRUE)
    expect_error(
        GenotypeHandle(ldMeta = f, block = "chr21:17513228-17592874"),
        "`block` must be a GRanges"
    )
    expect_error(GenotypeHandle(ldMeta = f), "requires `block`")
})

# === Tests migrated from test_h2ClassesSumstats.R (GenotypeHandle) ===

test_that("GenotypeHandle refuses the formats it no longer reads", {
    for (fmt in c("vcf", "gds", "plink1")) {
        expect_error(
            new(
                "GenotypeHandle",
                path = "/tmp/test",
                format = fmt,
                snpInfo = data.frame(),
                nSamples = 0L,
                sampleIds = character(),
                pgenPtr = NULL
            ),
            "'format' must be one of: plink2, bgen"
        )
    }
})

test_that("GenotypeHandle constructs and validates every format it reads", {
    for (fmt in c("plink2", "bgen")) {
        obj <- new(
            "GenotypeHandle",
            path = "/tmp/test",
            format = fmt,
            snpInfo = data.frame(),
            nSamples = 0L,
            sampleIds = character(),
            pgenPtr = NULL
        )
        expect_s4_class(obj, "GenotypeHandle")
        expect_equal(genotypeFormat(obj), fmt)
        expect_true(methods::validObject(obj))
    }
})

# ===========================================================================
# One-file-per-chromosome (sharded) handle: constructor + validation.
# genoMeta accepts a #chr,path meta file or a named chrom->path vector.
# The named-vector handle is built and checked with its extraction routing
# in test_genotypeIo.R.
# ===========================================================================
test_that("genoMeta meta-file resolves payloads relative to its own directory", {
    skip_if_not_installed("pgenlibr")
    # A meta file living in test_data referencing payloads by basename only.
    metafile <- file.path(test_data_dir, "tmp_chrom_meta_relative.tsv")
    on.exit(unlink(metafile), add = TRUE)
    writeLines(
        c(
            "#chr\tpath",
            "21\ttest_variants.pgen",
            "22\ttest_variants_chr22.pgen"
        ),
        metafile
    )
    h <- GenotypeHandle(genoMeta = metafile)
    expect_equal(sort(names(chromPaths(h))), c("chr21", "chr22"))
    expect_equal(nrow(snpInfo(h)), 2L * 349L)
})

test_that(".parseChromMeta matches chrom/path columns by name (order- and extra-column-tolerant)", {
    mk <- function(lines) {
        f <- tempfile(fileext = ".tsv")
        writeLines(lines, f)
        f
    }
    # extra legacy start/end columns + path not in position 2
    m1 <- pecotmr:::.parseChromMeta(
        mk(c(
            "#chrom\tstart\tend\tpath",
            "21\t1\t9\t/abs/a.bed",
            "22\t1\t9\t/abs/b.bed"
        ))
    )
    expect_equal(m1$chrom, c("chr21", "chr22"))
    expect_equal(m1$path, c("/abs/a.bed", "/abs/b.bed"))
    # reordered header (path first) + alias column names
    m2 <- pecotmr:::.parseChromMeta(
        mk(c("genotype\tchr", "/abs/a.bed\t21", "/abs/b.bed\t22"))
    )
    expect_equal(m2$chrom, c("chr21", "chr22"))
    expect_equal(m2$path, c("/abs/a.bed", "/abs/b.bed"))
})

test_that(".parseChromMeta falls back to the first two columns when no alias header matches", {
    mk <- function(lines) {
        f <- tempfile(fileext = ".tsv")
        writeLines(lines, f)
        f
    }
    m <- pecotmr:::.parseChromMeta(
        mk(c("foo\tbar", "21\t/abs/a.bed", "22\t/abs/b.bed"))
    )
    expect_equal(m$chrom, c("chr21", "chr22"))
    expect_equal(m$path, c("/abs/a.bed", "/abs/b.bed"))
})

test_that(".parseChromMeta errors when a chromosome maps to multiple payloads", {
    mk <- function(lines) {
        f <- tempfile(fileext = ".tsv")
        writeLines(lines, f)
        f
    }
    expect_error(
        pecotmr:::.parseChromMeta(
            mk(c("#chr\tpath", "21\t/abs/a.bed", "21\t/abs/b.bed"))
        ),
        "map to multiple genotype payloads"
    )
    # same guard on the named-vector path
    expect_error(
        pecotmr:::.parseChromMeta(c("21" = "/abs/a.bed", "21" = "/abs/b.bed")),
        "map to multiple genotype payloads"
    )
})

test_that("genoMeta errors on mismatched samples across shards", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("PlinkMatrix")
    # dummy_data and protocol_example.genotype are different sample panels:
    # one a sharded handle's check, the other a stacked PLINK1 panel's.
    expect_error(
        GenotypeHandle(
            genoMeta = c(
                "21" = file.path(test_data_dir, "test_variants.pgen"),
                "1" = file.path(test_data_dir, "dummy_data.pgen")
            )
        ),
        "identical sample IDs"
    )
    expect_error(
        readGenotypes(
            genoMeta = c(
                "21" = file.path(test_data_dir, "test_variants.bed"),
                "22" = file.path(test_data_dir, "protocol_example.genotype")
            )
        ),
        "identical sample IDs"
    )
})

test_that("genoMeta errors on mixed formats across shards", {
    # A panel's shards are told apart before any is read, so a mix is
    # refused for the mix, whichever reader each shard would take.
    expect_error(
        readGenotypes(
            genoMeta = c(
                "21" = paste0(plink_prefix, ".bed"),
                "22" = paste0(plink_prefix, ".pgen")
            )
        ),
        "share one"
    )
    expect_error(
        readGenotypes(
            genoMeta = c(
                "21" = file.path(test_data_dir, "test_variants.bed"),
                "22" = file.path(test_data_dir, "test_variants_chr22.gds")
            )
        ),
        "share one"
    )
})

test_that("genoMeta errors when a chromosome appears in two files", {
    skip_if_not_installed("pgenlibr")
    skip_if_not_installed("PlinkMatrix")
    expect_error(
        GenotypeHandle(
            genoMeta = c(
                "a" = file.path(test_data_dir, "test_variants.pgen"),
                "b" = file.path(test_data_dir, "test_variants.pgen")
            )
        ),
        "more than one"
    )
    expect_error(
        readGenotypes(
            genoMeta = c(
                "a" = file.path(test_data_dir, "test_variants.bed"),
                "b" = file.path(test_data_dir, "test_variants.bed")
            )
        ),
        "more than one"
    )
})

test_that("genoMeta sharded handle show() reports the layout", {
    skip_if_not_installed("pgenlibr")
    sh <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(test_data_dir, "test_variants.pgen"),
            "22" = file.path(test_data_dir, "test_variants_chr22.pgen")
        )
    )
    out <- paste(capture.output(show(sh)), collapse = "\n")
    expect_match(out, "per-chromosome files")
})

# ===========================================================================
# genoMeta `chroms`: read only the requested per-chromosome shards (I/O win).
# ===========================================================================
test_that("genoMeta chroms reads only the requested shard", {
    skip_if_not_installed("pgenlibr")
    meta <- c(
        "21" = file.path(test_data_dir, "test_variants.pgen"),
        "22" = file.path(test_data_dir, "test_variants_chr22.pgen")
    )
    full <- GenotypeHandle(genoMeta = meta)
    only21 <- GenotypeHandle(genoMeta = meta, chroms = "21")
    expect_equal(names(chromPaths(only21)), "chr21")
    expect_equal(nrow(snpInfo(only21)), 349L)
    # chr21 is the first shard, so the kept rows/fileIdx match the full handle's.
    expect_equal(snpInfo(only21)$SNP, snpInfo(full)$SNP[1:349])
    expect_equal(snpInfo(only21)$fileIdx, snpInfo(full)$fileIdx[1:349])
})

test_that("genoMeta chroms reads chromosome labels in any style", {
    skip_if_not_installed("pgenlibr")
    # "chr21" names the shard the panel declares as "21".
    h <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(test_data_dir, "test_variants.pgen"),
            "22" = file.path(test_data_dir, "test_variants_chr22.pgen")
        ),
        chroms = "chr21"
    )
    expect_equal(names(chromPaths(h)), "chr21")
})

test_that("genoMeta chroms falls back to all shards when none match", {
    skip_if_not_installed("pgenlibr")
    # chr9 is absent from the panel: rather than build an empty handle, read all
    # so the caller's own absence check can report the mismatch.
    h <- GenotypeHandle(
        genoMeta = c(
            "21" = file.path(test_data_dir, "test_variants.pgen"),
            "22" = file.path(test_data_dir, "test_variants_chr22.pgen")
        ),
        chroms = "9"
    )
    expect_equal(sort(names(chromPaths(h))), c("chr21", "chr22"))
})

test_that("genoMeta chroms skips a shard whose file does not exist", {
    skip_if_not_installed("pgenlibr")
    # A bogus chr21 payload proves the skip: restricting to chr22 must never open
    # it, while requesting chr21 too fails on the missing files.
    meta <- c(
        "22" = file.path(test_data_dir, "test_variants_chr22.pgen"),
        "21" = "/no/such/chr21/prefix"
    )
    h <- GenotypeHandle(genoMeta = meta, chroms = "22")
    expect_equal(names(chromPaths(h)), "chr22")
    expect_error(GenotypeHandle(genoMeta = meta, chroms = c("22", "21")))
})

test_that("chroms is rejected for a non-genoMeta source", {
    expect_error(
        readGenotypes(
            plink1Prefix = file.path(test_data_dir, "test_variants"),
            chroms = "21"
        ),
        "only supported with .genoMeta"
    )
})

# ===========================================================================
# .genotypeHandleFromLdMeta: row-resolution error branches.
# The real regionalLdMeta errors earlier for genuinely-uncovered regions
# (see the "region with no covering row" test above), so we mock it to drive
# the post-resolution branches (0 rows, >1 row, unsupported payload ext)
# without crafting whole on-disk meta layouts.
# ===========================================================================

test_that(".genotypeHandleFromLdMeta errors when no LD-meta row covers the region", {
    local_mocked_bindings(
        regionalLdMeta = function(...) {
            list(
                intersections = list(
                    LD_file_paths = character(0),
                    bimFilePaths = NULL
                )
            )
        },
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.genotypeHandleFromLdMeta("dummy.tsv", "chr1:1-100"),
        "no LD-meta row covers"
    )
})

test_that(".genotypeHandleFromLdMeta errors when the region spans multiple rows", {
    local_mocked_bindings(
        regionalLdMeta = function(...) {
            list(
                intersections = list(
                    LD_file_paths = c("/tmp/a.bed", "/tmp/b.bed"),
                    bimFilePaths = NULL
                )
            )
        },
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.genotypeHandleFromLdMeta("dummy.tsv", "chr1:1-100"),
        "spans multiple LD-meta"
    )
})

test_that(".genotypeHandleFromLdMeta errors on an unsupported payload extension", {
    local_mocked_bindings(
        regionalLdMeta = function(...) {
            list(
                intersections = list(
                    LD_file_paths = "/tmp/foo.txt",
                    bimFilePaths = NULL
                )
            )
        },
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.genotypeHandleFromLdMeta("dummy.tsv", "chr1:1-100"),
        "ldMeta: .* is neither a genotype file nor a PLINK prefix"
    )
})

# ===========================================================================
# .parseChromMeta: input-shape validation (no shard reads).
# ===========================================================================

test_that(".parseChromMeta errors on a meta file with fewer than 2 columns", {
    f <- tempfile(fileext = ".tsv")
    on.exit(unlink(f), add = TRUE)
    writeLines(c("chr", "21", "22"), f)
    expect_error(pecotmr:::.parseChromMeta(f), "must have at least 2 columns")
})

test_that(".parseChromMeta errors on an unrecognized genoMeta shape", {
    # An unnamed multi-element character vector (names forgotten) is neither a
    # meta-file path nor a chrom->path map.
    expect_error(
        pecotmr:::.parseChromMeta(c("21", "22")),
        "expected a path to a"
    )
})

# ===========================================================================
# .genotypeHandleFromChromMeta: empty meta input.
# ===========================================================================

test_that(".genotypeHandleFromChromMeta errors when the meta has no rows", {
    f <- tempfile(fileext = ".tsv")
    on.exit(unlink(f), add = TRUE)
    writeLines("#chr\tpath", f) # header only -> zero data rows
    expect_error(
        pecotmr:::.genotypeHandleFromChromMeta(f),
        "no chromosomes found"
    )
})

# ===========================================================================
# .resolveGenotypeShard: explicit-format dispatch, extension auto-detection,
# prefix sidecar probing, and the unresolvable-format error.
# ===========================================================================

test_that(".resolveGenotypeShard honours explicit plink2 format and .pgen extension", {
    skip_if_not_installed("pgenlibr")
    h1 <- pecotmr:::.resolveGenotypeShard(plink_prefix, format = "plink2")
    expect_s4_class(h1, "GenotypeHandle")
    expect_equal(genotypeFormat(h1), "plink2")
    h2 <- pecotmr:::.resolveGenotypeShard(paste0(plink_prefix, ".pgen"))
    expect_equal(genotypeFormat(h2), "plink2")
})

test_that(".shardFormat tells a shard's format without reading it", {
    f <- pecotmr:::.shardFormat
    expect_equal(f(gds_path), "gds")
    expect_equal(f(vcf_path), "vcf")
    expect_equal(f("/nonexistent/x.vcf.bgz"), "vcf")
    expect_equal(f("/nonexistent/x.bgen"), "bgen")
    # A prefix with both a .pgen and a .bed beside it is PLINK2.
    expect_equal(f(plink_prefix), "plink2")
    # A forced format wins over the extension.
    expect_equal(f(gds_path, format = "plink2"), "plink2")
})

test_that(".resolveGenotypeShard refuses a VCF, GDS or PLINK1 shard", {
    expect_error(
        pecotmr:::.resolveGenotypeShard(gds_path),
        "read as a genotype panel"
    )
    expect_error(
        pecotmr:::.resolveGenotypeShard(vcf_path, format = "vcf"),
        "read as a genotype panel"
    )
    expect_error(
        pecotmr:::.resolveGenotypeShard(paste0(plink_prefix, ".bed")),
        "read as a genotype panel"
    )
    expect_error(
        pecotmr:::.resolveGenotypeShard(plink_prefix, format = "plink1"),
        "read as a genotype panel"
    )
})

test_that(".resolveGenotypeShard probes the .pgen sidecar for an extension-less prefix", {
    # Only a .pgen sidecar exists (no .bed) so the prefix-probe routes to
    # plink2; mock the reader so we exercise the dispatch, not real IO.
    p <- tempfile()
    on.exit(unlink(paste0(p, ".pgen")), add = TRUE)
    file.create(paste0(p, ".pgen"))
    local_mocked_bindings(
        .makePlink2Handle = function(...) "STUB_PLINK2",
        .package = "pecotmr"
    )
    expect_identical(pecotmr:::.resolveGenotypeShard(p), "STUB_PLINK2")
})

test_that(".resolveGenotypeShard errors when the format cannot be determined", {
    # No recognized extension and no .bed/.pgen sidecar to probe.
    expect_error(
        pecotmr:::.resolveGenotypeShard(tempfile()),
        "genoMeta: .* is neither a genotype file nor a PLINK prefix"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(GenotypeHandle) does not error", {
    gh <- new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = makeTestSnpInfo(),
        nSamples = 100L,
        sampleIds = paste0("s", 1:100),
        pgenPtr = NULL
    )
    expect_output(show(gh), "GenotypeHandle")
})


# ---------------------------------------------------------------------------
# GhSeed / genotypeDelayedArray -- lazy dosages behind a Bioconductor assay.
#
# extractBlockGenotypes() already reads only the requested variants, so the
# seed is an adapter rather than a new reader. What needs pinning is the three
# ways an adapter like this goes quietly wrong: orientation, the meaning of a
# NULL index, and whether anything is read before it is asked for. The handle
# is test_variants' PLINK2 copy: 349 variants x 100 samples.
# ---------------------------------------------------------------------------

test_that("genotypeDelayedArray reports variants x samples", {
    # The Bioconductor assay orientation, which is the TRANSPOSE of
    # genotypes(). Backwards, this returns a plausible matrix and a
    # silently transposed LD matrix.
    gh <- .gh_pgenHandle()
    da <- genotypeDelayedArray(gh)
    expect_s4_class(da, "DelayedMatrix")
    expect_equal(dim(da), c(nrow(snpInfo(gh)), nSamples(gh)))
})

test_that("genotypeDelayedArray values match genotypes transposed", {
    gh <- .gh_pgenHandle()
    da <- genotypeDelayedArray(gh)
    truth <- t(pecotmr:::.dosageMatrix(gh, 1:5, meanImpute = FALSE))
    got <- as.matrix(da[1:5, 1:4])
    expect_equal(unname(got), unname(truth[, 1:4]))
    expect_identical(
        rownames(got),
        normalizeVariantId(snpInfo(gh)$SNP[1:5])
    )
    expect_identical(colnames(got), sampleIds(gh)[1:4])
})

test_that("a NULL index means every row or column, not none", {
    gh <- .gh_pgenHandle()
    seed <- DelayedArray::seed(genotypeDelayedArray(gh))
    full <- S4Arrays::extract_array(seed, list(NULL, NULL))
    expect_equal(dim(full), c(nrow(snpInfo(gh)), nSamples(gh)))
})

test_that("the seed answers an empty extraction on either axis", {
    gh <- .gh_pgenHandle()
    seed <- DelayedArray::seed(genotypeDelayedArray(gh))
    expect_equal(
        dim(S4Arrays::extract_array(seed, list(integer(0), NULL))),
        c(0L, nSamples(gh))
    )
    expect_equal(
        dim(S4Arrays::extract_array(seed, list(NULL, integer(0)))),
        c(nrow(snpInfo(gh)), 0L)
    )
})

test_that("genotypeDelayedArray stores the sample axis once, derives it on read", {
    # The sample ids are an anonymous projection dimension and would otherwise be
    # duplicated between GhSeed@dn[[2]] and handle@sampleIds. The seed stores only
    # the variant axis; dimnames() derives the sample axis from the handle.
    gh <- .gh_pgenHandle()
    seed <- genotypeDelayedArray(gh)@seed
    expect_null(seed@dn[[2L]]) # not stored
    expect_identical(dimnames(seed)[[2L]], as.character(sampleIds(gh))) # derived
    expect_identical(dimnames(seed)[[1L]], seed@dn[[1L]]) # variant axis stored
    expect_equal(dim(seed)[[2L]], length(dimnames(seed)[[2L]])) # dm consistent
})

test_that(".emptyGenotypeHandle empties the sample axis as well as the variants", {
    gh <- .gh_pgenHandle()
    e <- pecotmr:::.emptyGenotypeHandle(gh)
    expect_equal(nrow(snpInfo(e)), 0L)
    expect_equal(length(sampleIds(e)), 0L)
    expect_equal(nSamples(e), 0L)
    # the rebuilt seed is a consistent 0 x 0 (no 0 x nSamples mismatch)
    seed <- genotypeDelayedArray(e)@seed
    expect_equal(dim(seed), c(0L, 0L))
    expect_equal(lengths(dimnames(seed)), c(0L, 0L))
})

test_that("nothing is read until a block is actually touched", {
    # The point of the seed: a panel can be described without being read.
    gh <- .gh_pgenHandle()
    reads <- 0L
    # Captured before mocking: calling extractBlockGenotypes by name inside
    # the mock would re-enter the mock, not the reader.
    realExtract <- extractBlockGenotypes
    local_mocked_bindings(
        extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
            reads <<- reads + 1L
            realExtract(handle, snpIdx = snpIdx, meanImpute = meanImpute)
        },
        .package = "pecotmr"
    )
    da <- genotypeDelayedArray(gh)
    expect_equal(reads, 0L)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = da)
    )
    expect_equal(reads, 0L)
    invisible(as.matrix(da[1:3, 1:2]))
    expect_gt(reads, 0L)
})

test_that("only the requested variants are read", {
    gh <- .gh_pgenHandle()
    seen <- integer(0)
    realExtract <- extractBlockGenotypes
    local_mocked_bindings(
        extractBlockGenotypes = function(handle, snpIdx, meanImpute = TRUE) {
            seen <<- c(seen, length(snpIdx))
            realExtract(handle, snpIdx = snpIdx, meanImpute = meanImpute)
        },
        .package = "pecotmr"
    )
    da <- genotypeDelayedArray(gh)
    invisible(as.matrix(da[1:5, 1:4]))
    expect_true(all(seen <= 5L))
    expect_lt(sum(seen), nrow(snpInfo(gh)))
})

test_that("genotypeDelayedArray rejects a non-handle", {
    expect_error(genotypeDelayedArray("not a handle"), "GenotypeHandle")
})


test_that("the resource-path resolver passes non-scalar and NA paths through", {
    # Only a length-1 non-NA string can carry a "pecotmr://" reference, so
    # anything else is returned untouched rather than matched against.
    f <- pecotmr:::.resolveGenotypeResourcePath
    expect_identical(f(NA_character_), NA_character_)
    expect_length(f(character(0)), 0L)
    expect_identical(f(c("a", "b")), c("a", "b"))
    # ...and an ordinary scalar path still comes back unchanged.
    expect_identical(f("/tmp/plink"), "/tmp/plink")
})


test_that("an unresolvable bundled resource reference is an error", {
    expect_error(
        pecotmr:::.resolveGenotypeResourcePath(
            "pecotmr://extdata/no_such_resource_anywhere"
        ),
        "cannot resolve bundled genotype resource"
    )
})

test_that(".chromMetaPaths answers character(0) for no shards", {
    expect_identical(pecotmr:::.chromMetaPaths(list()), character(0))
})

test_that("pgenPtr returns the handle's pgen pointer", {
    # NULL for every non-pgen format: the pointer is opened lazily and only
    # by the pgen reader, so the accessor's job is to report its absence.
    h <- new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = makeTestSnpInfo(),
        nSamples = 100L,
        sampleIds = paste0("sample_", 1:100),
        pgenPtr = NULL
    )
    expect_null(pgenPtr(h))
})
