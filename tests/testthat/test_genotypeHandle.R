context("GenotypeHandle constructor")

# ===========================================================================
# Fixtures: 100 samples x 349 variants on chr21:17513228-17592874.
# Available in all formats (plink1, plink2, vcf.gz, gds).
# ===========================================================================

test_data_dir <- test_path("test_data")
plink_prefix  <- file.path(test_data_dir, "test_variants")
gds_path      <- file.path(test_data_dir, "test_variants.gds")
vcf_path      <- file.path(test_data_dir, "test_variants.vcf.gz")

# ===========================================================================
# Source-counting & input-validation branches (no file reads)
# ===========================================================================

test_that("GenotypeHandle: zero sources errors", {
  expect_error(GenotypeHandle(),
               "Exactly one of `path`")
})

test_that("GenotypeHandle: two simultaneous sources errors", {
  expect_error(GenotypeHandle(path = gds_path, plink2Prefix = plink_prefix),
               "Exactly one of `path`")
})

test_that("GenotypeHandle: partial bed/bim/fam triplet errors", {
  expect_error(GenotypeHandle(bed = "x.bed", bim = "x.bim"),
               "all three must be provided")
})

test_that("GenotypeHandle: partial pgen/pvar/psam triplet errors", {
  expect_error(GenotypeHandle(pgen = "x.pgen", pvar = "x.pvar"),
               "all three must be provided")
})

test_that("GenotypeHandle: ldMeta without region errors", {
  expect_error(GenotypeHandle(ldMeta = "x.tsv"),
               "`ldMeta` requires a `region`")
})

test_that("GenotypeHandle: region without ldMeta errors", {
  expect_error(GenotypeHandle(path = gds_path, region = "chr1:1-100"),
               "`region` is only meaningful when `ldMeta`")
})

# ===========================================================================
# Format-specific constructor branches (real file reads)
# ===========================================================================

test_that("GenotypeHandle: path = .gds uses readGenotypes (gds backend)", {
  skip_if_not_installed("SNPRelate")
  h <- GenotypeHandle(path = gds_path)
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "gds")
  expect_equal(h@nSamples, 100L)
  expect_equal(nrow(h@snpInfo), 349L)
})

test_that("GenotypeHandle: path = .vcf.gz uses readGenotypes (vcf backend)", {
  skip_if_not_installed("VariantAnnotation")
  h <- GenotypeHandle(path = vcf_path)
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "vcf")
  expect_equal(h@nSamples, 100L)
})

test_that("GenotypeHandle: plink1Prefix builds a plink1 handle", {
  skip_if_not_installed("snpStats")
  h <- GenotypeHandle(plink1Prefix = plink_prefix)
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink1")
  expect_equal(h@nSamples, 100L)
  expect_equal(nrow(h@snpInfo), 349L)
})

test_that("GenotypeHandle: plink2Prefix builds a plink2 handle", {
  skip_if_not_installed("pgenlibr")
  h <- GenotypeHandle(plink2Prefix = plink_prefix)
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink2")
  expect_equal(h@nSamples, 100L)
})

# ===========================================================================
# Explicit triplet constructors
# ===========================================================================

test_that("GenotypeHandle: bed/bim/fam triplet with matching stems builds plink1 handle", {
  skip_if_not_installed("snpStats")
  h <- GenotypeHandle(
    bed = paste0(plink_prefix, ".bed"),
    bim = paste0(plink_prefix, ".bim"),
    fam = paste0(plink_prefix, ".fam"))
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink1")
  expect_equal(h@nSamples, 100L)
})

test_that(".genotypeHandleFromPlink1Triplet: errors when stems disagree", {
  expect_error(
    pecotmr:::.genotypeHandleFromPlink1Triplet(
      bed = "a/x.bed", bim = "b/y.bim", fam = "c/z.fam"),
    "must share a common path stem"
  )
})

test_that(".genotypeHandleFromPlink1Triplet: errors on non-character input", {
  expect_error(
    pecotmr:::.genotypeHandleFromPlink1Triplet(bed = 1L, bim = "x.bim", fam = "x.fam"),
    "must be a single file path"
  )
})

test_that("GenotypeHandle: pgen/pvar/psam triplet with matching stems builds plink2 handle", {
  skip_if_not_installed("pgenlibr")
  h <- GenotypeHandle(
    pgen = paste0(plink_prefix, ".pgen"),
    pvar = paste0(plink_prefix, ".pvar"),
    psam = paste0(plink_prefix, ".psam"))
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink2")
})

test_that(".genotypeHandleFromPlink2Triplet: accepts .pvar.zst by stripping the .zst", {
  skip_if_not_installed("pgenlibr")
  # We only check that the stem-validation accepts the zst-suffixed pvar;
  # the actual file does not exist, so the downstream reader will error.
  expect_error(
    pecotmr:::.genotypeHandleFromPlink2Triplet(
      pgen = "/tmp/x.pgen", pvar = "/tmp/x.pvar.zst", psam = "/tmp/x.psam"),
    NA  # stems agree, so the stem-check should not fire
  ) |> suppressWarnings() |> tryCatch(error = function(e) {
    # The downstream reader will fail because /tmp/x.pgen doesn't exist —
    # but the failure mode is what we want to confirm is *not* the stem
    # mismatch error.
    expect_false(grepl("must share a common path stem", conditionMessage(e)))
  })
})

test_that(".genotypeHandleFromPlink2Triplet: errors when stems disagree", {
  expect_error(
    pecotmr:::.genotypeHandleFromPlink2Triplet(
      pgen = "a/x.pgen", pvar = "b/y.pvar", psam = "c/z.psam"),
    "must share a common path stem"
  )
})

test_that(".genotypeHandleFromPlink2Triplet: errors on non-character input", {
  expect_error(
    pecotmr:::.genotypeHandleFromPlink2Triplet(pgen = 1L, pvar = "x.pvar", psam = "x.psam"),
    "must be a single file path"
  )
})

# ===========================================================================
# .genotypeHandleFromLdMeta: dispatch by file extension
# ===========================================================================

.gh_makeLdMetaForGds <- function() {
  meta <- data.frame(
    chrom = "chr21",
    start = 17000000L,
    end   = 18000000L,
    path  = gds_path,
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

.gh_makeLdMetaForVcf <- function() {
  meta <- data.frame(
    chrom = "chr21",
    start = 17000000L,
    end   = 18000000L,
    path  = vcf_path,
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

.gh_makeLdMetaForBed <- function() {
  meta <- data.frame(
    chrom = "chr21",
    start = 17000000L,
    end   = 18000000L,
    path  = paste0(plink_prefix, ".bed"),
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

.gh_makeLdMetaForPgen <- function() {
  meta <- data.frame(
    chrom = "chr21",
    start = 17000000L,
    end   = 18000000L,
    path  = paste0(plink_prefix, ".pgen"),
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

.gh_makeLdMetaForCor <- function() {
  # findValidFilePath needs the target file to exist (otherwise it falls
  # back to returning the meta-TSV path itself, which would dead-end the
  # downstream format dispatch elsewhere). Materialize an empty .cor.xz so
  # we exercise the proper "cor.xz payload" rejection branch.
  corPath <- tempfile(fileext = ".cor.xz")
  file.create(corPath)
  meta <- data.frame(
    chrom = "chr21",
    start = 17000000L,
    end   = 18000000L,
    path  = corPath,
    stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".tsv")
  write.table(meta, f, sep = "\t", quote = FALSE, row.names = FALSE)
  list(meta = f, cor = corPath)
}

test_that("GenotypeHandle ldMeta: dispatches to gds reader for .gds path", {
  skip_if_not_installed("SNPRelate")
  f <- .gh_makeLdMetaForGds()
  on.exit(unlink(f), add = TRUE)
  h <- GenotypeHandle(ldMeta = f, region = "chr21:17513228-17592874")
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "gds")
})

test_that("GenotypeHandle ldMeta: dispatches to vcf reader for .vcf.gz path", {
  skip_if_not_installed("VariantAnnotation")
  f <- .gh_makeLdMetaForVcf()
  on.exit(unlink(f), add = TRUE)
  h <- GenotypeHandle(ldMeta = f, region = "chr21:17513228-17592874")
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "vcf")
})

test_that("GenotypeHandle ldMeta: dispatches to plink1 reader for .bed path", {
  skip_if_not_installed("snpStats")
  f <- .gh_makeLdMetaForBed()
  on.exit(unlink(f), add = TRUE)
  h <- GenotypeHandle(ldMeta = f, region = "chr21:17513228-17592874")
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink1")
})

test_that("GenotypeHandle ldMeta: dispatches to plink2 reader for .pgen path", {
  skip_if_not_installed("pgenlibr")
  f <- .gh_makeLdMetaForPgen()
  on.exit(unlink(f), add = TRUE)
  h <- GenotypeHandle(ldMeta = f, region = "chr21:17513228-17592874")
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink2")
})

test_that("GenotypeHandle ldMeta: .cor.xz payload is rejected (out of scope)", {
  paths <- .gh_makeLdMetaForCor()
  on.exit(unlink(c(paths$meta, paths$cor)), add = TRUE)
  expect_error(
    GenotypeHandle(ldMeta = paths$meta, region = "chr21:17513228-17592874"),
    "points at a pre-computed correlation matrix"
  )
})

test_that("GenotypeHandle ldMeta: region with no covering row errors", {
  f <- .gh_makeLdMetaForGds()
  on.exit(unlink(f), add = TRUE)
  # The upstream helper emits "No data for chromosome ..." when nothing on
  # the region's chromosome exists in the meta TSV; either phrasing
  # indicates the same "region uncovered" failure mode.
  expect_error(
    GenotypeHandle(ldMeta = f, region = "chr22:1-1000"),
    "No data for chromosome|no LD-meta row covers"
  )
})

# === Tests migrated from test_h2ClassesSumstats.R (GenotypeHandle) ===

test_that("GenotypeHandle constructs and validates correctly", {
  obj <- new("GenotypeHandle",
    path = "/tmp/test.gds",
    format = "gds",
    snpInfo = make_test_snp_info(),
    nSamples = 100L,
    sampleIds = paste0("sample_", 1:100),
    pgenPtr = NULL
  )
  expect_s4_class(obj, "GenotypeHandle")
  expect_equal(obj@format, "gds")
  expect_true(methods::validObject(obj))
})


test_that("GenotypeHandle accepts all valid formats", {
  for (fmt in c("gds", "vcf", "plink1", "plink2")) {
    obj <- new("GenotypeHandle",
      path = "/tmp/test",
      format = fmt,
      snpInfo = data.frame(),
      nSamples = 0L,
      sampleIds = character(),
      pgenPtr = NULL
    )
    expect_true(methods::validObject(obj))
  }
})


test_that("GenotypeHandle rejects invalid format", {
  expect_error(
    methods::validObject(
      new("GenotypeHandle",
        path = "/tmp/test",
        format = "bgen",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
      )
    ),
    "format.*must be one of"
  )
})

# ===========================================================================
# One-file-per-chromosome (sharded) handle: constructor + validation.
# genoMeta accepts a #chr,path meta file or a named chrom->path vector.
# Extraction routing for sharded handles is covered in test_genotypeIo.R.
# ===========================================================================
test_that("genoMeta (named vector) builds a sharded handle", {
  skip_if_not_installed("snpStats")
  h <- GenotypeHandle(genoMeta = c(
    "21" = file.path(test_data_dir, "test_variants"),
    "22" = file.path(test_data_dir, "test_variants_chr22")))
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink1")
  expect_equal(sort(names(h@chromPaths)), c("21", "22"))
  expect_equal(nrow(h@snpInfo), 2L * 349L)
})

test_that("genoMeta meta-file resolves payloads relative to its own directory", {
  skip_if_not_installed("snpStats")
  # A meta file living in test_data referencing payloads by basename only.
  metafile <- file.path(test_data_dir, "tmp_chrom_meta_relative.tsv")
  on.exit(unlink(metafile), add = TRUE)
  writeLines(c("#chr\tpath", "21\ttest_variants", "22\ttest_variants_chr22"),
             metafile)
  h <- GenotypeHandle(genoMeta = metafile)
  expect_equal(sort(names(h@chromPaths)), c("21", "22"))
  expect_equal(nrow(h@snpInfo), 2L * 349L)
})

test_that("genoMeta errors on mismatched samples across shards", {
  skip_if_not_installed("snpStats")
  # protocol_example.genotype is a different sample panel.
  expect_error(
    GenotypeHandle(genoMeta = c(
      "21" = file.path(test_data_dir, "test_variants"),
      "22" = file.path(test_data_dir, "protocol_example.genotype"))),
    "identical sample IDs")
})

test_that("genoMeta errors on mixed formats across shards", {
  skip_if_not_installed("snpStats")
  skip_if_not_installed("SNPRelate")
  expect_error(
    GenotypeHandle(genoMeta = c(
      "21" = file.path(test_data_dir, "test_variants"),
      "22" = file.path(test_data_dir, "test_variants_chr22.gds"))),
    "share one")
})

test_that("genoMeta errors when a chromosome appears in two files", {
  skip_if_not_installed("snpStats")
  expect_error(
    GenotypeHandle(genoMeta = c(
      "a" = file.path(test_data_dir, "test_variants"),
      "b" = file.path(test_data_dir, "test_variants"))),
    "more than one")
})

test_that("genoMeta sharded handle show() reports the layout", {
  skip_if_not_installed("snpStats")
  sh <- GenotypeHandle(genoMeta = c(
    "21" = file.path(test_data_dir, "test_variants"),
    "22" = file.path(test_data_dir, "test_variants_chr22")))
  out <- paste(capture.output(show(sh)), collapse = "\n")
  expect_match(out, "per-chromosome files")
})


