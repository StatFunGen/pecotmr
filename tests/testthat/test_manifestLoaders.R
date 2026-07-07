# Tests for R/manifestLoaders.R

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# The toy_ref PLINK1 panel (chr22) serves as the LD sketch. These variants
# are drawn verbatim from toy_ref.bim so LD-containment succeeds.
.toyRefPrefix <- function() {
  sub("\\.bed$", "", system.file("extdata", "toy_ref.bed", package = "pecotmr"))
}

.toyLdSketch <- function() {
  GenotypeHandle(plink1Prefix = .toyRefPrefix())
}

# A small GWAS sumstats data.frame over real toy_ref chr22 variants.
.toyGwasDf <- function(n = 5) {
  base <- data.frame(
    chrom      = rep("22", 5),
    pos        = c(14560203, 14564328, 14850625, 14870204, 14878387),
    variant_id = c("rs11089128", "rs7288972", "rs11167319", "rs8138488", "rs2186521"),
    A1         = c("G", "C", "G", "C", "A"),
    A2         = c("A", "T", "T", "T", "G"),
    z          = c(1.2, -0.4, 2.1, 0.05, -1.8),
    n_sample   = rep(10000L, 5),
    stringsAsFactors = FALSE)
  base[seq_len(n), , drop = FALSE]
}

# Write a data.frame to a TSV whose header line starts with "#" (so a tabix
# index can treat it as a comment/header line).
.writeSumstatsTsv <- function(df, path) {
  hdr <- names(df); hdr[[1L]] <- paste0("#", hdr[[1L]])
  con <- file(path, "w")
  writeLines(paste(hdr, collapse = "\t"), con)
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)
  close(con)
  path
}

# ---------------------------------------------------------------------------
# Column canonicalisation + scalar reconciliation
# ---------------------------------------------------------------------------

test_that(".canonManifestCols renames snake_case aliases to camelCase", {
  df <- data.frame(study_id = "s1", sumstats_path = "x.tsv",
                   stringsAsFactors = FALSE)
  out <- pecotmr:::.canonManifestCols(df, pecotmr:::.gwasSumStatsManifestAliases,
                                      required = c("study", "sumStatsPath"),
                                      label = "GwasSumStats")
  expect_true(all(c("study", "sumStatsPath") %in% names(out)))
  expect_false(any(c("study_id", "sumstats_path") %in% names(out)))
})

test_that(".canonManifestCols errors on missing required columns", {
  df <- data.frame(study = "s1", stringsAsFactors = FALSE)
  expect_error(
    pecotmr:::.canonManifestCols(df, pecotmr:::.gwasSumStatsManifestAliases,
                                 required = c("study", "sumStatsPath"),
                                 label = "GwasSumStats"),
    "missing required column")
})

test_that(".reconcileScalar resolves arg/column and flags conflicts", {
  expect_equal(pecotmr:::.reconcileScalar(NULL, "hg38", "genome"), "hg38")
  expect_equal(pecotmr:::.reconcileScalar(c("hg38", "hg38"), NULL, "genome"), "hg38")
  expect_equal(pecotmr:::.reconcileScalar(c("hg38", "hg38"), "hg38", "genome"), "hg38")
  expect_error(pecotmr:::.reconcileScalar(c("hg19", "hg38"), NULL, "genome"),
               "must be constant")
  expect_error(pecotmr:::.reconcileScalar("hg19", "hg38", "genome"), "disagrees")
  expect_error(pecotmr:::.reconcileScalar(NULL, NULL, "genome"), "must be provided")
})

# ---------------------------------------------------------------------------
# Genotype-format detection
# ---------------------------------------------------------------------------

test_that(".detectGenotypeFormat builds a PLINK1 handle from a prefix", {
  h <- pecotmr:::.detectGenotypeFormat(.toyRefPrefix())
  expect_s4_class(h, "GenotypeHandle")
  expect_equal(h@format, "plink1")
})

test_that("BCF sumstats are rejected", {
  expect_error(
    pecotmr:::.readSumStatsFile("x.bcf", NULL, NULL, NULL, NULL, "lbl"),
    "BCF sumstats are not supported")
})

# ---------------------------------------------------------------------------
# LD-sketch containment
# ---------------------------------------------------------------------------

test_that(".checkLdContainment passes for overlapping variants", {
  gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(5))
  expect_silent(pecotmr:::.checkLdContainment(.toyLdSketch(), gr, 0.5, "lbl"))
})

test_that(".checkLdContainment errors on a chromosome absent from the sketch", {
  df <- .toyGwasDf(3); df$chrom <- "1"
  gr <- pecotmr:::.dfToEntryGranges(df)
  expect_error(pecotmr:::.checkLdContainment(.toyLdSketch(), gr, 0.5, "lbl"),
               "absent from the LD sketch")
})

test_that(".checkLdContainment errors on zero overlap", {
  df <- .toyGwasDf(2); df$pos <- c(999999001, 999999002)
  gr <- pecotmr:::.dfToEntryGranges(df)
  expect_error(pecotmr:::.checkLdContainment(.toyLdSketch(), gr, 0.5, "lbl"),
               "none of the")
})

test_that(".checkLdContainment warns on low overlap", {
  df <- .toyGwasDf(4)
  df$pos[3:4] <- c(999999001, 999999002)  # 2/4 present -> 50%
  gr <- pecotmr:::.dfToEntryGranges(df)
  expect_warning(pecotmr:::.checkLdContainment(.toyLdSketch(), gr, 0.9, "lbl"),
                 "are present in the LD sketch")
})

# ---------------------------------------------------------------------------
# GwasSumStats loader
# ---------------------------------------------------------------------------

test_that("loadGwasSumStatsFromManifest builds from a data.frame manifest", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
  manifest <- data.frame(study = "study1", sumStatsPath = ssPath,
                         stringsAsFactors = FALSE)
  obj <- loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                      ldSketch = .toyLdSketch())
  expect_s4_class(obj, "GwasSumStats")
  expect_true(methods::validObject(obj))
  expect_equal(as.character(obj$study), "study1")
  expect_equal(getGenome(obj), "hg38")
  expect_equal(length(obj$entry[[1L]]), 5L)
  expect_true(all(c("SNP", "A1", "A2", "Z", "N") %in%
                  colnames(S4Vectors::mcols(obj$entry[[1L]]))))
  expect_equal(length(obj@qcInfo), 0L)  # loaders run no QC
})

test_that("loadGwasSumStatsFromManifest reads a manifest file and reconciles genome", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
  mfPath <- file.path(tmp, "gwas_manifest.tsv")
  readr::write_tsv(data.frame(study = "study1", sumStatsPath = ssPath,
                              genome = "hg38", stringsAsFactors = FALSE), mfPath)
  # genome column present; matching arg is fine.
  obj <- loadGwasSumStatsFromManifest(mfPath, genome = "hg38",
                                      ldSketch = .toyLdSketch())
  expect_equal(getGenome(obj), "hg38")
  # conflicting arg errors.
  expect_error(loadGwasSumStatsFromManifest(mfPath, genome = "hg19",
                                            ldSketch = .toyLdSketch()),
               "disagrees")
})

test_that("loadGwasSumStatsFromManifest rejects duplicate studies", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "s.tsv"))
  manifest <- data.frame(study = c("a", "a"),
                         sumStatsPath = c(ssPath, ssPath),
                         stringsAsFactors = FALSE)
  expect_error(loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                            ldSketch = .toyLdSketch()),
               "must be unique")
})

# ---------------------------------------------------------------------------
# Region reading (tabix vs plain text)
# ---------------------------------------------------------------------------

test_that("tabix-indexed text honours the region", {
  skip_if_not_installed("Rsamtools")
  tmp <- withr::local_tempdir()
  tsv <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
  gz <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(gz, seq = 1L, start = 2L, end = 2L, comment = "#")
  df <- pecotmr:::.readSumStatsText(gz, "22:14560000-14565000", NULL, "lbl")
  expect_equal(nrow(df), 2L)  # first two variants only
  expect_true(all(df$pos >= 14560000 & df$pos <= 14565000))
})

test_that("plain-text region is ignored with a warning", {
  tmp <- withr::local_tempdir()
  tsv <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
  expect_warning(
    df <- pecotmr:::.readSumStatsText(tsv, "22:14560000-14565000", NULL, "lbl"),
    "region ignored")
  expect_equal(nrow(df), 5L)
})

# ---------------------------------------------------------------------------
# YAML column mapping
# ---------------------------------------------------------------------------

test_that("loadGwasSumStatsFromManifest resolves a YAML column mapping", {
  tmp <- withr::local_tempdir()
  # Sumstats with non-standard source column names.
  df <- .toyGwasDf(5)
  names(df) <- c("CHR", "POS", "RSID", "EA", "OA", "ZSCORE", "SAMPLE_N")
  ssPath <- .writeSumstatsTsv(df, file.path(tmp, "study1.tsv"))
  mapPath <- file.path(tmp, "mapping.yaml")
  yaml::write_yaml(list(chrom = "CHR", pos = "POS", variant_id = "RSID",
                        A1 = "EA", A2 = "OA", z = "ZSCORE",
                        n_sample = "SAMPLE_N"), mapPath)
  manifest <- data.frame(study = "study1", sumStatsPath = ssPath,
                         columnMapping = mapPath, stringsAsFactors = FALSE)
  obj <- loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                      ldSketch = .toyLdSketch())
  expect_s4_class(obj, "GwasSumStats")
  expect_equal(length(obj$entry[[1L]]), 5L)
  expect_equal(as.character(S4Vectors::mcols(obj$entry[[1L]])$A1),
               c("G", "C", "G", "C", "A"))
})

test_that(".readColumnMapping accepts a named list and a YAML file alike", {
  tmp <- withr::local_tempdir()
  mp <- file.path(tmp, "m.yaml")
  yaml::write_yaml(list(chrom = "CHR", z = "ZSCORE"), mp)
  fromFile <- pecotmr:::.readColumnMapping(mp)
  fromList <- pecotmr:::.readColumnMapping(list(chrom = "CHR", z = "ZSCORE"))
  expect_equal(fromFile[["chrom"]], "CHR")
  expect_equal(fromFile[["z"]], "ZSCORE")
  expect_equal(fromList, fromFile)
})

# ---------------------------------------------------------------------------
# QtlSumStats loader
# ---------------------------------------------------------------------------

test_that("loadQtlSumStatsFromManifest builds a per-tuple collection", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "geneA.tsv"))
  manifest <- data.frame(study = "eqtl", context = "Whole_Blood",
                         trait = "geneA", sumStatsPath = ssPath,
                         stringsAsFactors = FALSE)
  obj <- loadQtlSumStatsFromManifest(manifest, genome = "hg38",
                                     ldSketch = .toyLdSketch())
  expect_s4_class(obj, "QtlSumStats")
  expect_true(methods::validObject(obj))
  expect_equal(as.character(obj$trait), "geneA")
  expect_equal(as.character(obj$context), "Whole_Blood")
})

# ---------------------------------------------------------------------------
# QtlDataset loader
# ---------------------------------------------------------------------------

# Build a tiny 2-trait phenotype BED for one context.
.writePhenoBed <- function(path, samples = paste0("S", 1:6)) {
  df <- data.frame(
    `#chr` = c("chr22", "chr22"),
    start  = c(14560000L, 14870000L),
    end    = c(14560001L, 14870001L),
    gene_id = c("geneA", "geneB"),
    check.names = FALSE, stringsAsFactors = FALSE)
  for (s in samples) df[[s]] <- stats::rnorm(2)
  readr::write_tsv(df, path)
  path
}

test_that("loadQtlDatasetFromManifest builds a single-study dataset", {
  tmp <- withr::local_tempdir()
  bed <- .writePhenoBed(file.path(tmp, "ctxA.bed"))
  manifest <- data.frame(context = "ctxA", phenotypePath = bed,
                         study = "studyX", genotypePath = .toyRefPrefix(),
                         stringsAsFactors = FALSE)
  qd <- loadQtlDatasetFromManifest(manifest)
  expect_s4_class(qd, "QtlDataset")
  expect_true(methods::validObject(qd))
  expect_equal(getStudy(qd), "studyX")
  expect_equal(getContexts(qd), "ctxA")
})

test_that("loadQtlDatasetFromManifest accepts study/genotypes as arguments", {
  tmp <- withr::local_tempdir()
  bed <- .writePhenoBed(file.path(tmp, "ctxA.bed"))
  manifest <- data.frame(context = "ctxA", phenotypePath = bed,
                         stringsAsFactors = FALSE)
  qd <- loadQtlDatasetFromManifest(manifest, study = "studyY",
                                   genotypes = .toyLdSketch())
  expect_equal(getStudy(qd), "studyY")
})

# ---------------------------------------------------------------------------
# MultiStudyQtlDataset loader
# ---------------------------------------------------------------------------

test_that("loadMultiStudyQtlDatasetFromManifest builds from >=2 studies", {
  tmp <- withr::local_tempdir()
  bedA <- .writePhenoBed(file.path(tmp, "a.bed"))
  bedB <- .writePhenoBed(file.path(tmp, "b.bed"))
  manifest <- data.frame(
    study = c("study1", "study2"),
    context = c("ctxA", "ctxA"),
    phenotypePath = c(bedA, bedB),
    genotypePath = rep(.toyRefPrefix(), 2),
    stringsAsFactors = FALSE)
  msd <- loadMultiStudyQtlDatasetFromManifest(manifest)
  expect_s4_class(msd, "MultiStudyQtlDataset")
  expect_true(methods::validObject(msd))
  expect_equal(sort(names(msd@qtlDatasets)), c("study1", "study2"))
})

test_that("loadMultiStudyQtlDatasetFromManifest attaches a summary-only study", {
  tmp <- withr::local_tempdir()
  bedA <- .writePhenoBed(file.path(tmp, "a.bed"))
  qdMan <- data.frame(study = "study1", context = "ctxA", phenotypePath = bedA,
                      genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "geneA.tsv"))
  ssMan <- data.frame(study = "eqtl2", context = "Blood", trait = "geneA",
                      sumStatsPath = ssPath, stringsAsFactors = FALSE)
  msd <- loadMultiStudyQtlDatasetFromManifest(
    qdMan, sumStatsManifest = ssMan, genome = "hg38", ldSketch = .toyLdSketch())
  expect_true(methods::validObject(msd))
  expect_s4_class(msd@sumStats, "QtlSumStats")
  expect_equal(names(msd@qtlDatasets), "study1")
})

# ===========================================================================
# Coverage: manifest ingestion, path resolution, scalar reconciliation
# ===========================================================================

test_that(".readManifest handles CSV, missing files, and bad input", {
  tmp <- withr::local_tempdir()
  csv <- file.path(tmp, "m.csv")
  ssPath <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "s.tsv"))
  readr::write_csv(data.frame(study = "s1", sumStatsPath = ssPath), csv)
  m <- pecotmr:::.readManifest(csv)
  expect_true(is.data.frame(m) && m$study == "s1")
  expect_error(pecotmr:::.readManifest("/no/such/manifest.tsv"), "not found")
  expect_error(pecotmr:::.readManifest(42L), "data.frame or a single file path")
})

test_that(".resolveRel resolves relative paths against a base only", {
  expect_equal(pecotmr:::.resolveRel("a.tsv", "/base"), "/base/a.tsv")
  expect_equal(pecotmr:::.resolveRel("/abs/a.tsv", "/base"), "/abs/a.tsv")
  expect_equal(pecotmr:::.resolveRel("a.tsv", NULL), "a.tsv")
  expect_equal(pecotmr:::.resolveRel("", "/base"), "")
})

test_that("loadGwasSumStatsFromManifest resolves manifest-relative paths", {
  tmp <- withr::local_tempdir()
  .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
  mfPath <- file.path(tmp, "manifest.tsv")
  # Relative sumStatsPath resolves against the manifest's own directory.
  readr::write_tsv(data.frame(study = "study1", sumStatsPath = "study1.tsv"),
                   mfPath)
  obj <- loadGwasSumStatsFromManifest(mfPath, genome = "hg38",
                                      ldSketch = .toyLdSketch())
  expect_equal(length(obj$entry[[1L]]), 5L)
})

test_that(".reconcileScalar returns NULL for an optional unset scalar", {
  expect_null(pecotmr:::.reconcileScalar(NULL, NULL, "x", required = FALSE))
})

# ===========================================================================
# Coverage: genotype-format detection and LD-sketch resolution
# ===========================================================================

test_that(".detectGenotypeFormat dispatches by extension", {
  expect_error(pecotmr:::.detectGenotypeFormat("x.vcf.gz"))   # vcf branch
  expect_error(pecotmr:::.detectGenotypeFormat("x.pgen"))     # plink2 branch
  expect_error(pecotmr:::.detectGenotypeFormat("x.unknown"),
               "Could not determine genotype format")
  h <- pecotmr:::.detectGenotypeFormat(paste0(.toyRefPrefix(), ".bed"))
  expect_equal(h@format, "plink1")
})

test_that(".resolveLdSketch accepts a genoMeta vector, a path, and rejects bad input", {
  sharded <- pecotmr:::.resolveLdSketch(c("22" = .toyRefPrefix()))
  expect_s4_class(sharded, "GenotypeHandle")
  expect_true("22" %in% names(sharded@chromPaths))
  expect_s4_class(pecotmr:::.resolveLdSketch(.toyRefPrefix()), "GenotypeHandle")
  expect_error(pecotmr:::.resolveLdSketch("nonexistent_meta.tsv"))
  expect_error(pecotmr:::.resolveLdSketch(42L), "must be a GenotypeHandle")
})

test_that("ldSketch is resolved from an ldSketchPath column and conflicts error", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
  manifest <- data.frame(study = "s1", sumStatsPath = ssPath,
                         ldSketchPath = .toyRefPrefix(), stringsAsFactors = FALSE)
  obj <- loadGwasSumStatsFromManifest(manifest, genome = "hg38")
  expect_s4_class(obj, "GwasSumStats")
  # Conflicting arg vs column.
  expect_error(
    loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                 ldSketch = "/some/other/prefix"),
    "disagrees")
  # Neither arg nor column.
  bare <- data.frame(study = "s1", sumStatsPath = ssPath, stringsAsFactors = FALSE)
  expect_error(loadGwasSumStatsFromManifest(bare, genome = "hg38"),
               "must be provided as an argument or an .ldSketchPath")
})

test_that("LD containment uses a sharded sketch's chromPaths", {
  sharded <- GenotypeHandle(genoMeta = c("22" = .toyRefPrefix()))
  gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(5))
  expect_silent(pecotmr:::.checkLdContainment(sharded, gr, 0.5, "lbl"))
})

test_that(".checkLdContainment short-circuits on an empty entry", {
  expect_null(pecotmr:::.checkLdContainment(.toyLdSketch(),
                                            GenomicRanges::GRanges(), 0.5, "lbl"))
})

test_that(".checkLdContainment warns when the sketch has no variant metadata", {
  h <- methods::new("GenotypeHandle", path = "meta", format = "plink1",
                    snpInfo = data.frame(), nSamples = 0L,
                    sampleIds = character(), pgenPtr = NULL,
                    chromPaths = c("22" = "/p"))
  gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(3))
  expect_warning(pecotmr:::.checkLdContainment(h, gr, 0.5, "lbl"),
                 "no variant metadata")
})

# ===========================================================================
# Coverage: region coercion and column-mapping errors
# ===========================================================================

test_that(".asGRegion coerces strings, GRanges, and data.frames", {
  gr <- GenomicRanges::GRanges("22", IRanges::IRanges(1, 100))
  expect_identical(pecotmr:::.asGRegion(gr), gr)
  fromStr <- pecotmr:::.asGRegion("22:1-100")
  expect_equal(GenomicRanges::start(fromStr), 1L)
  fromDf <- pecotmr:::.asGRegion(data.frame(chrom = "22", start = 1, end = 100))
  expect_equal(as.character(GenomicRanges::seqnames(fromDf)), "22")
  expect_error(pecotmr:::.asGRegion("not-a-region"), "chr:start-end")
  expect_error(pecotmr:::.asGRegion(42L), "chr:start-end")
})

test_that(".readColumnMapping errors on bad inputs", {
  tmp <- withr::local_tempdir()
  expect_error(pecotmr:::.readColumnMapping("/no/map.yaml"), "not found")
  bad <- file.path(tmp, "bad.yaml")
  yaml::write_yaml(list("a", "b"), bad)  # unnamed sequence
  expect_error(pecotmr:::.readColumnMapping(bad), "standardName: sourceName")
  expect_error(pecotmr:::.readColumnMapping(42L), "named list/vector or a path")
})

test_that(".resolveSumstatCols errors on bad mapping and missing fields", {
  df <- .toyGwasDf(3)
  names(df) <- c("chrom", "pos", "variant_id", "A1", "A2", "z", "n_sample")
  expect_error(
    pecotmr:::.resolveSumstatCols(df, list(chrom = "NOPE"), "lbl"),
    "is not a column")
  expect_error(
    pecotmr:::.resolveSumstatCols(df[, c("chrom", "pos")], NULL, "lbl"),
    "missing required field")
})

# ===========================================================================
# Coverage: nCase/nControl/varY columns
# ===========================================================================

test_that("loadGwasSumStatsFromManifest attaches nCase/nControl/varY", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
  manifest <- data.frame(study = "cc", sumStatsPath = ssPath,
                         nCase = 1200, nControl = 3400, varY = 0.19,
                         stringsAsFactors = FALSE)
  obj <- loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                      ldSketch = .toyLdSketch())
  expect_equal(as.numeric(obj$nCase), 1200)
  expect_equal(as.numeric(obj$nControl), 3400)
  expect_equal(as.numeric(obj$varY), 0.19)
})

test_that("loadQtlSumStatsFromManifest attaches varY and honours a region arg", {
  tmp <- withr::local_tempdir()
  ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "geneA.tsv"))
  manifest <- data.frame(study = "eqtl", context = "Blood", trait = "geneA",
                         sumStatsPath = ssPath, varY = 1.7,
                         stringsAsFactors = FALSE)
  obj <- suppressWarnings(loadQtlSumStatsFromManifest(
    manifest, genome = "hg38", ldSketch = .toyLdSketch(),
    region = data.frame(chrom = "22", start = 1, end = 5e7)))
  expect_equal(as.numeric(obj$varY), 1.7)
})

# ===========================================================================
# Coverage: QtlDataset covariates and error branches
# ===========================================================================

.writeCovTsv <- function(path, samples, npc = 3, transpose = FALSE) {
  set.seed(1)
  m <- matrix(stats::rnorm(length(samples) * npc), nrow = length(samples),
              dimnames = list(samples, paste0("PC", seq_len(npc))))
  if (transpose) {
    out <- data.frame(covariate = colnames(m), t(m), check.names = FALSE)
  } else {
    out <- data.frame(sample = rownames(m), m, check.names = FALSE)
  }
  readr::write_tsv(out, path)
  path
}

test_that("loadQtlDatasetFromManifest attaches per-context and genotype covariates", {
  tmp <- withr::local_tempdir()
  samples <- paste0("S", 1:6)
  bed  <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
  pcov <- .writeCovTsv(file.path(tmp, "pcov.tsv"), samples, npc = 3)
  gcov <- .writeCovTsv(file.path(tmp, "gcov.tsv"), samples, npc = 2)
  manifest <- data.frame(context = "ctx", phenotypePath = bed,
                         covariatePath = pcov, study = "s",
                         genotypePath = .toyRefPrefix(),
                         genotypeCovariatePath = gcov, stringsAsFactors = FALSE)
  qd <- loadQtlDatasetFromManifest(manifest)
  expect_equal(ncol(getPhenotypeCovariates(qd, "ctx")[["ctx"]]), 3L)
  expect_equal(ncol(getGenotypeCovariates(qd)), 2L)
})

test_that("loadQtlDatasetFromManifest reads transposed (QTLtools) covariates", {
  tmp <- withr::local_tempdir()
  samples <- paste0("S", 1:6)
  bed  <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
  pcov <- .writeCovTsv(file.path(tmp, "pcovT.tsv"), samples, npc = 4,
                       transpose = TRUE)
  manifest <- data.frame(context = "ctx", phenotypePath = bed,
                         covariatePath = pcov, study = "s",
                         genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  qd <- loadQtlDatasetFromManifest(manifest, transposeCovariates = TRUE)
  expect_equal(ncol(getPhenotypeCovariates(qd, "ctx")[["ctx"]]), 4L)
})

test_that("loadQtlDatasetFromManifest accepts a genotype prefix and matrix covariates", {
  tmp <- withr::local_tempdir()
  samples <- paste0("S", 1:6)
  bed <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
  manifest <- data.frame(context = "ctx", phenotypePath = bed, study = "s",
                         stringsAsFactors = FALSE)
  gmat <- matrix(0, nrow = 6, ncol = 2, dimnames = list(samples, c("A", "B")))
  qd <- loadQtlDatasetFromManifest(manifest, genotypes = .toyRefPrefix(),
                                   genotypeCovariates = gmat)
  expect_equal(ncol(getGenotypeCovariates(qd)), 2L)
})

test_that("QtlDataset builder errors on inconsistent per-context paths", {
  tmp <- withr::local_tempdir()
  bed1 <- .writePhenoBed(file.path(tmp, "c1.bed"))
  bed2 <- .writePhenoBed(file.path(tmp, "c2.bed"))
  # Same context, two different phenotype paths.
  manifest <- data.frame(context = c("ctx", "ctx"),
                         phenotypePath = c(bed1, bed2), study = "s",
                         genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  expect_error(loadQtlDatasetFromManifest(manifest),
               "exactly one phenotypePath")
})

# ===========================================================================
# Coverage: GWAS-VCF reading
# ===========================================================================

# Write a minimal GWAS-VCF (ES/SE/LP/SS/EAF FORMAT) over toy_ref chr22 sites.
# A1 (effect) = ALT, A2 = REF, matching the toy_ref alleles.
.writeGwasVcf <- function(path, samples = "study1", n = 5, dropSs = FALSE,
                          dropEs = FALSE) {
  pos <- c(14560203, 14564328, 14850625, 14870204, 14878387)[seq_len(n)]
  ids <- c("rs11089128", "rs7288972", "rs11167319", "rs8138488", "rs2186521")[seq_len(n)]
  ref <- c("A", "T", "T", "T", "G")[seq_len(n)]
  alt <- c("G", "C", "G", "C", "A")[seq_len(n)]
  es  <- c(0.12, -0.04, 0.21, 0.005, -0.18)[seq_len(n)]
  se  <- rep(0.1, n)
  lp  <- c(1.3, 0.2, 2.1, 0.05, 1.0)[seq_len(n)]
  eaf <- c(0.3, 0.4, 0.2, 0.5, 0.1)[seq_len(n)]
  # Build the FORMAT tag list and matching per-record cells dynamically.
  tags <- c(if (!dropEs) c("ES", "SE"), "LP", if (!dropSs) "SS", "EAF")
  hdrFor <- c(
    ES  = '##FORMAT=<ID=ES,Number=A,Type=Float,Description="effect">',
    SE  = '##FORMAT=<ID=SE,Number=A,Type=Float,Description="se">',
    LP  = '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10p">',
    SS  = '##FORMAT=<ID=SS,Number=A,Type=Float,Description="n">',
    EAF = '##FORMAT=<ID=EAF,Number=A,Type=Float,Description="eaf">')
  meta <- c("##fileformat=VCFv4.2", unname(hdrFor[tags]), "##contig=<ID=22>",
            paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "INFO", "FORMAT", samples), collapse = "\t"))
  rows <- vapply(seq_len(n), function(i) {
    vals <- c(ES = sprintf("%.4f", es[i]), SE = sprintf("%.4f", se[i]),
              LP = sprintf("%.4f", lp[i]), SS = "10000",
              EAF = sprintf("%.4f", eaf[i]))
    cell <- paste(vals[tags], collapse = ":")
    geno <- paste(rep(cell, length(samples)), collapse = "\t")
    paste(c("22", pos[i], ids[i], ref[i], alt[i], ".", "PASS", ".",
            paste(tags, collapse = ":"), geno), collapse = "\t")
  }, character(1))
  writeLines(c(meta, rows), path)
  path
}

test_that("loadGwasSumStatsFromManifest reads a GWAS-VCF (A1=ALT, Z=ES/SE)", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "study1.vcf"))
  manifest <- data.frame(study = "study1", sumStatsPath = vcf,
                         stringsAsFactors = FALSE)
  obj <- suppressWarnings(loadGwasSumStatsFromManifest(
    manifest, genome = "hg38", ldSketch = .toyLdSketch()))
  mc <- S4Vectors::mcols(obj$entry[[1L]])
  expect_equal(as.character(mc$A1), c("G", "C", "G", "C", "A"))  # ALT
  expect_equal(as.character(mc$A2), c("A", "T", "T", "T", "G"))  # REF
  expect_equal(mc$Z, c(0.12, -0.04, 0.21, 0.005, -0.18) / 0.1, tolerance = 1e-3)
  expect_equal(mc$N, rep(10000, 5))
})

test_that("GWAS-VCF sample selection is enforced for multi-sample files", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "multi.vcf"),
                       samples = c("study1", "study2"))
  manifest <- data.frame(study = "study1", sumStatsPath = vcf,
                         stringsAsFactors = FALSE)
  expect_error(
    suppressWarnings(loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                                  ldSketch = .toyLdSketch())),
    "pass .sampleSelect")
  obj <- suppressWarnings(loadGwasSumStatsFromManifest(
    manifest, genome = "hg38", ldSketch = .toyLdSketch(),
    sampleSelect = "study2"))
  expect_equal(length(obj$entry[[1L]]), 5L)
  expect_error(
    suppressWarnings(loadGwasSumStatsFromManifest(
      manifest, genome = "hg38", ldSketch = .toyLdSketch(),
      sampleSelect = "nope")),
    "not a sample")
})

test_that("GWAS-VCF without SS errors", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "noss.vcf"), dropSs = TRUE)
  manifest <- data.frame(study = "study1", sumStatsPath = vcf,
                         stringsAsFactors = FALSE)
  expect_error(
    suppressWarnings(loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                                  ldSketch = .toyLdSketch())),
    "SS")
})

test_that("GWAS-VCF region read: tabixed honours it, plain warns", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "study1.vcf"))
  gz <- Rsamtools::bgzip(vcf, dest = paste0(vcf, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(gz, format = "vcf")
  df <- suppressWarnings(pecotmr:::.readSumStatsVcf(
    gz, "22:14560000-14565000", NULL, NULL, "lbl"))
  expect_equal(nrow(df), 2L)
  # Plain (non-indexed) VCF ignores the region with a warning.
  expect_warning(
    pecotmr:::.readSumStatsVcf(vcf, "22:14560000-14565000", NULL, NULL, "lbl"),
    "region ignored")
})

test_that("GWAS-VCF without ES/SE errors", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "noes.vcf"), dropEs = TRUE)
  manifest <- data.frame(study = "study1", sumStatsPath = vcf,
                         stringsAsFactors = FALSE)
  expect_error(
    suppressWarnings(loadGwasSumStatsFromManifest(manifest, genome = "hg38",
                                                  ldSketch = .toyLdSketch())),
    "ES and SE")
})

test_that("GWAS-VCF region on an absent chromosome yields no rows", {
  skip_if_not_installed("VariantAnnotation")
  tmp <- withr::local_tempdir()
  vcf <- .writeGwasVcf(file.path(tmp, "s.vcf"))
  gz <- Rsamtools::bgzip(vcf, dest = paste0(vcf, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(gz, format = "vcf")
  df <- suppressWarnings(
    pecotmr:::.readSumStatsVcf(gz, "99:1-100", NULL, NULL, "lbl"))
  expect_equal(nrow(df), 0L)
})

# ===========================================================================
# Coverage: remaining error/edge branches
# ===========================================================================

test_that("optional sumstats columns (BETA/SE/P/MAF) are attached", {
  tmp <- withr::local_tempdir()
  df <- .toyGwasDf(5)
  df$beta <- df$z * 0.1; df$se <- rep(0.1, 5)
  df$p <- rep(0.05, 5); df$maf <- rep(0.2, 5)
  ssPath <- .writeSumstatsTsv(df, file.path(tmp, "s.tsv"))
  obj <- loadGwasSumStatsFromManifest(
    data.frame(study = "s1", sumStatsPath = ssPath, stringsAsFactors = FALSE),
    genome = "hg38", ldSketch = .toyLdSketch())
  mc <- S4Vectors::mcols(obj$entry[[1L]])
  expect_true(all(c("BETA", "SE", "P", "MAF") %in% colnames(mc)))
})

test_that(".detectGenotypeFormat probes a PLINK2 prefix", {
  tmp <- withr::local_tempdir()
  file.create(file.path(tmp, "g.pgen"))
  expect_error(pecotmr:::.detectGenotypeFormat(file.path(tmp, "g")))
})

test_that(".resolveLdSketch passes through NULL and a handle", {
  expect_null(pecotmr:::.resolveLdSketch(NULL))
  h <- .toyLdSketch()
  expect_identical(pecotmr:::.resolveLdSketch(h), h)
})

test_that("non-constant ldSketchPath column errors", {
  tmp <- withr::local_tempdir()
  ss1 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "a.tsv"))
  ss2 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "b.tsv"))
  manifest <- data.frame(study = c("a", "b"), sumStatsPath = c(ss1, ss2),
                         ldSketchPath = c(.toyRefPrefix(), "other"),
                         stringsAsFactors = FALSE)
  expect_error(loadGwasSumStatsFromManifest(manifest, genome = "hg38"),
               "must be constant")
})

test_that("tabix region with absent chrom or empty range returns no rows", {
  tmp <- withr::local_tempdir()
  tsv <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
  gz <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".gz"), overwrite = TRUE)
  Rsamtools::indexTabix(gz, seq = 1L, start = 2L, end = 2L, comment = "#")
  expect_equal(nrow(pecotmr:::.readSumStatsText(gz, "99:1-100", NULL, "l")), 0L)
  expect_equal(nrow(pecotmr:::.readSumStatsText(gz, "22:1-2", NULL, "l")), 0L)
})

test_that("phenotype BED missing required columns errors", {
  tmp <- withr::local_tempdir()
  bad <- file.path(tmp, "bad.bed")
  readr::write_tsv(
    data.frame(`#chr` = "chr22", start = 1L, end = 2L, check.names = FALSE), bad)
  manifest <- data.frame(context = "ctx", phenotypePath = bad, study = "s",
                         genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  expect_error(loadQtlDatasetFromManifest(manifest), "missing one of")
})

test_that("covariate file with no shared samples errors", {
  tmp <- withr::local_tempdir()
  bed <- .writePhenoBed(file.path(tmp, "c.bed"), paste0("S", 1:6))
  cov <- .writeCovTsv(file.path(tmp, "cov.tsv"), paste0("Z", 1:6))
  manifest <- data.frame(context = "ctx", phenotypePath = bed,
                         covariatePath = cov, study = "s",
                         genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  expect_error(loadQtlDatasetFromManifest(manifest), "No shared samples")
})

test_that("genotypeCovariates accepts a path string", {
  tmp <- withr::local_tempdir()
  samples <- paste0("S", 1:6)
  bed <- .writePhenoBed(file.path(tmp, "c.bed"), samples)
  gcov <- .writeCovTsv(file.path(tmp, "g.tsv"), samples, npc = 2)
  manifest <- data.frame(context = "ctx", phenotypePath = bed, study = "s",
                         genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  qd <- loadQtlDatasetFromManifest(manifest, genotypeCovariates = gcov)
  expect_equal(ncol(getGenotypeCovariates(qd)), 2L)
})

test_that("QtlSumStats manifest resolves a per-row columnMapping", {
  tmp <- withr::local_tempdir()
  df <- .toyGwasDf(5)
  names(df) <- c("CHR", "POS", "RSID", "EA", "OA", "ZSCORE", "SAMPLE_N")
  ssPath <- .writeSumstatsTsv(df, file.path(tmp, "g.tsv"))
  mapPath <- file.path(tmp, "m.yaml")
  yaml::write_yaml(list(chrom = "CHR", pos = "POS", variant_id = "RSID",
                        A1 = "EA", A2 = "OA", z = "ZSCORE",
                        n_sample = "SAMPLE_N"), mapPath)
  manifest <- data.frame(study = "eqtl", context = "Blood", trait = "geneA",
                         sumStatsPath = ssPath, columnMapping = mapPath,
                         stringsAsFactors = FALSE)
  obj <- loadQtlSumStatsFromManifest(manifest, genome = "hg38",
                                     ldSketch = .toyLdSketch())
  expect_equal(length(obj$entry[[1L]]), 5L)
})

test_that("QtlDataset builder errors on conflicting per-context/study paths", {
  tmp <- withr::local_tempdir()
  samples <- paste0("S", 1:6)
  bed <- .writePhenoBed(file.path(tmp, "c.bed"), samples)
  cov1 <- .writeCovTsv(file.path(tmp, "cov1.tsv"), samples)
  cov2 <- .writeCovTsv(file.path(tmp, "cov2.tsv"), samples)
  # Two rows, same context + phenotype, conflicting covariatePath.
  m1 <- data.frame(context = c("ctx", "ctx"), phenotypePath = c(bed, bed),
                   covariatePath = c(cov1, cov2), study = "s",
                   genotypePath = .toyRefPrefix(), stringsAsFactors = FALSE)
  expect_error(loadQtlDatasetFromManifest(m1), "multiple covariatePath")
  # Conflicting genotypeCovariatePath.
  m2 <- data.frame(context = c("ctx", "ctx"), phenotypePath = c(bed, bed),
                   study = "s", genotypePath = .toyRefPrefix(),
                   genotypeCovariatePath = c(cov1, cov2), stringsAsFactors = FALSE)
  expect_error(loadQtlDatasetFromManifest(m2), "multiple genotypeCovariatePath")
})

test_that("MultiStudy builder errors on conflicting per-study genotypePath", {
  tmp <- withr::local_tempdir()
  bedA <- .writePhenoBed(file.path(tmp, "a.bed"))
  bedB <- .writePhenoBed(file.path(tmp, "b.bed"))
  man <- data.frame(study = c("study1", "study1"), context = c("c1", "c2"),
                    phenotypePath = c(bedA, bedB),
                    genotypePath = c(.toyRefPrefix(), "other"),
                    stringsAsFactors = FALSE)
  expect_error(loadMultiStudyQtlDatasetFromManifest(man),
               "exactly one genotypePath")
})
