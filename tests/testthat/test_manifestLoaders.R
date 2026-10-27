# Tests for R/manifestLoaders.R

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# The toy_ref PLINK1 panel (chr22) serves as the LD sketch. These variants
# are drawn verbatim from toy_ref.bim so LD-containment succeeds.
.toyRefPrefix <- function() {
    sub(
        "\\.bed$",
        "",
        system.file("extdata", "toy_ref.bed", package = "pecotmr")
    )
}

# Read once and shared: the panel is a value, and some forty tests take one.
.toyLdRefPanel <- readGenotypes(plink1Prefix = .toyRefPrefix())

.toyLdRef <- function() {
    .toyLdRefPanel
}

# A small GWAS sumstats data.frame over real toy_ref chr22 variants.
.toyGwasDf <- function(n = 5) {
    base <- data.frame(
        chrom = rep("22", 5),
        pos = c(14560203, 14564328, 14850625, 14870204, 14878387),
        variant_id = c(
            "rs11089128",
            "rs7288972",
            "rs11167319",
            "rs8138488",
            "rs2186521"
        ),
        A1 = c("G", "C", "G", "C", "A"),
        A2 = c("A", "T", "T", "T", "G"),
        z = c(1.2, -0.4, 2.1, 0.05, -1.8),
        n_sample = rep(10000L, 5),
        # The A1 frequency, directional, which a column named `af` declares;
        # without one every load warns that top_loci$af will be NA.
        af = c(0.12, 0.33, 0.48, 0.71, 0.05),
        stringsAsFactors = FALSE
    )
    base[seq_len(n), , drop = FALSE]
}

# Write a data.frame to a TSV whose header line starts with "#" (so a tabix
# index can treat it as a comment/header line).
.writeSumstatsTsv <- function(df, path) {
    hdr <- names(df)
    hdr[[1L]] <- paste0("#", hdr[[1L]])
    con <- file(path, "w")
    writeLines(paste(hdr, collapse = "\t"), con)
    utils::write.table(
        df,
        con,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
    close(con)
    path
}

# ---------------------------------------------------------------------------
# Column canonicalisation + scalar reconciliation
# ---------------------------------------------------------------------------

test_that(".canonManifestCols renames snake_case aliases to camelCase", {
    df <- data.frame(
        study_id = "s1",
        sumstats_path = "x.tsv",
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.canonManifestCols(
        df,
        pecotmr:::.gwasSumStatsManifestAliases,
        required = c("study", "sumStatsPath"),
        label = "GwasSumStats"
    )
    expect_true(all(c("study", "sumStatsPath") %in% names(out)))
    expect_false(any(c("study_id", "sumstats_path") %in% names(out)))
})

test_that(".canonManifestCols errors on missing required columns", {
    df <- data.frame(study = "s1", stringsAsFactors = FALSE)
    expect_error(
        pecotmr:::.canonManifestCols(
            df,
            pecotmr:::.gwasSumStatsManifestAliases,
            required = c("study", "sumStatsPath"),
            label = "GwasSumStats"
        ),
        "manifest.*Colnames must include the elements"
    )
})

test_that(".reconcileScalar resolves arg/column and flags conflicts", {
    expect_equal(pecotmr:::.reconcileScalar(NULL, "hg38", "genome"), "hg38")
    expect_equal(
        pecotmr:::.reconcileScalar(c("hg38", "hg38"), NULL, "genome"),
        "hg38"
    )
    expect_equal(
        pecotmr:::.reconcileScalar(c("hg38", "hg38"), "hg38", "genome"),
        "hg38"
    )
    expect_error(
        pecotmr:::.reconcileScalar(c("hg19", "hg38"), NULL, "genome"),
        "must be constant"
    )
    expect_error(
        pecotmr:::.reconcileScalar("hg19", "hg38", "genome"),
        "disagrees"
    )
    expect_error(
        pecotmr:::.reconcileScalar(NULL, NULL, "genome"),
        "must be provided"
    )
})

# ---------------------------------------------------------------------------
# Genotype-format detection
# ---------------------------------------------------------------------------

test_that(".readGenotypePath reads a PLINK1 prefix as a panel", {
    panel <- pecotmr:::.readGenotypePath(.toyRefPrefix())
    expect_s4_class(panel, "RangedSummarizedExperiment")
    expect_equal(pecotmr:::.ldRefFormat(panel), "plink1")
})

test_that("BCF sumstats are rejected", {
    expect_error(
        pecotmr:::.readSumStatsFile("x.bcf", NULL, NULL, NULL, NULL, "lbl"),
        "BCF sumstats are not supported"
    )
})

# ---------------------------------------------------------------------------
# LD-sketch containment
# ---------------------------------------------------------------------------

test_that(".checkLdContainment passes for overlapping variants", {
    gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(5))
    expect_silent(pecotmr:::.checkLdContainment(.toyLdRef(), gr, 0.5, "lbl"))
})

test_that(".checkLdContainment errors on a chromosome absent from the sketch", {
    df <- .toyGwasDf(3)
    df$chrom <- "1"
    gr <- pecotmr:::.dfToEntryGranges(df)
    expect_error(
        pecotmr:::.checkLdContainment(.toyLdRef(), gr, 0.5, "lbl"),
        "absent from the LD reference"
    )
})

test_that(".checkLdContainment errors on zero overlap", {
    df <- .toyGwasDf(2)
    df$pos <- c(999999001, 999999002)
    gr <- pecotmr:::.dfToEntryGranges(df)
    expect_error(
        pecotmr:::.checkLdContainment(.toyLdRef(), gr, 0.5, "lbl"),
        "none of the"
    )
})

test_that(".checkLdContainment warns on low overlap", {
    df <- .toyGwasDf(4)
    df$pos[3:4] <- c(999999001, 999999002) # 2/4 present -> 50%
    gr <- pecotmr:::.dfToEntryGranges(df)
    expect_warning(
        pecotmr:::.checkLdContainment(.toyLdRef(), gr, 0.9, "lbl"),
        "are present in the LD reference"
    )
})

# ---------------------------------------------------------------------------
# GwasSumStats loader
# ---------------------------------------------------------------------------

test_that("loadGwasSumStatsFromManifest builds from a data.frame manifest", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
    manifest <- data.frame(
        study = "study1",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    expect_true(methods::validObject(obj))
    expect_equal(as.character(obj$study), "study1")
    expect_equal(unname(GenomeInfoDb::genome(obj)), "hg38")
    expect_equal(length(obj[[1L]]), 5L)
    expect_true(all(
        c("SNP", "A1", "A2", "Z", "N") %in%
            colnames(S4Vectors::mcols(obj[[1L]]))
    ))
    expect_equal(length(qcOptions(obj)), 0L) # loaders run no QC
})

test_that("loadGwasSumStatsFromManifest reads a file and reconciles genome", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
    mfPath <- file.path(tmp, "gwas_manifest.tsv")
    readr::write_tsv(
        data.frame(
            study = "study1",
            sumStatsPath = ssPath,
            genome = "hg38",
            stringsAsFactors = FALSE
        ),
        mfPath
    )
    # genome column present; matching arg is fine.
    obj <- loadGwasSumStatsFromManifest(
        mfPath,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_equal(unname(GenomeInfoDb::genome(obj)), "hg38")
    # conflicting arg errors.
    expect_error(
        loadGwasSumStatsFromManifest(
            mfPath,
            genome = "hg19",
            ldReference = .toyLdRef()
        ),
        "disagrees"
    )
})

test_that("loadGwasSumStatsFromManifest rejects duplicate studies", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "s.tsv"))
    manifest <- data.frame(
        study = c("a", "a"),
        sumStatsPath = c(ssPath, ssPath),
        stringsAsFactors = FALSE
    )
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = .toyLdRef()
        ),
        "must be unique"
    )
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
    df <- pecotmr:::.readSumStatsText(
        gz,
        GenomicRanges::GRanges("22:14560000-14565000"),
        NULL,
        "lbl"
    )
    expect_equal(nrow(df), 2L) # first two variants only
    expect_true(all(df$pos >= 14560000 & df$pos <= 14565000))
})

test_that("plain text cannot narrow a read to the region", {
    tmp <- withr::local_tempdir()
    tsv <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
    expect_warning(
        df <- pecotmr:::.readSumStatsText(
            tsv,
            GenomicRanges::GRanges("22:14560000-14565000"),
            NULL,
            "lbl"
        ),
        "the whole file is read"
    )
    expect_equal(nrow(df), 5L)
})

# ---------------------------------------------------------------------------
# YAML column mapping
# ---------------------------------------------------------------------------

test_that("loadGwasSumStatsFromManifest resolves a YAML column mapping", {
    tmp <- withr::local_tempdir()
    # Sumstats with non-standard source column names.
    df <- .toyGwasDf(5)
    names(df) <- c(
        "CHR",
        "POS",
        "RSID",
        "EA",
        "OA",
        "ZSCORE",
        "SAMPLE_N",
        "EA_FREQ"
    )
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "study1.tsv"))
    mapPath <- file.path(tmp, "mapping.yaml")
    yaml::write_yaml(
        list(
            chrom = "CHR",
            pos = "POS",
            variant_id = "RSID",
            A1 = "EA",
            A2 = "OA",
            z = "ZSCORE",
            n_sample = "SAMPLE_N",
            af = "EA_FREQ"
        ),
        mapPath
    )
    manifest <- data.frame(
        study = "study1",
        sumStatsPath = ssPath,
        columnMapping = mapPath,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    expect_equal(length(obj[[1L]]), 5L)
    expect_equal(
        as.character(S4Vectors::mcols(obj[[1L]])$A1),
        c("G", "C", "G", "C", "A")
    )
})

# ---------------------------------------------------------------------------
# Case/control counts -> effective sample size
# ---------------------------------------------------------------------------

test_that("loadGwasSumStatsFromManifest reads per-variant N_CASE/N_CONTROL", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL # no N column ...
    df$n_case <- rep(1000L, 5) # ... only per-variant case/control counts
    df$n_control <- rep(2000L, 5)
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "cc.tsv"))
    manifest <- data.frame(
        study = "cc",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_true(all(c("N_CASE", "N_CONTROL") %in% colnames(mc)))
    expect_false("N" %in% colnames(mc)) # effective N is derived later by summaryStatsQc
    expect_equal(as.integer(mc$N_CASE), rep(1000L, 5))
})

test_that(".resolveSumstatCols requires an N field or N_CASE + N_CONTROL", {
    df <- .toyGwasDf(3)
    df$n_sample <- NULL # drop the only sample-size source
    expect_error(
        pecotmr:::.resolveSumstatCols(df, NULL, "lbl"),
        "N_CASE + N_CONTROL",
        fixed = TRUE
    )
})

test_that(".resolveSumstatCols derives the Wald z from beta/se", {
    df <- .toyGwasDf(4)
    df$beta <- df$z * 0.5 # se = 0.5 -> beta/se reproduces z
    df$se <- rep(0.5, 4)
    df$z <- NULL # no z column
    out <- pecotmr:::.resolveSumstatCols(df, NULL, "lbl")
    expect_equal(out$Z, .toyGwasDf(4)$z) # Z == beta/se
    expect_true(all(c("BETA", "SE") %in% colnames(out))) # beta/se still attached
})

test_that("an af: key is directional AF; a maf: key is directionless", {
    # Same source column, different declaration -> different meaning.
    a <- suppressWarnings(.resolveSumstatCols(
        .afm_df("effect_allele_frequency", 0.8),
        .afm_cm(af = "effect_allele_frequency"),
        "af"
    ))
    expect_equal(a$AF, 0.8)
    expect_false("MAF" %in% names(a))
    b <- suppressWarnings(.resolveSumstatCols(
        .afm_df("effect_allele_frequency", 0.8),
        .afm_cm(maf = "effect_allele_frequency"),
        "maf"
    ))
    expect_false("AF" %in% names(b))
    expect_equal(b$MAF, 0.8)
})

test_that("an undeclared af warns distinctly and yields no AF column", {
    expect_warning(
        out <- .resolveSumstatCols(
            .afm_df("maf", 0.2),
            .afm_cm(maf = "maf"),
            "undeclared"
        ),
        "no effect-allele frequency declared"
    )
    expect_false("AF" %in% names(out))
})

test_that("a declared-but-unusable af warns with a different message", {
    # Distinct from "undeclared": this is a data problem, not a mapping one.
    expect_warning(
        .resolveSumstatCols(
            .afm_df("eaf", NA_real_),
            .afm_cm(af = "eaf"),
            "missing"
        ),
        "declared but its values are all missing"
    )
})

test_that(".resolveSumstatCols keeps a supplied z over beta/se", {
    df <- .toyGwasDf(4)
    df$beta <- rep(99, 4)
    df$se <- rep(1, 4) # beta/se = 99, must be ignored
    out <- pecotmr:::.resolveSumstatCols(df, NULL, "lbl")
    expect_equal(out$Z, .toyGwasDf(4)$z) # z column wins, not 99
})

test_that(".resolveSumstatCols errors when neither z nor beta+se is present", {
    df <- .toyGwasDf(3)
    df$z <- NULL # no z, no beta/se
    expect_error(
        pecotmr:::.resolveSumstatCols(df, NULL, "lbl"),
        "needs a z field",
        fixed = TRUE
    )
})

test_that("loadGwasSumStatsFromManifest builds from beta/se with no z column", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$beta <- df$z * 2
    df$se <- rep(2, 5) # beta/se reproduces z
    z_expected <- df$z
    df$z <- NULL
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "betase.tsv"))
    manifest <- data.frame(
        study = "bs",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_true("Z" %in% colnames(mc))
    expect_equal(as.numeric(mc$Z), z_expected) # derived Wald z
})

test_that(".resolveSumstatCols allows no N when allowNoN = TRUE", {
    df <- .toyGwasDf(3)
    df$n_sample <- NULL # no per-variant N, no counts
    out <- pecotmr:::.resolveSumstatCols(df, NULL, "lbl", allowNoN = TRUE)
    expect_false("N" %in% colnames(out)) # built without N ...
    expect_true(all(c("SNP", "A1", "A2", "Z") %in% colnames(out))) # ... core kept
})

test_that("loadGwasSumStatsFromManifest builds from a study nSample scalar", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL # sumstats has no per-variant N and no counts
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "studyn.tsv"))
    manifest <- data.frame(
        study = "sn",
        sumStatsPath = ssPath,
        n_sample = 487511,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    expect_false("N" %in% colnames(S4Vectors::mcols(obj[[1L]])))
    expect_equal(as.numeric(obj$nSample), 487511) # forwarded to the slot
})

test_that("loadGwasSumStatsFromManifest builds from case/control scalars", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "studycc.tsv"))
    manifest <- data.frame(
        study = "cc",
        sumStatsPath = ssPath,
        n_case = 5000,
        n_control = 15000,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "GwasSumStats")
    expect_false("N" %in% colnames(S4Vectors::mcols(obj[[1L]])))
    expect_equal(as.numeric(obj$nCase), 5000)
    expect_equal(as.numeric(obj$nControl), 15000)
})

test_that("loadGwasSumStatsFromManifest errors when a study has no N source", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL # no per-variant N, no counts ...
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "non.tsv"))
    manifest <- data.frame(
        study = "x",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    ) # ... and no study scalar either
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = .toyLdRef()
        ),
        "needs an N field"
    )
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
    manifest <- data.frame(
        study = "eqtl",
        context = "Whole_Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    obj <- loadQtlSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "QtlSumStats")
    expect_true(methods::validObject(obj))
    expect_equal(as.character(obj$trait), "geneA")
    expect_equal(as.character(obj$context), "Whole_Blood")
})

test_that("loadQtlSumStatsFromManifest derives the Wald z from beta/se", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    z_expected <- df$z
    df$beta <- df$z * 2
    df$se <- rep(2, 5) # beta/se reproduces z
    df$z <- NULL # ... but no z column
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "geneA.tsv"))
    manifest <- data.frame(
        study = "eqtl",
        context = "Whole_Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    obj <- loadQtlSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "QtlSumStats")
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_equal(as.numeric(mc$Z), z_expected) # derived Wald z
    expect_true(all(c("BETA", "SE") %in% colnames(mc))) # beta/se still attached
})

# A cis-QTL scan's nominal output: every gene in one table, effect columns named
# bhat/sebhat, no allele columns (the alleles live in the variant id), and the
# ids written in PLINK .bim order chr:pos:A1:A2.
.toyQtlCisDf <- function() {
    # positions + alleles of the first four toy_ref panel variants, so the ids
    # are consistent with the panel they will be matched against
    pos <- c(14560203L, 14564328L, 14850625L, 14870204L)
    a1 <- c("G", "C", "G", "C")
    a2 <- c("A", "T", "T", "T")
    data.frame(
        chrom = rep("22", 8L),
        pos = rep(pos, 2L),
        molecular_trait_id = rep(c("geneA", "geneB"), each = 4L),
        variant_id = rep(paste0("chr22:", pos, ":", a1, ":", a2), 2L),
        af = rep(0.3, 8L),
        pvalue = rep(0.1, 8L),
        bhat = c(1.2, -0.4, 2.1, 0.05, 0.9, -1.1, 0.2, 0.3),
        sebhat = rep(2, 8L),
        n = rep(49L, 8L),
        stringsAsFactors = FALSE
    )
}

.qtlCisManifest <- function(path, trait = "geneA") {
    data.frame(
        study = "eqtl",
        context = "context1",
        trait = trait,
        sumStatsPath = path,
        stringsAsFactors = FALSE
    )
}

test_that("loadQtlSumStatsFromManifest reads a cis-QTL nominal table", {
    # Covers the trait filter, id-derived alleles and bhat/sebhat -> Z.
    tmp <- withr::local_tempdir()
    ss <- .writeSumstatsTsv(.toyQtlCisDf(), file.path(tmp, "cis.tsv"))
    obj <- loadQtlSumStatsFromManifest(
        .qtlCisManifest(ss),
        genome = "hg38",
        ldReference = .toyLdRef(),
        traitColumn = "molecular_trait_id",
        variantIdAlleles = "A1A2"
    )
    expect_s4_class(obj, "QtlSumStats")
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_length(obj[[1L]], 4L) # geneB's rows dropped
    expect_equal(as.character(mc$A1), c("G", "C", "G", "C")) # id field 3 -> A1
    expect_equal(as.character(mc$A2), c("A", "T", "T", "T"))
    expect_equal(as.numeric(mc$Z), c(0.6, -0.2, 1.05, 0.025)) # bhat/sebhat
    expect_equal(as.integer(mc$N), rep(49L, 4L))
})

test_that("variantIdAlleles = 'A2A1' reads the id in canonical order", {
    tmp <- withr::local_tempdir()
    ss <- .writeSumstatsTsv(.toyQtlCisDf(), file.path(tmp, "cis.tsv"))
    obj <- loadQtlSumStatsFromManifest(
        .qtlCisManifest(ss),
        genome = "hg38",
        ldReference = .toyLdRef(),
        traitColumn = "molecular_trait_id",
        variantIdAlleles = "A2A1"
    )
    mc <- S4Vectors::mcols(obj[[1L]])
    # Alleles come off the id the other way round.
    expect_equal(as.character(mc$A1), c("A", "T", "T", "T"))
    expect_equal(as.character(mc$A2), c("G", "C", "G", "C"))
})

test_that("an explicit allele column beats variantIdAlleles", {
    tmp <- withr::local_tempdir()
    df <- .toyQtlCisDf()
    df$A1 <- "T"
    df$A2 <- "C"
    ss <- .writeSumstatsTsv(df, file.path(tmp, "cis.tsv"))
    obj <- loadQtlSumStatsFromManifest(
        .qtlCisManifest(ss),
        genome = "hg38",
        ldReference = .toyLdRef(),
        traitColumn = "molecular_trait_id",
        variantIdAlleles = "A1A2"
    )
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_equal(as.character(mc$A1), rep("T", 4L))
    expect_equal(as.character(mc$A2), rep("C", 4L))
})

test_that("the trait filter and the id-allele request report what went wrong", {
    tmp <- withr::local_tempdir()
    ss <- .writeSumstatsTsv(.toyQtlCisDf(), file.path(tmp, "cis.tsv"))
    expect_error(
        loadQtlSumStatsFromManifest(
            .qtlCisManifest(ss),
            genome = "hg38",
            ldReference = .toyLdRef(),
            traitColumn = "no_such_col"
        ),
        "not a column of the sumstats file"
    )
    expect_error(
        loadQtlSumStatsFromManifest(
            .qtlCisManifest(ss, trait = "geneZ"),
            genome = "hg38",
            ldReference = .toyLdRef(),
            traitColumn = "molecular_trait_id"
        ),
        "no rows with molecular_trait_id == 'geneZ'"
    )
    # rsIDs carry no alleles, so the derivation cannot be honoured
    rs <- .toyGwasDf(5)
    rs$A1 <- NULL
    rs$A2 <- NULL
    rsPath <- .writeSumstatsTsv(rs, file.path(tmp, "rs.tsv"))
    expect_error(
        loadQtlSumStatsFromManifest(
            .qtlCisManifest(rsPath),
            genome = "hg38",
            ldReference = .toyLdRef(),
            variantIdAlleles = "A2A1"
        ),
        "carry no allele pair"
    )
})

test_that("loadQtlSumStatsFromManifest builds from a tuple nSample scalar", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL # sumstats has no per-variant N
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "geneA.tsv"))
    manifest <- data.frame(
        study = "eqtl",
        context = "Whole_Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        n_sample = 838,
        stringsAsFactors = FALSE
    )
    obj <- loadQtlSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_s4_class(obj, "QtlSumStats")
    expect_false("N" %in% colnames(S4Vectors::mcols(obj[[1L]])))
    expect_equal(as.numeric(obj$nSample), 838) # forwarded to the slot
})

test_that("loadQtlSumStatsFromManifest errors when a tuple has no N source", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$n_sample <- NULL # no per-variant N ...
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "geneA.tsv"))
    manifest <- data.frame(
        study = "eqtl",
        context = "Whole_Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    ) # ... and no nSample scalar
    expect_error(
        loadQtlSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = .toyLdRef()
        ),
        "needs an N field"
    )
})

# ---------------------------------------------------------------------------
# QtlDataset loader
# ---------------------------------------------------------------------------

# Build a tiny 2-trait phenotype BED for one context.
.writePhenoBed <- function(path, samples = paste0("S", 1:6)) {
    df <- data.frame(
        `#chr` = c("chr22", "chr22"),
        start = c(14560000L, 14870000L),
        end = c(14560001L, 14870001L),
        gene_id = c("geneA", "geneB"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    for (s in samples) {
        df[[s]] <- stats::rnorm(2)
    }
    readr::write_tsv(df, path)
    path
}

test_that("a phenotype BED's chromosomes are read in UCSC style", {
    tmp <- withr::local_tempdir()
    path <- file.path(tmp, "ncbi.bed")
    readr::write_tsv(
        data.frame(
            `#chr` = c("22", "MT"),
            start = c(10L, 20L),
            end = c(11L, 21L),
            gene_id = c("gA", "gB"),
            S1 = c(1, 2),
            check.names = FALSE
        ),
        path
    )
    rr <- pecotmr:::.readPhenotypeBed(path)$rr
    expect_identical(
        as.character(GenomicRanges::seqnames(rr)),
        c("chr22", "chrM")
    )
    expect_identical(GenomeInfoDb::seqlevelsStyle(rr), "UCSC")
})

test_that("a manifest's BED-style ranges are read in UCSC style", {
    df <- data.frame(
        chrom = c("22", "chrX"),
        start = c(0L, 10L),
        end = c(5L, 20L)
    )
    gr <- pecotmr:::.manifestRanges(df, "start", "end", "QTL", "window")
    expect_identical(
        as.character(GenomicRanges::seqnames(gr)),
        c("chr22", "chrX")
    )
    expect_identical(GenomicRanges::start(gr), c(1L, 11L))
})

test_that("loadQtlDatasetFromManifest builds a single-study dataset", {
    tmp <- withr::local_tempdir()
    bed <- .writePhenoBed(file.path(tmp, "ctxA.bed"))
    manifest <- data.frame(
        context = "ctxA",
        phenotypePath = bed,
        study = "studyX",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    qd <- loadQtlDatasetFromManifest(manifest)
    expect_s4_class(qd, "QtlDataset")
    expect_true(methods::validObject(qd))
    expect_equal(studyName(qd), "studyX")
    expect_equal(contexts(qd), "ctxA")
})

test_that("loadQtlDatasetFromManifest accepts study/genotypes as arguments", {
    tmp <- withr::local_tempdir()
    bed <- .writePhenoBed(file.path(tmp, "ctxA.bed"))
    manifest <- data.frame(
        context = "ctxA",
        phenotypePath = bed,
        stringsAsFactors = FALSE
    )
    qd <- loadQtlDatasetFromManifest(
        manifest,
        studyName = "studyY",
        genotypes = .toyLdRef()
    )
    expect_equal(studyName(qd), "studyY")
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
        stringsAsFactors = FALSE
    )
    msd <- loadMultiStudyQtlDatasetFromManifest(manifest)
    expect_s4_class(msd, "MultiStudyQtlDataset")
    expect_true(methods::validObject(msd))
    expect_equal(sort(names(msd)), c("study1", "study2"))
})

test_that("loadMultiStudyQtlDatasetFromManifest attaches a summary study", {
    tmp <- withr::local_tempdir()
    bedA <- .writePhenoBed(file.path(tmp, "a.bed"))
    qdMan <- data.frame(
        study = "study1",
        context = "ctxA",
        phenotypePath = bedA,
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "geneA.tsv"))
    ssMan <- data.frame(
        study = "eqtl2",
        context = "Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    msd <- loadMultiStudyQtlDatasetFromManifest(
        qdMan,
        sumStatsManifest = ssMan,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_true(methods::validObject(msd))
    expect_s4_class(msd@sumStats, "QtlSumStats")
    expect_equal(names(msd), "study1")
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
    expect_error(
        pecotmr:::.readManifest("/no/such/manifest.tsv"),
        "File does not exist"
    )
    expect_error(
        pecotmr:::.readManifest(42L),
        "Must be of type 'string'"
    )
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
    readr::write_tsv(
        data.frame(study = "study1", sumStatsPath = "study1.tsv"),
        mfPath
    )
    obj <- loadGwasSumStatsFromManifest(
        mfPath,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_equal(length(obj[[1L]]), 5L)
})

test_that(".reconcileScalar returns NULL for an optional unset scalar", {
    expect_null(pecotmr:::.reconcileScalar(NULL, NULL, "x", required = FALSE))
})

# ===========================================================================
# Coverage: genotype-format detection and LD-sketch resolution
# ===========================================================================

test_that(".readGenotypePath dispatches by extension", {
    expect_error(pecotmr:::.readGenotypePath("x.vcf.gz")) # vcf branch
    expect_error(pecotmr:::.readGenotypePath("x.pgen")) # plink2 branch
    expect_error(
        pecotmr:::.readGenotypePath("x.unknown"),
        "manifest: .* is neither a genotype file nor a PLINK prefix"
    )
    panel <- pecotmr:::.readGenotypePath(paste0(.toyRefPrefix(), ".bed"))
    expect_equal(pecotmr:::.ldRefFormat(panel), "plink1")
})

test_that(".readGenotypePath reads a VCF as a panel", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), c("s1", "s2", "s3"))
    panel <- suppressWarnings(pecotmr:::.readGenotypePath(vcf))
    expect_equal(pecotmr:::.ldRefFormat(panel), "vcf")
})

test_that(".resolveLdRef accepts a genoMeta vector or path, rejects junk", {
    sharded <- pecotmr:::.resolveLdRef(c("22" = .toyRefPrefix()))
    expect_s4_class(sharded, "RangedSummarizedExperiment")
    expect_true("chr22" %in% names(pecotmr:::.ldRefChromPaths(sharded)))
    expect_s4_class(
        pecotmr:::.resolveLdRef(.toyRefPrefix()),
        "RangedSummarizedExperiment"
    )
    expect_error(pecotmr:::.resolveLdRef("nonexistent_meta.tsv"))
    expect_error(pecotmr:::.resolveLdRef(42L), "must be a genotype panel")
})

test_that(".resolveLdRef takes the panel readGenotypes returns", {
    # readGenotypes() is the only public way to open genotypes now, so the
    # loaders have to accept its return value wherever they accepted a
    # handle -- as it stands, since a panel need not read through a handle,
    # and one the caller narrowed stays narrowed.
    panel <- readGenotypes(plink1Prefix = .toyRefPrefix())
    expect_identical(pecotmr:::.resolveLdRef(panel), panel)
    narrowed <- panel[1:5, 1:10]
    expect_identical(pecotmr:::.resolveLdRef(narrowed), narrowed)
})

test_that(".materializeLdRef and .resolveGenotypes take a panel too", {
    panel <- readGenotypes(plink1Prefix = .toyRefPrefix())
    expect_identical(pecotmr:::.materializeLdRef(panel, "22"), panel)
    expect_identical(
        pecotmr:::.resolveGenotypes(NULL, "s1", NULL, panel),
        panel
    )
})

test_that(".entriesChroms collects the chromosomes across entries", {
    a <- pecotmr:::.dfToEntryGranges(.toyGwasDf(3)) # chr22
    b <- .toyGwasDf(2)
    b$chrom <- "1"
    b <- pecotmr:::.dfToEntryGranges(b)
    expect_setequal(pecotmr:::.entriesChroms(list(a, b)), c("chr22", "chr1"))
    expect_equal(pecotmr:::.entriesChroms(list()), character(0))
    expect_equal(pecotmr:::.entriesChroms(list(NULL)), character(0))
})

test_that(".materializeLdRef passes a handle through, restricts shards", {
    h <- .toyLdRef()
    # already in memory
    expect_identical(pecotmr:::.materializeLdRef(h, "22"), h)
    expect_null(pecotmr:::.materializeLdRef(NULL, "22"))
    # A chrom-sharded spec with a bogus chr21 payload: restricting to chr22 must
    # skip (never open) the chr21 shard; requesting chr21 too fails on it.
    spec <- c("22" = .toyRefPrefix(), "21" = "/no/such/chr21/prefix")
    restricted <- pecotmr:::.materializeLdRef(spec, "22")
    expect_s4_class(restricted, "RangedSummarizedExperiment")
    expect_equal(names(pecotmr:::.ldRefChromPaths(restricted)), "chr22")
    expect_error(pecotmr:::.materializeLdRef(spec, c("22", "21")))
})

test_that("loadGwasSumStatsFromManifest reads only the needed shards", {
    skip_if_not_installed("PlinkMatrix")
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "study1.tsv"))
    manifest <- data.frame(
        study = "study1",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    # Sharded ldReference: chr22 real (toy_ref), chr21 a bogus path. The
    # sumstats are all chr22, so the loader must skip the chr21 shard --
    # opening it would fail.
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = c("22" = .toyRefPrefix(), "21" = "/no/such/chr21/prefix")
    )
    expect_s4_class(obj, "GwasSumStats")
    expect_true(methods::validObject(obj))
})

test_that("ldReference resolves from ldReferencePath; conflicts error", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
    manifest <- data.frame(
        study = "s1",
        sumStatsPath = ssPath,
        ldReferencePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(manifest, genome = "hg38")
    expect_s4_class(obj, "GwasSumStats")
    # Conflicting arg vs column.
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = "/some/other/prefix"
        ),
        "disagrees"
    )
    # Neither arg nor column.
    bare <- data.frame(
        study = "s1",
        sumStatsPath = ssPath,
        stringsAsFactors = FALSE
    )
    expect_error(
        loadGwasSumStatsFromManifest(bare, genome = "hg38"),
        "must be provided as an argument or an .ldReferencePath"
    )
})

test_that("LD containment uses a sharded sketch's chromPaths", {
    skip_if_not_installed("PlinkMatrix")
    sharded <- readGenotypes(genoMeta = c("22" = .toyRefPrefix()))
    gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(5))
    expect_silent(pecotmr:::.checkLdContainment(sharded, gr, 0.5, "lbl"))
})

test_that(".checkLdContainment short-circuits on an empty entry", {
    expect_null(pecotmr:::.checkLdContainment(
        .toyLdRef(),
        GenomicRanges::GRanges(),
        0.5,
        "lbl"
    ))
})

test_that(".checkLdContainment warns when the sketch has no variant metadata", {
    h <- methods::new(
        "GenotypeHandle",
        path = "meta",
        format = "plink2",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL,
        chromPaths = c(chr22 = "/p")
    )
    gr <- pecotmr:::.dfToEntryGranges(.toyGwasDf(3))
    expect_warning(
        pecotmr:::.checkLdContainment(h, gr, 0.5, "lbl"),
        "no variant metadata"
    )
})

# ===========================================================================
# Coverage: region coercion and column-mapping errors
# ===========================================================================

test_that(".readColumnMapping errors on bad inputs", {
    tmp <- withr::local_tempdir()
    expect_error(
        pecotmr:::.readColumnMapping("/no/map.yaml"),
        "File does not exist"
    )
    bad <- file.path(tmp, "bad.yaml")
    yaml::write_yaml(list("a", "b"), bad) # unnamed sequence
    expect_error(pecotmr:::.readColumnMapping(bad), "standardName: sourceName")
    expect_error(
        pecotmr:::.readColumnMapping(42L),
        "named list/vector or a path"
    )
})

test_that(".resolveSumstatCols errors on bad mapping and missing fields", {
    df <- .toyGwasDf(3)
    names(df) <- c("chrom", "pos", "variant_id", "A1", "A2", "z", "n_sample")
    expect_error(
        pecotmr:::.resolveSumstatCols(df, list(chrom = "NOPE"), "lbl"),
        "is not a column"
    )
    expect_error(
        pecotmr:::.resolveSumstatCols(df[, c("chrom", "pos")], NULL, "lbl"),
        "missing required field"
    )
})

# ===========================================================================
# Coverage: nCase/nControl/varY columns
# ===========================================================================

test_that("loadGwasSumStatsFromManifest attaches nCase/nControl/varY", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
    manifest <- data.frame(
        study = "cc",
        sumStatsPath = ssPath,
        nCase = 1200,
        nControl = 3400,
        varY = 0.19,
        stringsAsFactors = FALSE
    )
    obj <- loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_equal(as.numeric(obj$nCase), 1200)
    expect_equal(as.numeric(obj$nControl), 3400)
    expect_equal(as.numeric(obj$varY), 0.19)
})

test_that("loadQtlSumStatsFromManifest attaches varY and honours region", {
    tmp <- withr::local_tempdir()
    ssPath <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "geneA.tsv"))
    manifest <- data.frame(
        study = "eqtl",
        context = "Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        varY = 1.7,
        stringsAsFactors = FALSE
    )
    obj <- suppressWarnings(loadQtlSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef(),
        which = GenomicRanges::GRanges("22:1-50000000")
    ))
    expect_equal(as.numeric(obj$varY), 1.7)
})

# ===========================================================================
# Coverage: QtlDataset covariates and error branches
# ===========================================================================

.writeCovTsv <- function(path, samples, npc = 3, transpose = FALSE) {
    set.seed(1)
    m <- matrix(
        stats::rnorm(length(samples) * npc),
        nrow = length(samples),
        dimnames = list(samples, paste0("PC", seq_len(npc)))
    )
    if (transpose) {
        out <- data.frame(covariate = colnames(m), t(m), check.names = FALSE)
    } else {
        out <- data.frame(sample = rownames(m), m, check.names = FALSE)
    }
    readr::write_tsv(out, path)
    path
}

test_that("loadQtlDatasetFromManifest attaches context/genotype covariates", {
    tmp <- withr::local_tempdir()
    samples <- paste0("S", 1:6)
    bed <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
    pcov <- .writeCovTsv(file.path(tmp, "pcov.tsv"), samples, npc = 3)
    gcov <- .writeCovTsv(file.path(tmp, "gcov.tsv"), samples, npc = 2)
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bed,
        covariatePath = pcov,
        study = "s",
        genotypePath = .toyRefPrefix(),
        genotypeCovariatePath = gcov,
        stringsAsFactors = FALSE
    )
    qd <- loadQtlDatasetFromManifest(manifest)
    expect_equal(ncol(phenotypeCovariates(qd, "ctx")[["ctx"]]), 3L)
    expect_equal(ncol(genotypeCovariates(qd)), 2L)
})

test_that("loadQtlDatasetFromManifest reads transposed (QTLtools) covariates", {
    tmp <- withr::local_tempdir()
    samples <- paste0("S", 1:6)
    bed <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
    pcov <- .writeCovTsv(
        file.path(tmp, "pcovT.tsv"),
        samples,
        npc = 4,
        transpose = TRUE
    )
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bed,
        covariatePath = pcov,
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    qd <- loadQtlDatasetFromManifest(manifest, transposeCovariates = TRUE)
    expect_equal(ncol(phenotypeCovariates(qd, "ctx")[["ctx"]]), 4L)
})

test_that("loadQtlDatasetFromManifest accepts a prefix and matrix covariates", {
    tmp <- withr::local_tempdir()
    samples <- paste0("S", 1:6)
    bed <- .writePhenoBed(file.path(tmp, "ctx.bed"), samples)
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bed,
        study = "s",
        stringsAsFactors = FALSE
    )
    gmat <- matrix(0, nrow = 6, ncol = 2, dimnames = list(samples, c("A", "B")))
    qd <- loadQtlDatasetFromManifest(
        manifest,
        genotypes = .toyRefPrefix(),
        genotypeCovariates = gmat
    )
    expect_equal(ncol(genotypeCovariates(qd)), 2L)
})

test_that("QtlDataset builder errors on inconsistent per-context paths", {
    tmp <- withr::local_tempdir()
    bed1 <- .writePhenoBed(file.path(tmp, "c1.bed"))
    bed2 <- .writePhenoBed(file.path(tmp, "c2.bed"))
    # Same context, two different phenotype paths.
    manifest <- data.frame(
        context = c("ctx", "ctx"),
        phenotypePath = c(bed1, bed2),
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    expect_error(
        loadQtlDatasetFromManifest(manifest),
        "exactly one phenotypePath"
    )
})

# ===========================================================================
# Coverage: GWAS-VCF reading
# ===========================================================================

# Write a minimal GWAS-VCF (ES/SE/LP/SS/EAF FORMAT) over toy_ref chr22 sites.
# A1 (effect) = ALT, A2 = REF, matching the toy_ref alleles.
.writeGwasVcf <- function(
    path,
    samples = "study1",
    n = 5,
    dropSs = FALSE,
    dropEs = FALSE,
    reference = NULL
) {
    pos <- c(14560203, 14564328, 14850625, 14870204, 14878387)[seq_len(n)]
    ids <- c(
        "rs11089128",
        "rs7288972",
        "rs11167319",
        "rs8138488",
        "rs2186521"
    )[seq_len(n)]
    ref <- c("A", "T", "T", "T", "G")[seq_len(n)]
    alt <- c("G", "C", "G", "C", "A")[seq_len(n)]
    es <- c(0.12, -0.04, 0.21, 0.005, -0.18)[seq_len(n)]
    se <- rep(0.1, n)
    lp <- c(1.3, 0.2, 2.1, 0.05, 1.0)[seq_len(n)]
    eaf <- c(0.3, 0.4, 0.2, 0.5, 0.1)[seq_len(n)]
    # Build the FORMAT tag list and matching per-record cells dynamically.
    tags <- c(if (!dropEs) c("ES", "SE"), "LP", if (!dropSs) "SS", "EAF")
    hdrFor <- c(
        ES = '##FORMAT=<ID=ES,Number=A,Type=Float,Description="effect">',
        SE = '##FORMAT=<ID=SE,Number=A,Type=Float,Description="se">',
        LP = '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10p">',
        SS = '##FORMAT=<ID=SS,Number=A,Type=Float,Description="n">',
        EAF = '##FORMAT=<ID=EAF,Number=A,Type=Float,Description="eaf">'
    )
    meta <- c(
        "##fileformat=VCFv4.2",
        if (!is.null(reference)) str_c("##reference=", reference),
        unname(hdrFor[tags]),
        "##contig=<ID=22>",
        paste(
            c(
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                samples
            ),
            collapse = "\t"
        )
    )
    rows <- map_chr(
        seq_len(n),
        function(i) {
            vals <- c(
                ES = sprintf("%.4f", es[i]),
                SE = sprintf("%.4f", se[i]),
                LP = sprintf("%.4f", lp[i]),
                SS = "10000",
                EAF = sprintf("%.4f", eaf[i])
            )
            cell <- paste(vals[tags], collapse = ":")
            geno <- paste(rep(cell, length(samples)), collapse = "\t")
            paste(
                c(
                    "22",
                    pos[i],
                    ids[i],
                    ref[i],
                    alt[i],
                    ".",
                    "PASS",
                    ".",
                    paste(tags, collapse = ":"),
                    geno
                ),
                collapse = "\t"
            )
        }
    )
    writeLines(c(meta, rows), path)
    path
}

test_that("loadGwasSumStatsFromManifest reads a GWAS-VCF (A1=ALT, Z=ES/SE)", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    vcf <- .writeGwasVcf(file.path(tmp, "study1.vcf"))
    manifest <- data.frame(
        study = "study1",
        sumStatsPath = vcf,
        stringsAsFactors = FALSE
    )
    obj <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    ))
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_equal(as.character(mc$A1), c("G", "C", "G", "C", "A")) # ALT
    expect_equal(as.character(mc$A2), c("A", "T", "T", "T", "G")) # REF
    expect_equal(
        mc$Z,
        c(0.12, -0.04, 0.21, 0.005, -0.18) / 0.1,
        tolerance = 1e-3
    )
    expect_equal(mc$N, rep(10000, 5))
})

test_that("a GWAS-VCF's recorded build stands in for a missing genome", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    vcf <- .writeGwasVcf(file.path(tmp, "study1.vcf"), reference = "hg38")
    manifest <- data.frame(study = "study1", sumStatsPath = vcf)
    obj <- suppressWarnings(
        loadGwasSumStatsFromManifest(manifest, ldReference = .toyLdRef())
    )
    expect_equal(unique(unname(GenomeInfoDb::genome(obj))), "hg38")
    # A build that contradicts the file's is refused, not taken.
    expect_error(
        suppressWarnings(loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg19",
            ldReference = .toyLdRef()
        )),
        "`genome` is 'hg19', but the GWAS-VCF files record 'hg38'"
    )
})

test_that("the given build and the VCFs' recorded ones reconcile", {
    f <- pecotmr:::.manifestBuild
    expect_equal(f("hg38", character(0), "f"), "hg38")
    expect_equal(f("hg38", "hg38", "f"), "hg38")
    expect_error(
        f(NULL, c("hg19", "hg38"), "f"),
        "record different genome builds \\(hg19, hg38\\)"
    )
    expect_error(f(NULL, character(0), "f"), "`genome` must be provided")
})

test_that("only a manifest's VCF paths are read for a build", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    df <- data.frame(
        sumStatsPath = c(
            .writeGwasVcf(file.path(tmp, "a.vcf"), reference = "hg38"),
            .writeGwasVcf(file.path(tmp, "b.vcf")),
            file.path(tmp, "c.tsv")
        )
    )
    expect_equal(pecotmr:::.manifestVcfBuilds(df, tmp), "hg38")
    expect_equal(
        pecotmr:::.manifestVcfBuilds(df[3L, , drop = FALSE], tmp),
        character(0)
    )
})

test_that("GWAS-VCF sample selection is enforced for multi-sample files", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    vcf <- .writeGwasVcf(
        file.path(tmp, "multi.vcf"),
        samples = c("study1", "study2")
    )
    manifest <- data.frame(
        study = "study1",
        sumStatsPath = vcf,
        stringsAsFactors = FALSE
    )
    expect_error(
        suppressWarnings(loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = .toyLdRef()
        )),
        "pass .sampleSelect"
    )
    obj <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef(),
        sampleSelect = "study2"
    ))
    expect_equal(length(obj[[1L]]), 5L)
    expect_error(
        suppressWarnings(loadGwasSumStatsFromManifest(
            manifest,
            genome = "hg38",
            ldReference = .toyLdRef(),
            sampleSelect = "nope"
        )),
        "not a sample"
    )
})

test_that("GWAS-VCF without SS, or without ES/SE, errors", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    # Each case drops fields the reader needs; the error names them.
    cases <- list(SS = list(dropSs = TRUE), `ES and SE` = list(dropEs = TRUE))
    for (field in names(cases)) {
        vcf <- exec(
            .writeGwasVcf,
            tempfile(tmpdir = tmp, fileext = ".vcf"),
            !!!cases[[field]]
        )
        manifest <- data.frame(
            study = "study1",
            sumStatsPath = vcf,
            stringsAsFactors = FALSE
        )
        expect_error(
            suppressWarnings(loadGwasSumStatsFromManifest(
                manifest,
                genome = "hg38",
                ldReference = .toyLdRef()
            )),
            field,
            info = field
        )
    }
})

# Drop the ES and SE subfields from one record (or leave a header line as-is).
# @noRd
.ezDropEsSe <- function(line) {
    if (!str_starts(line, "chr1")) {
        return(line)
    }
    fields <- str_split(line, "\t")[[1L]]
    keys <- str_split(fields[[9L]], ":")[[1L]]
    values <- str_split(fields[[10L]], ":")[[1L]]
    keep <- !is_in(keys, c("ES", "SE"))
    fields[[9L]] <- str_flatten(keys[keep], ":")
    fields[[10L]] <- str_flatten(values[keep], ":")
    str_flatten(fields, "\t")
}

test_that("a z-only GWAS-VCF reads without ES and SE", {
    # Reading EZ was missing altogether, so a file whose ES/SE were derived
    # upstream from z either lost that fact or, with no ES/SE at all, was
    # refused outright.
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(.ezss(TRUE), out)
    lines <- readLines(out)
    lines <- lines[!str_detect(lines, "^##FORMAT=<ID=(ES|SE),")]
    lines <- map_chr(lines, .ezDropEsSe)
    stripped <- withr::local_tempfile(fileext = ".vcf")
    writeLines(lines, stripped)
    back <- .readSumStatsVcf(stripped, NULL, NULL, NULL, "x")
    expect_equal(back$Z, c(2.5, -1.2))
    expect_true(all(back$z_supplied))
    expect_false(is_in("BETA", names(back)))
})

test_that("GWAS-VCF region read: tabixed honours it, plain warns", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    vcf <- .writeGwasVcf(file.path(tmp, "study1.vcf"))
    gz <- Rsamtools::bgzip(vcf, dest = paste0(vcf, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "vcf")
    df <- suppressWarnings(pecotmr:::.readSumStatsVcf(
        gz,
        GenomicRanges::GRanges("22:14560000-14565000"),
        NULL,
        NULL,
        "lbl"
    ))
    expect_equal(nrow(df), 2L)
    # A plain (non-indexed) VCF is read whole, with a warning.
    expect_warning(
        pecotmr:::.readSumStatsVcf(
            vcf,
            GenomicRanges::GRanges("22:14560000-14565000"),
            NULL,
            NULL,
            "lbl"
        ),
        "the whole file is read"
    )
})

test_that("GWAS-VCF region on an absent chromosome yields no rows", {
    skip_if_not_installed("VariantAnnotation")
    tmp <- withr::local_tempdir()
    vcf <- .writeGwasVcf(file.path(tmp, "s.vcf"))
    gz <- Rsamtools::bgzip(vcf, dest = paste0(vcf, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "vcf")
    df <- suppressWarnings(
        pecotmr:::.readSumStatsVcf(
            gz,
            GenomicRanges::GRanges("99:1-100"),
            NULL,
            NULL,
            "lbl"
        )
    )
    expect_equal(nrow(df), 0L)
})

# ===========================================================================
# Coverage: remaining error/edge branches
# ===========================================================================

test_that("optional sumstats columns (BETA/SE/P/MAF) are attached", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    df$beta <- df$z * 0.1
    df$se <- rep(0.1, 5)
    df$p <- rep(0.05, 5)
    df$maf <- rep(0.2, 5)
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "s.tsv"))
    obj <- loadGwasSumStatsFromManifest(
        data.frame(
            study = "s1",
            sumStatsPath = ssPath,
            stringsAsFactors = FALSE
        ),
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    mc <- S4Vectors::mcols(obj[[1L]])
    expect_true(all(c("BETA", "SE", "P", "MAF") %in% colnames(mc)))
})

test_that(".readGenotypePath probes a PLINK2 prefix", {
    tmp <- withr::local_tempdir()
    file.create(file.path(tmp, "g.pgen"))
    expect_error(pecotmr:::.readGenotypePath(file.path(tmp, "g")))
})

test_that(".resolveLdRef passes through NULL and a handle", {
    expect_null(pecotmr:::.resolveLdRef(NULL))
    h <- .toyLdRef()
    expect_identical(pecotmr:::.resolveLdRef(h), h)
})

# The toy panel's 200-variant copy keyed by canonical ids: the same samples
# and the toy GWAS variants, read from different files.
.toyCanonicalPrefix <- function() {
    sub(
        "\\.bed$",
        "",
        system.file("extdata", "toy_canonical.bed", package = "pecotmr")
    )
}

test_that("each study's ldReferencePath gives it its own LD reference", {
    skip_if_not_installed("PlinkMatrix")
    tmp <- withr::local_tempdir()
    ss1 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "a.tsv"))
    ss2 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "b.tsv"))
    ss3 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "c.tsv"))
    manifest <- data.frame(
        study = c("a", "b", "c"),
        sumStatsPath = c(ss1, ss2, ss3),
        ldReferencePath = c(
            .toyRefPrefix(),
            .toyCanonicalPrefix(),
            .toyRefPrefix()
        ),
        stringsAsFactors = FALSE
    )
    obj <- suppressWarnings(
        loadGwasSumStatsFromManifest(manifest, genome = "hg38")
    )
    # Two distinct panels, each read once and stored once; studies a and c
    # share the first.
    refs <- ldReference(obj)
    expect_length(refs, 2L)
    keys <- obj$ldReference[match(c("a", "b", "c"), obj$study)]
    expect_identical(keys[[1L]], keys[[3L]])
    expect_false(identical(keys[[1L]], keys[[2L]]))
    expect_match(
        pecotmr:::.ldRefSource(refs[[keys[[2L]]]])$path,
        "toy_canonical"
    )
})

test_that("a manifest row with no LD reference errors", {
    tmp <- withr::local_tempdir()
    ss1 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "a.tsv"))
    ss2 <- .writeSumstatsTsv(.toyGwasDf(3), file.path(tmp, "b.tsv"))
    manifest <- data.frame(
        study = c("a", "b"),
        sumStatsPath = c(ss1, ss2),
        ldReferencePath = c(.toyRefPrefix(), ""),
        stringsAsFactors = FALSE
    )
    expect_error(
        loadGwasSumStatsFromManifest(manifest, genome = "hg38"),
        "row\\(s\\) 2 have no `ldReferencePath`"
    )
})

test_that("tabix region with absent chrom or empty range returns no rows", {
    tmp <- withr::local_tempdir()
    tsv <- .writeSumstatsTsv(.toyGwasDf(5), file.path(tmp, "s.tsv"))
    gz <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(gz, seq = 1L, start = 2L, end = 2L, comment = "#")
    expect_equal(
        nrow(pecotmr:::.readSumStatsText(
            gz,
            GenomicRanges::GRanges("99:1-100"),
            NULL,
            "l"
        )),
        0L
    )
    expect_equal(
        nrow(pecotmr:::.readSumStatsText(
            gz,
            GenomicRanges::GRanges("22:1-2"),
            NULL,
            "l"
        )),
        0L
    )
})

test_that("phenotype BED missing required columns errors", {
    tmp <- withr::local_tempdir()
    bad <- file.path(tmp, "bad.bed")
    readr::write_tsv(
        data.frame(`#chr` = "chr22", start = 1L, end = 2L, check.names = FALSE),
        bad
    )
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bad,
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    expect_error(loadQtlDatasetFromManifest(manifest), "missing one of")
})

test_that("covariate file with no shared samples errors", {
    tmp <- withr::local_tempdir()
    bed <- .writePhenoBed(file.path(tmp, "c.bed"), paste0("S", 1:6))
    cov <- .writeCovTsv(file.path(tmp, "cov.tsv"), paste0("Z", 1:6))
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bed,
        covariatePath = cov,
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    expect_error(loadQtlDatasetFromManifest(manifest), "No shared samples")
})

test_that("genotypeCovariates accepts a path string", {
    tmp <- withr::local_tempdir()
    samples <- paste0("S", 1:6)
    bed <- .writePhenoBed(file.path(tmp, "c.bed"), samples)
    gcov <- .writeCovTsv(file.path(tmp, "g.tsv"), samples, npc = 2)
    manifest <- data.frame(
        context = "ctx",
        phenotypePath = bed,
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    qd <- loadQtlDatasetFromManifest(manifest, genotypeCovariates = gcov)
    expect_equal(ncol(genotypeCovariates(qd)), 2L)
})

test_that("QtlSumStats manifest resolves a per-row columnMapping", {
    tmp <- withr::local_tempdir()
    df <- .toyGwasDf(5)
    names(df) <- c(
        "CHR",
        "POS",
        "RSID",
        "EA",
        "OA",
        "ZSCORE",
        "SAMPLE_N",
        "EA_FREQ"
    )
    ssPath <- .writeSumstatsTsv(df, file.path(tmp, "g.tsv"))
    mapPath <- file.path(tmp, "m.yaml")
    yaml::write_yaml(
        list(
            chrom = "CHR",
            pos = "POS",
            variant_id = "RSID",
            A1 = "EA",
            A2 = "OA",
            z = "ZSCORE",
            n_sample = "SAMPLE_N",
            af = "EA_FREQ"
        ),
        mapPath
    )
    manifest <- data.frame(
        study = "eqtl",
        context = "Blood",
        trait = "geneA",
        sumStatsPath = ssPath,
        columnMapping = mapPath,
        stringsAsFactors = FALSE
    )
    obj <- loadQtlSumStatsFromManifest(
        manifest,
        genome = "hg38",
        ldReference = .toyLdRef()
    )
    expect_equal(length(obj[[1L]]), 5L)
})

test_that("QtlDataset builder errors on conflicting per-context/study paths", {
    tmp <- withr::local_tempdir()
    samples <- paste0("S", 1:6)
    bed <- .writePhenoBed(file.path(tmp, "c.bed"), samples)
    cov1 <- .writeCovTsv(file.path(tmp, "cov1.tsv"), samples)
    cov2 <- .writeCovTsv(file.path(tmp, "cov2.tsv"), samples)
    # Two rows, same context + phenotype, conflicting covariatePath.
    m1 <- data.frame(
        context = c("ctx", "ctx"),
        phenotypePath = c(bed, bed),
        covariatePath = c(cov1, cov2),
        study = "s",
        genotypePath = .toyRefPrefix(),
        stringsAsFactors = FALSE
    )
    expect_error(loadQtlDatasetFromManifest(m1), "multiple covariatePath")
    # Conflicting genotypeCovariatePath.
    m2 <- data.frame(
        context = c("ctx", "ctx"),
        phenotypePath = c(bed, bed),
        study = "s",
        genotypePath = .toyRefPrefix(),
        genotypeCovariatePath = c(cov1, cov2),
        stringsAsFactors = FALSE
    )
    expect_error(
        loadQtlDatasetFromManifest(m2),
        "multiple genotypeCovariatePath"
    )
})

test_that("MultiStudy builder errors on conflicting per-study genotypePath", {
    tmp <- withr::local_tempdir()
    bedA <- .writePhenoBed(file.path(tmp, "a.bed"))
    bedB <- .writePhenoBed(file.path(tmp, "b.bed"))
    man <- data.frame(
        study = c("study1", "study1"),
        context = c("c1", "c2"),
        phenotypePath = c(bedA, bedB),
        genotypePath = c(.toyRefPrefix(), "other"),
        stringsAsFactors = FALSE
    )
    expect_error(
        loadMultiStudyQtlDatasetFromManifest(man),
        "exactly one genotypePath"
    )
})


test_that("loadGwasSumStatsFromManifest splits by LD block when asked", {
    # Without a block manifest a single-chromosome file is ONE element; the
    # manifest is what gives cTWAS the per-block granularity its EM needs.
    gwasTsv <- system.file(
        "extdata",
        "manifests",
        "protocol_example.twas.gwas_sumstats.chr22.tsv.gz",
        package = "pecotmr"
    )
    skip_if(gwasTsv == "", "bundled GWAS manifest not available")
    ldStem <- file.path(
        system.file("extdata", "ld_reference", "chr22", package = "pecotmr"),
        "protocol_example.LD.chr22"
    )
    mf <- data.frame(study = "g1", sumStatsPath = gwasTsv)
    region <- GenomicRanges::GRanges("chr22:10000000-19000000")

    plain <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest = mf,
        genome = "hg38",
        ldReference = ldStem,
        which = region
    ))
    expect_equal(length(plain), 1L)
    # The window names the row it gives, as a window of a split would.
    expect_equal(plain$blockId, "chr22_10000000_19000000")
    expect_true(plain$regionGiven)

    blocks <- GenomicRanges::GRanges(
        "chr22",
        IRanges::IRanges(c(10000000L, 15000000L), c(14999999L, 19000000L))
    )
    names(blocks) <- c("blk1", "blk2")
    split <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest = mf,
        genome = "hg38",
        ldReference = ldStem,
        which = blocks
    ))
    expect_equal(length(split), 2L)
    expect_equal(split$blockId, c("blk1", "blk2"))
    expect_equal(split$study, c("g1", "g1"))
    # Same variants either way, just distributed across the blocks.
    expect_equal(sum(lengths(split)), sum(lengths(plain)))
    # The same blocks read from a BED file, whose start is 0-based.
    bedPath <- withr::local_tempfile(fileext = ".bed")
    writeLines(
        c("chr22\t9999999\t14999999\tblk1", "chr22\t14999999\t19000000\tblk2"),
        bedPath
    )
    fromBed <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest = mf,
        genome = "hg38",
        ldReference = ldStem,
        which = bedPath
    ))
    expect_equal(fromBed$blockId, c("blk1", "blk2"))
    expect_identical(
        as.character(fromBed$region),
        as.character(split$region)
    )
    # Given in R, a window is a GRanges; a string is taken for a path.
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest = mf,
            genome = "hg38",
            ldReference = ldStem,
            which = "chr22:10000000-19000000"
        ),
        "no file exists at 'chr22:10000000-19000000'"
    )
})


test_that("containment checking is skipped when there is no LD sketch", {
    # Nothing to contain the variants in, so the check is a no-op rather than
    # an error about a missing panel.
    expect_null(
        pecotmr:::.qtlSumStatsCheckContainment(NULL, NULL, NULL, 0.5)
    )
})

test_that(".deriveAllelesFromVariantId passes through with no id column", {
    df <- data.frame(foo = 1, bar = 2)
    # Neither an allele column nor a resolvable variant_id: the missing-field
    # error belongs to the caller, so nothing is derived or raised here.
    expect_identical(
        pecotmr:::.deriveAllelesFromVariantId(
            df,
            order = c("A1", "A2"),
            mapping = NULL,
            label = "lab"
        ),
        df
    )
})

test_that(".resolveGenoCov accepts a single genotypeCovariatePath", {
    dir <- withr::local_tempdir()
    covPath <- file.path(dir, "cov.tsv")
    cov <- data.frame(
        id = c("c1", "c2"),
        s1 = c(0.1, 0.3),
        s2 = c(0.2, 0.4),
        stringsAsFactors = FALSE
    )
    readr::write_tsv(cov, covPath, progress = FALSE)
    rows <- data.frame(
        genotypeCovariatePath = c("cov.tsv", "cov.tsv", NA),
        stringsAsFactors = FALSE
    )
    # Repeated identical paths collapse to one; only a genuine conflict errors.
    out <- pecotmr:::.resolveGenoCov(rows, "s1", dir, FALSE, NULL)
    expect_true(is.matrix(out))
    expect_equal(nrow(out), 2L)
    conflicting <- data.frame(
        genotypeCovariatePath = c("a.tsv", "b.tsv"),
        stringsAsFactors = FALSE
    )
    expect_error(
        pecotmr:::.resolveGenoCov(conflicting, "s1", dir, FALSE, NULL),
        "references multiple genotypeCovariatePath values"
    )
    # Blank / NA paths leave nothing to resolve, which yields the empty
    # covariate matrix rather than an error.
    none <- data.frame(
        genotypeCovariatePath = c(NA, ""),
        stringsAsFactors = FALSE
    )
    empty <- pecotmr:::.resolveGenoCov(none, "s1", dir, FALSE, NULL)
    expect_true(is.matrix(empty))
    expect_equal(dim(empty), c(0L, 0L))
})

test_that(".readTabixRegion returns a bare tibble for a headerless file", {
    skip_if_not_installed("Rsamtools")
    dir <- withr::local_tempdir()
    plain <- file.path(dir, "noheader.tsv")
    writeLines(c("chr1\t100\t200\tA", "chr1\t300\t400\tB"), plain)
    bgz <- Rsamtools::bgzip(
        plain,
        file.path(dir, "noheader.tsv.bgz"),
        overwrite = TRUE
    )
    Rsamtools::indexTabix(bgz, seq = 1L, start = 2L, end = 3L)
    # With no "#" header line there are no column names to build an empty
    # frame from, so the miss returns a column-less tibble.
    expect_length(Rsamtools::headerTabix(Rsamtools::TabixFile(bgz))$header, 0L)
    out <- pecotmr:::.readTabixRegion(bgz, GenomicRanges::GRanges("chr9:1-100"))
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 0L)
    expect_equal(ncol(out), 0L)
})

test_that("manifest loaders guard their arguments; QC is not one of them", {
    mf <- system.file("extdata", "manifests", package = "pecotmr")
    skip_if(mf == "", "manifest fixtures unavailable")
    # Loading applies no QC: the filter is a pipeline setting, and checks its
    # own values when it is built.
    expect_error(
        loadQtlDatasetFromManifest(
            data.frame(),
            genotypeFilterParam = GenotypeFilterParam()
        ),
        "unused argument"
    )
    expect_error(GenotypeFilterParam(mafCutoff = -1), "mafCutoff.*is not >= 0")
    expect_error(
        GenotypeFilterParam(keepIndel = NA),
        "keepIndel.*May not be NA"
    )
    expect_error(
        GenotypeFilterParam(xvarCutoff = -1),
        "xvarCutoff.*is not >= 0"
    )
    # A bare list cannot be checked, so a run refuses it outright.
    expect_error(
        pecotmr:::.qtlResolveFilter(list(mafCutoff = 0.01)),
        "must be built with GenotypeFilterParam"
    )
    expect_error(
        loadQtlDatasetFromManifest(data.frame(), studyName = 1L),
        "studyName.*Must be of type 'string'"
    )
    expect_error(
        loadMultiStudyQtlDatasetFromManifest(
            data.frame(),
            minLdOverlapWarn = 2
        ),
        "minLdOverlapWarn.*is not <= 1"
    )
})

# ---------------------------------------------------------------------------
# AF vs MAF on the GWAS-VCF path
#
# `AF` is the spec's ALTERNATIVE-allele frequency -- directional, and what the
# posterior view publishes as `af`. It used to be a fallback for the MAF slot
# and was then folded with pmin(af, 1 - af), so an AF of 0.8 came back as a
# MAF of 0.2 with no `AF` column at all: pecotmr's own VCF output did not
# round-trip, and the direction the standardized-to-raw weight conversion
# needs was gone. The mirror of the released writer bug where `AF` held MAF.
# ---------------------------------------------------------------------------

# @noRd
.afvcf_sumstats <- function(af = c(0.8, 0.3)) {
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
        SE = c(0.1, 0.1),
        z_supplied = c(TRUE, TRUE)
    )
    if (!is.null(af)) {
        cols$AF <- af
    }
    GwasSumStats(
        studyName = "s",
        entry = list(S4Vectors::`mcols<-`(gr, value = cols)),
        genome = "hg38"
    )
}

test_that("a directional AF survives the GWAS-VCF round trip", {
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(.afvcf_sumstats(), out)
    back <- .readSumStatsVcf(out, NULL, NULL, NULL, "x")
    # Unfolded, and under its own name -- 0.8 is an ALT frequency above 0.5,
    # which is exactly the case folding destroys.
    expect_equal(as.numeric(back$AF), c(0.8, 0.3))
    # A directionless MAF is still derived for QC, alongside AF rather than
    # instead of it.
    expect_equal(as.numeric(back$MAF), c(0.2, 0.3))
})

test_that("a GWAS-VCF with no AF has neither AF nor a folded one", {
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(.afvcf_sumstats(af = NULL), out)
    back <- .readSumStatsVcf(out, NULL, NULL, NULL, "x")
    expect_false(is_in("AF", names(back)))
    expect_false(is_in("MAF", names(back)))
})

test_that("an explicit formatMapping still redirects the AF slot", {
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(.afvcf_sumstats(), out)
    # Reading the spec's AF field into the MAF slot on purpose: the caller
    # says this file's AF is not to be trusted as directional.
    back <- .readSumStatsVcf(
        out,
        NULL,
        NULL,
        list(AF = "NOPE", MAF = "AF"),
        "x"
    )
    expect_false(is_in("AF", names(back)))
    expect_equal(as.numeric(back$MAF), c(0.2, 0.3))
})

test_that("a manifest's BED-style columns give each row its window", {
    gwasTsv <- system.file(
        "extdata",
        "manifests",
        "protocol_example.twas.gwas_sumstats.chr22.tsv.gz",
        package = "pecotmr"
    )
    ldStem <- file.path(
        system.file("extdata", "ld_reference", "chr22", package = "pecotmr"),
        "protocol_example.LD.chr22"
    )
    mf <- data.frame(
        study = "g1",
        sumStatsPath = gwasTsv,
        chrom = "chr22",
        start = 9999999L,
        end = 19000000L
    )
    g <- suppressWarnings(loadGwasSumStatsFromManifest(
        manifest = mf,
        genome = "hg38",
        ldReference = ldStem
    ))
    # The start is 0-based, as in BED.
    expect_identical(as.character(g$region), "chr22:10000000-19000000")
    expect_true(g$regionGiven)
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest = mf,
            genome = "hg38",
            ldReference = ldStem,
            which = GenomicRanges::GRanges("chr22:1-2")
        ),
        "cannot give them too"
    )
    # A QTL row may also give its trait's position, on its `chrom`.
    qmf <- data.frame(
        study = "s1",
        context = "c",
        trait = "t",
        sumStatsPath = gwasTsv,
        chrom = "chr22",
        start = 9999999L,
        end = 19000000L,
        traitStart = 15527999L,
        traitEnd = 15528000L
    )
    q <- suppressWarnings(loadQtlSumStatsFromManifest(
        manifest = qmf,
        genome = "hg38",
        ldReference = ldStem
    ))
    expect_identical(as.character(q$region), "chr22:10000000-19000000")
    expect_true(q$regionGiven)
    expect_identical(as.character(q$traitPos), "chr22:15528000")
    # A QTL tuple has one window.
    expect_error(
        loadQtlSumStatsFromManifest(
            manifest = qmf[, 1:4],
            genome = "hg38",
            ldReference = ldStem,
            which = GenomicRanges::GRanges(c("chr22:1-10", "chr22:20-30"))
        ),
        "must be one window"
    )
    # A window needs all three columns, on every row.
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest = mf[, c("study", "sumStatsPath", "start", "end")],
            genome = "hg38",
            ldReference = ldStem
        ),
        "missing chrom"
    )
    holes <- rbind(
        mf,
        data.frame(
            study = "g2",
            sumStatsPath = gwasTsv,
            chrom = "chr22",
            start = NA,
            end = 19000000L
        )
    )
    expect_error(
        loadGwasSumStatsFromManifest(
            manifest = holes,
            genome = "hg38",
            ldReference = ldStem
        ),
        "row(s) 2 leave their window empty",
        fixed = TRUE
    )
})

test_that(".matchRegionSeqnames reads the windows a file holds, once", {
    gr <- GenomicRanges::GRanges(
        c("chr22", "chr22", "chr3"),
        IRanges::IRanges(c(1L, 50L, 1L), c(100L, 200L, 10L))
    )
    m <- pecotmr:::.matchRegionSeqnames(gr, c("22", "21"))
    expect_identical(as.character(m), "22:1-200")
    expect_null(pecotmr:::.matchRegionSeqnames(gr[3L], "22"))
})
