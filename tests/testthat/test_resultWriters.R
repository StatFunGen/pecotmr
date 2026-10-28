# Tests for the published result writers (resultWriters.R): writeGwasVcf, the
# GWAS-VCF-style writer, and writeBed, the published long-format tables.

# =============================================================================
# Test data helpers
# =============================================================================

makeTestGenotypeHandle <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

makeTestGwasSumStats <- function(n = 5) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = seq(100, by = 100, length.out = n), width = 1)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("T", n),
        Z = seq(1.5, by = -0.5, length.out = n),
        N = rep(1000L, n)
    )
    GwasSumStats(
        studyName = "test_trait",
        entry = list(gr),
        genome = "hg38",
        ldReference = makeTestGenotypeHandle()
    )
}

makeTestFineMappingResult <- function(n = 5) {
    beta <- seq(0.5, by = -0.1, length.out = n)
    se <- rep(0.1, n)
    zv <- seq(5.0, by = -1.0, length.out = n)
    tl <- data.frame(
        variant_id = paste0(
            "chr1:",
            seq(100, by = 100, length.out = n),
            ":T:A"
        ),
        chrom = rep("1", n),
        pos = as.integer(seq(100, by = 100, length.out = n)),
        A1 = rep("A", n),
        A2 = rep("T", n),
        N = rep(1000, n),
        MAF = rep(0.1, n),
        marginal_beta = beta,
        marginal_se = se,
        marginal_z = zv,
        marginal_p = 2 * pnorm(-abs(zv)),
        pip = seq(0.9, by = -0.1, length.out = n),
        posterior_mean = beta * 0.5,
        posterior_sd = se * 0.5,
        cs_95 = paste0("susie_", c(1L, 1L, 0L, 2L, 0L)[seq_len(n)]),
        method = rep("susie", n),
        stringsAsFactors = FALSE
    )
    entry <- FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    GwasFineMappingResult(
        studyName = "test_study",
        method = "susie",
        entry = list(entry)
    )
}

# =============================================================================
# Fine-mapping result to VCF
# =============================================================================

test_that("writeGwasVcf(FineMapping): emits PIP + CS FORMAT fields", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")

    fm <- makeTestFineMappingResult(5)
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    expect_equal(writeGwasVcf(fm, out), out)
    ln <- readLines(out)
    # The posterior fields the legacy create_vcf carried must be present. The CS
    # field is DYNAMIC (CS<coverage>) -- the fixture carries cs_95 -> CS95.
    expect_true(any(grepl("ID=PIP,", ln, fixed = TRUE)))
    expect_true(any(grepl("ID=CS95,", ln, fixed = TRUE)))
    data1 <- ln[!grepl("^#", ln)][[1L]]
    fmtKeys <- strsplit(strsplit(data1, "\t")[[1L]][[9L]], ":")[[1L]]
    expect_true("PIP" %in% fmtKeys)
    expect_true("CS95" %in% fmtKeys)
    # cs_95 "susie_<k>" -> integer CS index; fixture row 1 is in CS 1.
    vals <- strsplit(strsplit(data1, "\t")[[1L]][[10L]], ":")[[1L]]
    expect_equal(vals[[which(fmtKeys == "CS95")]], "1")
})

# =============================================================================
# Fine-mapping result to BCF
# =============================================================================

test_that("writeGwasVcf writes a fine-mapping result to BCF", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("Rsamtools")

    # asBcf may be temporarily disabled in some Rsamtools versions
    asbcf_works <- tryCatch(
        {
            tmp_stem <- tempfile(fileext = ".vcf")
            tmp_bgz <- paste0(tmp_stem, ".bgz")
            gr <- GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(start = 100, width = 1, names = "chr1:100:A:G")
            )
            hdr <- VariantAnnotation::VCFHeader(
                header = IRanges::DataFrameList(
                    fileformat = S4Vectors::DataFrame(
                        Value = "VCFv4.2",
                        row.names = "fileformat"
                    )
                ),
                sample = "probe"
            )
            cd <- S4Vectors::DataFrame(Samples = "probe", row.names = "probe")
            v <- VariantAnnotation::VCF(
                rowRanges = gr,
                colData = cd,
                exptData = list(header = hdr)
            )
            VariantAnnotation::ref(v) <- Biostrings::DNAStringSet("A")
            VariantAnnotation::alt(v) <- Biostrings::DNAStringSetList(list("T"))
            VariantAnnotation::fixed(v)$FILTER <- "PASS"
            VariantAnnotation::writeVcf(v, tmp_stem, index = TRUE)
            bcf_stem <- tempfile()
            Rsamtools::asBcf(
                tmp_bgz,
                dictionary = "chr1",
                destination = bcf_stem
            )
            TRUE
        },
        error = function(e) FALSE
    )
    skip_if(!asbcf_works, "Rsamtools::asBcf is not functional in this build")

    fm <- makeTestFineMappingResult(5)
    out <- tempfile(fileext = ".bcf")
    on.exit(unlink(out), add = TRUE)

    result <- writeGwasVcf(fm, out)
    expect_equal(result, out)
    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 0)
})

# =============================================================================
# Empty fine-mapping result errors
# =============================================================================

test_that("writeGwasVcf errors on an empty fine-mapping result", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")

    empty_tl <- data.frame(
        variant_id = character(0),
        chrom = character(0),
        pos = integer(0),
        A1 = character(0),
        A2 = character(0),
        N = numeric(0),
        MAF = numeric(0),
        marginal_beta = numeric(0),
        marginal_se = numeric(0),
        marginal_z = numeric(0),
        marginal_p = numeric(0),
        pip = numeric(0),
        posterior_mean = numeric(0),
        posterior_sd = numeric(0),
        stringsAsFactors = FALSE
    )
    entry <- FineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = empty_tl
    )
    fm_empty <- GwasFineMappingResult(
        studyName = "test_study",
        method = "susie",
        entry = list(entry)
    )

    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)

    expect_error(writeGwasVcf(fm_empty, out), "no variants to write")
})

# =============================================================================
# splitByContext / splitByTrait — one file per tuple
# =============================================================================

.makeMultiTupleQtlFmr <- function() {
    contexts <- c("brain", "blood")
    traits <- c("ENSG_A", "ENSG_B")
    entries <- map(seq_along(contexts), function(i) {
        ids <- paste0("chr1:", 100 * (1:3), ":T:A")
        tl <- data.frame(
            variant_id = ids,
            chrom = rep("1", 3),
            pos = c(100L, 200L, 300L),
            A1 = rep("A", 3),
            A2 = rep("T", 3),
            N = rep(1000, 3),
            MAF = rep(0.1, 3),
            marginal_beta = c(0.3, 0.1, -0.2) + i / 100,
            marginal_se = rep(0.05, 3),
            marginal_z = c(6.0, 2.0, -4.0),
            marginal_p = c(1e-9, 0.045, 6e-5),
            pip = c(0.9, 0.5, 0.7),
            posterior_mean = rep(0.05, 3),
            posterior_sd = rep(0.02, 3),
            cs_95 = paste0("susie_", c(1L, 1L, 0L)),
            stringsAsFactors = FALSE
        )
        FineMappingRow(variantIds = ids, susieFit = list(), topLoci = tl)
    })
    QtlFineMappingResult(
        studyName = rep("study1", 2),
        context = contexts,
        trait = traits,
        method = rep("susie", 2),
        entry = entries
    )
}

test_that("writeGwasVcf(FineMapping): each split flag emits one VCF per tag", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fmr <- .makeMultiTupleQtlFmr()
    # Each path is decorated with the tags its flags select; both flags
    # combine them, context first.
    cases <- list(
        context = list(
            flags = list(splitByContext = TRUE),
            tags = c("brain", "blood")
        ),
        trait = list(
            flags = list(splitByTrait = TRUE),
            tags = c("ENSG_A", "ENSG_B")
        ),
        both = list(
            flags = list(splitByContext = TRUE, splitByTrait = TRUE),
            tags = c("brain.ENSG_A", "blood.ENSG_B")
        )
    )
    d <- withr::local_tempdir()
    for (case in names(cases)) {
        paths <- exec(
            writeGwasVcf,
            fmr,
            file.path(d, str_c(case, ".vcf")),
            !!!cases[[case]]$flags
        )
        expect_equal(length(paths), 2L, info = case)
        expect_true(
            setequal(
                basename(paths),
                str_c(case, ".", cases[[case]]$tags, ".vcf")
            ),
            info = case
        )
        expect_true(all(file.exists(paths)), info = case)
    }
})

test_that("writeGwasVcf(FineMapping): several fits make one multi-sample VCF", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fmr <- .makeMultiTupleQtlFmr()
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    # No selector and no split flag: every matching fit is a SAMPLE of one
    # file, which is what GWAS-VCF's per-sample design means.
    expect_identical(writeGwasVcf(fmr, out), out)
    lines <- readLines(out)
    header <- str_split(lines[str_starts(lines, "#CHROM")][[1L]], "\t")[[1L]]
    samples <- header[-seq_len(9L)]
    expect_length(samples, 2L)
    # each sample is named for the fit it came from, and described in a
    # ##SAMPLE header line
    expect_true(all(str_detect(samples, "\\|")))
    sampleLines <- keep(lines, str_starts, "##SAMPLE=")
    expect_length(sampleLines, 2L)
    expect_true(all(str_detect(sampleLines, "TotalVariants=")))
    # one record per variant, with a cell per sample
    body <- str_split(lines[!str_starts(lines, "#")][[1L]], "\t")[[1L]]
    expect_length(body, 9L + 2L)
})

test_that("writeGwasVcf(FineMapping): a fit missing a variant leaves a gap", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    # Two fits over DISJOINT variants: the records are the union, and each
    # sample is empty where its own fit had nothing.
    mk <- function(positions, ctx) {
        ids <- sprintf("chr1:%d:T:A", positions)
        tl <- data.frame(
            variant_id = ids,
            chrom = rep("1", length(ids)),
            pos = as.integer(positions),
            A1 = "A",
            A2 = "T",
            pip = seq(0.9, by = -0.1, length.out = length(ids)),
            cs_95 = "susie_1",
            method = "susie",
            stringsAsFactors = FALSE
        )
        QtlFineMappingResult(
            studyName = "s1",
            context = ctx,
            trait = "g1",
            method = "susie",
            entry = list(FineMappingRow(
                variantIds = ids,
                susieFit = list(),
                topLoci = tl
            ))
        )
    }
    fmr <- combineFineMappingResults(mk(c(100, 200), "a"), mk(c(300), "b"))
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    writeGwasVcf(fmr, out)
    lines <- readLines(out)
    cells <- str_split(keep(lines, str_starts, "chr1"), "\t")
    expect_identical(map_chr(cells, 2L), c("100", "200", "300"))
    keys <- str_split(cells[[1L]][[9L]], ":")[[1L]]
    pipOf <- function(cell) str_split(cell, ":")[[1L]][[which(keys == "PIP")]]
    # sample 1 covers 100/200, sample 2 only 300
    expect_identical(map_chr(map(cells, 10L), pipOf), c("0.9", "0.8", "."))
    expect_identical(map_chr(map(cells, 11L), pipOf), c(".", ".", "0.9"))
})

# =============================================================================
# GwasSumStats: multi-study selection, MAF/AF field, .gz rename
# =============================================================================

# Two studies over the same variants, told apart by their sample size: N is
# 1000 in studyA and 2000 in studyB.
.makeTwoStudyGwasSumStats <- function(n = 3) {
    mkGr <- function(nSample) {
        gr <- GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(
                start = seq(100, by = 100, length.out = n),
                width = 1
            )
        )
        mc <- S4Vectors::DataFrame(
            SNP = paste0("rs", seq_len(n)),
            A1 = rep("A", n),
            A2 = rep("T", n),
            Z = seq(1.5, by = -0.5, length.out = n),
            N = rep(nSample, n)
        )
        S4Vectors::mcols(gr) <- mc
        gr
    }
    GwasSumStats(
        studyName = c("studyA", "studyB"),
        entry = list(mkGr(1000L), mkGr(2000L)),
        genome = "hg38",
        ldReference = makeTestGenotypeHandle()
    )
}

test_that("writeGwasVcf(SumStats): multi-study needs a `studyName`", {
    skip_if_not_installed("VariantAnnotation")
    ss2 <- .makeTwoStudyGwasSumStats()
    out <- tempfile(fileext = ".vcf")
    expect_error(writeGwasVcf(ss2, out), "2 studies")
})

test_that("writeGwasVcf(SumStats): `studyName` picks the chosen study", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    ss2 <- .makeTwoStudyGwasSumStats()
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    res <- writeGwasVcf(ss2, out, studyName = "studyB")
    expect_equal(res, out)
    # The file holds studyB alone: one sample, named for it, with its N. A
    # Number=A field reads back as a list matrix, one value per cell.
    vcf <- VariantAnnotation::readVcf(out)
    expect_identical(colnames(vcf), "studyB")
    expect_identical(
        list_c(VariantAnnotation::geno(vcf)$NS[, 1L]),
        rep(2000L, 3L)
    )
})

test_that("writeGwasVcf(SumStats): emits AF when MAF is present", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    ss <- GwasSumStats(
        studyName = "t",
        entry = list(local({
            gr <- GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(start = c(100, 200, 300, 400), width = 1)
            )
            S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
                SNP = paste0("rs", 1:4),
                A1 = rep("A", 4),
                A2 = rep("T", 4),
                Z = c(1.5, 1.0, 0.5, 0.0),
                N = rep(1000L, 4),
                MAF = rep(0.2, 4)
            )
            gr
        })),
        genome = "hg38",
        ldReference = makeTestGenotypeHandle()
    )
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    writeGwasVcf(ss, out)
    expect_true(file.exists(out))
    # The AF FORMAT field is declared in the header when MAF is present.
    # A source carrying only MAF gets the MAF key, not AF: the spec's AF is
    # the ALTERNATIVE allele frequency, which a minor-allele frequency cannot
    # be converted to -- it has lost which allele it refers to.
    expect_true(any(grepl("ID=MAF", readLines(out))))
    expect_false(any(grepl("ID=AF,", readLines(out))))
})

test_that("writeGwasVcf(SumStats): .vcf.gz is renamed from .bgz", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    ss <- makeTestGwasSumStats(5)
    out <- tempfile(fileext = ".vcf.gz")
    on.exit(unlink(c(out, paste0(out, ".tbi")), force = TRUE), add = TRUE)
    res <- writeGwasVcf(ss, out)
    expect_equal(res, out)
    # VariantAnnotation::writeVcf() always emits .bgz; the writer renames it
    # to the requested .gz.
    expect_true(file.exists(out))
    expect_false(file.exists(sub("\\.gz$", ".bgz", out)))
})

# =============================================================================
# Fine-mapping result: explicit selectors, no-match error, composite/empty
# tags
# =============================================================================

test_that("writeGwasVcf(FineMapping): selectors pick a single row", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fmr <- .makeMultiTupleQtlFmr()
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    res <- writeGwasVcf(
        fmr,
        out,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    expect_equal(res, out)
    # One sample, the brain / ENSG_A fit, carrying that row's marginal
    # effects; the blood row's are 0.01 larger.
    vcf <- VariantAnnotation::readVcf(out)
    expect_identical(colnames(vcf), "study1|brain|ENSG_A|susie")
    expect_equal(
        list_c(VariantAnnotation::geno(vcf)$ES[, 1L]),
        c(0.31, 0.11, -0.19)
    )
})

test_that("writeGwasVcf(fine-mapping result): no matching rows errors", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fmr <- .makeMultiTupleQtlFmr()
    out <- tempfile(fileext = ".vcf")
    expect_error(
        writeGwasVcf(fmr, out, studyName = "does_not_exist"),
        "no rows match"
    )
})

test_that("writeGwasVcf(FineMapping): split decorates a .vcf.bgz path", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fmr <- .makeMultiTupleQtlFmr()
    baseOut <- tempfile(fileext = ".vcf.bgz")
    stem <- tools::file_path_sans_ext(tools::file_path_sans_ext(baseOut))
    on.exit(
        unlink(
            list.files(
                dirname(baseOut),
                pattern = basename(stem),
                full.names = TRUE
            ),
            force = TRUE
        ),
        add = TRUE
    )
    paths <- writeGwasVcf(fmr, baseOut, splitByContext = TRUE)
    expect_length(paths, 2L)
    # The composite extension survives; the tag is injected before it.
    expect_true(all(grepl("\\.vcf\\.bgz$", paths)))
    expect_true(any(grepl("\\.brain\\.vcf\\.bgz$", paths)))
    expect_true(any(grepl("\\.blood\\.vcf\\.bgz$", paths)))
})

test_that("writeGwasVcf(FineMapping): a context-less split keeps the path", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    # GwasFineMappingResult has no context/trait axes, so the split tag set is
    # empty and the original output path is used unchanged.
    fm <- makeTestFineMappingResult(5)
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    paths <- writeGwasVcf(fm, out, splitByContext = TRUE)
    expect_length(paths, 1L)
    expect_equal(paths[[1L]], out)
    expect_true(file.exists(out))
})

test_that("writeGwasVcf(FineMapping): AF comes from topLoci `af`", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    # The marginal view projects allele frequency to a column named `af`, so a
    # topLoci carrying `af` must surface as the AF genotype field.
    n <- 3
    tl <- data.frame(
        variant_id = paste0("chr1:", c(100, 200, 300), ":T:A"),
        chrom = rep("1", n),
        pos = c(100L, 200L, 300L),
        A1 = rep("A", n),
        A2 = rep("T", n),
        N = rep(1000, n),
        af = c(0.1, 0.2, 0.3),
        marginal_beta = c(0.3, -0.2, 0.1),
        marginal_se = rep(0.05, n),
        marginal_z = c(6.0, -4.0, 2.0),
        marginal_p = c(1e-9, 6e-5, 0.045),
        pip = c(0.9, 0.7, 0.5),
        posterior_mean = rep(0.05, n),
        posterior_sd = rep(0.02, n),
        cs_95 = paste0("susie_", c(1L, 0L, 2L)),
        stringsAsFactors = FALSE
    )
    entry <- FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    fm <- GwasFineMappingResult(
        studyName = "s",
        method = "susie",
        entry = list(entry)
    )
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    writeGwasVcf(fm, out)
    expect_true(file.exists(out))
    expect_true(any(grepl("ID=AF", readLines(out))))
})

test_that("writeGwasVcf(FineMapping): LBF / LFSR / PUR / fullFit fields", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    # A rich topLoci exercises every posterior FORMAT field: per-variant
    # logBF (LBF), local false sign rate (LFSR), per-CS purity (PUR<cov>),
    # and the wide fullFit
    # columns (within_cs_pip / cs_logbf_<lab> / cs_effect_<lab>).
    n <- 3
    tl <- data.frame(
        variant_id = paste0("chr1:", c(100, 200, 300), ":T:A"),
        chrom = rep("1", n),
        pos = c(100L, 200L, 300L),
        A1 = rep("A", n),
        A2 = rep("T", n),
        N = rep(1000, n),
        marginal_beta = c(0.3, -0.2, 0.1),
        marginal_se = rep(0.05, n),
        marginal_z = c(6, -4, 2),
        marginal_p = c(1e-9, 6e-5, 0.045),
        pip = c(0.9, 0.7, 0.5),
        posterior_mean = rep(0.05, n),
        posterior_sd = rep(0.02, n),
        logBF = c(3.1, 1.2, 0.4),
        lfsr = c(0.01, 0.2, 0.4),
        cs_95 = paste0("susie_", c(1L, 0L, 2L)),
        cs_95_purity = c(0.95, 0, 0.88),
        within_cs_pip = c(0.6, 0, 0.5),
        cs_logbf_susie_1 = c(2.2, NA, NA),
        cs_effect_susie_1 = c(0.4, NA, NA),
        stringsAsFactors = FALSE
    )
    entry <- FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    fm <- GwasFineMappingResult(
        studyName = "s",
        method = "susie",
        entry = list(entry)
    )
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    writeGwasVcf(fm, out)
    ln <- readLines(out)
    for (id in c(
        "ID=LBF",
        "ID=LFSR",
        "ID=PUR95",
        "ID=WITHIN_CS_PIP",
        "ID=CS_LOGBF_SUSIE_1",
        "ID=CS_EFFECT_SUSIE_1"
    )) {
        expect_true(any(grepl(id, ln, fixed = TRUE)), info = id)
    }
})

test_that("writeGwasVcf(FineMapping): falls back to marginal sumstats", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Biostrings")
    fm <- makeTestFineMappingResult(5)
    # Simulate an entry whose posterior table is unavailable (empty topLoci) so
    # the marginal univariate sumstats alone drive the body (the else-if).
    testthat::local_mocked_bindings(
        topLoci = function(x, ...) GenomicRanges::GRanges(),
        .package = "pecotmr"
    )
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    writeGwasVcf(fm, out)
    expect_true(file.exists(out))
    # marginal beta -> ES
    expect_true(any(grepl("ID=ES", readLines(out), fixed = TRUE)))
})


# ===========================================================================
# Output formats
# ===========================================================================

test_that("plain and bgzipped VCF are written, the latter with an index", {
    skip_if_not_installed("VariantAnnotation")
    data(gwasSumStatsS4Example)
    d <- withr::local_tempdir()
    # Each call returns the path it wrote.
    expect_equal(
        writeGwasVcf(gwasSumStatsS4Example, filename = file.path(d, "a.vcf")),
        file.path(d, "a.vcf")
    )
    expect_true(file.exists(file.path(d, "a.vcf")))
    expect_equal(
        writeGwasVcf(
            gwasSumStatsS4Example,
            filename = file.path(d, "b.vcf.bgz")
        ),
        file.path(d, "b.vcf.bgz")
    )
    expect_true(file.exists(file.path(d, "b.vcf.bgz")))
    # index = TRUE, so the tabix index is written alongside.
    expect_true(file.exists(file.path(d, "b.vcf.bgz.tbi")))
})

test_that("BCF output explains itself when Rsamtools cannot convert", {
    # Rsamtools >= 2.26 raises "asBcf() is temporarily disabled". Left bare
    # that surfaces as an opaque failure from a documented output format, so
    # it is translated into the cause plus the workaround.
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("Rsamtools")
    data(gwasSumStatsS4Example)
    d <- withr::local_tempdir()
    canConvert <- !inherits(
        tryCatch(
            Rsamtools::asBcf("nonexistent", "x", "y"),
            error = function(e) e
        ),
        "error"
    )
    skip_if(canConvert, "Rsamtools::asBcf() works here; nothing to translate")
    expect_error(
        writeGwasVcf(
            gwasSumStatsS4Example,
            filename = file.path(d, "c.bcf")
        ),
        "BCF output needs a working Rsamtools::asBcf"
    )
})

test_that(".vcfResolveBody uses posterior rows when no marginals exist", {
    local_mocked_bindings(
        topLoci = function(entry, signalCutoff) {
            # topLoci() returns ranges carrying the per-variant columns in
            # mcols(); these ids are deliberately not coordinate-shaped, so
            # the ranges are built directly rather than parsed from them.
            gr <- GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(start = c(1L, 2L), width = 1L)
            )
            S4Vectors::`mcols<-`(
                gr,
                value = S4Vectors::DataFrame(
                    variant_id = c("v1", "v2"),
                    pip = c(0.9, 0.1)
                )
            )
        },
        marginalEffects = function(entry) stop("none"),
        .package = "pecotmr"
    )
    out <- pecotmr:::.vcfResolveBody("e", "S1")
    expect_true(out$hasPost)
    expect_equal(nrow(out$base), 2L)
    # No marginal table to align against, so there is nothing to merge in.
    expect_null(out$m)
})

test_that(".writeVcfImpl writes UCSC chromosome names", {
    skip_if_not_installed("VariantAnnotation")
    out <- tempfile(fileext = ".vcf")
    on.exit(unlink(out), add = TRUE)
    geno <- matrix(c(1.5, -0.5), ncol = 1L, dimnames = list(NULL, "S1"))
    header <- S4Vectors::DataFrame(
        Number = "1",
        Type = "Float",
        Description = "Z score",
        row.names = "Z"
    )
    pecotmr:::.writeVcfImpl(
        chrom = c("1", "MT"),
        pos = c(100L, 200L),
        ref = c("A", "C"),
        alt = c("G", "T"),
        snpIds = c("rs1", "rs2"),
        geno = list(Z = geno),
        genoHeader = header,
        sampleName = "S1",
        outputPath = out
    )
    body <- readLines(out)
    body <- body[!startsWith(body, "#")]
    expect_setequal(stringr::str_extract(body, "^[^\t]+"), c("chr1", "chrM"))
})

# =============================================================================
# writeBed: the fine-mapping schemas
# =============================================================================

test_that("writeBed writes the published top-loci schema", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    expect_identical(
        writeBed(qtlFineMappingExample, out, signalCutoff = 0),
        out
    )
    tbl <- utils::read.delim(out, check.names = FALSE)
    # The whole schema, in order -- `context` is gone, `event_ID` already
    # carried it, as in the other two schemas.
    expect_identical(
        names(tbl)[1:6],
        c("#chr", "start", "end", "ID", "REF", "ALT")
    )
    # One vocabulary with the VCF: the spec's keys for what it defines, ours
    # alongside, snake_case for the per-fit fields.
    expect_true(all(c("ES", "SE", "LP", "NS") %in% names(tbl)))
    # EZ only when the source SUPPLIED z -- the spec defines it as the z used
    # to derive ES/SE, and this fixture recorded beta/se instead.
    expect_false("EZ" %in% names(tbl))
    expect_true(all(c("WGT", "WGTSE", "PIP", "LBF") %in% names(tbl)))
    # coverage levels are discovered, not hardcoded
    expect_true(all(c("CS95", "PUR95", "CS70", "CS50") %in% names(tbl)))
    expect_true(all(c("study", "context", "trait", "method") %in% names(tbl)))
    expect_true(all(
        c("region_id", "region_start", "region_end") %in% names(tbl)
    ))
    # and the columns this input has nothing to say about are absent, not
    # padded with NAs
    expect_false(any(c("MAF", "weight_standardized") %in% names(tbl)))
    expect_gt(nrow(tbl), 0L)
    # sorted by (chr, pos)
    first <- tbl[[1L]][[1L]]
    expect_false(is.unsorted(tbl$start[tbl[[1L]] == first]))
})

test_that("writeBed bgzips and tabix-indexes a .gz fsusie table", {
    data(fsusieFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed.gz")
    writeBed(fsusieFineMappingExample, out, schema = "fsusie")
    expect_true(file.exists(out))
    expect_true(file.exists(str_c(out, ".tbi")))
    hdr <- readLines(gzfile(out), n = 1L)
    cols <- str_split(hdr, "\t")[[1L]]
    expect_identical(cols[[1L]], "#chr")
    # MAF is absent from this fixture, so the file carries the schema minus
    # the columns that had nothing to say.
    expect_identical(
        cols,
        setdiff(
            .fmbFsusieColumns(
                c("wavelet", "trait"),
                c("CS95", "PUR95")
            ),
            "AF"
        )
    )
    expect_true(all(
        c(
            "grid_positions",
            "grid_effects",
            "grid_band_halfwidth",
            "joint_traits",
            "region_start",
            "region_end"
        ) %in%
            cols
    ))
})

test_that("writeBed emits the full schema for an empty result", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    # A cutoff no variant clears -> 0 data rows, but the schema must not
    # collapse. read.table() cannot read a header-only file, so read the line.
    writeBed(qtlFineMappingExample, out, signalCutoff = 2)
    lines <- readLines(out)
    expect_length(lines, 1L)
    # A header-only file documents the whole vocabulary rather than
    # collapsing to the required core.
    hdr <- str_split(lines[[1L]], "\t")[[1L]]
    expect_identical(hdr[1:6], c("#chr", "start", "end", "ID", "REF", "ALT"))
    expect_true(all(c("ES", "WGT", "PIP", "region_id") %in% hdr))
})

test_that("writeBed combines several results into one table", {
    data(qtlFineMappingExample)
    data(mvsusieFineMappingExample)
    one <- withr::local_tempfile(fileext = ".bed")
    both <- withr::local_tempfile(fileext = ".bed")
    writeBed(qtlFineMappingExample, one, signalCutoff = 0)
    writeBed(
        list(qtlFineMappingExample, mvsusieFineMappingExample),
        both,
        signalCutoff = 0
    )
    n1 <- length(readLines(one))
    n2 <- length(readLines(both))
    # One header either way, and the combined table carries strictly more rows.
    expect_gt(n2, n1)
    tbl <- utils::read.delim(both, check.names = FALSE)
    expect_setequal(
        unique(tbl$context),
        unique(c(
            as.character(contexts(qtlFineMappingExample)),
            as.character(contexts(mvsusieFineMappingExample))
        ))
    )
})

test_that("writeBed rejects something that is not a result collection", {
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(writeBed(data.frame(a = 1), out), "FineMappingResult")
    expect_error(writeBed(list(), out), "at least one")
})

test_that("the fsusie schema is still complete when no fsusie rows exist", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    # No fsusie rows at all -> a header-only table, not an error.
    writeBed(qtlFineMappingExample, out, schema = "fsusie")
    lines <- readLines(out)
    expect_length(lines, 1L)
    # a header-only file documents the schema, with no coverage discovered
    expect_identical(
        str_split(lines[[1L]], "\t")[[1L]],
        .fmbFsusieColumns(c("wavelet", "trait"))
    )
})

test_that("a merged multi-region fSuSiE table keeps each region's own window", {
    data(fsusieFineMappingExample)
    a <- fsusieFineMappingExample
    b <- fsusieFineMappingExample
    # Two regions whose (context, cs) labels collide -- the case a
    # context+cs key would silently merge onto the first region's window.
    S4Vectors::mcols(b)$study <- "studyB"
    wantA <- GenomicRanges::GRanges("chr22:1-9000000")
    wantB <- GenomicRanges::GRanges("chr22:2-8000000")
    withRegion <- function(x, want) {
        S4Vectors::mcols(x)$region <- rep(want, nrow(x))
        x
    }
    both <- combineFineMappingResults(
        withRegion(a, wantA),
        withRegion(b, wantB)
    )
    out <- withr::local_tempfile(fileext = ".bed.gz")
    writeBed(both, out, schema = "fsusie")
    tbl <- utils::read.delim(gzfile(out), check.names = FALSE)
    expect_setequal(
        unique(tbl$region_id),
        c("chr22_1_9000000", "chr22_2_8000000")
    )
    expect_setequal(unique(tbl$region_start), c(0, 1))
    expect_setequal(unique(tbl$region_end), c(9000000, 8000000))
})

test_that("extraColumns appends trailing fields the object cannot know", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        qtlFineMappingExample,
        out,
        signalCutoff = 0,
        extraColumns = list(source_file = "chr22_1_50000000")
    )
    tbl <- utils::read.delim(out, check.names = FALSE)
    expect_identical(names(tbl)[[ncol(tbl)]], "source_file")
    expect_true(all(tbl$source_file == "chr22_1_50000000"))

    # one value per row is equally valid
    perRow <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        qtlFineMappingExample,
        perRow,
        signalCutoff = 0,
        extraColumns = list(tag = seq_len(nrow(tbl)))
    )
    expect_identical(
        utils::read.delim(perRow, check.names = FALSE)$tag,
        seq_len(nrow(tbl))
    )
})

test_that("extraColumns refuses a clash, a bad length, an unnamed list", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(
        writeBed(
            qtlFineMappingExample,
            out,
            signalCutoff = 0,
            extraColumns = list(PIP = 1)
        ),
        "would overwrite"
    )
    expect_error(
        writeBed(
            qtlFineMappingExample,
            out,
            signalCutoff = 0,
            extraColumns = list(tag = c(1, 2))
        ),
        "one value per row"
    )
    expect_error(
        writeBed(
            qtlFineMappingExample,
            out,
            signalCutoff = 0,
            extraColumns = list(1)
        ),
        "NAMED list"
    )
})

test_that("extraColumns survive an empty table", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        qtlFineMappingExample,
        out,
        signalCutoff = 2,
        extraColumns = list(source_file = "b1")
    )
    lines <- readLines(out)
    expect_length(lines, 1L)
    # The header-only schema, with the extra column after it.
    plain <- withr::local_tempfile(fileext = ".bed")
    writeBed(qtlFineMappingExample, plain, signalCutoff = 2)
    expect_identical(
        str_split(lines[[1L]], "\t")[[1L]],
        c(str_split(readLines(plain)[[1L]], "\t")[[1L]], "source_file")
    )
})

test_that("genomicOrder sorts chromosomes naturally, not as text", {
    # 1..22 numerically, then the non-numeric contigs after all of them --
    # `sort -k1,1V` order, where 2 precedes 10.
    expect_identical(
        genomicOrder(c("10", "2", "X", "1"), c(1, 1, 1, 1)),
        c(4L, 2L, 1L, 3L)
    )
    expect_identical(
        genomicOrder(c("chr10", "chrX", "chr2", "chrY", "chr22"), rep(1, 5)),
        c(3L, 1L, 5L, 2L, 4L)
    )
    # a leading "chr" is cosmetic
    expect_identical(
        genomicOrder(c("chr2", "10"), c(1, 1)),
        genomicOrder(c("2", "chr10"), c(1, 1))
    )
    # positions sort numerically within a chromosome, not as text
    expect_identical(
        genomicOrder(rep("chr1", 3L), c("100", "9", "1000")),
        c(2L, 1L, 3L)
    )
    # a single chromosome is the common case and must be a no-op
    expect_identical(genomicOrder(rep("chr22", 3L), c(1, 2, 3)), 1:3)
})

test_that("genomicOrder uses `end` only to break ties, and only if numeric", {
    chrom <- rep("chr1", 3L)
    start <- c(10, 10, 5)
    expect_identical(
        genomicOrder(chrom, start, end = c(90, 20, 1)),
        c(3L, 2L, 1L)
    )
    # A third column that is not a position (an allele) is ignored rather than
    # coerced to NA and allowed to scramble the order.
    expect_identical(
        genomicOrder(chrom, start, end = c("A", "T", "G")),
        genomicOrder(chrom, start)
    )
    expect_identical(
        genomicOrder(chrom, start, end = NULL),
        genomicOrder(chrom, start)
    )
})

# =============================================================================
# writeBed: the TWAS-weight schemas
# =============================================================================

makeTestTwasWeights <- function(weights, method = "lasso", nVariants = 3L) {
    ids <- sprintf("chr1:%d:A:G", 100L * seq_len(nVariants))
    TwasWeights(
        studyName = "s1",
        context = "brain",
        trait = "g1",
        method = method,
        entry = list(TwasWeightsRow(
            variantIds = ids,
            weights = weights,
            cvResult = list(metrics = c(rsq = 0.5, pval = 0.01)),
            weightStandardized = TRUE
        ))
    )
}

test_that("writeBed writes the published TWAS-weight schema", {
    data(twasWeightsExample)
    out <- withr::local_tempfile(fileext = ".bed")
    expect_identical(writeBed(twasWeightsExample, out), out)
    tbl <- utils::read.delim(out, check.names = FALSE)
    expect_identical(
        names(tbl)[1:6],
        c("#chr", "start", "end", "ID", "REF", "ALT")
    )
    expect_true(all(
        c("WGT", "weight_standardized", "data_type") %in% names(tbl)
    ))
    # the posterior columns a weights table has nothing to say about are gone
    expect_false(any(c("PIP", "LBF", "CS95") %in% names(tbl)))
    # one row per (model, variant): 168 variants each
    # 2 of the 11 rows are fine-mapping-derived and are skipped
    expect_identical(nrow(tbl), 9L * 168L)
    expect_setequal(
        unique(tbl$method),
        setdiff(
            as.character(S4Vectors::mcols(twasWeightsExample)$method),
            .TWB_FINE_MAPPING_METHODS
        )
    )
    expect_true(all(tbl$end == tbl$start + 1L))
    expect_false(is.unsorted(tbl$start))
    expect_true(all(tbl$gene_ID == "ENSG00000130538"))
    # the analysis window, not the per-variant position
    expect_identical(unique(tbl$region_end), 18850784L)
})

test_that("writeBed writes one TWAS-weight row per matrix column", {
    w <- matrix(c(0.1, 0.2, 0.3, -0.4, -0.5, -0.6), nrow = 3L)
    colnames(w) <- c("mark_a", "mark_b")
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(makeTestTwasWeights(w), out)
    tbl <- utils::read.delim(out, check.names = FALSE)
    expect_identical(nrow(tbl), 6L)
    expect_setequal(unique(tbl$trait), c("mark_a", "mark_b"))
    expect_identical(
        tbl$WGT[tbl$trait == "mark_b"],
        c(-0.4, -0.5, -0.6)
    )
    # a plain weight vector has no feature axis, so the column is empty
    plain <- withr::local_tempfile(fileext = ".bed")
    writeBed(makeTestTwasWeights(c(0.1, 0.2, 0.3)), plain)
    plainTbl <- utils::read.delim(plain, check.names = FALSE)
    expect_identical(nrow(plainTbl), 3L)
    # a plain weight vector has no feature axis, so `trait` stays the row's
    expect_identical(unique(plainTbl$trait), "g1")
})

test_that("writeBed bgzips and tabix-indexes a .gz TWAS-weight table", {
    data(twasWeightsExample)
    out <- withr::local_tempfile(fileext = ".bed.gz")
    writeBed(twasWeightsExample, out)
    expect_true(file.exists(out))
    expect_true(file.exists(str_c(out, ".tbi")))
    cols <- str_split(readLines(gzfile(out), n = 1L), "\t")[[1L]]
    expect_true(all(c("#chr", "start", "end", "ID", "WGT") %in% cols))
})

test_that("exporting refuses a TwasWeights with no full-data fit", {
    # twasWeightsPipeline(fitFullData = FALSE) leaves each row an empty
    # element with no weights and only a cvResult. Fold-level results cannot
    # be published, so that is an error rather than an empty output file.
    cvOnly <- TwasWeights(
        studyName = "s1",
        context = "brain",
        trait = "g1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = character(0),
            weights = NULL,
            cvResult = list(metrics = c(rsq = 0.5))
        ))
    )
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(writeBed(cvOnly, out), "no full-data fit")
    expect_error(writeBed(cvOnly, out), "fitFullData")
    expect_error(
        writeGwasVcf(cvOnly, withr::local_tempfile(fileext = ".vcf")),
        "no full-data fit"
    )
    # the offending row is named, so a genome-wide run says which one
    expect_error(writeBed(cvOnly, out), "s1\\|brain\\|g1\\|susie")
    # and a weightless-but-populated row is refused for the same reason
    expect_error(
        writeBed(makeTestTwasWeights(NULL), out),
        "no full-data fit"
    )
})

test_that("exporting refuses a fine-mapping row with no fitted variants", {
    # Fine-mapping has no CV-only mode -- CV is attached on top of the fit --
    # but the export checks rather than assumes it.
    fmr <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(FineMappingRow(
            variantIds = character(0),
            susieFit = list(),
            topLoci = data.frame()
        ))
    )
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(writeBed(fmr, out), "no full-data fit")
})

test_that("writeBed combines a list of TwasWeights into one table", {
    one <- withr::local_tempfile(fileext = ".bed")
    both <- withr::local_tempfile(fileext = ".bed")
    a <- makeTestTwasWeights(c(0.1, 0.2, 0.3), method = "lasso")
    b <- makeTestTwasWeights(c(0.4, 0.5, 0.6), method = "enet")
    writeBed(a, one)
    writeBed(list(a, b), both)
    expect_identical(length(readLines(both)), 2L * length(readLines(one)) - 1L)
    tbl <- utils::read.delim(both, check.names = FALSE)
    expect_setequal(unique(tbl$method), c("lasso", "enet"))
})

test_that("writeBed rejects a list that holds neither result type", {
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(writeBed(list(), out), "at least one")
    expect_error(
        writeBed(list(1, 2), out),
        "TwasWeights, QtlFineMappingResult or GwasFineMappingResult"
    )
})

# =============================================================================
# writeGwasVcf: TWAS weights
# =============================================================================

test_that("writeGwasVcf writes TWAS weights with one sample per model", {
    data(twasWeightsExample)
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(twasWeightsExample, out)
    lines <- readLines(out)
    header <- str_split(
        lines[str_starts(lines, "#CHROM")][[1L]],
        "\t"
    )[[1L]]
    samples <- header[-seq_len(9L)]
    # Nine of the eleven rows: the two fine-mapping-derived ones are skipped,
    # the same as on the BED path.
    expect_length(samples, 9L)
    expect_true(all(str_detect(samples, "^s|protocol_example\\|bulk_rnaseq")))
    expect_true(any(str_detect(lines, "ID=WGT")))
    expect_identical(sum(!str_starts(lines, "#")), 168L)
    # the per-variant FORMAT column is the weight, one value per sample
    body <- str_split(lines[!str_starts(lines, "#")][[1L]], "\t")[[1L]]
    expect_identical(body[[9L]], "WGT")
    expect_length(body, 9L + 9L)
})

test_that("the TWAS-weight VCF's ##SAMPLE lines describe each model, not CV", {
    data(twasWeightsExample)
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(twasWeightsExample, out)
    sampleLines <- keep(readLines(out), str_starts, "##SAMPLE=")
    expect_length(sampleLines, 9L)
    expect_true(all(str_detect(sampleLines, "TotalVariants=")))
    susie <- keep(sampleLines, str_detect, fixed("Method=lasso,"))[[1L]]
    expect_true(str_detect(susie, "Trait=ENSG00000130538"))
    expect_true(str_detect(susie, "TotalVariants=168"))
    expect_true(str_detect(susie, "WeightStandardized=FALSE"))
    expect_true(str_detect(susie, "DataType=bulk_rnaseq"))
    # CV describes held-out folds; a published file has no data-secure way to
    # say anything fold-level, so the header stops at identity + counts.
    expect_false(any(str_detect(sampleLines, "rsq=|corr=|pval=|CvRegions=")))
})

test_that("writeGwasVcf narrows TWAS models with the identity selectors", {
    data(twasWeightsExample)
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(twasWeightsExample, out, method = "lasso")
    lines <- readLines(out)
    header <- str_split(
        lines[str_starts(lines, "#CHROM")][[1L]],
        "\t"
    )[[1L]]
    expect_length(header[-seq_len(9L)], 1L)
    expect_error(
        writeGwasVcf(twasWeightsExample, out, method = "nope"),
        "no rows match"
    )
    # Selecting a fine-mapping method leaves nothing exportable, which is a
    # different refusal from "no rows match".
    expect_error(
        suppressMessages(writeGwasVcf(
            twasWeightsExample,
            out,
            method = "susie"
        )),
        "every selected row is fine-mapping-derived"
    )
})

test_that("writeGwasVcf names TWAS samples, or refuses a mismatched count", {
    out <- withr::local_tempfile(fileext = ".vcf")
    tw <- makeTestTwasWeights(c(0.1, 0.2, 0.3))
    writeGwasVcf(tw, out, sampleName = "my_model")
    lines <- readLines(out)
    header <- str_split(
        lines[str_starts(lines, "#CHROM")][[1L]],
        "\t"
    )[[1L]]
    expect_identical(header[[10L]], "my_model")
    expect_error(
        writeGwasVcf(tw, out, sampleName = c("a", "b")),
        "one name per model"
    )
})

test_that("a TWAS-weight VCF unions the variants its models disagree on", {
    a <- makeTestTwasWeights(c(0.1, 0.2, 0.3), method = "lasso")
    b <- TwasWeights(
        studyName = "s1",
        context = "brain",
        trait = "g1",
        method = "enet",
        entry = list(TwasWeightsRow(
            variantIds = c("chr1:300:A:G", "chr1:400:A:G"),
            weights = c(0.7, 0.8)
        ))
    )
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(list(a, b), out)
    body <- keep(readLines(out), str_starts, "chr1")
    # union of {100,200,300} and {300,400}
    expect_length(body, 4L)
    cells <- str_split(body, "\t")
    pos <- map_chr(cells, 2L)
    expect_identical(pos, c("100", "200", "300", "400"))
    # the model that skipped a variant leaves its cell empty
    expect_identical(map_chr(cells, 10L), c("0.1", "0.2", "0.3", "."))
    expect_identical(map_chr(cells, 11L), c(".", ".", "0.7", "0.8"))
})

# =============================================================================
# Shared layer
# =============================================================================

test_that(".vcfAltList keeps a non-DNA ALT allele out of a DNAStringSet", {
    # VCF's "*" (a variant spanned by an upstream deletion) is a legal ALT but
    # not a DNA letter, and DNAStringSet rejects it outright.
    expect_s4_class(.vcfAltList(c("A", "GT")), "DNAStringSetList")
    expect_s4_class(.vcfAltList(c("A", "*")), "CharacterList")
    expect_s4_class(.vcfAltList(c("A", NA)), "CharacterList")
})

test_that(".pubIdentityKey keeps an absent axis from erasing the whole key", {
    # str_c() propagates NA, so a GWAS row's missing trait would make the
    # entire key NA -- and match() pairs NA with NA, joining unrelated rows.
    df <- data.frame(
        study = c("a", "b"),
        trait = c(NA_character_, NA_character_),
        variant_id = c("v1", "v1")
    )
    keys <- .pubIdentityKey(df, .PUB_VARIANT_COLS)
    expect_false(anyNA(keys))
    expect_false(keys[[1L]] == keys[[2L]])
})

test_that(".pubJoinMarginals gives an NA row where no marginal exists", {
    base <- data.frame(variant_id = c("v1", "v2"))
    marg <- data.frame(variant_id = "v1", beta = 0.5)
    joined <- .pubJoinMarginals(base, marg)
    expect_identical(nrow(joined), 2L)
    expect_identical(joined$beta, c(0.5, NA_real_))
})

# =============================================================================
# writeBed: variant-level column completeness
# =============================================================================

test_that("the top-loci schema identifies the fit each row came from", {
    data(qtlFineMappingExample)
    data(mvsusieFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        list(qtlFineMappingExample, mvsusieFineMappingExample),
        out,
        signalCutoff = 0
    )
    tbl <- utils::read.delim(out, check.names = FALSE)
    # Without study_ID / method a folded table cannot say which fit a row
    # came from, which is the whole point of folding several in.
    expect_true(all(c("study", "method") %in% names(tbl)))
    expect_setequal(
        unique(tbl$method),
        c(
            as.character(qtlFineMappingExample$method),
            as.character(mvsusieFineMappingExample$method)
        )
    )
    # the marginal statistics are complete: N and p, not only z, and they
    # are populated rather than present-but-empty. Only for the univariate
    # rows -- mvSuSiE has no marginal sumstats at all, so its rows carry no
    # betahat / sebetahat / z / N / p either, which is why these columns are
    # NA-able rather than required.
    expect_true(all(c("NS", "LP") %in% names(tbl)))
    univariate <- tbl[tbl$method == "susie", , drop = FALSE]
    expect_gt(nrow(univariate), 0L)
    expect_false(anyNA(univariate$NS))
    # LP is -log10 p of the marginal test, so it inverts to a probability
    expect_true(all(10^(-univariate$LP) > 0 & 10^(-univariate$LP) <= 1))
    expect_equal(
        10^(-univariate$LP),
        2 * pnorm(-abs(univariate$ES / univariate$SE)),
        tolerance = 1e-6
    )
    expect_true(all(is.na(tbl$LP[tbl$method == "mvsusie"])))
    # purity at every coverage the schema reports a credible set for
    expect_true(all(c("PUR95", "PUR70", "PUR50") %in% names(tbl)))
    # a true half-open BED like the other two schemas, so one reader rule
    # covers everything
    expect_identical(names(tbl)[1:3], c("#chr", "start", "end"))
    expect_true(all(tbl$end == tbl$start + 1L))
    # the analysis window travels with the row, so a read-back result does not
    # silently fall back to the span of the surviving variants
    expect_false(anyNA(tbl$region_ID))
    expect_true(all(tbl$region_start <= tbl$end))
})

test_that("columns = 'all' appends the rest of the source view", {
    data(qtlFineMappingExample)
    base <- withr::local_tempfile(fileext = ".bed")
    full <- withr::local_tempfile(fileext = ".bed")
    writeBed(qtlFineMappingExample, base, signalCutoff = 0)
    writeBed(qtlFineMappingExample, full, signalCutoff = 0, columns = "all")
    baseCols <- names(utils::read.delim(base, check.names = FALSE))
    fullCols <- names(utils::read.delim(full, check.names = FALSE))
    # additive: the published columns keep their place, the rest follow
    expect_identical(fullCols[seq_along(baseCols)], baseCols)
    # "all" is a request for whatever is there, so it stays subject to the
    # drop-empty rule: of the remaining view columns only `theta` holds values
    # in this fixture, so it alone is appended, while the empty `gene`,
    # `event` and `grange_*` stay out. An explicitly NAMED column is kept even
    # when empty -- that is what the next test checks.
    expect_identical(setdiff(fullCols, baseCols), "theta")
    # the range columns an as.data.frame() adds duplicate chr/pos, so they
    # must not come back under a second name
    expect_false(any(c("seqnames", "width", "strand") %in% fullCols))
})

test_that("columns takes an explicit set, and refuses a bad one", {
    data(qtlFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        qtlFineMappingExample,
        out,
        signalCutoff = 0,
        columns = "grange_start"
    )
    cols <- names(utils::read.delim(out, check.names = FALSE))
    expect_identical(cols[[length(cols)]], "grange_start")
    expect_error(
        writeBed(
            qtlFineMappingExample,
            out,
            signalCutoff = 0,
            columns = "not_a_column"
        ),
        "does not carry"
    )
    expect_error(
        writeBed(
            qtlFineMappingExample,
            out,
            signalCutoff = 0,
            columns = "pip"
        ),
        "already writes"
    )
})

test_that("the fsusie schema: both resolutions, its study, columns = 'all'", {
    data(fsusieFineMappingExample)
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(fsusieFineMappingExample, out, schema = "fsusie")
    tbl <- utils::read.delim(out, check.names = FALSE)
    # MAF is absent from this fixture, so the file carries the schema minus
    # the columns that had nothing to say.
    expect_identical(
        names(tbl),
        setdiff(
            .fmbFsusieColumns(
                c("wavelet", "trait"),
                c("CS95", "PUR95")
            ),
            "AF"
        )
    )
    # Both resolutions by default.
    expect_true(all(.FMB_FSUSIE_WAVELET_COLUMNS %in% names(tbl)))
    expect_true(all(.FMB_FSUSIE_TRAIT_COLUMNS %in% names(tbl)))
    expect_false(anyNA(tbl$study))

    full <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        fsusieFineMappingExample,
        full,
        schema = "fsusie",
        columns = "all"
    )
    fullCols <- names(utils::read.delim(full, check.names = FALSE))
    # "all" is additive: the published columns keep their place and the
    # non-empty remainder of the credible-set view follows.
    expect_identical(fullCols[seq_along(names(tbl))], names(tbl))
    expect_true(all(names(tbl) %in% fullCols))
    # an explicitly named column is written even though it is empty here
    named <- withr::local_tempfile(fileext = ".bed")
    writeBed(
        fsusieFineMappingExample,
        named,
        schema = "fsusie",
        columns = "cs_95_purity"
    )
    namedCols <- names(utils::read.delim(named, check.names = FALSE))
    expect_identical(namedCols[[length(namedCols)]], "cs_95_purity")
})

test_that("a TwasWeights element has no variant columns left to add", {
    data(twasWeightsExample)
    base <- withr::local_tempfile(fileext = ".bed")
    full <- withr::local_tempfile(fileext = ".bed")
    writeBed(twasWeightsExample, base)
    writeBed(twasWeightsExample, full, columns = "all")
    expect_identical(
        names(utils::read.delim(full, check.names = FALSE)),
        names(utils::read.delim(base, check.names = FALSE))
    )
})

test_that("fsusie `levels` drops one resolution, leaving a strict subset", {
    data(fsusieFineMappingExample)
    read <- function(lv) {
        out <- withr::local_tempfile(
            fileext = ".bed",
            .local_envir = parent.frame()
        )
        writeBed(
            fsusieFineMappingExample,
            out,
            schema = "fsusie",
            levels = lv
        )
        utils::read.delim(out, check.names = FALSE)
    }
    wav <- read("wavelet")
    tr <- read("trait")
    both <- read(c("wavelet", "trait"))

    # each level keeps its own group and drops the other's
    expect_true(all(.FMB_FSUSIE_WAVELET_COLUMNS %in% names(wav)))
    expect_false(any(.FMB_FSUSIE_TRAIT_COLUMNS %in% names(wav)))
    expect_true(all(.FMB_FSUSIE_TRAIT_COLUMNS %in% names(tr)))
    expect_false(any(
        .FMB_FSUSIE_WAVELET_COLUMNS %in%
            names(wav) &
            .FMB_FSUSIE_WAVELET_COLUMNS %in% names(tr)
    ))

    # a strict subset: same rows, same order, columns removed not reordered
    expect_identical(nrow(wav), nrow(both))
    expect_identical(
        names(wav),
        setdiff(names(both), .FMB_FSUSIE_TRAIT_COLUMNS)
    )
    expect_identical(
        names(tr),
        setdiff(names(both), .FMB_FSUSIE_WAVELET_COLUMNS)
    )
    expect_identical(wav$cs_id, both$cs_id)
    expect_identical(wav$grid_effects, both$grid_effects)

    expect_error(
        writeBed(
            fsusieFineMappingExample,
            withr::local_tempfile(),
            schema = "fsusie",
            levels = "nope"
        ),
        "wavelet"
    )
})

# ---------------------------------------------------------------------------
# EZ on the GwasSumStats path
#
# The spec defines EZ as the z "provided if it was used to derive the ES and
# SE fields". The fine-mapping writer already obeyed that through
# `.pubEz()`; this writer published any `Z` column it found, so a Wald z that
# pecotmr itself computed from a recorded beta/se went out as though the
# source had provided it -- which tells a reader the opposite of the truth
# about where ES and SE came from.
# ---------------------------------------------------------------------------

# @noRd
.ezFields <- function(x, name) {
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(x, out)
    lines <- readLines(out)
    record <- keep(lines, str_starts, "chr1")[[1L]]
    str_split(str_split(record, "\t")[[1L]][[9L]], ":")[[1L]]
}

test_that("writeGwasVcf(GwasSumStats) publishes EZ only for a supplied z", {
    expect_true(is_in("EZ", .ezFields(.ezss(TRUE))))
    expect_false(is_in("EZ", .ezFields(.ezss(FALSE))))
    # A collection built by hand that never recorded the provenance is
    # treated as not having supplied one: that is the safe direction, since
    # claiming ES/SE were derived from a z is the assertion that misleads.
    expect_false(is_in("EZ", .ezFields(.ezss(NULL))))
})

test_that("a supplied z survives the GwasSumStats VCF round trip", {
    skip_if_not_installed("VariantAnnotation")
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(.ezss(TRUE), out)
    back <- .readSumStatsVcf(out, NULL, NULL, NULL, "x")
    expect_equal(back$Z, c(2.5, -1.2))
    expect_true(all(back$z_supplied))
    # The derived case comes back derived, so a re-export does not promote it.
    writeGwasVcf(.ezss(FALSE), out)
    derived <- .readSumStatsVcf(out, NULL, NULL, NULL, "x")
    expect_false(any(derived$z_supplied))
    expect_equal(derived$Z, c(2.5, -1.2))
})

# ---------------------------------------------------------------------------
# `coverage =`: the dynamic credible-set levels
#
# Which coverages exist is a FIT-time decision (CredibleSetParam), and the
# writers READ the stored cs_<n> columns rather than recomputing, so a writer
# can restrict what it publishes but can never invent a level. The error says
# what IS available instead of emitting a silently empty column.
# ---------------------------------------------------------------------------

test_that("coverage = selects a subset of the levels the fit computed", {
    data(qtlFineMappingExample, envir = environment())
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(qtlFineMappingExample, out, coverage = 0.95)
    nms <- names(read.delim(out, check.names = FALSE))
    expect_true(is_in("CS95", nms))
    # The fit computed 70 and 50 as well; asking for 95 drops them.
    expect_false(any(is_in(c("CS70", "CS50"), nms)))
})

test_that("coverage = names what the fit did compute, when asked for more", {
    data(qtlFineMappingExample, envir = environment())
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(
        writeBed(qtlFineMappingExample, out, coverage = c(0.95, 0.8)),
        "asks for 80, which this fit did not compute.*Available: 95, 70, 50"
    )
})

test_that(".pubCsLevels reports the available set even when it is empty", {
    # A view with no cs_<n> column at all: the message has to say "<none>"
    # rather than interpolate an empty vector away.
    expect_error(
        .pubCsLevels(tibble(pip = 0.5), coverage = 0.95),
        "Available: <none>"
    )
    expect_equal(.pubCsLevels(tibble(cs_95 = "susie_1")), "95")
})

test_that("the fsusie schema refuses a result with no credible-set column", {
    expect_error(
        .fmbPrimaryCsLevel(tibble(pip = 0.5), NULL),
        "carries no credible-set column, so there are no rows to key on"
    )
})

# ---------------------------------------------------------------------------
# EZ on the fine-mapping path, and AF from the LD panel
# ---------------------------------------------------------------------------

test_that("a fine-mapping table publishes EZ when the source supplied z", {
    # `.pubEz()` is the fine-mapping counterpart of the GwasSumStats rule: the
    # z is published only where `marginal_z_supplied` says the source gave it,
    # never where pecotmr derived it from beta/se.
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
        chrom = "chr1",
        pos = c(100L, 200L),
        A1 = c("G", "T"),
        A2 = c("A", "C"),
        pip = c(0.9, 0.1),
        posterior_mean = c(0.2, -0.1),
        posterior_sd = c(0.05, 0.05),
        marginal_beta = c(0.3, -0.2),
        marginal_se = c(0.1, 0.1),
        marginal_z = c(5.5, -2.5),
        marginal_z_supplied = c(TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    row <- FineMappingRow(tl$variant_id, susieFit = list(), topLoci = tl)
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie",
        entry = list(row)
    )
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(res, out, signalCutoff = 0)
    tbl <- read.delim(out, check.names = FALSE)
    expect_equal(tbl$EZ, c(5.5, -2.5))
    # Flip the flag and the column goes away entirely rather than publishing a
    # z pecotmr computed itself.
    tl$marginal_z_supplied <- c(FALSE, FALSE)
    row2 <- FineMappingRow(tl$variant_id, susieFit = list(), topLoci = tl)
    res2 <- QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie",
        entry = list(row2)
    )
    writeBed(res2, out, signalCutoff = 0)
    expect_false(is_in("EZ", names(read.delim(out, check.names = FALSE))))
})

test_that("a weights export with no LD panel publishes no AF", {
    # A TwasWeights element records only alleles and a weight, so the only
    # frequency available is the panel's -- and AF is what a
    # standardized-to-raw conversion would need, so its absence is worth
    # seeing rather than faking.
    data(twasWeightsExample, envir = environment())
    expect_identical(ldReference(twasWeightsExample), list())
    out <- withr::local_tempfile(fileext = ".bed")
    suppressMessages(writeBed(twasWeightsExample, out))
    # An all-NA column is dropped rather than written as a column of ".",
    # so the absence of a frequency is visible in the file's shape.
    expect_false(is_in("AF", names(read.delim(out, check.names = FALSE))))
})

test_that("every weight row being fine-mapping-derived is an error", {
    # Those weights are the fine-mapping posterior re-expressed, so there is
    # nothing left to write once they are skipped.
    data(twasWeightsExample, envir = environment())
    fm <- twasWeightsExample[
        is_in(
            as.character(twasWeightsExample$method),
            .TWB_FINE_MAPPING_METHODS
        )
    ]
    expect_gt(nrow(fm), 0L)
    out <- withr::local_tempfile(fileext = ".bed")
    expect_error(
        suppressMessages(writeBed(fm, out)),
        "every row is fine-mapping-derived"
    )
})

test_that("the output path is left alone when no split flag is set", {
    expect_equal(
        .decorateOutputPath(
            "out.vcf.bgz",
            list(context = "brain"),
            FALSE,
            FALSE
        ),
        "out.vcf.bgz"
    )
})

test_that("an identity-less view keys on nothing; an empty range has no id", {
    # A bare single-entry view carries none of the identity axes, so the join
    # key degenerates to one empty string per row -- which is the right join
    # for one entry.
    expect_equal(
        .pubIdentityKey(tibble(pip = c(0.1, 0.2)), .PUB_ROW_COLS),
        c("", "")
    )
    expect_equal(.pubRegionIds(GenomicRanges::GRanges()), character(0))
    expect_equal(.pubRegionIds(NULL), character(0))
})

test_that("both formats skip fine-mapping-derived weight rows alike", {
    # The BED path skipped them and the VCF path published them, so the two
    # disagreed about what a weights file contains. The rule is the same in
    # both: those weights are the fine-mapping posterior re-expressed, and the
    # fine-mapping result says strictly more.
    skip_if_not_installed("VariantAnnotation")
    data(twasWeightsExample, envir = environment())
    derived <- is_in(
        as.character(twasWeightsExample$method),
        .TWB_FINE_MAPPING_METHODS
    )
    expect_gt(sum(derived), 0L)
    out <- withr::local_tempfile(fileext = ".vcf")
    expect_message(
        writeGwasVcf(twasWeightsExample, out),
        "writeGwasVcf: skipping \\d+ fine-mapping-derived"
    )
    header <- keep(readLines(out), str_starts, "#CHROM")[[1L]]
    samples <- str_split(header, "\t")[[1L]][-seq_len(9L)]
    expect_length(samples, sum(!derived))
    # None of the samples names a fine-mapping method.
    expect_false(any(str_detect(
        samples,
        str_c("\\|", .TWB_FINE_MAPPING_METHODS, "$", collapse = "|")
    )))
})

test_that("a weights export reads AF off the collection's LD panel", {
    # A TwasWeights element records only alleles and a weight, so the panel is
    # the only place a frequency can come from. The bundled sketch panel
    # records its variants' frequencies (from its .afreq), and they are what
    # is published; a panel recording none leaves the column absent rather
    # than faked.
    data(ctwasWeightsExample, twasWeightsExample, envir = environment())
    panel <- ldReference(ctwasWeightsExample)[[1L]]
    expect_false(is.null(panel))
    lasso <- twasWeightsExample[
        as.character(twasWeightsExample$method) == "lasso"
    ]
    expect_equal(nrow(lasso), 1L)
    exported <- function(ref) {
        withPanel <- TwasWeights(
            studyName = as.character(lasso$study),
            context = as.character(lasso$context),
            trait = as.character(lasso$trait),
            method = as.character(lasso$method),
            entry = list(.twrRowParts(lasso, 1L)),
            ldReference = ref
        )
        out <- withr::local_tempfile(fileext = ".bed")
        writeBed(withPanel, out)
        read.delim(out, check.names = FALSE)
    }
    tbl <- exported(panel)
    rr <- SummarizedExperiment::rowRanges(panel)
    expect_equal(tbl$AF, rr$AF[match(tbl$ID, .grVariantIds(rr))])
    expect_false(anyNA(tbl$AF))
    keep <- c("SNP", "A1", "A2")
    noAf <- SummarizedExperiment::`rowData<-`(
        panel,
        value = SummarizedExperiment::rowData(panel)[, keep]
    )
    expect_false(is_in("AF", names(exported(noAf))))
})

test_that("a sumstats collection with no Z publishes no EZ", {
    skip_if_not_installed("VariantAnnotation")
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100L, 200L), width = 1L)
    )
    gr <- S4Vectors::`mcols<-`(
        gr,
        value = S4Vectors::DataFrame(
            SNP = c("chr1:100:A:G", "chr1:200:C:T"),
            A1 = c("G", "T"),
            A2 = c("A", "C"),
            N = c(1000L, 1000L),
            BETA = c(0.25, -0.12),
            SE = c(0.1, 0.1)
        )
    )
    ss <- GwasSumStats(
        studyName = "s",
        entry = list(gr),
        genome = "hg38"
    )
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(ss, out)
    expect_false(any(str_detect(readLines(out), "ID=EZ,")))
})

# ---------------------------------------------------------------------------
# A joint (matrix-weight) model, and AF from the LD panel
# ---------------------------------------------------------------------------

# A multi-trait weight matrix on a panel that records allele frequency: the
# two shapes the single-trait fixtures never exercise. `mr_mash` because the
# matrix-weight fixture is a susie fit, which both exporters refuse as
# fine-mapping-derived -- a joint TWAS fit is the case that reaches the
# matrix paths at all.
# @noRd
.pbw_jointWeights <- function(nVariant = 6L) {
    data(ctwasWeightsExample, envir = environment())
    panel <- ldReference(ctwasWeightsExample)[[1L]]
    vids <- .grVariantIds(SummarizedExperiment::rowRanges(panel))
    SummarizedExperiment::rowRanges(panel)$AF <- seq(
        0.05,
        0.45,
        length.out = length(vids)
    )
    w <- matrix(
        c(0.1, -0.2, 0.3, 0, 0.05, -0.4, 0.2, 0.1, -0.1, 0.3, 0, 0.15),
        nrow = 6L,
        dimnames = list(NULL, c("tA", "tB"))
    )
    TwasWeights(
        studyName = "S1",
        context = "ctx",
        trait = "tA",
        method = "mr_mash",
        entry = list(TwasWeightsRow(
            variantIds = vids[seq_len(nVariant)],
            weights = w[seq_len(nVariant), , drop = FALSE]
        )),
        ldReference = panel
    )
}

test_that("AF is published from the LD panel when it records one", {
    tw <- .pbw_jointWeights()
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(tw, out)
    bed <- read.delim(out, check.names = FALSE)
    panel <- ldReference(tw)[[1L]]
    rr <- SummarizedExperiment::rowRanges(panel)
    expected <- mcols(rr)$AF[match(bed$ID, .grVariantIds(rr))]
    expect_equal(bed$AF, expected)
})

test_that("a joint fit's VCF samples carry the weight-column label", {
    skip_if_not_installed("VariantAnnotation")
    tw <- .pbw_jointWeights()
    out <- withr::local_tempfile(fileext = ".vcf")
    writeGwasVcf(tw, out)
    hdr <- readLines(out)
    chrom <- hdr[str_starts(hdr, "#CHROM")]
    # The label is appended to the row's identity, so two weight columns of
    # one row stay two distinguishable samples.
    expect_match(chrom, "S1\\|ctx\\|tA\\|mr_mash\\|tA")
    expect_match(chrom, "S1\\|ctx\\|tA\\|mr_mash\\|tB")
})

test_that("a recorded data_source wins over the derived one", {
    data(qtlFineMappingExample, envir = environment())
    x <- qtlFineMappingExample
    # These rows record where they came from. The collection carries no LD
    # panel, so the value the writer would DERIVE is "individual" -- and a
    # published file must report what the rows recorded instead of
    # re-deriving it from the object it happens to be held in.
    recorded <- QtlFineMappingResult(
        studyName = as.character(x$study),
        context = as.character(x$context),
        trait = as.character(x$trait),
        method = as.character(x$method),
        entry = map(seq_len(nrow(x)), .fmrRowParts, x = x),
        dataSource = rep("sumstat", nrow(x))
    )
    expect_equal(.pubDataSource(recorded), "individual")
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(recorded, out)
    bed <- read.delim(out, check.names = FALSE)
    expect_equal(unique(bed$data_source), "sumstat")
})

test_that("columns= appends a hand-built row's own variant columns", {
    # `TwasWeightsRow()` writes exactly A1 / A2 / weight, which is the whole
    # published source set -- so for a constructor-built collection there is
    # nothing left for `columns` to add, as ?writeBed says. A collection
    # assembled by hand is the case the argument exists for: its variants
    # carry fields pecotmr knows nothing about, and they travel into the file
    # AFTER the published schema rather than displacing any of it.
    vids <- sprintf("chr1:%d:A:G", 100L * (1:4))
    gr <- .variantIdsToGRanges(vids, "variantIds")
    mcols(gr)$weight <- c(0.1, -0.2, 0.3, 0)
    mcols(gr)$n_fold <- rep(5L, 4L)
    mcols(gr)$qc_note <- c("ok", "ok", "flipped", "ok")
    row <- methods::new(
        "TwasWeightsRow",
        variants = gr,
        weights = c(0.1, -0.2, 0.3, 0),
        weightStandardized = FALSE
    )
    tw <- TwasWeights(
        studyName = "S",
        context = "c",
        trait = "t",
        method = "lasso",
        entry = list(row)
    )
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(tw, out, columns = "all")
    bed <- read.delim(out, check.names = FALSE)
    expect_equal(bed$qc_note, c("ok", "ok", "flipped", "ok"))
    expect_equal(bed$n_fold, rep(5L, 4L))
    # The published columns still lead, and the weight is published as WGT
    # rather than appearing a second time under its source name.
    expect_equal(names(bed)[1:4], c("#chr", "start", "end", "ID"))
    expect_false(is_in("weight", names(bed)))
    # Naming one explicitly works the same way; naming a published one is an
    # error rather than a duplicate column.
    writeBed(tw, out, columns = "qc_note")
    expect_false(is_in("n_fold", names(read.delim(out, check.names = FALSE))))
    expect_error(
        writeBed(tw, out, columns = "weight"),
        "already writes"
    )
})

test_that("a weight matrix with no columns publishes a header-only file", {
    # A joint fit that retained no conditions leaves a 0-column weight
    # matrix, so the row contributes no (variant, trait) block at all. The
    # published file is then the schema with no rows -- the same shape an
    # empty result takes everywhere else -- rather than a file with no header
    # a reader could not parse.
    gr <- .variantIdsToGRanges(sprintf("chr1:%d:A:G", 100L * (1:4)), "v")
    w <- matrix(numeric(0), nrow = 4L, ncol = 0L)
    mcols(gr)$weight <- w
    tw <- TwasWeights(
        studyName = "S",
        context = "c",
        trait = "t",
        method = "lasso",
        entry = list(methods::new(
            "TwasWeightsRow",
            variants = gr,
            weights = w,
            weightStandardized = FALSE
        ))
    )
    expect_length(.twbWeightLabels(mcols(tw[[1L]])$weight), 0L)
    expect_null(.twbWeightsTable(tw, NULL, NULL))
    out <- withr::local_tempfile(fileext = ".bed")
    writeBed(tw, out)
    lines <- readLines(out)
    expect_length(lines, 1L)
    expect_equal(
        str_split_1(lines[[1L]], "\t")[1:4],
        c(
            "#chr",
            "start",
            "end",
            "ID"
        )
    )
})

# ---- THETA: the infinitesimal half of WGT -----------------------------------

.thetaFmr <- function(theta = NULL, nV = 12) {
    vids <- sprintf("chr1:%d:A:T", seq_len(nV) * 1000)
    fit <- .thetaFit(nV, theta)
    tl <- buildTopLoci(
        fit = fit,
        csTables = NULL,
        variantNames = vids,
        method = "susie_inf"
    )
    QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "g",
        method = "susie_inf",
        entry = list(FineMappingRow(
            variantIds = vids,
            susieFit = fit,
            topLoci = tl
        ))
    )
}

test_that("THETA round-trips through both published formats", {
    th <- rnorm(12, sd = 0.05)
    fmr <- .thetaFmr(th)
    want <- as.data.frame(topLoci(fmr, signalCutoff = 0))$theta
    expect_equal(sort(want), sort(th), tolerance = 1e-10)

    bed <- withr::local_tempfile(fileext = ".bed")
    writeBed(fmr, bed)
    expect_true(is_in("THETA", strsplit(readLines(bed, n = 1), "\t")[[1]]))
    gotBed <- as.data.frame(topLoci(
        readFineMappingResult(bed, weightStandardized = TRUE),
        signalCutoff = 0
    ))$theta
    expect_equal(sort(gotBed), sort(want), tolerance = 1e-5)

    vcf <- withr::local_tempfile(fileext = ".vcf.bgz")
    writeGwasVcf(fmr, vcf)
    hdr <- grep("^##FORMAT", readLines(gzfile(vcf)), value = TRUE)
    expect_true(any(grepl("ID=THETA,", hdr)))
    gotVcf <- as.data.frame(topLoci(
        readFineMappingResult(vcf, weightStandardized = TRUE),
        signalCutoff = 0
    ))$theta
    expect_equal(sort(gotVcf), sort(want), tolerance = 1e-5)
})

test_that("THETA is omitted for collections that never recorded one", {
    # Objects written before theta was published carry no posterior_theta, so
    # the view yields NA and the drop-empty rule keeps it out of the file --
    # rather than inventing a column of zeros that claims the split is known.
    # The example now records theta, so an older object is made from it.
    data(qtlFineMappingExample)
    old <- endoapply(qtlFineMappingExample, function(gr) {
        S4Vectors::mcols(gr)$posterior_theta <- NULL
        gr
    })
    d <- as.data.frame(topLoci(old, signalCutoff = 0))
    expect_true(all(is.na(d$theta)))
    bed <- withr::local_tempfile(fileext = ".bed")
    writeBed(old[1], bed)
    expect_false(is_in("THETA", strsplit(readLines(bed, n = 1), "\t")[[1]]))
})

test_that("data_source falls back per row, and to NA for a mixed collection", {
    rows <- map(1:2, function(i) {
        TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:3) + i),
            weights = c(0.1, 0.2, 0.3)
        )
    })
    ids <- sprintf("chr1:%d:A:G", 100L * (1:3) + 1L)
    ldm <- LdMatrix(`dimnames<-`(diag(3), list(ids, ids)))
    tw <- TwasWeights(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "lasso"),
        entry = rows,
        ldReference = list(ldm, NULL)
    )
    expect_identical(
        pecotmr:::.pubDataSourceOf(tw),
        c("sumstat", "individual")
    )
    expect_identical(pecotmr:::.pubDataSource(tw), NA_character_)
    expect_identical(pecotmr:::.pubDataSource(tw[1L]), "sumstat")
})

# =============================================================================
# A GWAS study fine-mapped block by block
# =============================================================================

test_that("writeGwasVcf joins a GWAS fit's blocks into one sample", {
    skip_if_not_installed("VariantAnnotation")
    out <- tempfile(fileext = ".vcf")
    writeGwasVcf(.rw_blockedGwasFmr(), out)
    sample <- .rw_sampleLines(out)
    expect_length(sample, 1L)
    expect_match(sample, "TotalVariants=6", fixed = TRUE)
    expect_match(sample, "BlockId=chr1_1_5000|chr1_5001_20000", fixed = TRUE)
    expect_match(sample, "V=0.1;0.2|0.3;0.4", fixed = TRUE)
    expect_match(sample, "CsIndex=1|2", fixed = TRUE)
    expect_match(sample, "RegionGiven=TRUE|TRUE", fixed = TRUE)
    # Credible sets are numbered within a block, so each variant says which.
    vcf <- VariantAnnotation::readVcf(out)
    bid <- unlist(VariantAnnotation::geno(vcf)$BID[, 1L])
    expect_identical(
        unname(bid),
        rep(c("chr1_1_5000", "chr1_5001_20000"), each = 3L)
    )
})

test_that("writeGwasVcf splitByBlock writes one file per LD block", {
    skip_if_not_installed("VariantAnnotation")
    out <- tempfile(fileext = ".vcf")
    paths <- writeGwasVcf(.rw_blockedGwasFmr(), out, splitByBlock = TRUE)
    expect_identical(
        basename(paths),
        str_c(
            basename(tools::file_path_sans_ext(out)),
            c(".chr1_1_5000", ".chr1_5001_20000"),
            ".vcf"
        )
    )
    second <- .rw_sampleLines(paths[[2L]])
    expect_match(second, "BlockId=chr1_5001_20000,", fixed = TRUE)
    expect_match(second, "V=0.3;0.4", fixed = TRUE)
    vcf <- VariantAnnotation::readVcf(paths[[2L]])
    expect_equal(nrow(vcf), 3L)
    expect_false("BID" %in% names(VariantAnnotation::geno(vcf)))
})

test_that("writeGwasVcf refuses to join blocks that disagree on a fit field", {
    skip_if_not_installed("VariantAnnotation")
    x <- .rw_blockedGwasFmr(dataType = c("gwas", "eqtl"))
    expect_error(
        writeGwasVcf(x, tempfile(fileext = ".vcf")),
        "disagree on `data_type`"
    )
    expect_length(
        writeGwasVcf(x, tempfile(fileext = ".vcf"), splitByBlock = TRUE),
        2L
    )
})

test_that("writeGwasVcf splitByBlock is for GWAS results only", {
    skip_if_not_installed("VariantAnnotation")
    expect_error(
        writeGwasVcf(
            .makeMultiTupleQtlFmr(),
            tempfile(fileext = ".vcf"),
            splitByBlock = TRUE
        ),
        "a QtlFineMappingResult fit is one row"
    )
})

test_that("writeBed gives a GWAS result's rows their own windows", {
    data(gwasFineMappingExample)
    out <- tempfile(fileext = ".bed")
    writeBed(gwasFineMappingExample, out, signalCutoff = 0)
    tbl <- utils::read.delim(out, check.names = FALSE)
    expect_true(all(!is.na(tbl$region_id)))
    # Two blocks of one fit each keep their own window.
    blocked <- .rw_blockedGwasFmr()
    out2 <- tempfile(fileext = ".bed")
    writeBed(blocked, out2, signalCutoff = 0)
    tbl2 <- utils::read.delim(out2, check.names = FALSE)
    windows <- unique(tbl2[, c("block_id", "region_id", "region_given")])
    expect_identical(windows$block_id, blocked$blockId)
    # The blocks were given their windows, so those are published, given.
    expect_identical(windows$region_id, c("chr1_1_5000", "chr1_5001_20000"))
    expect_identical(windows$region_given, c(TRUE, TRUE))
})

test_that("writeBed records whether each block's range was given", {
    x <- .rw_mixedGwasFmr()
    out <- tempfile(fileext = ".bed")
    writeBed(x, out, signalCutoff = 0)
    tbl <- utils::read.delim(out, check.names = FALSE)
    windows <- unique(tbl[, c("block_id", "region_id", "region_given")])
    # The second block was given none, so its window is its span.
    expect_identical(windows$region_id, c("chr1_1_5000", "chr1_10100_10300"))
    expect_identical(windows$region_given, c(TRUE, FALSE))
    back <- readFineMappingResult(out, weightStandardized = TRUE)
    expect_identical(back$regionGiven, c(TRUE, FALSE))
    expect_identical(
        as.character(genomicRegion(back)),
        as.character(genomicRegion(x))
    )
    # A file that does not record it reads back as computed: a window read
    # from a file says nothing of where it came from.
    old <- tempfile(fileext = ".bed")
    utils::write.table(
        tbl[, setdiff(names(tbl), "region_given")],
        old,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    expect_identical(
        readFineMappingResult(old, weightStandardized = TRUE)$regionGiven,
        c(FALSE, FALSE)
    )
})

test_that("writeGwasVcf(GwasSumStats) defaults to a study split across rows", {
    skip_if_not_installed("VariantAnnotation")
    data(gwasSumStatsS4Example)
    g <- gwasSumStatsS4Example
    split <- GwasSumStats(
        studyName = "trait1",
        entry = list(sumStats(g)),
        genome = "hg19",
        ldReference = ldReference(g)[[1L]],
        region = GenomicRanges::GRanges(
            "chr22",
            IRanges::IRanges(c(1L, 15200001L), c(15200000L, 51304566L))
        )
    )
    expect_equal(nrow(split), 2L)
    out <- tempfile(fileext = ".vcf")
    writeGwasVcf(split, out)
    expect_equal(nrow(VariantAnnotation::readVcf(out)), length(sumStats(g)))
})

test_that("writeGwasVcf splitByContext keeps each context's fits in one file", {
    skip_if_not_installed("VariantAnnotation")
    # Two traits in one context decorate to one path; both fits belong in
    # that file, rather than the second overwriting the first.
    fmr <- .makeMultiTupleQtlFmr()
    S4Vectors::mcols(fmr)$context <- c("brain", "brain")
    paths <- writeGwasVcf(
        fmr,
        tempfile(fileext = ".vcf"),
        splitByContext = TRUE
    )
    expect_length(paths, 1L)
    header <- VariantAnnotation::scanVcfHeader(paths[[1L]])
    expect_length(VariantAnnotation::samples(header), 2L)
})

test_that("writeGwasVcf joins blocks that carry no marginal effects", {
    skip_if_not_installed("VariantAnnotation")
    # The fits' top loci hold no marginal beta / se / z, so the joined
    # sample has no ES / SE / EZ to write.
    bare <- function(offset, V, csIndex) {
        tl <- .sc_makeTopLoci(3, offset = offset)
        tl <- tl[, setdiff(
            names(tl),
            c(
                "marginal_beta",
                "marginal_se",
                "marginal_z",
                "marginal_p"
            )
        )]
        FineMappingRow(
            variantIds = tl$variant_id,
            susieFit = list(V = V, sets = list(cs_index = csIndex)),
            topLoci = tl
        )
    }
    x <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(bare(0L, c(0.1, 0.2), 1L), bare(10000L, c(0.3, 0.4), 2L))
    )
    out <- tempfile(fileext = ".vcf")
    writeGwasVcf(x, out)
    vcf <- VariantAnnotation::readVcf(out)
    expect_equal(nrow(vcf), 6L)
    expect_false("ES" %in% names(VariantAnnotation::geno(vcf)))
    expect_true("PIP" %in% names(VariantAnnotation::geno(vcf)))
})

# =============================================================================
# The genome build in the file header, and fit coverage
# =============================================================================

test_that("a published table names its genome build above the header", {
    data(qtlFineMappingExample, envir = environment())
    x <- qtlFineMappingExample
    GenomeInfoDb::genome(x) <- "hg38"
    path <- tempfile(fileext = ".bed")
    on.exit(unlink(path), add = TRUE)
    writeBed(x, path)
    lines <- readLines(path, n = 2L)
    expect_equal(lines[[1L]], "##genome=hg38")
    expect_match(lines[[2L]], "^#chr\t")
    # With no build named there is no line to write.
    GenomeInfoDb::genome(x) <- NA_character_
    plain <- tempfile(fileext = ".bed")
    on.exit(unlink(plain), add = TRUE)
    writeBed(x, plain)
    expect_match(readLines(plain, n = 1L), "^#chr\t")
})

test_that("a published VCF names its build as reference and contig assembly", {
    skip_if_not_installed("VariantAnnotation")
    data(qtlFineMappingExample, envir = environment())
    x <- qtlFineMappingExample
    GenomeInfoDb::genome(x) <- "hg38"
    path <- tempfile(fileext = ".vcf")
    on.exit(unlink(path), add = TRUE)
    writeGwasVcf(x, path, context = "context_1")
    lines <- readLines(path)
    expect_true(is_in("##reference=hg38", lines))
    contigs <- lines[str_detect(lines, "^##contig=")]
    expect_gt(length(contigs), 0L)
    expect_true(all(str_detect(contigs, "assembly=hg38")))
})

test_that("writeGwasVcf() warns about the narrowed rows it writes, only", {
    skip_if_not_installed("VariantAnnotation")
    data(qtlFineMappingLbfExample, envir = environment())
    x <- qtlFineMappingLbfExample
    narrowed <- endoapply(x, utils::head, 50L)
    path <- tempfile(fileext = ".vcf")
    on.exit(unlink(path), add = TRUE)
    expect_warning(
        writeGwasVcf(narrowed, path),
        "writeGwasVcf\\(\\): 2 row\\(s\\) hold fewer variants"
    )
    # Rows a selector leaves out are not written, so not reported.
    mixed <- narrowed
    mixed[[2L]] <- x[[2L]]
    expect_no_warning(
        writeGwasVcf(mixed, path, context = mixed$context[[2L]])
    )
})
