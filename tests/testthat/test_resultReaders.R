# ---------------------------------------------------------------------------
# Published result readers: readFineMappingResult() / readTwasWeights()
#
# The reader and the writer are each other's only real test. Each direction
# had only ever been checked against a stored fixture before, which is how the
# two came to disagree about a reserved key while both "passed" -- so most of
# what is here is a round trip, and the comparisons are against the ORIGINAL
# object rather than against a file.
# ---------------------------------------------------------------------------

# @noRd
.rdr_tmp <- function(name) file.path(tempdir(), name)

# The ways a round trip legitimately differs from its input, so a test can say
# "everything else is identical" and mean it:
#   dataSource / standardized  ADDED on the way back (one is derived by the
#     writer, the other supplied to the reader)
#   cs_<n>  the published form is an integer index, and "susie_0" -- the
#     not-captured sentinel -- publishes as NA, which is the canonical form
#     every reader of the column already treats as "in no credible set"
#   z_supplied  resolved from NA to FALSE: a file publishing ES/SE and no EZ
#     says the z was NOT supplied by the source, which is what the fixture's
#     NA meant before the flag existed
# @noRd
.RDR_NORMALIZED <- c(
    "dataSource",
    "standardized",
    "cs_95",
    "cs_70",
    "cs_50",
    "z_supplied"
)

# Which columns of two per-variant tables disagree, aligned on variant id.
#
# Aligned rather than compared in place because a published file is sorted
# genomically and the stored table is not, and restricted to the variants
# `b` carries because an export is topLoci-FILTERED by default -- the
# question is whether what was published came back intact, not whether a
# filtered file holds every variant of the fit.
# @noRd
.rdr_diff <- function(a, b, tolerance = 1e-6) {
    a <- a[is_in(a$variant_id, b$variant_id), ]
    b <- b[match(a$variant_id, b$variant_id), ]
    shared <- intersect(names(a), names(b))
    bad <- character(0)
    for (nm in shared) {
        va <- a[[nm]]
        vb <- b[[nm]]
        ok <- if (is.numeric(va) && is.numeric(vb)) {
            isTRUE(all.equal(va, vb, tolerance = tolerance))
        } else {
            identical(as.character(va), as.character(vb))
        }
        if (!ok) {
            bad <- c(bad, nm)
        }
    }
    setdiff(bad, .RDR_NORMALIZED)
}

# A minimal published table written by hand, as a user converting another
# tool's output would produce it.
# @noRd
.rdr_writeTable <- function(df, name) {
    path <- .rdr_tmp(name)
    utils::write.table(
        df,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    path
}

# The bundled fine-mapping example, published once as a VCF and read back
# once. Reading a published VCF is the slow step of a round trip, so the tests
# that need only "the example after a VCF round trip" share this copy instead
# of each writing and reading the same file again.
# @noRd
.rdr_fmVcf <- .rdr_tmp("shared_fm.vcf")
if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
    data(qtlFineMappingExample, envir = environment())
    writeGwasVcf(qtlFineMappingExample, .rdr_fmVcf, context = "context_1")
    .rdr_fmVcfBack <- readFineMappingResult(
        .rdr_fmVcf,
        weightStandardized = TRUE
    )
}

# One three-variant fit per context: small enough to publish and read back in
# a fraction of the example's time, for the tests that edit a published VCF's
# header or need only some VCF to read.
# @noRd
.rdr_smallFmr <- function(contexts = "c1") {
    n <- length(contexts)
    QtlFineMappingResult(
        studyName = rep("s1", n),
        context = contexts,
        trait = rep("t1", n),
        method = rep("susie", n),
        entry = rep(list(.sc_makeFineMappingRow(3)), n)
    )
}

# ---- round trips ----------------------------------------------------------

test_that("a fine-mapping result round-trips through a published BED", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("rt_fm.bed")
    writeBed(qtlFineMappingExample, path, columns = "all")
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(nrow(back), nrow(qtlFineMappingExample))
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(qtlFineMappingExample)),
            as.data.frame(topLoci(back))
        ),
        character(0)
    )
    # The marginal view is the other half of the same table and has its own
    # published fields (ES/SE/EZ/LP), so it is checked separately.
    expect_equal(
        .rdr_diff(
            as.data.frame(marginalEffects(qtlFineMappingExample)),
            as.data.frame(marginalEffects(back))
        ),
        character(0)
    )
    # The analysis WINDOW, without which genomicRegion() silently degrades to
    # the span of the variants that survived QC.
    expect_equal(
        genomicRegion(back),
        genomicRegion(qtlFineMappingExample)
    )
})

test_that("a fine-mapping result round-trips through a published VCF", {
    skip_if_not_installed("VariantAnnotation")
    data(qtlFineMappingExample, envir = environment())
    back <- .rdr_fmVcfBack
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(qtlFineMappingExample)),
            as.data.frame(topLoci(back)),
            tolerance = 1e-5
        ),
        character(0)
    )
    expect_equal(genomicRegion(back), genomicRegion(qtlFineMappingExample))
})

test_that("the two formats agree with each other", {
    skip_if_not_installed("VariantAnnotation")
    # The point of one shared vocabulary: whichever format a result was
    # published in, it reads back the same. A `region_start` that is 0-based
    # in the table and 1-based in the VCF header was caught by exactly this.
    data(qtlFineMappingExample, envir = environment())
    bed <- .rdr_tmp("agree.bed")
    writeBed(qtlFineMappingExample, bed, columns = "all")
    fromBed <- readFineMappingResult(bed, weightStandardized = TRUE)
    fromVcf <- .rdr_fmVcfBack
    expect_equal(genomicRegion(fromBed), genomicRegion(fromVcf))
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(fromBed)),
            as.data.frame(topLoci(fromVcf)),
            tolerance = 1e-5
        ),
        character(0)
    )
})

test_that("a published VCF restores the per-effect fit the header carried", {
    # `V` and `sets$cs_index` are per EFFECT, so no per-variant column can
    # hold them -- the `##SAMPLE` block is why they survive at all, and
    # colocalization's credible-set filter reads them off the fit.
    skip_if_not_installed("VariantAnnotation")
    data(qtlFineMappingExample, envir = environment())
    fit <- susieFit(.rdr_fmVcfBack)
    original <- susieFit(qtlFineMappingExample)
    expect_equal(as.numeric(fit$V), as.numeric(original$V), tolerance = 1e-8)
    expect_equal(
        as.integer(fit$sets$cs_index),
        as.integer(original$sets$cs_index)
    )
    # A BED has no header to carry them, so there the fit is honestly absent
    # rather than half-invented.
    bed <- .rdr_tmp("fit.bed")
    writeBed(qtlFineMappingExample, bed)
    expect_null(susieFit(readFineMappingResult(bed, weightStandardized = TRUE)))
})

test_that("TWAS weights round-trip through a published BED", {
    data(twasWeightsExample, envir = environment())
    path <- .rdr_tmp("rt_w.bed")
    expect_message(
        writeBed(twasWeightsExample, path),
        "fine-mapping-derived"
    )
    back <- readTwasWeights(path)
    # Fine-mapping-derived rows are skipped on export, so fewer rows come
    # back than went out -- that is the writer's rule, not a reader loss.
    expect_lt(nrow(back), nrow(twasWeightsExample))
    md <- S4Vectors::mcols(back)
    for (i in seq_len(nrow(back))) {
        key <- list(
            studyName = as.character(md$study[[i]]),
            context = as.character(md$context[[i]]),
            trait = as.character(md$trait[[i]]),
            method = as.character(md$method[[i]])
        )
        was <- exec(resolveWeights, twasWeightsExample, !!!key)
        now <- exec(resolveWeights, back, !!!key)
        expect_setequal(now$variantIds, was$variantIds)
        expect_equal(
            as.numeric(now$weights)[order(now$variantIds)],
            as.numeric(was$weights)[order(was$variantIds)],
            tolerance = 1e-10
        )
    }
    expect_equal(
        genomicRegion(back),
        genomicRegion(twasWeightsExample)[
            match(
                as.character(md$method),
                as.character(S4Vectors::mcols(twasWeightsExample)$method)
            )
        ]
    )
})

# ---- multi-row tables -----------------------------------------------------

test_that("a merged table reads back as a multi-row collection", {
    data(qtlFineMappingExample, qtlFineMappingLbfExample, envir = environment())
    combined <- combineFineMappingResults(
        qtlFineMappingExample,
        qtlFineMappingLbfExample
    )
    path <- .rdr_tmp("merged.bed")
    writeBed(combined, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_equal(nrow(back), nrow(combined))
    md <- S4Vectors::mcols(back)
    expect_setequal(
        as.character(md$trait),
        as.character(S4Vectors::mcols(combined)$trait)
    )
})

test_that("a per-result value that disagrees within one tuple is an error", {
    # Two results concatenated under one identity. Picking either window
    # would silently publish one fit's region as the other's.
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.1),
        study = c("s", "s"),
        context = c("c", "c"),
        trait = c("t", "t"),
        method = c("susie", "susie"),
        region_start = c(0L, 500L),
        region_end = c(1000L, 1500L),
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "clash.bed")
    expect_error(
        readFineMappingResult(path),
        "`region_start` takes 2 different values"
    )
})

# ---- supplying what the file does not record -------------------------------

test_that("identity arguments fill in a file that did not record them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "bare.bed")
    back <- readFineMappingResult(
        path,
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie"
    )
    expect_s4_class(back, "QtlFineMappingResult")
    md <- S4Vectors::mcols(back)
    expect_equal(as.character(md$study), "s")
    expect_equal(as.character(md$method), "susie")
    expect_equal(as.numeric(topLoci(back)$pip), 0.9)
})

test_that("an identity argument that disagrees with the file is an error", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "recorded",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "disagree.bed")
    expect_error(
        readFineMappingResult(path, studyName = "supplied"),
        "was given as 'supplied' but the file records 'recorded'"
    )
    # Dropping the argument uses the file's value.
    expect_equal(
        as.character(S4Vectors::mcols(readFineMappingResult(path))$study),
        "recorded"
    )
})

test_that("a result keyed by axes the file lacks names them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        context = "c",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noaxes.bed")
    expect_error(
        readFineMappingResult(path),
        "does not record study, trait, method.*`studyName`, `trait`, `method`"
    )
})

# ---- permissiveness -------------------------------------------------------

test_that("a table carrying only PIP is enough for a fine-mapping result", {
    # TWAS needs only the weight, MR beta/se/PIP, colocalization the LBF. A
    # file converted from another tool has some and not others, so the reader
    # requires only what the class cannot exist without.
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.05),
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "pip_only.bed")
    back <- readFineMappingResult(path)
    tl <- topLoci(back, signalCutoff = 0)
    expect_equal(length(tl), 2L)
    expect_equal(as.numeric(S4Vectors::mcols(tl)$pip), c(0.9, 0.05))
    # Absent columns are absent, not NA-filled: nothing invents a beta.
    expect_true(all(is.na(marginalEffects(back)$beta)))
})

test_that("a table carrying only a weight is enough for TwasWeights", {
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        WGT = c(0.3, -0.2),
        study = "s",
        context = "c",
        trait = "t",
        method = "lasso",
        weight_standardized = TRUE,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wgt_only.bed")
    back <- readTwasWeights(path)
    w <- resolveWeights(back)
    expect_equal(as.numeric(w$weights), c(0.3, -0.2))
    expect_true(weightStandardized(back))
})

test_that("ids are rebuilt from coordinates when ID does not encode them", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "rs123",
        REF = "A",
        ALT = "G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "rsid.bed")
    back <- readFineMappingResult(path)
    expect_equal(variantIds(back), "chr1:100:A:G")
})

test_that("a table with no way to identify its variants says so", {
    df <- data.frame(
        ID = "rs123",
        PIP = 0.9,
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noid.bed")
    expect_error(
        readFineMappingResult(path, studyName = "s"),
        "must identify its variants"
    )
})

test_that("a weights file read as a fine-mapping result points back", {
    data(twasWeightsExample, envir = environment())
    wt <- .rdr_tmp("x_wt.bed")
    suppressMessages(writeBed(twasWeightsExample, wt))
    expect_error(
        readFineMappingResult(wt, weightStandardized = TRUE),
        "no `PIP` column.*readTwasWeights"
    )
})

test_that("a fine-mapping table reads as weights, because WGT is the weight", {
    # Not a category error in this direction: a susie posterior mean IS the
    # TWAS coefficient, which is why it publishes as WGT in both schemas.
    data(qtlFineMappingExample, envir = environment())
    fm <- .rdr_tmp("x_fm.bed")
    writeBed(qtlFineMappingExample, fm)
    asWeights <- readTwasWeights(fm, weightStandardized = TRUE)
    expect_s4_class(asWeights, "TwasWeights")
    expect_equal(
        as.numeric(resolveWeights(asWeights)$weights),
        as.numeric(
            S4Vectors::mcols(
                topLoci(qtlFineMappingExample)
            )$posterior_effect
        )[order(variantIds(
            readFineMappingResult(fm, weightStandardized = TRUE)
        ))],
        tolerance = 1e-8
    )
})

# ---- the weight scale -----------------------------------------------------

test_that("a published weight without its scale has to be told", {
    data(twasWeightsExample, envir = environment())
    path <- .rdr_tmp("scale.bed")
    suppressMessages(writeBed(twasWeightsExample, path))
    # This fixture records the flag, so it reads without one.
    expect_s4_class(readTwasWeights(path), "TwasWeights")
    df <- utils::read.delim(path, check.names = FALSE)
    df$weight_standardized <- NA
    stripped <- .rdr_writeTable(df, "scale_stripped.bed")
    expect_error(
        readTwasWeights(stripped),
        "publishes WGT but not `weight_standardized`"
    )
    expect_s4_class(
        readTwasWeights(stripped, weightStandardized = FALSE),
        "TwasWeights"
    )
})

# ---- opt-in Wakefield log Bayes factors ------------------------------------

test_that("priorVariance derives log Bayes factors only when asked", {
    df <- data.frame(
        `#chr` = c("chr1", "chr1"),
        start = c(99L, 199L),
        end = c(100L, 200L),
        ID = c("chr1:100:A:G", "chr1:200:C:T"),
        PIP = c(0.9, 0.05),
        ES = c(0.3, -0.1),
        SE = c(0.05, 0.05),
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wakefield.bed")
    plain <- readFineMappingResult(path)
    expect_true(all(is.na(
        S4Vectors::mcols(
            topLoci(plain, signalCutoff = 0)
        )$logBF
    )))
    derived <- readFineMappingResult(path, priorVariance = 0.04)
    lbf <- S4Vectors::mcols(topLoci(derived, signalCutoff = 0))$logBF
    expect_equal(
        lbf,
        .approxBayesFactor(df$ES / df$SE, df$SE, W = 0.04)$lbf,
        tolerance = 1e-10
    )
    # Never silently over published values.
    df$LBF <- c(5, 0)
    withLbf <- .rdr_writeTable(df, "wakefield_has_lbf.bed")
    expect_error(
        readFineMappingResult(withLbf, priorVariance = 0.04),
        "already carries LBF"
    )
})

test_that("priorVariance without a z or a standard error says so", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "wakefield_bare.bed")
    expect_error(
        readFineMappingResult(path, priorVariance = 0.04),
        "publishes no z \\(EZ, or ES/SE\\)"
    )
})

# ---- EZ, and what it means --------------------------------------------------

test_that("EZ round-trips as a SUPPLIED z and is not invented", {
    # The spec defines EZ as the z "provided if it was used to derive ES and
    # SE". A z pecotmr computed as ES/SE must not come back claiming to be
    # the source's, or a reader would believe ES/SE were derived.
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        ES = 0.3,
        SE = 0.05,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    noEz <- readFineMappingResult(.rdr_writeTable(df, "ez_no.bed"))
    me <- marginalEffects(noEz)
    expect_equal(as.numeric(me$z), 6, tolerance = 1e-8)
    expect_false(any(as.logical(me$z_supplied)))
    df$EZ <- 5.5
    withEz <- readFineMappingResult(.rdr_writeTable(df, "ez_yes.bed"))
    me <- marginalEffects(withEz)
    expect_equal(as.numeric(me$z), 5.5)
    expect_true(all(as.logical(me$z_supplied)))
})

# ---- class choice, format choice, region ------------------------------------

test_that("a result with no context or trait is a GWAS result", {
    data(gwasFineMappingExample, envir = environment())
    path <- .rdr_tmp("gwas.bed")
    writeBed(gwasFineMappingExample, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_s4_class(back, "GwasFineMappingResult")
    expect_equal(
        .rdr_diff(
            as.data.frame(topLoci(gwasFineMappingExample)),
            as.data.frame(topLoci(back))
        ),
        character(0)
    )
})

test_that("the format follows the extension unless it is named", {
    expect_equal(pecotmr:::.prdFormat("x/y.bed.gz", "auto"), "bed")
    expect_equal(pecotmr:::.prdFormat("x/y.tsv", "auto"), "bed")
    expect_equal(pecotmr:::.prdFormat("x/y.vcf.gz", "auto"), "vcf")
    expect_equal(pecotmr:::.prdFormat("x/y.bcf", "auto"), "vcf")
    # A file named something else is read as whatever it is said to be.
    expect_equal(pecotmr:::.prdFormat("x/y.txt", "vcf"), "vcf")
})

test_that("a tabix-indexed table can be read by region", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("region.bed.gz")
    writeBed(qtlFineMappingExample, path)
    all <- readFineMappingResult(path, weightStandardized = TRUE)
    full <- topLoci(all, signalCutoff = 0)
    window <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(full))[[1L]],
        IRanges::IRanges(
            start = GenomicRanges::start(full)[[1L]],
            end = GenomicRanges::start(full)[[3L]]
        )
    )
    part <- readFineMappingResult(
        path,
        which = window,
        weightStandardized = TRUE
    )
    expect_lt(length(variantIds(part)), length(variantIds(all)))
    expect_true(all(is_in(variantIds(part), variantIds(all))))
})

test_that("an un-indexed file ignores a region loudly", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    path <- .rdr_writeTable(df, "noindex.bed")
    expect_warning(
        readFineMappingResult(
            path,
            which = GenomicRanges::GRanges("chr1:1-1000")
        ),
        "not bgzipped \\+ tabix-indexed"
    )
    # A region given in R is a GRanges; a string is refused, not read.
    expect_error(
        readFineMappingResult(path, which = "chr1:1-1000"),
        "`which` must be a GRanges"
    )
})

test_that("a missing file names itself", {
    expect_error(
        readFineMappingResult(.rdr_tmp("absent_xyz.bed")),
        "does not exist"
    )
})

# ---- the ##SAMPLE key table -------------------------------------------------

test_that("every ##SAMPLE key is alphanumeric", {
    # Rsamtools' header parser splits the `<...>` on `,(?=[[:alnum:]]+=)`, so
    # an underscored key is not recognised as a key at all: it glues onto the
    # previous value and every field after it on the line is lost. Six of
    # these were unreadable for that reason.
    keys <- pecotmr:::.PUB_SAMPLE_KEYS
    expect_true(all(str_detect(keys, "^[A-Za-z0-9]+$")))
    expect_equal(length(unique(keys)), length(keys))
})

test_that("a per-fit field with no ##SAMPLE key stops the writer", {
    expect_error(
        pecotmr:::.pubSampleKeys(c("study", "brand_new_field")),
        "brand_new_field.*no\\s+##SAMPLE key"
    )
})

test_that("the ##SAMPLE block survives a VariantAnnotation read", {
    skip_if_not_installed("VariantAnnotation")
    meta <- pecotmr:::.prdVcfSampleMeta(VariantAnnotation::readVcf(.rdr_fmVcf))
    # Named in the shared vocabulary, not the header's CamelCase: the reader
    # inverts the key table so everything downstream sees one set of names.
    expect_true(all(is_in(
        c("study", "context", "trait", "method", "region_start", "cs_index"),
        colnames(meta)
    )))
    expect_equal(meta$method, "susie")
})

test_that("a header-only file reads back as an empty collection", {
    # An export of a result whose every entry was screened out writes its
    # schema and no rows. That is a legitimate file, so it reads as an empty
    # collection rather than erroring about identity axes it has no rows to
    # record.
    header <- str_flatten(pecotmr:::.pubSchema(character(0)), "\t")
    path <- .rdr_tmp("empty.bed")
    writeLines(header, path)
    back <- readFineMappingResult(path)
    expect_s4_class(back, "QtlFineMappingResult")
    expect_equal(nrow(back), 0L)
    weights <- readTwasWeights(path)
    expect_s4_class(weights, "TwasWeights")
    expect_equal(nrow(weights), 0L)
})

test_that("an empty result survives an export and a read", {
    data(qtlFineMappingExample, envir = environment())
    path <- .rdr_tmp("empty_rt.bed")
    # A signal cutoff above 1 retains nothing, so the table is vacuous.
    writeBed(qtlFineMappingExample, path, signalCutoff = 2)
    expect_equal(nrow(readFineMappingResult(path)), 0L)
})

# ---------------------------------------------------------------------------
# The VCF reader's edge paths
# ---------------------------------------------------------------------------

# A published VCF's bgzipped, tabix-indexed copy, so it can be read by region.
# @noRd
.rdr_indexedVcf <- function(plain) {
    gz <- Rsamtools::bgzip(plain, dest = str_c(plain, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "vcf")
    gz
}

test_that("a tabix-indexed VCF can be read by region", {
    skip_if_not_installed("VariantAnnotation")
    gz <- .rdr_indexedVcf(.rdr_fmVcf)
    # Read whole, the indexed copy holds what the shared round trip read.
    all <- .rdr_fmVcfBack
    full <- topLoci(all, signalCutoff = 0)
    window <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(full))[[1L]],
        IRanges::IRanges(
            start = GenomicRanges::start(full)[[1L]],
            end = GenomicRanges::start(full)[[5L]]
        )
    )
    part <- readFineMappingResult(
        gz,
        which = window,
        weightStandardized = TRUE
    )
    expect_lt(length(variantIds(part)), length(variantIds(all)))
    expect_true(all(is_in(variantIds(part), variantIds(all))))
    # The window still carries the analysis region from the ##SAMPLE header,
    # which is per-fit and does not shrink with the variants read.
    expect_equal(genomicRegion(part), genomicRegion(all))
})

test_that("a region naming no seqname in the VCF reads as empty", {
    skip_if_not_installed("VariantAnnotation")
    gz <- .rdr_indexedVcf(.rdr_fmVcf)
    back <- readFineMappingResult(
        gz,
        which = GenomicRanges::GRanges("chrZZ:1-1000"),
        weightStandardized = TRUE,
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie"
    )
    expect_equal(nrow(back), 0L)
})

test_that("an un-indexed VCF ignores a region loudly", {
    skip_if_not_installed("VariantAnnotation")
    path <- .rdr_tmp("plain.vcf")
    writeGwasVcf(.rdr_smallFmr(), path)
    expect_warning(
        readFineMappingResult(
            path,
            which = GenomicRanges::GRanges("chr1:1-1000"),
            weightStandardized = TRUE
        ),
        "not bgzipped \\+ tabix-indexed"
    )
})

test_that("a missing VCF names itself before VariantAnnotation is asked", {
    skip_if_not_installed("VariantAnnotation")
    expect_error(
        readFineMappingResult(.rdr_tmp("absent_xyz.vcf")),
        "does not exist"
    )
})

test_that("sampleSelect picks fits out of a multi-sample VCF", {
    skip_if_not_installed("VariantAnnotation")
    # Two rows -> two samples in one file, which is GWAS-VCF's own shape.
    path <- .rdr_tmp("multi.vcf")
    writeGwasVcf(.rdr_smallFmr(c("c1", "c2")), path)
    both <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_equal(nrow(both), 2L)
    samples <- VariantAnnotation::samples(
        VariantAnnotation::scanVcfHeader(path)
    )
    one <- readFineMappingResult(
        path,
        sampleSelect = samples[[1L]],
        weightStandardized = TRUE
    )
    expect_equal(nrow(one), 1L)
    expect_error(
        readFineMappingResult(path, sampleSelect = "not-a-sample"),
        "sampleSelect names not-a-sample, which the VCF does not have"
    )
})

test_that("a VCF with no ##SAMPLE block falls back to the sample name", {
    skip_if_not_installed("VariantAnnotation")
    path <- .rdr_tmp("nosample.vcf")
    writeGwasVcf(.rdr_smallFmr(), path)
    stripped <- .rdr_tmp("nosample_stripped.vcf")
    lines <- readLines(path)
    writeLines(lines[!str_starts(lines, "##SAMPLE=")], stripped)
    # Without the header the file records no method either, so the reader asks
    # for the one axis it cannot invent.
    expect_error(
        readFineMappingResult(stripped, weightStandardized = TRUE),
        "does not record method"
    )
    back <- readFineMappingResult(
        stripped,
        method = "susie",
        weightStandardized = TRUE
    )
    # Splitting the sample name into a tuple would be guesswork; taking it
    # whole keeps each sample a result of its own.
    expect_s4_class(back, "GwasFineMappingResult")
    expect_match(studyName(back), "\\|")
    # No header means no partial fit to restore either.
    expect_null(susieFit(back))
})

# ---------------------------------------------------------------------------
# Published-table edge paths
# ---------------------------------------------------------------------------

test_that("a table with start but no end still locates its variants", {
    # `end` is the 1-based position of a half-open single-base record, and
    # `start + 1` says the same thing -- a table carrying only `start` is
    # still readable.
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        check.names = FALSE
    )
    back <- readFineMappingResult(.rdr_writeTable(df, "startonly.bed"))
    expect_equal(variantIds(back), "chr1:100:A:G")
})

test_that("an identity-less header-only file reads as an empty GWAS result", {
    # No context or trait column at all, so the GWAS branch of the empty path.
    writeLines(
        str_flatten(c("#chr", "start", "end", "ID", "PIP"), "\t"),
        .rdr_tmp("empty_gwas.bed")
    )
    back <- readFineMappingResult(.rdr_tmp("empty_gwas.bed"))
    expect_s4_class(back, "GwasFineMappingResult")
    expect_equal(nrow(back), 0L)
})

test_that("a weights file read with no WGT column points at the other reader", {
    df <- data.frame(
        `#chr` = "chr1",
        start = 99L,
        end = 100L,
        ID = "chr1:100:A:G",
        PIP = 0.9,
        check.names = FALSE
    )
    expect_error(
        readTwasWeights(.rdr_writeTable(df, "nowgt.bed")),
        "no `WGT` column.*readFineMappingResult"
    )
})

test_that("a published window is read on UCSC sequence names", {
    df <- data.frame(
        region_id = "22_100_200",
        region_start = 99L,
        region_end = 200L
    )
    gr <- .prdRowRegion(df, 1L, "x")
    expect_identical(as.character(GenomicRanges::seqnames(gr)), "chr22")
    expect_identical(GenomicRanges::start(gr), 100L)
})

test_that("a region id with no underscores falls back to the variants' chrom", {
    expect_equal(
        .prdRegionChrom(
            data.frame(
                region_id = "noUnderscores",
                `#chr` = "chr7",
                check.names = FALSE
            ),
            1L,
            "x"
        ),
        "chr7"
    )
    # Nothing to go on at all.
    expect_true(is.na(.prdRegionChrom(data.frame(a = 1), 1L, "x")))
})

test_that("a packed header value unpacks, and an absent one is NULL", {
    expect_equal(.prdUnpackNumeric(list("1;2;3")), c(1, 2, 3))
    expect_null(.prdUnpackNumeric(NULL))
    expect_null(.prdUnpackNumeric(list(NA_character_)))
    expect_null(.prdUnpackNumeric(list("")))
})

# ---------------------------------------------------------------------------
# Headers sparser than the writer's own
# ---------------------------------------------------------------------------

# A published VCF whose `##SAMPLE` block is replaced wholesale, keeping the
# sample `ID` so the block still belongs to the sample. A file written by an
# older pecotmr -- or by anything else following the GWAS-VCF spec -- carries
# fewer keys than today's writer, and the reader has to work from whatever is
# there rather than requiring the full set.
# @noRd
.rdr_reheaderedVcf <- function(name, keys) {
    path <- .rdr_tmp(str_c(name, "_src.vcf"))
    writeGwasVcf(.rdr_smallFmr(), path)
    lines <- readLines(path)
    at <- which(str_starts(lines, "##SAMPLE="))
    expect_length(at, 1L)
    id <- str_match(lines[[at]], "ID=([^,>]+)")[, 2L]
    lines[[at]] <- str_c("##SAMPLE=<ID=", id, ",", keys, ">")
    out <- .rdr_tmp(str_c(name, ".vcf"))
    writeLines(lines, out)
    list(path = out, sample = id)
}

test_that("a ##SAMPLE block with no identity keys reads as one study", {
    skip_if_not_installed("VariantAnnotation")
    f <- .rdr_reheaderedVcf("noident", "TotalVariants=3")
    # Nothing in the header maps to an identity axis, so the sample name is
    # the only study there is -- and method, which cannot be invented, has to
    # be supplied.
    expect_error(
        readFineMappingResult(f$path, weightStandardized = TRUE),
        "does not record method"
    )
    back <- readFineMappingResult(
        f$path,
        method = "susie",
        weightStandardized = TRUE
    )
    expect_s4_class(back, "GwasFineMappingResult")
    expect_equal(studyName(back), f$sample)
    # A header carrying neither V nor CsIndex restores no partial fit.
    expect_null(susieFit(back))
})

test_that("a ##SAMPLE block with no RegionStart still reads", {
    skip_if_not_installed("VariantAnnotation")
    f <- .rdr_reheaderedVcf(
        "noregion",
        "Study=s1,Context=c1,Trait=t1,Method=susie"
    )
    back <- readFineMappingResult(f$path, weightStandardized = TRUE)
    expect_equal(as.character(back$study), "s1")
    expect_equal(as.character(back$method), "susie")
    # The window is the span of the variants the file holds, because the
    # header recorded no analysis region to prefer over it.
    rg <- genomicRegion(back)
    expect_equal(start(rg), min(IRanges::start(back[[1L]])))
    expect_null(susieFit(back))
})

# A header or data line trimmed to the eight fixed VCF columns.
# @noRd
.rdr_dropSampleColumns <- function(line) {
    if (str_starts(line, "##")) {
        return(line)
    }
    str_flatten(str_split_1(line, "\t")[1:8], "\t")
}

test_that("a VCF with no samples reads as an empty result", {
    skip_if_not_installed("VariantAnnotation")
    lines <- readLines(.rdr_fmVcf)
    keep <- !str_starts(lines, "##SAMPLE=") & !str_starts(lines, "##FORMAT=")
    # A sites-only VCF: the eight fixed columns and no sample at all, which
    # is a legal VCF and carries no results.
    sites <- map_chr(lines[keep], .rdr_dropSampleColumns)
    out <- .rdr_tmp("sites.vcf")
    writeLines(sites, out)
    # Naming the real problem: the earlier behaviour was to hand back an
    # empty table, which failed further in with "no `PIP` column".
    expect_error(
        readFineMappingResult(out, method = "susie"),
        "the VCF has no samples"
    )
    expect_error(readTwasWeights(out), "the VCF has no samples")
})

test_that("each sample of a VCF reads back only its own variants", {
    skip_if_not_installed("VariantAnnotation")
    # Two fits over different variants: the records are their union, and a
    # fit leaves the cells of the other's variants empty.
    fmr <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(2, offset = 1000L)
        )
    )
    path <- .rdr_tmp("disjoint.vcf")
    writeGwasVcf(fmr, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_equal(unname(lengths(back)), c(3L, 2L))
    expect_false(anyNA(pip(back[1L])))
    expect_false(anyNA(pip(back[2L])))
})

test_that("a BID naming a block the header does not list is refused", {
    row <- data.frame(block_id = "b1|b2")
    expect_identical(
        pecotmr:::.prdBlockPositions(row, c("b2", NA, "b1")),
        c(2L, NA, 1L)
    )
    expect_null(pecotmr:::.prdBlockPositions(data.frame(block_id = "b1"), "b1"))
    expect_error(
        pecotmr:::.prdBlockPositions(row, "b3"),
        "names block\\(s\\) b3 that the sample's header does not list"
    )
})

test_that("a header field one sample lacks reads back as NA for it", {
    skip_if_not_installed("VariantAnnotation")
    fmr <- GwasFineMappingResult(
        studyName = c("g1", "g2"),
        method = c("susie", "susie"),
        entry = list(
            .sc_makeFineMappingRow(3),
            .sc_makeFineMappingRow(2, offset = 1000L)
        ),
        dataType = c("gwas", NA)
    )
    path <- .rdr_tmp("onefield.vcf")
    writeGwasVcf(fmr, path)
    back <- readFineMappingResult(path, weightStandardized = TRUE)
    expect_identical(back$dataType, c("gwas", NA))
})

test_that(".prdBlockPiece reads one block's part of a header cell", {
    expect_null(pecotmr:::.prdBlockPiece(NULL, 1L))
    expect_null(pecotmr:::.prdBlockPiece(".", 1L))
    expect_null(pecotmr:::.prdBlockPiece("a|b", NA_integer_))
    expect_identical(pecotmr:::.prdBlockPiece("a|b", 2L), "b")
    # A single-block sample's cell is its one value whole.
    expect_identical(pecotmr:::.prdBlockPiece("a", 2L), "a")
})

# ---- the genome build -------------------------------------------------------

test_that("a result read back carries the build its file recorded", {
    x <- .rdr_smallFmr()
    GenomeInfoDb::genome(x) <- "hg38"
    bed <- .rdr_tmp("build_fm.bed")
    writeBed(x, bed, columns = "all")
    expect_equal(
        unique(unname(GenomeInfoDb::genome(
            readFineMappingResult(bed, weightStandardized = TRUE)
        ))),
        "hg38"
    )
    skip_if_not_installed("VariantAnnotation")
    vcf <- .rdr_tmp("build_fm.vcf")
    writeGwasVcf(x, vcf)
    expect_equal(
        unique(unname(GenomeInfoDb::genome(
            readFineMappingResult(vcf, weightStandardized = TRUE)
        ))),
        "hg38"
    )
})

test_that("a file that names no build reads back with none", {
    data(qtlFineMappingExample, envir = environment())
    x <- qtlFineMappingExample
    GenomeInfoDb::genome(x) <- NA_character_
    bed <- .rdr_tmp("nobuild_fm.bed")
    writeBed(x, bed, columns = "all")
    expect_true(all(is.na(GenomeInfoDb::genome(
        readFineMappingResult(bed, weightStandardized = TRUE)
    ))))
    expect_true(is.na(pecotmr:::.prdBedBuild(bed)))
})

test_that("TWAS weights read back with the build their file recorded", {
    data(twasWeightsExample, envir = environment())
    x <- twasWeightsExample
    GenomeInfoDb::genome(x) <- "hg19"
    bed <- .rdr_tmp("build_tw.bed")
    writeBed(x, bed)
    expect_equal(
        unique(unname(GenomeInfoDb::genome(readTwasWeights(bed)))),
        "hg19"
    )
})

test_that("a file of no lines, or of meta lines alone, still reads a build", {
    # Every line a `##` line is the count's other branch: there is no first
    # non-meta line to stop at.
    empty <- .rdr_tmp("empty_build.bed")
    writeLines(character(0), empty)
    expect_true(is.na(pecotmr:::.prdBedBuild(empty)))
    metaOnly <- .rdr_tmp("meta_only_build.bed")
    writeLines(c("##genome=hg38", "##source=pecotmr"), metaOnly)
    expect_equal(pecotmr:::.prdBedBuild(metaOnly), "hg38")
    expect_equal(pecotmr:::.prdLeadingCount(c(TRUE, FALSE, TRUE)), 1L)
})

# ===========================================================================
# readFineMappingResult: a joined sample splits back into its blocks
# ===========================================================================

test_that("readFineMappingResult splits a joined sample back into its blocks", {
    skip_if_not_installed("VariantAnnotation")
    x <- .rw_blockedGwasFmr()
    out <- tempfile(fileext = ".vcf")
    writeGwasVcf(x, out)
    back <- readFineMappingResult(out, weightStandardized = TRUE)
    expect_s4_class(back, "GwasFineMappingResult")
    expect_identical(back$blockId, x$blockId)
    expect_equal(unname(lengths(back)), c(3L, 3L))
    expect_equal(susieFit(back[2L])$V, c(0.3, 0.4))
    expect_identical(susieFit(back[2L])$sets$cs_index, 2L)
    expect_equal(pip(back[2L]), pip(x[2L]))
    expect_identical(
        as.character(genomicRegion(back)),
        as.character(genomicRegion(x))
    )
    expect_identical(back$regionGiven, c(TRUE, TRUE))
    # A block whose range is its variants' span reads back as computed.
    mixed <- .rw_mixedGwasFmr()
    out2 <- tempfile(fileext = ".vcf")
    writeGwasVcf(mixed, out2)
    expect_match(
        .rw_sampleLines(out2),
        "RegionGiven=TRUE|FALSE",
        fixed = TRUE
    )
    back2 <- readFineMappingResult(out2, weightStandardized = TRUE)
    expect_identical(back2$regionGiven, c(TRUE, FALSE))
})
