#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges DataFrameList
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @importFrom Rsamtools asBcf
#' @importFrom tools file_ext
#' @importFrom purrr map map_chr keep set_names list_flatten
NULL

#' @rdname writeSumstatsVcf
#' @export
setMethod(
    "writeSumstatsVcf",
    signature("GwasSumStats"),
    function(x, outputPath, sampleName = NULL, study = NULL, ...) {
        # nocov start
        if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
            stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")
        }
        # nocov end
        study <- .vcfResolveStudy(x, study)
        ss <- getSumStats(x, study = study)
        mc <- mcols(ss)
        .writeVcfImpl(
            chrom = as.character(seqnames(ss)),
            pos = start(ss),
            ref = mc$A2,
            alt = mc$A1,
            snpIds = mc$SNP,
            geno = .vcfSumstatsGeno(mc, length(ss)),
            genoHeader = .vcfSumstatsGenoHeader(),
            sampleName = sampleName %||% study,
            outputPath = outputPath
        )
    }
)

# Select which study to write (a GwasSumStats can hold many); a single-study
# collection defaults to its one study.
# @noRd
.vcfResolveStudy <- function(x, study) {
    if (!is.null(study)) {
        return(study)
    }
    if (nrow(x) != 1L) {
        stop(
            "This GwasSumStats has ",
            nrow(x),
            " studies. Pass `study = <name>` to select one."
        )
    }
    as.character(x$study[[1L]])
}

# Per-sample geno matrices (ES / SS / AF) for the fields present in the mcols.
# @noRd
.vcfSumstatsGeno <- function(mc, nSnps) {
    geno <- list()
    if ("Z" %in% colnames(mc)) {
        geno[["ES"]] <- matrix(mc$Z, nSnps)
    }
    if ("N" %in% colnames(mc)) {
        geno[["SS"]] <- matrix(as.integer(mc$N), nSnps)
    }
    if ("MAF" %in% colnames(mc)) {
        geno[["AF"]] <- matrix(mc$MAF, nSnps)
    }
    geno
}

# The fixed FORMAT header for the sumstats geno fields (ES / SS / AF).
# @noRd
.vcfSumstatsGenoHeader <- function() {
    DataFrame(
        Number = c("A", "A", "A"),
        Type = c("Float", "Integer", "Float"),
        Description = c(
            "Z-score of effect size estimate",
            "Sample size",
            "Minor allele frequency"
        ),
        row.names = c("ES", "SS", "AF")
    )
}

#' @rdname writeSumstatsVcf
#' @export
setMethod(
    "writeSumstatsVcf",
    signature("FineMappingResultBase"),
    function(
        x,
        outputPath,
        sampleName = NULL,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        splitByContext = FALSE,
        splitByTrait = FALSE,
        ...
    ) {
        # nocov start
        if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
            stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")
        }
        # nocov end

        # Resolve the set of rows to write. With both selectors NULL and no
        # split flags, the collection must have exactly one row. Splitting
        # iterates over the unique values of the requested axis.
        rowSpecs <- .resolveFineMappingRows(
            x,
            study = study,
            context = context,
            trait = trait,
            method = method,
            splitByContext = splitByContext,
            splitByTrait = splitByTrait
        )
        out <- character(length(rowSpecs))
        for (i in seq_along(rowSpecs)) {
            spec <- rowSpecs[[i]]
            out[[i]] <- .writeFineMappingVcf(
                x,
                spec,
                outputPath = outputPath,
                sampleName = sampleName,
                splitByContext = splitByContext,
                splitByTrait = splitByTrait
            )
        }
        invisible(out)
    }
)

# Resolve which (study, context, trait, method) rows to write. Without
# the split flags this returns a single spec; with `splitByContext` or
# `splitByTrait` the collection's rows are walked and one spec is emitted
# per row (after applying any explicit selector filters).
# @noRd
.resolveFineMappingRows <- function(
    x,
    study,
    context,
    trait,
    method,
    splitByContext,
    splitByTrait
) {
    hasContextSlot <- "context" %in% names(x)
    hasTraitSlot <- "trait" %in% names(x)
    rows <- seq_len(nrow(x))
    if (!is.null(study)) {
        rows <- rows[as.character(x$study)[rows] == study]
    }
    if (hasContextSlot && !is.null(context)) {
        rows <- rows[as.character(x$context)[rows] == context]
    }
    if (hasTraitSlot && !is.null(trait)) {
        rows <- rows[as.character(x$trait)[rows] == trait]
    }
    if (!is.null(method)) {
        rows <- rows[as.character(x$method)[rows] == method]
    }
    if (length(rows) == 0L) {
        stop("writeSumstatsVcf: no rows match the supplied selectors.")
    }
    if (!isTRUE(splitByContext) && !isTRUE(splitByTrait)) {
        if (length(rows) != 1L) {
            stop(
                "This FineMappingResult has ",
                length(rows),
                " matching rows. ",
                "Pass `study`/`context`/`trait`/`method` to select one, or ",
                "set `splitByContext = TRUE` / `splitByTrait = TRUE` to emit ",
                "one file per row."
            )
        }
        return(list(.rowSpec(x, rows[[1L]])))
    }
    lapply(rows, function(r) .rowSpec(x, r))
}

# Build a (study, context, trait, method) spec list for one row index.
# @noRd
.rowSpec <- function(x, r) {
    list(
        study = as.character(x$study)[r],
        context = if ("context" %in% names(x)) {
            as.character(x$context)[r]
        } else {
            NA_character_
        },
        trait = if ("trait" %in% names(x)) {
            as.character(x$trait)[r]
        } else {
            NA_character_
        },
        method = as.character(x$method)[r]
    )
}

# Extract column `nm` from `df`, or an all-NA vector of length nSnps when
# absent.
# @noRd
.vcfCol <- function(df, nm, nSnps) {
    if (!is.null(df) && nm %in% names(df)) df[[nm]] else rep(NA, nSnps)
}

# One VCF FORMAT field spec: name, per-variant `values`, VCF `type`
# (Float/Integer), and header `desc`. Number is always "A" (per-allele).
# @noRd
.vcfSpec <- function(name, values, type, desc) {
    list(name = name, values = values, type = type, desc = desc)
}

# Core FORMAT fields shared by univariate + posterior VCFs: effect size
# (posterior conditional effect, else marginal beta), standard error, -log10
# p-value, sample size, and effect-allele frequency. Returns a list of specs;
# empty (all-NA) ones are dropped downstream.
# @noRd
.vcfCoreSpecs <- function(base, m, nSnps) {
    es <- .vcfCol(base, "conditional_effect", nSnps)
    if (all(is.na(es))) {
        es <- .vcfCol(m, "beta", nSnps)
    }
    af <- .vcfCol(base, "af", nSnps)
    if (all(is.na(af))) {
        af <- .vcfCol(m, "af", nSnps)
    }
    p <- .vcfCol(m, "p", nSnps)
    lp <- ifelse(is.na(p) | p <= 0, NA_real_, -log10(p))
    list(
        .vcfSpec(
            "ES",
            es,
            "Float",
            paste0(
                "Effect size (posterior conditional effect, ",
                "else marginal beta), effect allele"
            )
        ),
        .vcfSpec(
            "SE",
            .vcfCol(m, "se", nSnps),
            "Float",
            "Standard error of the marginal effect-size estimate"
        ),
        .vcfSpec(
            "LP",
            lp,
            "Float",
            "-log10 p-value of the marginal univariate effect"
        ),
        .vcfSpec(
            "SS",
            as.integer(.vcfCol(m, "N", nSnps)),
            "Integer",
            "Sample size"
        ),
        .vcfSpec("AF", af, "Float", "Allele frequency (effect allele)")
    )
}

# Credible-set FORMAT fields. pecotmr emits CS<coverage> (+ PUR<coverage> when
# a purity column exists) for EVERY cs_<coverage> column present, so coverage
# levels are data-driven rather than the 50/70/95 xqtl-protocol default.
# @noRd
.vcfCsSpecs <- function(base) {
    csCols <- grep("^cs_[0-9.]+$", names(base), value = TRUE)
    perCol <- map(csCols, function(cc) {
        cov <- sub("^cs_", "", cc)
        idx <- suppressWarnings(as.integer(sub(
            ".*_",
            "",
            as.character(base[[cc]])
        )))
        idx[is.na(idx)] <- 0L
        csSpec <- .vcfSpec(
            paste0("CS", cov),
            idx,
            "Integer",
            sprintf(
                "Credible-set index at %s%% coverage (0 = not captured)",
                cov
            )
        )
        pc <- paste0(cc, "_purity")
        if (!pc %in% names(base)) {
            return(list(csSpec))
        }
        pur <- suppressWarnings(as.numeric(base[[pc]]))
        purSpec <- .vcfSpec(
            paste0("PUR", cov),
            pur,
            "Float",
            sprintf(
                "Purity (min abs corr) of the %s%% credible set",
                cov
            )
        )
        list(csSpec, purSpec)
    })
    list_flatten(perCol)
}

# Per-CS variant-level fullFit FORMAT fields (within_cs_pip scalar plus the
# wide within_cs_pip_<lab> / cs_logbf_<lab> / cs_effect_<lab> sets emitted when
# the topLoci was built with fullFit = TRUE). Each column becomes an uppercase
# field.
# @noRd
.vcfFullFitSpecs <- function(base) {
    cols <- grep(
        "^(within_cs_pip|cs_logbf_|cs_effect_)",
        names(base),
        value = TRUE
    )
    map(cols, function(cc) {
        v <- suppressWarnings(as.numeric(base[[cc]]))
        .vcfSpec(
            toupper(cc),
            v,
            "Float",
            sprintf(
                "Per-credible-set variant statistic (%s)",
                cc
            )
        )
    })
}

# Posterior FORMAT fields, only when a fine-mapping posterior table exists:
# PIP, per-variant log Bayes factor, local false sign rate, plus the dynamic
# credible-set and fullFit column sets.
# @noRd
.vcfPosteriorSpecs <- function(base, nSnps) {
    fixed <- list(
        .vcfSpec(
            "PIP",
            .vcfCol(base, "pip", nSnps),
            "Float",
            "Posterior inclusion probability"
        ),
        .vcfSpec(
            "LBF",
            .vcfCol(base, "logBF", nSnps),
            "Float",
            "Per-variant log Bayes factor (max single effect)"
        ),
        .vcfSpec(
            "LFSR",
            .vcfCol(base, "lfsr", nSnps),
            "Float",
            "Local false sign rate (per-condition posterior)"
        )
    )
    c(fixed, .vcfCsSpecs(base), .vcfFullFitSpecs(base))
}

# Resolve the per-variant body: the fine-mapping POSTERIOR (getTopLoci at
# signalCutoff 0 = every fitted variant, carrying pip / CS membership /
# conditional effect) joined to the MARGINAL univariate sumstats
# (getMarginalEffects) where those exist. mvSuSiE / fSuSiE have no marginal
# sumstats, so the posterior drives the variant set; univariate susie adds
# marginal ES=beta / SE / LP / AF on top. Returns list(base, m, hasPost).
# @noRd
.vcfResolveBody <- function(entry, sn) {
    post <- tryCatch(
        as.data.frame(getTopLoci(entry, signalCutoff = 0)),
        error = function(e) NULL
    )
    marg <- tryCatch(
        as.data.frame(getMarginalEffects(entry)),
        error = function(e) NULL
    )
    hasPost <- !is.null(post) && nrow(post) > 0L
    hasMarg <- !is.null(marg) && nrow(marg) > 0L
    if (hasPost) {
        base <- post
        m <- if (hasMarg) {
            marg[match(base$variant_id, marg$variant_id), , drop = FALSE]
        } else {
            NULL
        }
    } else if (hasMarg) {
        base <- marg
        m <- marg
    } else {
        stop("writeSumstatsVcf: entry [", sn, "] has no variants to write")
    }
    list(base = base, m = m, hasPost = hasPost)
}

# Materialise the collected FORMAT specs into the geno matrix list + header
# DataFrame (dropping all-NA specs) and write the VCF via .writeVcfImpl.
# @noRd
.vcfWriteSpecs <- function(base, specs, sn, finalPath, nSnps) {
    specs <- keep(specs, function(s) any(!is.na(s$values)))
    geno <- set_names(
        map(specs, function(s) matrix(s$values, nSnps)),
        map_chr(specs, "name")
    )
    genoHeader <- DataFrame(
        Number = rep("A", length(specs)),
        Type = map_chr(specs, "type"),
        Description = map_chr(specs, "desc"),
        row.names = map_chr(specs, "name")
    )
    .writeVcfImpl(
        chrom = .vcfCol(base, "chrom", nSnps),
        pos = .vcfCol(base, "pos", nSnps),
        ref = .vcfCol(base, "A2", nSnps),
        alt = .vcfCol(base, "A1", nSnps),
        snpIds = .vcfCol(base, "variant_id", nSnps),
        geno = geno,
        genoHeader = genoHeader,
        sampleName = sn,
        outputPath = finalPath
    )
}

# Internal worker: write one (study, context, trait, method) tuple to a
# single VCF. When `splitByContext` / `splitByTrait` is in play the
# output path is decorated with the corresponding tag(s) so multiple
# files don't collide.
# @noRd
.writeFineMappingVcf <- function(
    x,
    spec,
    outputPath,
    sampleName,
    splitByContext,
    splitByTrait
) {
    entry <- getFineMappingResult(
        x,
        spec$study,
        spec$context,
        spec$trait,
        spec$method
    )
    finalPath <- .decorateOutputPath(
        outputPath,
        spec,
        splitByContext,
        splitByTrait
    )
    sn <- sampleName %||%
        sprintf(
            "%s|%s|%s|%s",
            spec$study,
            spec$context %||% "_",
            spec$trait %||% "_",
            spec$method
        )
    body <- .vcfResolveBody(entry, sn)
    nSnps <- nrow(body$base)
    specs <- c(
        .vcfCoreSpecs(body$base, body$m, nSnps),
        if (body$hasPost) .vcfPosteriorSpecs(body$base, nSnps) else list()
    )
    .vcfWriteSpecs(body$base, specs, sn, finalPath, nSnps)
    finalPath
}

# Decorate `outputPath` with the spec's context / trait tags when split flags
# are set. Preserves the file extension. Examples: "out.vcf" + (context="brain")
# -> "out.brain.vcf" "out.vcf.bgz" + (context="brain", trait="ENSG1") ->
# "out.brain.ENSG1.vcf.bgz"
# @noRd
.decorateOutputPath <- function(
    outputPath,
    spec,
    splitByContext,
    splitByTrait
) {
    if (!isTRUE(splitByContext) && !isTRUE(splitByTrait)) {
        return(outputPath)
    }
    ext <- tolower(tools::file_ext(outputPath))
    composite <- ext == "bgz" || ext == "gz"
    base <- if (composite) {
        sub("\\.[^.]+\\.(bgz|gz)$", "", outputPath, ignore.case = TRUE)
    } else {
        tools::file_path_sans_ext(outputPath)
    }
    ext_keep <- substr(outputPath, nchar(base) + 1L, nchar(outputPath))
    tags <- character(0)
    if (
        isTRUE(splitByContext) &&
            !is.null(spec$context) &&
            !is.na(spec$context) &&
            nzchar(spec$context)
    ) {
        tags <- c(tags, spec$context)
    }
    if (
        isTRUE(splitByTrait) &&
            !is.null(spec$trait) &&
            !is.na(spec$trait) &&
            nzchar(spec$trait)
    ) {
        tags <- c(tags, spec$trait)
    }
    if (length(tags) == 0L) {
        return(outputPath)
    }
    paste0(base, ".", paste(tags, collapse = "."), ext_keep)
}

# Internal implementation shared by all methods
# @noRd
.writeVcfImpl <- function(
    chrom,
    pos,
    ref,
    alt,
    snpIds,
    geno,
    genoHeader,
    sampleName,
    outputPath
) {
    # Ensure chromosome names have the "chr" prefix.
    if (!all(grepl("^chr", chrom))) {
        chrom <- paste0("chr", chrom)
    }
    gr <- .vcfBuildRanges(chrom, pos, ref, alt, snpIds)
    hdr <- .vcfBuildHeader(sampleName, genoHeader, geno)
    vcf <- .vcfBuildObject(gr, sampleName, hdr, geno, ref, alt)
    .vcfWriteByFormat(vcf, outputPath, chrom)
    invisible(outputPath)
}

# Per-variant rowRanges (end extended to the wider of ref / alt for indels).
# @noRd
.vcfBuildRanges <- function(chrom, pos, ref, alt, snpIds) {
    GRanges(
        chrom,
        IRanges(
            start = as.integer(pos),
            end = as.integer(pos) + pmax(nchar(ref), nchar(alt)) - 1L,
            names = snpIds
        )
    )
}

# VCF header (VCFv4.2 + the geno header subset to the fields actually present).
# @noRd
.vcfBuildHeader <- function(sampleName, genoHeader, geno) {
    hdr <- VariantAnnotation::VCFHeader(
        header = DataFrameList(
            fileformat = DataFrame(Value = "VCFv4.2", row.names = "fileformat")
        ),
        sample = sampleName
    )
    VariantAnnotation::geno(hdr) <- genoHeader[
        rownames(genoHeader) %in% names(geno),
        ,
        drop = FALSE
    ]
    hdr
}

# Assemble + finalize the VCF object (ref / alt / FILTER, sorted).
# @noRd
.vcfBuildObject <- function(gr, sampleName, hdr, geno, ref, alt) {
    vcf <- VariantAnnotation::VCF(
        rowRanges = gr,
        colData = DataFrame(Samples = sampleName, row.names = sampleName),
        exptData = list(header = hdr),
        geno = SimpleList(geno)
    )
    VariantAnnotation::ref(vcf) <- DNAStringSet(ref)
    VariantAnnotation::alt(vcf) <- DNAStringSetList(as.list(alt))
    VariantAnnotation::fixed(vcf)$FILTER <- "PASS"
    sort(vcf)
}

# Write the VCF in the format implied by the output extension. writeVcf appends
# ".bgz" when index = TRUE, so paths are passed WITHOUT the .bgz/.gz suffix.
# @noRd
.vcfWriteByFormat <- function(vcf, outputPath, chrom) {
    ext <- file_ext(outputPath)
    if (ext == "bcf") {
        .vcfWriteBcf(vcf, outputPath, chrom)
    } else if (ext == "gz" || ext == "bgz") {
        .vcfWriteBgz(vcf, outputPath)
    } else {
        VariantAnnotation::writeVcf(vcf, outputPath)
    }
    invisible(NULL)
}

# BCF path: write a temporary bgzipped VCF, then convert to BCF via asBcf.
# @noRd
.vcfWriteBcf <- function(vcf, outputPath, chrom) {
    tmpVcfStem <- tempfile(fileext = ".vcf")
    tmpVcfBgz <- paste0(tmpVcfStem, ".bgz")
    on.exit(
        unlink(c(tmpVcfBgz, paste0(tmpVcfBgz, ".tbi")), force = TRUE),
        add = TRUE
    )
    VariantAnnotation::writeVcf(vcf, tmpVcfStem, index = TRUE)
    # asBcf appends ".bcf" to destination, so strip the extension.
    asBcf(
        tmpVcfBgz,
        dictionary = unique(chrom),
        destination = sub("\\.bcf$", "", outputPath)
    )
}

# bgzip / gzip path: writeVcf always creates .bgz; rename to .gz on request.
# @noRd
.vcfWriteBgz <- function(vcf, outputPath) {
    vcfStem <- sub("\\.(bgz|gz)$", "", outputPath)
    VariantAnnotation::writeVcf(vcf, vcfStem, index = TRUE)
    actualPath <- paste0(vcfStem, ".bgz")
    if (actualPath != outputPath && file.exists(actualPath)) {
        file.rename(actualPath, outputPath)
        tbiActual <- paste0(actualPath, ".tbi")
        if (file.exists(tbiActual)) {
            file.rename(tbiActual, paste0(outputPath, ".tbi"))
        }
    }
}
