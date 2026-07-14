#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame SimpleList mcols
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges DataFrameList
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @importFrom Rsamtools asBcf
#' @importFrom tools file_ext
NULL

#' @rdname writeSumstatsVcf
#' @export
setMethod("writeSumstatsVcf", signature("GwasSumStats"),
  function(x, outputPath, sampleName = NULL, study = NULL, ...) {
    # nocov start
    if (!requireNamespace("VariantAnnotation", quietly = TRUE))
      stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")
    # nocov end

    # Select which study to write (the new GwasSumStats can hold many).
    if (is.null(study)) {
      if (nrow(x) != 1L) {
        stop("This GwasSumStats has ", nrow(x),
             " studies. Pass `study = <name>` to select one.")
      }
      study <- as.character(x$study[[1L]])
    }
    ss <- getSumStats(x, study = study)
    mc <- mcols(ss)
    sampleName <- sampleName %||% study

    nSnps <- length(ss)
    geno <- list()
    if ("Z" %in% colnames(mc))
      geno[["ES"]] <- matrix(mc$Z, nSnps)
    if ("N" %in% colnames(mc))
      geno[["SS"]] <- matrix(as.integer(mc$N), nSnps)
    if ("MAF" %in% colnames(mc))
      geno[["AF"]] <- matrix(mc$MAF, nSnps)

    genoHeader <- DataFrame(
      Number = c("A", "A", "A"),
      Type = c("Float", "Integer", "Float"),
      Description = c(
        "Z-score of effect size estimate",
        "Sample size",
        "Minor allele frequency"),
      row.names = c("ES", "SS", "AF"))

    .writeVcfImpl(
      chrom = as.character(seqnames(ss)),
      pos = start(ss),
      ref = mc$A2,
      alt = mc$A1,
      snpIds = mc$SNP,
      geno = geno,
      genoHeader = genoHeader,
      sampleName = sampleName,
      outputPath = outputPath)
  })

#' @rdname writeSumstatsVcf
#' @export
setMethod("writeSumstatsVcf", signature("FineMappingResultBase"),
  function(x, outputPath, sampleName = NULL,
           study = NULL, context = NULL, trait = NULL, method = NULL,
           splitByContext = FALSE, splitByTrait = FALSE,
           ...) {
  # nocov start
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")
  # nocov end

  # Resolve the set of rows to write. With both selectors NULL and no
  # split flags, the collection must have exactly one row. Splitting
  # iterates over the unique values of the requested axis.
  rowSpecs <- .resolveFineMappingRows(
    x, study = study, context = context, trait = trait, method = method,
    splitByContext = splitByContext, splitByTrait = splitByTrait)
  out <- character(length(rowSpecs))
  for (i in seq_along(rowSpecs)) {
    spec <- rowSpecs[[i]]
    out[[i]] <- .writeFineMappingVcf(x, spec,
                                     outputPath = outputPath,
                                     sampleName = sampleName,
                                     splitByContext = splitByContext,
                                     splitByTrait   = splitByTrait)
  }
  invisible(out)
})

# Resolve which (study, context, trait, method) rows to write. Without
# the split flags this returns a single spec; with `splitByContext` or
# `splitByTrait` the collection's rows are walked and one spec is emitted
# per row (after applying any explicit selector filters).
# @noRd
.resolveFineMappingRows <- function(x, study, context, trait, method,
                                    splitByContext, splitByTrait) {
  hasContextSlot <- "context" %in% names(x)
  hasTraitSlot   <- "trait"   %in% names(x)
  rows <- seq_len(nrow(x))
  if (!is.null(study))   rows <- rows[as.character(x$study)[rows]   == study]
  if (hasContextSlot && !is.null(context))
    rows <- rows[as.character(x$context)[rows] == context]
  if (hasTraitSlot && !is.null(trait))
    rows <- rows[as.character(x$trait)[rows]   == trait]
  if (!is.null(method))  rows <- rows[as.character(x$method)[rows]  == method]
  if (length(rows) == 0L)
    stop("writeSumstatsVcf: no rows match the supplied selectors.")
  if (!isTRUE(splitByContext) && !isTRUE(splitByTrait)) {
    if (length(rows) != 1L)
      stop("This FineMappingResult has ", length(rows), " matching rows. ",
           "Pass `study`/`context`/`trait`/`method` to select one, or ",
           "set `splitByContext = TRUE` / `splitByTrait = TRUE` to emit ",
           "one file per row.")
    return(list(.rowSpec(x, rows[[1L]])))
  }
  lapply(rows, function(r) .rowSpec(x, r))
}

# Build a (study, context, trait, method) spec list for one row index.
# @noRd
.rowSpec <- function(x, r) {
  list(
    study   = as.character(x$study)[r],
    context = if ("context" %in% names(x)) as.character(x$context)[r]
              else NA_character_,
    trait   = if ("trait"   %in% names(x)) as.character(x$trait)[r]
              else NA_character_,
    method  = as.character(x$method)[r])
}

# Internal worker: write one (study, context, trait, method) tuple to a
# single VCF. When `splitByContext` / `splitByTrait` is in play the
# output path is decorated with the corresponding tag(s) so multiple
# files don't collide.
# @noRd
.writeFineMappingVcf <- function(x, spec, outputPath, sampleName,
                                 splitByContext, splitByTrait) {
  entry <- getFineMappingResult(x, spec$study, spec$context, spec$trait,
                                spec$method)
  finalPath <- .decorateOutputPath(outputPath, spec, splitByContext,
                                   splitByTrait)
  sn <- sampleName %||% sprintf("%s|%s|%s|%s",
                                 spec$study, spec$context %||% "_",
                                 spec$trait %||% "_", spec$method)

  # Per-variant body: the fine-mapping POSTERIOR (getTopLoci at signalCutoff 0 =
  # every fitted variant, carrying pip / CS membership / conditional effect),
  # joined to the MARGINAL univariate sumstats (getMarginalEffects) where those
  # exist. mvSuSiE / fSuSiE have no marginal sumstats, so the posterior drives
  # the variant set and supplies ES=conditional_effect + PIP + CS; univariate
  # susie adds marginal ES=beta / SE / LP / AF on top. This restores the legacy
  # create_vcf fields (ES / CS / PIP) the mv_susie cell wrote.
  post <- tryCatch(as.data.frame(getTopLoci(entry, signalCutoff = 0)),
                   error = function(e) NULL)
  marg <- tryCatch(as.data.frame(getMarginalEffects(entry)),
                   error = function(e) NULL)
  hasPost <- !is.null(post) && nrow(post) > 0L
  hasMarg <- !is.null(marg) && nrow(marg) > 0L
  if (hasPost) {
    base <- post
    m <- if (hasMarg) marg[match(base$variant_id, marg$variant_id), , drop = FALSE]
         else NULL
  } else if (hasMarg) {
    base <- marg; m <- marg
  } else {
    stop("writeSumstatsVcf: entry [", sn, "] has no variants to write")
  }
  nSnps <- nrow(base)
  col <- function(df, nm) if (!is.null(df) && nm %in% names(df)) df[[nm]]
                          else rep(NA, nSnps)

  geno <- list()
  hdrRows <- character(0); hdrNum <- character(0)
  hdrType <- character(0); hdrDesc <- character(0)
  addGeno <- function(name, vec, type, desc) {
    geno[[name]] <<- matrix(vec, nSnps)
    hdrRows <<- c(hdrRows, name); hdrNum <<- c(hdrNum, "A")
    hdrType <<- c(hdrType, type); hdrDesc <<- c(hdrDesc, desc)
  }
  # ES: posterior conditional effect (mvSuSiE/fSuSiE) when present, else the
  # marginal univariate beta (univariate susie).
  es <- col(base, "conditional_effect")
  if (all(is.na(es))) es <- col(m, "beta")
  if (any(!is.na(es)))
    addGeno("ES", es, "Float",
            "Effect size (posterior conditional effect, else marginal beta), effect allele")
  se <- col(m, "se")
  if (any(!is.na(se)))
    addGeno("SE", se, "Float", "Standard error of the marginal effect-size estimate")
  p <- col(m, "p")
  if (any(!is.na(p))) {
    lp <- ifelse(is.na(p) | p <= 0, NA_real_, -log10(p))
    addGeno("LP", lp, "Float", "-log10 p-value of the marginal univariate effect")
  }
  if (any(!is.na(col(m, "N"))))
    addGeno("SS", as.integer(col(m, "N")), "Integer", "Sample size")
  af <- col(base, "af"); if (all(is.na(af))) af <- col(m, "af")
  if (any(!is.na(af)))
    addGeno("AF", af, "Float", "Allele frequency (effect allele)")
  # Posterior fields, only when a posterior table is available.
  if (hasPost) {
    pip <- col(base, "pip")
    if (any(!is.na(pip)))
      addGeno("PIP", pip, "Float", "Posterior inclusion probability")
    lbf <- col(base, "logBF")
    if (any(!is.na(lbf)))
      addGeno("LBF", lbf, "Float", "Per-variant log Bayes factor (max single effect)")
    lfsr <- col(base, "lfsr")
    if (any(!is.na(lfsr)))
      addGeno("LFSR", lfsr, "Float", "Local false sign rate (per-condition posterior)")
    # Credible sets are DYNAMIC: pecotmr does not assume any fixed coverage, so
    # we emit a CS<coverage> (+ PUR<coverage>) field for every cs_<coverage>
    # column the pipeline actually produced (e.g. cs_95 -> CS95 / PUR95). The
    # 50/70/95 defaults are an xqtl-protocol pipeline choice, not baked in here.
    csCols <- grep("^cs_[0-9.]+$", names(base), value = TRUE)
    for (cc in csCols) {
      cov <- sub("^cs_", "", cc)
      idx <- suppressWarnings(as.integer(sub(".*_", "", as.character(base[[cc]]))))
      idx[is.na(idx)] <- 0L
      addGeno(paste0("CS", cov), idx, "Integer",
              sprintf("Credible-set index at %s%% coverage (0 = not captured)", cov))
      pc <- paste0(cc, "_purity")
      if (pc %in% names(base)) {
        pur <- suppressWarnings(as.numeric(base[[pc]]))
        if (any(!is.na(pur)))
          addGeno(paste0("PUR", cov), pur, "Float",
                  sprintf("Purity (min abs corr) of the %s%% credible set", cov))
      }
    }
    # Per-CS variant-level fullFit columns (within_cs_pip default scalar, and the
    # wide within_cs_pip_<lab> / cs_logbf_<lab> / cs_effect_<lab> / cs_effect_var_<lab>
    # sets when the topLoci was built with fullFit=TRUE). Emit each as an
    # uppercase FORMAT field.
    for (cc in grep("^(within_cs_pip|cs_logbf_|cs_effect_)", names(base), value = TRUE)) {
      v <- suppressWarnings(as.numeric(base[[cc]]))
      if (any(!is.na(v)))
        addGeno(toupper(cc), v, "Float",
                sprintf("Per-credible-set variant statistic (%s)", cc))
    }
  }

  genoHeader <- DataFrame(
    Number = hdrNum, Type = hdrType, Description = hdrDesc,
    row.names = hdrRows)

  .writeVcfImpl(
    chrom = col(base, "chrom"),
    pos = col(base, "pos"),
    ref = col(base, "A2"),
    alt = col(base, "A1"),
    snpIds = col(base, "variant_id"),
    geno = geno,
    genoHeader = genoHeader,
    sampleName = sn,
    outputPath = finalPath)
  finalPath
}

# Decorate `outputPath` with the spec's context / trait tags when split
# flags are set. Preserves the file extension. Examples:
#   "out.vcf" + (context="brain") -> "out.brain.vcf"
#   "out.vcf.bgz" + (context="brain", trait="ENSG1") -> "out.brain.ENSG1.vcf.bgz"
# @noRd
.decorateOutputPath <- function(outputPath, spec, splitByContext,
                                splitByTrait) {
  if (!isTRUE(splitByContext) && !isTRUE(splitByTrait)) return(outputPath)
  ext <- tolower(tools::file_ext(outputPath))
  composite <- ext == "bgz" || ext == "gz"
  base <- if (composite) {
    sub("\\.[^.]+\\.(bgz|gz)$", "", outputPath, ignore.case = TRUE)
  } else {
    tools::file_path_sans_ext(outputPath)
  }
  ext_keep <- substr(outputPath, nchar(base) + 1L, nchar(outputPath))
  tags <- character(0)
  if (isTRUE(splitByContext) &&
      !is.null(spec$context) && !is.na(spec$context) && nzchar(spec$context))
    tags <- c(tags, spec$context)
  if (isTRUE(splitByTrait) &&
      !is.null(spec$trait) && !is.na(spec$trait) && nzchar(spec$trait))
    tags <- c(tags, spec$trait)
  if (length(tags) == 0L) return(outputPath)
  paste0(base, ".", paste(tags, collapse = "."), ext_keep)
}

# Internal implementation shared by all methods
# @noRd
.writeVcfImpl <- function(chrom, pos, ref, alt, snpIds, geno, genoHeader,
                          sampleName, outputPath) {
  nSnps <- length(chrom)

  # Ensure chromosome names have "chr" prefix
  if (!all(grepl("^chr", chrom)))
    chrom <- paste0("chr", chrom)

  # Build GRanges for row ranges
  gr <- GRanges(
    chrom,
    IRanges(
      start = as.integer(pos),
      end = as.integer(pos) + pmax(nchar(ref), nchar(alt)) - 1L,
      names = snpIds))

  # Build VCF header
  coldata <- DataFrame(Samples = sampleName, row.names = sampleName)

  hdr <- VariantAnnotation::VCFHeader(
    header = DataFrameList(
      fileformat = DataFrame(
        Value = "VCFv4.2", row.names = "fileformat")),
    sample = sampleName)

  # Subset geno header to only fields present in geno
  genoHeader <- genoHeader[rownames(genoHeader) %in% names(geno), , drop = FALSE]
  VariantAnnotation::geno(hdr) <- genoHeader

  # Build VCF object
  genoSl <- SimpleList(geno)
  vcf <- VariantAnnotation::VCF(
    rowRanges = gr,
    colData = coldata,
    exptData = list(header = hdr),
    geno = genoSl)

  VariantAnnotation::ref(vcf) <- DNAStringSet(ref)
  VariantAnnotation::alt(vcf) <- DNAStringSetList(as.list(alt))
  VariantAnnotation::fixed(vcf)$FILTER <- "PASS"
  vcf <- sort(vcf)

  # Write based on output format
  # Note: VariantAnnotation::writeVcf appends ".bgz" to the path when
  # index = TRUE, so we must pass the path *without* the .bgz/.gz suffix.
  ext <- file_ext(outputPath)
  if (ext == "bcf") {
    # Write temporary bgzipped VCF, then convert to BCF
    tmpVcfStem <- tempfile(fileext = ".vcf")
    tmpVcfBgz <- paste0(tmpVcfStem, ".bgz")
    on.exit(unlink(c(tmpVcfBgz, paste0(tmpVcfBgz, ".tbi")),
                   force = TRUE), add = TRUE)
    VariantAnnotation::writeVcf(vcf, tmpVcfStem, index = TRUE)
    # asBcf appends ".bcf" to destination, so strip the extension
    bcfStem <- sub("\\.bcf$", "", outputPath)
    dict <- unique(chrom)
    asBcf(tmpVcfBgz, dictionary = dict,
                     destination = bcfStem)
  } else if (ext == "gz" || ext == "bgz") {
    # writeVcf will append .bgz, so strip it from the path
    vcfStem <- sub("\\.(bgz|gz)$", "", outputPath)
    VariantAnnotation::writeVcf(vcf, vcfStem, index = TRUE)
    # writeVcf always creates .bgz; rename if the user requested .gz
    actualPath <- paste0(vcfStem, ".bgz")
    if (actualPath != outputPath && file.exists(actualPath)) {
      file.rename(actualPath, outputPath)
      tbiActual <- paste0(actualPath, ".tbi")
      if (file.exists(tbiActual))
        file.rename(tbiActual, paste0(outputPath, ".tbi"))
    }
  } else {
    VariantAnnotation::writeVcf(vcf, outputPath)
  }

  invisible(outputPath)
}
