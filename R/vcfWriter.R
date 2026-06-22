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
    if (!requireNamespace("VariantAnnotation", quietly = TRUE))
      stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")

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
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")

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

  # Body of the VCF is exclusively marginal univariate effects — no
  # posterior output. By design the fine-mapping write-out emits the
  # marginal sumstats so consumers can run their own downstream
  # analysis (coloc, TWAS, etc.) on a uniform per-variant table.
  marginal <- getMarginalEffects(entry)
  if (nrow(marginal) == 0)
    stop("writeSumstatsVcf: entry [", sn, "] has no variants to write")

  nSnps <- nrow(marginal)
  geno <- list()
  hdrRows <- character(0); hdrNum <- character(0)
  hdrType <- character(0); hdrDesc <- character(0)
  addGeno <- function(name, vec, type, desc) {
    geno[[name]] <<- matrix(vec, nSnps)
    hdrRows <<- c(hdrRows, name); hdrNum <<- c(hdrNum, "A")
    hdrType <<- c(hdrType, type); hdrDesc <<- c(hdrDesc, desc)
  }
  if (any(!is.na(marginal$beta)))
    addGeno("ES", marginal$beta, "Float",
            "Marginal univariate effect-size estimate (effect allele)")
  if (any(!is.na(marginal$se)))
    addGeno("SE", marginal$se, "Float",
            "Standard error of the marginal effect-size estimate")
  if (any(!is.na(marginal$p))) {
    lp <- ifelse(is.na(marginal$p) | marginal$p <= 0,
                 NA_real_, -log10(marginal$p))
    addGeno("LP", lp, "Float",
            "-log10 p-value of the marginal univariate effect")
  }
  if (any(!is.na(marginal$N)))
    addGeno("SS", as.integer(marginal$N), "Integer", "Sample size")
  if (any(!is.na(marginal$MAF)))
    addGeno("AF", marginal$MAF, "Float", "Minor allele frequency")

  genoHeader <- DataFrame(
    Number = hdrNum, Type = hdrType, Description = hdrDesc,
    row.names = hdrRows)

  .writeVcfImpl(
    chrom = marginal$chrom,
    pos = marginal$pos,
    ref = marginal$A2,
    alt = marginal$A1,
    snpIds = marginal$variant_id,
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
