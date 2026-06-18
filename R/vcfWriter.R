#' @include allGenerics.R
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
           ...) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")

  # Resolve the single (study, context, trait, method) row to write.
  if (is.null(study) || is.null(context) || is.null(trait) || is.null(method)) {
    if (nrow(x) != 1L) {
      stop("This FineMappingResult has ", nrow(x), " entries. ",
           "Pass `study`, `context`, `trait`, and `method` to select one.")
    }
    study   <- as.character(x$study)[1L]
    context <- as.character(x$context)[1L]
    trait   <- as.character(x$trait)[1L]
    method  <- as.character(x$method)[1L]
  }
  entry <- getFineMappingResult(x, study, context, trait, method)
  sampleName <- sampleName %||% sprintf("%s|%s|%s|%s",
                                       study, context, trait, method)
  tl <- getTopLoci(entry)
  if (nrow(tl) == 0) stop("FineMappingEntry has no topLoci to write")

  parsed <- parseVariantId(tl$variant_id)
  nSnps <- nrow(parsed)

  geno <- list()
  genoHeaderRows <- character(0)
  genoNumber <- character(0)
  genoType <- character(0)
  genoDesc <- character(0)

  # PIP
  pipCol <- resolvePipColumn(tl)
  if (!is.null(pipCol)) {
    geno[["PIP"]] <- matrix(tl[[pipCol]], nSnps)
    genoHeaderRows <- c(genoHeaderRows, "PIP")
    genoNumber <- c(genoNumber, "A")
    genoType <- c(genoType, "Float")
    genoDesc <- c(genoDesc, "Posterior inclusion probability")
  }

  # CS
  csCol <- grep("^cs_index", colnames(tl), value = TRUE)
  if (length(csCol) > 0) {
    geno[["CS"]] <- matrix(as.integer(tl[[csCol[1]]]), nSnps)
    genoHeaderRows <- c(genoHeaderRows, "CS")
    genoNumber <- c(genoNumber, "A")
    genoType <- c(genoType, "Integer")
    genoDesc <- c(genoDesc, "Credible set index (0 = not in any CS)")
  }

  # Effect size / SE if available
  if ("beta" %in% colnames(tl)) {
    geno[["ES"]] <- matrix(tl$beta, nSnps)
    genoHeaderRows <- c(genoHeaderRows, "ES")
    genoNumber <- c(genoNumber, "A")
    genoType <- c(genoType, "Float")
    genoDesc <- c(genoDesc, "Effect size estimate relative to the alternative allele")
  }
  if ("se" %in% colnames(tl)) {
    geno[["SE"]] <- matrix(tl$se, nSnps)
    genoHeaderRows <- c(genoHeaderRows, "SE")
    genoNumber <- c(genoNumber, "A")
    genoType <- c(genoType, "Float")
    genoDesc <- c(genoDesc, "Standard error of effect size estimate")
  }
  if ("z" %in% colnames(tl)) {
    pval <- 2 * pnorm(-abs(tl$z))
    geno[["LP"]] <- matrix(-log10(pval), nSnps)
    genoHeaderRows <- c(genoHeaderRows, "LP")
    genoNumber <- c(genoNumber, "A")
    genoType <- c(genoType, "Float")
    genoDesc <- c(genoDesc, "-log10 p-value for effect estimate")
  }

  genoHeader <- DataFrame(
    Number = genoNumber,
    Type = genoType,
    Description = genoDesc,
    row.names = genoHeaderRows)

  .writeVcfImpl(
    chrom = parsed$chrom,
    pos = parsed$pos,
    ref = parsed$A2,
    alt = parsed$A1,
    snpIds = tl$variant_id,
    geno = geno,
    genoHeader = genoHeader,
    sampleName = sampleName,
    outputPath = outputPath)
})

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
