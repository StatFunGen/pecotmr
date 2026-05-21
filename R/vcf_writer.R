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
setMethod("writeSumstatsVcf", signature("GWASSumStats"), function(x, output_path, sample_name = NULL, ...) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")

  ss <- x@sumstats
  mc <- mcols(ss)
  sample_name <- sample_name %||% x@trait_name

  # Build GENO fields from GWASSumStats metadata
  n_snps <- length(ss)
  geno <- list()
  if ("Z" %in% colnames(mc))
    geno[["ES"]] <- matrix(mc$Z, n_snps)
  if ("N" %in% colnames(mc))
    geno[["SS"]] <- matrix(as.integer(mc$N), n_snps)
  if ("MAF" %in% colnames(mc))
    geno[["AF"]] <- matrix(mc$MAF, n_snps)

  geno_header <- DataFrame(
    Number = c("A", "A", "A"),
    Type = c("Float", "Integer", "Float"),
    Description = c(
      "Z-score of effect size estimate",
      "Sample size",
      "Minor allele frequency"),
    row.names = c("ES", "SS", "AF"))

  .write_vcf_impl(
    chrom = as.character(seqnames(ss)),
    pos = start(ss),
    ref = mc$A2,
    alt = mc$A1,
    snp_ids = mc$SNP,
    geno = geno,
    geno_header = geno_header,
    sample_name = sample_name,
    output_path = output_path)
})

#' @rdname writeSumstatsVcf
#' @export
setMethod("writeSumstatsVcf", signature("FineMappingResult"), function(x, output_path, sample_name = NULL, ...) {
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
    stop("Package 'VariantAnnotation' is required for writeSumstatsVcf")

  sample_name <- sample_name %||% x@method
  tl <- x@top_loci
  if (nrow(tl) == 0) stop("FineMappingResult has no top_loci to write")

  parsed <- parse_variant_id(tl$variant_id)
  n_snps <- nrow(parsed)

  geno <- list()
  geno_header_rows <- character(0)
  geno_number <- character(0)
  geno_type <- character(0)
  geno_desc <- character(0)

  # PIP
  pip_col <- resolve_pip_column(tl)
  if (!is.null(pip_col)) {
    geno[["PIP"]] <- matrix(tl[[pip_col]], n_snps)
    geno_header_rows <- c(geno_header_rows, "PIP")
    geno_number <- c(geno_number, "A")
    geno_type <- c(geno_type, "Float")
    geno_desc <- c(geno_desc, "Posterior inclusion probability")
  }

  # CS
  cs_col <- grep("^cs_index", colnames(tl), value = TRUE)
  if (length(cs_col) > 0) {
    geno[["CS"]] <- matrix(as.integer(tl[[cs_col[1]]]), n_snps)
    geno_header_rows <- c(geno_header_rows, "CS")
    geno_number <- c(geno_number, "A")
    geno_type <- c(geno_type, "Integer")
    geno_desc <- c(geno_desc, "Credible set index (0 = not in any CS)")
  }

  # Effect size / SE if available
  if ("beta" %in% colnames(tl)) {
    geno[["ES"]] <- matrix(tl$beta, n_snps)
    geno_header_rows <- c(geno_header_rows, "ES")
    geno_number <- c(geno_number, "A")
    geno_type <- c(geno_type, "Float")
    geno_desc <- c(geno_desc, "Effect size estimate relative to the alternative allele")
  }
  if ("se" %in% colnames(tl)) {
    geno[["SE"]] <- matrix(tl$se, n_snps)
    geno_header_rows <- c(geno_header_rows, "SE")
    geno_number <- c(geno_number, "A")
    geno_type <- c(geno_type, "Float")
    geno_desc <- c(geno_desc, "Standard error of effect size estimate")
  }
  if ("z" %in% colnames(tl)) {
    pval <- 2 * pnorm(-abs(tl$z))
    geno[["LP"]] <- matrix(-log10(pval), n_snps)
    geno_header_rows <- c(geno_header_rows, "LP")
    geno_number <- c(geno_number, "A")
    geno_type <- c(geno_type, "Float")
    geno_desc <- c(geno_desc, "-log10 p-value for effect estimate")
  }

  geno_header <- DataFrame(
    Number = geno_number,
    Type = geno_type,
    Description = geno_desc,
    row.names = geno_header_rows)

  .write_vcf_impl(
    chrom = parsed$chrom,
    pos = parsed$pos,
    ref = parsed$A2,
    alt = parsed$A1,
    snp_ids = tl$variant_id,
    geno = geno,
    geno_header = geno_header,
    sample_name = sample_name,
    output_path = output_path)
})

# Internal implementation shared by all methods
# @noRd
.write_vcf_impl <- function(chrom, pos, ref, alt, snp_ids, geno, geno_header,
                             sample_name, output_path) {
  n_snps <- length(chrom)

  # Ensure chromosome names have "chr" prefix
  if (!all(grepl("^chr", chrom)))
    chrom <- paste0("chr", chrom)

  # Build GRanges for row ranges
  gr <- GRanges(
    chrom,
    IRanges(
      start = as.integer(pos),
      end = as.integer(pos) + pmax(nchar(ref), nchar(alt)) - 1L,
      names = snp_ids))

  # Build VCF header
  coldata <- DataFrame(Samples = sample_name, row.names = sample_name)

  hdr <- VariantAnnotation::VCFHeader(
    header = DataFrameList(
      fileformat = DataFrame(
        Value = "VCFv4.2", row.names = "fileformat")),
    sample = sample_name)

  # Subset geno header to only fields present in geno
  geno_header <- geno_header[rownames(geno_header) %in% names(geno), , drop = FALSE]
  VariantAnnotation::geno(hdr) <- geno_header

  # Build VCF object
  geno_sl <- SimpleList(geno)
  vcf <- VariantAnnotation::VCF(
    rowRanges = gr,
    colData = coldata,
    exptData = list(header = hdr),
    geno = geno_sl)

  VariantAnnotation::ref(vcf) <- DNAStringSet(ref)
  VariantAnnotation::alt(vcf) <- DNAStringSetList(as.list(alt))
  VariantAnnotation::fixed(vcf)$FILTER <- "PASS"
  vcf <- sort(vcf)

  # Write based on output format
  # Note: VariantAnnotation::writeVcf appends ".bgz" to the path when
  # index = TRUE, so we must pass the path *without* the .bgz/.gz suffix.
  ext <- file_ext(output_path)
  if (ext == "bcf") {
    # Write temporary bgzipped VCF, then convert to BCF
    tmp_vcf_stem <- tempfile(fileext = ".vcf")
    tmp_vcf_bgz <- paste0(tmp_vcf_stem, ".bgz")
    on.exit(unlink(c(tmp_vcf_bgz, paste0(tmp_vcf_bgz, ".tbi")),
                   force = TRUE), add = TRUE)
    VariantAnnotation::writeVcf(vcf, tmp_vcf_stem, index = TRUE)
    # asBcf appends ".bcf" to destination, so strip the extension
    bcf_stem <- sub("\\.bcf$", "", output_path)
    dict <- unique(chrom)
    asBcf(tmp_vcf_bgz, dictionary = dict,
                     destination = bcf_stem)
  } else if (ext == "gz" || ext == "bgz") {
    # writeVcf will append .bgz, so strip it from the path
    vcf_stem <- sub("\\.(bgz|gz)$", "", output_path)
    VariantAnnotation::writeVcf(vcf, vcf_stem, index = TRUE)
    # writeVcf always creates .bgz; rename if the user requested .gz
    actual_path <- paste0(vcf_stem, ".bgz")
    if (actual_path != output_path && file.exists(actual_path)) {
      file.rename(actual_path, output_path)
      tbi_actual <- paste0(actual_path, ".tbi")
      if (file.exists(tbi_actual))
        file.rename(tbi_actual, paste0(output_path, ".tbi"))
    }
  } else {
    VariantAnnotation::writeVcf(vcf, output_path)
  }

  invisible(output_path)
}
