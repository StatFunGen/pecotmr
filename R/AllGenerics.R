#' @title S4 Generic Function Definitions
#' @description All S4 generic function definitions for pecotmr.
#' @name pecotmr-generics
#' @keywords internal
#' @importFrom methods setGeneric
NULL

#' Show methods for pecotmr S4 classes
#'
#' Compact console display methods for the package's S4 objects.
#'
#' @param object The object to display.
#' @return \code{object}, invisibly.
#' @examples
#' data(qtlDatasetExample)
#' show(qtlDatasetExample)
#' @name show-methods
#' @rdname show-methods
NULL

# =============================================================================
# High-level estimation generic
# =============================================================================

#' @title Estimate SNP Heritability
#' @description Estimate SNP heritability from GWAS summary statistics using
#'   one of four methods: S-LDSC, g-LDSC, LDER, or HDL/sHDL. Supplying
#'   \code{annotations} makes any of them stratified.
#' @param sumstats A \code{GwasSumStats} object.
#' @param ldRef An \code{LdStatistic} object (method-appropriate subclass):
#'   an \code{\link{LdScore}} for \code{"sldsc"} / \code{"gldsc"}, an
#'   \code{\link{LdEigen}} for \code{"lder"} / \code{"hdl"}. Build either
#'   from loaded LD with \code{\link{buildLdScore}} /
#'   \code{\link{buildLdEigen}}.
#' @param method Character, one of \code{"lder"} (the default),
#'   \code{"sldsc"}, \code{"gldsc"} or \code{"hdl"}:
#'   \describe{
#'     \item{\code{"sldsc"}}{Stratified LD score regression (Finucane et al.
#'       2015): weighted least squares of \eqn{\chi^2} on LD scores. Needs
#'       only per-variant scores when unstratified, so it is the cheapest of
#'       the four and the only one that runs on a scores-only
#'       \code{LdScore}.}
#'     \item{\code{"gldsc"}}{Generalized LD score regression (Xiong et al.
#'       2024): GLS using the full per-block LD matrices.}
#'     \item{\code{"lder"}}{LD eigenvalue regression (Song et al. 2022).}
#'     \item{\code{"hdl"}}{High-definition likelihood (Ning et al. 2020),
#'       stratified as sHDL. Models reference-panel noise explicitly, so it
#'       shrinks hard when \code{nRef} is small.}
#'   }
#' @param annotations An \code{AnnotationMatrix} object, or NULL for
#'   unstratified estimation.
#' @param local Logical, whether to compute per-block local estimates.
#' @param estimatorArgs Optional named list of estimator-specific options
#'   (\code{lambda} for lder / gldsc / hdl, \code{nIter} for sldsc).
#' @param ... Additional method-specific arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @return An \code{H2Estimate} object.
#' @examples
#' data(ldScoreExample)
#' # The estimators index the z-scores by the reference's own variant order,
#' # so build the sumstats from the reference rather than alongside it.
#' gr <- as(ldScoreExample, "GRanges")
#' set.seed(1)
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'   SNP = names(ldScoreExample),
#'   A1 = S4Vectors::mcols(ldScoreExample)$A1,
#'   A2 = S4Vectors::mcols(ldScoreExample)$A2,
#'   Z = rnorm(length(ldScoreExample)),
#'   N = 100000L)
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr")
#' )
#' ss <- GwasSumStats(studyName = "trait1", entry = list(gr),
#'   genome = unique(unname(GenomeInfoDb::genome(ldScoreExample))),
#'   ldReference = panel)
#' estimateH2(ss, ldScoreExample, method = "sldsc")
#' @export
setGeneric(
    "estimateH2",
    function(
        sumstats,
        ldRef,
        method = "lder",
        annotations = NULL,
        local = FALSE,
        ...
    ) {
        standardGeneric("estimateH2")
    }
)

# =============================================================================
# LD score computation
# =============================================================================

#' @title Compute LD Scores
#' @description Compute LD scores from an LD reference, optionally stratified by
#'   annotations. Each r-squared is ldsc's small-sample unbiased estimate,
#'   \eqn{r^2 - (1 - r^2)/(n - 2)} for a reference of \eqn{n} samples, as
#'   \code{\link{buildLdScore}} stores, so the two agree.
#' @param ldRef An \code{LdStatistic} object.
#' @param annotations An \code{AnnotationMatrix} object, or NULL.
#' @param ... Additional arguments.
#' @return A numeric matrix of LD scores (SNPs x annotations+1).
#' @examples
#' data(ldEigenExample)
#' computeLdScores(ldEigenExample)
#' @export
setGeneric("computeLdScores", function(ldRef, annotations = NULL, ...) {
    standardGeneric("computeLdScores")
})

# =============================================================================
# I/O generics
# =============================================================================

#' @title Read a Genotype Panel
#' @description Read genotype data (VCF, plink1, plink2, GDS or BGEN) as a
#'   \code{RangedSummarizedExperiment}: variants in \code{rowRanges},
#'   samples in \code{colData}, and dosages in a \code{DelayedArray} assay
#'   that reads from the file only when something touches it. Nothing is read
#'   here beyond the variant and sample metadata.
#'
#'   The panel is the same shape a \code{QtlDataset} uses for its genotypes
#'   and the same shape a collection carries as its \code{ldReference}, so the
#'   Bioconductor and tidy surfaces apply throughout: subset it with
#'   \code{panel[i, j]}, ask \code{rowRanges()} for the variants,
#'   \code{colnames()} for the samples, and hand it to
#'   \code{\link{computeLd}}.
#'
#'   A panel that is not one self-describing file is named by keyword
#'   instead of by \code{path}: pass \code{plink1Prefix} or
#'   \code{plink2Prefix} for a triplet stem, \code{bed}/\code{bim}/\code{fam}
#'   or \code{pgen}/\code{pvar}/\code{psam} for explicit triplet paths,
#'   \code{genoMeta} for a one-file-per-chromosome panel (optionally narrowed
#'   with \code{chroms}), or \code{ldMeta} plus \code{block} to resolve a
#'   block out of an LD-meta table. Exactly one source may be given.
#' @param path Character, path to a single self-describing genotype file
#'   (\code{.vcf.gz}, \code{.vcf.bgz}, \code{.gds} or \code{.bgen}). Omit it
#'   when naming the panel by one of the keyword sources described above.
#'
#'   A VCF must be bgzip-compressed with its tabix index (\code{.tbi}) beside
#'   it, through which its records are read a block at a time; a BCF cannot
#'   be read. Its dosage is its FORMAT/DS when the header declares one, else
#'   the copies of the ALT allele in each GT call, a call with a missing
#'   allele being missing.
#'
#'   A GDS file is SNPRelate's, whose genotype codes count the first of each
#'   SNP's two alleles (A1), or SeqArray's, read like a VCF: its FORMAT/DS
#'   when it stores one, else its genotype calls.
#'
#'   A BGEN file of any layout (BGEN 1.1, 1.2 or 1.3), phased or unphased
#'   and of any ploidy, needs its \code{.bgi} index beside it and its sample
#'   ids in the file or in a \code{.sample} file beside it. Its dosage is the
#'   expected count of A1, from 0 to the sample's ploidy.
#'
#'   A PLINK1 fileset, a \code{.bed} with its \code{.bim} and \code{.fam}
#'   beside it, is read through the PlinkMatrix package: its dosage counts
#'   each variant's A1 allele, the \code{.bim}'s fifth column, a missing
#'   call being missing. Its samples are named by the \code{.fam}'s
#'   individual ids.
#'
#'   In a VCF, a SeqArray GDS or a BGEN file, A2 is each record's first (REF)
#'   allele, and a record with more than one other allele is read, with a
#'   warning, as one variant per other allele, each counting its own.
#' @param format Character, one of "vcf", "plink1", "plink2", "gds", "bgen".
#'   If NULL, inferred from file extension.
#' @param ... The keyword source arguments described above, plus any further
#'   arguments forwarded to the format-specific reader.
#' @param plink1Prefix,plink2Prefix Character. Stem of a plink1
#'   (\code{.bed}/\code{.bim}/\code{.fam}) or plink2
#'   (\code{.pgen}/\code{.pvar}/\code{.psam}) triplet, without the extension.
#' @param bed,bim,fam Character. Explicit paths to the three files of a plink1
#'   triplet, when they do not share a stem.
#' @param pgen,pvar,psam Character. Explicit paths to the three files of a
#'   plink2 triplet, when they do not share a stem.
#' @param genoMeta One-file-per-chromosome specification: a path to a
#'   \code{#chr,path} meta file, or a named character vector (names =
#'   chromosomes, values = payload paths / prefixes). Every file must be of
#'   one format and hold the same samples in the same order.
#' @param chroms Optional character vector of chromosomes. Only the
#'   \code{genoMeta} shards whose chromosome is listed are read; a label may
#'   be in any style (\code{"chr1"} or \code{"1"}). Meaningful only with
#'   \code{genoMeta}.
#' @param ldMeta Path to an LD-meta TSV (columns \code{chrom}, \code{start},
#'   \code{end}, \code{path}). Requires \code{block}: the row covering it is
#'   resolved and the file it names is read with the matching backend.
#' @param block The \code{GRanges} window whose LD block to read out of
#'   \code{ldMeta}. It chooses a file rather than filtering records: the
#'   panel is the whole LD block that covers it.
#' @param vcfArgs Optional named list of arguments forwarded to
#'   \code{VariantAnnotation::readVcf} when a VCF's records are read (its
#'   genotypes are read a block at a time, as they are needed).
#' @return A \code{RangedSummarizedExperiment} of variants x samples.
#' @seealso \code{\link{computeLd}}
#' @examples
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr")
#' )
#' dim(panel)
#' head(rownames(panel))
#'
#' # The same panel named by its PLINK1 stem rather than by a file.
#' stem <- sub("\\.bed$", "", system.file(
#'   "extdata", "toy_ref.bed",
#'   package = "pecotmr"
#' ))
#' dim(readGenotypes(plink1Prefix = stem))
#' @export
setGeneric("readGenotypes", function(path, format = NULL, ...) {
    standardGeneric("readGenotypes")
})

#' @title Read Annotations
#' @description Read genomic annotations from files (BED, BigWig, S-LDSC .annot
#'   format, or GRanges objects) and create an AnnotationMatrix.
#' @param path Named character vector of file paths, or a named list of GRanges
#'   objects. Names become annotation names.
#' @param snpRanges A \code{GRanges} object defining SNP positions.
#' @param annotationMeta A \code{data.frame} with annotation metadata (name,
#'   tier, type). If NULL, auto-detected from file format.
#' @param genome Character, genome build.
#' @param ... Additional arguments.
#' @return An \code{AnnotationMatrix} object.
#' @examples
#' bedFile <- tempfile(fileext = ".bed")
#' writeLines("chr1\t100\t500\tregion1", bedFile)
#' snpRanges <- GenomicRanges::GRanges(
#'   rep("chr1", 3), IRanges::IRanges(c(50, 200, 600), width = 1))
#' readAnnotations(c(enhancer = bedFile), snpRanges, genome = "hg38")
#' @export
setGeneric(
    "readAnnotations",
    function(path, snpRanges, annotationMeta = NULL, genome = "hg19", ...) {
        standardGeneric("readAnnotations")
    }
)

# =============================================================================
# Accessor generics
# =============================================================================

#' @title Get Local Estimates
#' @description Extract per-block local estimates from a result object.
#' @param object An \code{H2Estimate} object.
#' @return A \code{data.frame} of local estimates, or NULL.
#' @examples
#' data(h2EstimateExample)
#' localH2Blocks(h2EstimateExample)
#' @export
setGeneric("localH2Blocks", function(object) standardGeneric("localH2Blocks"))

#' @title Get Enrichment Estimates
#' @description Extract annotation enrichment estimates from a result object.
#' @param object An \code{H2Estimate} object.
#' @return A \code{data.frame} of enrichment estimates, or NULL.
#' @examples
#' data(h2EstimateExample)
#' stratifiedHeritabilityEnrichment(h2EstimateExample)
#' @export
setGeneric("stratifiedHeritabilityEnrichment", function(object) {
    standardGeneric("stratifiedHeritabilityEnrichment")
})

#' @title Get Score Statistics
#' @description Extract score statistics for candidate annotations.
#' @param object An \code{H2Estimate} object.
#' @return A list with \code{z} and \code{R}, or NULL.
#' @examples
#' data(h2EstimateExample)
#' scoreStats(h2EstimateExample)
#' @export
setGeneric("scoreStats", function(object) standardGeneric("scoreStats"))

# =============================================================================
# GwasSumStats accessor generics
# =============================================================================

#' @title Get Z-scores
#' @description Extract z-score vector from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments (e.g., \code{study} for
#'   \code{GwasSumStats}; \code{study}, \code{context}, \code{trait} for
#'   \code{QtlSumStats}).
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of z-scores.
#' @export
setGeneric("z", function(x, ...) standardGeneric("z"))

#' @title Get the Number of Reference Samples
#' @description The number of samples an LD reference was built from: the
#'   samples carried by a \code{GenotypeHandle}, or the panel size an
#'   \code{\link{LdMatrix}} records. The per-variant sample size of summary
#'   statistics is \code{\link{sampleSize}}.
#' @param x An \code{LdMatrix} or \code{GenotypeHandle} object.
#' @param ... Class-specific selection arguments.
#' @return An integer of length 1 (\code{NA} when an \code{LdMatrix} records
#'   no size).
#' @details The genotype handle is the seed layer behind a panel and is not
#'   part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' nSamples(LdMatrix(cor(X), nSamples = nrow(X)))
#' @export
setGeneric("nSamples", function(x, ...) standardGeneric("nSamples"))

#' @title Get Per-Variant Sample Sizes
#' @description Extract the per-variant sample-size vector (the \code{N}
#'   column) of a \code{GwasSumStats} or \code{QtlSumStats} entry, selected
#'   by its identity tuple. Part of the summary-statistic column set alongside
#'   \code{\link{z}}, \code{\link{se}} and \code{\link{pval}}; replace it
#'   with \code{\link[=variantColumnReplace]{sampleSize<-}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of per-variant sample sizes, or \code{NULL} when the
#'   entry carries none.
#' @export
setGeneric("sampleSize", function(x, ...) standardGeneric("sampleSize"))

#' @title Get Association P-values
#' @description Extract the association p-value vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple. Part of the first-class summary-statistic column set alongside
#'   \code{\link{z}} / \code{\link{coef}} / \code{\link{se}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of p-values, or \code{NULL} if not available.
#' @export
setGeneric("pval", function(x, ...) standardGeneric("pval"))

#' @title Get Marginal Effect Sizes
#' @description Extract the marginal effect-size (beta) vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple. A GWAS / QTL effect size is a regression coefficient, so this
#'   adopts \code{stats::coef} rather than defining a pecotmr generic.
#'   Promoting the existing function keeps it as the default method, so
#'   \code{coef()} on a model fit is unaffected.
#' @param object A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of effect sizes, or \code{NULL} if not available.
#' @importFrom stats coef
#' @rdname coef-methods
setGeneric("coef")

#' @title Get Effect-Size Standard Errors
#' @description Extract the effect-size standard-error vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of standard errors, or \code{NULL} if not available.
#' @export
setGeneric("se", function(x, ...) standardGeneric("se"))

#' @title Get Minor Allele Frequencies
#' @description Extract the minor-allele-frequency vector from one entry of
#'   a summary-statistic collection, or from a \code{QtlDataset} over the
#'   whole dataset (narrow it first, with \code{x[i, j, ]}, to read part of
#'   it).
#' @param x A \code{GwasSumStats}, \code{QtlSumStats} or \code{QtlDataset}
#'   object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric vector of MAFs, or NULL if not available.
#' @export
setGeneric("maf", function(x, ...) standardGeneric("maf"))

#' @title Replace a Per-Variant Column
#' @description Supply or overwrite one per-variant column of a
#'   \code{GwasSumStats} / \code{QtlSumStats} entry after construction --
#'   the usual case being an effect-allele frequency the source file did not
#'   carry. The value is positional: a vector in \code{\link{sumStats}}
#'   order for the selection -- a GWAS study split across rows is stitched,
#'   as \code{sumStats()} returns it -- with exactly one value per variant.
#'   Length is checked hard, because a recycled or truncated vector is
#'   otherwise silently wrong. To take the values from a reference keyed by
#'   variant instead, matched on \code{(chrom, pos, A1, A2)} and oriented to
#'   the entry's alleles, use \code{\link{fillFromReference}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param studyName,context,trait Selectors pinning one study (every row of a
#'   GWAS study) or one \code{QtlSumStats} entry; \code{context} and
#'   \code{trait} are \code{QtlSumStats}-only.
#' @param ... Unused.
#' @param value A \code{numeric} vector, one value per variant in
#'   \code{sumStats()} order.
#' @return The object, with the selection's column replaced. \code{sampleSize<-}
#'   writes the per-variant \code{N}.
#' @examples
#' data(qtlSumStatsExample)
#' x <- qtlSumStatsExample
#' af(x) <- runif(length(variantIds(x)))
#' head(af(x))
#' @rdname variantColumnReplace
#' @aliases variantColumnReplace
#' @export
setGeneric("af<-", function(x, ..., value) standardGeneric("af<-"))

#' @rdname variantColumnReplace
#' @export
setGeneric("maf<-", function(x, ..., value) standardGeneric("maf<-"))

#' @rdname variantColumnReplace
#' @export
setGeneric("sampleSize<-", function(x, ..., value) {
    standardGeneric("sampleSize<-")
})

#' @title Get Effect-Allele Frequencies
#' @description Extract the directional effect-allele (A1) frequency vector
#'   for a \code{QtlDataset}, an \code{\link{LdMatrix}}, or one entry of a
#'   \code{GwasSumStats} / \code{QtlSumStats} collection. Unlike
#'   \code{\link{maf}}, the value is \emph{not} folded to the minor allele: it
#'   is the frequency of the dosage-counted allele (A1, the effect allele),
#'   matching the allele the marginal effect sizes and the fine-mapping
#'   \code{af} column report. For a summary-statistic collection this is
#'   \code{NULL} when the source
#'   carried no frequency column -- see \code{\link{variantColumnReplace}}
#'   for supplying one afterwards. A \code{QtlDataset} answers over the whole
#'   dataset; narrow it first, with \code{x[i, j, ]}, to read part of it.
#' @param x A \code{QtlDataset}, \code{LdMatrix}, \code{GwasSumStats} or
#'   \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Named numeric vector of effect-allele frequencies (names are variant
#'   IDs), or an empty vector when no variants are selected.
#' @examples
#' data(qtlDatasetExample)
#' af(qtlDatasetExample)
#' @export
setGeneric("af", function(x, ...) standardGeneric("af"))

#' @title Get Number of SNPs
#' @description Number of SNPs in a \code{GwasSumStats} or \code{QtlSumStats}
#'   entry, selected by its identity tuple, or the number of variants an
#'   \code{H2Estimate} was computed over.
#' @param x A \code{GwasSumStats}, \code{QtlSumStats} or \code{H2Estimate}
#'   object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Integer.
#' @export
setGeneric("nSnps", function(x, ...) standardGeneric("nSnps"))

#' @title Hierarchical Multiple-Testing Correction for cis-QTL Association
#' @description Apply hierarchical (local per-gene + global across-gene)
#'   multiple-testing correction to a per-gene \code{QtlSumStats} of cis-QTL
#'   association statistics (one row per gene/trait), returning the same object
#'   enriched with the corrected-statistic columns. See
#'   \code{\link{qtlAssociationPostprocess}}.
#' @param x A \code{QtlSumStats}.
#' @param ... Correction arguments.
#' @return The input \code{QtlSumStats} with added correction columns.
#' @examples
#' pBeta <- c(10^(-c(8, 7, 6, 5, 4, 3)), stats::ppoints(54))
#' G <- length(pBeta)
#' entries <- lapply(seq_len(G), function(i) {
#'   gr <- GenomicRanges::GRanges("chr1",
#'     IRanges::IRanges(seq(1000L, by = 50L, length.out = 4L), width = 1L))
#'   S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'     SNP = paste0("g", i, "_v", 1:4), A1 = "A", A2 = "G",
#'     P = c(pBeta[i] / 5, 0.2, 0.5, 0.8),
#'     af = c(0.3, 0.2, 0.005, 0.4),
#'     tss_distance = c(0L, 500L, 900000L, 2000000L),
#'     tes_distance = c(0L, 500L, 900000L, 2000000L))
#'   gr
#' })
#' qss <- QtlSumStats(studyName = rep("s", G), context = rep("brain", G),
#'   trait = paste0("g", seq_len(G)), entry = entries, genome = "hg19",
#'   extraCols = list(
#'     n_variants = rep(50L, G), n_variants_filtered = rep(30L, G),
#'     p_beta = pBeta, beta_shape1 = rep(1, G), beta_shape2 = rep(200, G)))
#' qtlAssociationPostprocess(qss)
#' @export
setGeneric("qtlAssociationPostprocess", function(x, ...) {
    standardGeneric("qtlAssociationPostprocess")
})

#' @title Extract Significant cis-QTL Variants
#' @description Derive the significant variants under a given correction method
#'   and FDR threshold from a \code{QtlSumStats} enriched by
#'   \code{\link{qtlAssociationPostprocess}} (significance is computed on
#'   demand, never stored).
#' @param x A \code{QtlSumStats}.
#' @param ... Selection arguments (\code{method}, \code{threshold}).
#' @return A \code{GRanges} (or list of GRanges) of the significant variants.
#' @examples
#' pBeta <- c(10^(-c(8, 7, 6, 5, 4, 3)), stats::ppoints(54))
#' G <- length(pBeta)
#' entries <- lapply(seq_len(G), function(i) {
#'   gr <- GenomicRanges::GRanges("chr1",
#'     IRanges::IRanges(seq(1000L, by = 50L, length.out = 4L), width = 1L))
#'   S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
#'     SNP = paste0("g", i, "_v", 1:4), A1 = "A", A2 = "G",
#'     P = c(pBeta[i] / 5, 0.2, 0.5, 0.8),
#'     af = c(0.3, 0.2, 0.005, 0.4),
#'     tss_distance = c(0L, 500L, 900000L, 2000000L),
#'     tes_distance = c(0L, 500L, 900000L, 2000000L))
#'   gr
#' })
#' qss <- QtlSumStats(studyName = rep("s", G), context = rep("brain", G),
#'   trait = paste0("g", seq_len(G)), entry = entries, genome = "hg19",
#'   extraCols = list(
#'     n_variants = rep(50L, G), n_variants_filtered = rep(30L, G),
#'     p_beta = pBeta, beta_shape1 = rep(1, G), beta_shape2 = rep(200, G)))
#' pp <- qtlAssociationPostprocess(qss)
#' significantQtls(pp)
#' @export
setGeneric("significantQtls", function(x, ...) {
    standardGeneric("significantQtls")
})

#' @title Reconcile Two Collections to Their Shared Variants
#' @description Restrict two fine-mapping collections to the variants they have
#'   in common, renormalizing each retained single-effect posterior over the
#'   variants that survive.
#' @details This intersects two objects against each other. Coloc needs it
#'   symmetric --
#'   both sides must be scored on the same variant set -- while TWAS, MR and
#'   cTWAS need it one-sided, adjusting the QTL side to the GWAS variant set
#'   and leaving the GWAS side alone.
#'
#'   Reconciliation renormalizes the STORED alpha rather than rebuilding it
#'   from \code{lbf_variable}: alpha already carries the fit's prior, so
#'   restricting and renormalizing is exact for any prior weights, whereas the
#'   rebuild silently substitutes a uniform one.
#'
#'   Renormalizing matches a refit rather than approximating one -- a refit on
#'   the retained set normalizes the same single-effect regression over the
#'   same variants. What it cannot recover is the information the dropped
#'   variants carried: coverage falls as overlap shrinks, equally for both, and
#'   that degradation is irreducible rather than an artifact of the method.
#'   \code{retainedMass()} is what makes it visible.
#' @param x A fine-mapping collection to reconcile.
#' @param y The collection to reconcile it against.
#' @param oneSided Logical (length 1). \code{FALSE} (default) returns both
#'   sides restricted to the shared variants. \code{TRUE} returns only
#'   \code{x}, adjusted to \code{y}'s variants, leaving \code{y} untouched.
#' @param ... Additional arguments passed to methods.
#' @return With \code{oneSided = FALSE}, a list of the two reconciled
#'   collections named \code{x} and \code{y}. With \code{oneSided = TRUE},
#'   the reconciled \code{x}.
#' @examples
#' data(qtlFineMappingExample)
#' data(gwasFineMappingExample)
#' both <- intersectVariants(qtlFineMappingExample, gwasFineMappingExample)
#' names(both)
#' @export
setGeneric("intersectVariants", function(x, y, oneSided = FALSE, ...) {
    standardGeneric("intersectVariants")
})

#' @title Per-Effect Retained Posterior Mass
#' @description The share of each single effect's posterior that survived a
#'   variant subset -- the diagnostic that says how much a reconciliation cost.
#' @details Reported always, never used as a filter. The mass is
#'   \strong{continuous}, not bimodal, once the causal variant is untyped
#'   (median 0.78 at 30% missingness, 0.58 at 50%, with most loci between 0.1
#'   and 0.9), so there is no natural cut point and no default drop threshold
#'   is imposed.
#'
#'   susieR's own \code{V -> 0} null rule is not a substitute: across 135
#'   effects it fired zero times, including all 52 whose retained mass was
#'   ~0 -- retained variants in LD with a dropped causal variant still beat the
#'   null. An explicit diagnostic is the only thing that surfaces this.
#' @param x A fine-mapping collection.
#' @param ... Additional arguments passed to methods.
#' @return A \code{DataFrame}, one row per (element, effect), with the
#'   element's identity columns, the effect index, and its retained mass.
#' @examples
#' data(qtlFineMappingExample)
#' data(gwasFineMappingExample)
#' # Nothing has been reconciled yet, so there is no mass to report.
#' nrow(retainedMass(qtlFineMappingExample))
#' both <- intersectVariants(qtlFineMappingExample, gwasFineMappingExample)
#' retainedMass(both$x)
#' @export
setGeneric("retainedMass", function(x, ...) {
    standardGeneric("retainedMass")
})

#' @title Get Phenotype Variance
#' @description Extract phenotype variance from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple. Returns
#'   \code{NULL} when the entry has no \code{varY} recorded.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric phenotype variance, or NULL.
#' @examples
#' data(qtlSumStatsExample)
#' varY(qtlSumStatsExample)
#' @export
setGeneric("varY", function(x, ...) standardGeneric("varY"))

#' @title Get a Single Summary-Statistic Entry or Embedded Collection
#' @description Behavior depends on the class of \code{x}:
#'   \describe{
#'     \item{For \code{GwasSumStats} / \code{QtlSumStats}}{Returns the
#'       per-variant \code{GRanges} of summary statistics for one entry,
#'       selected by its identity tuple (\code{study} for GWAS;
#'       \code{study}, \code{context}, \code{trait} for QTL).}
#'     \item{For \code{MultiStudyQtlDataset}}{Returns the embedded
#'       \code{QtlSumStats} collection (the summary-statistic-only
#'       studies), or \code{NULL} when absent. No selection arguments
#'       are accepted in this case.}
#'   }
#' @param x A \code{GwasSumStats}, \code{QtlSumStats}, or
#'   \code{MultiStudyQtlDataset} object.
#' @param ... Class-specific selection arguments (see above).
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return A \code{GRanges}, a \code{QtlSumStats}, or \code{NULL}.
#' @examples
#' data(qtlSumStatsExample)
#' sumStats(qtlSumStatsExample)
#' @export
setGeneric("sumStats", function(x, ...) standardGeneric("sumStats"))

#' @title Get the Per-Row QC Audit
#' @description What \code{\link{summaryStatsQc}} did to each row of a
#'   summary-statistics collection: one row per tuple, its identity columns,
#'   then the variant counts each QC step recorded (\code{variantsIn},
#'   \code{variantsOut}, and so on; a nested count is named by its path, as
#'   in \code{contentFilters.nDropped}). A step that did not run on a row is
#'   \code{NA} there. Only the identity columns are present before QC has run.
#'   The options the QC ran with are \code{\link{qcOptions}}; the per-variant
#'   LD-mismatch tables are \code{\link{qcDiagnostics}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused.
#' @return A \code{DataFrame}, one row per tuple.
#' @export
setGeneric("qcInfo", function(x, ...) standardGeneric("qcInfo"))

#' @title Get the QC Options
#' @description The record of the QC that ran over a summary-statistics
#'   collection: the options \code{\link{summaryStatsQc}} was called with and
#'   when. An empty \code{list()} (the default on construction) means QC has
#'   not run; pipelines that need harmonized summary statistics
#'   (\code{fineMappingPipeline}, \code{twasWeightsPipeline} and others)
#'   refuse such inputs. What QC did to each row is \code{\link{qcInfo}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused.
#' @return A \code{list}, possibly empty.
#' @export
setGeneric("qcOptions", function(x, ...) standardGeneric("qcOptions"))

#' @title Get SLALOM / DENTIST Diagnostics
#' @description Return the per-variant LD-mismatch diagnostics (SLALOM or
#'   DENTIST output) that \code{\link{summaryStatsQc}} attached to the
#'   selected entry when \code{ldMismatchQcMethod} was not \code{"none"}.
#'   The selectors pin one study (every row of a GWAS study, whose rows'
#'   tables are stacked) or one \code{QtlSumStats} entry, as for
#'   \code{\link{sumStats}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param studyName,context,trait Selectors pinning one study or entry;
#'   \code{context} and \code{trait} are \code{QtlSumStats}-only.
#' @param ... Unused.
#' @return A \code{DataFrame} of per-variant diagnostics (columns include
#'   \code{variant_id} plus the SLALOM / DENTIST output columns), or
#'   \code{NULL} when no diagnostics were preserved for the selection.
#' @export
setGeneric("qcDiagnostics", function(x, ...) {
    standardGeneric("qcDiagnostics")
})

#' @title Get LD Reference
#' @description Return the LD references of this object, each a genotype
#'   panel (see \code{\link{readGenotypes}}), a \code{\link{GenotypeMixture}}
#'   or an \code{\link{LdMatrix}}.
#'
#'   The collections whose rows can sit on different references --
#'   \code{GwasSumStats}, \code{QtlSumStats}, \code{FineMappingResultBase},
#'   \code{TwasWeights}, \code{ColocResult} and \code{ColocBoostResult} --
#'   return a named list holding each distinct reference their rows use, once.
#'   The \code{ldReference} column of their \code{mcols} (of
#'   \code{outcomeInfo} for a \code{ColocBoostResult}) names the entry each
#'   row used, \code{NA} for one fit without LD (e.g. a \code{TwasWeights}
#'   row learned from individual-level data via \code{QtlDataset}). The list
#'   is empty when no row has a reference.
#'
#'   An \code{LdStatistic} returns the one panel it reads its blocks' LD from,
#'   or \code{NULL} when it was built from \code{LdMatrix} blocks.
#' @param x An S4 object that carries an \code{ldReference} slot.
#' @param ... Unused.
#' @return A named list of LD references; for an \code{LdStatistic}, a
#'   \code{RangedSummarizedExperiment} genotype panel, a
#'   \code{GenotypeMixture}, an \code{LdMatrix}, or \code{NULL}.
#' @export
setGeneric("ldReference", function(x, ...) standardGeneric("ldReference"))

# =============================================================================
# LD accessor generics
# =============================================================================

#' @title Get LD Correlation Matrix
#' @description Extract the LD correlation matrix an \code{\link{LdMatrix}}
#'   holds, labelled by its variant ids. \code{as.matrix(x)} and
#'   \code{as(x, "matrix")} return the same matrix.
#' @param x An \code{LdMatrix}.
#' @param ... Unused.
#' @return A correlation matrix.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' ld <- LdMatrix(cor(X), nSamples = nrow(X))
#' ldMatrix(ld)[1:3, 1:3]
#' identical(as.matrix(ld), ldMatrix(ld))
#' @export
setGeneric("ldMatrix", function(x) standardGeneric("ldMatrix"))

#' @title Get Genotype Matrix
#' @description Extract a genotype matrix from an object that carries genotype
#'   data. For a \code{QtlDataset}, the whole dataset's dosages (samples x
#'   variants), with no QC applied: to read part of it, narrow the dataset
#'   first, e.g. \code{genotypes(x[window, samples, ])}.
#' @param x The object to extract from.
#' @param ... Passed on to methods.
#' @return A numeric matrix, a list of matrices, or \code{NULL}.
#' @examples
#' data(qtlDatasetExample)
#' genotypes(qtlDatasetExample)
#' @export
setGeneric("genotypes", function(x, ...) standardGeneric("genotypes"))

#' @title Get Variant IDs
#' @description Extract variant ID vector from an object that carries one (e.g.,
#'   \code{LdMatrix}, \code{FineMappingRow}, \code{TwasWeightsRow}) or from
#'   one entry of a collection class selected by its identity tuple. Every
#'   class renders the same canonical \code{chrom:pos:A2:A1} id from each
#'   variant's range and alleles; labels a source gave (an rsID, a file's own
#'   spelling) stay in the \code{SNP} column of its ranges.
#' @param x The object.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @return Character vector of variant IDs.
#' @examples
#' data(qtlFineMappingExample)
#' variantIds(qtlFineMappingExample)
#' @export
setGeneric("variantIds", function(x, ...) standardGeneric("variantIds"))

#' @title Get Molecular Traits
#' @description Extract molecular-trait data from an object that carries it.
#'   For a \code{QtlDataset}, select contexts and, optionally, traits by
#'   name. The data come back as stored: missing values and outliers are a
#'   run's phenotype QC (\code{\link{PhenotypeFilterParam}}), applied by the
#'   pipelines, not by this accessor.
#' @param x The object to extract from.
#' @param ... Class-specific selection arguments (e.g., \code{context},
#'   \code{trait}).
#' @param context Character vector. Context(s) whose data to extract.
#' @param trait Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @return A named list of \code{SummarizedExperiment}s keyed by context,
#'   one per selected context. A \code{QtlDataset} is a
#'   \code{MultiAssayExperiment}, so these are its experiments; take the
#'   measurement matrix itself with
#'   \code{SummarizedExperiment::assay()} rather than asking this accessor
#'   for one.
#' @examples
#' data(qtlDatasetExample)
#' molecularTraits(qtlDatasetExample, context = "brain")
#' @export
setGeneric("molecularTraits", function(x, ...) {
    standardGeneric("molecularTraits")
})
# =============================================================================
# Fine-mapping result accessor generics
# =============================================================================

#' @title Get One Fine-Mapping Fit
#' @description Return the rows of one fit of a \code{QtlFineMappingResult}
#'   or \code{GwasFineMappingResult}, selected by its identity, as a
#'   collection of the same class. A QTL fit is one row. A GWAS study
#'   fine-mapped block by block holds one row per LD block, so its fit is
#'   every block of the \code{(study, method)} pair, unless \code{blockId}
#'   names one or \code{region} narrows them.
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param studyName,context,trait,method Single character identifiers.
#'   Required when the collection holds more than one fit; optional when it
#'   holds one. \code{context} and \code{trait} are QTL-only.
#' @param blockId GWAS only: the \code{blockId} of one LD block, to select
#'   that block alone. \code{NULL} (default) selects every block of the fit.
#' @return A collection of the same class as \code{x}: one row for a QTL fit,
#'   one row per selected block for a GWAS fit.
#' @examples
#' data(qtlFineMappingExample)
#' fineMappingResult(qtlFineMappingExample, studyName = "study_1",
#'   context = "context_1", trait = "gene_1", method = "susie")
#' @export
setGeneric(
    "fineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        standardGeneric("fineMappingResult")
    }
)

#' @title Renormalize Fine-Mapping PIPs to a Variant Subset
#' @description Re-derive a \code{FineMappingRow}'s PIPs (and the
#'   \code{topLoci} table) after restricting to a kept variant subset. For each
#'   effect the \code{lbf_variable} row is subset to the kept variants,
#'   renormalized via \code{lbfToAlpha()}, and the per-variant PIPs are
#'   recomputed as \code{1 - prod_l(1 - alpha[l, p])}.
#'
#'   The two scenarios this supports:
#'   \itemize{
#'     \item The user declined to impute missing variants in a GWAS
#'           \code{SumStats}, so a downstream fine-mapping result needs
#'           PIPs restricted to the GWAS-covered intersection.
#'     \item Colocalization between a \code{GwasFineMappingResult} and a
#'           \code{QtlFineMappingResult} computed on different variant
#'           sets -- the GWAS PIPs (or QTL PIPs) get renormalized to the
#'           common variant set.
#'   }
#'
#' @param x A \code{FineMappingRow} or \code{FineMappingResultBase}.
#' @param keepVariants Character vector of variant IDs to keep. Intersected with
#'   the entry's own \code{variantIds}; an empty intersection raises an error.
#' @param ... Future expansion.
#' @return The same flavour of object with PIPs renormalized on the kept subset.
#' @noRd
setGeneric("adjustPips", function(x, keepVariants, ...) {
    standardGeneric("adjustPips")
})

#' @title Get PIP Values
#' @description Extract posterior inclusion probabilities from a
#'   \code{QtlFineMappingResult} or \code{GwasFineMappingResult}, narrowed to
#'   the entries its identity selectors pick.
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @return A named numeric vector of PIPs.
#' @examples
#' data(qtlFineMappingExample)
#' pip(qtlFineMappingExample)
#' @export
setGeneric("pip", function(x, ...) standardGeneric("pip"))

#' @title Get SuSiE Fit
#' @description Extract the SuSiE fit object from a fine-mapping entry or
#'   result. The fit may be the trimmed view (when the pipeline ran with the
#'   default \code{trim = TRUE}) or the full untrimmed \code{susie()} return
#'   (when \code{trim = FALSE}).
#' @param x A \code{FineMappingRow}, \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @return A list (the SuSiE fit object).
#' @examples
#' data(qtlFineMappingExample)
#' susieFit(qtlFineMappingExample)
#' @export
setGeneric("susieFit", function(x, ...) standardGeneric("susieFit"))

#' @title Compute Between-Credible-Set Correlation On Demand
#' @description Compute the between-credible-set correlation matrix for a
#'   \code{\link{FineMappingRow}} on demand from an LD source. The correlation
#'   is a view over the fit-time LD, which lives on the data object --- a
#'   \code{QtlDataset}'s genotypes or a \code{QtlSumStats} / \code{GwasSumStats}
#'   LD reference --- and is deliberately never stored on the fit. An
#'   \code{ldSource} is therefore always required.
#' @param x A \code{\link{FineMappingRow}} carrying the fit (credible-set
#'   membership and PIP) and its variant ids.
#' @param ldSource The object supplying the LD: a \code{QtlDataset}
#'   (individual-level; genotypes via \code{\link{genotypes}}) or a
#'   \code{QtlSumStats} / \code{GwasSumStats} (summary statistics; LD via
#'   \code{\link{ldReference}}). Any other type is an error.
#' @param ... Passed to methods.
#' @return A numeric \code{m x m} between-credible-set correlation matrix
#'   (\code{m} = number of credible sets), or \code{NULL} when the fit has
#'   fewer than two credible sets.
#' @examples
#' data(qtlFineMappingPairedExample)
#' data(qtlSumStatsExample)
#' fe <- fineMappingResult(qtlFineMappingPairedExample)
#' computeCsCorrelation(fe, qtlSumStatsExample)
#' @export
setGeneric("computeCsCorrelation", function(x, ldSource, ...) {
    # The accepted `ldSource` classes are known, so this is a validity check
    # rather than an ANY method: running it here means the message survives
    # whatever else the call passed, and no method needs `...` to catch it.
    .csCorrAssertLdSource(ldSource)
    standardGeneric("computeCsCorrelation")
})

# @noRd
.csCorrAssertLdSource <- function(ldSource) {
    ok <- !missing(ldSource) &&
        any(map_lgl(c("QtlDataset", "SumStatsBase"), is, object = ldSource))
    if (ok) {
        return(invisible(NULL))
    }
    abort(glue(
        "computeCsCorrelation() requires a QtlDataset, QtlSumStats, or ",
        "GwasSumStats as `ldSource`: the between-credible-set correlation ",
        "is derived from that object's LD and is never stored on the fit."
    ))
}

#' @title Get Cross-Validation Result
#' @description Extract the cross-validation payload stored on a
#'   \code{FineMappingRow} (or the matching entry of a
#'   \code{QtlFineMappingResult}). The payload is a list with components
#'   \code{samplePartition} (a \code{data.frame} of \code{Sample}/\code{Fold}
#'   assignments), \code{predictions} (a named list of per-method out-of-fold
#'   prediction matrices), and \code{performance} (a named list of per-method
#'   metric matrices). \code{NULL} when fine-mapping was run without
#'   cross-validation.
#' @param x A \code{FineMappingRow} or \code{QtlFineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A list (the CV payload) or \code{NULL}.
#' @examples
#' data(qtlFineMappingExample)
#' cvResult(qtlFineMappingExample)
#' @export
setGeneric("cvResult", function(x, ...) standardGeneric("cvResult"))

#' @title Get Marginal Effects
#' @description Extract per-variant marginal univariate effects from a
#'   fine-mapping entry or result, one range per variant. Its \code{mcols()}
#'   hold the identity columns (\code{variant_id, chrom, pos, A1, A2}),
#'   context (\code{N, MAF}), and the marginal effect columns
#'   (\code{beta, se, z, p}). Populated uniformly across the
#'   individual-level and RSS paths.
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @return A \code{GRanges}, one range per variant; when the selectors match
#'   more than one entry the row-identity columns are carried in
#'   \code{mcols()} too.
#' @examples
#' data(qtlFineMappingExample)
#' marginalEffects(qtlFineMappingExample)
#' @export
setGeneric("marginalEffects", function(x, ...) {
    standardGeneric("marginalEffects")
})

#' @title Get Top Loci (posterior view)
#' @description Extract the per-variant posterior fine-mapping payload as either
#'   a \code{data.frame} (default) or a \code{GRanges}. Returns identity columns
#'   (\code{variant_id, chrom, pos, A1, A2}), context (\code{N, MAF}), the
#'   posterior effect columns (\code{beta = posterior_mean, se = posterior_sd}),
#'   \code{pip}, and credible-set membership columns (\code{cs_95}, etc.). Rows
#'   are filtered by PIP by default -- set \code{signalCutoff = 0} to return
#'   every variant.
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param signalCutoff Numeric (length 1). Drop rows where \code{pip <=
#'   signalCutoff}. Default \code{0.025}. Use \code{signalCutoff = 0} to keep
#'   every variant.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @param minPurity Numeric or \code{NULL}. Minimum credible-set purity to
#'   retain; \code{NULL} applies no purity filter.
#' @return A \code{GRanges}, one range per retained variant, carrying the
#'   per-variant columns in \code{mcols()}; when the selectors match more than
#'   one entry the row-identity columns are carried there too. Take a table
#'   with \code{as.data.frame()} (which adds the standard
#'   \code{seqnames}/\code{start}/... columns) or
#'   \code{as.data.frame(S4Vectors::mcols(.))} for the per-variant columns
#'   alone.
#' @examples
#' data(qtlFineMappingExample)
#' topLoci(qtlFineMappingExample)
#' @export
setGeneric(
    "topLoci",
    function(x, signalCutoff = 0.025, ...) standardGeneric("topLoci")
)

#' @title Get Credible Sets
#' @description Extract credible set assignments at the requested coverage.
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Class-specific selection arguments plus \code{coverage}.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block of a fit; \code{NULL} matches every block.
#' @param coverage Numeric in (0, 1). Credible-set coverage level. Default
#'   \code{0.95}.
#' @param minPurity Numeric or \code{NULL}. Minimum credible-set purity to
#'   retain; \code{NULL} applies no purity filter.
#' @return A \code{GRanges}, one range per credible-set member, carrying the
#'   per-variant columns in \code{mcols()}; when the selectors match more
#'   than one entry the row-identity columns are carried there too. Narrow it
#'   with \code{subsetByOverlaps()}; take a table with
#'   \code{as.data.frame()}.
#' @examples
#' data(qtlFineMappingExample)
#' credibleSets(qtlFineMappingExample)
#' @export
setGeneric("credibleSets", function(x, ...) standardGeneric("credibleSets"))

# =============================================================================
# TwasWeights accessor generics
# =============================================================================

#' @title Resolve Per-Variant Weights From a Weight Source
#' @description Return an aligned \code{(variantIds, weights)} pair from a
#'   single weight-source entry, so cTWAS and
#'   \code{\link{causalInferencePipeline}} extract weights identically whether
#'   the source is a \code{TwasWeights} collection (its learned weight vector)
#'   or a fine-mapping one (its topLoci posterior effect). The selectors pin
#'   the single row the weights come from.
#' @param x A \code{TwasWeights}, \code{QtlFineMappingResult} or
#'   \code{GwasFineMappingResult}.
#' @param studyName,context,trait,method Optional length-1 selectors pinning one
#'   row of a collection, as elsewhere; each \code{NULL} (default) leaves that
#'   part of the tuple unconstrained. Ignored when \code{x} is already a row.
#' @param blockId GWAS fine-mapping results only: the \code{blockId} naming
#'   one LD block; \code{NULL} (default) leaves it unconstrained.
#' @param ... Reserved for future use.
#' @return A \code{list} with \code{variantIds} (character) and \code{weights}
#'   (numeric) of equal length; both empty when no usable weights are present.
#' @examples
#' data(twasWeightsExample)
#' w <- resolveWeights(twasWeightsExample, studyName = "protocol_example",
#'   context = "bulk_rnaseq", trait = "ENSG00000130538", method = "susie")
#' length(w$variantIds)
#' @export
setGeneric("resolveWeights", function(x, ...) standardGeneric("resolveWeights"))

#' @title Get cTWAS Fine-mapping Posteriors
#' @description Return the per-gene (and per-SNP) fine-mapping posterior table
#'   from a \code{\link{CtwasResultEntry}} or a \code{\link{CtwasResult}}
#'   collection (aggregated across rows, tagged with run identity).
#' @param x A \code{CtwasResultEntry} or \code{CtwasResult}.
#' @param ... Class-specific selection arguments.
#' @return On a \code{CtwasResultEntry}, the stored posteriors: a
#'   \code{GRanges} when the pipeline placed its rows (each table column in
#'   \code{mcols()}), the table as given otherwise, or \code{NULL} when
#'   absent. On a \code{CtwasResult}, one \code{GRanges} of every row's
#'   posteriors with the run identity in \code{mcols()} (empty when no row
#'   has any); a row whose payload carries no coordinates is an error.
#' @examples
#' cre <- CtwasResultEntry(
#'   posteriors = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' ctwasPosteriors(cre)
#' @export
setGeneric("ctwasPosteriors", function(x, ...) {
    standardGeneric("ctwasPosteriors")
})

#' @title Get cTWAS Per-effect Susie Alpha Table
#' @description Return the full per-effect susie alpha table
#'   (\code{ctwas::finemap_regions} \code{susie_alpha_res} shape) from a
#'   \code{\link{CtwasResultEntry}} or a \code{\link{CtwasResult}} collection
#'   (aggregated across rows, tagged with run identity). This is the fuller
#'   cTWAS output retained so the raw run is reconstructable.
#' @param x A \code{CtwasResultEntry} or \code{CtwasResult}.
#' @param ... Class-specific selection arguments.
#' @return On a \code{CtwasResultEntry}, the stored alphas: a
#'   \code{GRanges} when the pipeline placed its rows, the table as given
#'   otherwise, or \code{NULL} when absent. On a \code{CtwasResult}, one
#'   \code{GRanges} of every row's alphas with the run identity in
#'   \code{mcols()} (empty when no row has any).
#' @examples
#' cre <- CtwasResultEntry(
#'   posteriors = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' susieAlpha(cre)
#' @export
setGeneric("susieAlpha", function(x, ...) standardGeneric("susieAlpha"))

#' @title Get cTWAS Group Prior Parameters
#' @description Return the estimated \code{group_prior} / \code{group_prior_var}
#'   of a \code{\link{CtwasResultEntry}}.
#' @param x A \code{CtwasResultEntry}.
#' @param ... Reserved for future use.
#' @return The stored parameter object (typically a list), or \code{NULL}.
#' @examples
#' cre <- CtwasResultEntry(
#'   posteriors = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' ctwasGroupPriors(cre)
#' @export
setGeneric("ctwasGroupPriors", function(x, ...) {
    standardGeneric("ctwasGroupPriors")
})

#' @title Get Standardized Flag
#' @description Check whether weights are on the standardized scale.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return Logical.
#' @examples
#' data(ctwasWeightsExample)
#' weightStandardized(ctwasWeightsExample)
#' @export
setGeneric("weightStandardized", function(x, ...) {
    standardGeneric("weightStandardized")
})

#' @title Get Model Fits
#' @description Extract fitted model objects.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return Method-specific (typically a list).
#' @examples
#' data(ctwasWeightsExample)
#' # NULL unless the pipeline kept fits (fitRetention != "none")
#' methodFits(ctwasWeightsExample)
#' @export
setGeneric("methodFits", function(x, ...) standardGeneric("methodFits"))

#' @title Get the Full-Data Prior from a MashPrior
#' @description Accessor for the \code{fullFit} slot (the full-data data-driven
#'   prior payload).
#' @param x A \code{MashPrior} object.
#' @param ... Unused.
#' @return The full-data prior payload, or \code{NULL}.
#' @examples
#' U <- list(shared = diag(3), singleton = matrix(0.3, 3, 3) + diag(0.7, 3))
#' mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))
#' fullFit(mp)
#' @export
setGeneric("fullFit", function(x, ...) standardGeneric("fullFit"))

#' @title Get the Per-Fold Priors from a MashPrior
#' @description Accessor for the \code{cvFits} slot (per-fold priors +
#'   \code{samplePartition}).
#' @param x A \code{MashPrior} object.
#' @param ... Unused.
#' @return The \code{cvFits} list, or \code{NULL}.
#' @examples
#' U <- list(shared = diag(3), singleton = matrix(0.3, 3, 3) + diag(0.7, 3))
#' mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))
#' cvFits(mp)
#' @export
setGeneric("cvFits", function(x, ...) standardGeneric("cvFits"))

#' @title Get Method Names
#' @description Extract the method name(s) an object was produced with: one
#'   per entry for a collection, a single name for one estimate.
#' @param x A \code{QtlFineMappingResult}, \code{GwasFineMappingResult},
#'   \code{TwasWeights}, \code{CtwasResult} or \code{H2Estimate} object.
#' @return Character vector.
#' @examples
#' data(h2EstimateExample)
#' methodNames(h2EstimateExample)
#' @export
setGeneric("methodNames", function(x) standardGeneric("methodNames"))

#' @title Get Data Type
#' @description Extract the data-type tag.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A character vector or NULL.
#' @examples
#' twe <- TwasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5),
#'   weightStandardized = FALSE)
#' weightsDataType(twe)
#' @export
setGeneric("weightsDataType", function(x, ...) {
    standardGeneric("weightsDataType")
})

# =============================================================================
# Published-output generics (VCF/BCF and BED)
# =============================================================================

#' Write a result to VCF/BCF
#'
#' Writes GWAS summary statistics, fine-mapping results or TWAS weights to a
#' GWAS-VCF-style file: one record per variant, one SAMPLE per study / fitted
#' model, and the per-variant statistics as FORMAT fields. The file is written
#' with VariantAnnotation's \code{\link[VariantAnnotation]{writeVcf}}, so that
#' package must be installed.
#' Supports bgzipped VCF (.vcf.gz/.vcf.bgz); BCF (.bcf) additionally needs a
#' working \code{Rsamtools::asBcf()}, which recent Rsamtools releases disable,
#' so prefer .vcf.bgz.
#'
#' Which FORMAT fields are written depends on the input:
#' \describe{
#'   \item{\code{GwasSumStats}}{\code{ES} / \code{SE} / \code{EZ} /
#'     \code{NS} / \code{AF} / \code{MAF} for one study, whichever its
#'     variants carry.}
#'   \item{a fine-mapping result}{the marginal fields (\code{ES},
#'     \code{SE}, \code{EZ}, \code{LP}, \code{NS}), the fit's own
#'     (\code{AF}, \code{WGT}, \code{WGTSE}, \code{THETA}) and the
#'     posterior ones (\code{PIP}, \code{LBF}, \code{LFSR}, and a
#'     \code{CS}/\code{PUR} pair per credible-set coverage present).}
#'   \item{\code{TwasWeights}}{\code{WGT}, the prediction weight.}
#' }
#'
#' Every fit the selectors match becomes a SAMPLE of ONE file: the records are
#' the union of their variants, the FORMAT schema is the union of their fields,
#' and a fit that did not include a variant -- or does not report a field --
#' leaves those cells empty. So two methods over a region, or a dozen weight
#' methods for one gene, are compared by reading one file rather than a dozen.
#' Each sample is described in its own \code{##SAMPLE} header line, carrying
#' the fit's identity tuple and how many variants it contributed.
#' Cross-validation results are not written -- see \code{\link{writeBed}} for
#' why -- and writing requires the full-data fit.
#'
#' A GWAS study fine-mapped block by block is ONE fit, so its blocks are
#' joined into one sample. A \code{BID} FORMAT field then names each
#' variant's LD block -- credible-set indices are numbered within a block,
#' so a set is the \code{(BID, CS)} pair -- and the per-block header values
#' (block id, analysis window and whether it was given, \code{V},
#' \code{cs_index}) are packed in block order, separated by \code{|}.
#' \code{\link{readFineMappingResult}} splits such a sample back into its
#' blocks.
#'
#' \code{splitByContext} / \code{splitByTrait} opt out of that and write one
#' file per context / trait, at a path decorated with the axis value, holding
#' the fits that share it. \code{splitByBlock} writes one file per LD block
#' of a GWAS result, each holding every fit's part of that block.
#'
#' @param obj A \code{GwasSumStats}, \code{QtlFineMappingResult},
#'   \code{GwasFineMappingResult} or \code{TwasWeights} object, or a
#'   \code{list} of fine-mapping results or of \code{TwasWeights} objects to
#'   combine first.
#' @param filename File path for output. Extension determines format:
#'   \code{.vcf.gz} or \code{.vcf.bgz} for bgzipped VCF, \code{.bcf} for BCF,
#'   \code{.vcf} for uncompressed VCF.
#' @param sampleName Name(s) for the VCF sample column(s). Default \code{NULL}
#'   derives them from the object's identity tuple.
#' @param studyName Character or \code{NULL}. Restrict the written records to
#'   this study; \code{NULL} includes all studies.
#' @param context Character or \code{NULL}. Restrict the written records to this
#'   context; \code{NULL} includes all contexts.
#' @param trait Character or \code{NULL}. Restrict the written records to this
#'   trait; \code{NULL} includes all traits.
#' @param method Character or \code{NULL}. Restrict the written records to this
#'   method; \code{NULL} includes all methods.
#' @param splitByContext Logical. If \code{TRUE}, write one VCF per context.
#'   Default \code{FALSE}.
#' @param splitByTrait Logical. If \code{TRUE}, write one VCF per trait. Default
#'   \code{FALSE}.
#' @param splitByBlock Logical. \code{GwasFineMappingResult} only: if
#'   \code{TRUE}, write one VCF per LD block instead of joining a study's
#'   blocks into one file. Default \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#' @return Invisible path(s) to the written file(s).
#' @seealso \code{\link{writeBed}}, \code{\link{readFineMappingResult}}
#' @examples
#' data(gwasSumStatsS4Example)
#' writeGwasVcf(gwasSumStatsS4Example, tempfile(fileext = ".vcf"))
#' @export
setGeneric("writeGwasVcf", function(obj, filename, ...) {
    standardGeneric("writeGwasVcf")
})

#' Write a result to a published BED table
#'
#' Writes fine-mapping results or TWAS weights to a published FunGen-xQTL BED
#' table, bgzipping and tabix-indexing the output when the path ends in
#' \code{.gz} / \code{.bgz}. The tabular counterpart of
#' \code{\link{writeGwasVcf}}: the VCF is one record per variant
#' with the models as samples, while these tables are long-format, one row per
#' published observation.
#'
#' A \code{TwasWeights} writes one table, a row per (variant, weight column).
#' A fine-mapping result offers two, selected with \code{schema}:
#' \describe{
#'   \item{\code{"topLoci"}}{(default) the posterior view joined to the
#'     marginal one, a row per (variant, context).}
#'   \item{\code{"fsusie"}}{the functional view, a row per (context,
#'     credible set), carrying the fitted curve at both resolutions
#'     \code{levels} selects.}
#' }
#'
#' Every table identifies the fit a row came from -- \code{study_ID},
#' \code{gene_ID}, \code{event_ID} and (where it spans methods)
#' \code{method} -- because these collections are keyed by that tuple and
#' several of them are routinely folded into one output.
#'
#' Cross-validation results are NOT written, in any schema. They describe how
#' a model scored on held-out FOLDS; nothing about them is per-variant, and a
#' published file has no data-secure way to express anything fold-level,
#' because fold membership is sample-level information and sample names can be
#' protected. They remain on the object, where \code{\link{cvResult}} returns
#' them and the caller tabulates whatever they need. For the same reason
#' writing REQUIRES the full-data fit: a cross-validation-only collection
#' (what \code{twasWeightsPipeline(fitFullData = FALSE)} produces) is refused
#' rather than written out empty.
#'
#' @param obj A \code{QtlFineMappingResult}, \code{GwasFineMappingResult} or
#'   \code{TwasWeights} object, or a \code{list} of them, to combine first.
#' @param filename Output path. A \code{.gz} / \code{.bgz} extension selects
#'   bgzip + tabix (via Rsamtools, so no htslib binary is needed); any other
#'   extension writes a plain table.
#' @param schema Fine-mapping results only: which published schema to
#'   write, \code{"topLoci"} (default) or \code{"fsusie"}.
#' @param coverage Which credible-set coverage levels to publish, as
#'   proportions (\code{c(0.95, 0.5)}). \code{NULL} (default) publishes every
#'   level the input carries. Which levels exist is fixed at FIT time by
#'   \code{\link{CredibleSetParam}} and the stored \code{cs_<n>} columns are
#'   read rather than recomputed, so this selects a subset and errors --
#'   naming what is available -- for a level the fit never computed.
#' @param levels \code{"fsusie"} schema only: which resolution(s) of the
#'   fitted effect to write. fSuSiE carries the same effect twice --
#'   \code{"wavelet"}, the curve on the fit's own grid
#'   (\code{\link{fsusieCredibleBand}}), and \code{"trait"}, that curve
#'   reinterpolated onto the positions of the traits the fit ran over
#'   (\code{\link{fsusieTraitEffects}}). They are complementary, so both are
#'   written by default; naming one drops the other's columns, leaving a strict
#'   subset of the full table.
#' @param signalCutoff PIP cutoff passed to \code{\link{topLoci}}
#'   (fine-mapping results, \code{"topLoci"} schema only). Default 0.025.
#' @param minPurity Optional credible-set purity cutoff passed to
#'   \code{\link{topLoci}} (fine-mapping results, \code{"topLoci"} schema
#'   only). Default \code{NULL} for no purity filter.
#' @param columns Which of the SOURCE view's own variant-level columns to
#'   append after the published schema. \code{NULL} (default) appends none;
#'   \code{"all"} appends every column the view carries that the schema does
#'   not already write; a character vector appends exactly those. The
#'   published columns always lead, so adding to them never moves an existing
#'   one. The fine-mapping views carry a long tail worth asking for
#'   (\code{grange_start} / \code{grange_end}, the legacy \code{gene} /
#'   \code{event} aliases, and whatever a custom build added); a
#'   \code{TwasWeights} element carries only its alleles and weight, so there
#'   is nothing further to add unless the collection was built by hand.
#' @param extraColumns Optional NAMED list of trailing columns to append after
#'   the published schema. Each element is either length 1 (recycled across
#'   every row) or one value per row. For fields no S4 object can know -- the
#'   legacy top-loci BED's \code{block_id}, parsed from the source filename,
#'   is the motivating case. Names that collide with a published column are an
#'   error rather than a silent overwrite.
#' @param ... Additional arguments passed to methods.
#' @return \code{filename}, invisibly.
#' @seealso \code{\link{writeGwasVcf}},
#'   \code{\link{readFineMappingResult}}, \code{\link{genomicOrder}}
#' @examples
#' data(qtlFineMappingExample)
#' writeBed(qtlFineMappingExample, tempfile(fileext = ".bed"))
#' writeBed(
#'     qtlFineMappingExample, tempfile(fileext = ".bed"),
#'     schema = "topLoci", columns = "all"
#' )
#' data(twasWeightsExample)
#' writeBed(twasWeightsExample, tempfile(fileext = ".bed"))
#' @export
setGeneric("writeBed", function(obj, filename, ...) {
    standardGeneric("writeBed")
})

# =============================================================================
# QtlDataset accessor generics
# =============================================================================

#' @title Get Study Identifier
#' @description Return the study identifier carried by a \code{QtlDataset}.
#' @param x A \code{QtlDataset} object.
#' @return Character (length 1).
#' @export
setGeneric("studyName", function(x) standardGeneric("studyName"))

#' @title Get Context Names
#' @description Return the names of all contexts carried by an object (e.g., the
#'   keys of the \code{phenotypes} list on a \code{QtlDataset}, or the unique
#'   \code{context} values of a \code{QtlSumStats}).
#' @param x The object.
#' @return Character vector of context names.
#' @examples
#' data(qtlDatasetExample)
#' contexts(qtlDatasetExample)
#' @export
setGeneric("contexts", function(x) standardGeneric("contexts"))

#' @title Get Unique Trait Names
#' @description Return the unique trait identifiers carried by a collection
#'   class (e.g., \code{QtlSumStats}), the trait names of an
#'   \code{\link{SldscData}}, or the single trait label an
#'   \code{H2Estimate} was computed for.
#' @param x The object.
#' @return Character vector of unique trait names; length 1 for an
#'   \code{H2Estimate}.
#' @examples
#' data(qtlSumStatsExample)
#' traitNames(qtlSumStatsExample)
#' @export
setGeneric("traitNames", function(x) standardGeneric("traitNames"))

#' @title Get Per-Row Analysis Windows
#' @description Return the genomic interval each row was analysed over, as a
#'   \code{GRanges} with one range per row -- the trait's own region for a
#'   \code{TwasWeights}; for a fine-mapping result the window the fit was run
#'   over (the requested region, or the trait's cis-window when the fit was
#'   driven by \code{traitId} + \code{cisWindow}); for summary statistics the
#'   LD block or window the row was loaded for (the \code{region} given to the
#'   constructor or loader), else the span of its variants before QC. This is
#'   location provenance (e.g. for cTWAS LD-block placement). One method
#'   serves every collection.
#'
#'   This is the window that was ANALYSED, not the span of the variants that
#'   survived into the fit: QC and screening drop variants first, so a fit over
#'   \code{chr22:10000000-18000000} may contain variants spanning only
#'   \code{chr22:10414272-16916492}. The window is stored because it cannot be
#'   recovered from the object; the span can, and is simply \code{range(x)}.
#'
#'   Being a record of what was analysed, it deliberately does NOT change under
#'   subsetting -- a subset of a fit's variants was
#'   still produced by a fit over the original window. Use \code{range(x)} when
#'   you want the current content instead.
#'
#'   Collections built before the window was recorded fall back to the derived
#'   span, so the accessor always returns a range.
#' @param x The object.
#' @param ... Reserved for future use.
#' @return A \code{GRanges} with one range per row of \code{x}, or an empty
#'   \code{GRanges} when the collection has no rows.
#' @seealso \code{\link{traitPosition}} for the trait's own coordinates.
#' @export
setGeneric("genomicRegion", function(x, ...) standardGeneric("genomicRegion"))

#' @title Get The LD Block Identifier
#' @description Return the LD-block identifier of each row of a collection
#'   that keys its rows by block --- a \code{GwasSumStats}, or a
#'   \code{GwasFineMappingResult} fine-mapped over one.
#'
#'   Per ROW, like \code{\link{genomicRegion}} and unlike
#'   \code{\link{studyName}} / \code{\link{methodNames}}: a block id is
#'   provenance for a particular fit, not an axis the collection is indexed
#'   on, so the returned vector is always the same length as \code{nrow(x)}
#'   and can be used to label or join rows directly.
#'
#'   It keys the external block manifest the blocks were defined in. The
#'   block's BOUNDARIES are the \code{region} column, recorded as given or
#'   computed in \code{regionGiven} (see \code{\linkS4class{GwasSumStats}})
#'   --- unlike \code{range(x)}, which is the span of the variants that
#'   survived QC. A collection with no block axis (every QTL flavour)
#'   returns \code{NA} for every row rather than erroring, so code that
#'   handles both flavours needs no class test.
#' @param x The object.
#' @param ... Reserved for future use.
#' @return Character, one value per row of \code{x} (\code{NA} where the row
#'   records none).
#' @seealso \code{\link{genomicRegion}} for the window that was analysed.
#' @examples
#' data(gwasFineMappingExample, qtlFineMappingExample)
#' blockId(gwasFineMappingExample)
#' # A QTL result has no block axis, so one NA per row rather than an error.
#' blockId(qtlFineMappingExample)
#' @export
setGeneric("blockId", function(x, ...) standardGeneric("blockId"))

#' @title Get Per-Row Trait Positions
#' @description Return the molecular feature's OWN genomic coordinates (the gene
#'   / peak range; TSS = \code{start()}) for each row of a per-tuple collection,
#'   as a \code{GRanges} with one range per row. Distinct from
#'   \code{\link{genomicRegion}} (the window a fine-mapping result was fit
#'   over): the true trait position cannot be inferred from
#'   summary statistics, so it is threaded from the QtlDataset \code{rowRanges}
#'   or an explicit QtlSumStats trait-position and carried as provenance onto
#'   \code{QtlFineMappingResult} and \code{TwasWeights}.
#' @param x The object.
#' @param ... Reserved for future use.
#' @param trait Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @return A \code{GRanges}: one range per row of a collection, or an empty
#'   \code{GRanges} when the collection carries no trait-position provenance
#'   (or, for a \code{QtlSumStats} with \code{trait} given, none of those
#'   traits). For a \code{QtlDataset}, one range per requested trait, named
#'   by trait; a trait in no context is an error.
#' @export
setGeneric("traitPosition", function(x, ...) {
    standardGeneric("traitPosition")
})

#' @title Get Residualized Genotypes
#' @description Residualize the genotype matrix against the per-context
#'   phenotype covariates and the genotype covariates, over the whole
#'   dataset: narrow it first, with \code{x[i, j, ]}, to residualize part
#'   of it.
#' @param x A \code{QtlDataset} object.
#' @param ... Passed on to methods.
#' @param context Character vector. Context(s) whose data to extract.
#' @param residualizationParam What the data is residualized on, built with
#'   \code{\link{ResidualizationParam}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} choose the covariate sets,
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each
#'   side on or off, and \code{covariateNaAction} says how missing
#'   covariate values are handled. The same record the pipelines take, so a
#'   setting means the same thing by either route.
#' @return A numeric matrix (samples x variants).
#' @examples
#' data(qtlDatasetExample)
#' residualizedGenotypes(qtlDatasetExample, context = "brain")
#' @export
setGeneric("residualizedGenotypes", function(x, ...) {
    standardGeneric("residualizedGenotypes")
})

#' @title Get Residualized Phenotypes
#' @description Residualize the per-context phenotype matrices against the
#'   per-context phenotype covariates and the genotype covariates, for one or
#'   more requested contexts. Missing values and outliers are left as they
#'   are: they are a run's phenotype QC (\code{\link{PhenotypeFilterParam}}),
#'   applied by the pipelines.
#' @param x A \code{QtlDataset} object.
#' @param ... Passed on to methods.
#' @param context Character vector. Context(s) whose data to extract.
#' @param trait Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param residualizationParam What the data is residualized on, built with
#'   \code{\link{ResidualizationParam}}: \code{phenotypeCovariates} and
#'   \code{genotypeCovariates} choose the covariate sets,
#'   \code{residualizePhenotype} / \code{residualizeGenotype} turn each
#'   side on or off, and \code{covariateNaAction} says how missing
#'   covariate values are handled. The same record the pipelines take, so a
#'   setting means the same thing by either route.
#' @return A named list of numeric matrices keyed by context.
#' @examples
#' data(qtlDatasetExample)
#' residualizedPhenotypes(qtlDatasetExample, context = "brain")
#' @export
setGeneric("residualizedPhenotypes", function(x, ...) {
    standardGeneric("residualizedPhenotypes")
})

#' @title Get Per-Context Phenotype Covariates
#' @description Return per-context phenotype covariate matrices, taken from the
#'   \code{colData} of each context's \code{SummarizedExperiment}.
#' @param x A \code{QtlDataset} object.
#' @param context Character vector of context names (subset of
#'   \code{names(molecularTraits(x))}).
#' @return A named list of matrices keyed by context.
#' @examples
#' data(qtlDatasetExample)
#' phenotypeCovariates(qtlDatasetExample, context = "brain")
#' @export
setGeneric("phenotypeCovariates", function(x, context) {
    standardGeneric("phenotypeCovariates")
})

#' @title Get Genotype Covariates
#' @description Return the single genotype-derived covariate matrix carried by a
#'   \code{QtlDataset} (e.g., ancestry PCs).
#' @param x A \code{QtlDataset} object.
#' @return Numeric matrix (samples x covariates).
#' @examples
#' data(qtlDatasetExample)
#' genotypeCovariates(qtlDatasetExample)
#' @export
setGeneric("genotypeCovariates", function(x) {
    standardGeneric("genotypeCovariates")
})

# =============================================================================
# GenotypeHandle / LD-statistic / LdMatrix / Annotation / H2Estimate accessors
# =============================================================================

#' @title Get SNP Info
#' @description Return the cached SNP metadata data.frame (columns: SNP, CHR,
#'   BP, A1, A2, fileIdx, optionally MAF).
#' @param x An object carrying cached SNP metadata.
#' @return A data.frame.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly. This table is the READ-PATH
#'   view: unlike the panel's \code{rowRanges()} it carries \code{fileIdx},
#'   the on-disk variant position that keeps reads correct after a handle has
#'   been row-subset. Callers that only need chrom/pos/alleles should use
#'   \code{rowRanges()} on the panel instead.
#' @keywords internal
setGeneric("snpInfo", function(x) standardGeneric("snpInfo"))

#' @title Get Genotype Storage Format
#' @description Return the detected genotype storage format.
#' @param x A \code{GenotypeHandle}.
#' @return Character (length 1): one of "plink1", "plink2", "bgen".
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("genotypeFormat", function(x) standardGeneric("genotypeFormat"))

#' @title Get Sample Identifiers
#' @description Return the sample-id vector.
#' @param x A \code{GenotypeHandle}.
#' @return Character vector.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("sampleIds", function(x) standardGeneric("sampleIds"))

#' @title Get plink2 pgen Pointer
#' @description Return the cached external pointer to the plink2 pgen handle
#'   (NULL when the handle is not pgen-backed).
#' @param x A \code{GenotypeHandle}.
#' @return An external pointer or NULL.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("pgenPtr", function(x) standardGeneric("pgenPtr"))

#' @title Get Per-Block Eigendecompositions
#' @description Return the per-block eigendecompositions of an
#'   \code{LdEigen}: the stored ones, or, for one built from a genotype
#'   panel (which stores none), every block's computed from the panel now.
#'   That holds all of them in memory at once, which the estimators avoid by
#'   decomposing one block at a time.
#' @param x An \code{LdEigen}.
#' @return A list with one element per block, each a list with
#'   \code{values}, \code{vectors} and \code{snpIdx}. One computed from a
#'   panel also carries \code{nNull}: the number of the block's directions
#'   whose eigenvalue is exactly zero, which it does not list.
#' @examples
#' data(ldEigenExample)
#' eigenList(ldEigenExample)
#' @export
setGeneric("eigenList", function(x) standardGeneric("eigenList"))

#' @title Get LD Reference Panel Size
#' @description Return the reference-panel sample size used to compute an
#'   \code{LdStatistic}.
#' @param x An \code{LdStatistic}.
#' @return Integer (length 1).
#' @examples
#' data(ldEigenExample)
#' ldPanelNSamples(ldEigenExample)
#' @export
setGeneric("ldPanelNSamples", function(x) standardGeneric("ldPanelNSamples"))

#' @title Get In-Sample Flag
#' @description Whether the LD reference panel is from the same cohort as the
#'   GWAS (affects bias correction).
#' @param x An \code{LdStatistic}.
#' @return Logical (length 1).
#' @examples
#' data(ldEigenExample)
#' inSample(ldEigenExample)
#' @export
setGeneric("inSample", function(x) standardGeneric("inSample"))

#' @title Get Per-Block LD Matrix List
#' @description Return the per-block LD correlation matrices g-LDSC forms its
#'   FGLS residual covariance from: the stored ones, or, for an
#'   \code{LdScore} built from a genotype panel (which stores none), every
#'   block's computed from the panel now. That holds all of them in memory at
#'   once, which the estimators avoid by reading one block at a time.
#' @param x An \code{LdScore}.
#' @return A list with one element per block, each a list with \code{R} and
#'   \code{snpIdx}; empty for a scores-only \code{LdScore}.
#' @examples
#' data(ldScoreExample)
#' ldMatrixList(ldScoreExample)
#' @export
setGeneric("ldMatrixList", function(x) standardGeneric("ldMatrixList"))

#' @title Get LD Block Container
#' @description Return the \code{GRanges} of LD blocks carried by an
#'   \code{LdStatistic}.
#' @param x An \code{LdStatistic}.
#' @return A \code{GRanges} of LD block intervals.
#' @examples
#' data(ldEigenExample)
#' ldBlocks(ldEigenExample)
#' @export
setGeneric("ldBlocks", function(x) standardGeneric("ldBlocks"))

#' @title Get GenotypeHandle
#' @description Return the \code{GenotypeHandle} behind a \code{QtlDataset}'s
#'   genotype experiment.
#' @param x A \code{QtlDataset}.
#' @return A \code{GenotypeHandle}, or NULL for genotypes read through no
#'   handle (a VCF, GDS or PLINK1 file, or held in memory);
#'   \code{genotypes()} extracts a dosage block instead.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("genotypeHandle", function(x) {
    standardGeneric("genotypeHandle")
})

#' @title Get Mixture Weights
#' @description Return the weights a mixture's components are combined
#'   with: an \code{\link{LdMatrix}} mixture's (\code{NULL} for a single
#'   panel's LD), or a \code{\link{GenotypeMixture}}'s, named after its
#'   panels.
#' @param x An \code{LdMatrix} or a \code{GenotypeMixture}.
#' @return Numeric vector or NULL.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' half <- seq_len(nrow(X) / 2)
#' mix <- LdMatrix(
#'   list(LdMatrix(cor(X[half, ])), LdMatrix(cor(X[-half, ]))),
#'   mixtureWeights = c(0.75, 0.25))
#' mixtureWeights(mix)
#' @export
setGeneric("mixtureWeights", function(x) {
    standardGeneric("mixtureWeights")
})

#' @title Get Per-Block tau Matrix
#' @description Return the per-block jackknife tau matrix carried by an
#'   \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return A numeric matrix or NULL.
#' @examples
#' data(h2EstimateExample)
#' annotationJackknifeCoefs(h2EstimateExample)
#' @export
setGeneric("annotationJackknifeCoefs", function(x) {
    standardGeneric("annotationJackknifeCoefs")
})

#' @title Get Global SNP Heritability
#' @description Return the global SNP heritability estimate carried by an
#'   \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @examples
#' data(h2EstimateExample)
#' heritabilityEstimate(h2EstimateExample)
#' @export
setGeneric("heritabilityEstimate", function(x) {
    standardGeneric("heritabilityEstimate")
})

#' @title Get Heritability Standard Error
#' @description Return the standard error of the global SNP heritability
#'   estimate carried by an \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @examples
#' data(h2EstimateExample)
#' heritabilityEstimateSe(h2EstimateExample)
#' @export
setGeneric("heritabilityEstimateSe", function(x) {
    standardGeneric("heritabilityEstimateSe")
})

#' @title Get LD-Score Regression Intercept
#' @description Return the regression intercept carried by an
#'   \code{H2Estimate}. Values above 1 indicate confounding or sample overlap
#'   rather than polygenic signal.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @examples
#' data(h2EstimateExample)
#' heritabilityIntercept(h2EstimateExample)
#' @export
setGeneric("heritabilityIntercept", function(x) {
    standardGeneric("heritabilityIntercept")
})

#' @title Get Intercept Standard Error
#' @description Return the standard error of the regression intercept carried
#'   by an \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @examples
#' data(h2EstimateExample)
#' heritabilityInterceptSe(h2EstimateExample)
#' @export
setGeneric("heritabilityInterceptSe", function(x) {
    standardGeneric("heritabilityInterceptSe")
})

# Internal generics for the unified joint-analysis engine (see R/JointGroup.R
# and dev/jointSpecification-s4-refactor.md). Not exported: the engine and its
# fitters are package-internal machinery.

# fitJointGroup(group, pipeline, token, args) -- multiple dispatch on
# (JointGroup subclass, JointPipeline subclass). The 4 irreducible joint fits
# (individual/sumstats x fm/twas). Returns one fit entry (FineMappingRow or
# TwasWeightsRow).
setGeneric("fitJointGroup", function(group, pipeline, token, args) {
    standardGeneric("fitJointGroup")
})

# construct(pipeline, records) -- assemble the per-pipeline result collection
# (QtlFineMappingResult vs TwasWeights) from the driver's list of per-row
# records. The joint row identity (which axes collapse to "joint" +
# jointStudies/Contexts/Traits) is carried on each record, derived from each
# group's `conditions`.
setGeneric("construct", function(pipeline, records, ...) {
    standardGeneric("construct")
})

# ---- SldscData accessors ----
#' @title Get the annotation table from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A \code{DataFrame} of annotations (CHR, SNP, annotation columns).
#' @rdname annotData
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     annotationJackknifeCoefs = matrix(1e-7, 10, n, dimnames = list(NULL,
#'   cats)), nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annotData = annot, frqData = frq, traits = traits)
#' annotData(sd)
#' @export
setGeneric("annotData", function(x) standardGeneric("annotData"))

#' @title Get the allele-frequency table from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A \code{DataFrame} of reference-panel frequencies (SNP, MAF).
#' @rdname frqData
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     annotationJackknifeCoefs = matrix(1e-7, 10, n, dimnames = list(NULL,
#'   cats)), nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annotData = annot, frqData = frq, traits = traits)
#' frqData(sd)
#' @export
setGeneric("frqData", function(x) standardGeneric("frqData"))

#' @title Get the annotation column names from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A character vector of annotation column names.
#' @rdname annotCols
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     annotationJackknifeCoefs = matrix(1e-7, 10, n, dimnames = list(NULL,
#'   cats)), nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annotData = annot, frqData = frq, traits = traits)
#' annotCols(sd)
#' @export
setGeneric("annotCols", function(x) standardGeneric("annotCols"))

#' @title Get Baseline Annotations
#' @description Extract only baseline-tier annotations from an
#'   \code{AnnotationMatrix}.
#' @param x An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only baseline annotations.
#' @export
setGeneric("baselineAnnotations", function(x) {
    standardGeneric("baselineAnnotations")
})

#' @title Get Tested Annotations
#' @description Extract only candidate-tier annotations from an
#'   \code{AnnotationMatrix} -- the annotations whose enrichment is being
#'   tested, as opposed to the baseline set they are tested against.
#' @param x An \code{AnnotationMatrix} object.
#' @return An \code{AnnotationMatrix} with only candidate annotations.
#' @export
setGeneric("testedAnnotations", function(x) {
    standardGeneric("testedAnnotations")
})

#' @title Get sLDSC Results from an SldscData
#' @description Return the stratified LD-score regression runs held by an
#'   \code{\link{SldscData}}, keyed by trait. Each trait's entry holds its
#'   \code{single} runs (one per annotation) and its \code{joint} run, read
#'   with \code{[[}: \code{sldscResults(x, "T")[["T"]]$joint}.
#' @param x An \code{\link{SldscData}} object.
#' @param trait Character trait name(s), or \code{NULL} (default) for every
#'   trait. A name the object does not hold is an error.
#' @param ... Unused.
#' @return A named list keyed by trait, each element a list with
#'   \code{single} (a list of runs) and \code{joint} (one run, or
#'   \code{NULL} when the trait has none).
#' @rdname sldscResults
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     annotationJackknifeCoefs = matrix(1e-7, 10, n, dimnames = list(NULL,
#'   cats)), nBlocks = 10L)
#' }
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annotData = annot, frqData = frq, traits = traits)
#' sldscResults(sd, "traitX")[["traitX"]]$joint$categories
#' @export
setGeneric("sldscResults", function(x, trait = NULL, ...) {
    standardGeneric("sldscResults")
})

#' @title Colocalization Views
#' @description The five views over a \code{\linkS4class{ColocResult}}. The
#'   object stores one element per tested (QTL credible set, GWAS credible set,
#'   block) pair; these accessors project it to the granularity a given
#'   analysis needs, rather than any one of them being the stored shape.
#'
#'   \describe{
#'     \item{\code{colocPairs()}}{One row per tested pair -- coloc's native
#'       \code{$summary} granularity. Every testable pair is reported verbatim,
#'       including both halves of a credible set that straddles a block
#'       boundary; interpreting them is the caller's decision.}
#'     \item{\code{colocVariants()}}{One range per (pair, variant), carrying
#'       \code{colocPp = PP.H4.abf * SNP.PP.H4} -- the posterior that this
#'       variant is the shared causal one.}
#'     \item{\code{colocPooledVariants()}}{One range per (gene, variant),
#'       pooled across pairs by the rule below. A
#'       \code{\linkS4class{ColocBoostResult}} pools per (analysis, variant)
#'       instead, with the noisy-OR stage only.}
#'     \item{\code{colocCredibleSets()}}{One element per coloc credible set:
#'       the smallest set of variants whose cumulative \code{SNP.PP.H4}
#'       reaches \code{coverage}, as ranges carrying their \code{SNP.PP.H4},
#'       with the set's pair columns, \code{csSize}, \code{csCoverage},
#'       \code{leadVariant}, \code{leadPp} and \code{purity} in
#'       \code{mcols()}. This is a THIRD variant set, not guaranteed to be a
#'       subset of either input credible set, so its purity is recomputed from
#'       LD rather than inherited.}
#'     \item{\code{colocGenes()}}{One row per gene, pooled across pairs.}
#'   }
#'
#' @section Pooling: For a fixed QTL credible set \emph{i}, the per-(block,
#'   GWAS credible set) results \emph{j} are mutually exclusive -- the shared
#'   causal variant is in one block or the other -- so they are summed:
#'   \code{min(1, sum_j PP.H4_ij)}, with a warning when the sum exceeds 1
#'   (several GWAS credible sets competing for one QTL signal, which is
#'   diagnostic). Distinct QTL credible sets are independent signals and
#'   combine as \code{1 - prod_i (1 - p_i)}.
#'
#' @param x A \code{ColocResult}.
#' @param coverage Cumulative \code{SNP.PP.H4} the credible set must reach.
#' @param minPp4 Optional lower bound on the pair's \code{PP.H4.abf}. Applied
#'   before any LD work, so a stricter threshold costs strictly less.
#' @param requireMaxH4 Logical (length 1), default \code{FALSE}. When
#'   \code{TRUE}, keep only pairs whose \code{PP.H4.abf} is the largest of
#'   the five posterior probabilities, i.e. colocalization is the favoured
#'   hypothesis rather than merely a probable one.
#' @param minAbsCorr Minimum absolute correlation between credible-set members
#'   (purity). Requires the object's \code{ldReference}; when that is absent,
#'   purity is reported as \code{NA} and no set is dropped.
#' @param ... Additional arguments passed on to methods.
#' @return \code{colocVariants()} and \code{colocPooledVariants()} a
#'   \code{GRanges}, its columns in \code{mcols()};
#'   \code{colocCredibleSets()} a \code{GRangesList}, one element per set;
#'   the other views a \code{DataFrame}.
#' @name colocViews
#' @rdname colocViews
#' @examples
#' data(qtlFineMappingLbfExample)
#' data(gwasFineMappingLbfExample)
#' res <- colocPipeline(
#'     qtlFineMappingResult = qtlFineMappingLbfExample,
#'     gwasInput = gwasFineMappingLbfExample
#' )
#' head(colocPairs(res))
#' head(colocVariants(res))
#' head(colocPooledVariants(res))
#' colocGenes(res)
#' head(colocCredibleSets(res, minPp4 = 0.001))
NULL

#' @rdname colocViews
#' @export
setGeneric("colocPairs", function(x, ...) {
    standardGeneric("colocPairs")
})

#' @rdname colocViews
#' @export
setGeneric("colocVariants", function(x, ...) {
    standardGeneric("colocVariants")
})

#' @rdname colocViews
#' @export
setGeneric("colocPooledVariants", function(x, ...) {
    standardGeneric("colocPooledVariants")
})

#' @rdname colocViews
#' @export
setGeneric(
    "colocCredibleSets",
    function(x, coverage = 0.95, minPp4 = NULL, minAbsCorr = 0.8, ...) {
        standardGeneric("colocCredibleSets")
    }
)

#' @rdname colocViews
#' @export
setGeneric("colocGenes", function(x, ...) {
    standardGeneric("colocGenes")
})

#' @title Get Per-Chromosome Payload Paths
#' @description The chromosome-to-path map of a sharded
#'   \code{\link{GenotypeHandle}}. Empty (\code{character(0)}) for a
#'   single-file handle, whose one path is \code{\link{path}} instead; a
#'   non-empty map marks a one-file-per-chromosome handle whose extraction is
#'   routed by chromosome.
#' @param x A \code{GenotypeHandle}.
#' @param ... Additional arguments passed on to methods.
#' @return A named character vector mapping each chromosome, in UCSC style
#'   (\code{"chr21"}), to its path or prefix, possibly empty.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("chromPaths", function(x, ...) standardGeneric("chromPaths"))

#' Per-variant per-effect log Bayes factors (wide)
#'
#' The variant x effect matrix of single-effect log Bayes factors
#' (\code{lbf_variable}, or fSuSiE's \code{lBF}) from the stored fit: one range
#' per variant, one \code{lbf_L<k>} column per effect. The per-variant scalar
#' summary (max across effects) is the \code{logBF} column of
#' \code{\link{topLoci}}; this accessor keeps the full per-effect breakdown.
#' Across entries with different effect counts the collection method NA-fills
#' the ragged columns.
#'
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Ignored.
#' @return A \code{GRanges}, one range per variant, whose \code{mcols()}
#'   hold \code{variant_id} + \code{lbf_L1..lbf_LL} (the collection method
#'   also carries the entry identity columns).
#' @seealso \code{\link{topLoci}} (the scalar \code{logBF} column)
#' @examples
#' data(qtlFineMappingExample)
#' lbf(qtlFineMappingExample)
#' @export
setGeneric("lbf", function(x, ...) standardGeneric("lbf"))

#' Per-credible-set summary of a fine-mapping result
#'
#' One row per credible set (at a given coverage) with its size, purity, prior
#' variance, log Bayes factor, and lead variant -- the per-CS complement to the
#' per-variant \code{\link{credibleSets}}. Replaces the legacy
#' per-effect `effect.tsv`.
#'
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param coverage Credible-set coverage to summarise. Default 0.95.
#' @param ... Ignored.
#' @return A \code{DataFrame}: \code{cs, effect_id, coverage, n_variants,
#'   purity_min, purity_mean, V, cs_log10bf} (strongest member logBF),
#'   \code{cs_log_bf} (true per-effect single-effect log Bayes factor),
#'   \code{cs_pip} (summed member PIP = inclusion mass captured),
#'   \code{cs_mean_effect} (mean posterior conditional effect),
#'   \code{lead_variant, lead_pip} (the collection method also carries the entry
#'   identity columns).
#' @seealso \code{\link{credibleSets}}, \code{\link{topLoci}}
#' @examples
#' data(qtlFineMappingExample)
#' credibleSetSummary(qtlFineMappingExample)
#' @export
setGeneric("credibleSetSummary", function(x, ...) {
    standardGeneric("credibleSetSummary")
})

#' fSuSiE credible band (fitted effect + uncertainty band)
#'
#' For a functional-SuSiE (\code{fsusieR::susiF}) fine-mapping entry, returns
#' the fitted effect curve and its credible band over the functional grid, one
#' range per (credible set, grid position). Wraps the upstream fSuSiE band
#' computation (fsusieR's internal \code{update_cal_credible_band.susiF} +
#' \code{get_fitted_effect}). Uses saved bands in trimmed fits, or computes
#' them from an untrimmed fit. Returns an empty \code{GRanges} when the
#' functional fields are absent.
#'
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Ignored.
#' @return A \code{GRanges}, one width-1 range per (credible set, grid
#'   position), whose \code{mcols()} hold \code{cs, chrom, pos, effect,
#'   lower, upper, halfwidth} -- \code{halfwidth} being
#'   \code{(upper - lower) / 2}, the symmetric band half-width publication
#'   tables carry (the collection method additionally carries the entry
#'   identity columns).
#' @seealso \code{\link{fsusieTraitEffects}},
#'   \code{\link{fsusieAffectedRegions}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieCredibleBand(fsusieFineMappingExample)
#' @export
setGeneric("fsusieCredibleBand", function(x, ...) {
    standardGeneric("fsusieCredibleBand")
})

#' fSuSiE effect at each joint trait
#'
#' The fitted effect curve evaluated AT the traits the functional fit was run
#' over, one row per (credible set, trait). The curve itself lives on a wavelet
#' grid (\code{\link{fsusieCredibleBand}}); the traits sit at arbitrary
#' positions between grid points, so the effect is linearly interpolated with
#' constant extrapolation past the ends.
#'
#' This is the view that pairs a trait's NAME with its position: the names come
#' from the row's \code{jointTraits} column and the positions from the fit,
#' index-aligned because both describe the same joint group. Returns zero rows
#' when the functional fields are absent.
#'
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Ignored.
#' @return A long \code{DataFrame}: \code{cs, chrom, jointTrait, traitPos,
#'   effect} (the collection method additionally carries the entry identity
#'   columns).
#' @seealso \code{\link{fsusieCredibleBand}},
#'   \code{\link{fsusieAffectedRegions}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieTraitEffects(fsusieFineMappingExample)
#' @export
setGeneric("fsusieTraitEffects", function(x, ...) {
    standardGeneric("fsusieTraitEffects")
})

#' fSuSiE affected genomic regions
#'
#' The sub-intervals of the functional grid where a credible set's band excludes
#' zero (the genomic footprint of each fSuSiE effect), from the upstream
#' \code{fsusieR::affected_reg}. Returns a \code{GRanges} with the credible-set
#' label, its purity, and the effect \code{direction}
#' (\code{"pos"}/\code{"neg"}; upstream discards it). Uses retained functional
#' summaries or an untrimmed fit; empty when those fields are absent.
#'
#' @param x A \code{QtlFineMappingResult} or \code{GwasFineMappingResult}.
#' @param ... Ignored.
#' @return A \code{GRanges} of affected intervals with mcols \code{cs},
#'   \code{purity}, \code{direction} (plus entry identity for the collection).
#' @seealso \code{\link{fsusieCredibleBand}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieAffectedRegions(fsusieFineMappingExample)
#' @export
setGeneric("fsusieAffectedRegions", function(x, ...) {
    standardGeneric("fsusieAffectedRegions")
})

# =============================================================================
# MethodParam accessor generics
# -----------------------------------------------------------------------------
# Every Param is read with `$` / `[[` and written with `$<-` (MethodParam.R),
# which rebuilds the record through its constructor. Two classes also carry
# named getters, because they are records the user HOLDS rather than builds
# and forgets:
#
#   CredibleSetParam     postprocessFinemappingFits() and buildTopLoci() take
#                        an already-computed fit, so a stored fit is
#                        re-summarized under changed settings without
#                        refitting. `L` / `Lgreedy` are deliberately fixed:
#                        they are seeded onto each SuSiE-family token's own
#                        arguments, so changing one needs a refit, and `$<-`
#                        refuses them.
#   GwasFineMappingParam the only Param nesting others.
#
# Bare field names are not available as getter names -- `coverage` is an
# IRanges generic for read depth and `methods` is utils::methods.
# =============================================================================

#' @rdname CredibleSetParam
#' @param x A \code{CredibleSetParam}.
#' @param ... Unused, present for generic consistency.
#' @export
setGeneric("csCoverage", function(x, ...) {
    standardGeneric("csCoverage")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("secondaryCoverage", function(x, ...) {
    standardGeneric("secondaryCoverage")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("signalCutoff", function(x, ...) {
    standardGeneric("signalCutoff")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("minAbsCorr", function(x, ...) {
    standardGeneric("minAbsCorr")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("medianAbsCorr", function(x, ...) {
    standardGeneric("medianAbsCorr")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("includeAllCs", function(x, ...) {
    standardGeneric("includeAllCs")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("perCsColumns", function(x, ...) {
    standardGeneric("perCsColumns")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("maxNumSingleEffects", function(x, ...) {
    standardGeneric("maxNumSingleEffects")
})

#' @rdname CredibleSetParam
#' @export
setGeneric("maxNumSingleEffectsGreedy", function(x, ...) {
    standardGeneric("maxNumSingleEffectsGreedy")
})

#' @rdname GwasFineMappingParam
#' @param x A \code{GwasFineMappingParam}.
#' @param ... Unused, present for generic consistency.
#' @export
setGeneric("fineMappingMethods", function(x, ...) {
    standardGeneric("fineMappingMethods")
})

#' @rdname GwasFineMappingParam
#' @export
setGeneric("credibleSetParam", function(x, ...) {
    standardGeneric("credibleSetParam")
})

#' @rdname GwasFineMappingParam
#' @export
setGeneric("susieRssParam", function(x, ...) {
    standardGeneric("susieRssParam")
})

#' @rdname GwasFineMappingParam
#' @export
setGeneric("panelFilterParam", function(x, ...) {
    standardGeneric("panelFilterParam")
})

#' @rdname GwasFineMappingParam
#' @export
setGeneric("initializeWithSusieInf", function(x, ...) {
    standardGeneric("initializeWithSusieInf")
})

#' @rdname GwasFineMappingParam
#' @export
setGeneric("fitRetention", function(x, ...) {
    standardGeneric("fitRetention")
})
