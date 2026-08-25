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
#' @description Estimate SNP heritability from GWAS summary statistics using one
#'   of three methods: LDER, g-LDSC, or HDL/sHDL.
#' @param sumstats A \code{GwasSumStats} object.
#' @param ldRef An \code{LdStatistic} object (method-appropriate subclass).
#' @param method Character, one of "lder", "gldsc", "hdl".
#' @param annotations An \code{AnnotationMatrix} object, or NULL for
#'   unstratified estimation.
#' @param local Logical, whether to compute per-block local estimates.
#' @param ... Additional method-specific arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @return An \code{H2Estimate} object.
#' @examples
#' data(ldEigenExample)
#' gr <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(seq(50, by = 100, length.out = 20), width = 1))
#' S4Vectors::mcols(gr) <- S4Vectors::DataFrame(SNP = paste0("rs", 1:20),
#'   A1 = "A", A2 = "G", Z = rnorm(20), N = 10000L)
#' panel <- readGenotypes(
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr")
#' )
#' ss <- GwasSumStats(study = "trait1", entry = list(gr),
#'   genome = "hg19", ldSketch = panel)
#' estimateH2(ss, ldEigenExample, method = "lder")
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
#'   annotations.
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
#' @description Read genotype data (VCF, plink1, plink2 or GDS) as a
#'   \code{RangedSummarizedExperiment}: variants in \code{rowRanges},
#'   samples in \code{colData}, and dosages in a \code{DelayedArray} assay
#'   that reads from the file only when something touches it. Nothing is read
#'   here beyond the variant and sample metadata.
#'
#'   The panel is the same shape a \code{QtlDataset} uses for its genotypes
#'   and the same shape a collection carries as its \code{ldSketch}, so the
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
#'   with \code{chroms}), or \code{ldMeta} plus \code{region} to resolve a
#'   block out of an LD-meta table. Exactly one source may be given.
#' @param path Character, path to a single self-describing genotype file
#'   (\code{.vcf}, \code{.vcf.gz}, \code{.vcf.bgz}, \code{.bcf} or
#'   \code{.gds}). Omit it when naming the panel by one of the keyword
#'   sources described above.
#' @param format Character, one of "vcf", "plink1", "plink2", "gds". If NULL,
#'   inferred from file extension.
#' @param ... The keyword source arguments described above, plus any further
#'   arguments forwarded to the format-specific reader.
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
#' @param paths Named character vector of file paths, or a named list of GRanges
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
    function(paths, snpRanges, annotationMeta = NULL, genome = "hg19", ...) {
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
#' getLocal(h2EstimateExample)
#' @export
setGeneric("getLocal", function(object) standardGeneric("getLocal"))

#' @title Get Enrichment Estimates
#' @description Extract annotation enrichment estimates from a result object.
#' @param object An \code{H2Estimate} object.
#' @return A \code{data.frame} of enrichment estimates, or NULL.
#' @examples
#' data(h2EstimateExample)
#' getEnrichment(h2EstimateExample)
#' @export
setGeneric("getEnrichment", function(object) standardGeneric("getEnrichment"))

#' @title Get Score Statistics
#' @description Extract score statistics for candidate annotations.
#' @param object An \code{H2Estimate} object.
#' @return A list with \code{z} and \code{R}, or NULL.
#' @examples
#' data(h2EstimateExample)
#' getScoreStats(h2EstimateExample)
#' @export
setGeneric("getScoreStats", function(object) standardGeneric("getScoreStats"))

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
#' @return Numeric vector of z-scores.
#' @export
setGeneric("getZ", function(x, ...) standardGeneric("getZ"))

#' @title Get Sample Sizes
#' @description Extract sample size vector from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric vector of sample sizes.
#' @export
setGeneric("getN", function(x, ...) standardGeneric("getN"))

#' @title Get Association P-values
#' @description Extract the association p-value vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple. Part of the first-class summary-statistic column set alongside
#'   \code{\link{getZ}} / \code{\link{getBeta}} / \code{\link{getSe}}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric vector of p-values, or \code{NULL} if not available.
#' @export
setGeneric("getP", function(x, ...) standardGeneric("getP"))

#' @title Get Marginal Effect Sizes
#' @description Extract the marginal effect-size (beta) vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric vector of effect sizes, or \code{NULL} if not available.
#' @export
setGeneric("getBeta", function(x, ...) standardGeneric("getBeta"))

#' @title Get Effect-Size Standard Errors
#' @description Extract the effect-size standard-error vector from a
#'   \code{GwasSumStats} or \code{QtlSumStats} entry, selected by its identity
#'   tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric vector of standard errors, or \code{NULL} if not available.
#' @export
setGeneric("getSe", function(x, ...) standardGeneric("getSe"))

#' @title Get Minor Allele Frequencies
#' @description Extract MAF vector from a GwasSumStats object.
#' @param x A \code{GwasSumStats} or \code{QtlDataset} object.
#' @param ... Class-specific selection arguments (e.g., \code{region},
#'   \code{cisWindow} for \code{QtlDataset}).
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param cisWindow Integer or \code{NULL}. Half-width (bp) of the cis window
#'   around the trait; \code{NULL} uses the dataset default.
#' @param samples Character vector or \code{NULL}. Restrict to these sample IDs;
#'   \code{NULL} uses all samples.
#' @return Numeric vector of MAFs, or NULL if not available.
#' @export
setGeneric("getMaf", function(x, ...) standardGeneric("getMaf"))

#' @title Get Effect-Allele Frequencies
#' @description Extract the directional effect-allele (A1) frequency vector for
#'   a \code{QtlDataset}. Unlike \code{\link{getMaf}}, the value is \emph{not}
#'   folded to the minor allele: it is the frequency of the dosage-counted
#'   allele (A1, the effect allele), matching the allele the marginal effect
#'   sizes and the fine-mapping \code{af} column report.
#' @param x A \code{QtlDataset} object.
#' @param ... Class-specific selection arguments (e.g., \code{traitId},
#'   \code{region}, \code{cisWindow}, \code{samples} for \code{QtlDataset}).
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param cisWindow Integer or \code{NULL}. Half-width (bp) of the cis window
#'   around the trait; \code{NULL} uses the dataset default.
#' @param samples Character vector or \code{NULL}. Restrict to these sample IDs;
#'   \code{NULL} uses all samples.
#' @return Named numeric vector of effect-allele frequencies (names are variant
#'   IDs), or an empty vector when no variants are selected.
#' @examples
#' data(qtlDatasetExample)
#' getAf(qtlDatasetExample)
#' @export
setGeneric("getAf", function(x, ...) standardGeneric("getAf"))

#' @title Get Number of SNPs
#' @description Number of SNPs in a \code{GwasSumStats} or \code{QtlSumStats}
#'   entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
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
#' qss <- QtlSumStats(study = rep("s", G), context = rep("brain", G),
#'   trait = paste0("g", seq_len(G)), entry = entries, genome = "hg19",
#'   n_variants = rep(50L, G), n_variants_filtered = rep(30L, G),
#'   p_beta = pBeta, beta_shape1 = rep(1, G), beta_shape2 = rep(200, G))
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
#' qss <- QtlSumStats(study = rep("s", G), context = rep("brain", G),
#'   trait = paste0("g", seq_len(G)), entry = entries, genome = "hg19",
#'   n_variants = rep(50L, G), n_variants_filtered = rep(30L, G),
#'   p_beta = pBeta, beta_shape1 = rep(1, G), beta_shape2 = rep(200, G))
#' pp <- qtlAssociationPostprocess(qss)
#' getSignificantQtls(pp)
#' @export
setGeneric("getSignificantQtls", function(x, ...) {
    standardGeneric("getSignificantQtls")
})

#' @title Subset by Chromosome
#' @description Extract a chromosome-specific subset of a GwasSumStats object.
#' @param x A \code{GwasSumStats} object.
#' @param chr Character, chromosome name (e.g., "1", "chr1").
#' @return A \code{GwasSumStats} object.
#' @examples
#' data(gwasSumStatsS4Example)
#' subsetChr(gwasSumStatsS4Example, "chr22")
#' @export
setGeneric("subsetChr", function(x, chr) standardGeneric("subsetChr"))

#' @title Reconcile Two Collections to Their Shared Variants
#' @description Restrict two fine-mapping collections to the variants they have
#'   in common, renormalizing each retained single-effect posterior over the
#'   variants that survive.
#' @details \code{subsetRegion()} restricts ONE object to a window; this
#'   intersects TWO objects against each other. Coloc needs it symmetric --
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
#'   \code{getRetainedMass()} is what makes it visible.
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
#' @return A \code{data.frame}, one row per (element, effect), with the
#'   element's identity columns, the effect index, and its retained mass.
#' @examples
#' data(qtlFineMappingExample)
#' data(gwasFineMappingExample)
#' # Nothing has been reconciled yet, so there is no mass to report.
#' nrow(getRetainedMass(qtlFineMappingExample))
#' both <- intersectVariants(qtlFineMappingExample, gwasFineMappingExample)
#' getRetainedMass(both$x)
#' @export
setGeneric("getRetainedMass", function(x, ...) {
    standardGeneric("getRetainedMass")
})

#' @title Restrict a Collection to a Region
#' @description Narrow a \code{RangedTupleList} collection to an analysis
#'   window: elements that do not overlap it are dropped, and the survivors are
#'   trimmed to the overlapping variants. Elements that trim to zero variants
#'   are dropped too.
#' @details This is deliberately a new verb rather than a redefinition of a
#'   Bioconductor generic, because three genuinely different operations are in
#'   play and all three are useful:
#'   \describe{
#'     \item{\code{subsetByOverlaps(x, win)}}{selects elements, keeps them
#'       whole.}
#'     \item{\code{restrict(x, start, end)}}{trims elements, keeps empty ones.}
#'     \item{\code{subsetRegion(x, win)}}{selects AND trims AND drops empties.}
#'   }
#'   \code{subsetChr()} is the whole-seqname special case of this verb.
#' @param x A \code{RangedTupleList} collection.
#' @param region A \code{GRanges}, a \code{"chr:start-end"} string, or a
#'   one-row \code{data.frame} with \code{chrom} / \code{start} / \code{end}.
#' @param ... Additional arguments passed to methods. Fine-mapping collections
#'   accept \code{adjust} to renormalize PIPs over the retained variants.
#' @return An object of the same class as \code{x}, with the collection-level
#'   slots carried over and \code{mcols} narrowed to the surviving elements.
#' @examples
#' data(qtlSumStatsMulticontextExample)
#' mc <- qtlSumStatsMulticontextExample
#' sum(lengths(mc))
#' # Every element is narrowed to the window; the collection keeps its shape,
#' # so an element with nothing left stays as an empty one rather than a
#' # dropped row.
#' sub <- subsetRegion(mc, "chr22:14560203-15000000")
#' sum(lengths(sub))
#' nrow(sub)
#' @export
setGeneric("subsetRegion", function(x, region, ...) {
    standardGeneric("subsetRegion")
})

#' @title Get Phenotype Variance
#' @description Extract phenotype variance from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple. Returns
#'   \code{NULL} when the entry has no \code{varY} recorded.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @return Numeric phenotype variance, or NULL.
#' @examples
#' data(qtlSumStatsExample)
#' getVarY(qtlSumStatsExample)
#' @export
setGeneric("getVarY", function(x, ...) standardGeneric("getVarY"))

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
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param ranges A \code{GRanges} or \code{NULL} (default). A selected tuple can
#'   own several elements -- one per chromosome, since a multi-seqname entry is
#'   split at construction -- and they are stitched back into a single
#'   \code{GRanges}. Pass \code{ranges} to pull only the variants overlapping a
#'   region instead of the whole span; \code{NULL} returns everything.
#' @return A \code{GRanges}, a \code{QtlSumStats}, or \code{NULL}.
#' @examples
#' data(qtlSumStatsExample)
#' getSumStats(qtlSumStatsExample)
#' @export
setGeneric("getSumStats", function(x, ...) standardGeneric("getSumStats"))

#' @title Get Standardized Sumstat Data Frame for One Tuple
#' @description Return a per-tuple summary-statistics \code{data.frame} in the
#'   standardized layout \code{variant_id, chrom, pos, A1, A2, z, beta, se, N,
#'   maf} (optional columns omitted when absent on the entry). Combines
#'   tuple-keyed row selection (\code{getSumStats}) with mcols unpacking;
#'   replaces the pre-S4 idiom of pulling \code{S4Vectors::mcols(entry)$<col>}
#'   directly inside pipelines.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selectors (\code{study} for \code{GwasSumStats};
#'   \code{study}, \code{context}, \code{trait} for \code{QtlSumStats}) plus
#'   pass-throughs \code{require}, \code{derive}, \code{keepChrPrefix} forwarded
#'   to the underlying unpacker.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param require Character vector. Columns that must be present (derived if
#'   necessary) in the returned summary-statistics data frame.
#' @param derive Logical. Whether to derive missing standard columns (e.g. Z,
#'   BETA, SE) from the available ones.
#' @param keepChrPrefix Logical. If \code{TRUE}, keep the \code{chr} prefix on
#'   chromosome names; otherwise strip it.
#' @return A \code{data.frame}.
#' @examples
#' data(qtlSumStatsExample)
#' getSumstatDf(qtlSumStatsExample)
#' @export
setGeneric("getSumstatDf", function(x, ...) standardGeneric("getSumstatDf"))

#' @title Get the Embedded QtlDataset List
#' @description Return the named list of \code{QtlDataset} objects carried by a
#'   \code{MultiStudyQtlDataset}.
#' @param x A \code{MultiStudyQtlDataset} object.
#' @return A named list of \code{QtlDataset} objects.
#' @examples
#' data(multiStudyQtlDatasetExample)
#' getQtlDatasets(multiStudyQtlDatasetExample)
#' @export
setGeneric("getQtlDatasets", function(x) standardGeneric("getQtlDatasets"))

#' @title Get the Genome Build
#' @description Return the genome build that the collection's LD sketch and
#'   every entry are aligned to. Because all entries in a \code{GwasSumStats} or
#'   \code{QtlSumStats} share the LD sketch, the genome build is a single value
#'   at the collection level.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused (present for method-signature compatibility).
#' @return Character (length 1).
#' @export
setGeneric("getGenome", function(x, ...) standardGeneric("getGenome"))

#' @title Get QC Audit Record
#' @description Return the audit record of QC steps applied to this collection.
#'   An empty \code{list()} (default on construction) means
#'   \code{\link{summaryStatsQc}} has not yet been run. Pipelines that require
#'   harmonized sumstats (\code{fineMappingPipeline},
#'   \code{twasWeightsPipeline}, and downstream consumers) reject inputs where
#'   \code{length(getQcInfo(x)) == 0L}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused.
#' @return A \code{list} (possibly empty).
#' @export
setGeneric("getQcInfo", function(x, ...) standardGeneric("getQcInfo"))

#' @title Get SLALOM / DENTIST Diagnostics
#' @description Return the per-variant LD-mismatch diagnostics table (SLALOM or
#'   DENTIST output) attached to a sumstats entry by
#'   \code{\link{summaryStatsQc}} when \code{zMismatchQc} was not \code{"none"}.
#'   Convenience accessor over
#'   \code{getQcInfo(x)$entryAudit[[entry]]$ldMismatchDiagnostics}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param entry Integer index (default 1) of the entry whose diagnostics table
#'   to return. When \code{NULL}, returns a named list of every entry's
#'   diagnostics keyed by entry index.
#' @param ... Unused.
#' @return A \code{data.frame} of per-variant diagnostics (columns include
#'   \code{variant_id} plus the SLALOM / DENTIST output columns), \code{NULL}
#'   when no diagnostics were preserved for that entry, or a named list of such
#'   data.frames when \code{entry = NULL}.
#' @export
setGeneric("getQcDiagnostics", function(x, entry = 1L, ...) {
    standardGeneric("getQcDiagnostics")
})

#' @title Get LD Sketch
#' @description Return the genotype panel carrying the LD reference for
#'   this collection. Defined on classes that embed an \code{ldSketch} slot:
#'   \code{GwasSumStats}, \code{QtlSumStats}, \code{FineMappingResult},
#'   \code{TwasWeights}. Returns \code{NULL} when the slot is unset (e.g. a
#'   \code{TwasWeights} fit from individual-level data via \code{QtlDataset}).
#' @param x An S4 object that carries an \code{ldSketch} slot.
#' @param ... Unused.
#' @return A \code{RangedSummarizedExperiment} genotype panel, or
#'   \code{NULL}.
#' @export
setGeneric("getLdSketch", function(x, ...) standardGeneric("getLdSketch"))

# =============================================================================
# LdData accessor generics
# =============================================================================

#' @title Get LD Correlation Matrix
#' @description Extract the LD correlation matrix from an \code{LdData} object.
#'   If only a genotype handle is available, recomputes R from genotypes on the
#'   fly.
#' @param x An \code{LdData} object.
#' @return A correlation matrix, or a list of per-block matrices.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getCorrelation(ld)
#' @export
setGeneric("getCorrelation", function(x) standardGeneric("getCorrelation"))

#' @title Get Genotype Matrix
#' @description Extract a genotype matrix from an object that carries genotype
#'   data. For an \code{LdData}, returns the underlying genotype matrix via its
#'   handle (or \code{NULL} if no handle is available). For a \code{QtlDataset},
#'   returns the genotype matrix for a selected set of traits or region (see
#'   method documentation for the per-class selection arguments).
#' @param x The object to extract from.
#' @param ... Class-specific selection arguments (e.g., \code{traitId},
#'   \code{region}, \code{cisWindow} for \code{QtlDataset}).
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param cisWindow Integer or \code{NULL}. Half-width (bp) of the cis window
#'   around the trait; \code{NULL} uses the dataset default.
#' @param samples Character vector or \code{NULL}. Restrict to these sample IDs;
#'   \code{NULL} uses all samples.
#' @return A numeric matrix, a list of matrices, or \code{NULL}.
#' @examples
#' data(qtlDatasetExample)
#' getGenotypes(qtlDatasetExample)
#' @export
setGeneric("getGenotypes", function(x, ...) standardGeneric("getGenotypes"))

#' @title Check Genotype Availability
#' @description Check whether an \code{LdData} object has a genotype handle for
#'   extracting raw genotypes.
#' @param x An \code{LdData} object.
#' @return Logical.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' hasGenotypes(ld)
#' @export
setGeneric("hasGenotypes", function(x) standardGeneric("hasGenotypes"))

#' @title Get Variant IDs
#' @description Extract variant ID vector from an object that carries one (e.g.,
#'   \code{LdData}, \code{FineMappingRow}, \code{TwasWeightsRow}) or from
#'   one entry of a collection class selected by its identity tuple.
#' @param x The object.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @return Character vector of variant IDs.
#' @examples
#' data(qtlFineMappingExample)
#' getVariantIds(qtlFineMappingExample)
#' @export
setGeneric("getVariantIds", function(x, ...) standardGeneric("getVariantIds"))

#' @title Get Phenotype List
#' @description Extract phenotype data from an object that carries it. For a
#'   \code{QtlDataset}, the user can optionally select specific contexts,
#'   traits, or a region (see method documentation for the per-class selection
#'   arguments).
#' @param x The object to extract from.
#' @param ... Class-specific selection arguments (e.g., \code{contexts},
#'   \code{traitId}, \code{region}).
#' @param contexts Character vector. Context(s) whose data to extract.
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param naAction Character. How to handle missing phenotype values: one of
#'   \code{"keep"}, \code{"drop"}, or \code{"impute"}.
#' @param outlierAction Character. How to handle phenotype outliers: one of
#'   \code{"keep"} or \code{"drop"}.
#' @param outlierPvalThreshold Numeric. Two-sided p-value threshold for flagging
#'   phenotype outliers. Default \code{1e-3}.
#' @return A named list of phenotype matrices or \code{SummarizedExperiment}
#'   objects.
#' @examples
#' data(qtlDatasetExample)
#' getPhenotypes(qtlDatasetExample, contexts = "brain")
#' @export
setGeneric("getPhenotypes", function(x, ...) standardGeneric("getPhenotypes"))
# =============================================================================
# FineMappingResult accessor generics
# =============================================================================

#' @title Get a Single Fine-Mapping Entry
#' @description Return the \code{FineMappingRow} for one \code{(study,
#'   context, trait, method)} row of a \code{FineMappingResult} collection.
#' @param x A \code{FineMappingResult} object.
#' @param study,context,trait,method Single character identifiers. All required
#'   when the collection has more than one row; optional when the collection has
#'   a single row.
#' @return A \code{FineMappingRow} object.
#' @examples
#' data(qtlFineMappingExample)
#' getFineMappingResult(qtlFineMappingExample, study = "study_1",
#'   context = "context_1", trait = "gene_1", method = "susie")
#' @export
setGeneric(
    "getFineMappingResult",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        standardGeneric("getFineMappingResult")
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
#'     \item Colocalization between a GWAS \code{FineMappingResult} and a
#'           QTL \code{FineMappingResult} computed on different variant
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
#' @description Extract posterior inclusion probabilities from a single
#'   \code{FineMappingRow} or from one entry of a \code{FineMappingResult}
#'   (selected by its identity tuple).
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param returnList Logical. If \code{TRUE}, return a per-entry list keyed by
#'   identity tuple instead of a single flattened vector.
#' @return A named numeric vector of PIPs.
#' @examples
#' data(qtlFineMappingExample)
#' getPip(qtlFineMappingExample)
#' @export
setGeneric("getPip", function(x, ...) standardGeneric("getPip"))

#' @title Get SuSiE Fit
#' @description Extract the SuSiE fit object from a fine-mapping entry or
#'   result. The fit may be the trimmed view (when the pipeline ran with the
#'   default \code{trim = TRUE}) or the full untrimmed \code{susie()} return
#'   (when \code{trim = FALSE}).
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @return A list (the SuSiE fit object).
#' @examples
#' data(qtlFineMappingExample)
#' getSusieFit(qtlFineMappingExample)
#' @export
setGeneric("getSusieFit", function(x, ...) standardGeneric("getSusieFit"))

#' @title Compute Between-Credible-Set Correlation On Demand
#' @description Compute the between-credible-set correlation matrix for a
#'   \code{\link{fineMappingRow}} on demand from an LD source. The correlation
#'   is a view over the fit-time LD, which lives on the data object --- a
#'   \code{QtlDataset}'s genotypes or a \code{QtlSumStats} / \code{GwasSumStats}
#'   LD sketch --- and is deliberately never stored on the fit. An
#'   \code{ldSource} is therefore always required.
#' @param x A \code{\link{fineMappingRow}} carrying the fit (credible-set
#'   membership and PIP) and its variant ids.
#' @param ldSource The object supplying the LD: a \code{QtlDataset}
#'   (individual-level; genotypes via \code{\link{getGenotypes}}) or a
#'   \code{QtlSumStats} / \code{GwasSumStats} (summary statistics; LD via
#'   \code{\link{getLdSketch}}). Any other type is an error.
#' @param ... Passed to methods.
#' @return A numeric \code{m x m} between-credible-set correlation matrix
#'   (\code{m} = number of credible sets), or \code{NULL} when the fit has
#'   fewer than two credible sets.
#' @examples
#' data(qtlFineMappingPairedExample)
#' data(qtlSumStatsExample)
#' fe <- getFineMappingResult(qtlFineMappingPairedExample)
#' computeCsCorrelation(fe, qtlSumStatsExample)
#' @export
setGeneric(
    "computeCsCorrelation",
    function(x, ldSource, ...) standardGeneric("computeCsCorrelation")
)

#' @title Get Cross-Validation Result
#' @description Extract the cross-validation payload stored on a
#'   \code{FineMappingRow} (or the matching entry of a
#'   \code{FineMappingResult}). The payload is a list with components
#'   \code{samplePartition} (a \code{data.frame} of \code{Sample}/\code{Fold}
#'   assignments), \code{predictions} (a named list of per-method out-of-fold
#'   prediction matrices), and \code{performance} (a named list of per-method
#'   metric matrices). \code{NULL} when fine-mapping was run without
#'   cross-validation (\code{cvFolds <= 1}).
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A list (the CV payload) or \code{NULL}.
#' @examples
#' data(qtlFineMappingExample)
#' getCvResult(qtlFineMappingExample)
#' @export
setGeneric("getCvResult", function(x, ...) standardGeneric("getCvResult"))

#' @title Get Marginal Effects
#' @description Extract per-variant marginal univariate effects from a
#'   fine-mapping entry or result. Returns a \code{data.frame} with identity
#'   columns (\code{variant_id, chrom, pos, A1, A2}), context (\code{N, MAF}),
#'   and the marginal effect columns (\code{beta, se, z, p}). Populated
#'   uniformly across the individual-level and RSS paths.
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param maxPval Optional numeric (length 1). When non-\code{NULL}, filter rows
#'   where \code{p > maxPval}. Default \code{NULL} (no filter).
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @return A \code{data.frame}.
#' @examples
#' data(qtlFineMappingExample)
#' getMarginalEffects(qtlFineMappingExample)
#' @export
setGeneric("getMarginalEffects", function(x, maxPval = NULL, ...) {
    standardGeneric("getMarginalEffects")
})

#' @title Get Top Loci (posterior view)
#' @description Extract the per-variant posterior fine-mapping payload as either
#'   a \code{data.frame} (default) or a \code{GRanges}. Returns identity columns
#'   (\code{variant_id, chrom, pos, A1, A2}), context (\code{N, MAF}), the
#'   posterior effect columns (\code{beta = posterior_mean, se = posterior_sd}),
#'   \code{pip}, and credible-set membership columns (\code{cs_95}, etc.). Rows
#'   are filtered by PIP by default -- set \code{signalCutoff = 0} to return
#'   every variant.
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param type One of \code{"data.frame"} (default) or \code{"GRanges"}.
#' @param signalCutoff Numeric (length 1). Drop rows where \code{pip <=
#'   signalCutoff}. Default \code{0.025}. Use \code{signalCutoff = 0} to keep
#'   every variant.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param minPurity Numeric or \code{NULL}. Minimum credible-set purity to
#'   retain; \code{NULL} applies no purity filter.
#' @return A \code{data.frame} or a \code{GRanges}.
#' @examples
#' data(qtlFineMappingExample)
#' getTopLoci(qtlFineMappingExample)
#' @export
setGeneric(
    "getTopLoci",
    function(x, type = c("data.frame", "GRanges"), signalCutoff = 0.025, ...) {
        standardGeneric("getTopLoci")
    }
)

#' @title Get Credible Sets
#' @description Extract credible set assignments at the requested coverage.
#' @param x A \code{FineMappingRow} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments plus \code{coverage}.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param coverage Numeric in (0, 1). Credible-set coverage level. Default
#'   \code{0.95}.
#' @param minPurity Numeric or \code{NULL}. Minimum credible-set purity to
#'   retain; \code{NULL} applies no purity filter.
#' @return A data.frame of credible set information.
#' @examples
#' data(qtlFineMappingExample)
#' getCs(qtlFineMappingExample)
#' @export
setGeneric("getCs", function(x, ...) standardGeneric("getCs"))

# =============================================================================
# TwasWeights accessor generics
# =============================================================================

#' @title Get TWAS Weights
#' @description Extract weights from a \code{TwasWeightsRow} or from one entry
#'   of a \code{TwasWeights} collection.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A numeric vector or matrix of weights.
#' @examples
#' twe <- twasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' getWeights(twe)
#' @export
setGeneric("getWeights", function(x, ...) standardGeneric("getWeights"))

#' @title Resolve Per-Variant Weights From a Weight Source
#' @description Return an aligned \code{(variantIds, weights)} pair from a
#'   single weight-source entry, so cTWAS and
#'   \code{\link{causalInferencePipeline}} extract weights identically whether
#'   the source is a \code{TwasWeights} collection (its learned weight vector)
#'   or a \code{FineMappingResult} one (its topLoci posterior effect). The
#'   selectors pin the single row the weights come from.
#' @param x A \code{TwasWeights} or a \code{FineMappingResult}.
#' @param study,context,trait,method Optional length-1 selectors pinning one
#'   row of a collection, as elsewhere; each \code{NULL} (default) leaves that
#'   part of the tuple unconstrained. Ignored when \code{x} is already a row.
#' @param ... Reserved for future use.
#' @return A \code{list} with \code{variantIds} (character) and \code{weights}
#'   (numeric) of equal length; both empty when no usable weights are present.
#' @examples
#' data(twasWeightsExample)
#' w <- resolveWeights(twasWeightsExample, study = "protocol_example",
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
#' @return A \code{data.frame} of posteriors (or \code{NULL} when absent).
#' @examples
#' cre <- CtwasResultEntry(
#'   finemap = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' getFinemap(cre)
#' @export
setGeneric("getFinemap", function(x, ...) standardGeneric("getFinemap"))

#' @title Get cTWAS Per-effect Susie Alpha Table
#' @description Return the full per-effect susie alpha table
#'   (\code{ctwas::finemap_regions} \code{susie_alpha_res} shape) from a
#'   \code{\link{CtwasResultEntry}} or a \code{\link{CtwasResult}} collection
#'   (aggregated across rows, tagged with run identity). This is the fuller
#'   cTWAS output retained so the raw run is reconstructable.
#' @param x A \code{CtwasResultEntry} or \code{CtwasResult}.
#' @param ... Class-specific selection arguments.
#' @return A \code{data.frame} of per-effect alphas (or \code{NULL} when
#'   absent).
#' @examples
#' cre <- CtwasResultEntry(
#'   finemap = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' getSusieAlpha(cre)
#' @export
setGeneric("getSusieAlpha", function(x, ...) standardGeneric("getSusieAlpha"))

#' @title Get cTWAS Group Prior Parameters
#' @description Return the estimated \code{group_prior} / \code{group_prior_var}
#'   of a \code{\link{CtwasResultEntry}}.
#' @param x A \code{CtwasResultEntry}.
#' @param ... Reserved for future use.
#' @return The stored parameter object (typically a list), or \code{NULL}.
#' @examples
#' cre <- CtwasResultEntry(
#'   finemap = data.frame(id = c("g1", "g2"), susie_pip = c(0.9, 0.1)),
#'   susieAlpha = data.frame(id = c("g1", "g2"), alpha = c(0.9, 0.1)))
#' getCtwasParam(cre)
#' @export
setGeneric("getCtwasParam", function(x, ...) standardGeneric("getCtwasParam"))

#' @title Get Standardized Flag
#' @description Check whether weights are on the standardized scale.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return Logical.
#' @examples
#' data(ctwasWeightsExample)
#' getStandardized(ctwasWeightsExample)
#' @export
setGeneric("getStandardized", function(x, ...) {
    standardGeneric("getStandardized")
})

#' @title Get Model Fits
#' @description Extract fitted model objects.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return Method-specific (typically a list).
#' @examples
#' data(ctwasWeightsExample)
#' # NULL unless the pipeline was run with `retainFit = TRUE`
#' getFits(ctwasWeightsExample)
#' @export
setGeneric("getFits", function(x, ...) standardGeneric("getFits"))

#' @title Get the Full-Data Prior from a MashPrior
#' @description Accessor for the \code{fullFit} slot (the full-data data-driven
#'   prior payload).
#' @param x A \code{MashPrior} object.
#' @param ... Unused.
#' @return The full-data prior payload, or \code{NULL}.
#' @examples
#' U <- list(shared = diag(3), singleton = matrix(0.3, 3, 3) + diag(0.7, 3))
#' mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))
#' getFullFit(mp)
#' @export
setGeneric("getFullFit", function(x, ...) standardGeneric("getFullFit"))

#' @title Get the Per-Fold Priors from a MashPrior
#' @description Accessor for the \code{cvFits} slot (per-fold priors +
#'   \code{samplePartition}).
#' @param x A \code{MashPrior} object.
#' @param ... Unused.
#' @return The \code{cvFits} list, or \code{NULL}.
#' @examples
#' U <- list(shared = diag(3), singleton = matrix(0.3, 3, 3) + diag(0.7, 3))
#' mp <- MashPrior(fullFit = list(U = U, w = c(0.5, 0.5)))
#' getCvFits(mp)
#' @export
setGeneric("getCvFits", function(x, ...) standardGeneric("getCvFits"))

#' @title Get Method Names
#' @description Extract method names from a collection class.
#' @param x A \code{FineMappingResult} or \code{TwasWeights} object.
#' @return Character vector.
#' @export
setGeneric("getMethodNames", function(x) standardGeneric("getMethodNames"))

#' @title Get Data Type
#' @description Extract the data-type tag.
#' @param x A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param study Character (length 1) or \code{NULL}. Restrict the selection to
#'   this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A character vector or NULL.
#' @examples
#' twe <- twasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' getDataType(twe)
#' @export
setGeneric("getDataType", function(x, ...) standardGeneric("getDataType"))

# =============================================================================
# AlleleQcResult accessor generics
# QcResult accessor generics
# =============================================================================
# VCF/BCF writer generic
# =============================================================================

#' Write summary statistics or fine-mapping results to VCF/BCF
#'
#' Creates a VCF object from GWAS summary statistics or fine-mapping results and
#' writes it to disk. Supports bgzipped VCF (.vcf.gz/.vcf.bgz) and BCF (.bcf)
#' output formats via VariantAnnotation and Rsamtools.
#'
#' @param x Input data: a \code{GwasSumStats} object, a \code{FineMappingResult}
#'   object, or a data.frame with columns \code{chrom}, \code{pos}, \code{ref},
#'   \code{alt}.
#' @param outputPath File path for output. Extension determines format:
#'   \code{.vcf.gz} or \code{.vcf.bgz} for bgzipped VCF, \code{.bcf} for BCF,
#'   \code{.vcf} for uncompressed VCF.
#' @param sampleName Name for the VCF sample column (default: trait name or
#'   method name from the S4 object).
#' @param study Character or \code{NULL}. Restrict the written records to this
#'   study; \code{NULL} includes all studies.
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
#' @param ... Additional arguments passed to methods.
#' @return Invisible path to the written file.
#' @examples
#' data(gwasSumStatsS4Example)
#' writeSumstatsVcf(
#'   gwasSumStatsS4Example, outputPath = tempfile(fileext = ".vcf"))
#' @export
setGeneric("writeSumstatsVcf", function(x, outputPath, sampleName = NULL, ...) {
    standardGeneric("writeSumstatsVcf")
})

# =============================================================================
# QtlDataset accessor generics
# =============================================================================

#' @title Get Study Identifier
#' @description Return the study identifier carried by a \code{QtlDataset}.
#' @param x A \code{QtlDataset} object.
#' @return Character (length 1).
#' @export
setGeneric("getStudy", function(x) standardGeneric("getStudy"))

#' @title Get Context Names
#' @description Return the names of all contexts carried by an object (e.g., the
#'   keys of the \code{phenotypes} list on a \code{QtlDataset}, or the unique
#'   \code{context} values of a \code{QtlSumStats}).
#' @param x The object.
#' @return Character vector of context names.
#' @examples
#' data(qtlDatasetExample)
#' getContexts(qtlDatasetExample)
#' @export
setGeneric("getContexts", function(x) standardGeneric("getContexts"))

#' @title Get Unique Trait Names
#' @description Return the unique trait identifiers carried by a collection
#'   class (e.g., \code{QtlSumStats}).
#' @param x The object.
#' @return Character vector of unique trait names.
#' @examples
#' data(qtlSumStatsExample)
#' getTraits(qtlSumStatsExample)
#' @export
setGeneric("getTraits", function(x) standardGeneric("getTraits"))

#' @title Get Per-Row Genomic Regions
#' @description Return the genomic anchor of each row of a per-tuple collection
#'   as a \code{GRanges} with one range per row -- the trait's own region for a
#'   \code{TwasWeights}, or the fine-mapping window for a
#'   \code{FineMappingResult}. This is location provenance (e.g. for cTWAS
#'   LD-block placement).
#' @param x The object.
#' @param ... Reserved for future use.
#' @return A \code{GRanges} with one range per row of \code{x}, or an empty
#'   \code{GRanges} when the collection carries no region provenance.
#' @export
setGeneric("getRegion", function(x, ...) standardGeneric("getRegion"))

#' @title Get Per-Row Trait Positions
#' @description Return the molecular feature's OWN genomic coordinates (the gene
#'   / peak range; TSS = \code{start()}) for each row of a per-tuple collection,
#'   as a \code{GRanges} with one range per row. Distinct from
#'   \code{\link{getRegion}} (the fine-mapping window for a
#'   \code{FineMappingResult}): the true trait position cannot be inferred from
#'   summary statistics, so it is threaded from the QtlDataset \code{rowRanges}
#'   or an explicit QtlSumStats trait-position and carried as provenance onto
#'   \code{QtlFineMappingResult} and \code{TwasWeights}.
#' @param x The object.
#' @param ... Reserved for future use.
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @return A \code{GRanges} with one range per row of \code{x}, or an empty
#'   \code{GRanges} when the collection carries no trait-position provenance.
#' @export
setGeneric("getTraitPosition", function(x, ...) {
    standardGeneric("getTraitPosition")
})

#' @title Get Residualized Genotypes
#' @description Residualize the genotype matrix against the per-context
#'   phenotype covariates and the genotype covariates, optionally subsetting
#'   variants to those falling within a trait's cis-window or an explicit
#'   region.
#' @param x A \code{QtlDataset} object.
#' @param ... Selection arguments: \code{traitId}, \code{region},
#'   \code{cisWindow}, \code{phenotypeCovariatesToRemove},
#'   \code{genotypeCovariatesToRemove}, and \code{covariateNaAction}
#'   (\code{"impute"}, the default, mean-imputes missing covariate cells;
#'   \code{"drop"} removes samples with any missing covariate).
#' @param contexts Character vector. Context(s) whose data to extract.
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param cisWindow Integer or \code{NULL}. Half-width (bp) of the cis window
#'   around the trait; \code{NULL} uses the dataset default.
#' @param samples Character vector or \code{NULL}. Restrict to these sample IDs;
#'   \code{NULL} uses all samples.
#' @param phenotypeCovariatesToResidualize Character vector or \code{NULL}.
#'   Phenotype covariate names to residualize; \code{NULL} uses the dataset
#'   default set.
#' @param genotypeCovariatesToResidualize Character vector or \code{NULL}.
#'   Genotype covariate names to residualize; \code{NULL} uses the dataset
#'   default set.
#' @param residualizePhenotypeCovariates Logical. Whether to residualize the
#'   phenotype covariates. Default \code{TRUE}.
#' @param residualizeGenotypeCovariates Logical. Whether to residualize the
#'   genotype covariates. Default \code{TRUE}.
#' @param residualizePhenotypeCovariatesFromGenotypes Logical or \code{NULL}.
#'   Whether to residualize phenotype covariates out of the genotype design;
#'   \code{NULL} uses the default.
#' @param residualizeGenotypeCovariatesFromGenotypes Logical or \code{NULL}.
#'   Whether to residualize genotype covariates out of the genotype design;
#'   \code{NULL} uses the default.
#' @param covariateNaAction Character. How to handle missing covariate values.
#' @return A numeric matrix (samples x variants).
#' @examples
#' data(qtlDatasetExample)
#' getResidualizedGenotypes(qtlDatasetExample, contexts = "brain")
#' @export
setGeneric("getResidualizedGenotypes", function(x, ...) {
    standardGeneric("getResidualizedGenotypes")
})

#' @title Get Residualized Phenotypes
#' @description Residualize the per-context phenotype matrices against the
#'   per-context phenotype covariates and the genotype covariates, for one or
#'   more requested contexts.
#' @param x A \code{QtlDataset} object.
#' @param ... Selection arguments: \code{contexts} (required), \code{traitId},
#'   \code{region}, \code{phenotypeCovariatesToRemove},
#'   \code{genotypeCovariatesToRemove}, and \code{covariateNaAction}
#'   (\code{"impute"}, the default, mean-imputes missing covariate cells;
#'   \code{"drop"} removes samples with any missing covariate).
#' @param contexts Character vector. Context(s) whose data to extract.
#' @param traitId Character or \code{NULL}. Molecular trait / feature identifier
#'   to select; \code{NULL} matches all traits in the dataset.
#' @param region Character (length 1, \code{"chr:start-end"}) or \code{NULL}.
#'   Restrict variants to this region; \code{NULL} uses the full cis window /
#'   all regions.
#' @param phenotypeCovariatesToResidualize Character vector or \code{NULL}.
#'   Phenotype covariate names to residualize; \code{NULL} uses the dataset
#'   default set.
#' @param genotypeCovariatesToResidualize Character vector or \code{NULL}.
#'   Genotype covariate names to residualize; \code{NULL} uses the dataset
#'   default set.
#' @param residualizePhenotypeCovariates Logical. Whether to residualize the
#'   phenotype covariates. Default \code{TRUE}.
#' @param residualizeGenotypeCovariates Logical. Whether to residualize the
#'   genotype covariates. Default \code{TRUE}.
#' @param residualizePhenotypeCovariatesFromPhenotypes Logical or \code{NULL}.
#'   Whether to residualize phenotype covariates out of the phenotype;
#'   \code{NULL} uses the default.
#' @param residualizeGenotypeCovariatesFromPhenotypes Logical or \code{NULL}.
#'   Whether to residualize genotype covariates out of the phenotype;
#'   \code{NULL} uses the default.
#' @param naAction Character. How to handle missing phenotype values: one of
#'   \code{"keep"}, \code{"drop"}, or \code{"impute"}.
#' @param covariateNaAction Character. How to handle missing covariate values.
#' @param outlierAction Character. How to handle phenotype outliers: one of
#'   \code{"keep"} or \code{"drop"}.
#' @param outlierPvalThreshold Numeric. Two-sided p-value threshold for flagging
#'   phenotype outliers. Default \code{1e-3}.
#' @return A named list of numeric matrices keyed by context.
#' @examples
#' data(qtlDatasetExample)
#' getResidualizedPhenotypes(qtlDatasetExample, contexts = "brain")
#' @export
setGeneric("getResidualizedPhenotypes", function(x, ...) {
    standardGeneric("getResidualizedPhenotypes")
})

#' @title Get Per-Context Phenotype Covariates
#' @description Return per-context phenotype covariate matrices, taken from the
#'   \code{colData} of each context's \code{SummarizedExperiment}.
#' @param x A \code{QtlDataset} object.
#' @param contexts Character vector of context names (subset of
#'   \code{names(getPhenotypes(x))}).
#' @return A named list of matrices keyed by context.
#' @examples
#' data(qtlDatasetExample)
#' getPhenotypeCovariates(qtlDatasetExample, contexts = "brain")
#' @export
setGeneric("getPhenotypeCovariates", function(x, contexts) {
    standardGeneric("getPhenotypeCovariates")
})

#' @title Get Genotype Covariates
#' @description Return the single genotype-derived covariate matrix carried by a
#'   \code{QtlDataset} (e.g., ancestry PCs).
#' @param x A \code{QtlDataset} object.
#' @return Numeric matrix (samples x covariates).
#' @examples
#' data(qtlDatasetExample)
#' getGenotypeCovariates(qtlDatasetExample)
#' @export
setGeneric("getGenotypeCovariates", function(x) {
    standardGeneric("getGenotypeCovariates")
})

#' @title Get scaleResiduals Flag
#' @description Whether residualization accessors scale residuals to unit
#'   variance.
#' @param x A \code{QtlDataset} object.
#' @return Logical (length 1).
#' @examples
#' data(qtlDatasetExample)
#' getScaleResiduals(qtlDatasetExample)
#' @export
setGeneric("getScaleResiduals", function(x) {
    standardGeneric("getScaleResiduals")
})

# =============================================================================
# GenotypeHandle / LD-statistic / Annotation / LdData / H2Estimate accessors
# =============================================================================

#' @title Get SNP Info
#' @description Return the cached SNP metadata data.frame (columns: SNP, CHR,
#'   BP, A1, A2, optionally MAF).
#' @param x An object carrying cached SNP metadata.
#' @return A data.frame.
#' @keywords internal
setGeneric("getSnpInfo", function(x) standardGeneric("getSnpInfo"))

#' @title Get Genotype Storage Format
#' @description Return the detected genotype storage format.
#' @param x A \code{GenotypeHandle}.
#' @return Character (length 1): one of "gds", "vcf", "plink1", "plink2".
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getFormat", function(x) standardGeneric("getFormat"))

#' @title Get File Path
#' @description Return the underlying genotype file path or stem.
#' @param x A \code{GenotypeHandle}.
#' @return Character (length 1).
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getPath", function(x) standardGeneric("getPath"))

#' @title Get Sample Identifiers
#' @description Return the sample-id vector.
#' @param x A \code{GenotypeHandle}.
#' @return Character vector.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getSampleIds", function(x) standardGeneric("getSampleIds"))

#' @title Get plink2 pgen Pointer
#' @description Return the cached external pointer to the plink2 pgen handle
#'   (NULL when the handle is not pgen-backed).
#' @param x A \code{GenotypeHandle}.
#' @return An external pointer or NULL.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getPgenPtr", function(x) standardGeneric("getPgenPtr"))

#' @title Get Sample Count
#' @description Return the number of samples carried by a \code{GenotypeHandle}.
#' @param x A \code{GenotypeHandle}.
#' @return Integer (length 1).
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getNSamples", function(x) standardGeneric("getNSamples"))

#' @title Get Per-Block Eigendecompositions
#' @description Return the per-block eigendecomposition list carried by an
#'   \code{LdEigen} object.
#' @param x An \code{LdEigen}.
#' @return List of per-block eigen decompositions.
#' @examples
#' data(ldEigenExample)
#' getEigenList(ldEigenExample)
#' @export
setGeneric("getEigenList", function(x) standardGeneric("getEigenList"))

#' @title Get LD Reference Panel Size
#' @description Return the reference-panel sample size used to compute an
#'   \code{LdStatistic} or carried by an \code{LdData}.
#' @param x An \code{LdStatistic} or \code{LdData}.
#' @return Integer (length 1).
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getNRef(ld)
#' @export
setGeneric("getNRef", function(x) standardGeneric("getNRef"))

#' @title Get In-Sample Flag
#' @description Whether the LD reference panel is from the same cohort as the
#'   GWAS (affects bias correction).
#' @param x An \code{LdStatistic}.
#' @return Logical (length 1).
#' @examples
#' data(ldEigenExample)
#' getInSample(ldEigenExample)
#' @export
setGeneric("getInSample", function(x) standardGeneric("getInSample"))

#' @title Get LD Scores
#' @description Return the per-SNP LD score matrix carried by an \code{LdScore}
#'   object.
#' @param x An \code{LdScore}.
#' @return Numeric matrix (SNPs x annotations+1).
#' @examples
#' data(ldScoreExample)
#' getLdScores(ldScoreExample)
#' @export
setGeneric("getLdScores", function(x) standardGeneric("getLdScores"))

#' @title Get LD-Score Regression Weights
#' @description Return the per-SNP regression weights vector carried by an
#'   \code{LdScore} object.
#' @param x An \code{LdScore}.
#' @return Numeric vector.
#' @examples
#' data(ldScoreExample)
#' getLdScoreWeights(ldScoreExample)
#' @export
setGeneric("getLdScoreWeights", function(x) {
    standardGeneric("getLdScoreWeights")
})

#' @title Get Per-Block LD Matrix List
#' @description Return the list of per-block LD (R^2) matrices used for the FGLS
#'   residual covariance in g-LDSC.
#' @param x An \code{LdScore}.
#' @return List of matrices (empty list for S-LDSC).
#' @examples
#' data(ldScoreExample)
#' getLdMatrixList(ldScoreExample)
#' @export
setGeneric("getLdMatrixList", function(x) standardGeneric("getLdMatrixList"))

#' @title Get LD Block Container
#' @description Return the \code{GRanges} of LD blocks carried by an
#'   \code{LdStatistic}.
#' @param x An \code{LdStatistic}.
#' @return A \code{GRanges} of LD block intervals.
#' @examples
#' data(ldEigenExample)
#' getLdBlocks(ldEigenExample)
#' @export
setGeneric("getLdBlocks", function(x) standardGeneric("getLdBlocks"))


#' @title Get GenotypeHandle from LdData
#' @description Return the \code{GenotypeHandle} (or list of handles for mixture
#'   panels) carried by an \code{LdData}.
#' @param x An \code{LdData}.
#' @return A \code{GenotypeHandle}, a list of them, or NULL. For a
#'   \code{QtlDataset} this is the handle itself; \code{getGenotypes()}
#'   extracts a dosage block from it instead.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getGenotypeHandle", function(x) {
    standardGeneric("getGenotypeHandle")
})

#' @title Get Mixture Weights
#' @description Return the per-panel mixing proportions carried by an
#'   \code{LdData} when its \code{genotypeHandle} slot is a list of panels. NULL
#'   for single-panel objects.
#' @param x An \code{LdData}.
#' @return Numeric vector or NULL.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getMixtureWeights(ld)
#' @export
setGeneric("getMixtureWeights", function(x) {
    standardGeneric("getMixtureWeights")
})

#' @title Get SNP Indices
#' @description Return the integer indices into the handle's snpInfo carried by
#'   an \code{LdData}.
#' @param x An \code{LdData}.
#' @return Integer vector or NULL.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getSnpIdx(ld)
#' @export
setGeneric("getSnpIdx", function(x) standardGeneric("getSnpIdx"))

#' @title Get Variant GRanges
#' @description Return the variant metadata \code{GRanges} of an \code{LdData}.
#' @param x An \code{LdData}.
#' @return A \code{GRanges}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getVariantInfo(ld)
#' @export
setGeneric("getVariantInfo", function(x) standardGeneric("getVariantInfo"))

#' @title Get Block Metadata
#' @description Return the block metadata (\code{GRanges} or \code{data.frame})
#'   carried by an \code{LdData}.
#' @param x An \code{LdData}.
#' @return A \code{GRanges} or \code{data.frame}.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getBlockMetadata(ld)
#' @export
setGeneric("getBlockMetadata", function(x) standardGeneric("getBlockMetadata"))

#' @title Get Reference Panel (data.frame)
#' @description Flatten the variant \code{GRanges} of an \code{LdData} into a
#'   reference-panel data.frame.
#' @param x An \code{LdData}.
#' @return A data.frame.
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:8]
#' gr <- GenomicRanges::GRanges("22",
#'   IRanges::IRanges(seq(1L, by = 100L, length.out = 8), width = 1L))
#' ld <- LdData(correlation = cor(X), variants = gr,
#'   blockMetadata = S4Vectors::DataFrame(
#'     chrom = "22", start = 1L, end = 1000L))
#' getRefPanel(ld)
#' @export
setGeneric("getRefPanel", function(x) standardGeneric("getRefPanel"))

#' @title Get Per-Block tau Matrix
#' @description Return the per-block jackknife tau matrix carried by an
#'   \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return A numeric matrix or NULL.
#' @examples
#' data(h2EstimateExample)
#' getTauBlocks(h2EstimateExample)
#' @export
setGeneric("getTauBlocks", function(x) standardGeneric("getTauBlocks"))

#' @title Get Global SNP Heritability
#' @description Return the global SNP heritability estimate carried by an
#'   \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @examples
#' data(h2EstimateExample)
#' getH2(h2EstimateExample)
#' @export
setGeneric("getH2", function(x) standardGeneric("getH2"))

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
#' @return A \code{data.frame} of annotations (CHR, SNP, annotation columns).
#' @rdname getAnnotData
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getAnnotData(sd)
#' @export
setGeneric("getAnnotData", function(x) standardGeneric("getAnnotData"))

#' @title Get the allele-frequency table from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A \code{data.frame} of reference-panel frequencies (SNP, MAF).
#' @rdname getFrqData
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getFrqData(sd)
#' @export
setGeneric("getFrqData", function(x) standardGeneric("getFrqData"))

#' @title Get the per-trait runs list from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return The named list of per-trait \code{single}/\code{joint} runs.
#' @rdname getTraitRuns
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getTraitRuns(sd)
#' @export
setGeneric("getTraitRuns", function(x) standardGeneric("getTraitRuns"))

#' @title Get the trait names from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A character vector of trait names.
#' @rdname getTraitNames
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getTraitNames(sd)
#' @export
setGeneric("getTraitNames", function(x) standardGeneric("getTraitNames"))

#' @title Get the annotation column names from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @return A character vector of annotation column names.
#' @rdname getAnnotCols
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getAnnotCols(sd)
#' @export
setGeneric("getAnnotCols", function(x) standardGeneric("getAnnotCols"))

#' @title Get one trait's run from an SldscData
#' @param x An \code{\link{SldscData}} object.
#' @param trait Character. Trait name.
#' @param ... Further arguments: \code{mode} (\code{"single"}/\code{"joint"})
#'   and \code{idx} (which single run).
#' @param mode Character. Trait-run selection mode.
#' @param idx Integer. Index of the trait run to select.
#' @return A single run list, the list of single runs, or \code{NULL}.
#' @rdname getTraitRun
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
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
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
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' getTraitRun(sd, "traitX")
#' @export
setGeneric("getTraitRun", function(x, trait, ...) {
    standardGeneric("getTraitRun")
})

#' @title Colocalization Views
#' @description The four views over a \code{\linkS4class{ColocResult}}. The
#'   object stores one element per tested (QTL credible set, GWAS credible set,
#'   block) pair; these accessors project it to the granularity a given
#'   analysis needs, rather than any one of them being the stored shape.
#'
#'   \describe{
#'     \item{\code{getColocPairs()}}{One row per tested pair -- coloc's native
#'       \code{$summary} granularity. Every testable pair is reported verbatim,
#'       including both halves of a credible set that straddles a block
#'       boundary; interpreting them is the caller's decision.}
#'     \item{\code{getColocVariants()}}{One row per (pair, variant), carrying
#'       \code{colocPp = PP.H4.abf * SNP.PP.H4} -- the posterior that this
#'       variant is the shared causal one. With \code{pooled = TRUE}, one row
#'       per (gene, variant) pooled by the rule below.}
#'     \item{\code{getColocCredibleSets()}}{One row per coloc credible set: the
#'       smallest set of variants whose cumulative \code{SNP.PP.H4} reaches
#'       \code{coverage}. This is a THIRD variant set, not guaranteed to be a
#'       subset of either input credible set, so its purity is recomputed from
#'       LD rather than inherited.}
#'     \item{\code{getColocGenes()}}{One row per gene, pooled across pairs.}
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
#' @param pooled Logical. \code{getColocVariants()} only: pool to one row per
#'   (gene, variant) instead of per (pair, variant).
#' @param coverage Cumulative \code{SNP.PP.H4} the credible set must reach.
#' @param minPp4 Optional lower bound on the pair's \code{PP.H4.abf}. Applied
#'   before any LD work, so a stricter threshold costs strictly less.
#' @param requireMaxH4 Logical (length 1), default \code{FALSE}. When
#'   \code{TRUE}, keep only pairs whose \code{PP.H4.abf} is the largest of
#'   the five posterior probabilities, i.e. colocalization is the favoured
#'   hypothesis rather than merely a probable one.
#' @param minAbsCorr Minimum absolute correlation between credible-set members
#'   (purity). Requires the object's \code{ldSketch}; when that is absent,
#'   purity is reported as \code{NA} and no set is dropped.
#' @param ... Additional arguments passed on to methods.
#' @return A tibble at the requested granularity.
#' @name colocViews
#' @rdname colocViews
#' @examples
#' data(qtlFineMappingLbfExample)
#' data(gwasFineMappingLbfExample)
#' res <- colocPipeline(
#'     qtlFineMappingResult = qtlFineMappingLbfExample,
#'     gwasInput = gwasFineMappingLbfExample
#' )
#' head(getColocPairs(res))
#' head(getColocVariants(res))
#' getColocGenes(res)
#' head(getColocCredibleSets(res, minPp4 = 0.001))
NULL

#' @rdname colocViews
#' @export
setGeneric("getColocPairs", function(x, ...) {
    standardGeneric("getColocPairs")
})

#' @rdname colocViews
#' @export
setGeneric("getColocVariants", function(x, pooled = FALSE, ...) {
    standardGeneric("getColocVariants")
})

#' @rdname colocViews
#' @export
setGeneric(
    "getColocCredibleSets",
    function(x, coverage = 0.95, minPp4 = NULL, minAbsCorr = 0.8, ...) {
        standardGeneric("getColocCredibleSets")
    }
)

#' @rdname colocViews
#' @export
setGeneric("getColocGenes", function(x, ...) {
    standardGeneric("getColocGenes")
})

#' @title QtlDataset Filter Settings
#' @description The variant- and sample-filter settings a
#'   \code{\link{QtlDataset}} was built with. They are recorded on the object
#'   rather than reapplied per call, so reading them back is how a caller finds
#'   out what a dataset already excludes.
#' @param x A \code{QtlDataset}.
#' @param ... Additional arguments passed on to methods.
#' @return \code{getMafCutoff()}, \code{getMacCutoff()},
#'   \code{getXvarCutoff()} and \code{getImissCutoff()} single numerics;
#'   \code{getKeepVariants()} a character vector (empty when no restriction
#'   was requested); \code{getKeepIndel()} a single logical. The sample set
#'   is not among them: it lives in \code{colData(x)} and \code{sampleMap(x)}.
#' @name qtlDatasetFilters
#' @rdname qtlDatasetFilters
#' @examples
#' data(qtlDatasetExample)
#' getMafCutoff(qtlDatasetExample)
#' getKeepIndel(qtlDatasetExample)
NULL

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getMafCutoff", function(x, ...) standardGeneric("getMafCutoff"))

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getMacCutoff", function(x, ...) standardGeneric("getMacCutoff"))

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getXvarCutoff", function(x, ...) standardGeneric("getXvarCutoff"))

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getImissCutoff", function(x, ...) {
    standardGeneric("getImissCutoff")
})

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getKeepVariants", function(x, ...) {
    standardGeneric("getKeepVariants")
})

#' @rdname qtlDatasetFilters
#' @export
setGeneric("getKeepIndel", function(x, ...) standardGeneric("getKeepIndel"))

#' @title Get Per-Chromosome Payload Paths
#' @description The chromosome-to-path map of a sharded
#'   \code{\link{GenotypeHandle}}. Empty (\code{character(0)}) for a
#'   single-file handle, whose one path is \code{\link{getPath}} instead; a
#'   non-empty map marks a one-file-per-chromosome handle whose extraction is
#'   routed by chromosome.
#' @param x A \code{GenotypeHandle}.
#' @param ... Additional arguments passed on to methods.
#' @return A named character vector mapping canonical chromosome to path or
#'   prefix, possibly empty.
#' @details The genotype handle is the seed layer behind a panel and
#'   is not part of the public interface; obtain a panel from
#'   \code{readGenotypes()} and ask it directly.
#' @keywords internal
setGeneric("getChromPaths", function(x, ...) standardGeneric("getChromPaths"))

#' Per-variant per-effect log Bayes factors (wide)
#'
#' The variant x effect matrix of single-effect log Bayes factors
#' (\code{lbf_variable}, or fSuSiE's \code{lBF}) from the stored fit: one row
#' per variant, one \code{lbf_L<k>} column per effect. The per-variant scalar
#' summary (max across effects) is the \code{logBF} column of
#' \code{\link{getTopLoci}}; this accessor keeps the full per-effect breakdown.
#' Across entries with different effect counts the collection method NA-fills
#' the ragged columns.
#'
#' @param x A \code{FineMappingRow} or \code{FineMappingResultBase}.
#' @param ... Ignored.
#' @return A \code{tibble}: \code{variant_id} + \code{lbf_L1..lbf_LL} (the
#'   collection method also carries the entry identity columns).
#' @seealso \code{\link{getTopLoci}} (the scalar \code{logBF} column)
#' @examples
#' data(qtlFineMappingExample)
#' getLbf(qtlFineMappingExample)
#' @export
setGeneric("getLbf", function(x, ...) standardGeneric("getLbf"))

#' Per-credible-set summary of a fine-mapping result
#'
#' One row per credible set (at a given coverage) with its size, purity, prior
#' variance, log Bayes factor, and lead variant -- the per-CS complement to the
#' per-variant \code{\link{getCs}}. Replaces the legacy per-effect `effect.tsv`.
#'
#' @param x A \code{FineMappingRow} or \code{FineMappingResultBase}.
#' @param coverage Credible-set coverage to summarise. Default 0.95.
#' @param ... Ignored.
#' @return A \code{tibble}: \code{cs, effect_id, coverage, n_variants,
#'   purity_min, purity_mean, V, cs_log10bf} (strongest member logBF),
#'   \code{cs_log_bf} (true per-effect single-effect log Bayes factor),
#'   \code{cs_pip} (summed member PIP = inclusion mass captured),
#'   \code{cs_mean_effect} (mean posterior conditional effect),
#'   \code{lead_variant, lead_pip} (the collection method also carries the entry
#'   identity columns).
#' @seealso \code{\link{getCs}}, \code{\link{getTopLoci}}
#' @examples
#' data(qtlFineMappingExample)
#' getCredibleSetSummary(qtlFineMappingExample)
#' @export
setGeneric("getCredibleSetSummary", function(x, ...) {
    standardGeneric("getCredibleSetSummary")
})

#' fSuSiE credible band (fitted effect + uncertainty band)
#'
#' For a functional-SuSiE (\code{fsusieR::susiF}) fine-mapping entry, returns
#' the fitted effect curve and its credible band over the functional grid, one
#' row per (credible set, grid position). Wraps the upstream fSuSiE band
#' computation (fsusieR's internal \code{update_cal_credible_band.susiF} +
#' \code{get_fitted_effect}). Requires an UNtrimmed fit (the wavelet slots are
#' dropped by trimming); degrades to zero rows for non-fSuSiE or trimmed fits.
#'
#' @param x A \code{FineMappingRow} or \code{FineMappingResultBase}.
#' @param ... Ignored.
#' @return A long \code{tibble}: \code{cs, chrom, pos, effect, lower, upper}
#'   (the collection method additionally carries the entry identity columns).
#' @seealso \code{\link{fsusieAffectedRegions}}
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieCredibleBand(fsusieFineMappingExample)
#' @export
setGeneric("fsusieCredibleBand", function(x, ...) {
    standardGeneric("fsusieCredibleBand")
})

#' fSuSiE affected genomic regions
#'
#' The sub-intervals of the functional grid where a credible set's band excludes
#' zero (the genomic footprint of each fSuSiE effect), from the upstream
#' \code{fsusieR::affected_reg}. Returns a \code{GRanges} with the credible-set
#' label, its purity, and the effect \code{direction}
#' (\code{"pos"}/\code{"neg"}; upstream discards it). Requires an UNtrimmed fit;
#' empty otherwise.
#'
#' @param x A \code{FineMappingRow} or \code{FineMappingResultBase}.
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
