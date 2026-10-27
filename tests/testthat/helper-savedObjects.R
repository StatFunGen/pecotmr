# Objects as an earlier pecotmr saved them, for the updateObject() tests of
# the collections (test_RangedTupleList.R, test_AllClasses.R,
# test_twasWeights.R, test_ColocBoostResult.R), test_LdStatistic.R,
# test_H2Estimate.R, test_QtlDataset.R and the panel tests in
# test_genotypeIo.R.

# `x` with each slot in `renames` (old name = current name) stored under its
# old name, as an object saved before the rename holds it: R stores a slot as
# an attribute, and a NULL one as a placeholder symbol, which `slot<-` writes.
.savedUnder <- function(x, renames) {
    reduce2(names(renames), unname(renames), .savedUnderName, .init = x)
}

.savedUnderName <- function(x, old, new) {
    value <- methods::slot(x, new)
    attr(x, new) <- NULL
    `slot<-`(x, old, check = FALSE, value = value)
}

# A collection as saved before its references became per-row: one reference
# (or NULL) in the slot and no row key.
.savedWithOneLdRef <- function(x, ref) {
    old <- `slot<-`(x, "ldReference", check = FALSE, value = ref)
    md <- attr(old, "elementMetadata")
    md[["ldReference"]] <- NULL
    attr(old, "elementMetadata") <- md
    old
}

# ...and as saved before the slot was renamed, holding that reference as
# `ldSketch`.
.savedWithLdSketch <- function(x, ref) {
    .savedUnder(.savedWithOneLdRef(x, ref), c(ldSketch = "ldReference"))
}

# A genotype panel as saved before panels recorded where their genotypes come
# from: its dosages read through a GenotypeHandle's GhSeed at `path` -- one
# claiming `format`, which may be a format a handle no longer reads -- with
# no metadata and only the SNP / A1 / A2 variant columns, and reading the
# per-chromosome files `chromPaths` when given. Its variants, samples and
# sample data are copied from `panel`, today's panel over the same files.
.savedHandlePanel <- function(panel, path, format, chromPaths = character(0)) {
    rr <- SummarizedExperiment::rowRanges(panel)
    handle <- methods::new(
        "GenotypeHandle",
        path = path,
        format = "plink2",
        snpInfo = data.frame(
            SNP = rownames(panel),
            CHR = as.character(GenomicRanges::seqnames(rr)),
            BP = GenomicRanges::start(rr),
            A1 = rr$A1,
            A2 = rr$A2
        ),
        nSamples = ncol(panel),
        sampleIds = colnames(panel),
        pgenPtr = NULL,
        chromPaths = chromPaths
    )
    old <- `slot<-`(handle, "format", check = FALSE, value = format)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = genotypeDelayedArray(old)),
        rowRanges = S4Vectors::`mcols<-`(
            rr,
            value = S4Vectors::mcols(rr)[c("SNP", "A1", "A2")]
        ),
        colData = SummarizedExperiment::colData(panel)
    )
}
