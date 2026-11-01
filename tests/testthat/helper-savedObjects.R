# Objects as an earlier pecotmr saved them, for the updateObject() tests of
# the collections (test_RangedTupleList.R, test_AllClasses.R,
# test_twasWeights.R, test_crossValidation.R, test_ColocBoostResult.R),
# test_LdStatistic.R, test_H2Estimate.R, test_QtlDataset.R and the panel
# tests in test_genotypeIo.R.

# `x` with each slot in `renames` (old name = current name) stored under its
# old name, as an object saved before the rename holds it: R stores a slot as
# an attribute, and a NULL one as a placeholder symbol, which `slot<-` writes.
.savedUnder <- function(x, renames) {
    reduce2(names(renames), unname(renames), .savedUnderName, .init = x)
}

.savedUnderName <- function(x, old, new) {
    value <- methods::slot(x, new)
    attr(x, new) <- NULL
    methods::slot(x, old, check = FALSE) <- value
    x
}

# A collection as saved before its references became per-row: one reference
# (or NULL) in the slot and no row key.
.savedWithOneLdRef <- function(x, ref) {
    old <- x
    methods::slot(old, "ldReference", check = FALSE) <- ref
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
# claiming `format`, which may be a format a handle no longer reads, its
# variant table whole and the seed naming its rows -- with no metadata and
# only the SNP / A1 / A2 variant columns, and reading the per-chromosome
# files `chromPaths` when given. Its variants, samples and sample data are
# copied from `panel`, today's panel over the same files.
.savedHandlePanel <- function(panel, path, format, chromPaths = character(0)) {
    rr <- .withPanelLabels(SummarizedExperiment::rowRanges(panel))
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
    seed <- methods::new(
        "GhSeed",
        handle = `slot<-`(handle, "format", check = FALSE, value = format),
        dm = dim(panel),
        dn = list(rownames(panel), NULL)
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = DelayedArray::DelayedArray(seed)),
        rowRanges = S4Vectors::`mcols<-`(
            rr,
            value = S4Vectors::mcols(rr)[c("SNP", "A1", "A2")]
        ),
        colData = SummarizedExperiment::colData(panel)
    )
}

# `panel` as saved before panels shed their copies of their variants'
# labels: an SNP column repeating its row names, and a dosage seed naming
# its rows -- over `handle`, its variant table whole, for a panel read
# through a GenotypeHandle. The panel holds every row of its seed, in order;
# it may hold some of its columns (a mixture's cohort).
.savedLabelledPanel <- function(panel, handle = NULL) {
    seed <- .panelSeed(panel)
    labelled <- if (methods::is(seed, "GhSeed")) {
        methods::new(
            "GhSeed",
            handle = handle,
            dm = dim(seed),
            dn = list(rownames(panel), NULL)
        )
    } else {
        methods::initialize(
            seed,
            dn = list(rownames(panel), dimnames(seed)[[2L]])
        )
    }
    index <- .panelSeedIndex(panel)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            dosage = DelayedArray::DelayedArray(labelled)[,
                index[[2L]],
                drop = FALSE
            ]
        ),
        rowRanges = .withPanelLabels(SummarizedExperiment::rowRanges(panel)),
        colData = SummarizedExperiment::colData(panel),
        metadata = S4Vectors::metadata(panel)
    )
}

# `x` as saved with one row per (study, context, trait, method): each
# method's variants an element carrying its `weight` column (named by
# variant id when `named`, as rows once stored them), its payload the row's
# columns, plus any `extra` columns. Attributes are written directly, as a
# saved object holds them: that layout no longer validates.
.savedTwasPerMethod <- function(x, named = FALSE, extra = list()) {
    cells <- .twasCells(x)
    entries <- .collectionEntries(x)
    grs <- map(entries, .savedTwasElement, named = named)
    flat <- suppressWarnings(exec(c, !!!grs))
    md <- exec(
        S4Vectors::DataFrame,
        !!!c(
            list(
                study = cells$study,
                context = cells$context,
                trait = cells$trait,
                method = cells$method,
                fits = S4Vectors::SimpleList(map(entries, methodFits)),
                cvResult = S4Vectors::SimpleList(map(entries, cvResult)),
                standardized = map_lgl(entries, weightStandardized),
                dataType = S4Vectors::SimpleList(map(entries, weightsDataType)),
                ldReference = rep(NA_character_, nrow(cells))
            ),
            extra,
            list(check.names = FALSE)
        )
    )
    saved <- x
    attr(saved, "unlistData") <- flat
    attr(saved, "partitioning") <- IRanges::PartitioningByEnd(
        cumsum(lengths(grs))
    )
    attr(saved, "elementMetadata") <- md
    saved
}

.savedTwasElement <- function(e, named) {
    gr <- variants(e)
    if (named) {
        S4Vectors::mcols(gr)$weight <- weights(e)
    }
    gr
}
