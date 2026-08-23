#!/usr/bin/env Rscript
#
# migrate_qtl_dataset_to_mae.R
#
# One-time migration for the QtlDataset -> MultiAssayExperiment change, kept
# as the record of what was rebuilt and re-runnable without harm.
#
# QtlDataset used to hold `phenotypes` (a named list of SummarizedExperiments)
# and `genotypeCovariates` (a matrix) as its own slots. It now extends
# MultiAssayExperiment: each context is an experiment, the genotype handle
# backs a further `genotype` experiment, and the covariates are that
# experiment's colData. `keepSamples` is gone -- the sample set lives in
# colData / sampleMap.
#
# The shipped fixtures were saved under the old definition, so their slots no
# longer line up with the class. They are read with attr() rather than slot(),
# which needs a definition that no longer describes them, and rebuilt through
# the constructor:
#
#   qtlDatasetExample            one study, one context
#   multiStudyQtlDatasetExample  two studies, each a QtlDataset
#
# Usage
#   GITHUB_WORKSPACE=$PWD pixi run --environment r45 Rscript \
#     inst/scripts/migrate_qtl_dataset_to_mae.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

dataDir <- "data"

# TRUE for an object saved under the pre-MAE definition: the old slots are
# still attributes, and the MAE ones were never written.
isLegacyQtlDataset <- function(x) {
    !is.null(attr(x, "phenotypes")) && is.null(attr(x, "ExperimentList"))
}

# Rebuild one QtlDataset through the constructor. Idempotent: an object
# already carrying the MAE slots is returned untouched.
migrateQtlDataset <- function(x) {
    if (!isLegacyQtlDataset(x)) {
        return(NULL)
    }
    QtlDataset(
        study = attr(x, "study"),
        genotypes = attr(x, "genotypes"),
        phenotypes = attr(x, "phenotypes"),
        genotypeCovariates = attr(x, "genotypeCovariates"),
        scaleResiduals = attr(x, "scaleResiduals"),
        mafCutoff = attr(x, "mafCutoff"),
        macCutoff = attr(x, "macCutoff"),
        xvarCutoff = attr(x, "xvarCutoff"),
        imissCutoff = attr(x, "imissCutoff"),
        keepSamples = attr(x, "keepSamples"),
        keepVariants = attr(x, "keepVariants"),
        keepIndel = attr(x, "keepIndel")
    )
}

# A MultiStudyQtlDataset is a container of QtlDatasets, so migrating it is
# migrating each member. Its own slots did not change.
migrateMultiStudy <- function(x) {
    members <- attr(x, "qtlDatasets")
    if (!any(map_lgl(members, isLegacyQtlDataset))) {
        return(NULL)
    }
    MultiStudyQtlDataset(
        qtlDatasets = map(members, migrateOrKeep),
        sumStats = attr(x, "sumStats")
    )
}

# @noRd
migrateOrKeep <- function(x) {
    out <- migrateQtlDataset(x)
    if (is.null(out)) x else out
}

targets <- list(
    qtlDatasetExample = migrateQtlDataset,
    multiStudyQtlDatasetExample = migrateMultiStudy
)

for (nm in names(targets)) {
    env <- new.env()
    data(list = nm, envir = env)
    before <- get(nm, envir = env)
    after <- targets[[nm]](before)
    if (is.null(after)) {
        cat(sprintf("  %-29s already migrated, skipped\n", nm))
        next
    }
    assign(nm, after)
    save(
        list = nm,
        file = file.path(dataDir, str_c(nm, ".rda")),
        compress = "xz"
    )
    cat(sprintf("  %-29s rebuilt\n", nm))
}

cat("\nDone.\n")
