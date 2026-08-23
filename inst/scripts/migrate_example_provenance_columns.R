#!/usr/bin/env Rscript
#
# migrate_example_provenance_columns.R
#
# One-time migration for the §4.4 column retirement, kept as the record of what
# was changed and re-runnable without harm.
#
# §4.4 retired two mcols columns: `region` (a per-row GRanges of provenance,
# replaced by `traitPos`) and `region_id` (replaced by `blockId`). The code was
# updated, but four shipped fixtures were built before that and still carry the
# old columns:
#
#   ctwasWeightsExample        region, no traitPos   <- functionally broken
#   qtlFineMappingLbfExample   region + traitPos     <- dead weight
#   gwasFineMappingLbfExample  region + region_id    <- dead weight
#   mvsusieFineMappingExample  region + traitPos     <- dead weight
#
# Only the first is a live defect: the cTWAS weight-placement path requires
# `traitPos`, so `assembleCtwasInputs()` rejects the fixture its own @examples
# pass in. The other three merely carry a column nothing reads any more.
#
# `region` and `traitPos` hold the same thing -- one GRanges per row naming
# where the trait sits -- so the migration is a rename, not a recomputation.
#
# Usage
#   pixi run --environment r45 Rscript \
#     inst/scripts/migrate_example_provenance_columns.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

dataDir <- "data"

# The joint-key columns are optional and only present on joint fits, so a
# constructor call that does not name them silently drops them -- which is how
# the first run of this script lost `jointContexts` from the mvsusie fixture.
jointArgs <- function(x) {
    cn <- colnames(x)
    out <- list()
    for (nm in c("jointStudies", "jointContexts", "jointTraits")) {
        if (is_in(nm, cn)) {
            out[[nm]] <- mcols(x)[[nm]]
        }
    }
    out
}

# Rebuild a TwasWeights with `region` promoted to `traitPos`. Idempotent: a
# fixture already carrying `traitPos` and no `region` is left untouched.
migrateTwasWeights <- function(x) {
    cn <- colnames(x)
    if (!is_in("region", cn)) {
        return(NULL)
    }
    traitPos <- if (is_in("traitPos", cn)) {
        mcols(x)$traitPos
    } else {
        mcols(x)$region
    }
    rows <- map(seq_len(nrow(x)), .twrRowParts, x = x)
    exec(
        TwasWeights,
        study = as.character(x$study),
        context = as.character(x$context),
        trait = as.character(x$trait),
        method = as.character(x$method),
        entry = rows,
        traitPos = traitPos,
        ldSketch = getLdSketch(x),
        !!!jointArgs(x)
    )
}

# The fine-mapping fixtures need the same treatment, minus the constructor
# difference: their `region` is dead weight beside an existing `traitPos`.
migrateFineMapping <- function(x) {
    cn <- colnames(x)
    if (!is_in("region", cn) && !is_in("region_id", cn)) {
        return(NULL)
    }
    traitPos <- if (is_in("traitPos", cn)) {
        mcols(x)$traitPos
    } else if (is_in("region", cn)) {
        mcols(x)$region
    } else {
        NULL
    }
    rows <- map(seq_len(nrow(x)), .fmrRowParts, x = x)
    if (methods::is(x, "GwasFineMappingResult")) {
        blockId <- if (is_in("blockId", cn)) {
            as.character(x$blockId)
        } else if (is_in("region_id", cn)) {
            as.character(mcols(x)$region_id)
        } else {
            NULL
        }
        return(GwasFineMappingResult(
            study = as.character(x$study),
            method = as.character(x$method),
            entry = rows,
            blockId = blockId,
            traitPos = traitPos,
            ldSketch = getLdSketch(x)
        ))
    }
    exec(
        QtlFineMappingResult,
        study = as.character(x$study),
        context = as.character(x$context),
        trait = as.character(x$trait),
        method = as.character(x$method),
        entry = rows,
        traitPos = traitPos,
        ldSketch = getLdSketch(x),
        !!!jointArgs(x)
    )
}

# A fixture predating the RangedTupleList migration is still a DFrame: columns
# in `listData`, payloads as `FineMappingEntry` objects whose class no longer
# exists. Its slots are all still there, so it can be rebuilt -- accessed with
# attr() rather than slot(), which needs a class definition that is gone.
isLegacyLayout <- function(x) {
    !is.null(attr(x, "listData")) && is.null(attr(x, "partitioning"))
}

migrateLegacyFineMapping <- function(x) {
    ld <- attr(x, "listData")
    entries <- attr(ld$entry, "listData")
    rows <- map(entries, legacyRow)
    args <- list(
        study = as.character(ld$study),
        context = as.character(ld$context),
        trait = as.character(ld$trait),
        method = as.character(ld$method),
        entry = rows,
        traitPos = ld$traitPos,
        ldSketch = attr(x, "ldSketch")
    )
    for (nm in c("jointStudies", "jointContexts", "jointTraits")) {
        if (!is.null(ld[[nm]])) {
            args[[nm]] <- ld[[nm]]
        }
    }
    exec(QtlFineMappingResult, !!!args)
}

# @noRd
legacyRow <- function(e) {
    fineMappingRow(
        variantIds = attr(e, "variantIds"),
        susieFit = attr(e, "susieFit"),
        topLoci = attr(e, "topLoci"),
        cvResult = attr(e, "cvResult")
    )
}

targets <- list(
    ctwasWeightsExample = migrateTwasWeights,
    qtlFineMappingLbfExample = migrateFineMapping,
    gwasFineMappingLbfExample = migrateFineMapping,
    mvsusieFineMappingExample = migrateFineMapping
)

for (nm in names(targets)) {
    env <- new.env()
    data(list = nm, envir = env)
    before <- get(nm, envir = env)
    after <- if (isLegacyLayout(before)) {
        migrateLegacyFineMapping(before)
    } else {
        targets[[nm]](before)
    }
    if (is.null(after)) {
        cat(sprintf("  %-27s already migrated, skipped\n", nm))
        next
    }
    beforeCols <- if (isLegacyLayout(before)) {
        setdiff(names(attr(before, "listData")), "entry")
    } else {
        colnames(before)
    }
    beforeRows <- if (isLegacyLayout(before)) {
        attr(before, "nrows")
    } else {
        nrow(before)
    }
    stopifnot(nrow(after) == beforeRows)
    retired <- c("region", "region_id")
    lost <- setdiff(setdiff(beforeCols, retired), colnames(after))
    if (length(lost) > 0L) {
        stop(
            nm,
            ": migration dropped column(s) ",
            str_flatten(lost, ", ")
        )
    }
    assign(nm, after)
    save(
        list = nm,
        file = file.path(dataDir, str_c(nm, ".rda")),
        compress = "xz"
    )
    cat(sprintf(
        "  %-27s %s -> %s\n",
        nm,
        str_flatten(beforeCols, ", "),
        str_flatten(colnames(after), ", ")
    ))
}

cat("\nDone.\n")
