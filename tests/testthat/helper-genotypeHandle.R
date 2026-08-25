# The internal handle reader.
#
# readGenotypes() is the exported entry point and returns a genotype panel: a
# RangedSummarizedExperiment whose dosage assay reads lazily through a handle.
# The handle is the seed layer beneath it and is not public.
#
# Tests that exercise the reader machinery itself -- format detection, snpInfo,
# sample ids, sharded routing -- assert on the handle, so they read one
# directly rather than unwrapping a panel on every line.
# @noRd
readGenotypeHandle <- function(path, format = NULL, ...) {
    pecotmr:::.readGenotypeHandle(path, format = format, ...)
}
