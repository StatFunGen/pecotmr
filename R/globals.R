# Symbols referenced bare inside dplyr data-masking verbs (filter(), mutate(),
# arrange(), distinct(), ...) and the magrittr `.` placeholder, declared here so
# R CMD check's codetools pass does not flag them as undefined global variables.
# They are NOT package objects -- they are non-standard-evaluation symbols
# resolved at run time against the data frame in scope.
#
# Kept deliberately SHORT: the allele-harmonization / variant-id machinery and
# most other data-masking code now use the explicit `.data$col` pronoun (see
# `@importFrom rlang .data`), which is typo-safe and needs no entry here. Only
# these few remain bare -- pervasive coordinate/id columns threaded through the
# LD / RAISS / window helpers -- plus the magrittr `.` placeholder.
#
# NOTE when editing: audit with a `codetools::checkUsage(..., all = TRUE)` sweep
# on the SOURCE (e.g. R CMD check, or load with JIT off) -- `all = FALSE` and a
# byte-compiled `load_all` namespace UNDER-report and will hide genuinely-used
# columns.
utils::globalVariables(c(
    ".",
    "chrom",
    "pos",
    "variant_id",
    "index_global",
    "index_within_window"
))
