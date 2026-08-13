# Column names referenced bare inside dplyr data-masking verbs (filter(),
# mutate(), arrange(), ...) and the magrittr `.` placeholder. Declaring them
# here keeps R CMD check's codetools pass from flagging them as undefined
# global variables. These are NOT package objects -- they are non-standard
# evaluation symbols resolved at run time against the data frame in scope.
utils::globalVariables(c(
    ".",
    "A1",
    "A1.ref",
    "A1.target",
    "A2",
    "A2.ref",
    "A2.target",
    "chrom",
    "pos",
    "gpos",
    "id",
    "path",
    "variant_id",
    "exact_match",
    "ID_match",
    "INDEL",
    "flip1.ref",
    "flip2.ref",
    "sign_flip",
    "strand_flip",
    "strand_unambiguous",
    "variants_id_original",
    "variants_id_qced",
    "imputed_z",
    "original_z",
    "outlier_stat",
    "rsq",
    "z_diff",
    "index_global",
    "index_within_window"
))
