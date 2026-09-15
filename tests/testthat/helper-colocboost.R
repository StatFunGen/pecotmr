# Raw colocboost run objects and the outcome-info lookup that
# ColocBoostResult and .cbToResultObject are assembled from.
#
# Shared because the class file and the pipeline file both need one: the
# class tests exercise the row expansion, the pipeline tests exercise how
# .cbToResultObject keys runs by analysis and gwasStudy.

# A minimal stand-in for a colocboost object, shaped like the real one:
# cos_details keyed by cos_id, purity as a SQUARE matrix over sets, and
# top variables as a data frame with set ids for rownames.
.cbr_fake <- function(
    ids = "cos1:y1_y2",
    outcomes = list(c("t1", "t2")),
    members = list(2L),
    variants = list("chr1:200:C:T"),
    npc = 0.9,
    nRegion = 4L,
    focal = FALSE
) {
    regionIds <- str_c("chr1:", seq_len(nRegion) * 100L, ":C:T")
    vcp <- set_names(seq_len(nRegion) / (nRegion * 2), regionIds)
    purity <- matrix(
        1,
        nrow = length(ids),
        ncol = length(ids),
        dimnames = list(ids, ids)
    )
    structure(
        list(
            cos_summary = tibble(
                cos_id = ids,
                focal_outcome = focal,
                top_variable = map_chr(variants, 1L),
                top_variable_vcp = rep(0.8, length(ids))
            ),
            vcp = vcp,
            cos_details = list(
                cos = list(
                    cos_index = set_names(members, ids),
                    cos_variables = set_names(variants, ids)
                ),
                cos_outcomes = list(
                    outcome_name = set_names(outcomes, ids)
                ),
                cos_vcp = set_names(
                    rep(list(as.numeric(vcp)), length(ids)),
                    ids
                ),
                cos_npc = set_names(rep(npc, length(ids)), ids),
                cos_min_npc_outcome = set_names(rep(npc, length(ids)), ids),
                cos_purity = list(min_abs_cor = purity),
                cos_top_variables = data.frame(
                    top_index = unlist(members),
                    top_variables = unlist(variants),
                    row.names = ids,
                    stringsAsFactors = FALSE
                )
            )
        ),
        class = "colocboost"
    )
}

.cbr_info <- function(names = c("t1", "t2")) {
    data.frame(
        name = names,
        context = str_c("ctx", seq_along(names)),
        trait = "GENE1",
        study = "study1",
        dataForm = "individual",
        stringsAsFactors = FALSE
    )
}
