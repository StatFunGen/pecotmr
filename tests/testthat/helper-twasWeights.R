# TwasWeights cut to some of their methods, for the tests that build one from
# twasWeightsExample (test_tupleSelectors.R, test_resultWriters.R). The package
# itself cuts weights to chosen (gene, method) cells, never to whole methods.

# `x` with only the methods `methods` (camelCase tokens).
.twasKeepMethods <- function(x, methods) {
    cells <- .twasCells(x)
    .twasKeepCells(x, cells[is_in(cells$method, methods), ])
}
