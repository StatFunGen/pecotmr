# Validity of the method-selection Params lives in R/MethodParam.R, shared
# by FineMappingMethodsParam and TwasWeightsMethodsParam.

test_that("naming a method selects it, in any of the three slots", {
    p <- FineMappingMethodsParam(methods = "susie")
    expect_s4_class(p, "FineMappingMethodsParam")
    expect_s4_class(p, "MethodsSelectionParam")
    expect_s4_class(p, "MethodParam")
    expect_named(p$methods, "susie")
    expect_equal(
        FineMappingMethodsParam(methods = list(susie = list()))$methods,
        p$methods
    )
    expect_null(p$qtlDatasetMethods)
    expect_setequal(names(p), "methods")
})

test_that("one method may appear in both path slots, with its own options", {
    # The point of having two: on a dual-type input the same method reaches
    # a different engine on each path, so it needs separate arguments.
    p <- FineMappingMethodsParam(
        qtlDatasetMethods = list(susie = list(L = 20)),
        qtlSumStatsMethods = list(susie = list(L = 10))
    )
    expect_equal(p$qtlDatasetMethods$susie$L, 20)
    expect_equal(p$qtlSumStatsMethods$susie$L, 10)
    # Entries are stored as written: a per-method bundle is legitimately
    # mixed -- the wrapper's own formals alongside the engine's arguments --
    # so coercing it into the engine's constructor would reject the wrapper
    # half. Only the NAMES are judged.
    expect_type(p$qtlDatasetMethods$susie, "list")
    expect_false(isS4(p$qtlDatasetMethods$susie))
    # A name belonging to the other path is refused at construction, which
    # the union-validated aggregator accepted.
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(susie = list(z = 1))),
        "unknown argument\\(s\\) z for susie on QtlDataset"
    )
    # ... and a wrapper-level formal is accepted, not just engine arguments.
    expect_s4_class(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(mrmash = list(fitRetention = "slim"))
        ),
        "TwasWeightsMethodsParam"
    )
})

test_that("a method in `methods` may not also name an input path", {
    expect_error(
        FineMappingMethodsParam(
            methods = list(susie = list()),
            qtlDatasetMethods = list(susie = list())
        ),
        "named in `methods` and also under an input path"
    )
})

test_that("a path slot refuses options built for the other path", {
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = SusieRssOptions())
        ),
        "but `qtlDatasetMethods` is forwarded to QtlDataset"
    )
    expect_error(
        FineMappingMethodsParam(
            qtlSumStatsMethods = list(susie = SusieOptions())
        ),
        "options are for QtlDataset"
    )
})

test_that("a path slot offers only the methods that run on it", {
    expect_error(
        FineMappingMethodsParam(qtlSumStatsMethods = list(fsusie = list())),
        "fsusie is not a method this pipeline runs on QtlSumStats"
    )
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(ser = list())),
        "ser is not a method this pipeline runs on QtlDataset"
    )
    # `methods` states no path, so it accepts either; the pipeline's own
    # input class decides whether it can be honoured.
    expect_s4_class(
        FineMappingMethodsParam(methods = list(ser = list())),
        "FineMappingMethodsParam"
    )
})

test_that("`methods` entries must share one input class", {
    expect_error(
        FineMappingMethodsParam(
            methods = list(susie = SusieOptions(), ser = SerOptions())
        ),
        "cannot all be used with one input class"
    )
    expect_s4_class(
        FineMappingMethodsParam(
            methods = list(susie = SusieRssOptions(), ser = SerOptions())
        ),
        "FineMappingMethodsParam"
    )
})

test_that("`methods` defers name checking to the pipeline", {
    # `methods` names no input path, so there is nothing to check the
    # argument names against here -- the pipeline does it once dispatch has
    # settled the input class (.fmCheckMethodArgsForInput).
    expect_s4_class(
        FineMappingMethodsParam(methods = list(susie = list(L = 20))),
        "FineMappingMethodsParam"
    )
    expect_equal(
        pecotmr:::.fmNormalizeMethods(
            list(susie = list(z = 1)),
            inputKind = "QtlSumStats"
        )$methodArgs$susie$z,
        1
    )
    expect_error(
        pecotmr:::.fmNormalizeMethods(
            list(susie = list(z = 1)),
            inputKind = "QtlDataset"
        ),
        "unknown argument"
    )
})

test_that("malformed slots are rejected with the caller's own error", {
    expect_error(
        FineMappingMethodsParam(methods = list(list())),
        "every `methods` entry must be named"
    )
    expect_error(
        FineMappingMethodsParam(methods = 42),
        "must be a named list of per-method options"
    )
    expect_error(
        FineMappingMethodsParam(methods = "nosuch"),
        "nosuch is not a method this pipeline runs"
    )
    dup <- list(list(), list())
    names(dup) <- c("susie", "susie")
    expect_error(FineMappingMethodsParam(methods = dup), "given more than once")
    # Not wrapped in purrr's "In index: 1" framing.
    expect_no_match(
        tryCatch(
            FineMappingMethodsParam(
                qtlDatasetMethods = list(susie = list(LL = 1))
            ),
            error = function(e) conditionMessage(e)
        ),
        "In index"
    )
})

test_that("TwasWeightsMethodsParam shares the rules with its own registry", {
    p <- TwasWeightsMethodsParam(
        methods = list(susie = list()),
        qtlDatasetMethods = list(lasso = GlmnetOptions(alpha = 0.5)),
        qtlSumStatsMethods = list(lasso = LassosumOptions())
    )
    expect_s4_class(p, "TwasWeightsMethodsParam")
    expect_equal(p$qtlDatasetMethods$lasso$alpha, 0.5)
    # The registries differ: prsCs is summary-statistics only, enet
    # individual only, and both are methods a fine-mapping Param never has.
    expect_error(
        TwasWeightsMethodsParam(qtlDatasetMethods = list(prsCs = list())),
        "prsCs is not a method this pipeline runs on QtlDataset"
    )
    expect_error(
        TwasWeightsMethodsParam(qtlSumStatsMethods = list(enet = list())),
        "enet is not a method this pipeline runs on QtlSumStats"
    )
    expect_error(
        FineMappingMethodsParam(methods = list(lasso = list())),
        "lasso is not a method this pipeline runs"
    )
    # Cross-package per-path pairing, which one bag per method could not do.
    expect_error(
        TwasWeightsMethodsParam(
            qtlSumStatsMethods = list(lasso = GlmnetOptions())
        ),
        "options are for QtlDataset"
    )
})

test_that("a pipeline resolves only the slots its input class reads", {
    p <- FineMappingMethodsParam(
        methods = list(ser = list()),
        qtlDatasetMethods = list(susie = list(L = 20)),
        qtlSumStatsMethods = list(susieInf = list(L = 10))
    )
    # `methods` applies to every class; the other path's slot is not this
    # run's business.
    ind <- pecotmr:::.methodsParamResolve(p, "QtlDataset")
    expect_setequal(ind$tokens, c("ser", "susie"))
    expect_equal(ind$methodArgs$susie$L, 20)

    rss <- pecotmr:::.methodsParamResolve(p, "QtlSumStats")
    expect_setequal(rss$tokens, c("ser", "susieInf"))
    expect_equal(rss$methodArgs$susieInf$L, 10)

    # GwasSumStats is a summary-statistics class, so it reads the same slot.
    expect_setequal(
        pecotmr:::.methodsParamResolve(p, "GwasSumStats")$tokens,
        c("ser", "susieInf")
    )
})

test_that("a pipeline refuses a Param built for the other pipeline", {
    expect_error(
        pecotmr:::.methodsParamFor(
            TwasWeightsMethodsParam(methods = "lasso"),
            "QtlDataset",
            "FineMappingMethodsParam",
            "fineMappingPipeline"
        ),
        "configures a different pipeline"
    )
    # Its own family passes through untouched.
    p <- FineMappingMethodsParam(methods = "susie")
    expect_identical(
        pecotmr:::.methodsParamFor(
            p,
            "QtlDataset",
            "FineMappingMethodsParam",
            "fineMappingPipeline"
        ),
        p
    )
})

test_that("an entry must be built by that method's own constructor", {
    skip_if_not_installed("susieR")
    # The input-type check cannot see this: lasso and bayesA both run on
    # individual-level data, so GlmnetOptions() and QggOptions() are
    # indistinguishable by path and only the engine tells them apart.
    expect_s4_class(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(lasso = GlmnetOptions())
        ),
        "TwasWeightsMethodsParam"
    )
    expect_error(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(lasso = QggOptions())
        ),
        "built with the constructor for 'qgg'.*lasso on QtlDataset"
    )
    expect_error(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(bayesA = GlmnetOptions())
        ),
        "built with the constructor for 'glmnet'"
    )
    # Each path wants that path's own constructor.
    expect_s4_class(
        TwasWeightsMethodsParam(
            qtlSumStatsMethods = list(lasso = LassosumOptions())
        ),
        "TwasWeightsMethodsParam"
    )

    # susie and susieInf BOTH run susieR::susie on individual data, so the
    # callee cannot separate them -- only the token's own constructor pair
    # can, which is what the expected engine is derived from.
    expect_s4_class(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = SusieOptions())
        ),
        "FineMappingMethodsParam"
    )
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = SusieInfOptions())
        ),
        "built with the constructor for 'susieInf'"
    )
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = MvsusieOptions())
        ),
        "built with the constructor for 'mvsusie'"
    )
    expect_s4_class(
        FineMappingMethodsParam(
            qtlSumStatsMethods = list(susie = SusieRssOptions())
        ),
        "FineMappingMethodsParam"
    )
})

test_that("a fit-derived method is selected but never configured", {
    # These extract weights from a fit a prior fineMappingPipeline() run
    # produced, so they have no arguments of their own. Naming one with
    # list() SELECTS it -- the only thing it can say -- and anything else is
    # a misunderstanding the error names.
    expect_s4_class(
        TwasWeightsMethodsParam(qtlDatasetMethods = list(mvsusie = list())),
        "TwasWeightsMethodsParam"
    )
    expect_error(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(mvsusie = list(L = 1))
        ),
        "extracts weights from a fit supplied via"
    )
    # Path-independent: the answer does not depend on the input class, so
    # `methods` is held to it too.
    expect_error(
        TwasWeightsMethodsParam(methods = list(susie = list(L = 1))),
        "extracts weights from a fit supplied via"
    )
    # A method that does have arguments is unaffected.
    expect_s4_class(
        TwasWeightsMethodsParam(
            qtlDatasetMethods = list(lasso = list(alpha = 0.5))
        ),
        "TwasWeightsMethodsParam"
    )
})

test_that("a MultiStudyQtlDataset selects per half with path keys", {
    # Bare names run on whatever the dataset holds; a path key names methods
    # for that half only. fsusie has no summary-statistics implementation
    # and ser no individual one, so this is the only way to ask for both.
    p <- pecotmr:::.methodsParamForMulti(
        list("susie", qtlDataset = "fsusie", qtlSumStats = "ser"),
        "FineMappingMethodsParam",
        "fineMappingPipeline",
        TRUE
    )
    expect_named(p$methods, "susie")
    expect_named(p$qtlDatasetMethods, "fsusie")
    expect_named(p$qtlSumStatsMethods, "ser")
})

test_that("shared arguments are distributed to both halves", {
    skip_if_not_installed("susieR")
    # A method named at the top level carries arguments both halves share.
    # They are copied into both path slots, so each is checked against the
    # engine that half actually runs -- which is what makes a shared
    # argument have to be valid on BOTH.
    p <- pecotmr:::.methodsParamForMulti(
        list(susie = list(L = 20), qtlSumStats = list(susie = list(z = 1))),
        "FineMappingMethodsParam",
        "fineMappingPipeline",
        TRUE
    )
    expect_equal(p$qtlDatasetMethods$susie$L, 20)
    expect_equal(p$qtlSumStatsMethods$susie$L, 20)
    # The per-path extra lands on that half only.
    expect_null(p$qtlDatasetMethods$susie$z)
    expect_equal(p$qtlSumStatsMethods$susie$z, 1)
    # A per-path value of the same name wins over the shared one.
    q <- pecotmr:::.methodsParamForMulti(
        list(susie = list(L = 20), qtlSumStats = list(susie = list(L = 5))),
        "FineMappingMethodsParam",
        "fineMappingPipeline",
        TRUE
    )
    expect_equal(q$qtlDatasetMethods$susie$L, 20)
    expect_equal(q$qtlSumStatsMethods$susie$L, 5)
    # `z` is susie_rss-only, so sharing it is a mistake: distributing it to
    # the individual half fails, which is the intersection rule in action.
    expect_error(
        pecotmr:::.methodsParamForMulti(
            list(susie = list(z = 1)),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            TRUE
        ),
        "unknown argument"
    )
})

test_that("the two MultiStudy list shapes may not be mixed", {
    expect_error(
        pecotmr:::.methodsParamForMulti(
            list("susie", susieInf = list(L = 2)),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            TRUE
        ),
        "mixes two shapes"
    )
    expect_error(
        pecotmr:::.methodsParamForMulti(
            list(qtlDataset = 42),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            TRUE
        ),
        "must be a character vector of methods or a named list"
    )
})

test_that("summary-statistics methods need a dataset that carries some", {
    expect_error(
        pecotmr:::.methodsParamForMulti(
            list(qtlSumStats = "ser"),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            FALSE
        ),
        "carries no summary statistics"
    )
    # With sumstats present it is fine, and the individual half never needs
    # them.
    expect_s4_class(
        pecotmr:::.methodsParamForMulti(
            list(qtlSumStats = "ser"),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            TRUE
        ),
        "FineMappingMethodsParam"
    )
    expect_s4_class(
        pecotmr:::.methodsParamForMulti(
            list(qtlDataset = "fsusie"),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            FALSE
        ),
        "FineMappingMethodsParam"
    )
})

test_that("a Param keeps its structure through the MultiStudy helpers", {
    # The joint phase runs mrmash itself, so the per-method loop must not
    # also run it -- and dropping it has to leave a Param, or the
    # per-component recursion loses which path each entry was for.
    p <- TwasWeightsMethodsParam(
        methods = list(mrmash = list()),
        qtlDatasetMethods = list(lasso = list())
    )
    expect_setequal(pecotmr:::.methodsParamTokens(p), c("mrmash", "lasso"))
    stripped <- pecotmr:::.twasMsStripMrmash(p)
    expect_s4_class(stripped, "TwasWeightsMethodsParam")
    expect_setequal(pecotmr:::.methodsParamTokens(stripped), "lasso")
    # Stripping the only method leaves a Param that reports itself empty --
    # which is how the joint phase knows there is no per-method work left.
    only <- TwasWeightsMethodsParam(methods = list(mrmash = list()))
    expect_true(pecotmr:::.twasMethodsEmpty(
        pecotmr:::.twasMsStripMrmash(only)
    ))
})

test_that("show renders a slot holding per-method entries", {
    skip_if_not_installed("susieR")
    # format() cannot render a list of S4 records -- it falls through to
    # as.character() and errors -- so a slot holding per-method entries is
    # shown by naming them. Printing a populated Param used to fail, and
    # only the examples caught it: nothing called show() on one.
    out <- paste(
        capture.output(show(
            FineMappingMethodsParam(
                methods = list(susie = SusieOptions(L = 20))
            )
        )),
        collapse = " "
    )
    expect_match(out, "FineMappingMethodsParam")
    expect_match(out, "methods\\s+susie")
    expect_match(out, "unset: qtlDatasetMethods, qtlSumStatsMethods")

    multi <- paste(
        capture.output(show(TwasWeightsMethodsParam(
            methods = list(susie = list()),
            qtlDatasetMethods = list(lasso = GlmnetOptions(alpha = 0.5)),
            qtlSumStatsMethods = list(lasso = LassosumOptions())
        ))),
        collapse = " "
    )
    expect_match(multi, "qtlDatasetMethods\\s+lasso")
    expect_match(multi, "qtlSumStatsMethods\\s+lasso")

    # A settings Param still prints its values, not its names.
    expect_match(
        paste(
            capture.output(show(CredibleSetParam(coverage = 0.9))),
            collapse = " "
        ),
        "coverage\\s+0.9"
    )
})

# ===========================================================================
# show() and the slot-rendering helpers
# ---------------------------------------------------------------------------
# A populated Param used to fail to print: format() falls through to
# as.character() on an S4 slot value and errors. Nothing printed one, so the
# whole suite stayed green while show() was broken. These pin each rendering
# branch -- a nested Options record, a nested Param, a named list, an
# unnamed list -- so the next one cannot regress silently.
# ===========================================================================

test_that("show renders a Param with nothing set", {
    out <- capture.output(show(MashComponentParam()))
    expect_true(any(str_detect(out, "MashComponentParam")))
    expect_true(any(str_detect(out, fixed("(nothing set)"))))
    expect_true(any(str_detect(out, "unset: components")))
})

test_that("show names a nested Options record by its engine", {
    # Defaults carry no settings, so the record is reported as such.
    bare <- capture.output(show(SusieRssParam(
        control = SusieRssControlOptions()
    )))
    expect_true(any(str_detect(bare, fixed("<susieRssControl: defaults>"))))
    # With settings, they are listed instead.
    set <- capture.output(show(SusieRssParam(
        control = SusieRssControlOptions(r_tol = 1e-08)
    )))
    expect_true(any(str_detect(set, "susieRssControl: r_tol")))
})

test_that("show names a nested Param by what it contains", {
    # Not by slot name: the reader wants the fields, not "credibleSetArgs".
    out <- capture.output(show(GwasFineMappingParam()))
    expect_true(any(str_detect(out, "credibleSetArgs")))
    expect_true(any(str_detect(out, "coverage, secondaryCoverage")))
})

test_that("show counts an unnamed list slot and names a named one", {
    unnamed <- capture.output(show(CrossValidationParam(
        samplePartition = list(1:2, 3:4)
    )))
    expect_true(any(str_detect(unnamed, fixed("2 entries"))))
    one <- capture.output(show(CrossValidationParam(
        samplePartition = list(1:2)
    )))
    expect_true(any(str_detect(one, fixed("1 entry"))))
    named <- capture.output(show(FineMappingMethodsParam(
        qtlDatasetMethods = list(susie = SusieOptions())
    )))
    expect_true(any(str_detect(named, "susie")))
})

test_that(".paramNestedText reports a nested Param that is empty", {
    expect_match(
        as.character(pecotmr:::.paramNestedText(MashComponentParam())),
        "nothing set"
    )
})

test_that("$ and [[ refuse names and positions a Param does not have", {
    cs <- CredibleSetParam()
    expect_error(cs$nosuch, "is not a setting of CredibleSetParam")
    expect_error(cs$nosuch, "It has:")
    expect_error(cs[[99L]], "subscript out of bounds")
    expect_equal(cs[["coverage"]], 0.95)
    expect_error(cs[["nosuch"]], "is not a setting of CredibleSetParam")
    expect_error(cs[[NA_integer_]], "subscript out of bounds")
})

test_that(".assetMethodParam names the record it was handed", {
    # A Param of the wrong class reports its own class, not "a bare list".
    expect_error(
        pecotmr:::.assertMethodParam(
            CredibleSetParam(),
            "ColocPriorParam",
            "priors"
        ),
        "a CredibleSetParam"
    )
    expect_error(
        pecotmr:::.assertMethodParam(list(a = 1), "ColocPriorParam", "priors"),
        "a bare list"
    )
    # An empty list is "no settings", as it is for MethodOptions.
    expect_silent(pecotmr:::.assertMethodParam(
        list(),
        "ColocPriorParam",
        "priors"
    ))
})

test_that("validity catches unnamed entries that bypass the constructor", {
    # Two different checks produce near-identical text, and only one is on
    # the constructor path: .methodsNormalizeSlot() aborts first, so the
    # validity guard in .methodsSlotProblems() is reached only by a direct
    # new() or a slot<- afterwards. Asserting on the phrase alone cannot
    # tell them apart -- the prefix is what distinguishes them.
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(list(L = 1))),
        "FineMappingMethodsParam: every `qtlDatasetMethods` entry must be"
    )
    expect_error(
        new(
            "FineMappingMethodsParam",
            qtlDatasetMethods = list(list(L = 1))
        ),
        "`qtlDatasetMethods`: every entry must be named for its method"
    )
    # An empty-string name is the same failure as a missing one.
    expect_error(
        new(
            "FineMappingMethodsParam",
            qtlDatasetMethods = set_names(list(list()), "")
        ),
        "`qtlDatasetMethods`: every entry must be named for its method"
    )
})

test_that(".methodsParamFor and ...ForMulti route each input shape", {
    f <- pecotmr:::.methodsParamFor
    m <- pecotmr:::.methodsParamForMulti
    # A record of the right class passes straight through.
    p <- FineMappingMethodsParam(methods = "susie")
    expect_identical(f(p, "QtlDataset", "FineMappingMethodsParam", "x"), p)
    expect_identical(m(p, "FineMappingMethodsParam", "x", TRUE), p)
    # A character vector names methods with their defaults.
    expect_equal(
        names(f("susie", "QtlDataset", "FineMappingMethodsParam", "x")$methods),
        "susie"
    )
    # A record for the OTHER pipeline is refused by name, not silently used.
    expect_error(
        f(
            TwasWeightsMethodsParam(methods = "lasso"),
            "QtlDataset",
            "FineMappingMethodsParam",
            "fineMappingPipeline"
        ),
        "configures a different pipeline"
    )
    expect_error(
        m(
            TwasWeightsMethodsParam(methods = "lasso"),
            "FineMappingMethodsParam",
            "fineMappingPipeline",
            TRUE
        ),
        "configures a different pipeline"
    )
    # Anything that is neither a record, a character vector nor a list.
    expect_error(
        m(42, "FineMappingMethodsParam", "fineMappingPipeline", TRUE),
        "must be a character vector, a named list"
    )
})
