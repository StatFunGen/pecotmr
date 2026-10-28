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

test_that("a method's names are checked against the path that runs it", {
    skip_if_not_installed("susieR")
    indiv <- pecotmr:::.methodsPathAccepted(
        "FineMappingMethodsParam",
        "susie",
        "QtlDataset"
    )
    sumst <- pecotmr:::.methodsPathAccepted(
        "FineMappingMethodsParam",
        "susie",
        "QtlSumStats"
    )
    # `z` is a susie_rss argument and not a susie one.
    expect_false(is_in("z", indiv))
    expect_true(is_in("z", sumst))
    # `y` is susie-only, so the mirror image holds.
    expect_true(is_in("y", indiv))
    expect_false(is_in("y", sumst))
    # A name both callees accept is fine either way.
    expect_true(is_in("L", indiv))
    expect_true(is_in("L", sumst))
    # Why this matters: validating against the UNION of both callees accepts
    # `z` even for an individual-level run, and the setting then reaches an
    # engine that never reads it.
    both <- pecotmr:::.engineAcceptedNames(
        c("susieR::susie", "susieR::susie_rss")
    )
    expect_true(is_in("z", both))
    expect_false(is_in("z", pecotmr:::.engineAcceptedNames("susieR::susie")))
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
        qtlDatasetMethods = list(lasso = GlmnetOptions(nfolds = 5)),
        qtlSumStatsMethods = list(lasso = LassosumOptions())
    )
    expect_s4_class(p, "TwasWeightsMethodsParam")
    expect_equal(p$qtlDatasetMethods$lasso$nfolds, 5)
    # The registries differ: prsCs is summary-statistics only, dprVb
    # individual only, and both are methods a fine-mapping Param never has.
    expect_error(
        TwasWeightsMethodsParam(qtlDatasetMethods = list(prsCs = list())),
        "prsCs is not a method this pipeline runs on QtlDataset"
    )
    expect_error(
        TwasWeightsMethodsParam(qtlSumStatsMethods = list(dprVb = list())),
        "dprVb is not a method this pipeline runs on QtlSumStats"
    )
    # enet and ridge run on both, glmnet's options on one path and the
    # penalized RSS solver's on the other.
    both <- TwasWeightsMethodsParam(
        qtlDatasetMethods = list(enet = GlmnetOptions(), ridge = list()),
        qtlSumStatsMethods = list(
            enet = PenalizedRssOptions(thr = 1e-6),
            ridge = list()
        )
    )
    expect_s4_class(both, "TwasWeightsMethodsParam")
    expect_error(
        TwasWeightsMethodsParam(
            qtlSumStatsMethods = list(enet = GlmnetOptions())
        ),
        "options are for QtlDataset"
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
                methods = list(susie = SusieOptions(max_iter = 20))
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
            qtlDatasetMethods = list(lasso = GlmnetOptions(nfolds = 5)),
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
    # Not by slot name: the reader wants the fields, not the bundle name.
    out <- capture.output(show(GwasFineMappingParam()))
    expect_true(any(str_detect(out, "credibleSetParam")))
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

test_that("show names a BPPARAM by its backend and worker count", {
    out <- capture.output(show(CrossValidationParam()))
    expect_true(any(str_detect(out, fixed("SerialParam (1 worker)"))))
    two <- capture.output(show(CrossValidationParam(
        BPPARAM = BiocParallel::MulticoreParam(2)
    )))
    expect_true(any(str_detect(two, fixed("MulticoreParam (2 workers)"))))
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

# ---------------------------------------------------------------------------
# $<-
# ---------------------------------------------------------------------------

# Every concrete MethodParam class whose constructor builds with defaults.
.mp_defaultParams <- function() {
    ns <- asNamespace("pecotmr")
    classes <- names(methods::getClass("MethodParam")@subclasses)
    concrete <- Filter(function(cl) !methods::getClass(cl)@virtual, classes)
    built <- lapply(concrete, function(cl) {
        if (!exists(cl, envir = ns, mode = "function")) {
            return(NULL)
        }
        tryCatch(get(cl, envir = ns)(), error = function(e) NULL)
    })
    names(built) <- concrete
    Filter(Negate(is.null), built)
}

test_that("$<- writing a setting back to itself changes nothing", {
    # The rebuild goes through each class's own constructor, so every
    # constructor must accept the values it stored -- normalised lists and
    # nested records included. A round trip that is not the identity would
    # mean `$<-` quietly reshapes a record.
    params <- .mp_defaultParams()
    expect_gt(length(params), 20L)
    for (cl in names(params)) {
        p <- params[[cl]]
        fixed <- names(pecotmr:::.PARAM_FIXED_FIELDS[[cl]])
        for (f in setdiff(methods::slotNames(p), fixed)) {
            q <- do.call(`$<-`, list(p, f, methods::slot(p, f)))
            expect_identical(q, p, info = paste0(cl, "$", f))
        }
    }
})

test_that("$<- normalises through the constructor and unsets with NULL", {
    # FineMappingMethodsParam() turns a method name into its options entry;
    # `$<-` stores exactly what the constructor would.
    p <- FineMappingMethodsParam(methods = "susie")
    p$methods <- c("susie", "susieInf")
    expect_identical(
        p,
        FineMappingMethodsParam(methods = c("susie", "susieInf"))
    )
    cs <- CredibleSetParam(medianAbsCorr = 0.5)
    expect_true("medianAbsCorr" %in% names(cs))
    cs$medianAbsCorr <- NULL
    expect_null(cs$medianAbsCorr)
    expect_false("medianAbsCorr" %in% names(cs))
})

test_that("$<- refuses unknown and fixed settings, naming why", {
    cs <- CredibleSetParam()
    expect_error(cs$nosuch <- 1, "is not a setting of CredibleSetParam")
    expect_error(
        cs$maxNumSingleEffects <- 3L,
        "needs a refit. Build a new CredibleSetParam\\(\\) instead"
    )
    mc <- MashComponentParam()
    expect_error(mc$components <- list(), "assembles it from")
    pv <- PvalMethodsParam(aspu = AspuOptions(n.perm = 5000))
    expect_error(pv$methods <- list(), "one argument per method")
})

test_that(".assertMethodParam names the record it was handed", {
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

# ---- argument naming -----------------------------------------------------

# Every closure pecotmr defines: the namespace's functions plus the methods of
# its own generics. A method whose signature differs from its generic's is
# stored as a wrapper around an inner `.local`, whose formals the caller sees.
.mpUnwrapLocal <- function(f) {
    b <- body(f)
    isLocal <- is.call(b) &&
        identical(b[[1]], as.name("{")) &&
        is.call(b[[2]]) &&
        identical(b[[2]][[2]], as.name(".local"))
    if (isLocal) eval(b[[2]][[3]]) else f
}

.mpMethodsOf <- function(g, ns) {
    defs <- methods::findMethods(g, where = ns)
    set_names(
        map(defs, .mpUnwrapLocal),
        str_c(g, "[", map_chr(defs, .mpSignatureLabel), "]")
    )
}

.mpSignatureLabel <- function(m) str_flatten(m@defined, ",")

.mpClosures <- function(exportedOnly = FALSE) {
    ns <- asNamespace("pecotmr")
    nms <- if (exportedOnly) {
        getNamespaceExports(ns)
    } else {
        ls(ns, all.names = TRUE)
    }
    fns <- keep(
        mget(nms, envir = ns, inherits = TRUE, ifnotfound = list(NULL)),
        is.function
    )
    gens <- intersect(as.character(methods::getGenerics(where = ns)), nms)
    c(fns, list_flatten(map(gens, .mpMethodsOf, ns = ns)))
}

# The symbols a call tree passes as a function VALUE: the `.f` / FUN slot of
# the iteration helpers and the first argument of exec() and friends.
.mpFunSlot <- c(
    exec = 1L,
    do.call = 1L,
    match.fun = 1L,
    Filter = 1L,
    Reduce = 1L,
    Map = 1L,
    mapply = 1L,
    map = 2L,
    map_lgl = 2L,
    map_int = 2L,
    map_dbl = 2L,
    map_chr = 2L,
    imap = 2L,
    walk = 2L,
    keep = 2L,
    discard = 2L,
    reduce = 2L,
    accumulate = 2L,
    detect = 2L,
    lapply = 2L,
    sapply = 2L,
    vapply = 2L,
    pmap = 2L,
    map2 = 3L,
    imap_chr = 2L,
    apply = 3L
)

.mpFunValueSymbols <- function(e) {
    if (!is.call(e)) {
        return(character(0))
    }
    head <- if (is.symbol(e[[1]])) as.character(e[[1]]) else ""
    here <- character(0)
    if (is_in(head, names(.mpFunSlot)) && length(e) > .mpFunSlot[[head]]) {
        named <- intersect(names(e), c(".f", ".fn", "FUN", "f", "what"))
        a <- if (length(named)) {
            e[[named[[1]]]]
        } else {
            e[[.mpFunSlot[[head]] + 1L]]
        }
        if (is.symbol(a)) here <- as.character(a)
    }
    rest <- map(as.list(e)[-1], .mpFunValueSymbolsSafe)
    unique(c(here, unlist(rest, use.names = FALSE)))
}

.mpFunValueSymbolsSafe <- function(a) {
    if (missing(a) || is.null(a)) character(0) else .mpFunValueSymbols(a)
}

.mpIsNsFunction <- function(n, ns) is.function(get(n, envir = ns))

.mpShadowedValues <- function(f, label, own) {
    shadowing <- intersect(names(formals(f)), own)
    if (length(shadowing) == 0L) {
        return(NULL)
    }
    # A formal is passed as a value only through one of the .mpFunSlot
    # helpers, and only if it appears in the body at all; all.names() rules
    # out most closures before their call tree is walked.
    syms <- all.names(body(f))
    if (!any(is_in(shadowing, syms)) || !any(is_in(names(.mpFunSlot), syms))) {
        return(NULL)
    }
    hit <- intersect(.mpFunValueSymbols(body(f)), shadowing)
    if (length(hit)) hit else NULL
}

test_that("no formal that shadows a pecotmr function is passed as a value", {
    # Arguments are named after their class (credibleSetParam =
    # CredibleSetParam()), so a few share a name with an accessor generic --
    # credibleSetParam(), susieRssParam(), panelFilterParam(). That is safe
    # while the function only CALLS the accessor, because R skips non-function
    # bindings when it looks up a call. It breaks the moment the accessor is
    # passed as a value -- map(objs, credibleSetParam) -- because the formal
    # is found instead. This pins the condition the shadowing relies on.
    ns <- asNamespace("pecotmr")
    own <- keep(ls(ns, all.names = TRUE), .mpIsNsFunction, ns = ns)
    # The detector is not vacuous: it flags exactly that shape.
    probe <- function(credibleSetParam, fits) map(fits, credibleSetParam)
    expect_equal(
        .mpShadowedValues(probe, "probe", own = own),
        "credibleSetParam"
    )
    closures <- .mpClosures()
    bad <- compact(imap(closures, .mpShadowedValues, own = own))
    expect_equal(
        str_c(names(bad), ": ", map_chr(bad, str_flatten, ", ")),
        character(0)
    )
})

.mpParamFormals <- function(f, label) {
    # formals() is a pairlist, which purrr deprecated as an input.
    fm <- as.list(formals(f))
    isParam <- map_lgl(fm, .mpIsParamCall)
    map(names(fm)[isParam], .mpFormalRecord, fm = fm, label = label)
}

.mpIsParamCall <- function(d) {
    is.call(d) &&
        is.symbol(d[[1]]) &&
        methods::extends(as.character(d[[1]]), "MethodParam")
}

.mpFormalRecord <- function(nm, fm, label) {
    list(fn = label, formal = nm, class = as.character(fm[[nm]][[1]]))
}

.mpDescribeFormal <- function(r) {
    str_c(r$fn, ": ", r$formal, " = ", r$class, "()")
}

.mpNamedAfterClass <- function(r, roleNamed) {
    expected <- str_c(tolower(substr(r$class, 1, 1)), substring(r$class, 2))
    identical(r$formal, expected) || is_in(r$formal, roleNamed)
}

test_that("an argument taking a *Param record is named after its class", {
    # Settings pecotmr reads travel as *Param records, passed under the
    # class's own name (crossValidationParam = CrossValidationParam()); engine
    # argument bundles (*Options) keep the *Args / methodArgs names. The four
    # below are named for their role instead, and are the only exceptions.
    roleNamed <- c("prior", "priors", "columns", "pruning")
    expect_false(.mpNamedAfterClass(
        list(formal = "crossValidationArgs", class = "CrossValidationParam"),
        roleNamed = roleNamed
    ))
    found <- list_flatten(imap(
        .mpClosures(exportedOnly = TRUE),
        .mpParamFormals
    ))
    misnamed <- discard(found, .mpNamedAfterClass, roleNamed = roleNamed)
    expect_gt(length(found), 0L)
    expect_equal(unname(map_chr(misnamed, .mpDescribeFormal)), character(0))
})
