test_that(".engineCallee resolves pkg::fn and tolerates a missing package", {
    expect_identical(
        pecotmr:::.engineCallee("stats::var"),
        stats::var
    )
    # A package that cannot be there is not an error: building an argument
    # list is not running the engine.
    expect_null(pecotmr:::.engineCallee("notARealPackage::nope"))
    # Present package, absent function.
    expect_null(pecotmr:::.engineCallee("stats::notARealFunction"))
})

test_that(".engineCallee passes a function through unchanged", {
    # A caller that already holds the engine need not spell its package.
    expect_identical(pecotmr:::.engineCallee(stats::var), stats::var)
})

test_that(".newMethodOptions records NA for a callee given as a function", {
    res <- pecotmr:::.newMethodOptions(
        stats::var,
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs"
    )
    expect_true(is.na(S4Vectors::metadata(res)$callee))
    # Validation still happens -- the formals come from the function itself.
    expect_error(
        pecotmr:::.newMethodOptions(
            stats::var,
            defaults = list(),
            extra = list(bogus = 1),
            label = "varArgs"
        ),
        "unknown argument"
    )
})

test_that(".engineCallee rejects a malformed callee", {
    expect_error(
        pecotmr:::.engineCallee("susie_rss"),
        "must be 'pkg::fn'"
    )
})

test_that(".engineAcceptedNames returns formals, or NULL when unknowable", {
    expect_true(is_in("na.rm", pecotmr:::.engineAcceptedNames("stats::var")))
    # `paste` takes `...`, so it accepts anything and there is nothing to
    # validate against. So does the `mean` generic, whose `na.rm` lives on
    # mean.default rather than on the generic itself.
    expect_null(pecotmr:::.engineAcceptedNames("base::paste"))
    expect_null(pecotmr:::.engineAcceptedNames("base::mean"))
    expect_null(pecotmr:::.engineAcceptedNames("notARealPackage::nope"))
})

test_that(".newMethodOptions merges defaults under extras, dropping NULLs", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(na.rm = TRUE, use = NULL),
        extra = list(y = 1),
        label = "varArgs"
    )
    expect_true(pecotmr:::.isMethodOptions(res))
    # `use = NULL` means "leave it to the engine", so it is not forwarded.
    expect_false(is_in("use", names(res)))
    expect_true(res$na.rm)
    expect_equal(res$y, 1)
})

test_that(".newMethodOptions rejects a name the engine cannot accept", {
    expect_error(
        pecotmr:::.newMethodOptions(
            "stats::var",
            defaults = list(),
            extra = list(naRM = TRUE),
            label = "varArgs"
        ),
        "unknown argument\\(s\\) naRM"
    )
})

test_that(".newMethodOptions names the engine's accepted arguments in the error", {
    expect_error(
        pecotmr:::.newMethodOptions(
            "stats::var",
            defaults = list(),
            extra = list(nope = 1),
            label = "varArgs"
        ),
        "stats::var. accepts:.*na\\.rm"
    )
})

test_that(".newMethodOptions skips validation when the engine takes dots", {
    # base::paste has `...`, so an arbitrary name is legitimate.
    res <- pecotmr:::.newMethodOptions(
        "base::paste",
        defaults = list(sep = "-"),
        extra = list(anythingAtAll = 1),
        label = "pasteArgs"
    )
    expect_equal(res$anythingAtAll, 1)
})

test_that(".newMethodOptions skips validation when the package is absent", {
    res <- pecotmr:::.newMethodOptions(
        "notARealPackage::nope",
        defaults = list(),
        extra = list(whatever = 2),
        label = "nopeArgs"
    )
    expect_equal(res$whatever, 2)
})

test_that(".newMethodOptions requires every extra to be named", {
    expect_error(
        pecotmr:::.newMethodOptions(
            "base::paste",
            defaults = list(),
            extra = list(1, 2),
            label = "pasteArgs"
        ),
        "must be named"
    )
})

test_that(".newMethodOptions rejects an argument supplied twice", {
    # Two entries in `extra` -- what .methodOptionsWith() produces when a
    # wrapper's own contribution meets the caller's record -- have no winner.
    expect_error(
        pecotmr:::.newMethodOptions(
            "base::paste",
            defaults = list(),
            extra = list(sep = "-", sep = "+"),
            label = "pasteArgs"
        ),
        "given twice: sep"
    )
})

test_that(".newMethodOptions lets an explicit value override a default", {
    res <- pecotmr:::.newMethodOptions(
        "base::paste",
        defaults = list(sep = "-", collapse = ","),
        extra = list(sep = "+"),
        label = "pasteArgs"
    )
    expect_equal(res$sep, "+")
    expect_equal(res$collapse, ",")
    # The default keeps its position, so the record reads the same either way.
    expect_equal(names(res), c("sep", "collapse"))
})

test_that(".newMethodOptions records the engine and callee for reporting", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_equal(S4Vectors::metadata(res)$engine, "var")
    expect_equal(S4Vectors::metadata(res)$callee, "stats::var")
})

test_that("MethodOptions keeps its class through pecotmr's list operations", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(na.rm = TRUE, y = 1),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    # .fmNarrowAfterJoint / .fmQssJointPhase subset these records by name, and
    # a plain list carrying an S3 class loses the class right here.
    expect_true(pecotmr:::.isMethodOptions(res["na.rm"]))
    expect_true(pecotmr:::.isMethodOptions(res[1]))
    expect_equal(S4Vectors::metadata(res["na.rm"])$engine, "var")
})

test_that(".newMethodOptions splices into a call exactly like a bare list", {
    # The whole point of returning a plain list: every existing
    # exec(fn, !!!args) call site keeps working unchanged.
    res <- pecotmr:::.newMethodOptions(
        "base::paste",
        defaults = list(sep = "-"),
        extra = list(),
        label = "pasteArgs"
    )
    expect_equal(exec(paste, "a", "b", !!!res), "a-b")
})

test_that(".assertMethodOptions accepts a constructed record and empty list", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(),
        extra = list(),
        label = "varArgs"
    )
    expect_silent(pecotmr:::.assertMethodOptions(res, "varArgs", "methodArgs"))
    # `list()` carries no settings, so there is nothing a constructor would
    # have added and nothing to validate.
    expect_silent(pecotmr:::.assertMethodOptions(
        list(),
        "varArgs",
        "methodArgs"
    ))
})

test_that(".assertMethodOptions rejects a populated bare list", {
    expect_error(
        pecotmr:::.assertMethodOptions(
            list(na.rm = TRUE),
            "varArgs",
            "methodArgs"
        ),
        "`methodArgs` must be built with varArgs\\(\\)"
    )
})

test_that("show(MethodOptions) reports the engine and its settings", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_output(show(res), "MethodOptions: var")
    expect_output(show(res), "forwarded to stats::var")
    expect_output(show(res), "na\\.rm")
})

test_that("show(MethodOptions) says so when nothing is set", {
    res <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(),
        extra = list(),
        label = "varArgs"
    )
    expect_output(show(res), "no options set")
})

# --- character-or-constructor selection ------------------------------------

.ma_fooArgs <- function(a = 1, b = 2) {
    pecotmr:::.newMethodOptions(
        NULL,
        defaults = list(a = a, b = b),
        extra = list(),
        label = "fooArgs",
        engine = "foo"
    )
}

.ma_barArgs <- function(z = 9) {
    pecotmr:::.newMethodOptions(
        NULL,
        defaults = list(z = z),
        extra = list(),
        label = "barArgs",
        engine = "bar"
    )
}

.ma_ctors <- function() list(foo = .ma_fooArgs, bar = .ma_barArgs)

test_that(".engineOf reads the engine from a name or a record", {
    expect_equal(pecotmr:::.engineOf("foo"), "foo")
    expect_equal(pecotmr:::.engineOf(.ma_fooArgs()), "foo")
    expect_null(pecotmr:::.engineOf(list(a = 1)))
    expect_null(pecotmr:::.engineOf(c("foo", "bar")))
})

test_that(".resolveEngineChoice accepts a bare name with no options", {
    res <- pecotmr:::.resolveEngineChoice("foo", c("foo", "bar"), "engine")
    expect_equal(res$engine, "foo")
    expect_length(res$args, 0L)
    expect_true(pecotmr:::.isMethodOptions(res$args))
})

test_that(".resolveEngineChoice carries a constructor's options through", {
    res <- pecotmr:::.resolveEngineChoice(
        .ma_fooArgs(a = 3),
        c("foo", "bar"),
        "engine"
    )
    expect_equal(res$engine, "foo")
    expect_equal(res$args$a, 3)
})

test_that(".resolveEngineChoice rejects an engine outside the vocabulary", {
    expect_error(
        pecotmr:::.resolveEngineChoice("zzz", c("foo", "bar"), "engine"),
        "unknown engine 'zzz'"
    )
    # A record built for an engine this argument does not accept is caught by
    # name rather than left to fail further down.
    expect_error(
        pecotmr:::.resolveEngineChoice(.ma_barArgs(), c("foo"), "engine"),
        "unknown engine 'bar'"
    )
})

test_that(".resolveEngineChoice rejects a value that names no engine", {
    expect_error(
        pecotmr:::.resolveEngineChoice(list(a = 1), c("foo"), "engine"),
        "must be one of foo, or the matching constructor"
    )
})

# --- nested per-engine bundles ---------------------------------------------

test_that("a record says whether its names were checked", {
    checked <- pecotmr:::.newMethodOptions(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_output(show(checked), "argument names checked against stats::var")
    # An engine taking `...` accepts anything, so a misspelling reaches it.
    dots <- pecotmr:::.newMethodOptions(
        "base::paste",
        defaults = list(sep = "-"),
        extra = list(),
        label = "pasteArgs",
        engine = "paste"
    )
    expect_output(show(dots), "NOT checked")
    expect_output(show(dots), "takes `\\.\\.\\.`")
})

test_that("a record says when the engine's package is absent", {
    absent <- pecotmr:::.newMethodOptions(
        "notARealPackage::nope",
        defaults = list(a = 1),
        extra = list(),
        label = "nopeArgs",
        engine = "nope"
    )
    expect_output(show(absent), "is not installed")
})

test_that("a zero-length default survives; only NULL is dropped", {
    # compact() drops zero-length values as well as NULLs, which silently lost
    # keepSamples = character(0) and CovFlashOptions(greedy_args = list()).
    # Both are settings, not absent arguments.
    # An EXPLICIT character(0) is a real instruction ("no restriction") and
    # must survive; GenotypeFilterParam() defaults it to NULL, which means
    # "not set" and is correctly dropped.
    a <- GenotypeFilterParam(keepSamples = character(0))
    expect_true(is_in("keepSamples", names(a)))
    expect_identical(a$keepSamples, character(0))
    expect_false(is_in("keepSamples", names(GenotypeFilterParam())))
    expect_true(is_in("greedy_args", names(CovFlashOptions())))
    # A NULL default still means "leave it to the engine".
    expect_false(is_in("subset", names(CovFlashOptions())))
})

# TRUE when the package behind a "pkg::fn" string is installed.
.maPkgAvailable <- function(callee) {
    requireNamespace(sub("::.*$", "", callee), quietly = TRUE)
}

test_that("a constructor is unchecked only when its engine takes dots", {
    # The package-wide invariant behind splitting a method's individual and
    # summary-statistic paths into separate constructors: an argument name
    # goes unchecked ONLY where the engine's own signature ends in `...`, so
    # every name really is legal. Any other unchecked constructor is one
    # whose engine was never enumerated, and a misspelled option there is
    # dropped in silence -- the failure this whole design exists to prevent.
    ns <- asNamespace("pecotmr")
    unchecked <- character(0)
    # *Options() only: a *Param() carries no callee, so `show` prints no
    # check note for it and it could never be counted here anyway.
    for (n in sort(getNamespaceExports("pecotmr"))) {
        if (!str_detect(n, "Options$")) {
            next
        }
        f <- get(n, envir = ns)
        if (!is.function(f)) {
            next
        }
        rec <- tryCatch(f(), error = function(cnd) NULL)
        if (is.null(rec) || !is(rec, "MethodOptions")) {
            next
        }
        out <- paste(capture.output(show(rec)), collapse = " ")
        if (grepl("NOT checked", out)) unchecked <- c(unchecked, n)
    }
    # Each of these reaches an engine whose formals include `...`; where a
    # constructor serves several engines (poolr's three tests), one of them
    # taking dots is enough to make the whole record uncheckable.
    dotEngines <- c(
        CovEdOptions = "mashr::cov_ed",
        DprOptions = "RcppDPR::fit_model",
        GlmnetOptions = "glmnet::cv.glmnet",
        MashCorEmOptions = "mashr::mash_estimate_corr_em",
        MvsusieRssOptions = "mvsusieR::mvsusie_rss",
        NcvregOptions = "ncvreg::cv.ncvreg",
        PoolrOptions = "poolr::fisher",
        QvalueOptions = "qvalue::qvalue",
        RmaOptions = "metafor::rma"
    )
    # Only count the ones whose package is actually installed: an absent
    # package is a third, legitimate reason a record cannot be checked.
    installed <- keep(dotEngines, .maPkgAvailable)
    expect_setequal(intersect(unchecked, names(dotEngines)), names(installed))
    extra <- setdiff(unchecked, names(dotEngines))
    expect_equal(
        extra,
        character(0),
        info = paste(
            "unchecked for a reason other than the engine taking dots:",
            paste(extra, collapse = ", ")
        )
    )
})

test_that("an Options record knows which input classes it is for", {
    # A method on both data paths reaches a different engine on each -- often
    # in a different package -- so a record built for one path is not usable
    # on the other. Recording it lets a pipeline ask the record instead of
    # joining against a registry at the point of use.
    expect_equal(pecotmr:::.optionsInputType(SusieOptions()), "QtlDataset")
    expect_setequal(
        pecotmr:::.optionsInputType(SusieRssOptions()),
        c("QtlSumStats", "GwasSumStats")
    )
    # Cross-package pairs: the whole reason one constructor per method could
    # not carry this.
    expect_equal(pecotmr:::.optionsInputType(GlmnetOptions()), "QtlDataset")
    expect_equal(pecotmr:::.optionsInputType(LassosumOptions()), "QtlSumStats")
    expect_equal(pecotmr:::.optionsInputType(DprOptions()), "QtlDataset")
    expect_equal(pecotmr:::.optionsInputType(SdprOptions()), "QtlSumStats")

    # GwasSumStats is a summary-statistics class the multivariate methods do
    # not accept, which is what gwasAllowed already records -- so the
    # derivation must not hand it to mvsusie.
    expect_setequal(
        pecotmr:::.optionsInputType(SerOptions()),
        c("QtlSumStats", "GwasSumStats")
    )
    expect_equal(
        pecotmr:::.optionsInputType(MvsusieRssOptions()),
        "QtlSumStats"
    )
    expect_false(
        isTRUE(pecotmr:::.fineMappingMethodCapabilities$mvsusie$gwasAllowed)
    )

    # Most engines have no input-path dimension at all: a mash prior
    # component or a p-value test is reached the same way whatever the data
    # looked like.
    expect_null(pecotmr:::.optionsInputType(CovPcaOptions()))
    expect_null(pecotmr:::.optionsInputType(PoolrOptions()))
    expect_null(pecotmr:::.optionsInputType(CtwasOptions()))
    expect_null(pecotmr:::.optionsInputType(UdFitOptions()))
    # And it is reported where a user can see it.
    expect_match(
        paste(capture.output(show(SusieRssOptions())), collapse = " "),
        "for input class"
    )
    expect_no_match(
        paste(capture.output(show(CovPcaOptions())), collapse = " "),
        "for input class"
    )
})

test_that("inputType is derived, not declared per constructor", {
    # The index is built from the two dispatch tables, so a method moving
    # between paths needs no edit here. Checking the shape of that
    # dependency rather than the values: every Options whose callee is a
    # fine-mapping implementation must report that path.
    caps <- pecotmr:::.fineMappingMethodCapabilities
    expect_equal(
        pecotmr:::.fmCalleeInputType(caps$susie$individualImpl),
        "QtlDataset"
    )
    expect_setequal(
        pecotmr:::.fmCalleeInputType(caps$susie$sumstatImpl),
        c("QtlSumStats", "GwasSumStats")
    )
    expect_length(pecotmr:::.fmCalleeInputType("nosuch::fn"), 0L)
    expect_length(pecotmr:::.fmCalleeInputType(NULL), 0L)
    # The TWAS half keys on the constructor name, because its table stores
    # pecotmr implementations rather than callees.
    expect_equal(pecotmr:::.twasLabelInputType("GlmnetOptions"), "QtlDataset")
    expect_equal(
        pecotmr:::.twasLabelInputType("LassosumOptions"),
        "QtlSumStats"
    )
    expect_length(pecotmr:::.twasLabelInputType("NoSuchOptions"), 0L)
})

# ===========================================================================
# Why a name was or was not checked, and the nested-entry helpers
# ===========================================================================

test_that(".engineCheckNote says WHICH reason a name went unchecked", {
    note <- pecotmr:::.engineCheckNote
    # A package that is not installed cannot have its formals read.
    expect_match(
        as.character(note("nosuchpkg::nosuchfn")),
        "package 'nosuchpkg' is not installed"
    )
    # A primitive has no formals at all, which is a different reason.
    expect_match(
        as.character(note("base::sum")),
        "the engine's formals are unknown"
    )
    # An engine taking `...` accepts any name, so checking is meaningless.
    expect_match(
        as.character(note("base::paste")),
        "takes `...`, so any name is accepted"
    )
    # The normal case names the callee it checked against.
    expect_match(
        as.character(note("stats::rnorm")),
        "checked against stats::rnorm"
    )
})

test_that(".twasLabelInputType answers nothing for a label it cannot use", {
    f <- pecotmr:::.twasLabelInputType
    expect_equal(f(42), character(0))
    expect_equal(f(c("a", "b")), character(0))
    expect_equal(f(NULL), character(0))
})

test_that(".nestedExpectedEngine falls back to the key it was given", {
    # The engine is read off a record the constructor builds. When the
    # constructor cannot run, the key is the best name available.
    expect_equal(
        pecotmr:::.nestedExpectedEngine(
            "broken",
            list(broken = function(...) stop("no"))
        ),
        "broken"
    )
})

test_that(".nestedElement refuses an entry that is neither list nor record", {
    expect_error(
        pecotmr:::.nestedElement(42, "pca", list(pca = CovPcaOptions), "lbl"),
        "the `pca` entry must be a list or a constructor result"
    )
    # An already-built record passes straight through.
    rec <- CovPcaOptions()
    expect_identical(
        pecotmr:::.nestedElement(rec, "pca", list(pca = CovPcaOptions), "lbl"),
        rec
    )
})
