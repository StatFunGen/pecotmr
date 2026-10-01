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

test_that(".newMethodConfig records NA for a callee given as a function", {
    res <- pecotmr:::.newMethodConfig(
        stats::var,
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs"
    )
    expect_true(is.na(S4Vectors::metadata(res)$callee))
    # Validation still happens -- the formals come from the function itself.
    expect_error(
        pecotmr:::.newMethodConfig(
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

test_that(".newMethodConfig merges defaults under extras, dropping NULLs", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(na.rm = TRUE, use = NULL),
        extra = list(y = 1),
        label = "varArgs"
    )
    expect_true(pecotmr:::.isMethodConfig(res))
    # `use = NULL` means "leave it to the engine", so it is not forwarded.
    expect_false(is_in("use", names(res)))
    expect_true(res$na.rm)
    expect_equal(res$y, 1)
})

test_that(".newMethodConfig rejects a name the engine cannot accept", {
    expect_error(
        pecotmr:::.newMethodConfig(
            "stats::var",
            defaults = list(),
            extra = list(naRM = TRUE),
            label = "varArgs"
        ),
        "unknown argument\\(s\\) naRM"
    )
})

test_that(".newMethodConfig names the engine's accepted arguments in the error", {
    expect_error(
        pecotmr:::.newMethodConfig(
            "stats::var",
            defaults = list(),
            extra = list(nope = 1),
            label = "varArgs"
        ),
        "stats::var. accepts:.*na\\.rm"
    )
})

test_that(".newMethodConfig skips validation when the engine takes dots", {
    # base::paste has `...`, so an arbitrary name is legitimate.
    res <- pecotmr:::.newMethodConfig(
        "base::paste",
        defaults = list(sep = "-"),
        extra = list(anythingAtAll = 1),
        label = "pasteArgs"
    )
    expect_equal(res$anythingAtAll, 1)
})

test_that(".newMethodConfig skips validation when the package is absent", {
    res <- pecotmr:::.newMethodConfig(
        "notARealPackage::nope",
        defaults = list(),
        extra = list(whatever = 2),
        label = "nopeArgs"
    )
    expect_equal(res$whatever, 2)
})

test_that(".newMethodConfig requires every extra to be named", {
    expect_error(
        pecotmr:::.newMethodConfig(
            "base::paste",
            defaults = list(),
            extra = list(1, 2),
            label = "pasteArgs"
        ),
        "must be named"
    )
})

test_that(".newMethodConfig rejects an argument supplied twice", {
    # Two entries in `extra` -- what .methodConfigWith() produces when a
    # wrapper's own contribution meets the caller's record -- have no winner.
    expect_error(
        pecotmr:::.newMethodConfig(
            "base::paste",
            defaults = list(),
            extra = list(sep = "-", sep = "+"),
            label = "pasteArgs"
        ),
        "given twice: sep"
    )
})

test_that(".newMethodConfig lets an explicit value override a default", {
    res <- pecotmr:::.newMethodConfig(
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

test_that(".newMethodConfig records the engine and callee for reporting", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_equal(S4Vectors::metadata(res)$engine, "var")
    expect_equal(S4Vectors::metadata(res)$callee, "stats::var")
})

test_that("MethodConfig keeps its class through pecotmr's list operations", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(na.rm = TRUE, y = 1),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    # .fmNarrowAfterJoint / .fmQssJointPhase subset these records by name, and
    # a plain list carrying an S3 class loses the class right here.
    expect_true(pecotmr:::.isMethodConfig(res["na.rm"]))
    expect_true(pecotmr:::.isMethodConfig(res[1]))
    expect_equal(S4Vectors::metadata(res["na.rm"])$engine, "var")
})

test_that(".newMethodConfig splices into a call exactly like a bare list", {
    # The whole point of returning a plain list: every existing
    # exec(fn, !!!args) call site keeps working unchanged.
    res <- pecotmr:::.newMethodConfig(
        "base::paste",
        defaults = list(sep = "-"),
        extra = list(),
        label = "pasteArgs"
    )
    expect_equal(exec(paste, "a", "b", !!!res), "a-b")
})

test_that(".assertMethodConfig accepts a constructed record and empty list", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(),
        extra = list(),
        label = "varArgs"
    )
    expect_silent(pecotmr:::.assertMethodConfig(res, "varArgs", "methodArgs"))
    # `list()` carries no settings, so there is nothing a constructor would
    # have added and nothing to validate.
    expect_silent(pecotmr:::.assertMethodConfig(
        list(),
        "varArgs",
        "methodArgs"
    ))
})

test_that(".assertMethodConfig rejects a populated bare list", {
    expect_error(
        pecotmr:::.assertMethodConfig(
            list(na.rm = TRUE),
            "varArgs",
            "methodArgs"
        ),
        "`methodArgs` must be built with varArgs\\(\\)"
    )
})

test_that("show(MethodConfig) reports the engine and its settings", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_output(show(res), "MethodConfig: var")
    expect_output(show(res), "forwarded to stats::var")
    expect_output(show(res), "na\\.rm")
})

test_that("show(MethodConfig) says so when nothing is set", {
    res <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(),
        extra = list(),
        label = "varArgs"
    )
    expect_output(show(res), "no options set")
})

# --- character-or-constructor selection ------------------------------------

.ma_fooArgs <- function(a = 1, b = 2) {
    pecotmr:::.newMethodConfig(
        NULL,
        defaults = list(a = a, b = b),
        extra = list(),
        label = "fooArgs",
        engine = "foo"
    )
}

.ma_barArgs <- function(z = 9) {
    pecotmr:::.newMethodConfig(
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
    expect_true(pecotmr:::.isMethodConfig(res$args))
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

test_that(".newNestedConfig accepts plain lists and records interchangeably", {
    res <- pecotmr:::.newNestedConfig(
        list(foo = list(a = 5), bar = .ma_barArgs(z = 1)),
        .ma_ctors(),
        "demoArgs"
    )
    expect_true(pecotmr:::.isMethodConfig(res))
    expect_equal(res$foo$a, 5)
    expect_equal(res$bar$z, 1)
    # A plain list is the same validation written more briefly: it is spliced
    # into the constructor, so the defaults apply either way.
    expect_equal(res$foo$b, 2)
})

test_that(".newNestedConfig validates a plain list through its constructor", {
    expect_error(
        pecotmr:::.newNestedConfig(
            list(foo = list(aa = 5)),
            .ma_ctors(),
            "demoArgs"
        ),
        "unused argument"
    )
})

test_that(".newNestedConfig refuses an entry built by the wrong constructor", {
    expect_error(
        pecotmr:::.newNestedConfig(
            list(foo = .ma_barArgs()),
            .ma_ctors(),
            "demoArgs"
        ),
        "was built with the constructor for 'bar'"
    )
})

test_that(".newNestedConfig rejects unknown and unnamed entries", {
    expect_error(
        pecotmr:::.newNestedConfig(
            list(nope = list()),
            .ma_ctors(),
            "demoArgs"
        ),
        "unknown method\\(s\\) nope"
    )
    expect_error(
        pecotmr:::.newNestedConfig(list(list()), .ma_ctors(), "demoArgs"),
        "must be named for its method"
    )
})

test_that(".newNestedConfig rejects a non-list entry", {
    expect_error(
        pecotmr:::.newNestedConfig(list(foo = 5), .ma_ctors(), "demoArgs"),
        "must be a list or a constructor result"
    )
})

test_that(".newNestedConfig with no entries is an empty record", {
    res <- pecotmr:::.newNestedConfig(list(), .ma_ctors(), "demoArgs")
    expect_true(pecotmr:::.isMethodConfig(res))
    expect_length(res, 0L)
})

test_that("a record says whether its names were checked", {
    checked <- pecotmr:::.newMethodConfig(
        "stats::var",
        defaults = list(na.rm = TRUE),
        extra = list(),
        label = "varArgs",
        engine = "var"
    )
    expect_output(show(checked), "argument names checked against stats::var")
    # An engine taking `...` accepts anything, so a misspelling reaches it.
    dots <- pecotmr:::.newMethodConfig(
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
    absent <- pecotmr:::.newMethodConfig(
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
    # keepSamples = character(0) and covFlashConfig(greedy_args = list()).
    # Both are settings, not absent arguments.
    # An EXPLICIT character(0) is a real instruction ("no restriction") and
    # must survive; genotypeFilterConfig() defaults it to NULL, which means
    # "not set" and is correctly dropped.
    a <- genotypeFilterConfig(keepSamples = character(0))
    expect_true(is_in("keepSamples", names(a)))
    expect_identical(a$keepSamples, character(0))
    expect_false(is_in("keepSamples", names(genotypeFilterConfig())))
    expect_true(is_in("greedy_args", names(covFlashConfig())))
    # A NULL default still means "leave it to the engine".
    expect_false(is_in("subset", names(covFlashConfig())))
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
    for (n in sort(getNamespaceExports("pecotmr"))) {
        if (!grepl("Config$", n)) {
            next
        }
        f <- get(n, envir = ns)
        if (!is.function(f)) {
            next
        }
        rec <- tryCatch(f(), error = function(cnd) NULL)
        if (is.null(rec) || !is(rec, "MethodConfig")) {
            next
        }
        out <- paste(capture.output(show(rec)), collapse = " ")
        if (grepl("NOT checked", out)) unchecked <- c(unchecked, n)
    }
    # Each of these reaches an engine whose formals include `...`; where a
    # constructor serves several engines (poolr's three tests), one of them
    # taking dots is enough to make the whole record uncheckable.
    dotEngines <- c(
        covEdConfig = "mashr::cov_ed",
        dprConfig = "RcppDPR::fit_model",
        glmnetConfig = "glmnet::cv.glmnet",
        mashCorEmConfig = "mashr::mash_estimate_corr_em",
        mvsusieRssConfig = "mvsusieR::mvsusie_rss",
        ncvregConfig = "ncvreg::cv.ncvreg",
        poolrConfig = "poolr::fisher",
        qvalueConfig = "qvalue::qvalue",
        rmaConfig = "metafor::rma"
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
