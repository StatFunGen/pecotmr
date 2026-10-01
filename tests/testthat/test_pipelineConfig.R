test_that("credibleSetConfig owns L / Lgreedy", {
    # They bound how many credible sets can exist. They reach the engine by a
    # different route from the rest of the bundle -- seeded onto each
    # SuSiE-family token rather than read by postprocessFinemappingFits() --
    # but that is a plumbing difference, not a different setting.
    expect_equal(credibleSetConfig()$L, 10L)
    expect_false("Lgreedy" %in% names(credibleSetConfig()))
    expect_equal(credibleSetConfig(Lgreedy = 3L)$Lgreedy, 3L)
    # addSusieInf stands alone: it selects a chained initialisation between
    # methods rather than describing credible sets.
    expect_false("addSusieInf" %in% names(formals(credibleSetConfig)))
    expect_false("fitStructureArgs" %in% getNamespaceExports("pecotmr"))
})

test_that("credibleSetConfig owns the topLoci per-CS columns", {
    # perCsColumns used to be two booleans on fitRetentionArgs -- `fullFit`
    # and `fullFitAlphaOnly`, the second a documented no-op unless the first.
    # They govern the topLoci TABLE, next to includeAllCs which decides the
    # column labels, not the stored fit.
    expect_equal(credibleSetConfig()$perCsColumns, "none")
    expect_equal(credibleSetConfig(perCsColumns = "full")$perCsColumns, "full")
    expect_error(credibleSetConfig(perCsColumns = "yes"), "must be one of")
    # fitRetentionArgs() is retired: one axis is not a bundle.
    expect_false("fitRetentionArgs" %in% getNamespaceExports("pecotmr"))
})

test_that("signalScreenConfig resolves each metric to a pipeline spec", {
    # pip keeps the legacy bare-numeric spelling; the other three resolve to
    # list(metric, cutoff). .asScreen() canonicalizes either.
    expect_null(pecotmr:::.screenResolve(signalScreenConfig()))
    expect_equal(pecotmr:::.screenResolve(signalScreenConfig(pip = 0.5)), 0.5)
    expect_equal(
        pecotmr:::.screenResolve(signalScreenConfig(absZ = 5)),
        list(metric = "absZ", cutoff = 5)
    )
    expect_equal(
        pecotmr:::.screenResolve(signalScreenConfig(bf = 100)),
        list(metric = "bf", cutoff = 100)
    )
    expect_equal(
        pecotmr:::.screenResolve(signalScreenConfig(logBf = 3)),
        list(metric = "logBf", cutoff = 3)
    )
})

test_that("signalScreenConfig rejects two metrics and negative absZ / bf", {
    expect_error(
        signalScreenConfig(pip = 0.5, absZ = 5),
        "only one screening metric"
    )
    # 0 is the long-standing "off" spelling, so it is not a second metric.
    expect_equal(
        pecotmr:::.screenResolve(signalScreenConfig(pip = 0, absZ = 5)),
        list(metric = "absZ", cutoff = 5)
    )
    expect_error(signalScreenConfig(absZ = -1), "must be > 0")
    expect_error(signalScreenConfig(bf = -1), "must be > 0")
    # pip < 0 is the adaptive 3 / nVariants convention, and logBf is a log.
    expect_equal(pecotmr:::.screenResolve(signalScreenConfig(pip = -1)), -1)
    expect_equal(
        pecotmr:::.screenResolve(signalScreenConfig(logBf = -2)),
        list(metric = "logBf", cutoff = -2)
    )
})

test_that("sumstats inputs ignore residualization but refuse crossValidation", {
    # The asymmetry follows from the defaults, not from taste. CV is off
    # unless asked for (folds = 0), so a non-default value on a
    # summary-statistics run is an explicit request for something the input
    # cannot do -> error. Residualization is ON by default, so refusing a
    # non-default would reject the DEFAULT bundle and force every sumstats
    # caller to unset it -> ignore.
    expect_equal(crossValidationConfig()$folds, 0)
    expect_true(residualizationConfig()$residualizePhenotype)
    expect_true(residualizationConfig()$residualizeGenotype)
    # CV is refused on a summary-statistics input ...
    expect_error(
        pecotmr:::.cvRefuseOnSumstats(
            crossValidationConfig(folds = 5),
            "fineMappingPipeline",
            "QtlSumStats"
        ),
        "cross-validation"
    )
    # ... while the DEFAULT CV bundle passes, so one bundle still travels to
    # either input kind.
    expect_silent(pecotmr:::.cvRefuseOnSumstats(
        crossValidationConfig(),
        "fineMappingPipeline",
        "QtlSumStats"
    ))
    # Residualization has no equivalent refusal: the sumstats workers do not
    # take it at all, so it cannot reach an engine and cannot error.
    for (fn in c(
        ".fmPipelineQtlSumStats",
        ".fmPipelineGwas",
        ".twasPipelineQtlSumStats"
    )) {
        f <- get(fn, envir = asNamespace("pecotmr"))
        expect_false("residualization" %in% names(formals(f)), label = fn)
    }
})

# =============================================================================
# Filtering bundles (folded in from filterConfig.R)
# =============================================================================

test_that("the three filter bundles are separate on purpose", {
    # imissCutoff genuinely differs: it resolves to 0 on the genotype path
    # (admit no missingness) and 1 on the panel path (filter nothing, so the
    # allele-frequency sidecar can be read instead of materializing dosage).
    expect_equal(
        pecotmr:::.qtlResolveFilter(genotypeFilterConfig())$imissCutoff,
        0
    )
    expect_equal(panelFilterConfig()$imissCutoff, 1)
})

test_that("an unset genotype filter field is absent, not NULL", {
    # Absence is what distinguishes "pin this value" from "leave it alone":
    # QtlDataset resolves an unset field to its own default, a pipeline
    # leaves the dataset's construct-time slot untouched.
    expect_length(genotypeFilterConfig(), 0L)
    expect_equal(names(genotypeFilterConfig(mafCutoff = 0)), "mafCutoff")
    # A pinned zero is NOT the same as unset.
    expect_equal(genotypeFilterConfig(mafCutoff = 0)$mafCutoff, 0)
    expect_null(genotypeFilterConfig()$mafCutoff)
})

test_that("every genotypeFilterConfig field is one QtlDataset resolves", {
    expect_setequal(
        names(formals(genotypeFilterConfig)),
        names(pecotmr:::.qtlResolveFilter(genotypeFilterConfig()))
    )
})

test_that("each bundle carries only fields its consumers read", {
    expect_setequal(
        names(panelFilterConfig()),
        c("mafCutoff", "macCutoff", "imissCutoff")
    )
    expect_setequal(
        names(sumstatsFilterConfig()),
        c("removeIndels", "removeStrandAmbiguous", "infoCutoff", "nCutoff")
    )
})

test_that("a field belonging to another bundle is refused", {
    # The whole point of three constructors: a sumstats-row option passed to
    # the genotype filter would be silently ignored if they were one.
    expect_error(genotypeFilterConfig(removeIndels = TRUE), "unused argument")
    expect_error(panelFilterConfig(xvarCutoff = 0.1), "unused argument")
    expect_error(sumstatsFilterConfig(mafCutoff = 0.01), "unused argument")
})

test_that("all three are MethodConfig and splice like a list", {
    bundles <- list(
        genotypeFilterConfig(),
        panelFilterConfig(),
        sumstatsFilterConfig()
    )
    for (a in bundles) {
        expect_s4_class(a, "MethodConfig")
        expect_type(as.list(a), "list")
    }
})

test_that(".panelCutoffs answers NULL for a filter that keeps everything", {
    expect_null(pecotmr:::.panelCutoffs(panelFilterConfig()))
    expect_equal(
        pecotmr:::.panelCutoffs(panelFilterConfig(mafCutoff = 0.01))$mafCutoff,
        0.01
    )
})

test_that("sumstatsCleaningConfig is separate from sumstatsFilterConfig", {
    # Cleaning decides whether a row is a well-formed record; filtering
    # applies quality thresholds to rows that already are. No field overlaps.
    expect_length(
        intersect(
            names(formals(sumstatsCleaningConfig)),
            names(formals(sumstatsFilterConfig))
        ),
        0L
    )
})

test_that("every sumstatsCleaningConfig field is one the resolver hands on", {
    expect_setequal(
        names(formals(sumstatsCleaningConfig)),
        names(pecotmr:::.sumstatsCleaningResolve(sumstatsCleaningConfig()))
    )
    # And every one of them is an .applySanityChecks argument, since the
    # resolved record is spliced straight into it.
    expect_true(all(
        names(pecotmr:::.sumstatsCleaningResolve(sumstatsCleaningConfig())) %in%
            names(formals(pecotmr:::.applySanityChecks))
    ))
})

test_that(".sumstatsCleaningResolve fills defaults an empty bundle omits", {
    # .applySanityChecks tests these with `if (!flag)`, so a missing field
    # must not arrive as NULL.
    res <- pecotmr:::.sumstatsCleaningResolve(list())
    expect_equal(
        res,
        pecotmr:::.sumstatsCleaningResolve(sumstatsCleaningConfig())
    )
    expect_true(res$coerceNumeric)
    expect_equal(res$smallPFloor, 5e-324)
    expect_error(
        pecotmr:::.sumstatsCleaningResolve(list(coerceNumeric = FALSE)),
        "must be built with sumstatsCleaningConfig"
    )
})
