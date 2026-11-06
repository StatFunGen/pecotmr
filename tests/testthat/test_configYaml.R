# Settings records to and from YAML, R/configYaml.R.

# The message a writer's refusal carries, its causes included.
.cfg_refusal <- function(obj, path) {
    conditionMessage(tryCatch(writeConfigToYaml(obj, path), error = identity))
}

test_that("a record loads from its constructor's mapping", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(
        c(
            "CredibleSetParam:",
            "  coverage: 0.9",
            "  perCsColumns: alpha",
            "SusieOptions:",
            "  max_iter: 200"
        ),
        path
    )
    expect_identical(
        loadConfigFromYaml(path, "CredibleSetParam"),
        CredibleSetParam(coverage = 0.9, perCsColumns = "alpha")
    )
    expect_identical(
        loadConfigFromYaml(path, "SusieOptions"),
        SusieOptions(max_iter = 200L)
    )
})

test_that("a constructor with nothing set loads its defaults", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(c("PanelFilterParam:", "SusieOptions: {}"), path)
    expect_identical(
        loadConfigFromYaml(path, "PanelFilterParam"),
        PanelFilterParam()
    )
    expect_identical(loadConfigFromYaml(path, "SusieOptions"), SusieOptions())
})

test_that("scalars follow the YAML 1.2 core schema", {
    # Under 1.1 the key `n` reads as FALSE, `yes` as TRUE, and 1e-5 as a
    # string because a 1.1 float needs a dot.
    parsed <- pecotmr:::.configParse(str_c(
        "n: 5\ny: yes\nflag: true\np: 1e-5\nbig: 3000000000\n",
        "v: [1, 2.5]\ntiny: 4.94065645841247e-324\n"
    ))
    expect_named(parsed, c("n", "y", "flag", "p", "big", "v", "tiny"))
    expect_identical(parsed$y, "yes")
    expect_true(parsed$flag)
    expect_identical(parsed$p, 1e-5)
    expect_identical(parsed$big, 3e9)
    expect_identical(parsed$v, c(1, 2.5))
    expect_identical(parsed$tiny, 4.94065645841247e-324)
})

test_that("an empty sequence takes the type of the slot it fills", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(
        c(
            "GenotypeFilterParam:",
            "  keepSamples: []",
            "MashComponentParam:",
            "  pca: []"
        ),
        path
    )
    expect_identical(
        loadConfigFromYaml(path, "GenotypeFilterParam")$keepSamples,
        character(0)
    )
    # A slot that holds a record takes an empty list, which builds the
    # record with nothing set.
    expect_identical(
        loadConfigFromYaml(path, "MashComponentParam"),
        MashComponentParam(pca = CovPcaOptions())
    )
})

test_that("a nested mapping builds the record its default names", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(
        c(
            "GwasFineMappingParam:",
            "  credibleSetParam:",
            "    coverage: 0.9",
            "CovUdrOptions:",
            "  fit:",
            "    maxiter: 50"
        ),
        path
    )
    expect_identical(
        loadConfigFromYaml(path, "GwasFineMappingParam")$credibleSetParam,
        CredibleSetParam(coverage = 0.9)
    )
    expect_identical(
        loadConfigFromYaml(path, "CovUdrOptions"),
        CovUdrOptions(fit = UdFitOptions(maxiter = 50L))
    )
})

test_that("a mapping keyed by a constructor builds that record anywhere", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(
        c(
            "MashPriorParam:",
            "  engine:",
            "    CovUdrOptions: {}",
            "FineMappingMethodsParam:",
            "  qtlSumStatsMethods:",
            "    susie:",
            "      SusieRssOptions:",
            "        max_iter: 50"
        ),
        path
    )
    expect_identical(
        loadConfigFromYaml(path, "MashPriorParam"),
        MashPriorParam(engine = CovUdrOptions())
    )
    expect_identical(
        loadConfigFromYaml(path, "FineMappingMethodsParam"),
        FineMappingMethodsParam(
            qtlSumStatsMethods = list(susie = SusieRssOptions(max_iter = 50L))
        )
    )
})

test_that("loading names what is wrong", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeLines(c("CredibleSetParam:", "  coverage: 2"), path)
    expect_error(
        loadConfigFromYaml(path, "SusieOptions"),
        "has no `SusieOptions` record. It has: CredibleSetParam"
    )
    expect_error(
        loadConfigFromYaml(path, "qcOptions"),
        "not a pecotmr \\*Options\\(\\) or \\*Param\\(\\) constructor"
    )
    expect_error(
        loadConfigFromYaml(path, "CredibleSetParam"),
        "Could not build `CredibleSetParam`"
    )
    # A record inside another is named by the setting that holds it.
    writeLines(
        c("GwasFineMappingParam:", "  credibleSetParam:", "    coverage: 2"),
        path
    )
    nested <- conditionMessage(tryCatch(
        loadConfigFromYaml(path, "GwasFineMappingParam"),
        error = identity
    ))
    expect_match(nested, "Could not build `GwasFineMappingParam`")
    expect_match(nested, "In setting `credibleSetParam`.")
    expect_no_match(nested, "In index")
    writeLines(c("SusieOptions:", "  - 1", "  - 2"), path)
    expect_error(
        loadConfigFromYaml(path, "SusieOptions"),
        "must be a mapping of settings to values"
    )
    writeLines("- CredibleSetParam", path)
    expect_error(
        loadConfigFromYaml(path, "CredibleSetParam"),
        "must be a YAML mapping keyed by constructor name"
    )
})

test_that("every constructor's defaults survive a round trip", {
    path <- withr::local_tempfile(fileext = ".yaml")
    ctors <- stringr::str_subset(
        getNamespaceExports("pecotmr"),
        "^[A-Z][A-Za-z0-9]*(Options|Param)$"
    )
    records <- map(ctors, exec)
    expect_identical(writeConfigToYaml(records, path), path)
    # BPPARAM is not written, so the rebuilt record carries a fresh backend.
    backTo <- map(ctors, loadConfigFromYaml, filename = path)
    expect_identical(
        map2(backTo, records, pecotmr:::.configCarryBackends),
        records
    )
})

test_that("set values, nested records and awkward numbers round-trip", {
    path <- withr::local_tempfile(fileext = ".yaml")
    records <- list(
        GwasFineMappingParam(
            credibleSetParam = CredibleSetParam(coverage = 0.9)
        ),
        MashPriorParam(
            engine = CovUdrOptions(fit = UdFitOptions(maxiter = 50))
        ),
        GenotypeFilterParam(keepSamples = character(0), mafCutoff = 1e-5),
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = SusieOptions(max_iter = 200))
        ),
        # A method selected by name has nothing set: `susie: []`.
        TwasWeightsMethodsParam(methods = c("susie", "lasso")),
        PvalMethodsParam(aspu = AspuOptions(n.perm = 5000)),
        MashComponentParam(pca = list(subset = 1:5)),
        # No 15-17 digit decimal of this one reads back exactly in R, so it
        # is written as a hex float.
        CredibleSetParam(coverage = 0.123456789012345678),
        QtlEnrichmentParam(lambda = Inf, numGwas = 3e9),
        # NaN and NA_real_ each read back as themselves, not as the other.
        LdPruningOptions(maf = NaN, slide.max.n = NA_real_)
    )
    writeConfigToYaml(records, path)
    expect_match(read_file(path), "coverage: 0x", fixed = TRUE)
    ctorNames <- map_chr(records, pecotmr:::.configCtorName)
    expect_identical(
        map(ctorNames, loadConfigFromYaml, filename = path),
        records
    )
})

test_that("BPPARAM is left out of the file", {
    path <- withr::local_tempfile(fileext = ".yaml")
    writeConfigToYaml(
        CrossValidationParam(BPPARAM = BiocParallel::SnowParam(2)),
        path
    )
    expect_no_match(read_file(path), "BPPARAM")
})

test_that("writing refuses a record it cannot write, naming where", {
    path <- withr::local_tempfile(fileext = ".yaml")
    expect_error(
        writeConfigToYaml(list(CredibleSetParam(), CredibleSetParam()), path),
        "one record per constructor; given more than one of: CredibleSetParam"
    )
    expect_error(
        writeConfigToYaml(list(CredibleSetParam(), list()), path),
        "must be a \\*Options\\(\\) / \\*Param\\(\\) record or a list of them"
    )
    # The record and the setting are named above the cause, without purrr's
    # framing for each element it mapped over.
    named <- .cfg_refusal(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(
                susie = SusieOptions(prior_weights = c(a = 1))
            )
        ),
        path
    )
    expect_match(named, "Could not write `FineMappingMethodsParam` to")
    expect_match(
        named,
        "In setting `qtlDatasetMethods$susie$prior_weights`.",
        fixed = TRUE
    )
    expect_match(named, "holds a named numeric, which has no YAML form")
    expect_no_match(named, "In index")
    expect_match(
        .cfg_refusal(MashPriorParam(priorCovariances = list(diag(2))), path),
        "In setting `priorCovariances[[1]]`.",
        fixed = TRUE
    )
    # One that names no constructor is named by its position.
    stripped <- SusieOptions()
    S4Vectors::metadata(stripped)$constructor <- NULL
    unnamed <- .cfg_refusal(list(CredibleSetParam(), stripped), path)
    expect_match(unnamed, "Could not write record 2 to")
    expect_match(unnamed, "does not say which constructor built it")
    expect_false(file.exists(path))
})

test_that("writing refuses what would not read back", {
    path <- withr::local_tempfile(fileext = ".yaml")
    # A string that reads as an exponent number cannot be written. Where the
    # constructor takes the number, the settings that would change are named;
    # where it refuses the number, its own error is the cause.
    changed <- .cfg_refusal(LdPruningOptions(sample.id = c("1e5", "S2")), path)
    expect_match(changed, "Could not write `LdPruningOptions` to")
    expect_match(changed, "It would not read back as written.")
    expect_match(changed, "These settings would change: sample.id.")
    expect_match(
        .cfg_refusal(
            ResidualizationParam(phenotypeCovariates = c("age", "1e5")),
            path
        ),
        "These settings would change: phenotypeCovariates."
    )
    refused <- .cfg_refusal(RelatednessColumnsParam(iid1 = "1e5"), path)
    expect_match(refused, "It would not read back as written.")
    expect_match(refused, "invalid object for slot \"iid1\"")
    # metadata() is not written: the constructor sets it again.
    noted <- SusieOptions(max_iter = 5L)
    S4Vectors::metadata(noted)$note <- "x"
    expect_match(
        .cfg_refusal(noted, path),
        "its metadata() would come back as its constructor sets it",
        fixed = TRUE
    )
    expect_false(file.exists(path))
})
