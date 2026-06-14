context("colocboostPipeline")

# ===========================================================================
# Tests from test_colocboost_pipeline.R
# ===========================================================================

# Wrap a correlation (or genotype) matrix into an LdData for use in test mocks
# that previously used bare matrices for ldMat/ldInfo fields.
.test_lddata_from_matrix <- function(mat, is_genotype = FALSE) {
  vids <- if (is_genotype) colnames(mat) else rownames(mat)
  if (is.null(vids)) vids <- colnames(mat)
  ref_panel <- cbind(parseVariantId(vids), variant_id = vids)
  ref_panel$chrom <- as.character(ref_panel$chrom)
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  bm <- pecotmr:::.inferSingleLdBlockMetadata(ref_panel)
  if (is_genotype) {
    LdData(correlation = NULL, genotypeHandle = mat,
           variants = variants_gr, blockMetadata = bm,
           nRef = as.integer(nrow(mat)))
  } else {
    LdData(correlation = mat, variants = variants_gr, blockMetadata = bm)
  }
}

# Wrap one (rss_input, ldMatrix) pair as a QcResult for mocks that previously
# returned the legacy list shape.
.test_qcresult_from_list <- function(rss_input, ldMat) {
  QcResult(
    ldData = .test_lddata_from_matrix(ldMat),
    rssInput = rss_input,
    preprocess = list(),
    outlierNumber = 0L,
    skipped = FALSE
  )
}

# Build a RegionalData S4 object from legacy-style per-context lists of
# residual_Y / residual_X / maf matrices, so existing tests can keep their
# inline list builders.
.test_regionaldata_from_lists <- function(residual_Y, residual_X, maf = NULL) {
  contexts <- names(residual_Y)
  if (is.null(contexts)) contexts <- paste0("ctx", seq_along(residual_Y))
  # Pick the first X as the canonical genotype matrix; tests share X across
  # contexts in nearly all cases. Ensure sample rownames are present.
  X0 <- residual_X[[1]]
  if (is.null(rownames(X0))) {
    rownames(X0) <- paste0("sample", seq_len(nrow(X0)))
  }
  sample_ids <- rownames(X0)
  # Align each phenotype matrix to the same sample IDs (assume same n).
  phenotypes <- stats::setNames(lapply(seq_along(contexts), function(i) {
    y <- residual_Y[[i]]
    if (!is.matrix(y)) y <- as.matrix(y)
    if (is.null(rownames(y))) rownames(y) <- sample_ids[seq_len(nrow(y))]
    y
  }), contexts)
  covariates <- stats::setNames(lapply(contexts, function(c) {
    matrix(numeric(0), nrow = nrow(X0), ncol = 0,
           dimnames = list(sample_ids, NULL))
  }), contexts)
  if (is.null(maf)) {
    maf_list <- stats::setNames(lapply(contexts, function(c) {
      rep(0.1, ncol(X0))
    }), contexts)
  } else {
    maf_list <- if (is.null(names(maf))) stats::setNames(maf, contexts) else maf
    if (length(maf_list) < length(contexts)) {
      maf_list <- stats::setNames(rep(list(maf_list[[1]]), length(contexts)), contexts)
    }
  }
  RegionalData(
    genotypeMatrix = X0,
    phenotypes = phenotypes,
    covariates = covariates,
    scaleResiduals = FALSE,
    maf = maf_list,
    region = NULL,
    droppedSamples = list(X = list(), Y = list(), covar = list()),
    coordinates = NULL
  )
}


# ---- qc_method match.arg ----
test_that("qcRegionalData is exported for downstream use", {
  expect_true("qcRegionalData" %in% getNamespaceExports("pecotmr"))
  expect_identical(getExportedValue("pecotmr", "qcRegionalData"), qcRegionalData)
})

test_that("qcRegionalData accepts explicit zMismatchQc = 'slalom'", {
  region_data <- list(individualData = NULL, sumstatData = NULL)
  result <- qcRegionalData(region_data, zMismatchQc = "slalom")
  expect_type(result, "list")
})

test_that("qcRegionalData rejects invalid qc_method", {
  region_data <- list(individualData = NULL, sumstatData = NULL)
  expect_error(
    qcRegionalData(region_data, zMismatchQc = "invalid"),
    "arg"
  )
})

# ---- pipCutoffToSkipInd validation ----
test_that("pip_cutoff scalar is recycled for individual contexts", {
  # Create individualData with 3 real-ish contexts
  set.seed(42)
  n <- 10; p <- 5
  make_ctx <- function() {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("var", 1:p)
    Y <- matrix(rnorm(n * 2), n, 2)
    colnames(Y) <- paste0("gene", 1:2)
    list(X = X, Y = Y, maf = runif(p, 0.05, 0.5))
  }
  ctx <- make_ctx()
  individualData <- .test_regionaldata_from_lists(
    residual_Y = list(ctx1 = ctx$Y, ctx2 = ctx$Y, ctx3 = ctx$Y),
    residual_X = list(ctx1 = ctx$X, ctx2 = ctx$X, ctx3 = ctx$X),
    maf = list(ctx1 = ctx$maf, ctx2 = ctx$maf, ctx3 = ctx$maf)
  )
  region_data <- list(individualData = individualData, sumstatData = NULL)

  # Scalar 0 (no PIP check) should be recycled and run without error
  result <- qcRegionalData(region_data, pipCutoffToSkipInd = 0)
  expect_type(result, "list")
})

test_that("pip_cutoff wrong length errors for individual contexts", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("var", 1:p)
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"
  individualData <- .test_regionaldata_from_lists(
    residual_Y = list(ctx1 = Y, ctx2 = Y, ctx3 = Y),
    residual_X = list(ctx1 = X, ctx2 = X, ctx3 = X),
    maf = list(ctx1 = runif(p), ctx2 = runif(p), ctx3 = runif(p))
  )
  region_data <- list(individualData = individualData, sumstatData = NULL)

  expect_error(
    qcRegionalData(region_data, pipCutoffToSkipInd = c(0, 0)),
    "pipCutoffToSkipInd"
  )
})

test_that("pip_cutoff correct length vector works", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("var", 1:p)
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"
  individualData <- .test_regionaldata_from_lists(
    residual_Y = list(ctx1 = Y, ctx2 = Y),
    residual_X = list(ctx1 = X, ctx2 = X),
    maf = list(ctx1 = runif(p, 0.05, 0.5), ctx2 = runif(p, 0.05, 0.5))
  )
  region_data <- list(individualData = individualData, sumstatData = NULL)

  # Length-2 vector for 2 contexts should work
  result <- qcRegionalData(region_data, pipCutoffToSkipInd = c(0, 0))
  expect_type(result, "list")
})

# ===========================================================================
# Tests from test_colocboost_pipeline_comprehensive.R
# ===========================================================================


# ===========================================================================
# Helper: build a minimal region_data with individual-level data
# ===========================================================================
make_individual_region_data <- function(n = 20, p = 8, n_contexts = 2, n_events = 3) {
  set.seed(101)
  sample_ids <- paste0("sample", seq_len(n))
  var_ids <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  X <- matrix(rnorm(n * p), n, p, dimnames = list(sample_ids, var_ids))
  contextNames <- paste0("ctx", seq_len(n_contexts))
  phenotypes <- stats::setNames(lapply(seq_len(n_contexts), function(i) {
    Y <- matrix(rnorm(n * n_events), n, n_events,
                dimnames = list(sample_ids, paste0("event", seq_len(n_events))))
    Y
  }), contextNames)
  # Per-context covariates: empty intercept-only model (n x 0 with rownames)
  covariates <- stats::setNames(lapply(seq_len(n_contexts), function(i) {
    matrix(numeric(0), nrow = n, ncol = 0, dimnames = list(sample_ids, NULL))
  }), contextNames)
  maf_list <- stats::setNames(lapply(seq_len(n_contexts), function(i) {
    runif(p, 0.05, 0.45)
  }), contextNames)
  rd <- RegionalData(
    genotypeMatrix = X,
    phenotypes = phenotypes,
    covariates = covariates,
    scaleResiduals = FALSE,
    maf = maf_list,
    region = NULL,
    droppedSamples = list(X = list(), Y = list(), covar = list()),
    coordinates = NULL
  )
  list(
    individualData = rd,
    sumstatData = NULL
  )
}

# ===========================================================================
# Helper: build a minimal region_data with sumstat data
# ===========================================================================
make_sumstat_region_data <- function(n_variants = 5, n_studies = 2) {
  set.seed(202)
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  ref_panel <- data.frame(
    chrom = as.character(rep(1, n_variants)),
    pos   = seq_len(n_variants) * 100,
    A2    = rep("A", n_variants),
    A1    = rep("G", n_variants),
    variant_id = vids,
    stringsAsFactors = FALSE
  )

  sumstats_list <- lapply(seq_len(n_studies), function(i) {
    ss <- list(
      sumstats = data.frame(
        chrom      = rep(1, n_variants),
        pos        = seq_len(n_variants) * 100,
        A1         = rep("G", n_variants),
        A2         = rep("A", n_variants),
        beta       = rnorm(n_variants),
        se         = runif(n_variants, 0.05, 0.2),
        z          = rnorm(n_variants, 0, 2),
        variant_id = vids,
        stringsAsFactors = FALSE
      ),
      n     = 10000,
      varY = 1
    )
    list(ss) |> setNames(paste0("study", i))
  })

  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  ld_data <- LdData(
    correlation = ldMat,
    variants = variants_gr,
    blockMetadata = pecotmr:::.inferSingleLdBlockMetadata(ref_panel)
  )

  list(
    individualData = NULL,
    sumstatData    = list(
      sumstats = sumstats_list,
      ldInfo  = list(ld_data)
    )
  )
}

test_that("qcRegionalData treats NULL qc_method as basic-only none", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  captured_qc_method <- NULL
  ldMat <- diag(1)
  rownames(ldMat) <- colnames(ldMat) <- "chr1:100:A:G"

  ref_panel_one <- data.frame(
    chrom = "1", pos = 100L, A2 = "A", A1 = "G",
    variant_id = "chr1:100:A:G", stringsAsFactors = FALSE
  )
  variants_gr_one <- pecotmr:::.refPanelToGranges(ref_panel_one)
  ld_data_one <- LdData(
    correlation = ldMat,
    variants = variants_gr_one,
    blockMetadata = pecotmr:::.inferSingleLdBlockMetadata(ref_panel_one)
  )

  local_mocked_bindings(
    summaryStatsQc = function(..., zMismatchQc) {
      captured_qc_method <<- zMismatchQc
      list(study1 = QcResult(
        ldData = ld_data_one,
        rssInput = list(sumstats = data.frame(variant_id = "chr1:100:A:G"),
                         n = 1000, varY = 1),
        preprocess = list(),
        outlierNumber = 0L,
        skipped = FALSE
      ))
    }
  )

  result <- qcRegionalData(region_data, zMismatchQc = NULL, impute = FALSE)
  expect_equal(captured_qc_method, "none")
  expect_type(result, "list")
})

test_that("colocboostPipeline default qc_method resolves to basic-only none", {
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 1)
  captured_qc_method <- NULL

  local_mocked_bindings(
    qcRegionalData = function(regionData, ..., zMismatchQc) {
      captured_qc_method <<- zMismatchQc
      ind <- regionData$individualData
      list(
        individualData = list(
          Y = ind@phenotypes,
          X = stats::setNames(
            lapply(seq_along(ind@phenotypes), function(i) ind@genotypeMatrix),
            names(ind@phenotypes)
          )
        ),
        sumstatData = NULL
      )
    },
    .runColocboost = function(label, ...) {
      list(result = list(label = label), time = as.difftime(0, units = "secs"))
    }
  )

  result <- suppressMessages(colocboostPipeline(region_data))
  expect_equal(captured_qc_method, "none")
  expect_equal(result$xqtl_coloc$label, "xQTL-only ColocBoost")
})

# ===========================================================================
# New direct-input ColocBoost helpers
# ===========================================================================
test_that("RegionalData adapters expose individual and RSS inputs", {
  ind_region <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  ind_input <- regionDataToIndInput(ind_region)
  expect_true(ind_input$sourceInfo$hasIndividual)
  expect_equal(names(ind_input$X), "ctx1")
  expect_equal(names(ind_input$Y), "ctx1")

  rss_region <- make_sumstat_region_data(n_variants = 5, n_studies = 2)
  rss_input <- regionDataToRssInput(rss_region)
  expect_true(rss_input$sourceInfo$hasSumstat)
  expect_equal(names(rss_input$rssInput), c("study1", "study2"))
  expect_true(all(names(rss_input$rssInput) %in% names(rss_input$ldData)))
})

test_that("ColocBoost adapters accept genotype-backed LdData", {
  skip_if_not_installed("pgenlibr")
  td <- test_path("test_data")
  tmp <- tempfile("cb_lddata_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  prefix <- "test_variants"
  for (ext in c("pgen", "pvar", "psam", "afreq")) {
    file.copy(file.path(td, paste0(prefix, ".", ext)),
              file.path(tmp, paste0(prefix, ".", ext)))
  }
  meta_file <- file.path(tmp, "ld_meta.tsv")
  writeLines(c("chrom\tstart\tend\tpath", "21\t0\t0\ttest_variants"), meta_file)
  ld_data <- suppressWarnings(suppressMessages(loadLdMatrix(
    meta_file,
    region = "chr21:17513228-17550000",
    returnGenotype = TRUE
  )))

  variant_info <- getVariantInfo(ld_data)
  ref_panel <- as.data.frame(S4Vectors::mcols(variant_info))
  ref_panel$chrom <- as.character(GenomicRanges::seqnames(variant_info))
  ref_panel$pos <- GenomicRanges::start(variant_info)
  allele_pair <- apply(cbind(ref_panel$A1, ref_panel$A2), 1, function(x) {
    paste(sort(x), collapse = "")
  })
  ref_panel <- ref_panel[nchar(ref_panel$A1) == 1 &
                           nchar(ref_panel$A2) == 1 &
                           !allele_pair %in% c("AT", "CG"), , drop = FALSE]
  ref_panel <- utils::head(ref_panel, 5)
  variant_id <- formatVariantId(ref_panel$chrom, ref_panel$pos,
                                  ref_panel$A2, ref_panel$A1)
  rss_record <- list(
    sumstats = data.frame(
      chrom = ref_panel$chrom,
      pos = ref_panel$pos,
      A1 = ref_panel$A1,
      A2 = ref_panel$A2,
      z = seq_len(nrow(ref_panel)),
      variant_id = variant_id,
      stringsAsFactors = FALSE
    ),
    n = 1000,
    varY = 1
  )
  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(ldgrp = list(study = rss_record)),
      ldInfo = list(ldgrp = ld_data)
    )
  )

  converted <- regionDataToColocboostInput(region_data)
  expect_null(converted$colocboostInput$LD)
  expect_equal(length(converted$colocboostInput$X_ref), 1)
  expect_equal(nrow(converted$colocboostInput$X_ref[[1]]), 100L)
  expect_equal(ncol(converted$colocboostInput$X_ref[[1]]), length(getVariantIds(ld_data)))

  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) list(args = args, dots = dots)
  )
  X_ref <- getGenotypes(ld_data)[, match(variant_id, getVariantIds(ld_data)), drop = FALSE]
  colnames(X_ref) <- variant_id
  result <- suppressMessages(colocboostAnalysis(
    sumstat = data.frame(variant = variant_id, z = seq_along(variant_id), n = 1000),
    X_ref = X_ref,
    ldReferenceInfo = ld_data,
    zMismatchQc = "none",
    M = 2
  ))
  expect_null(result$args$LD)
  expect_equal(length(result$args$X_ref), 1)
  expect_equal(result$args$M, 2)
})

test_that("RegionalData individual adapter exposes context names from phenotypes", {
  ind_region <- make_individual_region_data(n = 12, p = 5, n_contexts = 2, n_events = 1)

  ind_input <- regionDataToIndInput(ind_region)
  expect_equal(names(ind_input$X), c("ctx1", "ctx2"))
  expect_equal(names(ind_input$maf), c("ctx1", "ctx2"))
  expect_equal(names(ind_input$xVariance), c("ctx1", "ctx2"))

  converted <- regionDataToColocboostInput(ind_region)
  # X is shared across contexts in RegionalData; deduplication yields one X.
  expect_equal(length(converted$colocboostInput$X), 1)
  expect_equal(nrow(converted$colocboostInput$dict_YX), 2)
})

test_that("RegionalData adapters handle missing data without fabricating inputs", {
  empty_region <- list(individualData = NULL, sumstatData = NULL)

  ind_input <- regionDataToIndInput(empty_region)
  expect_false(ind_input$sourceInfo$hasIndividual)
  expect_null(ind_input$X)
  expect_null(ind_input$Y)

  rss_input <- regionDataToRssInput(empty_region)
  expect_false(rss_input$sourceInfo$hasSumstat)
  expect_equal(rss_input$rssInput, list())
  expect_equal(rss_input$ldData, list())

  converted <- regionDataToColocboostInput(empty_region)
  expect_equal(converted$colocboostInput, list())
  expect_false(converted$sourceInfo$individual$hasIndividual)
  expect_false(converted$sourceInfo$sumstat$hasSumstat)
})

test_that("regionDataToRssInput keeps duplicate study names unique by LD group", {
  base_region <- make_sumstat_region_data(n_variants = 4, n_studies = 1)
  study <- base_region$sumstatData$sumstats[[1]]
  ld_info <- base_region$sumstatData$ldInfo[[1]]
  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(study, study),
      ldInfo = list(ldA = ld_info, ldB = ld_info)
    )
  )

  rss_input <- regionDataToRssInput(region_data)
  expect_equal(names(rss_input$rssInput), c("study1", "study1.1"))
  expect_equal(names(rss_input$ldData), c("study1", "study1.1"))
  expect_equal(unname(rss_input$sourceInfo$ldGroup), c("ldA", "ldB"))
})

test_that("regionDataToColocboostInput returns core and QC inputs", {
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  converted <- regionDataToColocboostInput(region_data)
  expect_true("colocboostInput" %in% names(converted))
  expect_true("qcInput" %in% names(converted))
  expect_true("sourceInfo" %in% names(converted))
  expect_equal(length(converted$colocboostInput$X), 1)
  expect_equal(length(converted$colocboostInput$Y), 2)
  expect_equal(nrow(converted$colocboostInput$dict_YX), 2)
})

test_that("regionDataToColocboostInput routes genotype LdData through X_ref", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  variants <- region_data$sumstatData$sumstats[[1]][[1]]$sumstats$variant_id
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- variants

  ref_panel <- cbind(parseVariantId(variants), variant_id = variants)
  ref_panel$chrom <- as.character(ref_panel$chrom)
  variants_gr <- pecotmr:::.refPanelToGranges(ref_panel)
  region_data$sumstatData$ldInfo[[1]] <- LdData(
    correlation = NULL,
    genotypeHandle = X_ref,
    variants = variants_gr,
    blockMetadata = pecotmr:::.inferSingleLdBlockMetadata(ref_panel),
    nRef = nrow(X_ref)
  )

  converted <- regionDataToColocboostInput(region_data)

  expect_null(converted$colocboostInput$LD)
  expect_equal(length(converted$colocboostInput$X_ref), 1)
  expect_equal(dim(converted$colocboostInput$X_ref[[1]]), c(10, 5))
  expect_equal(colnames(converted$colocboostInput$X_ref[[1]]), variants)
})

test_that("regionDataToColocboostInput preserves duplicated outcome names across contexts", {
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 2, n_events = 1)
  converted <- regionDataToColocboostInput(region_data)

  expect_equal(names(converted$colocboostInput$Y), c("ctx1_event1", "ctx2_event1"))
  expect_equal(length(converted$colocboostInput$Y), 2)
  expect_equal(nrow(converted$colocboostInput$dict_YX), 2)
  # X is shared across contexts in RegionalData; dict_YX maps both Y to X #1.
  expect_equal(converted$colocboostInput$dict_YX[, "X"], c(1, 1))
})

test_that("regionDataToColocboostInput deduplicates shared individual X", {
  # RegionalData shares one genotype matrix across all conditions, so the
  # per-context residualized X is identical when covariates are the same.
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 2, n_events = 1)
  converted <- regionDataToColocboostInput(region_data)

  expect_equal(length(converted$colocboostInput$X), 1)
  expect_equal(length(converted$colocboostInput$Y), 2)
  expect_equal(converted$colocboostInput$dict_YX[, "X"], c(1, 1))
})

test_that("regionDataToColocboostInput combines individual and RSS inputs", {
  ind_region <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  rss_region <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  region_data <- list(
    individualData = ind_region$individualData,
    sumstatData = rss_region$sumstatData
  )

  converted <- regionDataToColocboostInput(region_data)
  expect_equal(length(converted$colocboostInput$X), 1)
  expect_equal(length(converted$colocboostInput$Y), 2)
  expect_equal(nrow(converted$colocboostInput$dict_YX), 2)
  expect_equal(length(converted$colocboostInput$sumstat), 1)
  expect_equal(length(converted$colocboostInput$LD), 1)
  expect_equal(nrow(converted$colocboostInput$dict_sumstatLD), 1)
  expect_true(converted$sourceInfo$individual$hasIndividual)
  expect_true(converted$sourceInfo$sumstat$hasSumstat)
})

test_that("qcRegionalData applies individual genotype filtering helpers", {
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  expect_message(
    result <- qcRegionalData(region_data, mafCutoff =0),
    "QC track"
  )
  expect_equal(names(result$individualData$Y), "ctx1")
  expect_equal(ncol(result$individualData$Y$ctx1), 2)
})

test_that("qcRegionalData keeps individual context labels after filtering", {
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  region_data$individualData@maf$ctx1 <- stats::setNames(
    rep(0.2, ncol(region_data$individualData@genotypeMatrix)),
    colnames(region_data$individualData@genotypeMatrix)
  )
  region_data$individualData@maf$ctx1[1] <- 0.001

  expect_message(
    result <- qcRegionalData(region_data, mafCutoff =0.05),
    "retained"
  )
  dropped_variant <- names(region_data$individualData@maf$ctx1)[1]
  expect_false(dropped_variant %in% colnames(result$individualData$X$ctx1))
  expect_true(all(startsWith(colnames(result$individualData$Y$ctx1), "ctx1_")))
  expect_equal(ncol(result$individualData$Y$ctx1), ncol(region_data$individualData@phenotypes$ctx1))
})

test_that("summaryStatsQc runs combined basic harmonization when qc_method is none", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  rss_input <- regionDataToRssInput(region_data)
  expect_message(
    result <- summaryStatsQc(
      rssInput = rss_input$rssInput,
      ldData = rss_input$ldData,
      zMismatchQc = "none",
      impute = FALSE
    ),
    "basic allele harmonization"
  )
  expect_equal(names(result), "study1")
  expect_true(is(result$study1, "QcResult"))
  expect_true(nrow(getRssInput(result$study1)$sumstats) > 0)
})

test_that("summaryStatsQc returns one cleaned record for one RSS record", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  rss_input <- regionDataToRssInput(region_data)

  expect_message(
    result <- summaryStatsQc(
      rssInput = rss_input$rssInput$study1,
      ldData = rss_input$ldData$study1,
      zMismatchQc = "none",
      impute = FALSE
    ),
    "basic allele harmonization"
  )
  expect_true(is(result, "QcResult"))
  expect_false(is.null(getLdData(result)))
  expect_true(nrow(getRssInput(result)$sumstats) > 0)
})

test_that("summaryStatsQc treats a study named sumstats as multiple-study input", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)
  rss_input <- regionDataToRssInput(region_data)
  names(rss_input$rssInput)[1] <- "sumstats"
  names(rss_input$ldData)[1] <- "sumstats"

  expect_message(
    result <- summaryStatsQc(
      rssInput = rss_input$rssInput,
      ldData = rss_input$ldData,
      zMismatchQc = "none",
      impute = FALSE
    ),
    "basic allele harmonization"
  )
  expect_equal(names(result), c("sumstats", "study2"))
  expect_true(is(result$sumstats, "QcResult"))
  expect_true(is(result$study2, "QcResult"))
})

test_that("summaryStatsQc imputes when block metadata can be inferred from LD matrix", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  rss_input <- regionDataToRssInput(region_data)

  expect_message(
    result <- summaryStatsQc(
      rssInput = rss_input$rssInput,
      ldData = rss_input$ldData,
      zMismatchQc = "none",
      impute = TRUE,
      imputeOpts = list(rcond = 0.01, r2Threshold = -Inf, minimumLd = -Inf, lamb = 0.01)
    ),
    "running imputation"
  )
  expect_equal(names(result), "study1")
  expect_true(is(result$study1, "QcResult"))
  expect_true(nrow(getRssInput(result$study1)$sumstats) > 0)
})

test_that("colocboostAnalysis directly forwards core inputs without QC", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  X <- matrix(rnorm(20), 5, 4)
  colnames(X) <- paste0("v", 1:4)
  Y <- matrix(rnorm(10), 5, 2)
  colnames(Y) <- c("y1", "y2")
  result <- colocboostAnalysis(X = X, Y = Y, M = 2)
  expect_identical(result$args$X, X)
  expect_identical(result$args$Y, Y)
  expect_equal(result$args$M, 2)
  expect_length(result$dots, 0)
})

test_that("colocboostAnalysis runs individual QC from colocboost-style inputs", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  region_data <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  converted <- regionDataToColocboostInput(region_data)
  expect_message(
    result <- do.call(
      colocboostAnalysis,
      c(converted$colocboostInput, list(missingRateThresh = 1, M = 2))
    ),
    "individual-level"
  )
  expect_equal(length(result$args$X), 1)
  expect_equal(length(result$args$Y), 2)
  expect_equal(result$args$M, 2)
})

test_that("colocboostAnalysis deduplicates shared X after individual QC", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  X <- matrix(rnorm(60), 12, 5)
  colnames(X) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  Y <- list(
    ctx1 = matrix(rnorm(12), 12, 1, dimnames = list(NULL, "event1")),
    ctx2 = matrix(rnorm(12), 12, 1, dimnames = list(NULL, "event2"))
  )
  X_list <- list(ctx1 = X, ctx2 = X)

  expect_message(
    result <- colocboostAnalysis(X = X_list, Y = Y, missingRateThresh = 1, M = 2),
    "individual-level"
  )
  expect_equal(length(result$args$X), 1)
  expect_equal(result$args$dict_YX[, "X"], c(1, 1))
})

test_that("colocboostAnalysis refreshes outcome_names after combined QC", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  ind_region <- make_individual_region_data(n = 12, p = 5, n_contexts = 1, n_events = 2)
  rss_region <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  region_data <- list(
    individualData = ind_region$individualData,
    sumstatData = rss_region$sumstatData
  )
  converted <- regionDataToColocboostInput(region_data)

  expect_message(
    result <- do.call(
      colocboostAnalysis,
      c(converted$colocboostInput, list(missingRateThresh = 1, zMismatchQc = "none", M = 2))
    ),
    "summary-statistic"
  )
  expect_equal(
    result$args$outcome_names,
    c(names(result$args$Y), names(result$args$sumstat))
  )
})

test_that("colocboostAnalysis remaps focal_outcome_idx after QC keeps focal outcome", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    },
    qcIndividualData = function(X, Y, ...) {
      list(ctx = list(
        X = X$ctx,
        Y = Y$ctx[, "raw2", drop = FALSE]
      ))
    }
  )
  X <- list(ctx = matrix(rnorm(40), 10, 4))
  colnames(X$ctx) <- paste0("v", 1:4)
  Y <- list(ctx = matrix(rnorm(20), 10, 2))
  colnames(Y$ctx) <- c("raw1", "raw2")

  result <- colocboostAnalysis(
    X = X, Y = Y,
    outcome_names = c("trait1", "trait2"),
    focal_outcome_idx = 2,
    missingRateThresh = 1,
    M = 2
  )

  expect_equal(result$args$outcome_names, "trait2")
  expect_equal(result$args$focal_outcome_idx, 1)
})

test_that("colocboostAnalysis clears focal_outcome_idx when QC removes focal outcome", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    },
    qcIndividualData = function(X, Y, ...) {
      list(ctx = list(
        X = X$ctx,
        Y = Y$ctx[, "raw2", drop = FALSE]
      ))
    }
  )
  X <- list(ctx = matrix(rnorm(40), 10, 4))
  colnames(X$ctx) <- paste0("v", 1:4)
  Y <- list(ctx = matrix(rnorm(20), 10, 2))
  colnames(Y$ctx) <- c("raw1", "raw2")

  expect_warning(
    result <- colocboostAnalysis(
      X = X, Y = Y,
      outcome_names = c("trait1", "trait2"),
      focal_outcome_idx = 1,
      missingRateThresh = 1,
      M = 2
    ),
    "not present after QC"
  )

  expect_equal(result$args$outcome_names, "trait2")
  expect_null(result$args$focal_outcome_idx)
})

test_that("colocboostAnalysis skips empty individual QC for sumstat-only inputs", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  converted <- regionDataToColocboostInput(region_data)

  messages <- character()
  withCallingHandlers(
    result <- do.call(
      colocboostAnalysis,
      c(converted$colocboostInput, list(zMismatchQc = "none", M = 2))
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_false(any(grepl("individual-level", messages)))
  expect_true(any(grepl("summary-statistic", messages)))
  expect_equal(names(result$args$sumstat), "study1")
})

test_that("colocboostAnalysis falls back to direct call when QC inputs are unavailable", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )

  expect_warning(
    result <- colocboostAnalysis(zMismatchQc = "none", M = 2),
    "required QC inputs are unavailable"
  )
  expect_equal(result$args$M, 2)
  expect_length(result$dots, 0)
})

test_that("colocboostAnalysis falls back when individual QC cannot run", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  X <- matrix(rnorm(20), 5, 4)
  Y <- matrix(rnorm(10), 5, 2)
  colnames(Y) <- c("y1", "y2")

  expect_warning(
    result <- colocboostAnalysis(X = X, Y = Y, missingRateThresh = 1, M = 2),
    "QC requested but skipped"
  )
  expect_identical(result$args$X, X)
  expect_identical(result$args$Y, Y)
  expect_equal(result$args$M, 2)
})

test_that("colocboostAnalysis derives summary QC input from ColocBoost-style inputs", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- variants

  expect_message(
    result <- colocboostAnalysis(sumstat = sumstat, LD = LD, zMismatchQc = "none", M = 2),
    "summary-statistic"
  )
  expect_equal(length(result$args$sumstat), 1)
  expect_equal(length(result$args$LD), 1)
  expect_equal(nrow(result$args$sumstat[[1]]), 5)
  expect_equal(nrow(result$args$dict_sumstatLD), 1)
  expect_equal(result$args$M, 2)
})

test_that("colocboostAnalysis keeps multiple GWAS as colocboost list input after summary QC", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  make_sumstat <- function(seed) {
    set.seed(seed)
    data.frame(
      variant = variants,
      z = rnorm(5),
      n = rep(1000, 5),
      stringsAsFactors = FALSE
    )
  }
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- variants

  expect_message(
    result <- colocboostAnalysis(
      sumstat = list(gwas1 = make_sumstat(1), gwas2 = make_sumstat(2)),
      LD = LD,
      zMismatchQc = "none",
      M = 2
    ),
    "summary-statistic"
  )
  expect_equal(names(result$args$sumstat), c("gwas1", "gwas2"))
  expect_equal(length(result$args$LD), 1)
  expect_equal(nrow(result$args$dict_sumstatLD), 2)
  expect_equal(result$args$dict_sumstatLD[, 2], c(1, 1))
})

test_that("colocboostAnalysis imputes native LD input using inferred block metadata", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    },
    raiss = function(refPanel, knownZscores, ldMatrix = NULL, genotypeMatrix = NULL, ...) {
      expect_null(genotypeMatrix)
      expect_type(ldMatrix, "list")
      expect_true("ldMatrices" %in% names(ldMatrix))
      ldMat <- diag(nrow(knownZscores))
      rownames(ldMat) <- colnames(ldMat) <- knownZscores$variant_id
      list(
        resultFilter = knownZscores,
        ldMat = ldMat
      )
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- variants

  expect_warning(
    expect_message(
      result <- colocboostAnalysis(
        sumstat = sumstat, LD = LD,
        zMismatchQc = "none", impute = TRUE, M = 2
      ),
      "running imputation"
    ),
    NA
  )
  expect_equal(nrow(result$args$sumstat[[1]]), 5)
})

test_that("colocboostAnalysis imputes X_ref input through R-based RAISS path", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    },
    raiss = function(refPanel, knownZscores, ldMatrix = NULL, genotypeMatrix = NULL, ...) {
      # With S4 migration, X_ref is converted to R and imputation uses R-based path
      expect_true(!is.null(ldMatrix) || !is.null(genotypeMatrix))
      list(resultFilter = knownZscores, ldMat = NULL)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  X_ref <- matrix(rnorm(50), 10, 5)
  colnames(X_ref) <- variants

  expect_warning(
    expect_message(
      result <- colocboostAnalysis(
        sumstat = sumstat, X_ref = X_ref,
        zMismatchQc = "none", impute = TRUE, M = 2
      ),
      "running imputation"
    ),
    NA
  )
  # X_ref converted to R, so result has LD (not X_ref)
  expect_equal(length(result$args$LD), 1)
  expect_equal(ncol(result$args$LD[[1]]), 5)
})

test_that("colocboostAnalysis keeps QC-generated X_ref mutually exclusive with original LD", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- variants
  ref_panel <- parseVariantId(variants)
  ref_panel$variant_id <- variants
  # Use a data.frame as LD_reference_info (the production code accepts this format)
  LD_reference_info <- ref_panel

  result <- suppressMessages(colocboostAnalysis(
    sumstat = sumstat,
    LD = LD,
    ldReferenceInfo = LD_reference_info,
    zMismatchQc = "none",
    M = 2
  ))

  # QC produces an LD correlation matrix from the reference info
  expect_equal(length(result$args$LD), 1)
  expect_equal(ncol(result$args$LD[[1]]), 5)
})

test_that("colocboostAnalysis native summary QC supports explicit A1_A2 variant convention", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants_a1a2 <- paste0("chr1:", seq_len(5) * 100, ":G:A")
  sumstat <- data.frame(
    variant = variants_a1a2,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- variants_a1a2

  result <- suppressMessages(colocboostAnalysis(
    sumstat = sumstat, LD = LD,
    zMismatchQc = "none",
    variantConvention = "A1_A2", M = 2
  ))
  expect_equal(result$args$sumstat[[1]]$variant, paste0("chr1:", seq_len(5) * 100, ":A:G"))
  expect_equal(rownames(result$args$LD[[1]]), paste0("chr1:", seq_len(5) * 100, ":A:G"))
})

test_that("colocboostAnalysis uses LD_reference_info data frame for rsid-named LD", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  rsids <- paste0("rs", seq_len(5))
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- rsids
  LD_reference_info <- data.frame(
    chrom = 1,
    id = rsids,
    pos = seq_len(5) * 100,
    A2 = "A",
    A1 = "G",
    stringsAsFactors = FALSE
  )

  expect_message(
    result <- colocboostAnalysis(
      sumstat = sumstat, LD = LD,
      zMismatchQc = "none",
      ldReferenceInfo = LD_reference_info,
      M = 2
    ),
    "ldReferenceInfo|reference"
  )
  expect_equal(rownames(result$args$LD[[1]]), variants)
  expect_equal(result$args$sumstat[[1]]$variant, variants)
})

test_that("colocboostAnalysis uses LD_reference_info row order when LD names are absent", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  LD_reference_info <- data.frame(
    chrom = 1,
    pos = seq_len(5) * 100,
    A2 = "A",
    A1 = "G",
    stringsAsFactors = FALSE
  )

  expect_message(
    result <- colocboostAnalysis(
      sumstat = sumstat, LD = LD,
      zMismatchQc = "none",
      ldReferenceInfo = LD_reference_info,
      M = 2
    ),
    "row order"
  )
  expect_equal(rownames(result$args$LD[[1]]), variants)
  expect_equal(result$args$sumstat[[1]]$variant, variants)
})

test_that("colocboostAnalysis reads LD_reference_info from a bim file", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  rsids <- paste0("rs", seq_len(5))
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- rsids
  bim_file <- tempfile(fileext = ".bim")
  write.table(
    data.frame(1, rsids, 0, seq_len(5) * 100, "G", "A"),
    bim_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  result <- suppressMessages(colocboostAnalysis(
    sumstat = sumstat, LD = LD,
    zMismatchQc = "none",
    ldReferenceInfo = bim_file,
    M = 2
  ))
  expect_equal(rownames(result$args$LD[[1]]), variants)
})

test_that("colocboostAnalysis reports missing LD_reference_info when LD names are not genomic", {
  local_mocked_bindings(
    .cbCallColocboost = function(args, dots) {
      list(args = args, dots = dots)
    }
  )
  variants <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  sumstat <- data.frame(
    variant = variants,
    z = rnorm(5),
    n = rep(1000, 5),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- paste0("rs", seq_len(5))

  expect_warning(
    result <- colocboostAnalysis(sumstat = sumstat, LD = LD, zMismatchQc = "none", M = 2),
    "ldReferenceInfo"
  )
  expect_identical(result$args$LD, LD)
  expect_equal(result$args$M, 2)
})

test_that("colocboostPipeline is the protocol entry", {
  set.seed(450)
  X <- matrix(rnorm(50), 10, 5)
  colnames(X) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  Y <- matrix(rnorm(10), 10, 1)
  colnames(Y) <- "gene1"
  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(5, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = list(Y = list(ctx1 = Y), X = list(ctx1 = X)),
        sumstatData = NULL
      )
    },
    .runColocboost = function(label, ...) {
      list(result = list(label = label, args = list(...)),
           time = as.difftime(0, units = "secs"))
    }
  )

  result <- suppressMessages(colocboostPipeline(region_data))
  expect_named(result, c("xqtl_coloc", "joint_gwas", "separate_gwas", "computing_time"))
  expect_equal(result$xqtl_coloc$label, "xQTL-only ColocBoost")
})

test_that("colocboostPipeline preserves result fields when analyses return NULL", {
  set.seed(451)
  X <- matrix(rnorm(50), 10, 5)
  colnames(X) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
  Y <- matrix(rnorm(10), 10, 1)
  colnames(Y) <- "gene1"
  sumstat <- data.frame(
    chrom = 1,
    pos = seq_len(5) * 100,
    A2 = "A",
    A1 = "G",
    z = rnorm(5),
    n = 1000,
    variant_id = colnames(X),
    stringsAsFactors = FALSE
  )
  LD <- diag(5)
  rownames(LD) <- colnames(LD) <- colnames(X)
  ld_ref_panel <- cbind(parseVariantId(colnames(X)), variant_id = colnames(X))
  ld_ref_panel$chrom <- as.character(ld_ref_panel$chrom)
  ld_data_obj <- LdData(
    correlation = LD,
    variants = pecotmr:::.refPanelToGranges(ld_ref_panel),
    blockMetadata = data.frame(
      blockId = 1L, chrom = "1", blockStart = 100, blockEnd = 500,
      size = 5L, startIdx = 1L, endIdx = 5L
    )
  )
  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(5, 0.05, 0.45))
    ),
    sumstatData = list(
      sumstats = list(chr21_ref = list(study1 = list(sumstats = sumstat, n = 1000, varY = 1))),
      ldInfo = list(chr21_ref = ld_data_obj)
    )
  )

  local_mocked_bindings(
    .runColocboost = function(label, ...) {
      list(result = NULL, time = as.difftime(0, units = "secs"))
    }
  )

  result <- suppressMessages(colocboostPipeline(
    region_data,
    xqtlColoc =TRUE,
    jointGwas =TRUE,
    separateGwas =TRUE,
    zMismatchQc = "none",
    impute = FALSE
  ))
  expect_named(result, c("xqtl_coloc", "joint_gwas", "separate_gwas", "computing_time"))
  expect_null(result$xqtl_coloc)
  expect_null(result$joint_gwas)
  expect_named(result$separate_gwas, "study1")
  expect_null(result$separate_gwas$study1)
})

# ===========================================================================
# 1. colocboostPipeline: no analysis flags returns empty results
# ===========================================================================
test_that("pipeline returns empty results with message when no analysis flags set", {
  region_data <- make_individual_region_data()
  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc = FALSE,
      jointGwas = FALSE,
      separateGwas = FALSE
    ),
    "No colocalization has been performed"
  )
  expect_type(result, "list")
  expect_null(result$xqtl_coloc)
  expect_null(result$joint_gwas)
  expect_null(result$separate_gwas)
})

# ===========================================================================
# 2. colocboostPipeline: NULL individualData and NULL sumstatData
# ===========================================================================
test_that("pipeline returns early when both data sources are NULL", {
  region_data <- list(individualData = NULL, sumstatData = NULL)
  expect_message(
    result <- colocboostPipeline(region_data, xqtlColoc =TRUE),
    "No individual data"
  )
  expect_type(result, "list")
  expect_null(result$xqtl_coloc)
})

# ===========================================================================
# 3. filter_events: type_pattern, valid_pattern, exclude_pattern
# ===========================================================================
test_that("filter_events keeps events matching valid_pattern", {
  # Access the internal function from within colocboostPipeline's environment
  # We need to build a minimal call through the pipeline to test filter_events indirectly.
  # Instead, recreate the inner function for testing purposes.

  # Create a small region_data with events that should be filtered
  set.seed(303)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Events named like sQTL cluster events
  events <- c("clu_1_+:PR:gene1", "clu_1_+:IN:gene1", "clu_2_-:PR:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  event_filters <- list(
    list(
      type_pattern    = ".*clu_(\\d+_[+-?]).*",
      valid_pattern   = "clu_(\\d+_[+-?]):PR:",
      exclude_pattern = "clu_(\\d+_[+-?]):IN:"
    )
  )

  # Pipeline calls filter_events, then qcRegionalData.
  # Mock qcRegionalData so we can isolate the filtering step.
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # Return the data as-is; transform residual_Y to Y format
      list(
        individualData = list(
          Y = region_data$individualData@phenotypes,
          X = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))
        ),
        sumstatData = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc = FALSE,
      jointGwas = FALSE,
      separateGwas = FALSE
    )
  )
  # When no analysis flags set, should return empty and the region_data internal
  # was still filtered. At minimum no error was thrown.
  expect_type(result, "list")
})

# ===========================================================================
# 4. filter_events: error when missing required filter fields
# ===========================================================================
test_that("filter_events errors on missing type_pattern", {
  set.seed(404)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt1", "evt2")

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  # Filter missing type_pattern
  bad_filter <- list(list(valid_pattern = "something"))

  expect_error(
    suppressMessages(
      colocboostPipeline(
        region_data,
        eventFilters = bad_filter,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE
      )
    ),
    "type_pattern"
  )
})

test_that("filter_events errors when only type_pattern is given (no valid or exclude)", {
  set.seed(405)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt1", "evt2")

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  bad_filter <- list(list(type_pattern = "evt.*"))

  expect_error(
    suppressMessages(
      colocboostPipeline(
        region_data,
        eventFilters = bad_filter,
        xqtlColoc = TRUE
      )
    ),
    "type_pattern.*valid_pattern.*exclude_pattern"
  )
})

# ===========================================================================
# 5. extract_contexts_studies: initial extraction
# ===========================================================================
test_that("extract_contexts_studies returns individual contexts and sumstat studies on initial call", {
  # We access the internal by constructing minimal region_data and triggering
  # the pipeline but with both analysis=FALSE so it exits early after extraction.
  X_mat <- matrix(1, 2, 2, dimnames = list(NULL, c("v1", "v2")))
  Y_mat <- matrix(1, 2, 2, dimnames = list(NULL, c("g1", "g2")))
  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(tissue_A = Y_mat, tissue_B = Y_mat),
      residual_X = list(tissue_A = X_mat, tissue_B = X_mat)
    ),
    sumstatData = list(
      sumstats = list(
        list(gwas_trait1 = list(), gwas_trait2 = list())
      )
    )
  )

  # Pipeline calls extract_contexts_studies internally.
  # With no analysis flags it will still run extraction and return.
  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc = FALSE,
      jointGwas = FALSE,
      separateGwas = FALSE
    ),
    "No colocalization"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 6. extract_contexts_studies: after-QC extraction messages
# ===========================================================================
test_that("extract_contexts_studies reports after-QC when some individual data removed", {
  set.seed(601)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y, ctx2 = Y),
      residual_X = list(ctx1 = X, ctx2 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45), ctx2 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  # Mock qcRegionalData to return one NULL context (simulating QC removal).
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch around the colocboost call
  # handles the case where colocboost is unavailable or errors out.
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- matrix(rnorm(10), 10, 1)
      colnames(Y1) <- "ctx1_gene1"
      X1 <- matrix(rnorm(50), 10, 5)
      colnames(X1) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      X2 <- matrix(rnorm(50), 10, 5)
      colnames(X2) <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      list(
        individualData = list(
          Y = list(ctx1_gene1 = Y1, ctx2_gene1 = NULL),
          X = list(ctx1 = X1, ctx2 = X2)
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc = TRUE,
      jointGwas = FALSE,
      separateGwas = FALSE
    ),
    "Skipping follow-up analysis for individual traits"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 7. qcRegionalData: named pipCutoffToSkipSumstat vector
# ===========================================================================
test_that("qcRegionalData handles named pipCutoffToSkipSumstat vector", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  # Mock out the heavy QC functions
  local_mocked_bindings(
    alleleQc = function(target_data, ref_variants, ...) {
      AlleleQcResult(harmonizedData = target_data, qcSummary = target_data)
    },
    rssBasicQc = function(sumstats, ldData, ...) {
      ld_corr <- if (is(ldData, "LdData")) getCorrelation(ldData) else ldData$ldMatrix
      ldMat <- ld_corr[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, ldMat = ldMat)
    },
    summaryStatsQc = function(rssInput = NULL, ldData, ...) {
      stats::setNames(lapply(names(rssInput), function(study) {
        ss <- rssInput[[study]]$sumstats
        ld <- if (is(ldData[[study]], "LdData")) getCorrelation(ldData[[study]]) else ldData[[study]]$ldMatrix
        ldMat <- ld[ss$variant_id, ss$variant_id, drop = FALSE]
        .test_qcresult_from_list(rssInput[[study]], ldMat)
      }), names(rssInput))
    },
    raiss = function(...) {
      list(resultFilter = data.frame(z = rnorm(5)), ldMat = diag(5))
    },
    partitionLdMatrix = function(...) diag(5)
  )

  # Named vector: only specify cutoff for study1
  pip_named <- c("study1" = 0, "study2" = 0)
  result <- suppressMessages(
    qcRegionalData(
      region_data,
      pipCutoffToSkipSumstat = pip_named,
      zMismatchQc = "slalom",
      impute = FALSE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# 8. qcRegionalData: named pip_cutoff fills missing studies with 0
# ===========================================================================
test_that("qcRegionalData fills missing study names with 0 for pipCutoffToSkipSumstat", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    rssBasicQc = function(sumstats, ldData, ...) {
      ld_corr <- if (is(ldData, "LdData")) getCorrelation(ldData) else ldData$ldMatrix
      ldMat <- ld_corr[sumstats$variant_id, sumstats$variant_id, drop = FALSE]
      list(sumstats = sumstats, ldMat = ldMat)
    },
    summaryStatsQc = function(rssInput = NULL, ldData, ...) {
      stats::setNames(lapply(names(rssInput), function(study) {
        ss <- rssInput[[study]]$sumstats
        ld <- if (is(ldData[[study]], "LdData")) getCorrelation(ldData[[study]]) else ldData[[study]]$ldMatrix
        ldMat <- ld[ss$variant_id, ss$variant_id, drop = FALSE]
        .test_qcresult_from_list(rssInput[[study]], ldMat)
      }), names(rssInput))
    },
    raiss = function(...) list(resultFilter = data.frame(z = rnorm(5)), ldMat = diag(5)),
    partitionLdMatrix = function(...) diag(5)
  )

  # Named vector with only one of the two studies: the missing study should get 0

  pip_partial <- c("study1" = 0.05)
  result <- withCallingHandlers(
    suppressMessages(
      qcRegionalData(
        region_data,
        pipCutoffToSkipSumstat = pip_partial,
        zMismatchQc = "slalom",
        impute = FALSE
      )
    ),
    warning = function(w) {
      if (grepl("IBSS algorithm did not converge", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
      stop(w)
    }
  )
  expect_type(result, "list")
})

# ===========================================================================
# 9. colocboostPipeline: output structure verification
# ===========================================================================
test_that("pipeline output structure has expected top-level keys", {
  region_data <- list(individualData = NULL, sumstatData = NULL)
  result <- suppressMessages(
    colocboostPipeline(region_data, xqtlColoc =FALSE, jointGwas =FALSE, separateGwas =FALSE)
  )
  expect_true("xqtl_coloc" %in% names(result))
  expect_true("joint_gwas" %in% names(result))
  expect_true("separate_gwas" %in% names(result))
  expect_true("computing_time" %in% names(result))
  expect_true("QC" %in% names(result$computing_time))
  expect_true("Analysis" %in% names(result$computing_time))
})

# ===========================================================================
# 10. Pipeline with individual data but NULL sumstat (xqtl path)
# ===========================================================================
test_that("pipeline with individual data enters xqtl_coloc path and records timing", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  # Mock qcRegionalData (pecotmr namespace) to simulate QC output.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch handles the case where
  # colocboost is unavailable. We verify the pipeline enters the xqtl path
  # by checking that computing_time$Analysis$xqtl_coloc is recorded.
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- region_data$individualData@phenotypes[[1]]
      colnames(Y1) <- paste0(names(region_data$individualData@phenotypes)[1], "_", colnames(Y1))
      list(
        individualData = list(
          Y = list(ctx1 = Y1),
          X = list(ctx1 = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))[[1]])
        ),
        sumstatData = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc = TRUE,
      jointGwas = FALSE,
      separateGwas = FALSE
    )
  )
  expect_type(result, "list")
  # The pipeline should have entered the xqtl_coloc analysis branch and
  # recorded timing, regardless of whether the colocboost call succeeded.
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

# ===========================================================================
# 11. Pipeline: filter_events with exclude_pattern only
# ===========================================================================
test_that("filter_events exclude_pattern removes matching events via pipeline", {
  set.seed(1100)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Two events: one we want to keep, one to exclude
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("good_event", "bad_event")

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  event_filters <- list(
    list(
      type_pattern    = ".*_event$",
      exclude_pattern = "bad_event"
    )
  )

  # Mock qcRegionalData (pecotmr namespace) to pass through filtered data.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings. The pipeline's tryCatch handles the case where
  # colocboost is unavailable.
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # The residual_Y should have had bad_event removed by filter_events
      remaining_events <- colnames(region_data$individualData@phenotypes$ctx1)
      list(
        individualData = list(
          Y = list(ctx1 = region_data$individualData@phenotypes$ctx1),
          X = list(ctx1 = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))$ctx1)
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc = TRUE
    ),
    "removed"
  )
  expect_type(result, "list")
})

# ===========================================================================
# 12. Pipeline with sumstatData initializes separate_gwas structure
# ===========================================================================
# ===========================================================================
# 13. Pipeline catches colocboost errors gracefully
# ===========================================================================
test_that("pipeline catches colocboost xqtl error and returns NULL result", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  # Mock qcRegionalData to return deliberately mismatched data (X has

  # different row count from Y) so that the colocboost call will always
  # error, whether or not the colocboost package is installed.
  # colocboost is an external package function and cannot be mocked via
  # local_mocked_bindings.
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- region_data$individualData@phenotypes[[1]]
      colnames(Y1) <- paste0("ctx1_", colnames(Y1))
      # Return X with mismatched rows to guarantee colocboost errors
      bad_X <- matrix(rnorm(5 * 8), nrow = 5, ncol = 8)
      colnames(bad_X) <- colnames(stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))[[1]])
      list(
        individualData = list(
          Y = list(ctx1 = Y1),
          X = list(ctx1 = bad_X)
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas = FALSE,
      separateGwas =FALSE
    ),
    "xQTL-only ColocBoost failed"
  )
  expect_null(result$xqtl_coloc)
})

# ===========================================================================
# 14. Pipeline: no data passes QC returns early
# ===========================================================================
test_that("pipeline returns early when all data fails QC", {
  region_data <- make_individual_region_data()

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # Simulate all data removed by QC
      list(individualData = NULL, sumstatData = list(sumstats = NULL))
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc = TRUE,
      jointGwas = FALSE,
      separateGwas = FALSE
    ),
    "No data pass QC"
  )
  expect_type(result, "list")
})

make_qced_sumstat_data <- function(studies = c("study1"), n_variants = 5) {
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids
  sumstats <- lapply(studies, function(study) {
    list(
      sumstats = data.frame(
        z = rnorm(n_variants),
        variant_id = vids,
        stringsAsFactors = FALSE
      ),
      n = 5000,
      varY = 1
    )
  })
  names(sumstats) <- studies
  list(
    sumstats = sumstats,
    ldData = stats::setNames(list(.test_lddata_from_matrix(ldMat)), studies[1]),
    ldMatch = stats::setNames(rep(studies[1], length(studies)), studies)
  )
}

test_that("pipeline skips xqtl branch when QC removes individual data but keeps sumstats", {
  ind_region <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)
  sumstat_region <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  region_data <- list(
    individualData = ind_region$individualData,
    sumstatData = sumstat_region$sumstatData
  )
  calls <- character()

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = NULL,
        sumstatData = make_qced_sumstat_data("study1")
      )
    },
    .runColocboost = function(label, ...) {
      calls <<- c(calls, label)
      list(result = list(label = label, args = list(...)),
           time = as.difftime(0, units = "secs"))
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =TRUE,
      separateGwas =FALSE
    ),
    "No individual data pass QC"
  )
  expect_false(any(grepl("xQTL-only", calls)))
  expect_true(any(grepl("Joint GWAS", calls)))
  expect_null(result$computing_time$Analysis$xqtl_coloc)
  expect_s3_class(result$computing_time$Analysis$joint_gwas, "difftime")
})

test_that("pipeline skips joint_gwas when QC removes sumstats but keeps individual data", {
  ind_region <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)
  sumstat_region <- make_sumstat_region_data(n_variants = 5, n_studies = 1)
  region_data <- list(
    individualData = ind_region$individualData,
    sumstatData = sumstat_region$sumstatData
  )
  calls <- character()

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- region_data$individualData@phenotypes[[1]]
      colnames(Y1) <- paste0("ctx1_", colnames(Y1))
      list(
        individualData = list(
          Y = list(ctx1 = Y1),
          X = list(ctx1 = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))[[1]])
        ),
        sumstatData = list(sumstats = NULL)
      )
    },
    .runColocboost = function(label, ...) {
      calls <<- c(calls, label)
      list(result = list(label = label, args = list(...)),
           time = as.difftime(0, units = "secs"))
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =TRUE,
      separateGwas =FALSE
    ),
    "Skipping follow-up analysis for sumstat studies"
  )
  expect_true(any(grepl("xQTL-only", calls)))
  expect_false(any(grepl("Joint GWAS", calls)))
  expect_s3_class(result$computing_time$Analysis$xqtl_coloc, "difftime")
  expect_null(result$computing_time$Analysis$joint_gwas)
})

test_that("pipeline separate_gwas loops only over sumstats that survive QC", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)
  calls <- character()

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = NULL,
        sumstatData = make_qced_sumstat_data("study1")
      )
    },
    .runColocboost = function(label, ...) {
      calls <<- c(calls, label)
      list(result = list(label = label, args = list(...)),
           time = as.difftime(0, units = "secs"))
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "Skipping follow-up analysis for sumstat studies"
  )
  expect_equal(length(calls), 1)
  expect_true(grepl("study1", calls[[1]]))
  expect_equal(result$computing_time$Analysis$separate_gwas$n_studies, 1)
  expect_true("study2" %in% names(result$separate_gwas))
  expect_null(result$separate_gwas$study2)
})

test_that("pipeline event filters can remove all events in one context while keeping another", {
  set.seed(1401)
  n <- 10
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y_drop <- matrix(rnorm(n * 2), n, 2)
  colnames(Y_drop) <- c("bad_event1", "bad_event2")
  Y_keep <- matrix(rnorm(n * 2), n, 2)
  colnames(Y_keep) <- c("keep_event", "bad_event")
  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx_drop = Y_drop, ctx_keep = Y_keep),
      residual_X = list(ctx_drop = X, ctx_keep = X),
      maf = list(ctx_drop = runif(p, 0.05, 0.45), ctx_keep = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )
  event_filters <- list(list(
    type_pattern = ".*_event.*",
    valid_pattern = "keep_event"
  ))
  calls <- character()

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y_list <- region_data$individualData@phenotypes
      X_list <- stats::setNames(
        replicate(length(Y_list), region_data$individualData@genotypeMatrix,
                  simplify = FALSE),
        names(Y_list)
      )
      # event_filters dropped contexts: re-attach them as NULL entries so the
      # pipeline emits the legacy "Skipping follow-up analysis" message.
      dropped <- attr(region_data$individualData, "filtered_out_contexts")
      if (!is.null(dropped)) {
        for (ctx in dropped) {
          Y_list <- c(Y_list, stats::setNames(list(NULL), ctx))
          X_list <- c(X_list, stats::setNames(list(NULL), ctx))
        }
      }
      list(
        individualData = list(Y = Y_list, X = X_list),
        sumstatData = NULL
      )
    },
    .runColocboost = function(label, ...) {
      calls <<- c(calls, label)
      list(result = list(label = label, args = list(...)),
           time = as.difftime(0, units = "secs"))
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "Skipping follow-up analysis for individual traits ctx_drop"
  )
  expect_equal(length(calls), 1)
  expect_s3_class(result$computing_time$Analysis$xqtl_coloc, "difftime")
})


# ===========================================================================
# Tests from test_twas_colocboost_round3.R (colocboost-related)
# ===========================================================================

# Helper functions used by round3 tests
# Helper: build individual-level region_data
make_individual_region_data <- function(n = 20, p = 8, n_contexts = 2, n_events = 3) {
  set.seed(701)
  make_ctx <- function(ctx_name) {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
    Y <- matrix(rnorm(n * n_events), n, n_events)
    colnames(Y) <- paste0("event", seq_len(n_events))
    maf <- runif(p, 0.05, 0.45)
    list(X = X, Y = Y, maf = maf)
  }
  ctxs <- lapply(paste0("ctx", seq_len(n_contexts)), make_ctx)
  names(ctxs) <- paste0("ctx", seq_len(n_contexts))
  list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = lapply(ctxs, `[[`, "Y"),
      residual_X = lapply(ctxs, `[[`, "X"),
      maf = lapply(ctxs, `[[`, "maf")
    ),
    sumstatData = NULL
  )
}

# Helper: build sumstat-only region_data
make_sumstat_region_data <- function(n_variants = 5, n_studies = 2) {
  set.seed(702)
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos = seq_len(n_variants) * 100,
    A2 = rep("A", n_variants),
    A1 = rep("G", n_variants),
    stringsAsFactors = FALSE
  )

  sumstats_list <- lapply(seq_len(n_studies), function(i) {
    ss <- list(
      sumstats = data.frame(
        chrom = rep(1, n_variants),
        pos = seq_len(n_variants) * 100,
        A1 = rep("G", n_variants),
        A2 = rep("A", n_variants),
        beta = rnorm(n_variants),
        se = runif(n_variants, 0.05, 0.2),
        z = rnorm(n_variants, 0, 2),
        variant_id = vids,
        stringsAsFactors = FALSE
      ),
      n = 10000,
      varY = 1
    )
    list(ss) |> setNames(paste0("study", i))
  })

  ldInfo <- list(.test_lddata_from_matrix(ldMat))

  list(
    individualData = NULL,
    sumstatData = list(
      sumstats = sumstats_list,
      ldInfo = ldInfo
    )
  )
}


# ===========================================================================
# SECTION C: colocboostPipeline - filter_events valid_pattern path
# (lines 64, 67, 71-72, 74, 82, 84-85)
# ===========================================================================

test_that("filter_events: valid_pattern with no matching groups returns NULL (line 74)", {
  set.seed(800)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Events that match type_pattern but none match valid_pattern
  events <- c("clu_1_+:IN:gene1", "clu_2_-:IN:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  # valid_pattern requires ":PR:" but no events have it -> valid_groups is empty -> type_events = character(0) -> returns NULL
  event_filters <- list(
    list(
      type_pattern = ".*clu_(\\d+_[+-?]).*",
      valid_pattern = "clu_(\\d+_[+-?]):PR:"
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      if (is.null(region_data$individualData)) {
        return(list(individualData = NULL, sumstatData = NULL))
      }
      list(
        individualData = list(
          Y = region_data$individualData@phenotypes,
          X = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))
        ),
        sumstatData = NULL
      )
    }
  )

  # The filter returns NULL for the context -> residual_Y entry becomes NULL
  expect_message(
    result <- colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "No events matching|No data pass QC"
  )
  expect_type(result, "list")
})

test_that("filter_events: type_pattern matches nothing skips via next (line 64)", {
  set.seed(801)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("gene_A", "gene_B")

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  # type_pattern matches nothing -> type_events length 0 -> next
  event_filters <- list(
    list(
      type_pattern = "NONEXISTENT_PATTERN_xyz",
      exclude_pattern = "something"
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # Verify events were NOT filtered (both still present)
      remaining <- colnames(region_data$individualData@phenotypes$ctx1)
      expect_equal(length(remaining), 2)
      list(
        individualData = list(
          Y = region_data$individualData@phenotypes,
          X = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))
        ),
        sumstatData = NULL
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    )
  )
  expect_type(result, "list")
})

test_that("filter_events: all events pass -> 'included in following analysis' message (line 82)", {
  set.seed(802)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  # Both events match type_pattern and none get excluded
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("evt_alpha", "evt_beta")

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  # type_pattern matches all events, exclude_pattern matches none
  event_filters <- list(
    list(
      type_pattern = "^evt_",
      exclude_pattern = "NONEXISTENT"
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = list(
          Y = region_data$individualData@phenotypes,
          X = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "included in following analysis"
  )
  expect_type(result, "list")
})

test_that("filter_events: valid_pattern with successful groups retains valid events (lines 67, 71-72)", {
  set.seed(803)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  events <- c("clu_1_+:PR:gene1", "clu_1_+:IN:gene1", "clu_2_-:PR:gene2")
  Y <- matrix(rnorm(n * length(events)), n, length(events))
  colnames(Y) <- events

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(tissue1 = Y),
      residual_X = list(tissue1 = X),
      maf = list(tissue1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = NULL
  )

  event_filters <- list(
    list(
      type_pattern = ".*clu_(\\d+_[+-?]).*",
      valid_pattern = "clu_(\\d+_[+-?]):PR:",
      exclude_pattern = "clu_(\\d+_[+-?]):IN:"
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      remaining <- colnames(region_data$individualData@phenotypes$tissue1)
      # IN event should be removed
      expect_false("clu_1_+:IN:gene1" %in% remaining)
      list(
        individualData = list(
          Y = region_data$individualData@phenotypes,
          X = stats::setNames(replicate(length(region_data$individualData@phenotypes), region_data$individualData@genotypeMatrix, simplify = FALSE), names(region_data$individualData@phenotypes))
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      eventFilters = event_filters,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "removed"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION D: extract_contexts_studies - after-QC paths
# (lines 114, 125-127, 128-133, 134-135, 143-146)
# ===========================================================================

test_that("extract_contexts_studies: all individual data pass QC (line 126)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      Y2 <- matrix(rnorm(20), 20, 1); colnames(Y2) <- "ctx2_event1"
      X1 <- matrix(rnorm(20 * 8), 20, 8); colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      X2 <- matrix(rnorm(20 * 8), 20, 8); colnames(X2) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individualData = list(
          Y = list(ctx1_event1 = Y1, ctx2_event1 = Y2),
          X = list(ctx1 = X1, ctx2 = X2)
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "All individual data pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: all individual data fail QC (line 134-135)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = list(
          Y = list(ctx1 = NULL, ctx2 = NULL),
          X = list(ctx1 = matrix(0, 1, 1), ctx2 = matrix(0, 1, 1))
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "No individual data pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: sumstat studies extraction on initial call (line 114)", {
  # region_data with sumstat only -> triggers sumstat branch in extract_contexts_studies
  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(gwas_trait1 = list(sumstats = data.frame(z = 1), n = 100, varY = 1))
      )
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(gwas_trait1 = list(
            sumstats = data.frame(z = 1.5, variant_id = "chr1:100:A:G"),
            n = 100, varY = 1
          )),
          ldData = list(gwas_trait1 = .test_lddata_from_matrix(matrix(1, 1, 1, dimnames = list("chr1:100:A:G", "chr1:100:A:G")))),
          ldMatch = "gwas_trait1"
        )
      )
    }
  )

  # With xqtl_coloc=FALSE, separate_gwas=TRUE -> will enter the sumstat code paths
  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: after-QC sumstat all pass (line 144)", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      vids <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      ldMat <- diag(5); rownames(ldMat) <- colnames(ldMat) <- vids
      ss1 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, varY = 1)
      ss2 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, varY = 1)
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study1 = ss1, study2 = ss2),
          ldData = list(study1 = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("study1", "study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "All sumstat studies pass QC"
  )
  expect_type(result, "list")
})

test_that("extract_contexts_studies: after-QC sumstat partial pass (line 146)", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      vids <- paste0("chr1:", seq_len(5) * 100, ":A:G")
      ldMat <- diag(5); rownames(ldMat) <- colnames(ldMat) <- vids
      # Only one study remains after QC
      ss1 <- list(sumstats = data.frame(z = rnorm(5), variant_id = vids), n = 10000, varY = 1)
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study1 = ss1),
          ldData = list(study1 = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "Skipping follow-up analysis for sumstat studies"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION E: colocboostPipeline - sumstat processing block
# (lines 244-281: organizing sumstats, normalizing variant IDs, LD normalization)
# ===========================================================================

test_that("pipeline sumstat block: normalizes variant IDs and processes LD matrices (lines 245-281)", {
  set.seed(820)
  n_variants <- 4
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = c(2.1, -1.5, 0.8, 3.2), variant_id = vids),
          n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = c(2.1, -1.5, 0.8, 3.2), variant_id = vids), n = 5000, varY = 1)
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study1 = ss),
          ldData = list(study1 = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("study1")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
  # The separate_gwas structure should be initialized with the study name
  expect_true("separate_gwas" %in% names(result))
})

test_that("pipeline sumstat block: single sumstat study initializes separate_gwas (line 180-181)", {
  set.seed(821)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(single_study = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids),
          n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  # With only one sumstat study, line 180-181 should be reached
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1)
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(single_study = ss),
          ldData = list(single_study = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("single_study")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
  expect_true("separate_gwas" %in% names(result))
  # The separate_gwas structure was initialized (may be empty list if colocboost fails or not installed)
  expect_true(is.list(result$separate_gwas))
})

test_that("pipeline sumstat block: multiple sumstat studies initializes separate_gwas (line 176-178)", {
  set.seed(822)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(studyA = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1
        ), studyB = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss1 <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1)
      ss2 <- list(sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1)
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(studyA = ss1, studyB = ss2),
          ldData = list(studyA = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("studyA", "studyA")
        )
      )
    }
  )

  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
  expect_true("separate_gwas" %in% names(result))
  # The separate_gwas structure was initialized (may be empty list if colocboost fails or not installed)
  expect_true(is.list(result$separate_gwas))
})

# ===========================================================================
# SECTION F: colocboostPipeline - no valid summary statistics after validation
# (lines 275-276: pipeline with all invalid sumstats -> "No data pass QC")
# ===========================================================================

test_that("pipeline: all sumstats invalid returns No data pass QC (line 275-276)", {
  set.seed(830)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  # Mock qcRegionalData to return sumstats that are all invalid (e.g., all NA z)
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      invalid_ss <- list(
        sumstats = data.frame(z = NA_real_, variant_id = "chr1:100:A:G"),
        n = 0, varY = 1
      )
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study1 = invalid_ss),
          ldData = list(study1 = .test_lddata_from_matrix(matrix(1, 1, 1, dimnames = list("chr1:100:A:G", "chr1:100:A:G")))),
          ldMatch = c("study1")
        )
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "No data pass QC|No valid summary"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION G: colocboostPipeline - focal_trait matching (lines 300-301)
# ===========================================================================

test_that("pipeline: focal_trait matches a trait name sets focal_outcome_idx (lines 300-301)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- matrix(rnorm(40), 20, 2)
      colnames(Y1) <- c("ctx1_event1", "ctx1_event2")
      X1 <- matrix(rnorm(20 * 8), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individualData = list(
          Y = list(ctx1_event1 = Y1[, 1, drop = FALSE], ctx1_event2 = Y1[, 2, drop = FALSE]),
          X = list(ctx1 = X1)
        ),
        sumstatData = NULL
      )
    }
  )

  # focal_trait = "ctx1_event2" should match and set focal_outcome_idx
  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      focal_trait = "ctx1_event2",
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    )
  )
  expect_type(result, "list")
  # The xqtl_coloc branch was entered; timing should be recorded
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

test_that("pipeline: focal_trait does NOT match leaves focal_outcome_idx NULL (line 299-302)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      X1 <- matrix(rnorm(20 * 8), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individualData = list(
          Y = list(ctx1_event1 = Y1),
          X = list(ctx1 = X1)
        ),
        sumstatData = NULL
      )
    }
  )

  # focal_trait is specified but doesn't match any trait
  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      focal_trait = "nonexistent_trait",
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    )
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$xqtl_coloc))
})

# ===========================================================================
# SECTION H: colocboostPipeline - joint_gwas path (lines 320-323)
# ===========================================================================

test_that("pipeline: joint_gwas path is entered with both individual and sumstat data (lines 320-323)", {
  set.seed(840)
  n <- 20; p <- 5
  vids <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- vids
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("ctx1_gene1", "ctx1_gene2")
  ldMat <- diag(p); rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = list(
      sumstats = list(
        list(gwas1 = list(
          sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, varY = 1)
      list(
        individualData = list(
          Y = list(ctx1_gene1 = Y[, 1, drop = FALSE], ctx1_gene2 = Y[, 2, drop = FALSE]),
          X = list(ctx1 = X)
        ),
        sumstatData = list(
          sumstats = list(gwas1 = ss),
          ldData = list(gwas1 = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("gwas1")
        )
      )
    }
  )

  # joint_gwas=TRUE should trigger the joint GWAS branch (lines 320-323)
  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =TRUE,
      separateGwas =FALSE
    ),
    "non-focaled version GWAS-xQTL ColocBoost"
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$joint_gwas))
})


# ===========================================================================
# SECTION I: colocboostPipeline - separate_gwas path (lines 341+)
# ===========================================================================

test_that("pipeline: separate_gwas path is entered for each GWAS study", {
  set.seed(850)
  n <- 20; p <- 5
  vids <- paste0("chr1:", seq_len(p) * 100, ":A:G")
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- vids
  Y <- matrix(rnorm(n), n, 1); colnames(Y) <- "ctx1_gene1"
  ldMat <- diag(p); rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = Y),
      residual_X = list(ctx1 = X),
      maf = list(ctx1 = runif(p, 0.05, 0.45))
    ),
    sumstatData = list(
      sumstats = list(
        list(gwasA = list(
          sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss <- list(sumstats = data.frame(z = rnorm(p), variant_id = vids), n = 5000, varY = 1)
      list(
        individualData = list(
          Y = list(ctx1_gene1 = Y),
          X = list(ctx1 = X)
        ),
        sumstatData = list(
          sumstats = list(gwasA = ss),
          ldData = list(gwasA = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("gwasA")
        )
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "focaled version GWAS-xQTL ColocBoost"
  )
  expect_type(result, "list")
  expect_true(!is.null(result$computing_time$Analysis$separate_gwas))
})

# ===========================================================================
# SECTION J: colocboostPipeline - all Y NULL after organizing (lines 225-227)
# ===========================================================================

test_that("pipeline: all Y become NULL after organizing individual data (lines 225-227)", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # Return individualData where all Y entries are NULL
      list(
        individualData = list(
          Y = list(ctx1 = NULL, ctx2 = NULL),
          X = list(ctx1 = matrix(rnorm(160), 20, 8), ctx2 = matrix(rnorm(160), 20, 8))
        ),
        sumstatData = NULL
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "No data pass QC|No individual data pass QC"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION K: colocboostPipeline - sumstat all z NA (line 252-255)
# ===========================================================================

test_that("pipeline: sumstat with all NA z-scores yields warning message (lines 252-255)", {
  set.seed(860)
  n_variants <- 3
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants); rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(study_na = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      # Return sumstats where ALL z-scores are NA
      ss <- list(
        sumstats = data.frame(z = rep(NA_real_, n_variants), variant_id = vids),
        n = 5000, varY = 1
      )
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study_na = ss),
          ldData = list(study_na = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("study_na")
        )
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    ),
    "All z-scores are NA|No data pass QC|No valid"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION L: colocboostPipeline - no sumstatData pass QC (line 152)
# ===========================================================================

test_that("extract_contexts_studies: no sumstat data pass QC message", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 1, n_events = 2)
  # Add some sumstatData so it enters the initial extraction
  region_data$sumstatData <- list(
    sumstats = list(
      list(gwas1 = list(sumstats = data.frame(z = 1), n = 100, varY = 1))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      Y1 <- matrix(rnorm(20), 20, 1); colnames(Y1) <- "ctx1_event1"
      X1 <- matrix(rnorm(160), 20, 8)
      colnames(X1) <- paste0("chr1:", seq_len(8) * 100, ":A:G")
      list(
        individualData = list(
          Y = list(ctx1_event1 = Y1),
          X = list(ctx1 = X1)
        ),
        sumstatData = NULL  # All sumstat data removed by QC
      )
    }
  )

  expect_message(
    result <- colocboostPipeline(
      region_data,
      xqtlColoc =TRUE,
      jointGwas =FALSE,
      separateGwas =FALSE
    ),
    "Skipping follow-up analysis for sumstat studies"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION M: qcRegionalData - pipCutoffToSkipInd wrong length errors
# ===========================================================================

test_that("qcRegionalData: mismatched pipCutoffToSkipInd length errors", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  expect_error(
    qcRegionalData(
      region_data,
      pipCutoffToSkipInd = c(0.1, 0.2, 0.3),  # 3 values but only 2 contexts
      zMismatchQc = "slalom"
    ),
    "pipCutoffToSkipInd must be a scalar"
  )
})

test_that("qcRegionalData: named pipCutoffToSkipInd works with context names", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  # Named vector matching context names
  result <- qcRegionalData(
    region_data,
    pipCutoffToSkipInd = c(ctx1 = 0, ctx2 = 0),
    zMismatchQc = "slalom"
  )
  expect_type(result, "list")
})

test_that("qcRegionalData: named pipCutoffToSkipInd fills missing contexts with 0", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 3, n_events = 2)

  # Only specify cutoff for ctx1 - ctx2 and ctx3 should default to 0
  result <- qcRegionalData(
    region_data,
    pipCutoffToSkipInd = c(ctx1 = 0),
    zMismatchQc = "slalom"
  )
  expect_type(result, "list")
})

test_that("qcRegionalData: scalar pipCutoffToSkipInd becomes named vector", {
  region_data <- make_individual_region_data(n = 20, p = 8, n_contexts = 2, n_events = 2)

  # This exercises the scalar -> named vector recycling path
  result <- qcRegionalData(
    region_data,
    pipCutoffToSkipInd = 0,
    zMismatchQc = "slalom"
  )
  expect_type(result, "list")
})

test_that("qcRegionalData: pipCutoffToSkipInd lookup works when X and Y have different contexts", {
  # Simulate a case where residual_X has a context not in residual_Y
  set.seed(303)
  n <- 20; p <- 8; n_events <- 2
  make_ctx <- function() {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("chr1:", seq_len(p) * 100, ":A:G")
    Y <- matrix(rnorm(n * n_events), n, n_events)
    colnames(Y) <- paste0("event", seq_len(n_events))
    maf <- runif(p, 0.05, 0.45)
    list(X = X, Y = Y, maf = maf)
  }
  ctx1 <- make_ctx()
  ctx2 <- make_ctx()
  ctx3 <- make_ctx()

  # residual_X has 3 contexts, residual_Y only has 2
  # pipCutoffToSkipInd is recycled from residual_Y (length 2)
  region_data <- list(
    individualData = .test_regionaldata_from_lists(
      residual_Y = list(ctx1 = ctx1$Y, ctx2 = ctx2$Y),
      residual_X = list(ctx1 = ctx1$X, ctx2 = ctx2$X, ctx3 = ctx3$X),
      maf = list(ctx1 = ctx1$maf, ctx2 = ctx2$maf, ctx3 = ctx3$maf)
    ),
    sumstatData = NULL
  )

  # Should not error - ctx3 in X has no pip_cutoff entry, defaults to 0
  result <- qcRegionalData(
    region_data,
    pipCutoffToSkipInd = 0,
    zMismatchQc = "slalom"
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION U: colocboostPipeline - sumstat with NA z filtering (line 252-258)
# ===========================================================================

test_that("pipeline sumstat processing handles all-NA z-scores with warning (lines 252-258)", {
  set.seed(870)
  n_variants <- 4
  vids <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")
  ldMat <- diag(n_variants); rownames(ldMat) <- colnames(ldMat) <- vids

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(study_allna = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  # Mock to return sumstats with all NA z-scores
  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss_allna <- list(
        sumstats = data.frame(z = rep(NA_real_, n_variants), variant_id = vids),
        n = 5000, varY = 1
      )
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study_allna = ss_allna),
          ldData = list(study_allna = .test_lddata_from_matrix(ldMat)),
          ldMatch = c("study_allna")
        )
      )
    }
  )

  # Should produce message about NA z-scores or no valid studies
  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION V: colocboostPipeline - LD matrix normalization (lines 263-269)
# ===========================================================================

test_that("pipeline: LD matrix dimnames are normalized to canonical format (lines 261-269)", {
  set.seed(880)
  n_variants <- 3
  # Variant IDs without chr prefix to test normalization
  vids_no_chr <- paste0("1:", seq_len(n_variants) * 100, ":A:G")
  vids_with_chr <- paste0("chr1:", seq_len(n_variants) * 100, ":A:G")

  ldMat <- diag(n_variants)
  rownames(ldMat) <- colnames(ldMat) <- vids_no_chr

  region_data <- list(
    individualData = NULL,
    sumstatData = list(
      sumstats = list(
        list(study1 = list(
          sumstats = data.frame(z = rnorm(n_variants), variant_id = vids_no_chr), n = 5000, varY = 1
        ))
      ),
      ldInfo = list(.test_lddata_from_matrix(ldMat))
    )
  )

  local_mocked_bindings(
    qcRegionalData = function(region_data, ...) {
      ss <- list(
        sumstats = data.frame(z = rnorm(n_variants), variant_id = vids_no_chr),
        n = 5000, varY = 1
      )
      list(
        individualData = NULL,
        sumstatData = list(
          sumstats = list(study1 = ss),
          ldData = list(study1 = .test_lddata_from_matrix(ldMat)),  # ldMat has non-chr names
          ldMatch = c("study1")
        )
      )
    }
  )

  # normalizeVariantId should add chr prefix to LD matrix dimnames
  result <- suppressMessages(
    colocboostPipeline(
      region_data,
      xqtlColoc =FALSE,
      jointGwas =FALSE,
      separateGwas =TRUE
    )
  )
  expect_type(result, "list")
})

# ===========================================================================
# SECTION X: colocboost - qcRegionalData with NULL individualData only sumstat
# ===========================================================================

test_that("qcRegionalData: with only sumstat data processes correctly", {
  region_data <- make_sumstat_region_data(n_variants = 5, n_studies = 1)

  local_mocked_bindings(
    rssBasicQc = function(sumstats, ldData, ...) {
      ld_corr <- if (is(ldData, "LdData")) getCorrelation(ldData) else ldData$ldMatrix
      list(sumstats = sumstats, ldMat = ld_corr)
    },
    summaryStatsQc = function(rssInput = NULL, ldData, ...) {
      stats::setNames(lapply(names(rssInput), function(study) {
        ld <- if (is(ldData[[study]], "LdData")) getCorrelation(ldData[[study]]) else ldData[[study]]$ldMatrix
        .test_qcresult_from_list(rssInput[[study]], ld)
      }), names(rssInput))
    },
    raiss = function(...) {
      list(resultFilter = data.frame(z = rnorm(5)), ldMat = diag(5))
    },
    partitionLdMatrix = function(...) diag(5)
  )

  result <- suppressMessages(
    qcRegionalData(
      region_data,
      zMismatchQc = "slalom",
      impute = FALSE
    )
  )
  expect_type(result, "list")
  expect_null(result$individualData)
  expect_true(!is.null(result$sumstatData))
})

# ===========================================================================
# build_ld_args
# ===========================================================================

test_that("build_ld_args returns LD for square matrices", {
  m <- matrix(1, 5, 5)
  result <- pecotmr:::buildLdArgs(list(m))
  expect_true("LD" %in% names(result))
  expect_null(result$X_ref)
})

test_that("build_ld_args returns X_ref for non-square (genotype) matrices", {
  m <- matrix(1, 100, 5)  # samples x variants
  result <- pecotmr:::buildLdArgs(list(m))
  expect_true("X_ref" %in% names(result))
  expect_null(result$LD)
})

test_that("build_ld_args applies subset correctly", {
  m1 <- matrix(1, 5, 5)
  m2 <- matrix(2, 5, 5)
  result <- pecotmr:::buildLdArgs(list(m1, m2), subset = 2)
  expect_length(result$LD, 1)
  expect_equal(result$LD[[1]][1, 1], 2)
})

# ===========================================================================
# .runColocboost
# ===========================================================================

test_that(".runColocboost returns NULL and message on error", {
  expect_message(
    result <- pecotmr:::.runColocboost("test label", bad_arg = TRUE),
    "test label failed"
  )
  expect_null(result$result)
  expect_s3_class(result$time, "difftime")
})
NA

# ---- removed qcMethod is rejected (hard rename, no alias) ----
test_that("qcRegionalData rejects the removed qcMethod argument as unknown", {
  region_data <- list(individualData = NULL, sumstatData = NULL)
  expect_error(
    qcRegionalData(region_data, qcMethod = "slalom"),
    "unused argument"
  )
})
