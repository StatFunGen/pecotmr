# Tests for R/jointEngine.R — the unified joint-analysis engine. A joint fit over
# N conditions is SLICED into N per-context rows (real study/context/trait) that
# carry the ";"-joined co-fit members in jointStudies/jointContexts/jointTraits
# as provenance — exactly like running the univariate method per context, except
# the per-context rows share PIP/CS (fm) / the joint fit (twas). The fits are
# mocked so these assert the engine wiring (per-context expansion), not the fit.

.je_mkGroup <- function(tid, n = 10L) {
  X <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("v1", "v2")))
  Y <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("c1", "c2")))
  new("IndividualJointGroup",
      conditions = data.frame(study = "S", context = c("c1", "c2"),
                              trait = tid, stringsAsFactors = FALSE),
      X = X, Y = Y)
}

.je_synthCell <- function(groups) {
  new("JointDispatchCell", pattern = "context", dataForm = "individual",
      enumerate = function(data, scope, args) groups, minGroup = 2L)
}

# Mock postprocess: called once per condition (the fitter passes conditionIdx);
# returns one FineMappingEntry, so the fitter yields one entry per condition.
.je_mockPostprocess <- function(fit, method, dataX, dataY, coverage,
                                secondaryCoverage, signalCutoff, minAbsCorr,
                                csInput = NULL, af = NULL, region = NULL,
                                conditionIdx = NULL) {
  vids <- colnames(dataX)
  FineMappingEntry(
    variantIds = vids,
    susieFit   = list(method = method, cond = conditionIdx),
    topLoci    = data.frame(variant_id = vids,
                            pip = seq(0.9, by = -0.1, length.out = length(vids)),
                            stringsAsFactors = FALSE))
}

test_that(".runJointCell: cross-context FM expands to per-context rows", {
  set.seed(1)
  cell <- .je_synthCell(list(.je_mkGroup("G1"), .je_mkGroup("G2")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, cvFolds = 0))
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(
    fitMvsusie        = function(...) list(),
    .fmPostprocessOne = .je_mockPostprocess,
    .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mvsusie")
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 4L)                                   # 2 genes x 2 ctx
  expect_equal(as.character(res$context), c("c1", "c2", "c1", "c2"))  # REAL
  expect_equal(as.character(res$trait), c("G1", "G1", "G2", "G2"))
  expect_equal(as.character(res$method), rep("mvsusie", 4L))
  expect_equal(as.character(res$jointContexts), rep("c1;c2", 4L))     # provenance
})

test_that(".runJointCell: cross-context FM uses the per-fold mr.mash CV prior", {
  set.seed(2)
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, cvFolds = 2))

  # Prior mr.mash CV payload. Looked up by the FIXED axes (study=S, trait=G1)
  # with context match-any, so a per-context OR a legacy "joint" row both match.
  U  <- list(K = diag(2))
  fp <- list(dataDrivenPriorMatrices = list(U = U, w = c(K = 1)),
             w0 = c(K_grid1 = 1), V = diag(2))
  sp <- data.frame(Sample = paste0("s", 1:10), Fold = rep(1:2, each = 5),
                   stringsAsFactors = FALSE)
  cvPayload <- list(samplePartition = sp,
                    foldFits = list(fold_1 = fp, fold_2 = fp))
  twEntry <- TwasWeightsEntry(variantIds = c("v1", "v2"), weights = c(0.1, 0.2),
                              fits = fp, cvResult = cvPayload)
  tw <- TwasWeights(study = "S", context = "c1", trait = "G1",
                    method = "mrmash", entry = list(twEntry))

  captured <- NULL
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(
    rescaleCovW0      = function(w0) c(K = 1),
    fitMvsusie        = function(...) list(),
    .fmPostprocessOne = .je_mockPostprocess,
    .fmCrossValidate  = function(X, Y, tokens, methodArgs, fold,
                                 samplePartition = NULL, coverage = 0.95,
                                 pos = NULL, verbose = 1, mvPrior = NULL,
                                 mvPriorCv = NULL) {
      captured <<- list(mvPriorCv = mvPriorCv, samplePartition = samplePartition)
      list(samplePartition = samplePartition, prediction = list(),
           performance = list())
    },
    .fmSliceCv  = function(cv, token) cv,
    .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mvsusie", args = list(twasWeights = tw))
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)                            # per-context rows
  # Per-fold priors built (one per fold) + threaded into CV, sharing the folds.
  expect_false(is.null(captured$mvPriorCv))
  expect_equal(names(captured$mvPriorCv), c("1", "2"))
  expect_identical(captured$samplePartition, sp)
})

# ---- twas column (mr.mash) --------------------------------------------------

.je_fakeMrmashFit <- function() {
  list(dataDrivenPriorMatrices = list(U = list(K = diag(2)), w = c(K = 1)),
       w0 = c(K_grid1 = 1), V = diag(2))
}

.je_mockLearnTwas <- function(X, Y, weightMethods, study, context, trait,
                              retainFits, retainFitDetail, standardized,
                              dataType, verbose, ...) {
  W <- matrix(0.1, ncol(X), ncol(Y),
              dimnames = list(colnames(X), colnames(Y)))
  e <- TwasWeightsEntry(variantIds = colnames(X), weights = W,
                        fits = .je_fakeMrmashFit())
  TwasWeights(study = study, context = context, trait = trait,
              method = "mrmash", entry = list(e))
}

.je_mockTwasCv <- function(X, Y, fold, samplePartitions = NULL, weightMethods,
                           retainFits, ..., verbose) {
  sp <- data.frame(Sample = rownames(X),
                   Fold = rep(1:2, length.out = nrow(X)),
                   stringsAsFactors = FALSE)
  metric6 <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
  list(samplePartition = sp,
       prediction = list(mrmash_predicted = matrix(
         0, nrow(X), ncol(Y), dimnames = list(rownames(X), colnames(Y)))),
       performance = list(mrmash_performance = matrix(
         0, ncol(Y), 6, dimnames = list(colnames(Y), metric6))),
       foldFits = list(fold_1 = list(mrmash_weights = .je_fakeMrmashFit()),
                       fold_2 = list(mrmash_weights = .je_fakeMrmashFit())))
}

test_that(".runJointCell: cross-context twas expands to per-context weight vectors", {
  set.seed(3)
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 0))
  local_mocked_bindings(learnTwasWeights = .je_mockLearnTwas, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$context), c("c1", "c2"))
  expect_equal(as.character(res$trait), c("G1", "G1"))
  expect_equal(as.character(res$jointContexts), c("c1;c2", "c1;c2"))
  # Each row carries that context's weight VECTOR (the matrix column).
  w1 <- getWeights(res$entry[[1L]])
  expect_false(is.matrix(w1)); expect_length(w1, 2L)
  expect_false(is.null(getFits(res$entry[[1L]])))    # shared joint fit on each
})

test_that(".runJointCell: cross-context twas attaches per-condition CV slices", {
  set.seed(4)
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 2))
  local_mocked_bindings(learnTwasWeights = .je_mockLearnTwas,
                        twasWeightsCv = .je_mockTwasCv, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash",
                                 args = list(dataDrivenPriorMatricesCv = list(1, 2)))
  expect_equal(nrow(res), 2L)
  cv <- getCvResult(res$entry[[1L]])
  expect_equal(names(cv$foldFits), c("fold_1", "fold_2"))  # shared per-fold fits
  expect_false(is.null(cv$predictions))                    # this context's slice
  expect_false(is.null(cv$samplePartition))
  expect_false(is.matrix(getWeights(res$entry[[1L]])))     # per-context vector
})

test_that(".runJointCell: cross-context twas CV-only rows (fitFullData=FALSE)", {
  set.seed(5)
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 2, fitFullData = FALSE))
  local_mocked_bindings(twasWeightsCv = .je_mockTwasCv, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_equal(nrow(res), 2L)
  e <- res$entry[[1L]]
  expect_length(getVariantIds(e), 0L)                      # placeholder weights
  expect_null(getWeights(e))
  expect_equal(names(getCvResult(e)$foldFits), c("fold_1", "fold_2"))
})

# ---- sumstats column (RSS; no sample folds) ---------------------------------

.je_mkSsGroup <- function(tid, p = 3L, k = 2L) {
  Z <- matrix(rnorm(p * k), p, k,
              dimnames = list(paste0("v", seq_len(p)), paste0("c", seq_len(k))))
  R <- diag(p); dimnames(R) <- list(paste0("v", seq_len(p)), paste0("v", seq_len(p)))
  new("SumStatsJointGroup",
      conditions = data.frame(study = "S", context = paste0("c", seq_len(k)),
                              trait = tid, stringsAsFactors = FALSE),
      Z = Z, R = R, N = c(100, 120))
}

.je_ssCell <- function(groups) {
  new("JointDispatchCell", pattern = "context", dataForm = "sumstats",
      enumerate = function(data, scope, args) groups, minGroup = 2L)
}

test_that(".runJointCell: cross-context FM sumstats (mvsusie_rss) -> per-context", {
  set.seed(6)
  cell <- .je_ssCell(list(.je_mkSsGroup("G1")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  captured <- NULL
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(
    fitMvsusieRss     = function(Z, R, N, prior_variance, coverage, ...) {
      captured <<- list(N = N); list() },
    .fmPostprocessOne = .je_mockPostprocess,
    .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mvsusie")
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$context), c("c1", "c2"))
  expect_equal(as.character(res$jointContexts), c("c1;c2", "c1;c2"))
  expect_equal(captured$N, 110)   # median(c(100, 120)) passed once to mvsusie_rss
})

test_that(".runJointCell: cross-context twas sumstats (mr.mash.rss) -> per-context", {
  set.seed(7)
  cell <- .je_ssCell(list(.je_mkSsGroup("G1")))
  pipe <- new("TwasJointPipeline", config = list())
  local_mocked_bindings(
    mrmashRssWeights = function(stat, LD, retainFit, fitDetail) {
      W <- matrix(0.2, nrow(LD), ncol(stat$z),
                  dimnames = list(rownames(LD), colnames(stat$z)))
      attr(W, "fit") <- .je_fakeMrmashFit()
      W
    },
    .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$context), c("c1", "c2"))
  expect_equal(as.character(res$jointContexts), c("c1;c2", "c1;c2"))
  expect_false(is.matrix(getWeights(res$entry[[1L]])))
  expect_false(is.null(getFits(res$entry[[1L]])))
})

test_that(".lookupJointCell: present cells resolve, absent cells error", {
  for (df in c("individual", "sumstats")) {
    expect_s4_class(pecotmr:::.lookupJointCell("context", df), "JointDispatchCell")
    expect_s4_class(pecotmr:::.lookupJointCell("trait", df), "JointDispatchCell")
  }
  expect_s4_class(pecotmr:::.lookupJointCell("study", "sumstats"),
                  "JointDispatchCell")
  expect_error(pecotmr:::.lookupJointCell("study", "individual"),
               "No joint dispatch cell")
})

# ---- cross-study pattern (study jointed; sumstats-only) ---------------------

test_that(".runJointCell: cross-study (twas sumstats) -> per-study rows + jointStudies", {
  set.seed(10)
  Z <- matrix(rnorm(6), 3, 2, dimnames = list(paste0("v", 1:3), c("S1", "S2")))
  R <- diag(3); dimnames(R) <- list(paste0("v", 1:3), paste0("v", 1:3))
  grp <- new("SumStatsJointGroup",
    conditions = data.frame(study = c("S1", "S2"), context = "brain",
                            trait = "G1", stringsAsFactors = FALSE),
    Z = Z, R = R, N = c(100, 120))
  cell <- new("JointDispatchCell", pattern = "study", dataForm = "sumstats",
              enumerate = function(data, scope, args) list(grp), minGroup = 2L)
  pipe <- new("TwasJointPipeline", config = list())
  local_mocked_bindings(
    mrmashRssWeights = function(stat, LD, retainFit, fitDetail) {
      W <- matrix(0.2, nrow(LD), ncol(stat$z),
                  dimnames = list(rownames(LD), colnames(stat$z)))
      attr(W, "fit") <- .je_fakeMrmashFit(); W
    }, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$study), c("S1", "S2"))   # study is the jointed axis
  expect_equal(as.character(res$context), c("brain", "brain"))
  expect_equal(as.character(res$trait), c("G1", "G1"))
  expect_equal(as.character(res$jointStudies), c("S1;S2", "S1;S2"))
})

# ---- composed pattern (>1 axis varies) --------------------------------------

test_that(".runJointCell: composed (context + trait vary) -> per-tuple rows", {
  set.seed(11)
  n <- 10L
  X <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("v1", "v2")))
  Y <- matrix(rnorm(n * 3), n, 3,
              dimnames = list(paste0("s", seq_len(n)), c("c1:gA", "c1:gB", "c2:gA")))
  conds <- data.frame(study = "S", context = c("c1", "c1", "c2"),
                      trait = c("gA", "gB", "gA"), stringsAsFactors = FALSE)
  grp <- new("IndividualJointGroup", conditions = conds, X = X, Y = Y)
  cell <- new("JointDispatchCell", pattern = "composed", dataForm = "individual",
              enumerate = function(data, scope, args) list(grp), minGroup = 2L)
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, cvFolds = 0))
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(fitMvsusie = function(...) list(),
                        .fmPostprocessOne = .je_mockPostprocess, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mvsusie")
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 3L)                                  # one row per tuple
  expect_equal(as.character(res$study), rep("S", 3L))          # constant -> fixed
  expect_equal(as.character(res$context), c("c1", "c1", "c2"))  # real per-tuple
  expect_equal(as.character(res$trait), c("gA", "gB", "gA"))
  # Both varying axes carry their provenance member list on every row.
  expect_equal(as.character(res$jointContexts), rep("c1;c2", 3L))
  expect_equal(as.character(res$jointTraits), rep("gA;gB", 3L))
})

# ---- cross-trait pattern (trait jointed, context fixed) ---------------------

.je_mkTraitGroup <- function(cx, n = 10L) {
  X <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("v1", "v2")))
  Y <- matrix(rnorm(n * 2), n, 2,    # conditions are traits here
              dimnames = list(paste0("s", seq_len(n)), c("G1", "G2")))
  new("IndividualJointGroup",
      conditions = data.frame(study = "S", context = cx,
                              trait = c("G1", "G2"), stringsAsFactors = FALSE),
      X = X, Y = Y)
}

.je_traitCell <- function(groups) {
  new("JointDispatchCell", pattern = "trait", dataForm = "individual",
      enumerate = function(data, scope, args) groups, minGroup = 2L)
}

test_that(".runJointCell: cross-trait FM -> per-trait rows (context fixed)", {
  set.seed(8)
  cell <- .je_traitCell(list(.je_mkTraitGroup("brain")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, cvFolds = 0))
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(fitMvsusie = function(...) list(),
                        .fmPostprocessOne = .je_mockPostprocess, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mvsusie")
  expect_s4_class(res, "QtlFineMappingResult")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$context), c("brain", "brain"))   # fixed
  expect_equal(as.character(res$trait), c("G1", "G2"))            # real per-trait
  expect_equal(as.character(res$jointTraits), c("G1;G2", "G1;G2"))
})

test_that(".runJointCell: cross-trait twas -> per-trait weight vectors", {
  set.seed(9)
  cell <- .je_traitCell(list(.je_mkTraitGroup("brain")))
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 0))
  local_mocked_bindings(learnTwasWeights = .je_mockLearnTwas, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_s4_class(res, "TwasWeights")
  expect_equal(nrow(res), 2L)
  expect_equal(as.character(res$context), c("brain", "brain"))
  expect_equal(as.character(res$trait), c("G1", "G2"))
  expect_equal(as.character(res$jointTraits), c("G1;G2", "G1;G2"))
  expect_false(is.matrix(getWeights(res$entry[[1L]])))
})

test_that("fitJointGroup(Individual, Fm): fsusie returns one entry per trait", {
  set.seed(12)
  n <- 10L
  X <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("v1", "v2")))
  Y <- matrix(rnorm(n * 2), n, 2,
              dimnames = list(paste0("s", seq_len(n)), c("G1", "G2")))
  grp <- new("IndividualJointGroup",
    conditions = data.frame(study = "S", context = "brain",
                            trait = c("G1", "G2"), stringsAsFactors = FALSE),
    X = X, Y = Y, pos = c(100, 200))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  captured <- NULL
  local_mocked_bindings(
    fitFsusie = function(...) { captured <<- list(...); list() },
    .fmPostprocessOne = .je_mockPostprocess, .package = "pecotmr")

  entries <- pecotmr:::fitJointGroup(grp, pipe, "fsusie", list())
  expect_type(entries, "list")
  expect_length(entries, 2L)                        # one per trait
  expect_s4_class(entries[[1L]], "FineMappingEntry")
  expect_equal(captured$pos, c(100, 200))           # functional domain threaded
})

test_that("fitJointGroup(Individual, Fm): fsusie without pos errors; unknown token errors", {
  X <- matrix(0, 6, 2, dimnames = list(paste0("s", 1:6), c("v1", "v2")))
  Y <- matrix(0, 6, 2, dimnames = list(paste0("s", 1:6), c("G1", "G2")))
  grp <- new("IndividualJointGroup",
    conditions = data.frame(study = "S", context = "brain",
                            trait = c("G1", "G2"), stringsAsFactors = FALSE),
    X = X, Y = Y)  # no pos
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  expect_error(pecotmr:::fitJointGroup(grp, pipe, "fsusie", list()), "pos")
  expect_error(pecotmr:::fitJointGroup(grp, pipe, "bogus", list()),
               "unsupported token")
})

test_that(".runJointCell: composed/sumstats (context+trait vary) -> per-tuple rows", {
  set.seed(13)
  Z <- matrix(rnorm(9), 3, 3,
              dimnames = list(paste0("v", 1:3), c("c1:gA", "c1:gB", "c2:gA")))
  R <- diag(3); dimnames(R) <- list(paste0("v", 1:3), paste0("v", 1:3))
  conds <- data.frame(study = "S", context = c("c1", "c1", "c2"),
                      trait = c("gA", "gB", "gA"), stringsAsFactors = FALSE)
  grp <- new("SumStatsJointGroup", conditions = conds, Z = Z, R = R,
             N = c(100, 100, 120))
  cell <- new("JointDispatchCell", pattern = "composed", dataForm = "sumstats",
              enumerate = function(data, scope, args) list(grp), minGroup = 2L)
  pipe <- new("TwasJointPipeline", config = list())
  local_mocked_bindings(
    mrmashRssWeights = function(stat, LD, retainFit, fitDetail) {
      W <- matrix(0.2, nrow(LD), ncol(stat$z),
                  dimnames = list(rownames(LD), colnames(stat$z)))
      attr(W, "fit") <- .je_fakeMrmashFit(); W
    }, .package = "pecotmr")

  res <- pecotmr:::.runJointCell(cell, pipe, data = NULL, scope = NULL,
                                 tokens = "mrmash")
  expect_equal(nrow(res), 3L)
  expect_equal(as.character(res$context), c("c1", "c1", "c2"))
  expect_equal(as.character(res$trait), c("gA", "gB", "gA"))
  expect_equal(as.character(res$jointContexts), rep("c1;c2", 3L))
  expect_equal(as.character(res$jointTraits), rep("gA;gB", 3L))
  expect_false(is.matrix(getWeights(res$entry[[1L]])))
})

# ---- SR-TWAS ensemble layer (.twasEnsembleLayer) ----------------------------
# The ensemble orchestration formerly inside .twasWeightsPipelineMatrix is now a
# LAYER ON TOP of per-method fitting: per condition, read each method's retained
# out-of-fold CV predictions + weights + R^2, drop methods below the cutoff
# (needs >= 2), and stack via `ensembleWeights` per context.

.je_ensGroup <- function(n = 40L, nCond = 2L) {
  samp <- paste0("s", seq_len(n))
  new("IndividualJointGroup",
      conditions = data.frame(study = "S", context = paste0("c", seq_len(nCond)),
                              trait = "g", stringsAsFactors = FALSE),
      X = matrix(0, n, 3, dimnames = list(samp, paste0("v", 1:3))),
      Y = matrix(rnorm(n * nCond), n, nCond,
                 dimnames = list(samp, paste0("c", seq_len(nCond)))))
}
# One method's per-condition entries; CV predictions correlate with Y by predCor
# so the R^2 the layer reads is controllable.
.je_ensEntries <- function(group, predCor) {
  Y <- group@Y; vars <- colnames(group@X)
  lapply(seq_len(ncol(Y)), function(r) {
    pr <- predCor * Y[, r] + rnorm(nrow(Y), sd = 0.3); names(pr) <- rownames(Y)
    rsq <- stats::cor(Y[, r], pr)^2
    TwasWeightsEntry(variantIds = vars, weights = rnorm(length(vars)),
      cvResult = list(samplePartition = NULL, predictions = pr,
        metrics = c(corr = sqrt(rsq), rsq = rsq, adj_rsq = rsq, pval = 0.01,
                    RMSE = 1, MAE = 1), foldFits = NULL))
  })
}

test_that(".twasEnsembleLayer: >= 2 methods passing -> per-condition ensemble entries", {
  set.seed(1); g <- .je_ensGroup()
  pte <- list(lasso = .je_ensEntries(g, 0.85), enet = .je_ensEntries(g, 0.70))
  ens <- pecotmr:::.twasEnsembleLayer(g, pte, list(
    ensembleR2Threshold = 0.01, ensembleSolver = "quadprog",
    ensembleAlpha = 1, standardized = FALSE))
  expect_length(ens, 2L)
  expect_s4_class(ens[[1L]], "TwasWeightsEntry")
  expect_length(getWeights(ens[[1L]]), 3L)
  coef <- getCvResult(ens[[1L]])$methodCoef
  expect_true(all(coef >= -1e-8)); expect_equal(sum(coef), 1, tolerance = 1e-6)
})

test_that(".twasEnsembleLayer: < 2 methods pass the R^2 cutoff -> NULL (skip)", {
  set.seed(2); g <- .je_ensGroup()
  pte <- list(lasso = .je_ensEntries(g, 0.85), enet = .je_ensEntries(g, 0.70))
  ens <- pecotmr:::.twasEnsembleLayer(g, pte, list(
    ensembleR2Threshold = 0.999, ensembleSolver = "quadprog",
    ensembleAlpha = 1, standardized = FALSE))
  expect_true(all(vapply(ens, is.null, logical(1))))
})

# ---- engine twas fitter: orchestration absorbed from .twasWeightsPipelineMatrix

test_that("fitJointGroup(twas): leakage warning when a full-data mr.mash prior is reused across folds", {
  set.seed(3)
  g <- .je_mkGroup("G1")
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 2L, ensemble = FALSE))
  local_mocked_bindings(learnTwasWeights = .je_mockLearnTwas,
                        twasWeightsCv = .je_mockTwasCv, .package = "pecotmr")
  expect_warning(
    pecotmr:::fitJointGroup(g, pipe, "mrmash", list(methodList = list(
      mrmash_weights = list(dataDrivenPriorMatrices = list(U = diag(2)))))),
    "information leakage")
})

test_that("fitJointGroup(twas): spike-and-slab pi is estimated from an internal mr.ash fit", {
  set.seed(4); n <- 30L
  X <- matrix(rnorm(n * 3), n, 3,
              dimnames = list(paste0("s", 1:n), paste0("v", 1:3)))
  Y <- matrix(rnorm(n), n, 1, dimnames = list(paste0("s", 1:n), "c1"))
  g <- new("IndividualJointGroup",
           conditions = data.frame(study = "S", context = "c1", trait = "g",
                                   stringsAsFactors = FALSE), X = X, Y = Y)
  pipe <- new("TwasJointPipeline",
              config = list(cvFolds = 0L, ensemble = FALSE, estimatePi = TRUE))
  capturedPi <- NULL
  local_mocked_bindings(
    mrashWeights = function(X, y, ...) {
      out <- matrix(0.05, ncol(X), 1L, dimnames = list(colnames(X), NULL))
      attr(out, "fit") <- list(pi = c(0.8, 0.1, 0.1)); out
    },
    bayesCWeights = function(X, y, pi, ...) {
      capturedPi <<- pi
      matrix(0, ncol(X), 1L, dimnames = list(colnames(X), NULL))
    },
    .package = "pecotmr")
  pecotmr:::fitJointGroup(g, pipe, "bayes_c",
                          list(methodList = list(bayes_c_weights = list())))
  expect_false(is.null(capturedPi))
  expect_equal(as.numeric(capturedPi), 1 - 0.8, tolerance = 1e-8)
})

test_that("fitJointGroup(twas): FM-derived method reuses fine-mapping's CV (handoff, no re-CV)", {
  set.seed(5)
  g <- .je_mkGroup("G1")                       # 2 conditions (c1, c2)
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 2L, ensemble = FALSE))
  samp <- rownames(g@X)
  fmCv <- list(
    samplePartition = data.frame(Sample = samp,
                                 Fold = rep(1:2, length.out = length(samp))),
    prediction = list(mvsusie_predicted = matrix(
      rnorm(length(samp) * 2), length(samp), 2,
      dimnames = list(samp, c("c1", "c2")))),
    performance = list(mvsusie_performance = matrix(
      0.5, 2, 6, dimnames = list(c("c1", "c2"),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")))))
  cvCalled <- FALSE
  local_mocked_bindings(
    twasWeightsCv = function(...) { cvCalled <<- TRUE; .je_mockTwasCv(...) },
    .package = "pecotmr")
  local_mocked_bindings(
    mvsusieWeights = function(X, Y, mvsusieFit = NULL, ...)
      matrix(0.1, ncol(X), ncol(Y), dimnames = list(colnames(X), colnames(Y))),
    .package = "pecotmr")
  entries <- pecotmr:::fitJointGroup(g, pipe, "mvsusie", list(
    methodList = list(mvsusie_weights = list()),
    fittedModels = list(mvsusie = list(dummy = TRUE)),
    fineMappingCv = fmCv))
  expect_false(cvCalled)                        # handoff used, no re-CV
  expect_false(is.null(getCvResult(entries[[1L]])$predictions))
})

# =============================================================================
# Enumerators (pattern x dataForm -> list<JointGroup>)
# -----------------------------------------------------------------------------
# The dispatch table stores each enumerator in a JointDispatchCell@enumerate
# slot; the pipeline reaches them via `cell@enumerate(...)`. covr cannot trace a
# function invoked through such a stored reference, so these call the enumerators
# DIRECTLY (pecotmr:::.enum*) with the per-group X/Y/Z builders mocked, asserting
# the enumeration wiring (scope gating + per-group conditions).
# =============================================================================

# A QtlSumStats-shaped data.frame: the sumstat enumerators only touch
# data$study / data$context / data$trait and nrow(data).
.je_ssDf <- function(studies = "S", contexts = c("c1", "c2"), traits = "t1") {
  expand.grid(study = studies, context = contexts, trait = traits,
              stringsAsFactors = FALSE)
}
# Mock .buildJointSumstatZMatrix: a (p x k) Z plus n vector, keyed by colLabels.
.je_mockJointZ <- function(data, tupleRows, colLabels, errorLabel) {
  p <- 3L
  list(Z = matrix(seq_len(p * length(colLabels)), p, length(colLabels),
                  dimnames = list(paste0("v", seq_len(p)), colLabels)),
       nVec = rep(100, length(colLabels)),
       variantIds = paste0("v", seq_len(p)))
}
.je_mockLd <- function(sketch, vids)
  matrix(0, length(vids), length(vids), dimnames = list(vids, vids))

# A real SE (rowRanges carry the trait coordinates fsusie's `pos` needs).
.je_mkSe <- function(traits = c("G1", "G2"), n = 6L) {
  rng <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = seq(100L, by = 100L, length.out = length(traits)),
                     width = 1L))
  names(rng) <- traits
  SummarizedExperiment::SummarizedExperiment(
    assays = list(e = matrix(0, length(traits), n,
                  dimnames = list(traits, paste0("s", seq_len(n))))),
    rowRanges = rng)
}

test_that(".enumCrossContextIndividual: one group per trait in >= 2 contexts", {
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = c("G1", "G2")))
  local_mocked_bindings(
    getStudy = function(data) "S",
    .buildIndividualCrossContextXY = function(data, tid, scopedContexts,
                                              cisWindow, verbose, label,
                                              region = NULL) {
      if (tid == "G2") return(NULL)                       # skip branch (461)
      X <- matrix(0, 4, 2, dimnames = list(paste0("s", 1:4), c("v1", "v2")))
      Y <- matrix(0, 4, 2, dimnames = list(paste0("s", 1:4), c("c1", "c2")))
      list(X = X, Y = Y, perTraitContexts = c("c1", "c2"))
    },
    .package = "pecotmr")
  g <- pecotmr:::.enumCrossContextIndividual(NULL, scope)
  expect_length(g, 1L)                                    # only G1 survives
  expect_s4_class(g[[1L]], "IndividualJointGroup")
  expect_equal(as.character(g[[1L]]@conditions$context), c("c1", "c2"))
  expect_equal(as.character(g[[1L]]@conditions$trait), c("G1", "G1"))
})

test_that(".enumCrossContextIndividual: study not in scope / < 2 contexts -> empty", {
  local_mocked_bindings(getStudy = function(data) "S", .package = "pecotmr")
  expect_length(pecotmr:::.enumCrossContextIndividual(
    NULL, list(studies = "OTHER", contexts = list(), traits = list())), 0L)
  expect_length(pecotmr:::.enumCrossContextIndividual(
    NULL, list(studies = "S", contexts = list(S = "c1"),
               traits = list(S = "G1"))), 0L)
})

test_that(".enumCrossContextSumstats: groups per (study, trait) with >= 2 contexts", {
  df <- .je_ssDf(studies = "S", contexts = c("c1", "c2"), traits = "t1")
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = "t1"))
  local_mocked_bindings(
    getLdSketch = function(x) "SKETCH",
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  g <- pecotmr:::.enumCrossContextSumstats(df, scope)
  expect_length(g, 1L)
  expect_s4_class(g[[1L]], "SumStatsJointGroup")
  expect_equal(as.character(g[[1L]]@conditions$context), c("c1", "c2"))
})

test_that(".enumCrossContextSumstats: < 2 contexts and < 2 tuple rows skip", {
  local_mocked_bindings(
    getLdSketch = function(x) "SKETCH",
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  # study scoped to two contexts but only one row present -> < 2 tupleRows skip
  df <- .je_ssDf(studies = "S", contexts = "c1", traits = "t1")
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = "t1"))
  expect_length(pecotmr:::.enumCrossContextSumstats(df, scope), 0L)
  # < 2 scoped contexts -> skip
  scope2 <- list(studies = "S", contexts = list(S = "c1"),
                 traits = list(S = "t1"))
  expect_length(pecotmr:::.enumCrossContextSumstats(df, scope2), 0L)
})

test_that(".enumCrossTraitIndividual: one group per context with >= 2 traits + pos", {
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = c("G1", "G2")))
  local_mocked_bindings(
    getStudy = function(data) "S",
    .buildIndividualCrossTraitXY = function(data, cx, scopedTraits, cisWindow,
                                            verbose, label, study,
                                            region = NULL) {
      if (cx == "c2") return(NULL)                        # skip branch (511)
      X <- matrix(0, 4, 2, dimnames = list(paste0("s", 1:4), c("v1", "v2")))
      Y <- matrix(0, 4, 2, dimnames = list(paste0("s", 1:4), c("G1", "G2")))
      list(X = X, Y = Y, traitsHere = c("G1", "G2"), se = .je_mkSe())
    },
    .package = "pecotmr")
  g <- pecotmr:::.enumCrossTraitIndividual(NULL, scope)
  expect_length(g, 1L)
  expect_equal(as.character(g[[1L]]@conditions$context), c("c1", "c1"))
  expect_equal(as.character(g[[1L]]@conditions$trait), c("G1", "G2"))
  expect_equal(g[[1L]]@pos, c(100, 200))                  # rowRanges midpoints
})

test_that(".enumCrossTraitIndividual: study not in scope -> empty", {
  local_mocked_bindings(getStudy = function(data) "S", .package = "pecotmr")
  expect_length(pecotmr:::.enumCrossTraitIndividual(
    NULL, list(studies = "X", contexts = list(), traits = list())), 0L)
})

test_that(".enumCrossTraitSumstats: groups per (study, context) with >= 2 traits", {
  df <- .je_ssDf(studies = "S", contexts = "c1", traits = c("t1", "t2"))
  scope <- list(studies = "S", contexts = list(S = "c1"),
                traits = list(S = c("t1", "t2")))
  local_mocked_bindings(
    getLdSketch = function(x) "SKETCH",
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  g <- pecotmr:::.enumCrossTraitSumstats(df, scope)
  expect_length(g, 1L)
  expect_equal(as.character(g[[1L]]@conditions$trait), c("t1", "t2"))
  # < 2 traits present -> skip
  df1 <- .je_ssDf(studies = "S", contexts = "c1", traits = "t1")
  expect_length(pecotmr:::.enumCrossTraitSumstats(df1, scope), 0L)
})

test_that(".enumCrossStudySumstats: group per (context, trait) in >= 2 studies", {
  df <- .je_ssDf(studies = c("S1", "S2"), contexts = "c1", traits = "t1")
  scope <- list(studies = c("S1", "S2"),
                contexts = list(S1 = "c1", S2 = "c1"),
                traits = list(S1 = "t1", S2 = "t1"))
  local_mocked_bindings(
    getLdSketch = function(x) "SKETCH",
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  g <- pecotmr:::.enumCrossStudySumstats(df, scope)
  expect_length(g, 1L)
  expect_equal(as.character(g[[1L]]@conditions$study), c("S1", "S2"))
  # only one study in scope for the tuple -> filtered to < 2 -> skip
  scope1 <- list(studies = c("S1", "S2"),
                 contexts = list(S1 = "c1", S2 = "other"),
                 traits = list(S1 = "t1", S2 = "t1"))
  expect_length(pecotmr:::.enumCrossStudySumstats(df, scope1), 0L)
})

test_that(".enumComposedIndividual: one group joining every (context, trait) tuple", {
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = c("gA", "gB")))
  local_mocked_bindings(
    getStudy = function(data) "S",
    .buildComposedIndividualXY = function(data, scope, study, cisWindow,
                                          verbose, label, region = NULL) {
      Y <- matrix(0, 4, 3, dimnames = list(paste0("s", 1:4),
                                           c("c1:gA", "c1:gB", "c2:gA")))
      X <- matrix(0, 4, 2, dimnames = list(paste0("s", 1:4), c("v1", "v2")))
      list(X = X, Y = Y, tuples = list())
    },
    .package = "pecotmr")
  g <- pecotmr:::.enumComposedIndividual(NULL, scope)
  expect_length(g, 1L)
  expect_equal(as.character(g[[1L]]@conditions$context), c("c1", "c1", "c2"))
  expect_equal(as.character(g[[1L]]@conditions$trait), c("gA", "gB", "gA"))
})

test_that(".enumComposedIndividual: study not in scope / NULL xy -> empty", {
  local_mocked_bindings(getStudy = function(data) "S",
    .buildComposedIndividualXY = function(...) NULL, .package = "pecotmr")
  expect_length(pecotmr:::.enumComposedIndividual(
    NULL, list(studies = "X")), 0L)                       # study not in scope
  expect_length(pecotmr:::.enumComposedIndividual(
    NULL, list(studies = "S", contexts = list(S = "c1"),
               traits = list(S = "gA"))), 0L)             # xy NULL
})

test_that(".enumUnivariateIndividual: one 1-condition group per (context, trait)", {
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = c("G1", "G2")))
  samp <- paste0("s", 1:5)
  local_mocked_bindings(
    getStudy = function(data) "S",
    getPhenotypes = function(data, contexts) .je_mkSe(c("G1", "G2")),
    .fmResidPheno = function(data, contexts, traitId, naAction = "drop")
      matrix(0, 5, 1, dimnames = list(samp, traitId)),
    .fmResidGeno = function(data, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 5, 2, dimnames = list(samp, c("v1", "v2"))),
    .package = "pecotmr")
  g <- pecotmr:::.enumUnivariateIndividual(NULL, scope)
  expect_length(g, 4L)                                    # 2 ctx x 2 traits
  expect_true(all(vapply(g, function(x) nrow(x@conditions), integer(1)) == 1L))
})

test_that(".enumUnivariateIndividual: too few shared samples skips the tuple", {
  scope <- list(studies = "S", contexts = list(S = "c1"),
                traits = list(S = "G1"))
  local_mocked_bindings(
    getStudy = function(data) "S",
    getPhenotypes = function(data, contexts) .je_mkSe("G1"),
    .fmResidPheno = function(data, contexts, traitId, naAction = "drop")
      matrix(0, 1, 1, dimnames = list("s1", traitId)),     # one sample
    .fmResidGeno = function(data, contexts, traitId = NULL, cisWindow = NULL,
                            region = NULL)
      matrix(0, 1, 2, dimnames = list("s1", c("v1", "v2"))),
    .package = "pecotmr")
  expect_length(pecotmr:::.enumUnivariateIndividual(NULL, scope), 0L)
})

test_that(".enumUnivariateIndividual: study not in scope -> empty", {
  local_mocked_bindings(getStudy = function(data) "S", .package = "pecotmr")
  expect_length(pecotmr:::.enumUnivariateIndividual(
    NULL, list(studies = "OTHER")), 0L)
})

test_that(".enumComposedSumstats: one group per fixed-axis row block", {
  df <- .je_ssDf(studies = "S", contexts = c("c1", "c2"), traits = "t1")
  scope <- list(studies = "S", contexts = list(S = c("c1", "c2")),
                traits = list(S = "t1"))
  local_mocked_bindings(
    getLdSketch = function(x) "SKETCH",
    .enumerateComposedSumstatGroups = function(spec, data, scope)
      list(groups = list(c(1L, 2L)),
           studyCol = c("S", "S"), contextCol = c("c1", "c2"),
           traitCol = c("t1", "t1")),
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  g <- pecotmr:::.enumComposedSumstats(df, scope, args = list(axes = c("context", "trait")))
  expect_length(g, 1L)
  expect_equal(as.character(g[[1L]]@conditions$context), c("c1", "c2"))
})

test_that(".enumComposedSumstats: NULL group index and singleton blocks skip", {
  local_mocked_bindings(getLdSketch = function(x) "SKETCH",
    .enumerateComposedSumstatGroups = function(spec, data, scope) NULL,
    .package = "pecotmr")
  expect_length(pecotmr:::.enumComposedSumstats(
    .je_ssDf(), list(studies = "S")), 0L)                 # gi NULL (646)
  local_mocked_bindings(getLdSketch = function(x) "SKETCH",
    .enumerateComposedSumstatGroups = function(spec, data, scope)
      list(groups = list(1L), studyCol = "S", contextCol = "c1",
           traitCol = "t1"),
    .buildJointSumstatZMatrix = .je_mockJointZ, .fmLdFromSketch = .je_mockLd,
    .package = "pecotmr")
  expect_length(pecotmr:::.enumComposedSumstats(
    .je_ssDf(), list(studies = "S"),
    args = list(axes = c("context", "trait"))), 0L)       # < 2 gIdx (649)
})

# =============================================================================
# fitJointGroup branches not covered by the happy-path tests above
# =============================================================================

test_that("fitJointGroup(Individual, Fm): fsusie honest per-fold CV is attached", {
  set.seed(20); n <- 8L
  X <- matrix(rnorm(n * 2), n, 2, dimnames = list(paste0("s", 1:n), c("v1", "v2")))
  Y <- matrix(rnorm(n * 2), n, 2, dimnames = list(paste0("s", 1:n), c("G1", "G2")))
  grp <- new("IndividualJointGroup",
    conditions = data.frame(study = "S", context = "brain",
                            trait = c("G1", "G2"), stringsAsFactors = FALSE),
    X = X, Y = Y, pos = c(100, 200))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, cvFolds = 3))
  cvCalled <- FALSE
  local_mocked_bindings(
    fitFsusie         = function(...) list(),
    fsusieWeights     = function(fsusieFit, variantIds) NULL,
    .fmPostprocessOne = .je_mockPostprocess,
    .fmCrossValidate  = function(X, Y, token, methodArgs, fold, ...) {
      cvCalled <<- TRUE; list(samplePartition = NULL) },
    .fmSliceCv        = function(cv, token) list(prediction = NULL),
    .fmAttachCv       = function(e, cv) e,
    .package = "pecotmr")
  entries <- pecotmr:::fitJointGroup(grp, pipe, "fsusie", list())
  expect_true(cvCalled)                                   # CV path exercised
  expect_length(entries, 2L)
})

test_that("fitJointGroup(Individual, Fm): SER pre-screen skips when < 2 survivors", {
  g <- .je_mkGroup("G1")
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, verbose = 1))
  local_mocked_bindings(
    .fmScreenActive    = function(cut) TRUE,
    .fmSerScreenColumns = function(X, Y, cut) c(TRUE, FALSE),  # 1 survivor
    .package = "pecotmr")
  entries <- suppressMessages(
    pecotmr:::fitJointGroup(g, pipe, "mvsusie", list(pipCutoffToSkip = 0.8)))
  expect_length(entries, 2L)                              # one per ORIGINAL cond
  expect_true(all(vapply(entries, is.null, logical(1))))  # all-NULL (skipped)
})

test_that("fitJointGroup(Individual, Fm): SER pre-screen keeps a subset of conditions", {
  set.seed(21); n <- 10L
  X <- matrix(rnorm(n * 2), n, 2, dimnames = list(paste0("s", 1:n), c("v1", "v2")))
  Y <- matrix(rnorm(n * 3), n, 3,
              dimnames = list(paste0("s", 1:n), c("c1", "c2", "c3")))
  g <- new("IndividualJointGroup",
    conditions = data.frame(study = "S", context = c("c1", "c2", "c3"),
                            trait = "G1", stringsAsFactors = FALSE), X = X, Y = Y)
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95, verbose = 1))
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  local_mocked_bindings(
    .fmScreenActive     = function(cut) TRUE,
    .fmSerScreenColumns = function(X, Y, cut) c(TRUE, FALSE, TRUE),  # drop c2
    fitMvsusie          = function(...) list(),
    .fmPostprocessOne   = .je_mockPostprocess,
    .package = "pecotmr")
  entries <- suppressMessages(
    pecotmr:::fitJointGroup(g, pipe, "mvsusie", list(pipCutoffToSkip = 0.8)))
  expect_length(entries, 3L)
  expect_null(entries[[2L]])                              # screened-out condition
  expect_s4_class(entries[[1L]], "FineMappingEntry")
  expect_s4_class(entries[[3L]], "FineMappingEntry")
})

test_that("fitJointGroup(SumStats, Fm): fsusie and unknown tokens error", {
  grp <- .je_mkSsGroup("G1")
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  expect_error(pecotmr:::fitJointGroup(grp, pipe, "fsusie", list()),
               "no RSS variant")
  expect_error(pecotmr:::fitJointGroup(grp, pipe, "bogus", list()),
               "unsupported")
})

test_that("fitJointGroup(SumStats, Fm): a reweighted-prior residual variance is threaded", {
  grp <- .je_mkSsGroup("G1")
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  captured <- NULL
  local_mocked_bindings(
    .buildMvsusieReweightedPrior = function(fitParts, conditions, ddCut)
      list(priorVariance = "PV", residualVariance = diag(2)),
    fitMvsusieRss = function(Z, R, N, prior_variance, coverage,
                             residual_variance = NULL, ...) {
      captured <<- residual_variance; list() },
    .fmPostprocessOne = .je_mockPostprocess,
    .package = "pecotmr")
  pecotmr:::fitJointGroup(grp, pipe, "mvsusie", list())
  expect_false(is.null(captured))                         # residual_variance set
})

test_that("fitJointGroup(twas): spike-and-slab pi feeds bayes_b probIn", {
  set.seed(22); n <- 30L
  X <- matrix(rnorm(n * 3), n, 3, dimnames = list(paste0("s", 1:n), paste0("v", 1:3)))
  Y <- matrix(rnorm(n), n, 1, dimnames = list(paste0("s", 1:n), "c1"))
  g <- new("IndividualJointGroup",
           conditions = data.frame(study = "S", context = "c1", trait = "g",
                                   stringsAsFactors = FALSE), X = X, Y = Y)
  pipe <- new("TwasJointPipeline",
              config = list(cvFolds = 0L, ensemble = FALSE, estimatePi = TRUE))
  capturedProbIn <- NULL
  local_mocked_bindings(
    mrashWeights = function(X, y, ...) {
      out <- matrix(0.05, ncol(X), 1L, dimnames = list(colnames(X), NULL))
      attr(out, "fit") <- list(pi = c(0.7, 0.2, 0.1)); out
    },
    bayesBWeights = function(X, y, probIn, ...) {
      capturedProbIn <<- probIn
      matrix(0, ncol(X), 1L, dimnames = list(colnames(X), NULL))
    },
    .package = "pecotmr")
  pecotmr:::fitJointGroup(g, pipe, "bayes_b",
                          list(methodList = list(bayes_b_weights = list())))
  expect_equal(as.numeric(capturedProbIn), 1 - 0.7, tolerance = 1e-8)
})

test_that("fitJointGroup(SumStats, twas): a vector weight without rownames falls back to Z rows", {
  grp <- .je_mkSsGroup("G1", p = 3L, k = 2L)
  pipe <- new("TwasJointPipeline", config = list())
  local_mocked_bindings(
    mrmashRssWeights = function(stat, LD, retainFit, fitDetail) {
      # Return a bare numeric vector (one column collapsed) with no names.
      w <- as.numeric(rep(0.2, nrow(LD) * ncol(stat$z)))
      attr(w, "fit") <- .je_fakeMrmashFit(); w
    },
    .package = "pecotmr")
  res <- pecotmr:::.runJointCell(
    .je_ssCell(list(grp)), pipe, data = NULL, scope = NULL, tokens = "mrmash")
  expect_s4_class(res, "TwasWeights")
  expect_equal(getVariantIds(res$entry[[1L]]), rownames(grp@Z))  # fallback vids
})

# =============================================================================
# .runJointCell + ensemble-layer branches
# =============================================================================

test_that(".runJointCell: empty enumeration -> NULL", {
  emptyCell <- new("JointDispatchCell", pattern = "context",
                   dataForm = "individual",
                   enumerate = function(data, scope, args) list(), minGroup = 2L)
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  expect_null(pecotmr:::.runJointCell(emptyCell, pipe, NULL, NULL, "mvsusie"))
})

test_that(".runJointCell: a fitter returning all-NULL entries yields no rows -> NULL", {
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  local_mocked_bindings(create_mixture_prior = function(...) "PRIOR",
                        .package = "mvsusieR")
  # fitJointGroup returns a list of NULLs (every condition screened out).
  local_mocked_bindings(
    fitJointGroup = function(group, pipeline, token, args)
      vector("list", nrow(group@conditions)),
    .package = "pecotmr")
  expect_null(pecotmr:::.runJointCell(cell, pipe, NULL, NULL, "mvsusie"))
})

test_that(".runJointCell: twas ensemble layer adds 'ensemble' rows on top of >= 2 methods", {
  set.seed(30)
  cell <- .je_synthCell(list(.je_ensGroup(n = 40L, nCond = 2L)))
  pipe <- new("TwasJointPipeline", config = list(
    cvFolds = 2L, ensemble = TRUE, ensembleR2Threshold = 0.01,
    ensembleSolver = "quadprog", ensembleAlpha = 1, standardized = FALSE))
  # Two methods, each returning per-condition entries with CV predictions.
  local_mocked_bindings(
    fitJointGroup = function(group, pipeline, token, args)
      .je_ensEntries(group, if (token == "lasso") 0.85 else 0.7),
    .package = "pecotmr")
  res <- pecotmr:::.runJointCell(cell, pipe, NULL, NULL,
                                 tokens = c("lasso", "enet"))
  expect_s4_class(res, "TwasWeights")
  expect_true("ensemble" %in% as.character(res$method))
})

test_that(".twasEnsembleLayer: entries lacking CV predictions are skipped", {
  set.seed(31); g <- .je_ensGroup(nCond = 1L)
  good <- .je_ensEntries(g, 0.85)
  # A method whose entry carries no CV predictions -> contributes nothing.
  noCv <- list(TwasWeightsEntry(variantIds = colnames(g@X),
                                weights = rnorm(ncol(g@X)), cvResult = NULL))
  ens <- pecotmr:::.twasEnsembleLayer(
    g, list(a = good, b = noCv, c = NULL),
    list(ensembleR2Threshold = 0.01, ensembleSolver = "quadprog",
         ensembleAlpha = 1, standardized = FALSE))
  expect_true(all(vapply(ens, is.null, logical(1))))      # < 2 usable -> NULL
})

test_that(".twasEnsembleLayer: ensembleWeights returning NULL -> NULL entry", {
  set.seed(32); g <- .je_ensGroup(nCond = 1L)
  pte <- list(lasso = .je_ensEntries(g, 0.85), enet = .je_ensEntries(g, 0.70))
  local_mocked_bindings(ensembleWeights = function(...) NULL, .package = "pecotmr")
  ens <- pecotmr:::.twasEnsembleLayer(g, pte, list(
    ensembleR2Threshold = 0.01, ensembleSolver = "quadprog",
    ensembleAlpha = 1, standardized = FALSE))
  expect_true(all(vapply(ens, is.null, logical(1))))
})

test_that(".twasEnsembleLayer: unnamed ensemble weights fall back to a method's variant ids", {
  set.seed(33); g <- .je_ensGroup(nCond = 1L)
  pte <- list(lasso = .je_ensEntries(g, 0.85), enet = .je_ensEntries(g, 0.70))
  local_mocked_bindings(
    ensembleWeights = function(cvResults, Y, twasWeightList, contextIndex,
                               solver, alpha)
      list(ensembleTwasWeights = as.numeric(rep(0.1, ncol(g@X))),  # no names
           methodCoef = c(0.5, 0.5), methodPerformance = c(0.8, 0.7)),
    .package = "pecotmr")
  ens <- pecotmr:::.twasEnsembleLayer(g, pte, list(
    ensembleR2Threshold = 0.01, ensembleSolver = "quadprog",
    ensembleAlpha = 1, standardized = FALSE))
  expect_s4_class(ens[[1L]], "TwasWeightsEntry")
  expect_equal(getVariantIds(ens[[1L]]), colnames(g@X))   # fallback ids
})

# =============================================================================
# Small slicers, .jointTwasCvResult, .twasFmHandoffCv, construct, .runJointSpecs
# =============================================================================

test_that(".fmSliceCvCondition / .sliceTwasCvResultToCondition: NULL passes through", {
  expect_null(pecotmr:::.fmSliceCvCondition(NULL, 1L))
  expect_null(pecotmr:::.sliceTwasCvResultToCondition(NULL, 1L))
})

test_that(".jointTwasCvResult: NULL cv and empty/absent payloads degrade gracefully", {
  expect_null(pecotmr:::.jointTwasCvResult(NULL, "mrmash"))
  # Empty prediction/performance lists and all-NULL foldFits -> NULL components.
  cv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L),
             prediction = list(), performance = list(),
             foldFits = list(fold_1 = list(other = 1)))
  out <- pecotmr:::.jointTwasCvResult(cv, "mrmash")
  expect_null(out$predictions)                            # pickByBase empty (252)
  expect_null(out$foldFits)                               # all-NULL ffKey (260)
  # A method token absent from a non-empty prediction list -> NULL (255).
  cv2 <- list(samplePartition = NULL,
              prediction = list(lasso_predicted = matrix(0, 1, 1)),
              performance = list())
  expect_null(pecotmr:::.jointTwasCvResult(cv2, "mrmash")$predictions)
})

test_that(".twasFmHandoffCv: a token absent from the FM CV predictions -> NULL", {
  fmCv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L),
               prediction = list(susie_predicted = matrix(0, 1, 1)),
               performance = list())
  expect_null(pecotmr:::.twasFmHandoffCv(fmCv, "mvsusie"))
  expect_null(pecotmr:::.twasFmHandoffCv(NULL, "mvsusie"))
})

test_that("construct: empty rows -> NULL for both pipelines", {
  empty <- pecotmr:::.jointRows()
  expect_null(pecotmr:::construct(new("FmJointPipeline", config = list()), empty))
  expect_null(pecotmr:::construct(new("TwasJointPipeline", config = list()), empty))
})

test_that(".runJointSpecs: no methods or no specs -> NULL", {
  pipe <- new("FmJointPipeline", config = list(ldSketch = NULL))
  expect_null(pecotmr:::.runJointSpecs(list(), NULL, "individual", pipe,
                                       jointMethods = "mvsusie",
                                       contexts = NULL, traitIds = NULL))
  expect_null(pecotmr:::.runJointSpecs(list(list(axes = "context")), NULL,
                                       "individual", pipe,
                                       jointMethods = character(0),
                                       contexts = NULL, traitIds = NULL))
})

# =============================================================================
# Remaining branches: .twasGroupArgs CV-partition handoff, a token whose fit
# is NULL, and .runJointSpecs' region-mode trait restriction.
# =============================================================================

test_that(".twasGroupArgs: takes the CV partition from the fine-mapping CV when none is set", {
  g <- .je_mkGroup("G1")                                  # IndividualJointGroup
  pipe <- new("TwasJointPipeline", config = list(cvFolds = 2L))
  sp <- data.frame(Sample = rownames(g@X),
                   Fold = rep(1:2, length.out = nrow(g@X)),
                   stringsAsFactors = FALSE)
  local_mocked_bindings(
    .twasFineMappingFits = function(fineMappingResult, study, context, trait)
      list(),
    .twasCvResultFor = function(fmRes, s, c, t) list(samplePartition = sp),
    .package = "pecotmr")
  out <- pecotmr:::.twasGroupArgs(g, pipe, list(fineMappingResult = "FMR"))
  expect_identical(out$samplePartition, sp)               # line 692
})

test_that(".runJointCell: a token whose fitter returns NULL is skipped", {
  cell <- .je_synthCell(list(.je_mkGroup("G1")))
  pipe <- new("FmJointPipeline", config = list(coverage = 0.95))
  local_mocked_bindings(fitJointGroup = function(...) NULL, .package = "pecotmr")
  expect_null(pecotmr:::.runJointCell(cell, pipe, NULL, NULL, "mvsusie"))  # 755
})

test_that(".runJointSpecs: region mode without traitId restricts scoped traits to the locus", {
  pipe <- new("FmJointPipeline", config = list(ldSketch = NULL))
  region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  captured <- NULL
  local_mocked_bindings(
    # S2 has no scoped contexts -> the region-restriction loop skips it (873).
    .fmResolveSpecScope = function(spec, data, contexts, traitIds)
      list(studies = c("S", "S2"),
           contexts = list(S = c("c1", "c2"), S2 = character(0)),
           traits = list(S = c("g1", "g2"), S2 = "g9")),
    getPhenotypes = function(data, contexts) .je_mkSe(c("g1", "g2")),
    .fmTraitsInRegion = function(se, traits, region) { captured <<- traits; "g1" },
    .lookupJointCell = function(pattern, dataForm) .je_synthCell(list()),
    .package = "pecotmr")
  res <- pecotmr:::.runJointSpecs(
    list(list(axes = "context", scope = NULL)), data = NULL,
    dataForm = "individual", pipeline = pipe, jointMethods = "mvsusie",
    contexts = NULL, traitIds = NULL, args = list(region = region))
  expect_null(res)                                        # empty cell -> NULL
  expect_equal(captured, c("g1", "g2"))                   # 871-875 ran
})
