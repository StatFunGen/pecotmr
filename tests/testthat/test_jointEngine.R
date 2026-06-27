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
