context("twas: twasZ and harmonize deprecated wrappers")

# ===========================================================================
# Fixture builder: a small genotype + LD matrix.
# ===========================================================================

.tz_makeLd <- function(n = 100, p = 8, seed = 7) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, runif(p, 0.2, 0.8)), nrow = n, ncol = p)
  vid <- paste0("v", seq_len(p))
  colnames(X) <- vid
  af <- colMeans(X) / 2
  Xstd <- sweep(X, 2, 2 * af)
  Xstd <- sweep(Xstd, 2, sqrt(2 * af * (1 - af)), "/")
  R <- crossprod(Xstd) / (n - 1)
  rownames(R) <- colnames(R) <- vid
  list(X = X, Xstd = Xstd, R = R, vid = vid, n = n, p = p)
}

# ===========================================================================
# twasZ: input coercion and validity
# ===========================================================================

test_that("twasZ: numeric vector input is coerced to a single-method matrix", {
  d <- .tz_makeLd()
  w <- rnorm(d$p)
  names(w) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(w, z, R = d$R)
  expect_equal(dim(res$Z), c(1L, 2L))
  expect_equal(colnames(res$Z), c("Z", "pval"))
  expect_equal(rownames(res$Z), "method1")
})

test_that("twasZ: matrix without colnames gets method1/method2 names", {
  d <- .tz_makeLd()
  W <- cbind(rnorm(d$p), rnorm(d$p))
  rownames(W) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(W, z, R = d$R)
  expect_equal(rownames(res$Z), c("method1", "method2"))
})

test_that("twasZ: non-matrix non-numeric weights errors out", {
  expect_error(pecotmr:::twasZ(list(a = 1), z = 1, R = matrix(1)),
               "must be a numeric vector or a matrix")
})

test_that("twasZ: length mismatch between weights and z errors", {
  expect_error(
    pecotmr:::twasZ(c(0.1, 0.2, 0.3), z = c(1, 2)),
    "nrow\\(weights\\) must equal length\\(z\\)"
  )
})

test_that("twasZ: missing R/X/V triplet errors", {
  expect_error(
    pecotmr:::twasZ(c(0.1, 0.2), z = c(1, 2)),
    "provide R, X, or the .V, D, nSketch. SVD triplet"
  )
})

# ===========================================================================
# twasZ: R path
# ===========================================================================

test_that("twasZ: R path with named alignment matches the manual formula", {
  d <- .tz_makeLd()
  set.seed(1)
  z <- rnorm(d$p); names(z) <- d$vid
  w <- rnorm(d$p); names(w) <- d$vid
  res <- pecotmr:::twasZ(w, z, R = d$R)
  expected_stat <- sum(w * z)
  expected_denom <- sqrt(as.numeric(crossprod(w, d$R %*% w)))
  expect_equal(as.numeric(res$Z[, "Z"]),
               expected_stat / expected_denom,
               tolerance = 1e-10)
})

test_that("twasZ: R path realigns by rowname order", {
  d <- .tz_makeLd()
  z <- rnorm(d$p)
  w <- rnorm(d$p); names(w) <- d$vid
  # Permute R rows/columns; results must be identical to the un-permuted R.
  perm <- sample(seq_len(d$p))
  R_perm <- d$R[perm, perm]
  res_orig <- pecotmr:::twasZ(w, z, R = d$R)
  res_perm <- pecotmr:::twasZ(w, z, R = R_perm)
  expect_equal(res_perm$Z, res_orig$Z, tolerance = 1e-10)
})

test_that("twasZ: R path errors when R is missing rows named in weights", {
  d <- .tz_makeLd()
  w <- rnorm(d$p); names(w) <- d$vid
  R_short <- d$R[1:(d$p - 1), 1:(d$p - 1)]
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
    "R is missing rows for"
  )
})

test_that("twasZ: R path positional alignment errors on dim mismatch", {
  d <- .tz_makeLd()
  w <- rnorm(d$p)   # unnamed -> positional alignment
  R_short <- unname(d$R[1:(d$p - 1), 1:(d$p - 1)])
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
    "positional alignment requires nrow\\(R\\) == nrow\\(weights\\)"
  )
})

# ===========================================================================
# twasZ: X path (computes R via computeLd)
# ===========================================================================

test_that("twasZ: X path computes R via computeLd and matches that R explicitly", {
  d <- .tz_makeLd()
  set.seed(2)
  w <- rnorm(d$p); names(w) <- d$vid
  z <- rnorm(d$p)
  # The X path delegates to computeLd(X, method = "sample"), which uses
  # cor(X). Compare against the same correlation matrix passed explicitly.
  R_from_X <- cor(d$X)
  rownames(R_from_X) <- colnames(R_from_X) <- d$vid
  res_X   <- pecotmr:::twasZ(w, z, X = d$X)
  res_Rfx <- pecotmr:::twasZ(w, z, R = R_from_X)
  expect_equal(res_X$Z, res_Rfx$Z, tolerance = 1e-10)
})

# ===========================================================================
# twasZ: SVD path coverage (V-row missing error)
# ===========================================================================

test_that("twasZ: SVD path errors when V is missing rows named in weights", {
  d <- .tz_makeLd()
  s <- svd(d$Xstd)
  rownames(s$v) <- d$vid[seq_len(nrow(s$v))]
  w <- rnorm(d$p); names(w) <- c(d$vid[1:(d$p - 1)], "ghost")
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p),
                    V = s$v, D = s$d, nSketch = d$n),
    "V is missing rows for"
  )
})

test_that("twasZ: SVD path positional alignment errors on dim mismatch", {
  d <- .tz_makeLd()
  s <- svd(d$Xstd)
  # Pass an unnamed V with the wrong nrow.
  w <- rnorm(d$p)
  V_short <- s$v[1:(d$p - 1), , drop = FALSE]
  expect_error(
    pecotmr:::twasZ(w, rnorm(d$p),
                    V = V_short, D = s$d, nSketch = d$n),
    "positional alignment requires nrow\\(V\\) == nrow\\(weights\\)"
  )
})

# ===========================================================================
# twasZ: combineMethods integration
# ===========================================================================

test_that("twasZ: combineMethods K=1 returns the per-tuple p-value unchanged", {
  d <- .tz_makeLd()
  w <- rnorm(d$p); names(w) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(w, z, R = d$R, combineMethods = c("acat", "bonferroni"))
  expect_false(is.null(res$combined))
  expect_equal(res$combined$results$acat$pval, res$Z[1, "pval"])
  expect_equal(res$combined$results$bonferroni$pval, res$Z[1, "pval"])
  expect_equal(res$combined$input$nValid, 1L)
})

test_that("twasZ: combineMethods K>=2 forwards to combinePValues with correlation", {
  d <- .tz_makeLd()
  set.seed(11)
  W <- cbind(method_a = rnorm(d$p), method_b = rnorm(d$p))
  rownames(W) <- d$vid
  z <- rnorm(d$p)
  res <- pecotmr:::twasZ(W, z, R = d$R, combineMethods = "acat")
  expect_false(is.null(res$combined))
  # ACAT combines two p-values; result must be in (0, 1).
  expect_true(res$combined$results$acat$pval >= 0 &&
              res$combined$results$acat$pval <= 1)
  expect_equal(res$combined$input$nValid, 2L)
})

# ===========================================================================
# Deprecated wrappers
# ===========================================================================

test_that("harmonizeTwas is a deprecated no-op", {
  expect_warning(res <- harmonizeTwas(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})

test_that("harmonizeGwas is a deprecated no-op", {
  expect_warning(res <- harmonizeGwas(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})

test_that("twasPipeline is a deprecated no-op", {
  expect_warning(res <- twasPipeline(), "has been removed", ignore.case = TRUE)
  expect_null(res)
})
