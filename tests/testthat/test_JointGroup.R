# Tests for R/JointGroup.R — the uniform contract for the joint-analysis engine:
# the JointGroup hierarchy (conditions-table identity model), the
# JointDispatchCell wiring row, and the pipeline marker classes. Construction is
# validated, so a malformed group / mistyped cell fails loudly at the source.

.jg_cond <- function(study = "S", context = c("c1", "c2"), trait = "G") {
  data.frame(study = study, context = context, trait = trait,
             stringsAsFactors = FALSE)
}

test_that("JointGroup subclasses construct from a conditions table", {
  X <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("v1", "v2")))
  Y <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("c1", "c2")))
  g <- new("IndividualJointGroup", conditions = .jg_cond(), X = X, Y = Y)
  expect_s4_class(g, "JointGroup")
  expect_s4_class(g, "IndividualJointGroup")
  expect_equal(nrow(g@conditions), 2L)

  Z <- matrix(0, 3, 2); R <- diag(3)
  sg <- new("SumStatsJointGroup", conditions = .jg_cond(), Z = Z, R = R,
            N = c(100, 120))
  expect_s4_class(sg, "JointGroup")
  expect_s4_class(sg, "SumStatsJointGroup")
})

test_that("JointGroup validity rejects malformed groups", {
  X <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("v1", "v2")))
  Y <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("c1", "c2")))
  # 1 condition is valid (the univariate cell); 0 conditions is not.
  expect_s4_class(
    new("IndividualJointGroup", conditions = .jg_cond(context = "c1"),
        X = X[, 1, drop = FALSE], Y = Y[, 1, drop = FALSE]),
    "IndividualJointGroup")
  expect_error(new("IndividualJointGroup",
                   conditions = data.frame(study = character(0),
                                           context = character(0),
                                           trait = character(0)),
                   X = X, Y = Y[, 0, drop = FALSE]),
               ">= 1 condition")
  # Missing identity column.
  expect_error(new("IndividualJointGroup",
                   conditions = data.frame(study = "S", context = c("c1", "c2")),
                   X = X, Y = Y), "must have columns")
  # X/Y row mismatch.
  expect_error(new("IndividualJointGroup", conditions = .jg_cond(),
                   X = X, Y = Y[1:5, , drop = FALSE]), "dimension")
  # ncol(Y) must equal nrow(conditions).
  expect_error(new("IndividualJointGroup",
                   conditions = .jg_cond(context = c("c1", "c2", "c3")),
                   X = X, Y = Y), "ncol\\(Y\\)")
  # Non-square LD.
  expect_error(new("SumStatsJointGroup", conditions = .jg_cond(),
                   Z = matrix(0, 3, 2), R = matrix(0, 3, 2), N = 1), "square")
})

test_that("JointDispatchCell + pipeline markers validate at construction", {
  cell <- new("JointDispatchCell", pattern = "context", dataForm = "individual",
              enumerate = function(data, scope, args) list(), minGroup = 2L)
  expect_s4_class(cell, "JointDispatchCell")
  expect_error(new("JointDispatchCell", pattern = "context", dataForm = "bogus",
                   enumerate = function() NULL, minGroup = 2L), "dataForm")
  expect_error(new("JointDispatchCell", pattern = "context",
                   dataForm = "individual", enumerate = function() NULL,
                   minGroup = 0L), "minGroup")

  fm <- new("FmJointPipeline", config = list(coverage = 0.95))
  tw <- new("TwasJointPipeline", config = list(retainFit = TRUE))
  expect_s4_class(fm, "JointPipeline")
  expect_s4_class(tw, "JointPipeline")
  expect_equal(fm@config$coverage, 0.95)
})
