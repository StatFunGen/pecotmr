// penalized_rss.cpp — Coordinate descent for penalized regression on RSS objective
//
// The package's one RSS coordinate-descent solver, for the penalties
//   LASSO, MCP, SCAD, L0, L0L1, L0L2
// lassosum (Mak et al 2017) is its LASSO at alpha = 1: the same updates in the
// same order, so the same coefficients bit for bit.
//
// Objective (smooth part): min_beta  beta'R beta - 2 beta'z + penalty(beta)
// where R is the LD matrix shrunk to R_s = (1 - s) R + s I (s = `shrink`;
// 0 for an LD matrix already shrunk) and z = bhat / sqrt(n).
//
// LASSO, MCP and SCAD mix in a ridge term the way glmnet and ncvreg do:
// l1 = lambda * alpha, l2 = lambda * (1 - alpha), penalty = 2 P(beta; l1) +
// l2 beta'beta. LASSO with alpha < 1 is the elastic net, alpha = 0 ridge.
//
// The coordinate descent update is identical across penalties; only the
// proximal/thresholding operator differs.  MCP/SCAD operators are ported
// from ncvreg (Breheny & Huang 2011).  L0 thresholding follows L0Learn
// (Hazimeh & Mazumder 2020) with an optional swap phase.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <cmath>
#include <string>

#include "rss_ld.h"

using namespace cpp11;
using namespace arma;

// ---- Penalty enum --------------------------------------------------------

enum Penalty { PEN_LASSO = 0, PEN_MCP = 1, PEN_SCAD = 2, PEN_L0 = 3 };

static Penalty parse_penalty(const std::string& s) {
  if (s == "lasso") return PEN_LASSO;
  if (s == "MCP")   return PEN_MCP;
  if (s == "SCAD")  return PEN_SCAD;
  if (s == "L0" || s == "L0L1" || s == "L0L2") return PEN_L0;
  cpp11::stop("Unknown penalty: '%s'. Use lasso, MCP, SCAD, L0, L0L1, or L0L2.", s.c_str());
  return PEN_LASSO; // unreachable
}

// ---- Thresholding operators ----------------------------------------------
// v = R_jj (diagonal of LD matrix for coordinate j)
// l1, l2 = L1 and L2 penalty weights (for elastic net: l1 = lambda*alpha,
//          l2 = lambda*(1-alpha))
//
// Each operator is the exact minimizer of the one-coordinate problem
// v b^2 - 2 z b + 2 P(b; l1) + l2 b^2, whose curvature is a = v + l2. ncvreg
// writes the curvature v (1 + l2), which is the same only when v = 1: a
// correlation matrix's diagonal, but not a diagonal the shrinkage or a
// genotype factor leaves slightly off 1.

static inline double lasso_thresh(double z, double l1, double l2, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  return std::copysign(az - l1, z) / (v + l2);
}

// MCP (minimax concave penalty) — Breheny & Huang 2011
// gamma > 1 controls concavity; gamma -> inf recovers LASSO
static inline double mcp_thresh(double z, double l1, double l2,
                                double gamma, double v) {
  double az = std::abs(z);
  double a = v + l2;
  if (az <= l1) return 0.0;
  if (az <= gamma * l1 * a)
    return std::copysign(az - l1, z) / (a - 1.0 / gamma);
  return z / a;
}

// SCAD (smoothly clipped absolute deviation) — Fan & Li 2001
// gamma > 2 controls transition; default 3.7
static inline double scad_thresh(double z, double l1, double l2,
                                 double gamma, double v) {
  double az = std::abs(z);
  double a = v + l2;
  if (az <= l1) return 0.0;
  if (az <= l1 * (1.0 + a))
    return std::copysign(az - l1, z) / a;
  if (az <= gamma * l1 * a)
    return std::copysign(az - gamma * l1 / (gamma - 1.0), z) /
           (a - 1.0 / (gamma - 1.0));
  return z / a;
}

// L0 (+L1+L2) — Hazimeh & Mazumder 2020
// Soft-threshold for L1, then hard-threshold for L0.
// lambda0 = L0 penalty weight, l1 = L1 weight, l2 = L2 weight
static inline double l0_thresh(double z, double lambda0,
                               double l1, double l2, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  double denom = v + 2.0 * l2;
  double reg = (az - l1) / denom;
  double thr = std::sqrt(2.0 * lambda0 / denom);
  if (reg < thr) return 0.0;
  return std::copysign(reg, z);
}

// ---- Per-coordinate thresholding dispatch --------------------------------

static inline double apply_threshold(Penalty pen, double t, double v,
                                     double l1, double l2,
                                     double gamma, double lambda0) {
  switch (pen) {
    case PEN_LASSO: return lasso_thresh(t, l1, l2, v);
    case PEN_MCP:   return mcp_thresh(t, l1, l2, gamma, v);
    case PEN_SCAD:  return scad_thresh(t, l1, l2, gamma, v);
    case PEN_L0:    return l0_thresh(t, lambda0, l1, l2, v);
  }
  return 0.0; // unreachable
}

// ---- Full penalty value (for fbeta computation) --------------------------

static double mcp_penalty_single(double bj, double lambda, double gamma) {
  double ab = std::abs(bj);
  if (ab <= lambda * gamma)
    return lambda * ab - ab * ab / (2.0 * gamma);
  return 0.5 * lambda * lambda * gamma;
}

static double scad_penalty_single(double bj, double lambda, double gamma) {
  double ab = std::abs(bj);
  if (ab <= lambda)
    return lambda * ab;
  if (ab <= gamma * lambda)
    return (2.0 * gamma * lambda * ab - ab * ab - lambda * lambda) /
           (2.0 * (gamma - 1.0));
  return lambda * lambda * (gamma + 1.0) / 2.0;
}

// The penalty the thresholding operators minimize: 2 P(beta; l1) + l2
// beta'beta, so MCP and SCAD are evaluated at l1 = lambda * alpha, the
// threshold the operators use, not at lambda.
static double compute_penalty(Penalty pen, const vec& beta,
                               double lambda, double alpha, double gamma,
                               double lambda0, double lambda2) {
  double val = 0.0;
  double l1 = lambda * alpha;
  double l2 = lambda * (1.0 - alpha);
  int p = beta.n_elem;
  switch (pen) {
    case PEN_LASSO:
      val = 2.0 * l1 * sum(abs(beta)) + l2 * dot(beta, beta);
      break;
    case PEN_MCP:
      for (int j = 0; j < p; j++)
        val += 2.0 * mcp_penalty_single(beta(j), l1, gamma);
      val += l2 * dot(beta, beta);
      break;
    case PEN_SCAD:
      for (int j = 0; j < p; j++)
        val += 2.0 * scad_penalty_single(beta(j), l1, gamma);
      val += l2 * dot(beta, beta);
      break;
    case PEN_L0: {
      int nnz = 0;
      for (int j = 0; j < p; j++)
        if (beta(j) != 0.0) nnz++;
      val = 2.0 * lambda0 * nnz +
            2.0 * l1 * sum(abs(beta)) +
            2.0 * lambda2 * dot(beta, beta);
      break;
    }
  }
  return val;
}

// ---- Single-block coordinate descent -------------------------------------
// `ld` answers R's entries for the current beta (rss_ld.h): DenseLd tracks
// R * beta as a running sum, FactorLd tracks F * beta.

template <class Ld>
static int penalized_cd_rss(Penalty pen, double lambda, double gamma,
                            double alpha, double lambda0, double lambda2,
                            Ld& ld, const vec& z,
                            double thr, vec& beta, int maxiter) {
  int p = z.n_elem;
  double l1, l2;

  // For LASSO/MCP/SCAD the elastic net decomposition applies
  if (pen != PEN_L0) {
    l1 = lambda * alpha;
    l2 = lambda * (1.0 - alpha);
  } else {
    // For L0 variants, lambda controls the L1 part directly,
    // lambda2 controls L2 separately
    l1 = lambda;
    l2 = lambda2;
  }

  int conv = 0;
  for (int k = 0; k < maxiter; k++) {
    double dlx = 0.0;
    for (int j = 0; j < p; j++) {
      double bj = beta(j);
      double v = ld.diag(j);
      // Partial correlation: z_j - sum_{k!=j} R[j,k]*beta[k]
      double t = z(j) - ld.row(j, beta) + v * bj;

      // Apply penalty-specific thresholding
      beta(j) = apply_threshold(pen, t, v, l1, l2, gamma, lambda0);

      if (beta(j) == bj) continue;
      double del = beta(j) - bj;
      dlx = std::max(dlx, std::abs(del));
      // Update the running R * beta
      ld.update(j, del);
    }
    cpp11::check_user_interrupt();
    if (dlx < thr) {
      conv = 1;
      break;
    }
  }
  return conv;
}

// ---- L0 swap optimization (single block) ---------------------------------
// Attempts to improve the L0 solution by swapping non-zero coefficients
// with zero ones that have higher correlation with the residual.
// Ported from L0Learn CDL012Swaps.cpp.

template <class Ld>
static bool l0_swap_round(double lambda0, double l1, double l2,
                          Ld& ld, const vec& z,
                          double thr, vec& beta, int maxiter) {
  int p = beta.n_elem;
  bool found_better = false;

  // Collect non-zero indices
  std::vector<int> nnz_idx;
  for (int i = 0; i < p; i++)
    if (beta(i) != 0.0) nnz_idx.push_back(i);

  for (int ii = 0; ii < (int)nnz_idx.size(); ii++) {
    int i = nnz_idx[ii];
    // Compute partial correlation for all variables with beta[i] added back
    // riX_j = z_j - Rbeta_j + R_ji * beta_i
    vec riX = z - ld.product(beta) + ld.column(i) * beta(i);

    double max_corr = -1.0;
    int max_idx = -1;
    for (int j = 0; j < p; j++) {
      if (beta(j) == 0.0 && std::abs(riX(j)) > max_corr) {
        max_corr = std::abs(riX(j));
        max_idx = j;
      }
    }
    if (max_idx < 0) continue;

    // Check if swap is worthwhile
    double denom = ld.diag(max_idx) + 2.0 * l2;
    if (max_corr > denom * std::abs(beta(i)) + l1) {
      // Perform swap: zero out i, set j to optimal
      double old_bi = beta(i);
      beta(i) = 0.0;
      ld.update(i, -old_bi);

      double t_j = z(max_idx) - ld.row(max_idx, beta) + ld.diag(max_idx) * 0.0;
      beta(max_idx) = l0_thresh(t_j, lambda0, l1, l2, ld.diag(max_idx));
      if (beta(max_idx) != 0.0) {
        ld.update(max_idx, beta(max_idx));
      }

      // Re-run CD from swapped solution.
      // For PEN_L0: inside penalized_cd_rss, l1 = lambda param, l2 = lambda2.
      penalized_cd_rss(PEN_L0, l1, 0.0, 1.0, lambda0, l2,
                       ld, z, thr, beta, maxiter);
      found_better = true;
      break; // restart from scratch
    }
  }
  return found_better;
}

// ---- One block at one lambda: coordinate descent, then the L0 swaps ------

template <class Ld>
static int penalized_block(Penalty pen, double lambda, double gamma,
                           double alpha, double lambda0, double lambda2,
                           Ld& ld, const vec& z, double thr, vec& beta,
                           int maxiter, int maxSwaps) {
  int conv = penalized_cd_rss(pen, lambda, gamma, alpha, lambda0, lambda2,
                              ld, z, thr, beta, maxiter);
  // L0 swap phase (per block)
  if (pen == PEN_L0 && maxSwaps > 0) {
    for (int sw = 0; sw < maxSwaps; sw++) {
      bool improved = l0_swap_round(lambda0, lambda, lambda2, ld, z, thr,
                                    beta, maxiter);
      if (!improved) break;
    }
  }
  return conv;
}

// The lambda path's results, in the shape both entry points return.
static cpp11::writable::list penalized_path_result(const mat& beta_mat,
                                                   const vec& lambda,
                                                   const Col<int>& conv_vec,
                                                   const vec& loss_vec,
                                                   const vec& fbeta_vec) {
  using namespace cpp11::literals;
  writable::list result({
    "beta"_nm   = as_doubles_matrix(beta_mat),
    "lambda"_nm = as_doubles(lambda),
    "conv"_nm   = as_integers(conv_vec),
    "loss"_nm   = as_doubles(loss_vec),
    "fbeta"_nm  = as_doubles(fbeta_vec)
  });
  return result;
}

// ---- Registered entry points ---------------------------------------------

// `shrink` is s: each block is fitted against R_s = (1 - s) R + s I, which
// is never formed (DenseLd).
[[cpp11::register]]
cpp11::writable::list penalizedRssCpp(
    const doubles& zR,
    const list& LD,
    double shrink,
    const doubles& lambdaR,
    const std::string& penaltyStr,
    double gamma,
    double alpha,
    double lambda0,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps) {

  Penalty pen = parse_penalty(penaltyStr);
  vec z = as_Col(zR);
  vec lambda = as_Col(lambdaR);

  // Cache LD blocks
  int n_blocks = LD.size();
  std::vector<mat> ld_blocks(n_blocks);
  std::vector<int> block_start(n_blocks), block_end(n_blocks);
  int p = 0;
  for (int b = 0; b < n_blocks; b++) {
    ld_blocks[b] = as_Mat(doubles_matrix<>(LD[b]));
    block_start[b] = p;
    p += ld_blocks[b].n_rows;
    block_end[b] = p - 1;
  }

  if ((int)z.n_elem != p)
    cpp11::stop("Length of z must equal total rows across all LD blocks.");

  int nlambda = lambda.n_elem;
  mat beta_mat(p, nlambda, fill::zeros);
  Col<int> conv_vec(nlambda, fill::zeros);
  vec loss_vec(nlambda, fill::zeros);
  vec fbeta_vec(nlambda, fill::zeros);

  // Working beta — warm-started across lambda path
  vec beta(p, fill::zeros);

  for (int i = 0; i < nlambda; i++) {
    // Block-wise coordinate descent
    int out = 1;
    for (int b = 0; b < n_blocks; b++) {
      const mat& Rb = ld_blocks[b];
      int s = block_start[b];
      int e = block_end[b];
      vec z_blk = z.subvec(s, e);
      vec beta_blk = beta.subvec(s, e);
      DenseLd ld(Rb, shrink, beta_blk);

      int conv_blk = penalized_block(pen, lambda(i), gamma, alpha, lambda0,
                                     lambda2, ld, z_blk, thr, beta_blk,
                                     maxiter, maxSwaps);
      beta.subvec(s, e) = beta_blk;
      out = std::min(out, conv_blk);
    }

    beta_mat.col(i) = beta;
    conv_vec(i) = out;

    // Compute loss = beta'R_s beta - 2 z'beta (block-wise)
    double loss = -2.0 * dot(z, beta);
    for (int b = 0; b < n_blocks; b++) {
      const mat& Rb = ld_blocks[b];
      int s = block_start[b];
      int e = block_end[b];
      vec beta_blk = beta.subvec(s, e);
      loss += (1.0 - shrink) * as_scalar(beta_blk.t() * Rb * beta_blk)
              + shrink * dot(beta_blk, beta_blk);
    }
    loss_vec(i) = loss;
    fbeta_vec(i) = loss + compute_penalty(pen, beta, lambda(i), alpha,
                                           gamma, lambda0, lambda2);
  }

  return penalized_path_result(beta_mat, lambda, conv_vec, loss_vec, fbeta_vec);
}

// penalizedRssCpp() for one LD block given by a factor F of R = F'F / (m - 1)
// and shrunk to R_s = (1 - s) R + s I: the same lambda path, the same updates
// and swaps, without forming R or R_s.
[[cpp11::register]]
cpp11::writable::list penalizedRssFactorCpp(
    const doubles& zR,
    const doubles_matrix<>& factorR,
    double shrink,
    const doubles& lambdaR,
    const std::string& penaltyStr,
    double gamma,
    double alpha,
    double lambda0,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps) {

  Penalty pen = parse_penalty(penaltyStr);
  vec z = as_Col(zR);
  vec lambda = as_Col(lambdaR);
  mat F = as_Mat(factorR);
  if (F.n_cols != z.n_elem)
    cpp11::stop("The LD factor has %d columns but there are %d z-scores.",
                (int)F.n_cols, (int)z.n_elem);

  int p = z.n_elem;
  int nlambda = lambda.n_elem;
  mat beta_mat(p, nlambda, fill::zeros);
  Col<int> conv_vec(nlambda, fill::zeros);
  vec loss_vec(nlambda, fill::zeros);
  vec fbeta_vec(nlambda, fill::zeros);

  // Working beta — warm-started across lambda path
  vec beta(p, fill::zeros);

  for (int i = 0; i < nlambda; i++) {
    FactorLd ld(F, shrink, beta);
    conv_vec(i) = penalized_block(pen, lambda(i), gamma, alpha, lambda0,
                                  lambda2, ld, z, thr, beta, maxiter, maxSwaps);
    beta_mat.col(i) = beta;
    double loss = -2.0 * dot(z, beta) + ld.quad(beta);
    loss_vec(i) = loss;
    fbeta_vec(i) = loss + compute_penalty(pen, beta, lambda(i), alpha,
                                           gamma, lambda0, lambda2);
  }

  return penalized_path_result(beta_mat, lambda, conv_vec, loss_vec, fbeta_vec);
}
