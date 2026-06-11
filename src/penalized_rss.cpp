// penalized_rss.cpp — Coordinate descent for penalized regression on RSS objective
//
// Generalizes lassosum_rss.cpp to support multiple penalties:
//   LASSO, MCP, SCAD, L0, L0L1, L0L2
//
// Objective (smooth part): min_beta  beta'R beta - 2 beta'z + penalty(beta)
// where R is a (possibly pre-shrunk) LD matrix and z = bhat / sqrt(n).
//
// The coordinate descent update is identical across penalties; only the
// proximal/thresholding operator differs.  MCP/SCAD operators are ported
// from ncvreg (Breheny & Huang 2011).  L0 thresholding follows L0Learn
// (Hazimeh & Mazumder 2020) with an optional swap phase.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <cmath>
#include <string>

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

static inline double lasso_thresh(double z, double l1, double l2, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  return std::copysign(az - l1, z) / (v * (1.0 + l2));
}

// MCP (minimax concave penalty) — Breheny & Huang 2011
// gamma > 1 controls concavity; gamma -> inf recovers LASSO
static inline double mcp_thresh(double z, double l1, double l2,
                                double gamma, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  if (az <= gamma * l1 * (1.0 + l2))
    return std::copysign(az - l1, z) / (v * (1.0 + l2 - 1.0 / gamma));
  return z / (v * (1.0 + l2));
}

// SCAD (smoothly clipped absolute deviation) — Fan & Li 2001
// gamma > 2 controls transition; default 3.7
static inline double scad_thresh(double z, double l1, double l2,
                                 double gamma, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  if (az <= l1 * (2.0 + l2))
    return std::copysign(az - l1, z) / (v * (1.0 + l2));
  if (az <= gamma * l1 * (1.0 + l2))
    return std::copysign(az - gamma * l1 / (gamma - 1.0), z) /
           (v * (1.0 - 1.0 / (gamma - 1.0) + l2));
  return z / (v * (1.0 + l2));
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

static double compute_penalty(Penalty pen, const vec& beta,
                               double lambda, double alpha, double gamma,
                               double lambda0, double lambda2) {
  double val = 0.0;
  double l1 = lambda * alpha;
  int p = beta.n_elem;
  switch (pen) {
    case PEN_LASSO:
      val = 2.0 * l1 * sum(abs(beta));
      break;
    case PEN_MCP:
      for (int j = 0; j < p; j++)
        val += 2.0 * alpha * mcp_penalty_single(beta(j), lambda, gamma);
      // Add ridge part
      if (alpha < 1.0)
        val += lambda * (1.0 - alpha) * dot(beta, beta);
      break;
    case PEN_SCAD:
      for (int j = 0; j < p; j++)
        val += 2.0 * alpha * scad_penalty_single(beta(j), lambda, gamma);
      if (alpha < 1.0)
        val += lambda * (1.0 - alpha) * dot(beta, beta);
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

static int penalized_cd_rss(Penalty pen, double lambda, double gamma,
                            double alpha, double lambda0, double lambda2,
                            const vec& diag_R, const mat& R, const vec& z,
                            double thr, vec& beta, vec& Rbeta, int maxiter) {
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
      // Partial correlation: z_j - sum_{k!=j} R[j,k]*beta[k]
      double t = z(j) - Rbeta(j) + diag_R(j) * bj;

      // Apply penalty-specific thresholding
      beta(j) = apply_threshold(pen, t, diag_R(j), l1, l2, gamma, lambda0);

      if (beta(j) == bj) continue;
      double del = beta(j) - bj;
      dlx = std::max(dlx, std::abs(del));
      // Update running Rbeta
      Rbeta += del * R.col(j);
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

static bool l0_swap_round(double lambda0, double l1, double l2,
                          const vec& diag_R, const mat& R, const vec& z,
                          double thr, vec& beta, vec& Rbeta, int maxiter) {
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
    vec riX = z - Rbeta + R.col(i) * beta(i);

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
    double denom = diag_R(max_idx) + 2.0 * l2;
    if (max_corr > denom * std::abs(beta(i)) + l1) {
      // Perform swap: zero out i, set j to optimal
      double old_bi = beta(i);
      beta(i) = 0.0;
      Rbeta -= old_bi * R.col(i);

      double t_j = z(max_idx) - Rbeta(max_idx) + diag_R(max_idx) * 0.0;
      beta(max_idx) = l0_thresh(t_j, lambda0, l1, l2, diag_R(max_idx));
      if (beta(max_idx) != 0.0) {
        Rbeta += beta(max_idx) * R.col(max_idx);
      }

      // Re-run CD from swapped solution.
      // For PEN_L0: inside penalized_cd_rss, l1 = lambda param, l2 = lambda2.
      penalized_cd_rss(PEN_L0, l1, 0.0, 1.0, lambda0, l2,
                       diag_R, R, z, thr, beta, Rbeta, maxiter);
      found_better = true;
      break; // restart from scratch
    }
  }
  return found_better;
}

// ---- Registered entry point ----------------------------------------------

[[cpp11::register]]
cpp11::writable::list penalizedRssRcpp(
    const doubles& zR,
    const list& LD,
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
      vec diag_R = Rb.diag();
      vec z_blk = z.subvec(s, e);
      vec beta_blk = beta.subvec(s, e);
      vec Rbeta_blk = Rb * beta_blk;

      int conv_blk = penalized_cd_rss(pen, lambda(i), gamma, alpha,
                                       lambda0, lambda2,
                                       diag_R, Rb, z_blk, thr,
                                       beta_blk, Rbeta_blk, maxiter);

      // L0 swap phase (per block)
      if (pen == PEN_L0 && maxSwaps > 0) {
        for (int sw = 0; sw < maxSwaps; sw++) {
          bool improved = l0_swap_round(lambda0, lambda(i), lambda2,
                                        diag_R, Rb, z_blk, thr,
                                        beta_blk, Rbeta_blk, maxiter);
          if (!improved) break;
        }
      }

      beta.subvec(s, e) = beta_blk;
      out = std::min(out, conv_blk);
    }

    beta_mat.col(i) = beta;
    conv_vec(i) = out;

    // Compute loss = beta'R beta - 2 z'beta (block-wise)
    double loss = -2.0 * dot(z, beta);
    for (int b = 0; b < n_blocks; b++) {
      const mat& Rb = ld_blocks[b];
      int s = block_start[b];
      int e = block_end[b];
      vec beta_blk = beta.subvec(s, e);
      loss += as_scalar(beta_blk.t() * Rb * beta_blk);
    }
    loss_vec(i) = loss;
    fbeta_vec(i) = loss + compute_penalty(pen, beta, lambda(i), alpha,
                                           gamma, lambda0, lambda2);
  }

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
