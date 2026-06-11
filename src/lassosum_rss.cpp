// lassosum_rss.cpp — Coordinate descent for LASSO on RSS objective
// Ported from lassosum (Mak et al 2017) functions.cpp elnet()/repelnet()
//
// Objective: min_beta  beta'R beta - 2 beta'z + 2 lambda ||beta||_1
// where R is a (possibly pre-shrunk) LD matrix and z = bhat / sqrt(n).

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>

using namespace cpp11;
using namespace arma;

// Single-block coordinate descent — mirrors lassosum elnet()
// Tracks Rbeta = R * beta as a running sum (analogous to yhat = X * beta
// in the genotype version). This avoids recomputing the full dot product
// each iteration and matches the original algorithm exactly.
static int elnet_rss(double lambda1, const vec& diag_R,
                     const mat& R, const vec& z,
                     double thr, vec& beta, vec& Rbeta,
                     int maxiter) {
  int p = z.n_elem;
  double dlx, del, t, bj;

  int conv = 0;
  for (int k = 0; k < maxiter; k++) {
    dlx = 0.0;
    for (int j = 0; j < p; j++) {
      bj = beta(j);
      // t = z(j) - sum_{k != j} R[j,k]*beta[k]
      //   = z(j) - (Rbeta(j) - R[j,j]*beta(j))
      t = z(j) - Rbeta(j) + diag_R(j) * bj;
      // Soft-thresholding (same as lassosum line 438)
      beta(j) = 0.0;
      if (std::abs(t) - lambda1 > 0.0)
        beta(j) = std::copysign(std::abs(t) - lambda1, t) / diag_R(j);
      if (beta(j) == bj) continue;
      del = beta(j) - bj;
      dlx = std::max(dlx, std::abs(del));
      // Update running Rbeta (analogous to yhat += del * X.col(j))
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

[[cpp11::register]]
cpp11::writable::list lassosumRssRcpp(const doubles& zR,
                                        const list& LD,
                                        const doubles& lambdaR,
                                        double thr,
                                        int maxiter) {
  vec z = as_Col(zR);
  vec lambda = as_Col(lambdaR);

  // Cache LD blocks once (avoid re-copying from R on every lambda iteration)
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

  // Working beta vector — warm-started across lambda path
  vec beta(p, fill::zeros);

  for (int i = 0; i < nlambda; i++) {
    // Block-wise coordinate descent — mirrors lassosum repelnet()
    int out = 1;
    for (int b = 0; b < n_blocks; b++) {
      const mat& Rb = ld_blocks[b];
      int s = block_start[b];
      int e = block_end[b];
      vec diag_R = Rb.diag();
      vec z_blk = z.subvec(s, e);
      vec beta_blk = beta.subvec(s, e);
      vec Rbeta_blk = Rb * beta_blk;

      int conv_blk = elnet_rss(lambda(i), diag_R, Rb, z_blk,
                                thr, beta_blk, Rbeta_blk, maxiter);
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
    fbeta_vec(i) = loss + 2.0 * lambda(i) * sum(abs(beta));
  }

  using namespace cpp11::literals;
  writable::list result({
    "beta"_nm = as_doubles_matrix(beta_mat),
    "lambda"_nm = as_doubles(lambda),
    "conv"_nm = as_integers(conv_vec),
    "loss"_nm = as_doubles(loss_vec),
    "fbeta"_nm = as_doubles(fbeta_vec)
  });

  return result;
}
