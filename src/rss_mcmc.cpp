// rss_mcmc.cpp -- SDPR's likelihood and chain, and the cpp11 entry points of
// every summary-statistics sampler in rss_mcmc.h.
//
// Derived from PRS-CS, SDPR, qgg, hibayes and bigsnpr, each under its own
// license; rss_mcmc.h lists what each part comes from, and inst/COPYRIGHTS
// has the notices.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "rss_mcmc.h"

// ==== SDPR: likelihood and chain ================================================

// SDPR's likelihood and chain. The "Original:" notes point to SDPR's
// mcmc.cpp. The original used GSL and x86 SSE intrinsics. Here GSL's BLAS
// and Cholesky calls are Armadillo products, solves and chol(), its random
// variates are <random>'s (and beta_distribution, a ratio of gammas), and
// the SSE log/exp of sample_assignment() are Armadillo's vectorized float
// log/exp. The sampler reads each LD block through the shared block types:
// the dense matrices SDPR itself precomputes (sdpr_dense_blocks()) or their
// low-rank form from a factor of the LD (sdpr_factor_block()). The original
// ran each block's assignment step on a thread pool sharing one random
// number generator; the blocks run in turn here, as the original does with
// one thread.

namespace {

inline double square(double x) { return x * x; }

// The precision of eta's prior, N(1, 1 / eta_prior_precision) (sample_eta()).
constexpr double eta_prior_precision = 1.0;

// The sampler's state over `blocks`.
template <class LdBlock>
class SdprChain {
public:
    SdprChain(const std::vector<SdprBlock<LdBlock>> &blocks, size_t M,
              double a0k, double b0k, unsigned int seed,
              const ChainOptions &opts, double gap0);

    // Original: mcmc.cpp lines 434-513, one iteration of the main loop.
    void sweep() {
        sample_sigma2();
        for (size_t j = 0; j < blocks_.size(); j++) {
            calc_b(j);
        }
        for (size_t j = 0; j < blocks_.size(); j++) {
            sample_assignment(j);
        }
        update_suffstats();
        sample_V();
        update_p();
        sample_alpha();
        eta_beta = eta;
        for (size_t j = 0; j < blocks_.size(); j++) {
            sample_beta(j);
        }
        sample_eta();
    }

    // The reported effects are eta beta, and the heritability eta^2 beta'R beta.
    // Rao-Blackwellized, an effect is the conditional mean of beta given the
    // rest of the state when it was drawn, times eta at that time.
    void keep(int n_kept) {
        if (check_divergence(est_, opts_, eta * beta, gap0_)) {
            return;
        }
        compute_h2();
        est_.h2 += h2 * square(eta) / n_kept;
        keep_effects(est_, opts_, eta * beta, eta_beta * beta_mean, n_kept);
        if (opts_.track) {
            for (size_t j = 0; j < blocks_.size(); j++) {
                arma::span blk = layout_.span(j);
                est_.corr(blk) += blocks_[j].R.times(eta * beta(blk)) / n_kept;
            }
        }
    }

    void report(int itr, int /* n_iter */) {
        compute_h2();
        Rprintf("%d iter. h2: %g max beta: %g\n",
                itr, h2 * square(eta), arma::max(beta) * eta);
    }

    bool stopped() const { return est_.diverged; }

    const SdprEstimates &estimates() const { return est_; }

private:
    const std::vector<SdprBlock<LdBlock>> &blocks_;
    BlockLayout layout_;
    size_t M, n_snp;
    double a0k;
    double b0k;
    std::mt19937 r;

    double alpha;
    double eta;
    double h2;
    double N;
    arma::vec beta;
    arma::vec b;
    std::vector<int> cls_assgn;
    std::vector<double> V;
    std::vector<double> p;
    std::vector<double> log_p;
    std::vector<double> cluster_var;
    std::vector<unsigned> suff_stats;
    std::vector<double> sumsq;
    // Each block's terms of the eta update, from sample_beta().
    std::vector<double> num;
    std::vector<double> denom;
    ChainOptions opts_;
    double gap0_;
    arma::vec beta_mean;  // the conditional means of beta, for Rao-Blackwell
    double eta_beta = 1.0;  // eta when beta was drawn
    SdprEstimates est_;

    void sample_sigma2();
    void calc_b(size_t j);
    void sample_assignment(size_t j);
    void update_suffstats();
    void sample_V();
    void update_p();
    void sample_alpha();
    void sample_beta(size_t j);
    void compute_h2();
    void sample_eta();
};

template <class LdBlock>
SdprChain<LdBlock>::SdprChain(const std::vector<SdprBlock<LdBlock>> &blocks,
                              size_t M, double a0k, double b0k,
                              unsigned int seed, const ChainOptions &opts,
                              double gap0)
    : blocks_(blocks), layout_(block_sizes(blocks)), M(M), a0k(a0k),
      b0k(b0k), r(seed), opts_(opts), gap0_(gap0) {
    n_snp = 0;
    for (size_t j = 0; j < layout_.n_blocks(); j++) {
        n_snp += layout_.size(j);
    }
    // N (sz) is absorbed into A and B (SDPR, May 20 2021), so it is 1.
    N = 1.0;
    alpha = 1;
    eta = 1;
    h2 = 0;
    beta = arma::vec(n_snp, arma::fill::zeros);
    b = arma::vec(n_snp, arma::fill::zeros);
    p.assign(M, 1.0 / M);
    log_p.assign(M, 0);
    for (size_t i = 0; i < M; i++) {
        log_p[i] = std::log(static_cast<float>(p[i] + 1e-40));
    }

    cluster_var.assign(M, 0.0);
    suff_stats.assign(M, 0);
    sumsq.assign(M, 0.0);
    V.assign(M, 0.0);
    num.assign(blocks.size(), 0.0);
    denom.assign(blocks.size(), 0.0);
    // Every variant starts in the null component (k = 0). The original
    // started from a uniform draw over 0..M-1, which makes nearly every
    // variant causal, so the first sample_beta() factors a matrix the size
    // of the block; starting from null lets the sampler find the causal
    // assignments.
    cls_assgn.assign(n_snp, 0);
    est_.beta = arma::vec(n_snp, arma::fill::zeros);
    est_.corr = arma::vec(n_snp, arma::fill::zeros);
    beta_mean = arma::vec(n_snp, arma::fill::zeros);
    update_suffstats();
}

// ---------------------------------------------------------------------------
// sample_sigma2: Sample cluster variances from inverse-gamma posterior
//   var_k ~ InvGamma(suff_stats[k]/2 + a0k, sumsq[k]/2 + b0k)
// Original: mcmc.cpp lines 25-40
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_sigma2() {
    std::gamma_distribution<double> dist;
    for (size_t i = 1; i < M; i++) {
        double a = suff_stats[i] / 2.0 + a0k;
        double b = 1.0 / (sumsq[i] / 2.0 + b0k);
        dist = std::gamma_distribution<double>(a, b);
        cluster_var[i] = 1.0 / dist(r);
        if (std::isinf(cluster_var[i])) {
            cluster_var[i] = 1e5;
            REprintf("Cluster variance is infinite.\n");
        }
        else if (cluster_var[i] == 0) {
            cluster_var[i] = 1e-10;
            REprintf("Cluster variance is zero.\n");
        }
    }
}

// ---------------------------------------------------------------------------
// calc_b: Compute the b vector for block j
//   b_j = eta^2 * (diag(B)*beta - B*beta) + eta * A^T * beta_mrg
// Original: mcmc.cpp lines 42-62
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::calc_b(size_t j) {
    const SdprBlock<LdBlock> &blk = blocks_[j];
    arma::vec beta_j = beta(layout_.span(j));
    b(layout_.span(j)) = eta * eta * (blk.B.diag() % beta_j - blk.B.times(beta_j))
                         + eta * blk.calc_b_tmp;
}

// ---------------------------------------------------------------------------
// sample_assignment: Sample cluster assignments for each SNP in block j
//   For each SNP i, compute log P(z_i = k) for k = 0..M-1, then sample.
//
// Original: mcmc.cpp lines 64-165, which vectorized log/exp over the M
// clusters with SSE intrinsics in single precision; Armadillo's float
// vectors do the same.
//
// Math for cluster k >= 1:
//   C_k = eta^2 * N * B[i,i] * var_k + 1
//   log P(z_i=k) = -0.5*log(C_k) + log(p_k) + (N*b_i)^2 * var_k / (2*C_k)
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_assignment(size_t j) {
    size_t start_i   = layout_.start(j);
    size_t n_snp_blk = layout_.size(j);
    const arma::vec &diag_B = blocks_[j].B.diag();

    std::uniform_real_distribution<float> unif(0.0f, 1.0f);

    // N = 1.0 after May 21 2021 (absorbed into A, B in solve_ldmat)
    float C = static_cast<float>(eta * eta * N);

    arma::fvec cvar(M);
    arma::fvec log_p_fvec(M);
    for (size_t k = 0; k < M; k++) {
        cvar(k)       = static_cast<float>(cluster_var[k]);
        log_p_fvec(k) = static_cast<float>(log_p[k]);
    }

    for (size_t i = 0; i < n_snp_blk; i++) {
        float Bjj_i = static_cast<float>(diag_B(i));
        float bj_i  = static_cast<float>(b(start_i + i));
        float rnd_i = unif(r);

        arma::fvec prob(M);
        prob(0) = log_p_fvec(0);

        arma::fvec Ck = C * Bjj_i * cvar.subvec(1, M - 1) + 1.0f;

        prob.subvec(1, M - 1) = -0.5f * arma::log(Ck) + log_p_fvec.subvec(1, M - 1)
            + square(N * bj_i) * cvar.subvec(1, M - 1) / (2.0f * Ck);

        // Log-sum-exp for numerical stability
        float max_elem    = prob.max();
        float log_exp_sum = max_elem
            + std::log(arma::accu(arma::exp(prob - max_elem)));

        // Categorical sampling via inverse CDF
        // Original: mcmc.cpp lines 155-163
        cls_assgn[i + start_i] = M - 1;
        for (size_t k = 0; k < M - 1; k++) {
            rnd_i -= std::exp(prob(k) - log_exp_sum);
            if (rnd_i < 0) {
                cls_assgn[i + start_i] = k;
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// update_suffstats: Count SNPs per cluster and sum of squared effects
// Original: mcmc.cpp lines 167-175
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::update_suffstats() {
    std::fill(suff_stats.begin(), suff_stats.end(), 0);
    std::fill(sumsq.begin(), sumsq.end(), 0.0);
    for (size_t i = 0; i < n_snp; i++) {
        suff_stats[cls_assgn[i]]++;
        double tmp = beta(i);
        sumsq[cls_assgn[i]] += square(tmp);
    }
}

// ---------------------------------------------------------------------------
// sample_V: Sample stick-breaking fractions from Beta posterior
//   V[k] ~ Beta(1 + suff_stats[k], alpha + sum_{k'>k} suff_stats[k'])
// Original: mcmc.cpp lines 177-189
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_V() {
    std::vector<double> a(M - 1);
    a[M - 2] = suff_stats[M - 1];
    for (int i = static_cast<int>(M) - 3; i >= 0; i--) {
        a[i] = suff_stats[i + 1] + a[i + 1];
    }

    for (size_t i = 0; i < M - 1; i++) {
        beta_distribution dist(1 + suff_stats[i], alpha + a[i]);
        V[i] = dist(r);
    }
    V[M - 1] = 1;
}

// ---------------------------------------------------------------------------
// update_p: Convert stick-breaking fractions V to cluster probabilities p
//   p[0] = V[0]; p[k] = V[k] * prod_{j<k}(1-V[j])
// Original: mcmc.cpp lines 191-222
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::update_p() {
    std::vector<double> cumprod(M - 1);
    cumprod[0] = 1 - V[0];
    for (size_t i = 1; i < M - 1; i++) {
        cumprod[i] = cumprod[i - 1] * (1 - V[i]);
        if (V[i] == 1) {
            std::fill(cumprod.begin() + i + 1, cumprod.end(), 0.0);
            break;
        }
    }

    p[0] = V[0];
    for (size_t i = 1; i < M - 1; i++) {
        p[i] = cumprod[i - 1] * V[i];
    }

    double sum = std::accumulate(p.begin(), p.end() - 1, 0.0);
    p[M - 1] = (1 - sum > 0) ? (1 - sum) : 0;

    for (size_t i = 0; i < M; i++) {
        log_p[i] = std::log(static_cast<float>(p[i]) + 1e-40f);
    }
}

// ---------------------------------------------------------------------------
// sample_alpha: Sample Dirichlet process concentration parameter
//   alpha ~ Gamma(0.1 + m - 1, 1/(0.1 - sum(log(1-V))))
// Original: mcmc.cpp lines 224-237
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_alpha() {
    double sum = 0, m = 0;
    for (size_t i = 0; i < M; i++) {
        if (V[i] != 1) {
            sum += std::log(1 - V[i]);
            m++;
        }
    }
    if (m == 0) m = 1;

    std::gamma_distribution<double> dist(0.1 + m - 1, 1.0 / (0.1 - sum));
    alpha = dist(r);
}

// ---------------------------------------------------------------------------
// sample_beta: Sample effect sizes for causal SNPs in block j
//   Precision = eta^2*N*B_gamma + Sigma_0^{-1}, mean Precision^-1 A_vec,
//   with A_vec = N*eta*A^T*beta_mrg over the causal SNPs (gaussian_draw()).
//   The block's other effects are zero, so eta's terms are
//   num = N beta'A^T beta_mrg and denom = N beta'B beta.
//
// Original: mcmc.cpp lines 239-329. The original drew a lone causal SNP in
// closed form from its b_j instead, which carries the residuals of the
// previous effects of SNPs no longer causal, and took num = b_j beta_j,
// which is eta times the num above; sample_eta() then centres eta near
// eta^2, so on windows with one causal SNP eta, and with it the effects and
// h2, ran away. A lone causal SNP takes the joint draw here, which draws
// from the same random stream.
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_beta(size_t j) {
    const SdprBlock<LdBlock> &blk = blocks_[j];
    size_t start_i = layout_.start(j);
    size_t end_i   = start_i + layout_.size(j);

    // Build causal list (SNPs not in null cluster)
    std::vector<size_t> causal_list;
    for (size_t i = start_i; i < end_i; i++) {
        if (cls_assgn[i] != 0) {
            causal_list.push_back(i);
        }
    }

    // Zero out this block's betas
    beta(layout_.span(j)).zeros();
    beta_mean(layout_.span(j)).zeros();

    if (causal_list.empty()) {
        num[j]   = 0;
        denom[j] = 0;
        return;
    }

    size_t nc = causal_list.size();
    arma::uvec causal_idx(nc);
    arma::vec var_c(nc);
    for (size_t i = 0; i < nc; i++) {
        causal_idx(i) = causal_list[i] - start_i;
        var_c(i) = cluster_var[cls_assgn[causal_list[i]]];
    }

    // Original: mcmc.cpp lines 273-289
    double C_coeff = square(eta) * N;
    arma::vec A_vec = N * eta * blk.calc_b_tmp(causal_idx);
    arma::mat B_gamma = C_coeff * blk.B.sub(causal_idx);
    B_gamma.diag() += 1.0 / var_c;

    std::normal_distribution<double> dist(0.0, 1.0);
    arma::vec z = standard_normals(nc, dist, r);
    arma::vec mean_c;
    arma::vec beta_c = opts_.rao_blackwell ? gaussian_draw(B_gamma, A_vec, z, 1.0, mean_c)
                                           : gaussian_draw(B_gamma, A_vec, z);

    // Compute eta-related terms for sample_eta(): B_gamma's diagonal back to
    // eta^2 * N * B[i,i], without 1/var_k.
    // Original: mcmc.cpp lines 312-321
    B_gamma.diag() = C_coeff * blk.B.diag()(causal_idx);
    num[j] = arma::dot(A_vec, beta_c) / eta;
    arma::vec tmp = B_gamma * beta_c;
    denom[j] = arma::dot(beta_c, tmp) / square(eta);

    for (size_t i = 0; i < nc; i++) {
        beta(causal_list[i]) = beta_c(i);
        if (opts_.rao_blackwell) {
            beta_mean(causal_list[i]) = mean_c(i);
        }
    }
}

// ---------------------------------------------------------------------------
// compute_h2: Compute heritability h2 = sum_j beta_j^T * R_j * beta_j
// Original: mcmc.cpp lines 331-343
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::compute_h2() {
    h2 = 0;
    for (size_t j = 0; j < blocks_.size(); j++) {
        arma::vec beta_j = beta(layout_.span(j));
        h2 += blocks_[j].R.quad(beta_j);
    }
}

// ---------------------------------------------------------------------------
// sample_eta: Sample global scaling parameter
//   With the prior eta ~ N(1, 1 / eta_prior_precision),
//   eta ~ N((num_sum + prec) / (denom_sum + prec), 1 / (denom_sum + prec)).
// Original: mcmc.cpp lines 345-352, which had no prior on eta beyond the
// guard denom_sum += 1e-6, in effect N(0, 1e6). With no causal SNP, num and
// denom vanish and eta was drawn from that prior: an |eta| in the hundreds
// makes every SNP's inclusion unlikely, so a window could stay null for
// good.
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_eta() {
    double num_sum   = std::accumulate(num.begin(), num.end(), 0.0);
    double denom_sum = std::accumulate(denom.begin(), denom.end(), 0.0);
    double precision = denom_sum + eta_prior_precision;

    std::normal_distribution<double> dist(
        (num_sum + eta_prior_precision) / precision,
        std::sqrt(1.0 / precision));
    eta = dist(r);
}

// ---------------------------------------------------------------------------
// sdpr_llk2_matrix: The optLlk = 2 matrix B2 from one block's LD
//   Each pair's LD is scaled by min(N_i, N_j) / (1.1 * N_i * N_j) and zeroed
//   between arrays 1 and 2; the diagonal is then set to 1 / N_i, raised by
//   1.1 times the scaled matrix's most negative eigenvalue, if any, to keep
//   B2 positive definite.
//
// Original: mcmc.cpp lines 368-402
// ---------------------------------------------------------------------------
arma::mat sdpr_llk2_matrix(arma::mat B, const SdprData &dat, size_t start) {
    size_t blk_size = B.n_rows;
    for (size_t j = 0; j < blk_size; j++) {
        for (size_t k = 0; k < blk_size; k++) {
            size_t idx1 = j + start;
            size_t idx2 = k + start;
            if ((dat.array[idx1] == 1 && dat.array[idx2] == 2) ||
                (dat.array[idx1] == 2 && dat.array[idx2] == 1)) {
                B(j, k) = 0;
            } else {
                B(j, k) *= std::min(dat.sz[idx1], dat.sz[idx2])
                           / (1.1 * dat.sz[idx1] * dat.sz[idx2]);
            }
        }
    }

    // Force positive definite via eigenvalue correction. Only the smallest
    // eigenvalue is read, so no eigenvectors are computed.
    double eval_min = arma::eig_sym(B).min();

    // Restore exact symmetry before the diagonal is reset
    B = arma::symmatu(B);

    // Add to diagonal to ensure PD
    for (size_t j = 0; j < blk_size; j++) {
        double diag_val = 1.0 / dat.sz[j + start];
        if (eval_min < 0) diag_val -= 1.1 * eval_min;
        B(j, j) = diag_val;
    }
    return B;
}

// ---------------------------------------------------------------------------
// sdpr_factor_root: R = Z Z' for the factor F (m x p) of R = F'F / (m - 1)
//   Z (p x r) has orthogonal columns with Z'Z = diag(lambda), lambda the
//   nonzero eigenvalues of R, from the eigendecomposition of the smaller of
//   G G' and G'G (G = F / sqrt(m - 1)). Eigenvalues at the rounding level of
//   that decomposition are dropped.
// ---------------------------------------------------------------------------
arma::mat sdpr_factor_root(const arma::mat &F) {
    arma::mat G = F / std::sqrt(static_cast<double>(F.n_rows) - 1.0);
    bool wide = G.n_rows <= G.n_cols;
    arma::mat gram = wide ? arma::mat(G * G.t()) : arma::mat(G.t() * G);
    arma::vec lam;
    arma::mat Q;
    arma::eig_sym(lam, Q, gram);
    double tol = lam.is_empty() ? 0.0
        : lam.max() * gram.n_rows * arma::datum::eps;
    arma::uvec keep = arma::find(lam > tol);
    if (wide) {
        return G.t() * Q.cols(keep);
    }
    return Q.cols(keep) * arma::diagmat(arma::sqrt(lam(keep)));
}

// ---------------------------------------------------------------------------
// sdpr_llk2_solve: B2^{-1} rhs for the block R = Z Z', with B2 as
// sdpr_llk2_matrix() builds it, without forming B2 while it has a narrow
// low-rank part
//   B2 = T - diag(T) + D0 with T = R .* S .* Mask, D0 = diag(1 / N) + fix.
//   S_ij = 1 / (1.1 max(N_i, N_j)) is a sum over the K distinct sample sizes
//   n_1 < ... < n_K,
//     S = sum_l Delta_l 1_l 1_l' / 1.1,   1_l = [N_i <= n_l],
//     Delta_l = 1 / n_l - 1 / n_(l+1)  (1 / n_K for l = K),
//   so R .* S = U_S U_S' with U_S = [sqrt(Delta_l / 1.1) diag(1_l) Z]_l.
//   The array mask keeps the pairs within arrays {0, 1} and within {0, 2},
//   which counts the pairs within array 0 twice, so T = U Sigma U' with
//   U = [P01 U_S, P02 U_S, P0 U_S] and Sigma = diag(I, I, -I).
//   While U is narrower than p, T is singular, its most negative eigenvalue
//   is that of the small (U'U)^(1/2) Sigma (U'U)^(1/2), and the Woodbury
//   identity solves B2 = U Sigma U' + D through an inner matrix as wide as U.
//   Otherwise B2 is formed exactly as the dense route forms it.
// ---------------------------------------------------------------------------
arma::mat sdpr_llk2_solve(const SdprData &dat, const arma::mat &Z,
                          const arma::mat &rhs) {
    size_t p = Z.n_rows;
    size_t r = Z.n_cols;
    arma::vec sz(dat.sz);
    arma::ivec arr = arma::conv_to<arma::ivec>::from(dat.array);
    arma::vec levels = arma::unique(sz);
    size_t K = levels.n_elem;
    bool masked = arma::any(arr == 1) && arma::any(arr == 2);
    size_t width = (masked ? 3 : 1) * K * r;

    if (width >= p) {
        arma::mat L_chol = arma::chol(sdpr_llk2_matrix(Z * Z.t(), dat, 0),
                                      "lower");
        arma::mat Y = arma::solve(arma::trimatl(L_chol), rhs);
        return arma::solve(arma::trimatu(L_chol.t()), Y);
    }

    arma::mat U_S(p, K * r);
    for (size_t l = 0; l < K; l++) {
        double next = (l + 1 < K) ? 1.0 / levels(l + 1) : 0.0;
        double delta = 1.0 / levels(l) - next;
        arma::vec in_l = arma::conv_to<arma::vec>::from(sz <= levels(l));
        U_S.cols(l * r, (l + 1) * r - 1) =
            std::sqrt(delta / 1.1) * (Z.each_col() % in_l);
    }

    arma::mat U = U_S;
    arma::vec sigma(width, arma::fill::ones);
    // Unmasked, T = U_S U_S' is positive semidefinite (and singular), so its
    // smallest eigenvalue is 0 and B2 needs no fix.
    double eval_min = 0.0;
    if (masked) {
        arma::vec p01 = arma::conv_to<arma::vec>::from(arr != 2);
        arma::vec p02 = arma::conv_to<arma::vec>::from(arr != 1);
        arma::vec p0  = arma::conv_to<arma::vec>::from(arr == 0);
        U = arma::join_rows(U_S.each_col() % p01, U_S.each_col() % p02,
                            U_S.each_col() % p0);
        sigma.tail(K * r).fill(-1.0);
        arma::vec omega;
        arma::mat P;
        arma::eig_sym(omega, P, U.t() * U);
        arma::mat root = P
            * arma::diagmat(arma::sqrt(arma::clamp(omega, 0.0,
                                                   arma::datum::inf)))
            * P.t();
        arma::vec small = arma::eig_sym(root * arma::diagmat(sigma) * root);
        eval_min = std::min(0.0, small.min());
    }

    // D: 1 / N_i in place of T's own diagonal R_ii / (1.1 N_i), plus the fix.
    arma::vec d = 1.0 / sz - arma::sum(arma::square(Z), 1) / (1.1 * sz);
    if (eval_min < 0) d -= 1.1 * eval_min;

    arma::mat UD    = U.each_col() / d;
    arma::mat D_rhs = rhs.each_col() / d;
    arma::mat inner = arma::diagmat(sigma) + U.t() * UD;
    return D_rhs - UD * arma::solve(inner, U.t() * D_rhs);
}

}  // namespace

// ---------------------------------------------------------------------------
// sdpr_dense_blocks: Precompute LD-derived matrices B and A' beta_mrg for
// each block
//   opt_llk=1: B0 = R + aI; A = B0^{-1} R * sz; B = R * A
//   opt_llk=2: B0 = sdpr_llk2_matrix(R); A = B0^{-1} R; B = R * A
//
// Original: mcmc.cpp lines 354-431 (solve_ldmat). The original also kept A
// and R .* R per block, which nothing reads after this, so neither is kept
// here.
// ---------------------------------------------------------------------------
std::vector<SdprBlock<DenseLdBlock>> sdpr_dense_blocks(const SdprData &dat,
                                                       std::vector<arma::mat> ref_ld_mat,
                                                       double a, unsigned sz,
                                                       int opt_llk) {
    std::vector<SdprBlock<DenseLdBlock>> blocks;
    size_t start = 0;
    for (size_t i = 0; i < ref_ld_mat.size(); i++) {
        const arma::mat &R = ref_ld_mat[i];
        size_t blk_size = R.n_rows;

        arma::mat B0;
        if (opt_llk == 1) {
            // B0 = R + a*I  (diagonal shrinkage)
            B0 = R;
            B0.diag() += a;
        } else {
            B0 = sdpr_llk2_matrix(R, dat, start);
        }

        // A = B0^{-1} * R  via forward then back-solve on B0 = L L'
        arma::mat L_chol = arma::chol(B0, "lower");
        arma::mat A = arma::solve(arma::trimatl(L_chol), R);
        A = arma::solve(arma::trimatu(L_chol.t()), A);

        // Scale by sample size for opt_llk=1
        if (opt_llk == 1) {
            A *= sz;
        }

        arma::vec beta_mrg(blk_size);
        for (size_t j = 0; j < blk_size; j++) {
            beta_mrg(j) = dat.beta_mrg[j + start];
        }

        // B = R * A. The sampler reads it as symmetric (the original through
        // dsymv), so it is symmetrized once here rather than on every product.
        blocks.push_back(SdprBlock<DenseLdBlock>{
            DenseLdBlock(arma::symmatu(R * A)),
            DenseLdBlock(arma::symmatu(R)),
            A.t() * beta_mrg
        });
        start += blk_size;
    }
    return blocks;
}

// ---------------------------------------------------------------------------
// sdpr_factor_block: The one block's B and A' beta_mrg from a factor of R
//   With R = Z Z' = V diag(lambda) V':
//   opt_llk=1: A = n (R + aI)^{-1} R = n V diag(lambda / (lambda + a)) V' and
//              B = R A = W W', W = sqrt(n) V diag(lambda / sqrt(lambda + a)).
//   opt_llk=2: A = B2^{-1} R and B = R A = Z H Z', H = Z' B2^{-1} Z (r x r),
//              so W = Z H^(1/2) and B2^{-1} is needed only on Z and beta_mrg.
//   B and R are held through their factors W' and Z'.
// ---------------------------------------------------------------------------
SdprBlock<FactorLdBlock> sdpr_factor_block(const SdprData &dat, const arma::mat &F,
                                           double a, unsigned sz, int opt_llk) {
    arma::mat Z = sdpr_factor_root(F);
    size_t p = Z.n_rows;
    size_t r = Z.n_cols;
    arma::vec beta_mrg(dat.beta_mrg);

    arma::mat W;
    arma::vec calc_b_tmp;
    if (r == 0) {
        // No LD at all (every variant monomorphic): B and A vanish.
        W.zeros(p, 0);
        calc_b_tmp.zeros(p);
    } else if (opt_llk == 1) {
        double n = static_cast<double>(sz);
        arma::vec lam = arma::sum(arma::square(Z), 0).t();
        arma::mat V = Z.each_row() / arma::sqrt(lam).t();
        W = std::sqrt(n) * (V.each_row() % (lam / arma::sqrt(lam + a)).t());
        calc_b_tmp = n * (V * ((lam / (lam + a)) % (V.t() * beta_mrg)));
    } else {
        arma::mat Y = sdpr_llk2_solve(dat, Z, arma::join_rows(Z, beta_mrg));
        arma::mat H = Z.t() * Y.head_cols(r);
        H = 0.5 * (H + H.t());
        arma::vec h;
        arma::mat Q;
        arma::eig_sym(h, Q, H);
        W = Z * Q
            * arma::diagmat(arma::sqrt(arma::clamp(h, 0.0, arma::datum::inf)));
        calc_b_tmp = Z * (Z.t() * Y.col(r));
    }
    return SdprBlock<FactorLdBlock>{
        FactorLdBlock(W.t()), FactorLdBlock(Z.t()), calc_b_tmp
    };
}

template <class LdBlock>
SdprEstimates sdpr_mcmc(const std::vector<SdprBlock<LdBlock>> &blocks, size_t M,
                        double a0k, double b0k, const ChainSchedule &schedule,
                        bool verbose, unsigned int seed, const ChainOptions &opts,
                        double gap0) {
    SdprChain<LdBlock> chain(blocks, M, a0k, b0k, seed, opts, gap0);
    run_chain(chain, schedule, verbose);
    const SdprEstimates &est = chain.estimates();
    if (verbose) {
        Rprintf("h2: %g max: %g\n", est.h2, arma::max(est.beta));
    }
    return est;
}

template SdprEstimates sdpr_mcmc<DenseLdBlock>(
    const std::vector<SdprBlock<DenseLdBlock>> &, size_t, double, double,
    const ChainSchedule &, bool, unsigned int, const ChainOptions &, double);
template SdprEstimates sdpr_mcmc<FactorLdBlock>(
    const std::vector<SdprBlock<FactorLdBlock>> &, size_t, double, double,
    const ChainSchedule &, bool, unsigned int, const ChainOptions &, double);

// ==== Entry points ==============================================================

using namespace cpp11;
using namespace arma;

namespace {

// The seed R passes, or one drawn for NULL.
unsigned int seed_from(sexp seed) {
	if (seed != R_NilValue) {
		return cpp11::as_cpp<unsigned int>(seed);
	}
	return std::random_device{}();
}

ChainOptions chain_options(bool raoBlackwell, bool track) {
	ChainOptions opts;
	opts.rao_blackwell = raoBlackwell;
	opts.track = track;
	return opts;
}

// A factor of the LD as R passes it, refused unless it has a column per
// effect.
mat checked_factor(const doubles_matrix<>& factorR, R_xlen_t n_effects) {
	mat F = as_Mat(factorR);
	if (F.n_cols != static_cast<uword>(n_effects)) {
		cpp11::stop("The LD factor has %d columns but there are %d effects.",
		            (int)F.n_cols, (int)n_effects);
	}
	return F;
}

// What a chain reports to R: its effects, the sampler's own estimates
// (`own`, each a named field), and the fitted correlations and divergence
// every chain tracks.
template <class... Fields>
writable::list chain_report(const ChainEstimates& est, Fields... own) {
	using namespace cpp11::literals;
	return writable::list({
		"betaEst"_nm = as_doubles(est.beta),
		own...,
		"corrEst"_nm = as_doubles(est.corr),
		"diverged"_nm = cpp11::as_sexp(est.diverged)
	});
}

}  // namespace

// ==== PRS-CS ====================================================================

namespace {

// Run the sampler on `blocks` with the arguments as R passes them: a NULL
// `maf` is a vector of zeros, a NULL `phi` is learned from the data, and a
// NULL `seed` draws one.
template <class LdBlock>
cpp11::writable::list prs_cs_run(double a, double b, sexp phi, doubles bhat, sexp maf, int n,
                                 const std::vector<LdBlock>& blocks, int nIter, int nBurnin,
                                 int thin, bool verbose, sexp seed,
                                 const ChainOptions& opts) {
	vec beta_mrg = as_Col(bhat);
	vec maf_vec = maf == R_NilValue ? vec(beta_mrg.n_elem, fill::zeros)
	                                : as_Col(doubles(maf));
	bool learn_phi = phi == R_NilValue;
	double phi_val = learn_phi ? 0.0 : cpp11::as_cpp<double>(phi);
	ChainSchedule schedule(nIter, nBurnin, thin);
	PrsCsEstimates est = prs_cs_mcmc(a, b, learn_phi, phi_val, beta_mrg, maf_vec, n,
	                                 blocks, schedule, verbose, seed_from(seed), opts);
	using namespace cpp11::literals;
	return chain_report(est, "psiEst"_nm = as_doubles(est.psi),
	                    "sigmaEst"_nm = cpp11::as_sexp(est.sigma),
	                    "phiEst"_nm = cpp11::as_sexp(est.phi));
}

}  // namespace

/**
 * @brief PRS-CS from LD matrices.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If NULL, it will be estimated automatically.
 * @param bhat Vector of effect sizes.
 * @param maf Vector of minor allele frequencies. If NULL, it is assumed to be a vector of zeros.
 * @param n Sample size.
 * @param ldBlk List of LD blocks.
 * @param nIter Number of MCMC iterations.
 * @param nBurnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If NULL, one is drawn.
 * @param shrinkCorr The LD shrinkage c: the sampler reads c R + (1 - c) I.
 * @param raoBlackwell Report the blocks' conditional means rather than the
 *   draws' average.
 * @param track Stop a diverging chain and report the mean fitted
 *   correlations, for combining several chains.
 * @return A list containing the posterior estimates.
 */
[[cpp11::register]]
cpp11::writable::list prsCsCpp(double a, double b, sexp phi,
                       doubles bhat, sexp maf,
                       int n, list ldBlk,
                       int nIter, int nBurnin, int thin,
                       bool verbose, sexp seed, double shrinkCorr,
                       bool raoBlackwell, bool track) {
	std::vector<DenseLdBlock> blocks;
	for (int i = 0; i < ldBlk.size(); ++i) {
		mat S = as_Mat(doubles_matrix<>(ldBlk[i]));
		if (shrinkCorr != 1.0) {
			S *= shrinkCorr;
			S.diag() += 1.0 - shrinkCorr;
		}
		blocks.emplace_back(S);
	}
	return prs_cs_run(a, b, phi, bhat, maf, n, blocks, nIter, nBurnin, thin, verbose, seed,
	                  chain_options(raoBlackwell, track));
}

/**
 * @brief PRS-CS from a factor F of one LD block, R = F'F / (m - 1), never forming R.
 *
 * The arguments are prsCsCpp()'s with `factorR` in place of the LD blocks. The chain
 * targets the same posterior as prsCsCpp() on R but draws differently (FactorLdBlock).
 */
[[cpp11::register]]
cpp11::writable::list prsCsFactorCpp(double a, double b, sexp phi,
                             doubles bhat, sexp maf,
                             int n, const doubles_matrix<>& factorR,
                             int nIter, int nBurnin, int thin,
                             bool verbose, sexp seed, double shrinkCorr,
                             bool raoBlackwell, bool track) {
	mat F = checked_factor(factorR, bhat.size());
	mat W = F / std::sqrt(static_cast<double>(F.n_rows) - 1.0);
	// c R + (1 - c) I: the factor of c R, and the ridge 1 - c.
	if (shrinkCorr != 1.0) {
		W *= std::sqrt(shrinkCorr);
	}
	std::vector<FactorLdBlock> blocks;
	blocks.emplace_back(W, 1.0 - shrinkCorr);
	return prs_cs_run(a, b, phi, bhat, maf, n, blocks, nIter, nBurnin, thin, verbose, seed,
	                  chain_options(raoBlackwell, track));
}

// ==== SDPR ======================================================================

namespace {

// The per-variant data as R passes it: the marginal betas deflated by `c`, a
// NULL perVariantSampleSize is n for every variant, and a NULL array is
// array 1 for every variant.
SdprData sdpr_data(const doubles& bhatR, int n, sexp perVariantSampleSize, sexp array,
                   double c) {
	SdprData dat;
	dat.beta_mrg = cpp11::as_cpp<std::vector<double>>(bhatR);
	// Deflation correction
	for (size_t i = 0; i < dat.beta_mrg.size(); i++) {
		dat.beta_mrg[i] /= c;
	}
	if (perVariantSampleSize != R_NilValue) {
		dat.sz = cpp11::as_cpp<std::vector<double>>(perVariantSampleSize);
	} else {
		dat.sz = std::vector<double>(dat.beta_mrg.size(), n);
	}
	if (array != R_NilValue) {
		dat.array = cpp11::as_cpp<std::vector<int>>(array);
	} else {
		dat.array = std::vector<int>(dat.beta_mrg.size(), 1);
	}
	return dat;
}

cpp11::writable::list sdpr_output(const SdprEstimates& est) {
	using namespace cpp11::literals;
	return chain_report(est, "h2"_nm = cpp11::as_sexp(est.h2));
}

}  // namespace

// cpp11 interface function
[[cpp11::register]]
cpp11::writable::list sdprCpp(
	const doubles& bhatR,
	const list& LD,
	int n,
	sexp perVariantSampleSize = R_NilValue,
	sexp array = R_NilValue,
	double a = 0.1,
	double c = 1.0,
	int M = 1000,
	double a0k = 0.5,
	double b0k = 0.5,
	int iter = 1000,
	int burn = 200,
	int thin = 5,
	int optLlk = 1,
	bool verbose = true,
	sexp seed = R_NilValue,
	bool raoBlackwell = false,
	bool track = false
	) {
	ChainSchedule schedule(iter, burn, thin);
	std::vector<mat> ref_ld_mat;
	for (int i = 0; i < LD.size(); i++) {
		ref_ld_mat.push_back(as_Mat(doubles_matrix<>(LD[i])));
	}

	SdprData data = sdpr_data(bhatR, n, perVariantSampleSize, array, c);
	std::vector<SdprBlock<DenseLdBlock>> blocks =
		sdpr_dense_blocks(data, std::move(ref_ld_mat), a, n, optLlk);
	return sdpr_output(sdpr_mcmc(blocks, M, a0k, b0k, schedule, verbose, seed_from(seed),
	                             chain_options(raoBlackwell, track),
	                             divergence_gap(as_Col(bhatR))));
}

/**
 * @brief SDPR from a factor F of the LD, R = F'F / (m - 1), never forming R.
 *
 * The arguments are sdprCpp()'s with `factorR` (m x p) in place of the LD
 * blocks; the window is one block. The sampler is SDPR's own, run on the
 * low-rank form of its matrices (sdpr_factor_block()), so it draws as
 * sdprCpp() does on the matrix the factor gives.
 */
[[cpp11::register]]
cpp11::writable::list sdprFactorCpp(
	const doubles& bhatR,
	const doubles_matrix<>& factorR,
	int n,
	sexp perVariantSampleSize = R_NilValue,
	sexp array = R_NilValue,
	double a = 0.1,
	double c = 1.0,
	int M = 1000,
	double a0k = 0.5,
	double b0k = 0.5,
	int iter = 1000,
	int burn = 200,
	int thin = 5,
	int optLlk = 1,
	bool verbose = true,
	sexp seed = R_NilValue,
	bool raoBlackwell = false,
	bool track = false
	) {
	ChainSchedule schedule(iter, burn, thin);
	mat F = checked_factor(factorR, bhatR.size());

	SdprData data = sdpr_data(bhatR, n, perVariantSampleSize, array, c);
	std::vector<SdprBlock<FactorLdBlock>> blocks = {
		sdpr_factor_block(data, F, a, n, optLlk)
	};
	return sdpr_output(sdpr_mcmc(blocks, M, a0k, b0k, schedule, verbose, seed_from(seed),
	                             chain_options(raoBlackwell, track),
	                             divergence_gap(as_Col(bhatR))));
}

// ==== The Bayesian alphabet =====================================================

namespace {

template <class Chain>
cpp11::writable::list alphabet_run(Chain& chain, const ChainSchedule& schedule,
                                   bool verbose) {
	run_chain(chain, schedule, verbose);
	const AlphabetEstimates& est = chain.estimates();
	using namespace cpp11::literals;
	return chain_report(est, "pip"_nm = as_doubles(est.pip),
	                    "pi"_nm = cpp11::as_sexp(est.pi),
	                    "residualVariance"_nm = cpp11::as_sexp(est.ve));
}

// qgg's family of priors for every method but bayesB and Park and Casella's
// bayesL.
QggFormChain::Method qgg_method(const std::string& method) {
	if (method == "bayesN") return QggFormChain::N;
	if (method == "bayesA") return QggFormChain::A;
	if (method == "bayesL") return QggFormChain::L;
	if (method == "bayesC") return QggFormChain::C;
	if (method == "bayesR") return QggFormChain::R;
	cpp11::stop("Unknown Bayesian alphabet method '%s'.", method.c_str());
}

}  // namespace

/**
 * @brief The Bayesian alphabet's posterior means from a factor F of the LD,
 * R = F'F / (k - 1), never forming R.
 *
 * @param bhat The standardized marginal effects.
 * @param factorR The factor F (k x p).
 * @param weights The per-variant weights w: X'X = W^1/2 R W^1/2 and
 *   X'y = w % bhat.
 * @param n The sample size.
 * @param method "bayesA", "bayesB", "bayesC", "bayesL", "bayesN" or "bayesR".
 * @param legarra For bayesL, Legarra et al.'s form (qgg's) rather than Park
 *   and Casella's (hibayes').
 * @param pi For bayesB and bayesC, the starting proportion of non-null
 *   effects.
 * @param h2 The prior guess of the heritability that sets the starting
 *   variances.
 * @param nIter The chain's iterations, burn-in included.
 * @param nBurnin The burn-in.
 * @param thin The thinning interval.
 * @param verbose Whether to report progress.
 * @param seed Seed of the sampler's random number generator.
 * @param shrinkCorr The LD shrinkage c: X'X = W^1/2 (c R + (1 - c) I) W^1/2.
 * @param allowJumpSign Keep a drawn effect whose sign differs from the
 *   current one, rather than setting it to zero.
 * @param raoBlackwell Report the effects' conditional means rather than the
 *   draws' average.
 * @param track Stop a diverging chain and report the mean fitted
 *   correlations, for combining several chains.
 */
[[cpp11::register]]
cpp11::writable::list bayesAlphabetFactorCpp(
	const doubles& bhat, const doubles_matrix<>& factorR, const doubles& weights,
	int n, std::string method, bool legarra, double pi, double h2,
	int nIter, int nBurnin, int thin, bool verbose, int seed,
	double shrinkCorr, bool allowJumpSign, bool raoBlackwell, bool track) {
	ChainSchedule schedule(nIter, nBurnin, thin);
	mat F = checked_factor(factorR, bhat.size());
	AlphabetData data;
	data.w = as_Col(weights);
	data.bhat = as_Col(bhat);
	data.Fw = F.each_row() % arma::sqrt(data.w).t();
	data.wy = data.w % data.bhat;
	data.n = n;
	data.shrink = shrinkCorr;
	data.allow_jump_sign = allowJumpSign;
	data.opts.rao_blackwell = raoBlackwell;
	data.opts.track = track;
	unsigned int s = static_cast<unsigned int>(seed);
	if (method == "bayesB" || (method == "bayesL" && !legarra)) {
		HibayesFormChain chain(data,
		                       method == "bayesB" ? HibayesFormChain::B : HibayesFormChain::L,
		                       pi, h2, s);
		return alphabet_run(chain, schedule, verbose);
	}
	QggFormChain chain(data, qgg_method(method), pi, h2, s);
	return alphabet_run(chain, schedule, verbose);
}

// ==== LDpred2 ===================================================================

/**
 * @brief One LDpred2 chain from a factor F of the LD, R = F'F / (k - 1),
 * never forming R.
 *
 * @param bhat The standardized marginal effects.
 * @param factorR The factor F (k x p).
 * @param n The per-variant sample sizes.
 * @param automatic LDpred2-auto rather than LDpred2-grid.
 * @param p The proportion of causal variants: fixed for the grid, the
 *   starting value for auto.
 * @param h2 The heritability: fixed for the grid, the starting value for
 *   auto.
 * @param sparse For the grid, set an effect to zero when its inclusion
 *   probability falls below p.
 * @param shrinkCorr The LD shrinkage c in R_s = c R + (1 - c) I.
 * @param allowJumpSign Keep a drawn effect whose sign differs from the
 *   current one, rather than setting it to zero.
 * @param pLower,pUpper For auto, the bounds on p.
 * @param meanLd For auto, the mean LD score that divides the counts of p's
 *   Beta full conditional.
 * @param nIter The chain's iterations, burn-in included.
 * @param nBurnin The burn-in.
 * @param verbose Whether to report progress.
 * @param seed Seed of the sampler's random number generator.
 */
[[cpp11::register]]
cpp11::writable::list ldpred2FactorCpp(
	const doubles& bhat, const doubles_matrix<>& factorR, const doubles& n,
	bool automatic, double p, double h2, bool sparse, double shrinkCorr,
	bool allowJumpSign, double pLower, double pUpper, double meanLd,
	int nIter, int nBurnin, bool verbose, int seed) {
	ChainSchedule schedule(nIter, nBurnin, 1);
	mat F = checked_factor(factorR, bhat.size());
	vec b = as_Col(bhat);
	vec nv = as_Col(n);
	if (nv.n_elem != b.n_elem) {
		cpp11::stop("There are %d sample sizes but %d effects.", (int)nv.n_elem,
		            (int)b.n_elem);
	}
	Ldpred2Settings settings;
	settings.automatic = automatic;
	settings.p = p;
	settings.h2 = h2;
	settings.sparse = sparse;
	settings.shrink = shrinkCorr;
	settings.allow_jump_sign = allowJumpSign;
	settings.p_lower = pLower;
	settings.p_upper = pUpper;
	settings.mean_ld = meanLd;
	Ldpred2Chain chain(F, b, nv, settings, static_cast<unsigned int>(seed));
	run_chain(chain, schedule, verbose);
	const Ldpred2Estimates& est = chain.estimates();
	using namespace cpp11::literals;
	return chain_report(
		est, "postp"_nm = as_doubles(est.postp),
		"pathP"_nm = writable::doubles(est.path_p.begin(), est.path_p.end()),
		"pathH2"_nm = writable::doubles(est.path_h2.begin(), est.path_h2.end()));
}
