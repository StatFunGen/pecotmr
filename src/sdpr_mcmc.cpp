// sdpr_mcmc.cpp -- SDPR MCMC sampler using Armadillo
//
// Clean Armadillo port of the original SDPR by Zhou et al.
// (https://github.com/eldronzhou/SDPR), which used GSL + x86 SSE intrinsics.
//
// Translation strategy:
//   GSL gsl_blas_dsymv  -> arma::symmatu(M) * v
//   GSL gsl_blas_dgemv  -> A.t() * v
//   GSL gsl_blas_dgemm  -> A * B
//   GSL gsl_blas_dtrsv  -> arma::solve(trimatl(L), v)
//   GSL gsl_blas_dtrsm  -> arma::solve(trimatl(L), A)
//   GSL gsl_blas_ddot   -> arma::dot(x, y)
//   GSL gsl_blas_daxpy  -> y += alpha * x
//   GSL gsl_linalg_cholesky_decomp1 -> arma::chol(A, "lower")
//   GSL gsl_ran_gamma   -> std::gamma_distribution
//   GSL gsl_ran_beta    -> beta_distribution (ratio of gammas)
//   GSL gsl_ran_ugaussian -> std::normal_distribution(0,1)
//   GSL gsl_rng_uniform -> std::uniform_real_distribution(0,1)
//
// The original SSE intrinsics (log_ps, exp_ps, _mm_max_ps, _mm_hadd_ps)
// in sample_assignment() are replaced with Armadillo vectorized operations
// (arma::log, arma::exp, arma::max, arma::accu) which delegate to
// platform-optimal SIMD (NEON on ARM, SSE/AVX on x86) through the
// underlying BLAS/compiler auto-vectorization.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>
#include <fstream>
#include <numeric>
#include <random>
#include "function_pool.h"
#include "sdpr_mcmc.h"

using std::cout; using std::endl;
using std::thread; using std::ref;
using std::vector; using std::ofstream;
using std::string; using std::min;

#define square(x) ((x)*(x))

// ---------------------------------------------------------------------------
// sample_sigma2: Sample cluster variances from inverse-gamma posterior
//   var_k ~ InvGamma(suff_stats[k]/2 + a0k, sumsq[k]/2 + b0k)
// Original: mcmc.cpp lines 25-40
// ---------------------------------------------------------------------------
void MCMC_state::sample_sigma2() {
    std::gamma_distribution<double> dist;
    for (size_t i = 1; i < M; i++) {
        double a = suff_stats[i] / 2.0 + a0k;
        double b = 1.0 / (sumsq[i] / 2.0 + b0k);
        dist = std::gamma_distribution<double>(a, b);
        cluster_var[i] = 1.0 / dist(r);
        if (std::isinf(cluster_var[i])) {
            cluster_var[i] = 1e5;
            Rcpp::Rcerr << "Cluster variance is infinite." << std::endl;
        }
        else if (cluster_var[i] == 0) {
            cluster_var[i] = 1e-10;
            Rcpp::Rcerr << "Cluster variance is zero." << std::endl;
        }
    }
}

// ---------------------------------------------------------------------------
// calc_b: Compute the b vector for block j
//   b_j = eta^2 * (diag(B)*beta - B*beta) + eta * A^T * beta_mrg
// Original: mcmc.cpp lines 42-62
//   gsl_blas_dsymv(CblasUpper, -eta*eta, B, beta, eta*eta, b)
//   gsl_blas_daxpy(eta, calc_b_tmp, b)
// ---------------------------------------------------------------------------
void MCMC_state::calc_b(size_t j, const mcmc_data &dat,
                        const ldmat_data &ldmat_dat) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i   = dat.boundary[j].second;

    arma::vec beta_j = beta.subvec(start_i, end_i - 1);
    arma::vec diag_B = ldmat_dat.B[j].diag();

    // Original GSL: b = eta^2 * diag(B)*beta;
    //               b = eta^2*b - eta^2 * B*beta  (via dsymv with alpha=-eta^2, beta=eta^2)
    //               b += eta * calc_b_tmp
    arma::vec b_j = eta * eta * (diag_B % beta_j
                    - arma::symmatu(ldmat_dat.B[j]) * beta_j)
                    + eta * ldmat_dat.calc_b_tmp[j];

    b.subvec(start_i, end_i - 1) = b_j;
}

// ---------------------------------------------------------------------------
// sample_assignment: Sample cluster assignments for each SNP in block j
//   For each SNP i, compute log P(z_i = k) for k = 0..M-1, then sample.
//
// Original: mcmc.cpp lines 64-165
//   Used x86 SSE intrinsics (log_ps, exp_ps, _mm_max_ps, _mm_hadd_ps)
//   to vectorize log/exp over M clusters. We use Armadillo vectorized
//   ops instead, which auto-vectorize via NEON (ARM) or SSE/AVX (x86).
//
// Math for cluster k >= 1:
//   C_k = eta^2 * N * B[i,i] * var_k + 1
//   log P(z_i=k) = -0.5*log(C_k) + log(p_k) + (N*b_i)^2 * var_k / (2*C_k)
// ---------------------------------------------------------------------------
void MCMC_state::sample_assignment(size_t j, const mcmc_data &dat,
                                   const ldmat_data &ldmat_dat) {
    size_t start_i   = dat.boundary[j].first;
    size_t end_i     = dat.boundary[j].second;
    size_t n_snp_blk = end_i - start_i;

    std::uniform_real_distribution<float> unif(0.0f, 1.0f);

    // N = 1.0 after May 21 2021 (absorbed into A, B in solve_ldmat)
    float C = static_cast<float>(eta * eta * N);

    // Pre-convert cluster variances and log-probabilities to float vectors
    // for Armadillo vectorized ops (matching original's float precision)
    arma::fvec cvar(M);
    arma::fvec log_p_fvec(M);
    for (size_t k = 0; k < M; k++) {
        cvar(k)       = static_cast<float>(cluster_var[k]);
        log_p_fvec(k) = static_cast<float>(log_p[k]);
    }

    for (size_t i = 0; i < n_snp_blk; i++) {
        float Bjj_i = static_cast<float>(ldmat_dat.B[j](i, i));
        float bj_i  = static_cast<float>(b(start_i + i));
        float rnd_i = unif(r);

        // prob[0] = log_p[0] (null cluster)
        // prob[k] for k>=1: see math above
        arma::fvec prob(M);
        prob(0) = log_p_fvec(0);

        // Vectorized: Ck = C * Bjj * cvar[1..M-1] + 1
        arma::fvec Ck = C * Bjj_i * cvar.subvec(1, M - 1) + 1.0f;

        // prob[k] = -0.5*log(Ck) + log_p[k] + (N*bj)^2 * var_k / (2*Ck)
        // Original: tmp[k] = log(prob[k]); prob[k] = -0.5*tmp[k] + log_p[k] + ...
        prob.subvec(1, M - 1) = -0.5f * arma::log(Ck) + log_p_fvec.subvec(1, M - 1)
            + square(N * bj_i) * cvar.subvec(1, M - 1) / (2.0f * Ck);

        // Log-sum-exp for numerical stability (replaces SSE _mm_max_ps + exp_ps + _mm_hadd_ps)
        float max_elem    = prob.max();
        float log_exp_sum = max_elem
            + std::logf(arma::accu(arma::exp(prob - max_elem)));

        // Categorical sampling via inverse CDF
        // Original: mcmc.cpp lines 155-163
        cls_assgn[i + start_i] = M - 1;
        for (size_t k = 0; k < M - 1; k++) {
            rnd_i -= std::expf(prob(k) - log_exp_sum);
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
void MCMC_state::update_suffstats() {
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
void MCMC_state::sample_V() {
    vector<double> a(M - 1);
    a[M - 2] = suff_stats[M - 1];
    for (int i = M - 3; i >= 0; i--) {
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
void MCMC_state::update_p() {
    vector<double> cumprod(M - 1);
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
        log_p[i] = std::logf(static_cast<float>(p[i]) + 1e-40f);
    }
}

// ---------------------------------------------------------------------------
// sample_alpha: Sample Dirichlet process concentration parameter
//   alpha ~ Gamma(0.1 + m - 1, 1/(0.1 - sum(log(1-V))))
// Original: mcmc.cpp lines 224-237
// ---------------------------------------------------------------------------
void MCMC_state::sample_alpha() {
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
//   Precision = eta^2*N*B_gamma + Sigma_0^{-1}
//   Mean via Cholesky: L*L^T = Precision; mu = L^{-T}*L^{-1}*A_vec
//   Sample: z ~ N(0,I); beta = L^{-T}*(L^{-1}*A_vec + z)
//
// Original: mcmc.cpp lines 239-329
//   Uses gsl_linalg_cholesky_decomp1, gsl_blas_dtrsv for forward/back-solve
// ---------------------------------------------------------------------------
void MCMC_state::sample_beta(size_t j, const mcmc_data &dat,
                             ldmat_data &ldmat_dat) {
    size_t start_i = dat.boundary[j].first;
    size_t end_i   = dat.boundary[j].second;

    // Build causal list (SNPs not in null cluster)
    vector<size_t> causal_list;
    for (size_t i = start_i; i < end_i; i++) {
        if (cls_assgn[i] != 0) {
            causal_list.push_back(i);
        }
    }

    // Zero out this block's betas
    beta.subvec(start_i, end_i - 1).zeros();

    if (causal_list.empty()) {
        ldmat_dat.num[j]   = 0;
        ldmat_dat.denom[j] = 0;
        return;
    }

    // Single causal SNP: closed-form sampling
    // Original: mcmc.cpp lines 259-270
    if (causal_list.size() == 1) {
        double var_k = cluster_var[cls_assgn[causal_list[0]]];
        double bj    = b(causal_list[0]);
        double Bjj   = ldmat_dat.B[j](causal_list[0] - start_i,
                                       causal_list[0] - start_i);
        // C = var_k / (N * var_k * eta^2 * Bjj + 1)
        double C_val = var_k / (N * var_k * square(eta) * Bjj + 1.0);
        std::normal_distribution<double> dist(0.0, std::sqrt(C_val));
        double rv = dist(r) + C_val * N * bj;
        beta(causal_list[0]) = rv;
        ldmat_dat.num[j]   = bj * rv;
        ldmat_dat.denom[j] = square(rv) * Bjj;
        return;
    }

    // Multiple causal SNPs: multivariate normal sampling via Cholesky
    size_t nc = causal_list.size();

    // A_vec = N * eta * A^T * beta_mrg (restricted to causal indices)
    // Original: mcmc.cpp line 279: N*eta*gsl_vector_get(calc_b_tmp, idx)
    arma::vec A_vec(nc);
    double C_coeff = square(eta) * N;  // Original: C = square(eta)*N (line 273)

    // Build precision matrix B_gamma and A_vec
    // Original: mcmc.cpp lines 278-289
    arma::mat B_gamma(nc, nc);
    for (size_t i = 0; i < nc; i++) {
        size_t idx_i = causal_list[i] - start_i;
        A_vec(i) = N * eta * ldmat_dat.calc_b_tmp[j](idx_i);

        for (size_t k = 0; k < nc; k++) {
            size_t idx_k = causal_list[k] - start_i;
            if (i != k) {
                // Off-diagonal: eta^2 * N * B[i,k]
                B_gamma(i, k) = C_coeff * ldmat_dat.B[j](idx_i, idx_k);
            } else {
                // Diagonal: eta^2 * N * B[i,i] + 1/var_k
                B_gamma(i, k) = C_coeff * ldmat_dat.B[j](idx_i, idx_i)
                    + 1.0 / cluster_var[cls_assgn[causal_list[i]]];
            }
        }
    }

    arma::vec A_vec2 = A_vec;  // Save for eta computation

    // Sample z ~ N(0, I)
    arma::vec beta_c(nc);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (size_t i = 0; i < nc; i++) {
        beta_c(i) = dist(r);
    }

    // Cholesky: B_gamma = L * L^T
    // Original: gsl_linalg_cholesky_decomp1(&B.matrix)
    arma::mat L = arma::chol(B_gamma, "lower");

    // mu = L^{-1} * A_vec  (forward-solve)
    // Original: gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &B, A_vec)
    A_vec = arma::solve(arma::trimatl(L), A_vec);

    // beta_c = mu + z ~ N(mu, I)
    beta_c += A_vec;

    // beta_c = L^{-T} * beta_c ~ N(L^{-T}*mu, (L*L^T)^{-1})
    // Original: gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, &B, beta_c)
    beta_c = arma::solve(arma::trimatu(L.t()), beta_c);

    // Compute eta-related terms for sample_eta()
    // Original: mcmc.cpp lines 312-321
    // Restore diagonal of B_gamma to eta^2 * N * B[i,i] (without 1/var_k)
    for (size_t i = 0; i < nc; i++) {
        size_t idx_i = causal_list[i] - start_i;
        B_gamma(i, i) = C_coeff * ldmat_dat.B[j](idx_i, idx_i);
    }

    // num = A_vec2^T * beta_c / eta
    ldmat_dat.num[j] = arma::dot(A_vec2, beta_c) / eta;

    // denom = beta_c^T * B_gamma * beta_c / eta^2
    // Original: gsl_blas_dsymv(CblasUpper, 1.0, &B, beta_c, 0, A_vec)
    arma::vec tmp = arma::symmatu(B_gamma) * beta_c;
    ldmat_dat.denom[j] = arma::dot(beta_c, tmp) / square(eta);

    // Write sampled betas back
    for (size_t i = 0; i < nc; i++) {
        beta(causal_list[i]) = beta_c(i);
    }
}

// ---------------------------------------------------------------------------
// compute_h2: Compute heritability h2 = sum_j beta_j^T * R_j * beta_j
// Original: mcmc.cpp lines 331-343
// ---------------------------------------------------------------------------
void MCMC_state::compute_h2(const mcmc_data &dat) {
    h2 = 0;
    for (size_t j = 0; j < dat.ref_ld_mat.size(); j++) {
        size_t start_i = dat.boundary[j].first;
        size_t end_i   = dat.boundary[j].second;
        arma::vec beta_j = beta.subvec(start_i, end_i - 1);
        // Original: gsl_blas_dsymv(CblasUpper, 1.0, ref_ld_mat, beta, 0, tmp)
        arma::vec tmp = arma::symmatu(dat.ref_ld_mat[j]) * beta_j;
        h2 += arma::dot(tmp, beta_j);
    }
}

// ---------------------------------------------------------------------------
// sample_eta: Sample global scaling parameter
//   eta ~ N(num_sum/denom_sum, 1/denom_sum)
// Original: mcmc.cpp lines 345-352
// ---------------------------------------------------------------------------
void MCMC_state::sample_eta(const ldmat_data &ldmat_dat) {
    double num_sum   = std::accumulate(ldmat_dat.num.begin(),
                                       ldmat_dat.num.end(), 0.0);
    double denom_sum = std::accumulate(ldmat_dat.denom.begin(),
                                       ldmat_dat.denom.end(), 0.0);
    denom_sum += 1e-6;

    std::normal_distribution<double> dist(num_sum / denom_sum,
                                          std::sqrt(1.0 / denom_sum));
    eta = dist(r);
}

// ---------------------------------------------------------------------------
// solve_ldmat: Precompute LD-derived matrices A, B, L for each block
//   opt_llk=1: B = (R + aI); A = B^{-1} R * sz; B = R*A; L = R.*R
//   opt_llk=2: Multi-array with eigenvalue correction
//
// Original: mcmc.cpp lines 354-431
//   Uses gsl_linalg_cholesky_decomp1, gsl_blas_dtrsm for matrix solves
// ---------------------------------------------------------------------------
void solve_ldmat(const mcmc_data &dat, ldmat_data &ldmat_dat,
                 const double a, unsigned sz, int opt_llk) {
    for (size_t i = 0; i < dat.ref_ld_mat.size(); i++) {
        size_t blk_size = dat.boundary[i].second - dat.boundary[i].first;
        arma::mat A = dat.ref_ld_mat[i];
        arma::mat B = dat.ref_ld_mat[i];
        arma::mat L = dat.ref_ld_mat[i];  // Will become R .* R (element-wise square)

        if (opt_llk == 1) {
            // B = R + a*I  (diagonal shrinkage)
            // Original: gsl_vector_add_constant(&diag.vector, a)
            B.diag() += a;
        }
        else {
            // Multi-array: scale by min(N_i, N_j) / (1.1 * N_i * N_j)
            // Original: mcmc.cpp lines 368-402
            for (size_t j = 0; j < blk_size; j++) {
                for (size_t k = 0; k < blk_size; k++) {
                    size_t idx1 = j + dat.boundary[i].first;
                    size_t idx2 = k + dat.boundary[i].first;
                    if ((dat.array[idx1] == 1 && dat.array[idx2] == 2) ||
                        (dat.array[idx1] == 2 && dat.array[idx2] == 1)) {
                        B(j, k) = 0;
                    } else {
                        B(j, k) *= min(dat.sz[idx1], dat.sz[idx2])
                                   / (1.1 * dat.sz[idx1] * dat.sz[idx2]);
                    }
                }
            }

            // Force positive definite via eigenvalue correction
            arma::vec eval;
            arma::mat evec;
            arma::eig_sym(eval, evec, B);
            double eval_min = eval.min();

            // Restore symmetry (eig_sym may modify lower triangle)
            B = arma::symmatu(B);

            // Add to diagonal to ensure PD
            for (size_t j = 0; j < blk_size; j++) {
                double diag_val = 1.0 / dat.sz[j + dat.boundary[i].first];
                if (eval_min < 0) diag_val -= 1.1 * eval_min;
                B(j, j) = diag_val;
            }
        }

        // Cholesky: B = L_chol * L_chol^T
        // Original: gsl_linalg_cholesky_decomp1(B)
        arma::mat L_chol = arma::chol(B, "lower");

        // A = B^{-1} * R  via forward then back-solve
        // Original: gsl_blas_dtrsm(Left, Lower, NoTrans, NonUnit, 1.0, B, A)
        //           gsl_blas_dtrsm(Left, Lower, Trans,   NonUnit, 1.0, B, A)
        A = arma::solve(arma::trimatl(L_chol), A);
        A = arma::solve(arma::trimatu(L_chol.t()), A);

        // Scale by sample size for opt_llk=1
        if (opt_llk == 1) {
            A *= sz;
        }

        // B = R * A  (overwrite B with the product)
        // Original: gsl_blas_dgemm(NoTrans, NoTrans, 1.0, L, A, 0, B)
        // Note: L still holds the original R at this point
        B = L * A;

        // L = R .* R  (element-wise square)
        // Original: gsl_matrix_mul_elements(L, L)
        L %= L;

        // Precompute A^T * beta_mrg for calc_b()
        // Original: gsl_blas_dgemv(CblasTrans, 1.0, A, beta_mrg, 0, b_tmp)
        arma::vec beta_mrg(blk_size);
        for (size_t j = 0; j < blk_size; j++) {
            beta_mrg(j) = dat.beta_mrg[j + dat.boundary[i].first];
        }
        arma::vec b_tmp = A.t() * beta_mrg;

        ldmat_dat.A.push_back(A);
        ldmat_dat.B.push_back(B);
        ldmat_dat.L.push_back(L);
        ldmat_dat.calc_b_tmp.push_back(b_tmp);
        ldmat_dat.beta_mrg.push_back(beta_mrg);
        ldmat_dat.denom.push_back(0);
        ldmat_dat.num.push_back(0);
    }
}

// ---------------------------------------------------------------------------
// mcmc: Main MCMC loop
// Original: mcmc.cpp lines 434-513
// ---------------------------------------------------------------------------
std::unordered_map<std::string, arma::vec> mcmc(
    mcmc_data& data,
    unsigned   sz,
    double     a,
    double     c,
    size_t     M,
    double     a0k,
    double     b0k,
    int        iter,
    int        burn,
    int        thin,
    unsigned   n_threads,
    int        opt_llk,
    bool       verbose,
    unsigned int seed
    ) {

    int n_pst = (iter - burn) / thin;

    ldmat_data ldmat_dat;

    MCMC_state state(data.beta_mrg.size(), M, a0k, b0k, sz, seed);

    // Deflation correction
    for (size_t i = 0; i < data.beta_mrg.size(); i++) {
        data.beta_mrg[i] /= c;
    }

    MCMC_samples samples(data.beta_mrg.size());

    solve_ldmat(data, ldmat_dat, a, sz, opt_llk);
    state.update_suffstats();

    Function_pool func_pool(n_threads);

    for (int j = 1; j < iter + 1; j++) {
        state.sample_sigma2();

        for (size_t i = 0; i < data.ref_ld_mat.size(); i++) {
            state.calc_b(i, data, ldmat_dat);
        }

        // sample_assignment dispatched to thread pool across LD blocks.
        // NOTE: std::mt19937 is not thread-safe. With n_threads > 1 there
        // is a data race on MCMC_state::r. This matches the original SDPR
        // (eldronzhou/SDPR) which has the same issue with gsl_rng. The
        // default n_threads=1 avoids the race. For safe parallelism, each
        // block would need its own RNG seeded from the shared one.
        for (size_t i = 0; i < data.ref_ld_mat.size(); i++) {
            func_pool.push(std::bind(&MCMC_state::sample_assignment,
                                     &state, i, ref(data), ref(ldmat_dat)));
        }
        func_pool.waitFinished();

        state.update_suffstats();
        state.sample_V();
        state.update_p();
        state.sample_alpha();

        for (size_t i = 0; i < data.ref_ld_mat.size(); i++) {
            state.sample_beta(i, data, ldmat_dat);
        }

        state.sample_eta(ldmat_dat);

        // Collect posterior samples
        if ((j > burn) && (j % thin == 0)) {
            state.compute_h2(data);
            samples.h2   += state.h2 * square(state.eta) / n_pst;
            samples.beta += state.eta / n_pst * state.beta;
        }

        if (verbose && j % 100 == 0) {
            state.compute_h2(data);
            Rcpp::Rcout << j << " iter. h2: "
                        << state.h2 * square(state.eta)
                        << " max beta: "
                        << arma::max(state.beta) * state.eta << endl;
        }
    }

    if (verbose) {
        Rcpp::Rcout << "h2: " << samples.h2
                    << " max: " << arma::max(samples.beta) << endl;
    }

    std::unordered_map<std::string, arma::vec> results;
    results["beta"] = samples.beta;
    results["h2"]   = arma::vec(1, arma::fill::value(samples.h2));

    return results;
}
