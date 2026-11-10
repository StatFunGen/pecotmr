// rss_mcmc.cpp -- SDPR's likelihood and chain, for one study or across
// studies as SDPRX, and the cpp11 entry points of every summary-statistics
// sampler in rss_mcmc.h.
//
// Derived from PRS-CS, PRS-CSx, JointPRS, SDPR and SDPRX, each under its own
// license; rss_mcmc.h lists what each part comes from, and inst/COPYRIGHTS
// has the notices.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "rss_mcmc.h"

// ==== SDPR: likelihood and chain ================================================

// SDPR's likelihood and chain, across studies with SDPRX's kinds of
// components (rss_mcmc.h). The "Original:" notes point to SDPR's
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

// The grid a learned correlation is drawn on: 0, 0.01, ..., 0.99
// (sample_rho()).
constexpr int rho_grid_points = 100;
constexpr double rho_grid_step = 0.01;

// A set of studies, bit k for study k.
using StudySet = std::uint64_t;

inline int set_size(StudySet s) {
    return static_cast<int>(std::bitset<64>(s).count());
}

inline bool set_holds(StudySet s, arma::uword k) {
    return ((s >> k) & 1U) != 0;
}

// The first study of a nonempty `s`.
inline arma::uword set_first(StudySet s) {
    arma::uword k = 0;
    while (!set_holds(s, k)) {
        k++;
    }
    return k;
}

// The studies of `s`, in order.
inline arma::uvec set_studies(StudySet s) {
    std::vector<arma::uword> out;
    for (arma::uword k = 0; k < 64; k++) {
        if (set_holds(s, k)) {
            out.push_back(k);
        }
    }
    return arma::conv_to<arma::uvec>::from(out);
}

// The correlations of a set of two or more studies, Sigma (rho's submatrix):
// its symmetric square root and inverse.
struct SetCorr {
    arma::uvec studies;
    arma::mat root;
    arma::mat inv;
};

// The sampler's state over each study's `blocks`, its variants placed in the
// union by `union_idx`. A union variant's component is t: 0 for the null,
// else 1 + q (M - 1) + c for component c of kind q. Its effects in the
// studies `o` that both carry it and belong to t's kind are the ones t draws;
// it has none elsewhere.
template <class LdBlock>
class SdprChain {
public:
    SdprChain(const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks,
              const std::vector<arma::uvec> &union_idx, arma::uword n_union,
              const SdprSharing &sharing, size_t M, double a0k, double b0k,
              unsigned int seed, const ChainOptions &opts,
              const std::vector<double> &gap0);

    // Original: mcmc.cpp lines 434-513, one iteration of the main loop. The
    // kinds' shares and a learned rho are drawn only across studies.
    void sweep() {
        sample_sigma2();
        for (size_t j = 0; j < n_blocks_; j++) {
            calc_b(j);
        }
        for (size_t j = 0; j < n_blocks_; j++) {
            sample_assignment(j);
        }
        update_suffstats();
        sample_kind_weights();
        sample_V();
        update_p();
        sample_alpha();
        eta_beta = eta;
        for (size_t j = 0; j < n_blocks_; j++) {
            sample_beta(j);
        }
        if (learn_rho_) {
            sample_rho();
        }
        sample_eta();
    }

    // The reported effects are eta beta, and each study's heritability
    // eta^2 beta'R beta. Rao-Blackwellized, an effect is the conditional mean
    // of beta given the rest of the state when it was drawn, times eta at
    // that time. None is kept once any study's draws diverge.
    void keep(int n_kept) {
        bool diverged = false;
        for (size_t k = 0; k < K_; k++) {
            diverged = check_divergence(est_.study[k], opts_, eta * beta_[k],
                                        gap0_[k]) ||
                       diverged;
        }
        if (diverged) {
            return;
        }
        compute_h2();
        for (size_t k = 0; k < K_; k++) {
            est_.h2(k) += h2_[k] * square(eta) / n_kept;
            keep_effects(est_.study[k], opts_, eta * beta_[k],
                         eta_beta * beta_mean_[k], n_kept);
            if (opts_.track) {
                keep_fitted(k, n_kept);
            }
        }
        est_.rho += rho_ / n_kept;
    }

    void report(int itr, int /* n_iter */) {
        compute_h2();
        Rprintf("%d iter. h2:", itr);
        for (size_t k = 0; k < K_; k++) {
            Rprintf(" %g", h2_[k] * square(eta));
        }
        Rprintf(" max beta: %g\n", max_beta() * eta);
    }

    bool stopped() const { return est_.diverged(); }

    const SdprEstimates &estimates() const { return est_; }

private:
    const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks_;
    const std::vector<arma::uvec> &union_idx_;
    size_t K_;
    size_t n_blocks_;
    arma::uword n_union_;
    std::vector<BlockLayout> layouts_;
    size_t M;
    size_t n_kinds_;
    size_t T_;
    double a0k;
    double b0k;
    std::mt19937 r;

    // Each kind's studies; the distinct sets of studies carrying a variant.
    std::vector<StudySet> kind_set_;
    std::vector<StudySet> patterns_;
    arma::mat rho_;
    bool learn_rho_ = false;
    // Whether some kind holds both of a pair of studies.
    arma::umat paired_;
    // The correlations of each set of studies a variant's effects can span.
    std::map<StudySet, SetCorr> corr_;
    // Each union variant's studies, its place in each, and each block's union
    // variants in order.
    std::vector<StudySet> present_;
    arma::umat local_;
    std::vector<std::vector<arma::uword>> block_union_;

    double alpha;
    double eta;
    double N;
    std::vector<double> h2_;
    std::vector<arma::vec> beta_;
    std::vector<arma::vec> b_;
    std::vector<int> cls_assgn;
    std::vector<double> V;
    std::vector<double> p;
    std::vector<double> log_p;
    std::vector<double> cluster_var;
    std::vector<double> kind_weight_;
    std::vector<unsigned> suff_stats;
    // Each component's effects and their prior quadratic form,
    // beta' Sigma^-1 beta over each variant's effects (beta^2 for one).
    std::vector<unsigned> n_effects;
    std::vector<double> sumsq;
    // Each block's terms of the eta update, from sample_beta().
    std::vector<double> num;
    std::vector<double> denom;
    ChainOptions opts_;
    std::vector<double> gap0_;
    std::vector<arma::vec> beta_mean_;  // the conditional means, for Rao-Blackwell
    double eta_beta = 1.0;  // eta when beta was drawn
    SdprEstimates est_;

    size_t kind_of(size_t t) const { return (t - 1) / (M - 1); }

    // The set of studies in which component t draws variant u's effects.
    StudySet drawn_set(arma::uword u, size_t t) const {
        return kind_set_[kind_of(t)] & present_[u];
    }

    // Variant u's effects in `studies`.
    arma::vec effects(arma::uword u, const arma::uvec &studies) const {
        arma::vec x(studies.n_elem);
        for (arma::uword a = 0; a < studies.n_elem; a++) {
            x(a) = beta_[studies(a)](local_(u, studies(a)));
        }
        return x;
    }

    double max_beta() const {
        double top = -arma::datum::inf;
        for (const arma::vec &beta : beta_) {
            top = std::max(top, arma::max(beta));
        }
        return top;
    }

    void place_variants();
    void refresh_corr();
    void sample_sigma2();
    void calc_b(size_t j);
    void sample_assignment(size_t j);
    void spanning_probs(arma::uword u, size_t j, StudySet o, size_t lo, size_t hi,
                        const arma::fvec &cvar, const arma::fvec &log_p_fvec,
                        arma::fvec &prob) const;
    void update_suffstats();
    void sample_kind_weights();
    void sample_V();
    void update_p();
    void sample_alpha();
    void sample_beta(size_t j);
    void compute_h2();
    void sample_rho();
    double rho_log_target(const arma::mat &trial,
                          const std::vector<arma::uword> &spanning,
                          arma::uword k1, arma::uword k2) const;
    void sample_eta();
    void keep_fitted(size_t k, int n_kept);
};

template <class LdBlock>
SdprChain<LdBlock>::SdprChain(
    const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks,
    const std::vector<arma::uvec> &union_idx, arma::uword n_union,
    const SdprSharing &sharing, size_t M, double a0k, double b0k,
    unsigned int seed, const ChainOptions &opts, const std::vector<double> &gap0)
    : blocks_(blocks), union_idx_(union_idx), K_(blocks.size()),
      n_blocks_(blocks[0].size()), n_union_(n_union), M(M),
      n_kinds_(sharing.kinds.size()), T_(1 + sharing.kinds.size() * (M - 1)),
      a0k(a0k), b0k(b0k), r(seed), rho_(sharing.rho), opts_(opts), gap0_(gap0) {
    for (size_t k = 0; k < K_; k++) {
        layouts_.emplace_back(block_sizes(blocks[k]));
        arma::uword p_k = union_idx[k].n_elem;
        beta_.emplace_back(p_k, arma::fill::zeros);
        b_.emplace_back(p_k, arma::fill::zeros);
        beta_mean_.emplace_back(p_k, arma::fill::zeros);
        ChainEstimates est;
        est.beta.zeros(p_k);
        est.corr.zeros(p_k);
        est_.study.push_back(est);
    }
    place_variants();
    paired_.zeros(K_, K_);
    for (const arma::uvec &kind : sharing.kinds) {
        StudySet s = 0;
        for (arma::uword k : kind) {
            s |= StudySet(1) << k;
            for (arma::uword k2 : kind) {
                paired_(k, k2) = k != k2;
            }
        }
        kind_set_.push_back(s);
        learn_rho_ = learn_rho_ || (sharing.learn_rho && kind.n_elem >= 2);
    }
    refresh_corr();
    // N (sz) is absorbed into A and B (SDPR, May 20 2021), so it is 1.
    N = 1.0;
    alpha = 1;
    eta = 1;
    h2_.assign(K_, 0.0);
    p.assign(T_, 1.0 / T_);
    log_p.assign(T_, 0);
    for (size_t i = 0; i < T_; i++) {
        log_p[i] = std::log(static_cast<float>(p[i] + 1e-40));
    }

    cluster_var.assign(T_, 0.0);
    kind_weight_.assign(n_kinds_, 1.0 / n_kinds_);
    suff_stats.assign(T_, 0);
    n_effects.assign(T_, 0);
    sumsq.assign(T_, 0.0);
    V.assign(T_, 0.0);
    num.assign(n_blocks_, 0.0);
    denom.assign(n_blocks_, 0.0);
    // Every variant starts in the null component (k = 0). The original
    // started from a uniform draw over 0..M-1, which makes nearly every
    // variant causal, so the first sample_beta() factors a matrix the size
    // of the block; starting from null lets the sampler find the causal
    // assignments.
    cls_assgn.assign(n_union_, 0);
    est_.h2.zeros(K_);
    est_.rho.zeros(K_, K_);
    update_suffstats();
}

// Each union variant's studies and its place in each, and each block's union
// variants in order; a variant lies in the same block in every study
// carrying it.
template <class LdBlock>
void SdprChain<LdBlock>::place_variants() {
    present_.assign(n_union_, 0);
    local_.zeros(n_union_, K_);
    block_union_.assign(n_blocks_, std::vector<arma::uword>());
    std::vector<size_t> block_of(n_union_, n_blocks_);
    for (size_t k = 0; k < K_; k++) {
        for (size_t j = 0; j < n_blocks_; j++) {
            arma::uword start = layouts_[k].start(j);
            for (arma::uword i = start; i < start + layouts_[k].size(j); i++) {
                arma::uword u = union_idx_[k](i);
                present_[u] |= StudySet(1) << k;
                local_(u, k) = i;
                if (block_of[u] == n_blocks_) {
                    block_of[u] = j;
                    block_union_[j].push_back(u);
                } else if (block_of[u] != j) {
                    cpp11::stop("Union variant %d lies in different LD blocks "
                                "across the studies.", (int)u + 1);
                }
            }
        }
    }
    for (arma::uword u = 0; u < n_union_; u++) {
        if (present_[u] == 0) {
            cpp11::stop("Union variant %d is in no study.", (int)u + 1);
        }
    }
    for (std::vector<arma::uword> &members : block_union_) {
        std::sort(members.begin(), members.end());
    }
    patterns_ = present_;
    std::sort(patterns_.begin(), patterns_.end());
    patterns_.erase(std::unique(patterns_.begin(), patterns_.end()),
                    patterns_.end());
}

// The correlations at the current rho of every set of two or more studies
// a variant's effects can span: a kind's studies that carry it.
template <class LdBlock>
void SdprChain<LdBlock>::refresh_corr() {
    corr_.clear();
    for (StudySet kind : kind_set_) {
        for (StudySet pattern : patterns_) {
            StudySet o = kind & pattern;
            if (set_size(o) < 2 || corr_.count(o) > 0) {
                continue;
            }
            SetCorr c;
            c.studies = set_studies(o);
            arma::mat S = rho_.submat(c.studies, c.studies);
            c.root = arma::sqrtmat_sympd(S);
            c.inv = arma::inv_sympd(S);
            corr_.emplace(o, c);
        }
    }
}

// ---------------------------------------------------------------------------
// sample_sigma2: Sample cluster variances from inverse-gamma posterior
//   var_k ~ InvGamma(n_effects[k]/2 + a0k, sumsq[k]/2 + b0k)
// Original: mcmc.cpp lines 25-40, where a component's effects are its
// variants (n_effects = suff_stats); across studies, each variant counts
// the effects it has in the component's studies.
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_sigma2() {
    std::gamma_distribution<double> dist;
    for (size_t i = 1; i < T_; i++) {
        double a = n_effects[i] / 2.0 + a0k;
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
// calc_b: Compute the b vector for block j of each study
//   b_j = eta^2 * (diag(B)*beta - B*beta) + eta * A^T * beta_mrg
// Original: mcmc.cpp lines 42-62
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::calc_b(size_t j) {
    for (size_t k = 0; k < K_; k++) {
        if (layouts_[k].size(j) == 0) {
            continue;
        }
        const SdprBlock<LdBlock> &blk = blocks_[k][j];
        arma::span span = layouts_[k].span(j);
        arma::vec beta_j = beta_[k](span);
        b_[k](span) = eta * eta * (blk.B.diag() % beta_j - blk.B.times(beta_j))
                      + eta * blk.calc_b_tmp;
    }
}

// ---------------------------------------------------------------------------
// sample_assignment: Sample cluster assignments for each union variant in
// block j
//   For each variant i, compute log P(z_i = k) for k = 0..T-1, then sample.
//
// Original: mcmc.cpp lines 64-165, which vectorized log/exp over the M
// clusters with SSE intrinsics in single precision; Armadillo's float
// vectors do the same.
//
// Math for cluster k >= 1 whose effects on i lie in one study:
//   C_k = eta^2 * N * B[i,i] * var_k + 1
//   log P(z_i=k) = -0.5*log(C_k) + log(p_k) + (N*b_i)^2 * var_k / (2*C_k)
// In none, log P(z_i=k) = log(p_k); in several, spanning_probs().
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_assignment(size_t j) {
    std::uniform_real_distribution<float> unif(0.0f, 1.0f);

    // N = 1.0 after May 21 2021 (absorbed into A, B in solve_ldmat)
    float C = static_cast<float>(eta * eta * N);

    arma::fvec cvar(T_);
    arma::fvec log_p_fvec(T_);
    for (size_t k = 0; k < T_; k++) {
        cvar(k)       = static_cast<float>(cluster_var[k]);
        log_p_fvec(k) = static_cast<float>(log_p[k]);
    }

    for (arma::uword u : block_union_[j]) {
        float rnd_i = unif(r);

        arma::fvec prob(T_);
        prob(0) = log_p_fvec(0);

        for (size_t q = 0; q < n_kinds_; q++) {
            size_t lo = 1 + q * (M - 1);
            size_t hi = lo + M - 2;
            StudySet o = kind_set_[q] & present_[u];
            int n_o = set_size(o);
            if (n_o == 0) {
                prob.subvec(lo, hi) = log_p_fvec.subvec(lo, hi);
            } else if (n_o == 1) {
                arma::uword k = set_first(o);
                arma::uword i = local_(u, k);
                float Bjj_i = static_cast<float>(
                    blocks_[k][j].B.diag()(i - layouts_[k].start(j)));
                float bj_i  = static_cast<float>(b_[k](i));
                arma::fvec Ck = C * Bjj_i * cvar.subvec(lo, hi) + 1.0f;
                prob.subvec(lo, hi) = -0.5f * arma::log(Ck) + log_p_fvec.subvec(lo, hi)
                    + square(N * bj_i) * cvar.subvec(lo, hi) / (2.0f * Ck);
            } else {
                spanning_probs(u, j, o, lo, hi, cvar, log_p_fvec, prob);
            }
        }

        // Log-sum-exp for numerical stability
        float max_elem    = prob.max();
        float log_exp_sum = max_elem
            + std::log(arma::accu(arma::exp(prob - max_elem)));

        // Categorical sampling via inverse CDF
        // Original: mcmc.cpp lines 155-163
        cls_assgn[u] = T_ - 1;
        for (size_t k = 0; k < T_ - 1; k++) {
            rnd_i -= std::exp(prob(k) - log_exp_sum);
            if (rnd_i < 0) {
                cls_assgn[u] = k;
                break;
            }
        }
    }
}

// The log weights of components lo..hi, of a kind whose effects on union
// variant u lie in the studies `o`, two or more, drawn N(0, var_k Sigma):
// with D = diag(eta^2 N B_ii) and h = N b_i over those studies and
// Sigma^1/2 D Sigma^1/2 = Q diag(lambda) Q', g = Q' Sigma^1/2 h,
//   log P(z_i=k) = log(p_k) - 0.5 * sum_a log(1 + var_k lambda_a)
//                  + 0.5 * var_k * sum_a g_a^2 / (1 + var_k lambda_a),
// SDPRX's closed forms for its shared component, in any number of studies.
template <class LdBlock>
void SdprChain<LdBlock>::spanning_probs(arma::uword u, size_t j, StudySet o,
                                        size_t lo, size_t hi,
                                        const arma::fvec &cvar,
                                        const arma::fvec &log_p_fvec,
                                        arma::fvec &prob) const {
    const SetCorr &c = corr_.at(o);
    arma::uword n_o = c.studies.n_elem;
    arma::vec d(n_o);
    arma::vec h(n_o);
    for (arma::uword a = 0; a < n_o; a++) {
        arma::uword k = c.studies(a);
        arma::uword i = local_(u, k);
        d(a) = eta * eta * N * blocks_[k][j].B.diag()(i - layouts_[k].start(j));
        h(a) = N * b_[k](i);
    }
    arma::vec lambda;
    arma::mat Q;
    arma::eig_sym(lambda, Q,
                  arma::mat(arma::symmatu(c.root * arma::diagmat(d) * c.root)));
    arma::vec g = Q.t() * (c.root * h);
    arma::fvec cv = cvar.subvec(lo, hi);
    arma::fvec logdet(cv.n_elem, arma::fill::zeros);
    arma::fvec quad(cv.n_elem, arma::fill::zeros);
    for (arma::uword a = 0; a < n_o; a++) {
        arma::fvec t = 1.0f + static_cast<float>(std::max(lambda(a), 0.0)) * cv;
        logdet += arma::log(t);
        quad += static_cast<float>(square(g(a))) * cv / t;
    }
    prob.subvec(lo, hi) = log_p_fvec.subvec(lo, hi) - 0.5f * logdet + 0.5f * quad;
}

// ---------------------------------------------------------------------------
// update_suffstats: Count variants per cluster, and the effects each
// component draws and their quadratic form
// Original: mcmc.cpp lines 167-175
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::update_suffstats() {
    std::fill(suff_stats.begin(), suff_stats.end(), 0);
    std::fill(n_effects.begin(), n_effects.end(), 0);
    std::fill(sumsq.begin(), sumsq.end(), 0.0);
    for (arma::uword u = 0; u < n_union_; u++) {
        size_t t = cls_assgn[u];
        suff_stats[t]++;
        if (t == 0) {
            continue;
        }
        StudySet o = drawn_set(u, t);
        int n_o = set_size(o);
        if (n_o == 1) {
            arma::uword k = set_first(o);
            double tmp = beta_[k](local_(u, k));
            sumsq[t] += square(tmp);
        } else if (n_o >= 2) {
            const SetCorr &c = corr_.at(o);
            arma::vec x = effects(u, c.studies);
            sumsq[t] += arma::dot(x, c.inv * x);
        }
        n_effects[t] += static_cast<unsigned>(n_o);
    }
}

// ---------------------------------------------------------------------------
// sample_kind_weights: Sample the non-null kinds' shares of the weight
//   from the Dirichlet given the kinds' counts (SDPRX's sample_pi_pop(),
//   without the null, which is the first stick here). One kind has it all.
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_kind_weights() {
    if (n_kinds_ == 1) {
        return;
    }
    std::vector<double> counts(n_kinds_, 0.0);
    for (size_t t = 1; t < T_; t++) {
        counts[kind_of(t)] += suff_stats[t];
    }
    kind_weight_ = dirichlet_plus_one(counts, r);
}

// ---------------------------------------------------------------------------
// sample_V: Sample stick-breaking fractions from Beta posterior
//   V[k] ~ Beta(1 + suff_stats[k], alpha + sum_{k'>k} suff_stats[k'])
//   over the null's stick against every other component, then over each
//   kind's components, its last taking the remainder.
// Original: mcmc.cpp lines 177-189
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_V() {
    double rest = 0;
    for (size_t t = 1; t < T_; t++) {
        rest += suff_stats[t];
    }
    beta_distribution first(1 + suff_stats[0], alpha + rest);
    V[0] = first(r);

    for (size_t q = 0; q < n_kinds_; q++) {
        size_t lo = 1 + q * (M - 1);
        size_t hi = lo + M - 2;
        std::vector<double> a(M - 2);
        a[M - 3] = suff_stats[hi];
        for (int i = static_cast<int>(M) - 4; i >= 0; i--) {
            a[i] = suff_stats[lo + i + 1] + a[i + 1];
        }
        for (size_t i = 0; i < M - 2; i++) {
            beta_distribution dist(1 + suff_stats[lo + i], alpha + a[i]);
            V[lo + i] = dist(r);
        }
        V[hi] = 1;
    }
}

// ---------------------------------------------------------------------------
// update_p: Convert stick-breaking fractions V to cluster probabilities p
//   p[0] = V[0]; each kind's share of (1 - V[0]) broken over its components,
//   p[k] = share * V[k] * prod_{j<k}(1-V[j]) within the kind. The last
//   component takes what the others leave.
// Original: mcmc.cpp lines 191-222
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::update_p() {
    p[0] = V[0];
    double rest = 1 - V[0];
    for (size_t q = 0; q < n_kinds_; q++) {
        size_t lo = 1 + q * (M - 1);
        size_t hi = lo + M - 2;
        double cumprod = rest * kind_weight_[q];
        for (size_t t = lo; t < hi; t++) {
            p[t] = cumprod * V[t];
            cumprod *= (1 - V[t]);
        }
        p[hi] = cumprod;
    }

    double sum = std::accumulate(p.begin(), p.end() - 1, 0.0);
    p[T_ - 1] = (1 - sum > 0) ? (1 - sum) : 0;

    for (size_t i = 0; i < T_; i++) {
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
    for (size_t i = 0; i < T_; i++) {
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
//   Across studies the causal effects are each study's, stacked: B_gamma is
//   block diagonal over the studies, and a variant's effects in several
//   studies have the prior precision Sigma^-1 / var_k between them.
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
    // Each study's causal SNPs in the block: those whose cluster draws an
    // effect in that study. Zero out the block's betas.
    std::vector<arma::uvec> causal(K_);
    arma::uword nc = 0;
    for (size_t k = 0; k < K_; k++) {
        if (layouts_[k].size(j) == 0) {
            continue;
        }
        arma::span span = layouts_[k].span(j);
        beta_[k](span).zeros();
        beta_mean_[k](span).zeros();
        arma::uword start = layouts_[k].start(j);
        std::vector<arma::uword> causal_list;
        for (arma::uword i = 0; i < layouts_[k].size(j); i++) {
            int t = cls_assgn[union_idx_[k](start + i)];
            if (t != 0 && set_holds(kind_set_[kind_of(t)], k)) {
                causal_list.push_back(i);
            }
        }
        causal[k] = arma::conv_to<arma::uvec>::from(causal_list);
        nc += causal[k].n_elem;
    }

    if (nc == 0) {
        num[j]   = 0;
        denom[j] = 0;
        return;
    }

    // Original: mcmc.cpp lines 273-289. B_eta is eta^2 * N * B over the
    // causal SNPs, with B's diagonal as the original reads it for eta's
    // terms (lines 312-321).
    double C_coeff = square(eta) * N;
    arma::vec A_vec(nc);
    arma::mat B_gamma(nc, nc, arma::fill::zeros);
    arma::mat B_eta(nc, nc, arma::fill::zeros);
    std::map<arma::uword, std::vector<arma::uword>> spanning;
    arma::uword at = 0;
    for (size_t k = 0; k < K_; k++) {
        if (causal[k].is_empty()) {
            continue;
        }
        const SdprBlock<LdBlock> &blk = blocks_[k][j];
        arma::span stacked(at, at + causal[k].n_elem - 1);
        A_vec(stacked) = N * eta * blk.calc_b_tmp(causal[k]);
        arma::mat B_k = C_coeff * blk.B.sub(causal[k]);
        B_gamma(stacked, stacked) = B_k;
        B_k.diag() = C_coeff * blk.B.diag()(causal[k]);
        B_eta(stacked, stacked) = B_k;
        arma::uword start = layouts_[k].start(j);
        for (arma::uword i : causal[k]) {
            arma::uword u = union_idx_[k](start + i);
            int t = cls_assgn[u];
            if (set_size(drawn_set(u, t)) == 1) {
                B_gamma(at, at) += 1.0 / cluster_var[t];
            } else {
                spanning[u].push_back(at);
            }
            at++;
        }
    }
    for (const auto &entry : spanning) {
        int t = cls_assgn[entry.first];
        const SetCorr &c = corr_.at(drawn_set(entry.first, t));
        arma::uvec places = arma::conv_to<arma::uvec>::from(entry.second);
        B_gamma.submat(places, places) += c.inv / cluster_var[t];
    }

    std::normal_distribution<double> dist(0.0, 1.0);
    arma::vec z = standard_normals(nc, dist, r);
    arma::vec mean_c;
    arma::vec beta_c = opts_.rao_blackwell ? gaussian_draw(B_gamma, A_vec, z, 1.0, mean_c)
                                           : gaussian_draw(B_gamma, A_vec, z);

    num[j] = arma::dot(A_vec, beta_c) / eta;
    arma::vec tmp = B_eta * beta_c;
    denom[j] = arma::dot(beta_c, tmp) / square(eta);

    at = 0;
    for (size_t k = 0; k < K_; k++) {
        arma::uword start = layouts_[k].start(j);
        for (arma::uword i : causal[k]) {
            beta_[k](start + i) = beta_c(at);
            if (opts_.rao_blackwell) {
                beta_mean_[k](start + i) = mean_c(at);
            }
            at++;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_h2: Compute each study's heritability h2 = sum_j beta_j^T R_j beta_j
// Original: mcmc.cpp lines 331-343
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::compute_h2() {
    for (size_t k = 0; k < K_; k++) {
        h2_[k] = 0;
        for (size_t j = 0; j < n_blocks_; j++) {
            if (layouts_[k].size(j) == 0) {
                continue;
            }
            arma::vec beta_j = beta_[k](layouts_[k].span(j));
            h2_[k] += blocks_[k][j].R.quad(beta_j);
        }
    }
}

// ---------------------------------------------------------------------------
// sample_rho: Draw each pair of studies' correlation given the rest
//   On the grid 0, 0.01, ..., 0.99, the conditional is the density of the
//   effects of every variant whose cluster spans both studies,
//   N(0, var_k Sigma) over the studies it spans, at the grid values that
//   keep rho positive definite. A random walk from 0, as JointPRS takes, moves
//   little in a chain's length when few effects inform rho.
// ---------------------------------------------------------------------------
template <class LdBlock>
void SdprChain<LdBlock>::sample_rho() {
    std::vector<arma::uword> spanning;
    for (arma::uword u = 0; u < n_union_; u++) {
        int t = cls_assgn[u];
        if (t != 0 && set_size(drawn_set(u, t)) >= 2) {
            spanning.push_back(u);
        }
    }
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (arma::uword k1 = 0; k1 < K_; k1++) {
        for (arma::uword k2 = k1 + 1; k2 < K_; k2++) {
            if (!paired_(k1, k2)) {
                continue;
            }
            arma::vec log_w(rho_grid_points);
            for (int g = 0; g < rho_grid_points; g++) {
                arma::mat trial = rho_;
                trial(k1, k2) = g * rho_grid_step;
                trial(k2, k1) = g * rho_grid_step;
                log_w(g) = rho_log_target(trial, spanning, k1, k2);
            }
            double top = log_w.max();
            if (!std::isfinite(top)) {
                continue;
            }
            arma::vec w = arma::exp(log_w - top);
            double rnd = unif(r) * arma::accu(w);
            int pick = rho_grid_points - 1;
            for (int g = 0; g < rho_grid_points - 1; g++) {
                rnd -= w(g);
                if (rnd < 0) {
                    pick = g;
                    break;
                }
            }
            rho_(k1, k2) = pick * rho_grid_step;
            rho_(k2, k1) = pick * rho_grid_step;
        }
    }
    refresh_corr();
}

// The log density, up to a constant, of the effects of the variants in
// `spanning` whose clusters span studies k1 and k2, at the correlations
// `trial`; -inf where `trial` is not positive definite.
template <class LdBlock>
double SdprChain<LdBlock>::rho_log_target(const arma::mat &trial,
                                          const std::vector<arma::uword> &spanning,
                                          arma::uword k1, arma::uword k2) const {
    if (K_ > 2 && arma::eig_sym(trial).min() <= 0.0) {
        return -arma::datum::inf;
    }
    double total = 0.0;
    for (arma::uword u : spanning) {
        int t = cls_assgn[u];
        StudySet o = drawn_set(u, t);
        if (!set_holds(o, k1) || !set_holds(o, k2)) {
            continue;
        }
        arma::uvec studies = set_studies(o);
        arma::mat L;
        if (!arma::chol(L, arma::mat(trial.submat(studies, studies)), "lower")) {
            return -arma::datum::inf;
        }
        arma::vec w = arma::solve(arma::trimatl(L), effects(u, studies));
        total -= arma::accu(arma::log(L.diag())) + 0.5 * arma::dot(w, w) / cluster_var[t];
    }
    return total;
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

// Study k's fitted correlations, R eta beta, added block by block.
template <class LdBlock>
void SdprChain<LdBlock>::keep_fitted(size_t k, int n_kept) {
    for (size_t j = 0; j < n_blocks_; j++) {
        if (layouts_[k].size(j) == 0) {
            continue;
        }
        arma::span blk = layouts_[k].span(j);
        est_.study[k].corr(blk) += blocks_[k][j].R.times(eta * beta_[k](blk)) / n_kept;
    }
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

        // A block holding no variant has nothing to solve; the chain skips
        // it, as PRS-CS's does.
        if (blk_size == 0) {
            blocks.push_back(SdprBlock<DenseLdBlock>{
                DenseLdBlock(arma::mat()), DenseLdBlock(arma::mat()), arma::vec()
            });
            continue;
        }

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

namespace {

// The inputs sdpr_mcmc() reads, checked before the chain is built.
template <class LdBlock>
void check_sdpr_inputs(const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks,
                       const std::vector<arma::uvec> &union_idx,
                       const SdprSharing &sharing, size_t M,
                       const std::vector<double> &gap0) {
    size_t K = blocks.size();
    if (K == 0 || K > 64) {
        cpp11::stop("SDPR runs on 1 to 64 studies, not %d.", (int)K);
    }
    if (union_idx.size() != K || gap0.size() != K) {
        cpp11::stop("Each study needs its places in the union and its bound "
                    "of divergence.");
    }
    if (M < 4) {
        cpp11::stop("'M' must be at least 4.");
    }
    for (size_t k = 0; k < K; k++) {
        if (blocks[k].size() != blocks[0].size()) {
            cpp11::stop("Every study needs the same number of LD blocks.");
        }
        arma::uword p = 0;
        for (const SdprBlock<LdBlock> &blk : blocks[k]) {
            p += blk.size();
        }
        if (p != union_idx[k].n_elem) {
            cpp11::stop("Study %d has %d variants in its LD blocks but %d "
                        "places in the union.", (int)k + 1, (int)p,
                        (int)union_idx[k].n_elem);
        }
    }
    if (sharing.kinds.empty()) {
        cpp11::stop("SDPR needs at least one kind of component.");
    }
    for (const arma::uvec &kind : sharing.kinds) {
        if (kind.is_empty() || kind.max() >= K) {
            cpp11::stop("Each kind must hold some of the %d studies.", (int)K);
        }
    }
    if (sharing.rho.n_rows != K || sharing.rho.n_cols != K) {
        cpp11::stop("rho must be %d x %d.", (int)K, (int)K);
    }
}

}  // namespace

template <class LdBlock>
SdprEstimates sdpr_mcmc(const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks,
                        const std::vector<arma::uvec> &union_idx, arma::uword n_union,
                        const SdprSharing &sharing, size_t M, double a0k, double b0k,
                        const ChainSchedule &schedule, bool verbose, unsigned int seed,
                        const ChainOptions &opts, const std::vector<double> &gap0) {
    check_sdpr_inputs(blocks, union_idx, sharing, M, gap0);
    SdprChain<LdBlock> chain(blocks, union_idx, n_union, sharing, M, a0k, b0k, seed,
                             opts, gap0);
    run_chain(chain, schedule, verbose);
    const SdprEstimates &est = chain.estimates();
    if (verbose) {
        Rprintf("h2:");
        double top = -arma::datum::inf;
        for (size_t k = 0; k < est.study.size(); k++) {
            Rprintf(" %g", est.h2(k));
            top = std::max(top, arma::max(est.study[k].beta));
        }
        Rprintf(" max: %g\n", top);
    }
    return est;
}

template SdprEstimates sdpr_mcmc<DenseLdBlock>(
    const std::vector<std::vector<SdprBlock<DenseLdBlock>>> &,
    const std::vector<arma::uvec> &, arma::uword, const SdprSharing &, size_t, double,
    double, const ChainSchedule &, bool, unsigned int, const ChainOptions &,
    const std::vector<double> &);
template SdprEstimates sdpr_mcmc<FactorLdBlock>(
    const std::vector<std::vector<SdprBlock<FactorLdBlock>>> &,
    const std::vector<arma::uvec> &, arma::uword, const SdprSharing &, size_t, double,
    double, const ChainSchedule &, bool, unsigned int, const ChainOptions &,
    const std::vector<double> &);

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

// ==== PRS-CS, PRS-CSx and JointPRS ============================================

namespace {

// An LD matrix as the PRS-CS chain reads it: c R + (1 - c) I for the LD
// shrinkage c = `shrinkCorr`.
DenseLdBlock shrunk_dense_block(mat S, double shrinkCorr) {
	if (shrinkCorr != 1.0) {
		S *= shrinkCorr;
		S.diag() += 1.0 - shrinkCorr;
	}
	return DenseLdBlock(std::move(S));
}

// A factor F of R = F'F / (m - 1) as the chain reads it: the factor of c R,
// and the ridge 1 - c.
FactorLdBlock shrunk_factor_block(const mat& F, double shrinkCorr) {
	mat W = F / std::sqrt(static_cast<double>(F.n_rows) - 1.0);
	if (shrinkCorr != 1.0) {
		W *= std::sqrt(shrinkCorr);
	}
	return FactorLdBlock(std::move(W), 1.0 - shrinkCorr);
}

// One study's places in a union of its own variants: 0, ..., p - 1.
uvec own_places(uword p) {
	uvec idx(p);
	for (uword i = 0; i < p; ++i) {
		idx(i) = i;
	}
	return idx;
}

// One study's chain with the arguments as R passes them: a NULL `maf` keeps
// the effects standardized, else they are returned per allele; a NULL `phi`
// is learned from the data, and a NULL `seed` draws one.
template <class LdBlock>
cpp11::writable::list prs_cs_run(double a, double b, sexp phi, doubles bhat, sexp maf, int n,
                                 const std::vector<std::vector<LdBlock>>& blocks,
                                 int nIter, int nBurnin, int thin, bool verbose,
                                 sexp seed, const ChainOptions& opts) {
	vec beta_mrg = as_Col(bhat);
	vec maf_vec = maf == R_NilValue ? vec(beta_mrg.n_elem, fill::zeros)
	                                : as_Col(doubles(maf));
	bool learn_phi = phi == R_NilValue;
	double phi_val = learn_phi ? 0.0 : cpp11::as_cpp<double>(phi);
	ChainSchedule schedule(nIter, nBurnin, thin);
	vec one_constraint = {1.0};
	vec one_n = {static_cast<double>(n)};
	PrsCsEstimates est = prs_cs_mcmc(a, b, learn_phi, phi_val, one_constraint, {beta_mrg},
	                                 one_n, blocks, {own_places(beta_mrg.n_elem)},
	                                 beta_mrg.n_elem, schedule, verbose, seed_from(seed),
	                                 opts);
	ChainEstimates one = est.study[0];
	if (arma::max(maf_vec) > 0) {
		one.beta /= arma::sqrt(2.0 * maf_vec % (1.0 - maf_vec));
	}
	using namespace cpp11::literals;
	return chain_report(one, "psiEst"_nm = as_doubles(est.psi),
	                    "sigmaEst"_nm = cpp11::as_sexp(est.sigma(0)),
	                    "phiEst"_nm = cpp11::as_sexp(est.phi));
}

// The studies' marginal effects.
std::vector<vec> study_effects(const list& bhatList) {
	std::vector<vec> out;
	for (R_xlen_t k = 0; k < bhatList.size(); ++k) {
		out.push_back(as_Col(doubles(bhatList[k])));
	}
	return out;
}

// Several studies' chain with the arguments as R passes them, each study's
// effects and fitted correlations reported over the union, zero where it
// lacks a variant.
template <class LdBlock>
cpp11::writable::list joint_prs_run(double a, double b, sexp phi,
                                    const doubles& constraint,
                                    const std::vector<vec>& beta_mrg,
                                    const doubles& n,
                                    const std::vector<std::vector<LdBlock>>& blocks,
                                    const std::vector<uvec>& union_idx, int nUnion,
                                    int nIter, int nBurnin, int thin, bool verbose,
                                    sexp seed, const ChainOptions& opts) {
	if (constraint.size() != static_cast<R_xlen_t>(beta_mrg.size()) ||
	    n.size() != static_cast<R_xlen_t>(beta_mrg.size())) {
		cpp11::stop("Each study needs one sample size and one constraint.");
	}
	bool learn_phi = phi == R_NilValue;
	double phi_val = learn_phi ? 0.0 : cpp11::as_cpp<double>(phi);
	ChainSchedule schedule(nIter, nBurnin, thin);
	PrsCsEstimates est = prs_cs_mcmc(
	    a, b, learn_phi, phi_val, as_Col(constraint), beta_mrg, as_Col(n), blocks,
	    union_idx, static_cast<uword>(nUnion), schedule, verbose, seed_from(seed), opts);
	mat beta(nUnion, beta_mrg.size(), fill::zeros);
	mat corr(nUnion, beta_mrg.size(), fill::zeros);
	for (uword k = 0; k < beta_mrg.size(); ++k) {
		uvec col = {k};
		beta.submat(union_idx[k], col) = est.study[k].beta;
		corr.submat(union_idx[k], col) = est.study[k].corr;
	}
	using namespace cpp11::literals;
	return writable::list({
		"betaEst"_nm = as_doubles_matrix(beta),
		"psiEst"_nm = as_doubles(est.psi),
		"sigmaEst"_nm = as_doubles(est.sigma),
		"phiEst"_nm = cpp11::as_sexp(est.phi),
		"rhoEst"_nm = as_doubles_matrix(est.rho),
		"corrEst"_nm = as_doubles_matrix(corr),
		"diverged"_nm = cpp11::as_sexp(est.diverged())
	});
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
	std::vector<std::vector<DenseLdBlock>> blocks(1);
	for (int i = 0; i < ldBlk.size(); ++i) {
		blocks[0].push_back(shrunk_dense_block(as_Mat(doubles_matrix<>(ldBlk[i])),
		                                       shrinkCorr));
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
	std::vector<std::vector<FactorLdBlock>> blocks(1);
	blocks[0].push_back(shrunk_factor_block(checked_factor(factorR, bhat.size()),
	                                        shrinkCorr));
	return prs_cs_run(a, b, phi, bhat, maf, n, blocks, nIter, nBurnin, thin, verbose, seed,
	                  chain_options(raoBlackwell, track));
}

/**
 * @brief PRS-CSx / JointPRS from each study's LD matrix.
 *
 * K studies' PRS-CS chain sharing the local and global scales, their effects
 * correlated by rho (PrsCsChain). `bhatList` holds each study's marginal
 * effects, `n` its sample size, `ldList` its LD matrix (one block, shrunk to
 * c R + (1 - c) I by `shrinkCorr` as prsCsCpp() shrinks it) and
 * `unionIdxList` its variants' 1-based places in the union of the studies'
 * variants, of size `nUnion`, in one allele frame. `constraint` bounds each
 * study's correlations: 1 learns them on [0, 0.99], 0 holds them at 0, which is
 * PRS-CSx. The rest are prsCsCpp()'s. Returns the posterior means: each
 * study's effects over the union (`betaEst`, union x studies, zero where a
 * study lacks a variant), the shared `psiEst` and `phiEst`, each study's
 * `sigmaEst`, the correlations `rhoEst`, and the fitted correlations and
 * divergence every chain tracks.
 */
[[cpp11::register]]
cpp11::writable::list jointPrsCpp(double a, double b, sexp phi, doubles constraint,
                                  list bhatList, doubles n, list ldList,
                                  list unionIdxList, int nUnion, int nIter,
                                  int nBurnin, int thin, bool verbose, sexp seed,
                                  double shrinkCorr, bool raoBlackwell, bool track) {
	std::vector<vec> beta_mrg = study_effects(bhatList);
	std::vector<uvec> union_idx = union_index(unionIdxList, beta_mrg, nUnion, "effects");
	if (ldList.size() != static_cast<R_xlen_t>(beta_mrg.size())) {
		cpp11::stop("There are %d LD matrices for %d studies.", (int)ldList.size(),
		            (int)beta_mrg.size());
	}
	std::vector<std::vector<DenseLdBlock>> blocks(beta_mrg.size());
	for (uword k = 0; k < beta_mrg.size(); ++k) {
		mat S = as_Mat(doubles_matrix<>(ldList[static_cast<R_xlen_t>(k)]));
		if (S.n_rows != beta_mrg[k].n_elem || S.n_cols != beta_mrg[k].n_elem) {
			cpp11::stop("Study %d's LD is %d x %d but it has %d effects.", (int)k + 1,
			            (int)S.n_rows, (int)S.n_cols, (int)beta_mrg[k].n_elem);
		}
		blocks[k].push_back(shrunk_dense_block(std::move(S), shrinkCorr));
	}
	return joint_prs_run(a, b, phi, constraint, beta_mrg, n, blocks, union_idx, nUnion,
	                     nIter, nBurnin, thin, verbose, seed,
	                     chain_options(raoBlackwell, track));
}

/**
 * @brief PRS-CSx / JointPRS from a factor F of each study's LD, R = F'F / (m - 1),
 * never forming R.
 *
 * The arguments are jointPrsCpp()'s with `factorList` in place of the LD
 * matrices; each study's chain draws as prsCsFactorCpp()'s does.
 */
[[cpp11::register]]
cpp11::writable::list jointPrsFactorCpp(double a, double b, sexp phi,
                                        doubles constraint, list bhatList, doubles n,
                                        list factorList, list unionIdxList, int nUnion,
                                        int nIter, int nBurnin, int thin, bool verbose,
                                        sexp seed, double shrinkCorr, bool raoBlackwell,
                                        bool track) {
	std::vector<vec> beta_mrg = study_effects(bhatList);
	std::vector<uvec> union_idx = union_index(unionIdxList, beta_mrg, nUnion, "effects");
	if (factorList.size() != static_cast<R_xlen_t>(beta_mrg.size())) {
		cpp11::stop("There are %d LD factors for %d studies.", (int)factorList.size(),
		            (int)beta_mrg.size());
	}
	std::vector<std::vector<FactorLdBlock>> blocks(beta_mrg.size());
	for (uword k = 0; k < beta_mrg.size(); ++k) {
		blocks[k].push_back(shrunk_factor_block(
		    checked_factor(doubles_matrix<>(factorList[static_cast<R_xlen_t>(k)]),
		                   static_cast<R_xlen_t>(beta_mrg[k].n_elem)),
		    shrinkCorr));
	}
	return joint_prs_run(a, b, phi, constraint, beta_mrg, n, blocks, union_idx, nUnion,
	                     nIter, nBurnin, thin, verbose, seed,
	                     chain_options(raoBlackwell, track));
}

// ==== SDPR and SDPRX ============================================================

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

// One study's chain over its blocks: SDPR, with the one kind {0}.
template <class LdBlock>
cpp11::writable::list sdpr_one(const std::vector<std::vector<SdprBlock<LdBlock>>>& blocks,
                               const doubles& bhatR, int M, double a0k, double b0k,
                               const ChainSchedule& schedule, bool verbose, sexp seed,
                               const ChainOptions& opts) {
	SdprSharing sharing;
	sharing.kinds = {uvec{0}};
	sharing.rho = mat(1, 1, fill::eye);
	SdprEstimates est = sdpr_mcmc(blocks, {own_places(bhatR.size())}, bhatR.size(),
	                              sharing, M, a0k, b0k, schedule, verbose,
	                              seed_from(seed), opts, {divergence_gap(as_Col(bhatR))});
	using namespace cpp11::literals;
	return chain_report(est.study[0], "h2"_nm = cpp11::as_sexp(est.h2(0)));
}

// The kinds and rho as R passes them: each kind's 1-based studies, and rho a
// K x K correlation matrix, or NULL to learn it starting from the identity.
SdprSharing sdprx_sharing(const list& kinds, sexp rho, uword K) {
	SdprSharing sharing;
	std::vector<std::vector<uword>> seen;
	for (R_xlen_t q = 0; q < kinds.size(); ++q) {
		integers studies(kinds[q]);
		std::vector<uword> kind;
		for (int v : studies) {
			if (v < 1 || v > static_cast<int>(K)) {
				cpp11::stop("Kind %d's studies must lie in 1..%d.", (int)q + 1, (int)K);
			}
			kind.push_back(static_cast<uword>(v - 1));
		}
		std::sort(kind.begin(), kind.end());
		if (kind.empty() || std::adjacent_find(kind.begin(), kind.end()) != kind.end() ||
		    std::find(seen.begin(), seen.end(), kind) != seen.end()) {
			cpp11::stop("Kind %d must name distinct studies, and differ from the other "
			            "kinds.", (int)q + 1);
		}
		seen.push_back(kind);
		sharing.kinds.push_back(conv_to<uvec>::from(kind));
	}
	if (rho == R_NilValue) {
		sharing.rho = mat(K, K, fill::eye);
		sharing.learn_rho = true;
		return sharing;
	}
	sharing.rho = as_Mat(doubles_matrix<>(rho));
	mat L;
	if (sharing.rho.n_rows != K || sharing.rho.n_cols != K ||
	    !sharing.rho.is_symmetric() || arma::any(sharing.rho.diag() != 1.0) ||
	    !arma::chol(L, sharing.rho)) {
		cpp11::stop("rho must be a %d x %d positive-definite correlation matrix.",
		            (int)K, (int)K);
	}
	return sharing;
}

// Several studies' chain with the arguments as R passes them, each study's
// effects and fitted correlations reported over the union, zero where it
// lacks a variant.
template <class LdBlock>
cpp11::writable::list sdprx_run(const std::vector<std::vector<SdprBlock<LdBlock>>>& blocks,
                                const std::vector<vec>& beta_mrg,
                                const std::vector<uvec>& union_idx, int nUnion,
                                const SdprSharing& sharing, int M, double a0k, double b0k,
                                const ChainSchedule& schedule, bool verbose, sexp seed,
                                const ChainOptions& opts) {
	std::vector<double> gap0;
	for (const vec& b : beta_mrg) {
		gap0.push_back(divergence_gap(b));
	}
	SdprEstimates est = sdpr_mcmc(blocks, union_idx, static_cast<uword>(nUnion), sharing,
	                              M, a0k, b0k, schedule, verbose, seed_from(seed), opts,
	                              gap0);
	mat beta(nUnion, beta_mrg.size(), fill::zeros);
	mat corr(nUnion, beta_mrg.size(), fill::zeros);
	for (uword k = 0; k < beta_mrg.size(); ++k) {
		uvec col = {k};
		beta.submat(union_idx[k], col) = est.study[k].beta;
		corr.submat(union_idx[k], col) = est.study[k].corr;
	}
	using namespace cpp11::literals;
	return writable::list({
		"betaEst"_nm = as_doubles_matrix(beta),
		"h2"_nm = as_doubles(est.h2),
		"rhoEst"_nm = as_doubles_matrix(est.rho),
		"corrEst"_nm = as_doubles_matrix(corr),
		"diverged"_nm = cpp11::as_sexp(est.diverged())
	});
}

// The studies' marginal effects, union places, sample sizes and kinds as R
// passes them, checked against each other.
struct SdprxInputs {
	std::vector<vec> beta_mrg;
	std::vector<uvec> union_idx;
	SdprSharing sharing;
};

SdprxInputs sdprx_inputs(const list& bhatList, const list& ldList, const integers& n,
                         const list& unionIdxList, int nUnion, const list& kinds,
                         sexp rho) {
	SdprxInputs in;
	in.beta_mrg = study_effects(bhatList);
	in.union_idx = union_index(unionIdxList, in.beta_mrg, nUnion, "effects");
	uword K = in.beta_mrg.size();
	if (ldList.size() != static_cast<R_xlen_t>(K) || n.size() != static_cast<R_xlen_t>(K)) {
		cpp11::stop("Each of the %d studies needs one LD and one sample size.", (int)K);
	}
	for (int nk : n) {
		if (nk < 1) {
			cpp11::stop("Each study's sample size must be positive.");
		}
	}
	in.sharing = sdprx_sharing(kinds, rho, K);
	return in;
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
	std::vector<std::vector<SdprBlock<DenseLdBlock>>> blocks = {
		sdpr_dense_blocks(data, std::move(ref_ld_mat), a, n, optLlk)
	};
	return sdpr_one(blocks, bhatR, M, a0k, b0k, schedule, verbose, seed,
	                chain_options(raoBlackwell, track));
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
	std::vector<std::vector<SdprBlock<FactorLdBlock>>> blocks = {
		{sdpr_factor_block(data, F, a, n, optLlk)}
	};
	return sdpr_one(blocks, bhatR, M, a0k, b0k, schedule, verbose, seed,
	                chain_options(raoBlackwell, track));
}

/**
 * @brief SDPRX: SDPR across studies from each study's LD matrix.
 *
 * K studies' SDPR chain whose non-null components each belong to a kind, a
 * set of the studies (sdpr_mcmc()). `bhatList` holds each study's marginal
 * effects, `ldList` its LD matrix (one block), `n` its sample size and
 * `unionIdxList` its variants' 1-based places in the union of the studies'
 * variants, of size `nUnion`, in one allele frame. `kinds` lists each kind's
 * 1-based studies, and `rho` is the studies' K x K correlation matrix, or
 * NULL to learn it. The likelihood is SDPR's optLlk = 1; the rest are
 * sdprCpp()'s. Returns the posterior means: each study's effects over the
 * union (`betaEst`, union x studies, zero where a study lacks a variant), its
 * heritability `h2`, the correlations `rhoEst`, and the fitted correlations
 * and divergence every chain tracks.
 */
[[cpp11::register]]
cpp11::writable::list sdprxCpp(list bhatList, list ldList, integers n, list unionIdxList,
                               int nUnion, list kinds, sexp rho, double a, double c,
                               int M, double a0k, double b0k, int iter, int burn,
                               int thin, bool verbose, sexp seed, bool raoBlackwell,
                               bool track) {
	ChainSchedule schedule(iter, burn, thin);
	SdprxInputs in = sdprx_inputs(bhatList, ldList, n, unionIdxList, nUnion, kinds, rho);
	std::vector<std::vector<SdprBlock<DenseLdBlock>>> blocks;
	for (uword k = 0; k < in.beta_mrg.size(); ++k) {
		R_xlen_t at = static_cast<R_xlen_t>(k);
		mat R = as_Mat(doubles_matrix<>(ldList[at]));
		if (R.n_rows != in.beta_mrg[k].n_elem || R.n_cols != in.beta_mrg[k].n_elem) {
			cpp11::stop("Study %d's LD is %d x %d but it has %d effects.", (int)k + 1,
			            (int)R.n_rows, (int)R.n_cols, (int)in.beta_mrg[k].n_elem);
		}
		std::vector<mat> one;
		one.push_back(std::move(R));
		SdprData data = sdpr_data(doubles(bhatList[at]), n[at], R_NilValue, R_NilValue, c);
		blocks.push_back(sdpr_dense_blocks(data, std::move(one), a, n[at], 1));
	}
	return sdprx_run(blocks, in.beta_mrg, in.union_idx, nUnion, in.sharing, M, a0k, b0k,
	                 schedule, verbose, seed, chain_options(raoBlackwell, track));
}

/**
 * @brief SDPRX from a factor F of each study's LD, R = F'F / (m - 1), never
 * forming R.
 *
 * The arguments are sdprxCpp()'s with `factorList` in place of the LD
 * matrices; each study's matrices take the low-rank form sdprFactorCpp()
 * reads.
 */
[[cpp11::register]]
cpp11::writable::list sdprxFactorCpp(list bhatList, list factorList, integers n,
                                     list unionIdxList, int nUnion, list kinds, sexp rho,
                                     double a, double c, int M, double a0k, double b0k,
                                     int iter, int burn, int thin, bool verbose,
                                     sexp seed, bool raoBlackwell, bool track) {
	ChainSchedule schedule(iter, burn, thin);
	SdprxInputs in = sdprx_inputs(bhatList, factorList, n, unionIdxList, nUnion, kinds,
	                              rho);
	std::vector<std::vector<SdprBlock<FactorLdBlock>>> blocks;
	for (uword k = 0; k < in.beta_mrg.size(); ++k) {
		R_xlen_t at = static_cast<R_xlen_t>(k);
		mat F = checked_factor(doubles_matrix<>(factorList[at]),
		                       static_cast<R_xlen_t>(in.beta_mrg[k].n_elem));
		SdprData data = sdpr_data(doubles(bhatList[at]), n[at], R_NilValue, R_NilValue, c);
		blocks.push_back({sdpr_factor_block(data, F, a, n[at], 1)});
	}
	return sdprx_run(blocks, in.beta_mrg, in.union_idx, nUnion, in.sharing, M, a0k, b0k,
	                 schedule, verbose, seed, chain_options(raoBlackwell, track));
}
