// rss_mcmc.h -- the summary-statistics samplers: what they share, and each
// prior on it.
//
// The samplers draw the effects beta of a window's variants from their
// posterior under the RSS likelihood of the marginal effects,
// bhat ~ N(R beta, R / n), each with its own prior. What does not depend on
// the prior comes first, so a prior is one Gibbs sweep and what a kept draw
// adds to the posterior means:
//   * run_chain() runs the sweeps, keeps the draws past the burn-in at the
//     thinning interval, reports progress, lets the user interrupt and ends
//     a chain that has stopped;
//   * ChainOptions are what every sampler can do beyond its prior:
//     Rao-Blackwellized means, and tracking for combining several chains;
//   * ChainEstimates are what every chain reports, which keep_effects()
//     and check_divergence() update;
//   * BlockLayout places the variants in contiguous LD blocks;
//   * DenseLdBlock and FactorLdBlock hold a block's matrix S -- the LD R, or
//     a matrix a likelihood derives from it -- whole or as a factor W with
//     S = W'W, and answer what the samplers ask of it: its diagonal,
//     products, quadratic forms, principal submatrices and the joint
//     Gaussian draw of the block's effects;
//   * gaussian_draw() draws from a Gaussian given its precision, and
//     standard_normals(), gigrnd(), dirichlet_plus_one() and
//     beta_distribution add the variates <random> lacks.
// Then the priors, each a chain on the above:
//   * PRS-CS's continuous shrinkage prior (PrsCsChain), and across studies
//     as PRS-CSx and JointPRS share it;
//   * SDPR's Dirichlet process prior (SdprData and SdprBlock here; its
//     chain is in rss_mcmc.cpp), and across studies as SDPRX shares it.
// rss_mcmc.cpp has every sampler's cpp11 entry points.
//
// Derived from these projects, each under its own license (the notices are
// in inst/COPYRIGHTS):
//   * PRS-CS (https://github.com/getian107/PRScs), Copyright (c) 2018 Tian
//     Ge, distributed under the MIT License: PrsCsChain, gigrnd(), and
//     DenseLdBlock::draw(), which follows PRS-CS's draw of a block's
//     effects;
//   * PRS-CSx (https://github.com/getian107/PRScsx), Copyright (c) 2019 Tian
//     Ge, and JointPRS (https://github.com/LeqiXu/JointPRS), Copyright (c)
//     2024 Leqi Xu, both distributed under the MIT License: PrsCsChain's
//     coupling across studies;
//   * SDPR (https://github.com/eldronzhou/SDPR), Copyright (C) the SDPR
//     authors (Geyu Zhou and Hongyu Zhao), distributed under the GNU General
//     Public License, version 3: the SDPR section, and gaussian_draw(),
//     which follows the joint draw of SDPR's sample_beta();
//   * SDPRX (https://github.com/eldronzhou/SDPRX), Copyright (c) 2022 Geyu
//     Zhou, distributed under the MIT License: the SDPR chain's kinds of
//     components across studies.
#pragma once

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include "rss_ld.h"

// ==== Shared ====================================================================

// ---- Random variates --------------------------------------------------------

// n standard normal draws from `normal`, in order.
inline arma::vec standard_normals(arma::uword n,
                                  std::normal_distribution<double>& normal,
                                  std::mt19937& rng) {
	arma::vec z(n);
	for (arma::uword i = 0; i < n; ++i) {
		z(i) = normal(rng);
	}
	return z;
}

// The functions of Devroye's (2014) generator of the generalized inverse
// Gaussian: the log density of the transformed variate, its derivative, and
// the envelope.
namespace gig_detail {

inline double psi(double x, double alpha, double lambda) {
	return -alpha * (std::cosh(x) - 1.0) - lambda * (std::exp(x) - x - 1.0);
}

inline double dpsi(double x, double alpha, double lambda) {
	return -alpha * std::sinh(x) - lambda * (std::exp(x) - 1.0);
}

inline double envelope(double x, double sd, double td, double f1, double f2) {
	if (x >= -sd && x <= td) {
		return 1.0;
	} else if (x > td) {
		return f1;
	}
	return f2;
}

}  // namespace gig_detail

// A draw from the generalized inverse Gaussian GIG(p, a, b), with density
// proportional to x^(p - 1) exp(-(a x + b / x) / 2), by Devroye's (2014)
// rejection sampler.
inline double gigrnd(double p, double a, double b, std::mt19937& rng) {
	using gig_detail::psi;
	using gig_detail::dpsi;
	double lambda = p;
	double omega = std::sqrt(a * b);
	bool swap = false;
	if (lambda < 0) {
		lambda = -lambda;
		swap = true;
	}

	double alpha = std::sqrt(std::pow(omega, 2) + std::pow(lambda, 2)) - lambda;

	double x = -psi(1.0, alpha, lambda);
	double t;
	if (x >= 0.5 && x <= 2.0) {
		t = 1.0;
	} else if (x > 2.0) {
		if (alpha == 0 && lambda == 0) {
			t = 1.0;
		} else {
			t = std::sqrt(2.0 / (alpha + lambda));
		}
	} else {
		if (alpha == 0 && lambda == 0) {
			t = 1.0;
		} else {
			t = std::log(4.0 / (alpha + 2.0 * lambda));
		}
	}

	x = -psi(-1.0, alpha, lambda);
	double s;
	if (x >= 0.5 && x <= 2.0) {
		s = 1.0;
	} else if (x > 2.0) {
		if (alpha == 0 && lambda == 0) {
			s = 1.0;
		} else {
			s = std::sqrt(4.0 / (alpha * std::cosh(1.0) + lambda));
		}
	} else {
		if (alpha == 0 && lambda == 0) {
			s = 1.0;
		} else if (alpha == 0) {
			s = 1.0 / lambda;
		} else if (lambda == 0) {
			s = std::log(1.0 + 1.0 / alpha + std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha));
		} else {
			s = std::min(1.0 / lambda, std::log(1.0 + 1.0 / alpha +
			                                    std::sqrt(1.0 / std::pow(alpha, 2) + 2.0 / alpha)));
		}
	}

	double eta = -psi(t, alpha, lambda);
	double zeta = -dpsi(t, alpha, lambda);
	double theta = -psi(-s, alpha, lambda);
	double xi = dpsi(-s, alpha, lambda);

	double p_r = 1.0 / xi;
	double r = 1.0 / zeta;

	double td = t - r * eta;
	double sd = s - p_r * theta;
	double q = td + sd;

	std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
	double rnd = 0.0;
	while (true) {
		double U = unif_dist(rng);
		double V = unif_dist(rng);
		double W = unif_dist(rng);

		if (U < q / (p_r + q + r)) {
			rnd = -sd + q * V;
		} else if (U < (q + r) / (p_r + q + r)) {
			rnd = td - r * std::log(V);
		} else {
			rnd = -sd + p_r * std::log(V);
		}

		double f1 = std::exp(-eta - zeta * (rnd - t));
		double f2 = std::exp(-theta + xi * (rnd + s));

		if (W * gig_detail::envelope(rnd, sd, td, f1, f2) <=
		    std::exp(psi(rnd, alpha, lambda))) {
			break;
		}
	}

	rnd = std::exp(rnd) * (lambda / omega + std::sqrt(1.0 + std::pow(lambda / omega, 2)));

	if (swap) {
		rnd = 1.0 / rnd;
	}

	double result = rnd / std::sqrt(a / b);

	// An exact zero (possible at extreme parameters) would be an infinite
	// precision downstream.
	if (result == 0.0) {
		result = std::numeric_limits<double>::min();
	}

	return result;
}

// A draw from the Dirichlet with parameters counts + 1, through gamma draws.
inline std::vector<double> dirichlet_plus_one(const std::vector<double>& counts,
                                              std::mt19937& rng) {
	std::vector<double> draw(counts.size());
	double total = 0.0;
	for (std::size_t j = 0; j < counts.size(); ++j) {
		std::gamma_distribution<double> gamma(counts[j] + 1.0, 1.0);
		draw[j] = gamma(rng);
		total += draw[j];
	}
	for (double& d : draw) {
		d /= total;
	}
	return draw;
}

// Beta(alpha, beta) as the ratio of two gamma draws, x / (x + y).
class beta_distribution {
public:
	beta_distribution(double alpha, double beta)
	    : a_distribution(alpha, 1.0), b_distribution(beta, 1.0) {}

	template <typename RNG>
	double operator()(RNG& rng) {
		double x = a_distribution(rng);
		double y = b_distribution(rng);
		return x / (x + y);
	}

private:
	std::gamma_distribution<double> a_distribution;
	std::gamma_distribution<double> b_distribution;
};

// A draw from N(P^-1 a, scale P^-1) for a symmetric positive definite P,
// from the standard normals z: with P = L L', L'^-1 (L^-1 a + sqrt(scale) z).
inline arma::vec gaussian_draw(const arma::mat& P, const arma::vec& a,
                               const arma::vec& z, double scale = 1.0) {
	arma::mat L = arma::chol(P, "lower");
	return arma::solve(arma::trimatu(L.t()),
	                   arma::solve(arma::trimatl(L), a) + std::sqrt(scale) * z);
}

// The same draw, with its mean P^-1 a as well: the draw is computed as
// gaussian_draw() computes it, so the two agree bit for bit.
inline arma::vec gaussian_draw(const arma::mat& P, const arma::vec& a,
                               const arma::vec& z, double scale, arma::vec& mean) {
	arma::mat L = arma::chol(P, "lower");
	arma::vec x = arma::solve(arma::trimatl(L), a);
	mean = arma::solve(arma::trimatu(L.t()), x);
	return arma::solve(arma::trimatu(L.t()), x + std::sqrt(scale) * z);
}

// What every sampler can do beyond its prior: report Rao-Blackwellized
// means, the average over the kept iterations of each effect's conditional
// posterior mean, instead of the draws' average; and when several chains run
// (`track`), stop a chain once it diverges and average the fitted
// correlations over the kept iterations, for the filter that combines the
// chains.
struct ChainOptions {
	bool rao_blackwell = false;
	bool track = false;
};

// What every chain reports: its posterior mean effects and, when tracked,
// its mean fitted correlations and whether it diverged.
struct ChainEstimates {
	arma::vec beta;
	arma::vec corr;
	bool diverged = false;
};

// The bound of LDpred2's test of divergence: twice the marginal effects' sum
// of squares.
inline double divergence_gap(const arma::vec& bhat) {
	return 2.0 * arma::dot(bhat, bhat);
}

// A tracked chain's kept draw tested for divergence: its effects' squares
// summing past `gap0`, bigsnpr's test for LDpred2. bigsnpr tests every sweep
// from the first; here the test starts past the burn-in, since some priors
// start far from the posterior (SDPR's first sweeps overshoot the bound many
// times over, then settle well below it). Returns whether the chain has
// diverged, which stops it.
inline bool check_divergence(ChainEstimates& est, const ChainOptions& opts,
                             const arma::vec& beta, double gap0) {
	if (opts.track && arma::dot(beta, beta) > gap0) {
		est.diverged = true;
	}
	return est.diverged;
}

// A kept iteration's effects added to the posterior means: the conditional
// means when Rao-Blackwellized, else the draws.
inline void keep_effects(ChainEstimates& est, const ChainOptions& opts,
                         const arma::vec& draws, const arma::vec& cond_means,
                         int n_kept) {
	est.beta += (opts.rao_blackwell ? cond_means : draws) / n_kept;
}

// The sizes of `blocks`, for a BlockLayout.
template <class Block>
std::vector<arma::uword> block_sizes(const std::vector<Block>& blocks) {
	std::vector<arma::uword> sizes;
	for (const Block& block : blocks) {
		sizes.push_back(block.size());
	}
	return sizes;
}

// ---- The chain --------------------------------------------------------------

// Iterations 1, ..., n_iter, of which the multiples of thin past the burn-in
// are kept for the posterior means.
class ChainSchedule {
public:
	ChainSchedule(int n_iter, int n_burnin, int thin)
	    : n_iter_(n_iter), n_burnin_(n_burnin), thin_(thin) {
		if (n_burnin < 0 || thin < 1) {
			cpp11::stop("The burn-in must be non-negative and the thinning "
			            "interval positive.");
		}
		if (n_kept() < 1) {
			cpp11::stop("The chain keeps no draws: it must run at least one "
			            "thinning interval past the burn-in.");
		}
	}

	int n_iter() const { return n_iter_; }

	// How many iterations are kept: the multiples of thin in (n_burnin, n_iter].
	int n_kept() const { return n_iter_ / thin_ - n_burnin_ / thin_; }

	bool keeps(int itr) const { return itr > n_burnin_ && itr % thin_ == 0; }

private:
	int n_iter_;
	int n_burnin_;
	int thin_;
};

// Runs `chain` through `schedule`: each iteration is one chain.sweep(), then
// chain.keep(n_kept) when the schedule keeps it and, when verbose,
// chain.report(itr, n_iter) every 100th iteration. A chain that has stopped
// (chain.stopped(), as a tracked chain that diverged has) ends the run. The
// user can interrupt between iterations.
template <class Chain>
void run_chain(Chain& chain, const ChainSchedule& schedule, bool verbose) {
	const int n_kept = schedule.n_kept();
	for (int itr = 1; itr <= schedule.n_iter(); ++itr) {
		cpp11::check_user_interrupt();
		chain.sweep();
		if (schedule.keeps(itr)) {
			chain.keep(n_kept);
		}
		if (chain.stopped()) {
			break;
		}
		if (verbose && itr % 100 == 0) {
			chain.report(itr, schedule.n_iter());
		}
	}
}

// ---- LD blocks --------------------------------------------------------------

// The variants laid out in contiguous blocks: block k holds size(k) of them
// from start(k) on.
class BlockLayout {
public:
	explicit BlockLayout(const std::vector<arma::uword>& sizes) : size_(sizes) {
		arma::uword start = 0;
		for (arma::uword size : sizes) {
			start_.push_back(start);
			start += size;
		}
	}

	arma::uword n_blocks() const { return size_.size(); }
	arma::uword start(arma::uword k) const { return start_[k]; }
	arma::uword size(arma::uword k) const { return size_[k]; }

	// Block k's variants; the block must not be empty.
	arma::span span(arma::uword k) const {
		return arma::span(start_[k], start_[k] + size_[k] - 1);
	}

private:
	std::vector<arma::uword> size_;
	std::vector<arma::uword> start_;
};

// A block's symmetric positive semidefinite p x p matrix S held whole.
class DenseLdBlock {
public:
	explicit DenseLdBlock(arma::mat S) : S_(std::move(S)), diag_(S_.diag()) {}

	arma::uword size() const { return S_.n_rows; }
	const arma::vec& diag() const { return diag_; }
	arma::vec times(const arma::vec& x) const { return S_ * x; }
	double quad(const arma::vec& x) const { return arma::dot(S_ * x, x); }
	arma::mat sub(const arma::uvec& idx) const { return S_.submat(idx, idx); }

	// Draws the block's effects, beta ~ N(P^-1 a, scale P^-1) with
	// P = S + diag(1 / v) for prior variances v (gaussian_draw()), and
	// returns beta'P beta. With `mean`, the mean P^-1 a too.
	double draw(const arma::vec& v, const arma::vec& a, double scale,
	            std::normal_distribution<double>& normal, std::mt19937& rng,
	            arma::vec& beta, arma::vec* mean = nullptr) const {
		arma::mat P = S_ + arma::diagmat(1.0 / v);
		arma::vec z = standard_normals(S_.n_rows, normal, rng);
		beta = mean ? gaussian_draw(P, a, z, scale, *mean) : gaussian_draw(P, a, z, scale);
		return arma::as_scalar(beta.t() * P * beta);
	}

private:
	arma::mat S_;
	arma::vec diag_;
};

// A block's symmetric positive semidefinite p x p matrix S = W'W + e I held
// as its factor W (r x p) and the ridge e, so that nothing p x p is formed: a
// product costs O(rp). The ridge is what an LD shrunk toward the identity,
// c R + (1 - c) I, adds to the factor of c R; it is 0 otherwise.
class FactorLdBlock {
public:
	explicit FactorLdBlock(arma::mat W, double extra = 0.0)
	    : W_(std::move(W)), extra_(extra),
	      diag_(arma::sum(arma::square(W_), 0).t() + extra) {}

	arma::uword size() const { return W_.n_cols; }
	const arma::vec& diag() const { return diag_; }
	arma::vec times(const arma::vec& x) const {
		arma::vec Sx = W_.t() * (W_ * x);
		if (extra_ != 0.0) {
			Sx += extra_ * x;
		}
		return Sx;
	}

	double quad(const arma::vec& x) const {
		arma::vec Wx = W_ * x;
		return arma::dot(Wx, Wx) + extra_ * arma::dot(x, x);
	}

	arma::mat sub(const arma::uvec& idx) const {
		arma::mat W_idx = W_.cols(idx);
		arma::mat S = W_idx.t() * W_idx;
		if (extra_ != 0.0) {
			S.diag() += extra_;
		}
		return S;
	}

	// The same draw as DenseLdBlock's without forming S or P. With D = diag(v)
	// and M = W D W' + I (r x r), Woodbury gives P^-1 = D - D W' M^-1 W D, so
	// the mean P^-1 a costs one solve against M, and the noise comes from the
	// sampler of Bhattacharya, Chakraborty and Mallick (2016): u ~ N(0, D),
	// delta ~ N(0, I), theta = u - D W' M^-1 (W u + delta) is N(0, P^-1). One
	// Cholesky factor of M serves both, at O(r^2 p + r^3) per draw instead of
	// the O(p^3) of factoring P. The draws consume the random stream
	// differently from DenseLdBlock's, so the two agree in distribution, not
	// draw for draw.
	// The ridge joins the prior's diagonal: P = W'W + diag(1 / v'), with
	// 1 / v' = 1 / v + e. With `mean`, the mean P^-1 a is returned too.
	double draw(const arma::vec& v, const arma::vec& a, double scale,
	            std::normal_distribution<double>& normal, std::mt19937& rng,
	            arma::vec& beta, arma::vec* mean = nullptr) const {
		const arma::vec vp = extra_ == 0.0 ? v : 1.0 / (1.0 / v + extra_);
		arma::uword p = W_.n_cols;
		arma::uword r = W_.n_rows;
		arma::mat WD = W_.each_row() % vp.t();
		arma::mat M = WD * W_.t();
		M.diag() += 1.0;
		arma::mat L = arma::chol(M, "lower");
		arma::vec centre = vp % a - WD.t() * solve_chol(L, WD * a);
		arma::vec u = arma::sqrt(vp) % standard_normals(p, normal, rng);
		arma::vec delta = standard_normals(r, normal, rng);
		arma::vec theta = u - WD.t() * solve_chol(L, W_ * u + delta);
		beta = centre + std::sqrt(scale) * theta;
		if (mean) {
			*mean = centre;
		}
		arma::vec Wbeta = W_ * beta;
		return arma::dot(Wbeta, Wbeta) + arma::sum(arma::square(beta) / vp);
	}

private:
	arma::mat W_;
	double extra_;
	arma::vec diag_;

	// M^-1 x from the lower Cholesky factor L of M.
	static arma::vec solve_chol(const arma::mat& L, const arma::vec& x) {
		return arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), x));
	}
};

// ==== PRS-CS, PRS-CSx and JointPRS ============================================

// PRS-CS's prior is beta_j ~ N(0, sigma / n psi_j) with the
// Strawderman-Berger local scales psi_j ~ G(a, delta_j), delta_j ~ G(b, phi),
// and the global scale phi either fixed or learned under a half-Cauchy
// prior. One chain runs it for one study, or for K studies coupled as JointPRS
// couples them: the effects b_u = (sqrt(n_k) beta_uk) of the studies carrying
// variant u are drawn from N(0, psi_u C), C's diagonal the studies' residual
// variances sigma_k and its off-diagonal rho_kk' sqrt(sigma_k sigma_k'). The
// local scale psi_u, over the union of the studies' variants, and the global
// scale phi are shared, as PRS-CSx shares them; at rho = 0 the chain is
// PRS-CSx's, and for one study PRS-CS's. Each study's prior on a variant is
// read through the inverse of the correlation matrix of the studies carrying
// it, R_g^-1, the identity at rho = 0, so that no term the coupling adds moves
// a value there.
// A sweep draws each study's blocks jointly given psi and the others' current
// effects (the block's draw(); the prior's conditional is a precision and a
// shift of the marginal effects) and its sigma, as JointPRS draws them; then
// rho by JointPRS's random-walk Metropolis step on [0, min(c_k c_k', 0.99)]
// for the studies' constraints c, keeping the correlation matrix positive
// definite (a constraint of 0 holds a study's correlations at 0); then the
// shared delta, psi and, when it is learned, phi.

// The posterior means: each study's effects over its own variants, with the
// fitted correlations and divergence every chain reports; the shared local
// scales over the union; the studies' residual variances and correlation; and
// the global scale.
struct PrsCsEstimates {
	std::vector<ChainEstimates> study;
	arma::vec psi;
	arma::vec sigma;
	arma::mat rho;
	double phi = 0.0;

	bool diverged() const {
		for (const ChainEstimates& s : study) {
			if (s.diverged) {
				return true;
			}
		}
		return false;
	}
};

// The chain on study k's marginal effects `beta_mrg[k]` of the variants in
// `blocks[k]` (each a DenseLdBlock of R or each a FactorLdBlock with
// R = W'W), placed in the union by `union_idx[k]`.
template <class LdBlock>
class PrsCsChain {
public:
	PrsCsChain(double a, double b, bool learn_phi, double phi,
	           const arma::vec& constraint,
	           const std::vector<arma::vec>& beta_mrg, const arma::vec& n,
	           const std::vector<std::vector<LdBlock>>& blocks,
	           const std::vector<arma::uvec>& union_idx, arma::uword n_union,
	           unsigned int seed, const ChainOptions& opts = ChainOptions())
	    : a_(a), b_(b), learn_phi_(learn_phi), constraint_(constraint),
	      beta_mrg_(beta_mrg), n_(n), blocks_(blocks), union_idx_(union_idx),
	      n_union_(n_union), k_(beta_mrg.size()), rng_(seed),
	      normal_(0.0, 1.0), psi_(n_union, arma::fill::ones),
	      phi_(learn_phi ? 1.0 : phi), sigma_(k_, arma::fill::ones),
	      rho_(k_, k_, arma::fill::eye),
	      beta_union_(n_union, k_, arma::fill::zeros), opts_(opts) {
		for (arma::uword k = 0; k < k_; ++k) {
			arma::uword p = beta_mrg[k].n_elem;
			layouts_.emplace_back(block_sizes(blocks[k]));
			beta_.emplace_back(p, arma::fill::zeros);
			beta_mean_.emplace_back(p, arma::fill::zeros);
			gamma_sigma_.emplace_back((n(k) + p) / 2.0, 1.0);
			gap0_.push_back(divergence_gap(beta_mrg[k]));
			ChainEstimates est;
			est.beta.zeros(p);
			est.corr.zeros(p);
			est_.study.push_back(est);
		}
		build_patterns();
		refresh_rinv();
		est_.psi.zeros(n_union);
		est_.sigma.zeros(k_);
		est_.rho.zeros(k_, k_);
	}

	void sweep() {
		for (arma::uword k = 0; k < k_; ++k) {
			draw_study(k);
		}
		if (update_rho()) {
			refresh_rinv();
		}
		update_scales();
	}

	// A kept iteration: none is kept once any study's draws diverge, which
	// stops the chain.
	void keep(int n_kept) {
		bool diverged = false;
		for (arma::uword k = 0; k < k_; ++k) {
			diverged = check_divergence(est_.study[k], opts_, beta_[k], gap0_[k]) ||
			           diverged;
		}
		if (diverged) {
			return;
		}
		for (arma::uword k = 0; k < k_; ++k) {
			keep_effects(est_.study[k], opts_, beta_[k], beta_mean_[k], n_kept);
			if (opts_.track) {
				keep_fitted(k, n_kept);
			}
		}
		est_.psi += psi_ / n_kept;
		est_.sigma += sigma_ / n_kept;
		est_.rho += rho_ / n_kept;
		est_.phi += phi_ / n_kept;
	}

	void report(int itr, int n_iter) const {
		Rprintf("Iteration %4d of %d\n", itr, n_iter);
	}

	bool stopped() const { return est_.diverged(); }

	const PrsCsEstimates& estimates() const { return est_; }

private:
	double a_;
	double b_;
	bool learn_phi_;
	arma::vec constraint_;
	std::vector<arma::vec> beta_mrg_;
	arma::vec n_;
	const std::vector<std::vector<LdBlock>>& blocks_;
	std::vector<arma::uvec> union_idx_;
	arma::uword n_union_;
	arma::uword k_;
	std::mt19937 rng_;
	std::normal_distribution<double> normal_;
	arma::vec psi_;
	double phi_;
	arma::vec sigma_;
	arma::mat rho_;
	arma::mat beta_union_;  // each study's current effects over the union
	ChainOptions opts_;
	std::vector<BlockLayout> layouts_;
	std::vector<arma::vec> beta_;
	std::vector<arma::vec> beta_mean_;  // the conditional means, for Rao-Blackwell
	std::vector<std::gamma_distribution<double>> gamma_sigma_;
	std::vector<double> gap0_;
	// Each union variant's pattern: the studies carrying it, in order.
	std::vector<arma::uvec> pattern_studies_;
	arma::uvec pattern_of_;
	// Each study's place in the pattern of each of its variants.
	std::vector<arma::uvec> self_;
	// Each pattern's inverse correlation matrix, the identity at rho = 0.
	std::vector<arma::mat> rinv_;
	// The pairs of studies sharing variants, and the variants they share.
	std::vector<std::pair<arma::uword, arma::uword>> pairs_;
	std::vector<arma::uvec> shared_;
	PrsCsEstimates est_;

	// JointPRS's step for the random walk on a correlation.
	static constexpr double rho_step = 0.02;

	void build_patterns() {
		std::vector<std::vector<arma::uword>> carriers(n_union_);
		for (arma::uword k = 0; k < k_; ++k) {
			for (arma::uword u : union_idx_[k]) {
				carriers[u].push_back(k);
			}
		}
		std::map<std::vector<arma::uword>, arma::uword> seen;
		pattern_of_.set_size(n_union_);
		for (arma::uword u = 0; u < n_union_; ++u) {
			auto it = seen.find(carriers[u]);
			if (it == seen.end()) {
				it = seen.emplace(carriers[u], pattern_studies_.size()).first;
				pattern_studies_.push_back(
				    arma::conv_to<arma::uvec>::from(carriers[u]));
			}
			pattern_of_(u) = it->second;
		}
		for (arma::uword k = 0; k < k_; ++k) {
			arma::uvec self(union_idx_[k].n_elem);
			for (arma::uword i = 0; i < self.n_elem; ++i) {
				const arma::uvec& s = pattern_studies_[pattern_of_(union_idx_[k](i))];
				self(i) = arma::as_scalar(arma::find(s == k, 1));
			}
			self_.push_back(self);
		}
		for (arma::uword k1 = 0; k1 < k_; ++k1) {
			for (arma::uword k2 = k1 + 1; k2 < k_; ++k2) {
				std::vector<arma::uword> both;
				for (arma::uword u = 0; u < n_union_; ++u) {
					const std::vector<arma::uword>& c = carriers[u];
					bool has1 = std::find(c.begin(), c.end(), k1) != c.end();
					bool has2 = std::find(c.begin(), c.end(), k2) != c.end();
					if (has1 && has2) {
						both.push_back(u);
					}
				}
				if (!both.empty()) {
					pairs_.emplace_back(k1, k2);
					shared_.push_back(arma::conv_to<arma::uvec>::from(both));
				}
			}
		}
	}

	// Each pattern's R_g^-1 from the current rho.
	void refresh_rinv() {
		rinv_.clear();
		for (const arma::uvec& s : pattern_studies_) {
			arma::mat R = rho_.submat(s, s);
			rinv_.push_back(s.n_elem == 1 ? arma::mat(1, 1, arma::fill::ones)
			                              : arma::mat(arma::inv_sympd(R)));
		}
	}

	// Study k's blocks drawn given the others' current effects, then its
	// sigma. The prior's conditional for variant u has precision
	// (R_g^-1)_kk / psi_u in PRS-CS's units, and shifts the marginal effect
	// by sqrt(sigma_k / n_k) sum_t (R_g^-1)_kt sqrt(n_t / sigma_t) beta_ut /
	// psi_u over the other studies t carrying u.
	void draw_study(arma::uword k) {
		const arma::uvec& u = union_idx_[k];
		arma::uword p = u.n_elem;
		arma::vec v(p);
		arma::vec shifted(p);
		for (arma::uword i = 0; i < p; ++i) {
			arma::uword g = pattern_of_(u(i));
			const arma::uvec& s = pattern_studies_[g];
			const arma::mat& ri = rinv_[g];
			double shift = 0.0;
			for (arma::uword t = 0; t < s.n_elem; ++t) {
				if (t != self_[k](i)) {
					shift += ri(self_[k](i), t) * std::sqrt(n_(s(t)) / sigma_(s(t))) *
					         beta_union_(u(i), s(t));
				}
			}
			v(i) = psi_(u(i)) / ri(self_[k](i), self_[k](i));
			shifted(i) = beta_mrg_[k](i) -
			             std::sqrt(sigma_(k) / n_(k)) * shift / psi_(u(i));
		}
		double quad = 0.0;
		const BlockLayout& layout = layouts_[k];
		for (arma::uword b = 0; b < layout.n_blocks(); ++b) {
			if (layout.size(b) == 0) {
				continue;
			}
			arma::span blk = layout.span(b);
			arma::vec beta_blk;
			arma::vec mean_blk;
			quad += blocks_[k][b].draw(v(blk), shifted(blk), sigma_(k) / n_(k),
			                           normal_, rng_, beta_blk,
			                           opts_.rao_blackwell ? &mean_blk : nullptr);
			beta_[k](blk) = beta_blk;
			if (opts_.rao_blackwell) {
				beta_mean_[k](blk) = mean_blk;
			}
		}
		arma::uvec col = {k};
		beta_union_.submat(u, col) = beta_[k];
		arma::vec psi_k = psi_(u);
		double prior_psi = arma::sum(arma::pow(beta_[k], 2) / psi_k);
		double prior_v = arma::sum(arma::pow(beta_[k], 2) / v);
		// JointPRS reads sigma's quadratic with the prior's diagonal beta'
		// beta / psi; the draws' was beta'beta / v. Their difference is
		// exactly 0 at rho = 0, where this is PRS-CS's update.
		double err = std::max(
		    n_(k) / 2.0 *
		        (1.0 - 2.0 * arma::dot(beta_[k], beta_mrg_[k]) +
		         (quad + (prior_psi - prior_v))),
		    n_(k) / 2.0 * prior_psi);
		sigma_(k) = err / gamma_sigma_[k](rng_);
	}

	// JointPRS's Metropolis step on each pair's correlation. Whether any
	// moved.
	bool update_rho() {
		bool moved = false;
		std::uniform_real_distribution<double> unif(0.0, 1.0);
		for (std::size_t r = 0; r < pairs_.size(); ++r) {
			arma::uword k1 = pairs_[r].first;
			arma::uword k2 = pairs_[r].second;
			double current = rho_(k1, k2);
			std::uniform_real_distribution<double> proposal(current - rho_step,
			                                                current + rho_step);
			double star = proposal(rng_);
			double bound = std::min(constraint_(k1) * constraint_(k2), 0.99);
			if (star < 0.0 || star > bound) {
				continue;
			}
			arma::mat trial = rho_;
			trial(k1, k2) = star;
			trial(k2, k1) = star;
			if (arma::eig_sym(trial).min() <= 0.0) {
				continue;
			}
			double log_ratio = rho_log_target(r, star) - rho_log_target(r, current);
			if (log_ratio >= 0.0 || std::exp(log_ratio) >= unif(rng_)) {
				rho_ = trial;
				moved = true;
			}
		}
		return moved;
	}

	// The log density, up to a constant, of the effects pair r's studies share
	// at correlation `value`: their scaled effects bivariate normal given psi.
	double rho_log_target(std::size_t r, double value) const {
		arma::uword k1 = pairs_[r].first;
		arma::uword k2 = pairs_[r].second;
		const arma::uvec& sh = shared_[r];
		double s1 = sigma_(k1);
		double s2 = sigma_(k2);
		double n1 = n_(k1);
		double n2 = n_(k2);
		double c = value * std::sqrt(s1 * s2);
		double det = s1 * s2 - c * c;
		arma::vec col1 = beta_union_.col(k1);
		arma::vec col2 = beta_union_.col(k2);
		arma::vec b1 = col1(sh);
		arma::vec b2 = col2(sh);
		arma::vec ps = psi_(sh);
		double q11 = arma::sum(arma::square(b1) / ps);
		double q12 = arma::sum(b1 % b2 / ps);
		double q22 = arma::sum(arma::square(b2) / ps);
		return -0.5 * sh.n_elem * std::log(det) -
		       0.5 * n1 * n2 / det *
		           (q11 * s2 / n2 - 2.0 * q12 * c / std::sqrt(n1 * n2) + q22 * s1 / n1);
	}

	// The shared scales: delta, then psi from each variant's effects in the
	// studies carrying it, b_u'(psi_u C)^-1 b_u's quadratic, then phi.
	void update_scales() {
		arma::vec delta(n_union_);
		for (arma::uword u = 0; u < n_union_; ++u) {
			std::gamma_distribution<double> gamma_delta(a_ + b_, 1.0 / (psi_(u) + phi_));
			delta(u) = gamma_delta(rng_);
		}
		for (arma::uword u = 0; u < n_union_; ++u) {
			const arma::uvec& s = pattern_studies_[pattern_of_(u)];
			psi_(u) = gigrnd(a_ - 0.5 * s.n_elem, 2.0 * delta(u), scaled_quad(u), rng_);
		}
		// psi is capped at 1, as in PRS-CS, against explosive local shrinkage.
		psi_ = arma::clamp(psi_, 0.0, 1.0);
		if (learn_phi_) {
			std::gamma_distribution<double> gamma_w(1.0, 1.0 / (phi_ + 1.0));
			double w = gamma_w(rng_);
			std::gamma_distribution<double> gamma_phi(n_union_ * b_ + 0.5,
			                                          1.0 / (arma::sum(delta) + w));
			phi_ = gamma_phi(rng_);
		}
	}

	// sum_st (R_g^-1)_st sqrt(n_s n_t / (sigma_s sigma_t)) beta_us beta_ut
	// over the studies carrying u, its diagonal terms written as PRS-CS
	// writes n beta^2 / sigma.
	double scaled_quad(arma::uword u) const {
		const arma::uvec& s = pattern_studies_[pattern_of_(u)];
		const arma::mat& ri = rinv_[pattern_of_(u)];
		double xx = 0.0;
		for (arma::uword i = 0; i < s.n_elem; ++i) {
			double bi = beta_union_(u, s(i));
			xx += ri(i, i) * (n_(s(i)) * std::pow(bi, 2) / sigma_(s(i)));
			for (arma::uword j = 0; j < s.n_elem; ++j) {
				if (j != i) {
					xx += ri(i, j) *
					      std::sqrt(n_(s(i)) * n_(s(j)) / (sigma_(s(i)) * sigma_(s(j)))) *
					      bi * beta_union_(u, s(j));
				}
			}
		}
		return xx;
	}

	// Study k's fitted correlations, R beta, added block by block.
	void keep_fitted(arma::uword k, int n_kept) {
		const BlockLayout& layout = layouts_[k];
		for (arma::uword b = 0; b < layout.n_blocks(); ++b) {
			if (layout.size(b) > 0) {
				arma::span blk = layout.span(b);
				est_.study[k].corr(blk) += blocks_[k][b].times(beta_[k](blk)) / n_kept;
			}
		}
	}
};

/**
 * @brief The PRS-CS posterior means of one or more studies' effects.
 *
 * @param a,b The shapes of the local scales' prior.
 * @param learn_phi Whether the global scale is learned; if not, it is `phi`.
 * @param constraint Each study's bound on its correlations: 1 lets them be
 *   learned on [0, 0.99], 0 holds them at 0 (PRS-CSx).
 * @param beta_mrg Each study's marginal effects.
 * @param n Each study's sample size.
 * @param blocks Each study's LD blocks, each a DenseLdBlock or each a
 *   FactorLdBlock.
 * @param union_idx Each study's variants' 0-based places in the union.
 */
template <class LdBlock>
PrsCsEstimates prs_cs_mcmc(double a, double b, bool learn_phi, double phi,
                           const arma::vec& constraint,
                           const std::vector<arma::vec>& beta_mrg,
                           const arma::vec& n,
                           const std::vector<std::vector<LdBlock>>& blocks,
                           const std::vector<arma::uvec>& union_idx,
                           arma::uword n_union, const ChainSchedule& schedule,
                           bool verbose, unsigned int seed,
                           const ChainOptions& opts = ChainOptions()) {
	if (verbose) {
		Rprintf("Running Markov Chain Monte Carlo (MCMC) sampler...\n");
	}
	PrsCsChain<LdBlock> chain(a, b, learn_phi, phi, constraint, beta_mrg, n, blocks,
	                          union_idx, n_union, seed, opts);
	run_chain(chain, schedule, verbose);
	PrsCsEstimates est = chain.estimates();
	if (verbose && learn_phi) {
		Rprintf("Estimated global shrinkage parameter: %g\n", est.phi);
	}
	if (verbose) {
		Rprintf("MCMC sampling completed.\n");
	}
	return est;
}

// ==== SDPR and SDPRX ============================================================

// SDPR's prior on each effect is a Dirichlet process mixture of normals, the
// first component a point mass at zero, and its likelihood is SDPR's
// approximation to the RSS likelihood through the matrices B and A of each
// block. A sweep draws the components' variances, each variant's component,
// the stick-breaking weights and the concentration, the causal variants'
// effects jointly per block, and the effects' global scale eta, which has an
// N(1, 1) prior (the original's is flat).
//
// One chain runs it for one study, or for K studies as SDPRX runs two. Each
// non-null component then belongs to a kind, a set of the studies, and a
// variant in it has effects in those studies drawn from N(0, var Sigma), with
// Sigma the studies' correlations rho restricted to the kind: SDPRX's
// population-specific and shared components. Each kind holds M - 1
// components. The null is the process's first stick, as in SDPR; the rest of
// the weight is split over the kinds by a Dirichlet (SDPRX's over the kinds),
// and each kind's share is broken over its components with the one
// concentration. A variant a study lacks may take any kind, its effect in
// that study integrated out. The studies share eta. rho is given, or each
// pair's is drawn from its conditional on a grid over [0, 0.99]. For one
// study there is one kind, and the chain is SDPR's draw for draw.

// The per-variant data SDPR's likelihood reads.
//
// - beta_mrg: the marginal effects, deflated by `c`.
// - sz: the sample size of each variant.
// - array: the genotyping array of each variant (0 = every cohort, 1 or 2 =
//   the cohorts genotyped on that array). Pairs of variants on arrays 1 and
//   2 share no samples, so the optLlk = 2 likelihood zeroes their LD.
struct SdprData {
	std::vector<double> beta_mrg;
	std::vector<double> sz;
	std::vector<int> array;
};

// What the sampler reads from one LD block, R: B, the matrix of the
// likelihood, n R (R + aI)^-1 R for optLlk = 1 or R B2^-1 R for optLlk = 2
// (B2 the scaled R of sdpr_llk2_matrix()), calc_b_tmp, the matching
// A' beta_mrg, and R itself for beta'R beta (the heritability). B and R are
// both held whole (DenseLdBlock) or both as factors (FactorLdBlock).
template <class LdBlock>
struct SdprBlock {
	LdBlock B;
	LdBlock R;
	arma::vec calc_b_tmp;

	arma::uword size() const { return B.size(); }
};

// The posterior means: each study's effects over its own variants, with the
// fitted correlations and divergence every chain reports, and its
// heritability; and the studies' correlations.
struct SdprEstimates {
	std::vector<ChainEstimates> study;
	arma::vec h2;
	arma::mat rho;

	bool diverged() const {
		for (const ChainEstimates& s : study) {
			if (s.diverged) {
				return true;
			}
		}
		return false;
	}
};

// How K studies' effects are shared: each non-null kind's studies (0-based,
// in order), their correlations rho (K x K), and whether rho is learned,
// starting from the given matrix. One study has the one kind {0}.
struct SdprSharing {
	std::vector<arma::uvec> kinds;
	arma::mat rho;
	bool learn_rho = false;
};

/**
 * @brief The dense LD blocks of the SDPR likelihood (solve_ldmat).
 *
 * @param dat The per-variant data, with beta_mrg already deflated by `c`.
 * @param ref_ld_mat The LD matrix of each block, taken over.
 * @param a Factor to shrink the reference LD matrix (optLlk = 1).
 * @param sz The sample size of the GWAS (optLlk = 1).
 * @param opt_llk 1 for equation 6 (slightly shrink the correlation of
 *   variants) and 2 for equation 5 (variants genotyped on different arrays in
 *   a separate cohort).
 */
std::vector<SdprBlock<DenseLdBlock>> sdpr_dense_blocks(const SdprData &dat,
                                                       std::vector<arma::mat> ref_ld_mat,
                                                       double a, unsigned sz,
                                                       int opt_llk);

/**
 * @brief The one LD block of the SDPR likelihood from a factor of its LD,
 * R = F'F / (m - 1), never forming R or B.
 *
 * The arguments are sdpr_dense_blocks()'s with the factor F (m x p) in place
 * of the LD matrices. B and A' beta_mrg agree with the dense block's to
 * rounding: optLlk = 1 has them in closed form from the eigenvalues of R;
 * optLlk = 2 applies B2^-1 to r + 1 vectors, through the Woodbury identity
 * while B2's low-rank part is narrower than p, and otherwise by forming B2.
 */
SdprBlock<FactorLdBlock> sdpr_factor_block(const SdprData &dat, const arma::mat &F,
                                           double a, unsigned sz, int opt_llk);

/**
 * @brief The SDPR posterior means of one or more studies' effects and
 * heritabilities.
 *
 * @param blocks Each study's LD blocks (sdpr_dense_blocks() or
 *   sdpr_factor_block()), the same number for every study; a variant the
 *   studies share lies in the same block in each.
 * @param union_idx Each study's variants' 0-based places in the union.
 * @param n_union The size of the union.
 * @param sharing The kinds of the non-null components and rho.
 * @param M Max number of variance components of each kind, the null's
 *   among them: SDPR's M for one study.
 * @param a0k,b0k The inverse gamma prior of the components' variances.
 * @param schedule The chain's iterations, burn-in and thinning.
 * @param verbose Whether to report progress.
 * @param seed Seed of the sampler's random number generator.
 * @param opts Rao-Blackwellized means, and tracking for several chains.
 * @param gap0 Each study's twice its squared marginal effects' sum, past
 *   which the effects' squares mark a tracked chain diverged.
 */
template <class LdBlock>
SdprEstimates sdpr_mcmc(const std::vector<std::vector<SdprBlock<LdBlock>>> &blocks,
                        const std::vector<arma::uvec> &union_idx, arma::uword n_union,
                        const SdprSharing &sharing, size_t M, double a0k, double b0k,
                        const ChainSchedule &schedule, bool verbose, unsigned int seed,
                        const ChainOptions &opts, const std::vector<double> &gap0);
