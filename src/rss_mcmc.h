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
//     and check_divergence() update; ChainRandom is the generator of the
//     single-site samplers;
//   * BlockLayout places the variants in contiguous LD blocks;
//   * DenseLdBlock and FactorLdBlock hold a block's matrix S -- the LD R, or
//     a matrix a likelihood derives from it -- whole or as a factor W with
//     S = W'W, and answer what the samplers ask of it: its diagonal,
//     products, quadratic forms, principal submatrices and the joint
//     Gaussian draw of the block's effects;
//   * gaussian_draw() draws from a Gaussian given its precision, and
//     standard_normals(), gigrnd(), rinvgauss(), dirichlet_plus_one() and
//     beta_distribution add the variates <random> lacks.
// Then the priors, each a chain on the above:
//   * PRS-CS's continuous shrinkage prior (PrsCsChain);
//   * SDPR's Dirichlet process prior (SdprData and SdprBlock here; its
//     chain is in rss_mcmc.cpp);
//   * the Bayesian alphabet in qgg's and hibayes' forms, on a factor of the
//     LD (QggFormChain, HibayesFormChain);
//   * LDpred2's point-normal prior, on a factor of the LD (Ldpred2Chain).
// rss_mcmc.cpp has every sampler's cpp11 entry points.
//
// Derived from these projects, each under its own license (the notices are
// in inst/COPYRIGHTS):
//   * PRS-CS (https://github.com/getian107/PRScs), Copyright (c) 2018 Tian
//     Ge, distributed under the MIT License: PrsCsChain, gigrnd(), and
//     DenseLdBlock::draw(), which follows PRS-CS's draw of a block's
//     effects;
//   * SDPR (https://github.com/eldronzhou/SDPR), Copyright (C) the SDPR
//     authors (Geyu Zhou and Hongyu Zhao), distributed under the GNU General
//     Public License, version 3: the SDPR section, and gaussian_draw(),
//     which follows the joint draw of SDPR's sample_beta();
//   * qgg (https://github.com/psoerensen/qgg) and hibayes
//     (https://github.com/YinLiLin/hibayes), both distributed under the GNU
//     General Public License, version 3: the Bayesian alphabet section;
//   * bigsnpr (https://github.com/privefl/bigsnpr), Copyright (c) Florian
//     Prive, distributed under the GNU General Public License, version 3:
//     the LDpred2 section, from its ldpred2.cpp and ldpred2-auto.cpp.
#pragma once

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
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

// A draw from the inverse Gaussian with mean mu and shape lambda, by Michael,
// Schucany and Haas's (1976) transformation, as qgg and hibayes draw it.
inline double rinvgauss(double mu, double lambda, std::mt19937& rng) {
	std::normal_distribution<double> normal(0.0, 1.0);
	double z = normal(rng);
	double z2 = z * z;
	double x = mu + 0.5 * mu * mu * z2 / lambda -
	           0.5 * (mu / lambda) * std::sqrt(4.0 * mu * lambda * z2 + mu * mu * z2 * z2);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	if (unif(rng) <= mu / (mu + x)) {
		return x;
	}
	return mu * mu / x;
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

// The single-site samplers' generator, each variate from a fresh
// distribution.
class ChainRandom {
public:
	explicit ChainRandom(unsigned int seed) : rng_(seed) {}

	double normal(double mean, double sd) {
		std::normal_distribution<double> dist(mean, sd);
		return dist(rng_);
	}

	double chisq(double df) {
		std::chi_squared_distribution<double> dist(df);
		return dist(rng_);
	}

	double unif() {
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		return dist(rng_);
	}

	std::mt19937& engine() { return rng_; }

private:
	std::mt19937 rng_;
};

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

// ==== PRS-CS ====================================================================

// PRS-CS's prior is beta_j ~ N(0, sigma / n psi_j) with the
// Strawderman-Berger local scales psi_j ~ G(a, delta_j), delta_j ~ G(b, phi),
// and the global scale phi either fixed or learned under a half-Cauchy
// prior. A sweep draws each block's beta jointly given psi (the block's
// draw()), then sigma, delta, psi and, when it is learned, phi.

// The posterior means, besides the effects every chain reports.
struct PrsCsEstimates : ChainEstimates {
	arma::vec psi;
	double sigma = 0.0;
	double phi = 0.0;
};

// The sampler's state on the marginal effects `beta_mrg` of the variants in
// `blocks` (each a DenseLdBlock of R or each a FactorLdBlock with R = W'W).
template <class LdBlock>
class PrsCsChain {
public:
	PrsCsChain(double a, double b, bool learn_phi, double phi,
	           const arma::vec& beta_mrg, int n,
	           const std::vector<LdBlock>& blocks, unsigned int seed,
	           const ChainOptions& opts = ChainOptions())
	    : a_(a), b_(b), learn_phi_(learn_phi), n_(n),
	      p_(static_cast<int>(beta_mrg.n_elem)), beta_mrg_(beta_mrg),
	      blocks_(blocks), layout_(block_sizes(blocks)), rng_(seed),
	      normal_(0.0, 1.0), gamma_sigma_((n + p_) / 2.0, 1.0),
	      beta_(p_, arma::fill::zeros), psi_(p_, arma::fill::ones),
	      sigma_(1.0), phi_(learn_phi ? 1.0 : phi), opts_(opts),
	      beta_mean_(p_, arma::fill::zeros),
	      gap0_(divergence_gap(beta_mrg)) {
		est_.beta.zeros(p_);
		est_.psi.zeros(p_);
		est_.corr.zeros(p_);
	}

	void sweep() {
		double quad = 0.0;
		for (arma::uword k = 0; k < layout_.n_blocks(); ++k) {
			if (layout_.size(k) == 0) {
				continue;
			}
			arma::span blk = layout_.span(k);
			arma::vec beta_blk;
			arma::vec mean_blk;
			quad += blocks_[k].draw(psi_(blk), beta_mrg_(blk), sigma_ / n_,
			                        normal_, rng_, beta_blk,
			                        opts_.rao_blackwell ? &mean_blk : nullptr);
			beta_(blk) = beta_blk;
			if (opts_.rao_blackwell) {
				beta_mean_(blk) = mean_blk;
			}
		}
		double err = std::max(n_ / 2.0 * (1.0 - 2.0 * arma::dot(beta_, beta_mrg_) + quad),
		                      n_ / 2.0 * arma::sum(arma::pow(beta_, 2) / psi_));

		// sigma ~ InvGamma((n + p) / 2, rate = err), as err / Gamma((n + p) / 2, 1).
		sigma_ = err / gamma_sigma_(rng_);

		arma::vec delta(p_);
		for (int j = 0; j < p_; ++j) {
			std::gamma_distribution<double> gamma_delta(a_ + b_, 1.0 / (psi_(j) + phi_));
			delta(j) = gamma_delta(rng_);
		}

		for (int j = 0; j < p_; ++j) {
			psi_(j) = gigrnd(a_ - 0.5, 2.0 * delta(j), n_ * std::pow(beta_(j), 2) / sigma_,
			                 rng_);
		}
		// psi is capped at 1, as in PRS-CS, against explosive local shrinkage.
		psi_ = arma::clamp(psi_, 0.0, 1.0);

		if (learn_phi_) {
			std::gamma_distribution<double> gamma_w(1.0, 1.0 / (phi_ + 1.0));
			double w = gamma_w(rng_);
			std::gamma_distribution<double> gamma_phi(p_ * b_ + 0.5,
			                                          1.0 / (arma::sum(delta) + w));
			phi_ = gamma_phi(rng_);
		}
	}

	void keep(int n_kept) {
		if (check_divergence(est_, opts_, beta_, gap0_)) {
			return;
		}
		keep_effects(est_, opts_, beta_, beta_mean_, n_kept);
		est_.psi += psi_ / n_kept;
		est_.sigma += sigma_ / n_kept;
		est_.phi += phi_ / n_kept;
		if (opts_.track) {
			for (arma::uword k = 0; k < layout_.n_blocks(); ++k) {
				if (layout_.size(k) > 0) {
					arma::span blk = layout_.span(k);
					est_.corr(blk) += blocks_[k].times(beta_(blk)) / n_kept;
				}
			}
		}
	}

	void report(int itr, int n_iter) const {
		Rprintf("Iteration %4d of %d\n", itr, n_iter);
	}

	bool stopped() const { return est_.diverged; }

	const PrsCsEstimates& estimates() const { return est_; }

private:
	double a_;
	double b_;
	bool learn_phi_;
	int n_;
	int p_;
	arma::vec beta_mrg_;
	const std::vector<LdBlock>& blocks_;
	BlockLayout layout_;
	std::mt19937 rng_;
	std::normal_distribution<double> normal_;
	std::gamma_distribution<double> gamma_sigma_;
	arma::vec beta_;
	arma::vec psi_;
	double sigma_;
	double phi_;
	ChainOptions opts_;
	arma::vec beta_mean_;  // the blocks' conditional means, for Rao-Blackwell
	double gap0_;
	PrsCsEstimates est_;
};

/**
 * @brief The PRS-CS posterior means of the effects.
 *
 * @param a,b The shapes of the local scales' prior.
 * @param learn_phi Whether the global scale is learned; if not, it is `phi`.
 * @param beta_mrg The marginal effects.
 * @param maf The minor allele frequencies, all zero to keep the effects
 *   standardized; otherwise beta is returned per allele.
 * @param n The sample size.
 * @param blocks The LD blocks, each a DenseLdBlock or each a FactorLdBlock.
 */
template <class LdBlock>
PrsCsEstimates prs_cs_mcmc(double a, double b, bool learn_phi, double phi,
                           const arma::vec& beta_mrg, const arma::vec& maf, int n,
                           const std::vector<LdBlock>& blocks,
                           const ChainSchedule& schedule, bool verbose,
                           unsigned int seed,
                           const ChainOptions& opts = ChainOptions()) {
	if (verbose) {
		Rprintf("Running Markov Chain Monte Carlo (MCMC) sampler...\n");
	}
	PrsCsChain<LdBlock> chain(a, b, learn_phi, phi, beta_mrg, n, blocks, seed, opts);
	run_chain(chain, schedule, verbose);
	PrsCsEstimates est = chain.estimates();

	if (arma::max(maf) > 0) {
		est.beta /= arma::sqrt(2.0 * maf % (1.0 - maf));
	}
	if (verbose && learn_phi) {
		Rprintf("Estimated global shrinkage parameter: %g\n", est.phi);
	}
	if (verbose) {
		Rprintf("MCMC sampling completed.\n");
	}
	return est;
}

// ==== SDPR ======================================================================

// SDPR's prior on each effect is a Dirichlet process mixture of normals, the
// first component a point mass at zero, and its likelihood is SDPR's
// approximation to the RSS likelihood through the matrices B and A of each
// block. A sweep draws the components' variances, each variant's component,
// the stick-breaking weights and the concentration, the causal variants'
// effects jointly per block, and the effects' global scale eta, which has an
// N(1, 1) prior (the original's is flat).

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

// The posterior mean of the heritability, besides the effects every chain
// reports.
struct SdprEstimates : ChainEstimates {
	double h2 = 0.0;
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
 * @brief The SDPR posterior means of the effects and the heritability.
 *
 * @param blocks The LD blocks (sdpr_dense_blocks() or sdpr_factor_block()).
 * @param M Max number of variance components.
 * @param a0k,b0k The inverse gamma prior of the components' variances.
 * @param schedule The chain's iterations, burn-in and thinning.
 * @param verbose Whether to report progress.
 * @param seed Seed of the sampler's random number generator.
 * @param opts Rao-Blackwellized means, and tracking for several chains.
 * @param gap0 Twice the squared marginal effects' sum, past which the
 *   effects' squares mark a tracked chain diverged.
 */
template <class LdBlock>
SdprEstimates sdpr_mcmc(const std::vector<SdprBlock<LdBlock>> &blocks, size_t M,
                        double a0k, double b0k, const ChainSchedule &schedule,
                        bool verbose, unsigned int seed, const ChainOptions &opts,
                        double gap0);

// ==== The Bayesian alphabet =====================================================

// qgg's and hibayes' summary-statistics samplers read the LD as a matrix.
// This one reads it as a factor F (k x p) of R = F'F / (k - 1), keeping
// u = F_w b for the current effects b through rss_ld.h's FactorLd, so a
// single-site update costs O(k) and nothing p x p is formed. The likelihood
// is SBayes's: with per-variant weights w (the sample size n, or qgg's
// w = 1 / (se^2 + b^2 / n)), X'X = W^1/2 R W^1/2 and X'y = w % bhat, both
// held through F_w = F diag(sqrt(w)).
//
// Two families of priors, each sampled as its package samples it:
//   * qgg's (CRAN 1.1.6): bayesA, bayesC, bayesN and bayesR as its
//     summary-statistics sampler sbayes() has them, and bayesL in Legarra et
//     al.'s (2011) form, as its individual-level sampler has it. Scaled
//     inverse chi-square priors with nu = 4 on the marker and residual
//     variances; pi drawn from a Dirichlet.
//   * hibayes': bayesB as BayesBpi, and bayesL in Park and Casella's (2008)
//     form, as its SBayesD() has them. A scaled inverse chi-square prior on
//     the genetic variance and an improper one on the residual variance.

// What every alphabet sampler reads, and the options every sampler shares:
// the LD shrinkage c, under which X'X is W^1/2 (c R + (1 - c) I) W^1/2, so
// only its off-diagonal part is scaled; whether a drawn effect may change
// sign (without, it is zero instead, and leaves its class); and
// ChainOptions.
struct AlphabetData {
	arma::mat Fw;  // F diag(sqrt(w)), k x p: X'X = Fw'Fw / (k - 1)
	arma::vec w;  // the per-variant weights
	arma::vec wy;  // X'y = w % bhat
	arma::vec bhat;
	int n;
	double shrink = 1.0;
	bool allow_jump_sign = true;
	ChainOptions opts;
};

// The posterior means, and when the chain is tracked its mean fitted
// correlations and whether it diverged.
struct AlphabetEstimates : ChainEstimates {
	arma::vec pip;  // share of kept draws with a non-null effect
	double pi = 0.0;  // proportion of non-null effects
	double ve = 0.0;  // residual variance
};

// What the two families share: the effects, the factor that keeps X'X b, the
// generator and the posterior means.
class AlphabetChainBase : protected ChainRandom {
public:
	const AlphabetEstimates& estimates() const { return est_; }

	void report(int itr, int n_iter) const {
		Rprintf("Iteration %4d of %d: residual variance %g\n", itr, n_iter, ve_);
	}

	bool stopped() const { return est_.diverged; }

protected:
	AlphabetChainBase(const AlphabetData& data, unsigned int seed)
	    : ChainRandom(seed), data_(data), p_(data.wy.n_elem),
	      b_(p_, arma::fill::zeros), ld_(data.Fw, 0.0, b_), ve_(0.0),
	      cond_mean_(p_, arma::fill::zeros), gap0_(divergence_gap(data.bhat)) {
		est_.beta.zeros(p_);
		est_.pip.zeros(p_);
		est_.corr.zeros(p_);
	}

	// (X'X)_jj under the LD shrinkage.
	double xx(arma::uword j) const {
		const double c = data_.shrink;
		return c * ld_.diag(j) + (1.0 - c) * data_.w(j);
	}

	// X'y_j - (X'X b)_j + (X'X)_jj b_j: variant j's statistic net of the
	// other effects, whose LD the shrinkage scales.
	double rhs(arma::uword j) const {
		const double c = data_.shrink;
		return data_.wy(j) - c * ld_.row(j, b_) + c * ld_.diag(j) * b_(j);
	}

	// Sets b_j, which without sign jumps is zero where the draw's sign
	// differs from the current effect's; returns what was set.
	double set(arma::uword j, double bn) {
		if (!data_.allow_jump_sign && bn * b_(j) < 0.0) {
			bn = 0.0;
		}
		double diff = bn - b_(j);
		if (diff != 0.0) {
			ld_.update(j, diff);
			b_(j) = bn;
		}
		return bn;
	}

	// X'y - X'X b.
	arma::vec residuals() const {
		arma::vec r = data_.wy - data_.shrink * ld_.product(b_);
		if (data_.shrink != 1.0) {
			r -= (1.0 - data_.shrink) * (data_.w % b_);
		}
		return r;
	}

	// Kept draws add to the means of the effects (Rao-Blackwellized: their
	// conditional means), the inclusions, pi and ve, and when tracked to the
	// fitted correlations X'X b / w.
	void keep_common(int n_kept, const std::vector<int>& d, double pi) {
		if (check_divergence(est_, data_.opts, b_, gap0_)) {
			return;
		}
		keep_effects(est_, data_.opts, b_, cond_mean_, n_kept);
		if (data_.opts.track) {
			est_.corr += (data_.wy - residuals()) / data_.w / n_kept;
		}
		for (arma::uword j = 0; j < d.size(); ++j) {
			if (d[j] != 0) {
				est_.pip(j) += 1.0 / n_kept;
			}
		}
		est_.pi += pi / n_kept;
		est_.ve += ve_ / n_kept;
	}

	const AlphabetData& data_;
	arma::uword p_;
	arma::vec b_;
	FactorLd ld_;
	double ve_;
	arma::vec cond_mean_;  // each effect's conditional mean, for Rao-Blackwell
	double gap0_;
	AlphabetEstimates est_;
};

// qgg's priors: bayesN, bayesA, bayesC, bayesR and Legarra's bayesL.
class QggFormChain : public AlphabetChainBase {
public:
	enum Method { N, A, L, C, R };

	QggFormChain(const AlphabetData& data, Method method, double pi, double h2,
	             unsigned int seed)
	    : AlphabetChainBase(data, seed), method_(method), d_(p_, 1) {
		init_priors(pi, h2);
		// qgg enters the variants with the largest marginal effects first.
		order_.resize(p_);
		std::iota(order_.begin(), order_.end(), 0);
		const arma::vec& bhat = data_.bhat;
		std::stable_sort(order_.begin(), order_.end(), [&](arma::uword i, arma::uword j) {
			return bhat(i) * bhat(i) > bhat(j) * bhat(j);
		});
	}

	void sweep() {
		for (arma::uword j : order_) {
			sweep_one(j);
		}
		if (method_ == L) {
			update_lasso_precisions();
		} else {
			update_marker_variance();
		}
		if (method_ == C || method_ == R) {
			pi_ = dirichlet_plus_one(class_counts(), engine());
		}
		update_residual_variance();
	}

	void keep(int n_kept) {
		keep_common(n_kept, d_, 1.0 - pi_[0]);
	}

private:
	Method method_;
	std::vector<arma::uword> order_;
	std::vector<int> d_;  // each effect's variance class; 0 is null
	std::vector<double> pi_;
	std::vector<double> gamma_;
	double vb_ = 0.0;
	double ssb_prior_ = 0.0;
	double sse_prior_ = 0.0;
	const double nub_ = 4.0;
	const double nue_ = 4.0;
	arma::vec vbi_;  // bayesA: each effect's variance
	arma::vec lambda_;  // bayesL: each effect's prior precision
	double lambda2_ = 0.0;

	// qgg's sbayes_sparse() defaults; bayesR's classes ignore `pi`, which
	// then sets the slab variance at qgg's default of 0.001. qgg's y'y,
	// median((b^2 + (n - 2) se^2) w), is median(w) for standardized
	// statistics.
	void init_priors(double pi, double h2) {
		double vy = arma::median(data_.w) / (data_.n - 1.0);
		ve_ = vy * (1.0 - h2);
		double vg = vy * h2;
		double m = static_cast<double>(p_);
		pi_ = {1.0 - pi, pi};
		gamma_ = {0.0, 1.0};
		if (method_ == R) {
			pi = 0.001;
			pi_ = {0.95, 0.02, 0.02, 0.01};
			gamma_ = {0.0, 0.01, 0.1, 1.0};
		}
		vb_ = (method_ == C || method_ == R) ? vg / (m * pi) : vg / m;
		ssb_prior_ = ((nub_ - 2.0) / nub_) * vb_;
		sse_prior_ = ((nue_ - 2.0) / nue_) * ve_;
		vbi_ = arma::vec(p_, arma::fill::value(vb_));
		lambda_ = arma::vec(p_, arma::fill::value(ve_ / vb_));
		lambda2_ = m;
		if (method_ == N || method_ == A || method_ == L) {
			std::fill(d_.begin(), d_.end(), 1);
		} else {
			std::fill(d_.begin(), d_.end(), 0);
		}
	}

	void sweep_one(arma::uword j) {
		double dj = xx(j);
		double rhs_j = rhs(j);
		switch (method_) {
		case N:
			cond_mean_(j) = slab_mean(rhs_j, dj, vb_);
			set(j, slab_draw(rhs_j, dj, vb_));
			break;
		case A:
			vbi_(j) = (b_(j) * b_(j) + ssb_prior_ * nub_) / chisq(1.0 + nub_);
			cond_mean_(j) = slab_mean(rhs_j, dj, vbi_(j));
			set(j, slab_draw(rhs_j, dj, vbi_(j)));
			break;
		case L: {
			double lhs = dj / ve_ + lambda_(j);
			cond_mean_(j) = (rhs_j / ve_) / lhs;
			set(j, normal((rhs_j / ve_) / lhs, std::sqrt(1.0 / lhs)));
			break;
		}
		case C:
		case R:
			sweep_mixture(j, rhs_j, dj);
			break;
		}
	}

	// b_j's conditional under a normal prior of variance v, and its mean.
	double slab_draw(double rhs_j, double dj, double v) {
		double lhs = dj + ve_ / v;
		return normal(rhs_j / lhs, std::sqrt(ve_ / lhs));
	}

	double slab_mean(double rhs_j, double dj, double v) const {
		return rhs_j / (dj + ve_ / v);
	}

	// b_j's conditional mean over the variance classes, weighted by their
	// posterior probabilities.
	double mixture_mean(const std::vector<double>& loglik, double rhs_j, double dj) const {
		double top = *std::max_element(loglik.begin(), loglik.end());
		double total = 0.0;
		double mean = 0.0;
		for (std::size_t c = 0; c < loglik.size(); ++c) {
			double weight = std::exp(loglik[c] - top);
			total += weight;
			if (c > 0) {
				mean += weight * slab_mean(rhs_j, dj, vb_ * gamma_[c]);
			}
		}
		return mean / total;
	}

	// bayesC and bayesR: the variance class from its marginal likelihood,
	// then the effect within it.
	void sweep_mixture(arma::uword j, double rhs_j, double dj) {
		std::size_t nc = pi_.size();
		std::vector<double> loglik(nc);
		double v0 = dj * ve_;
		loglik[0] = -0.5 * std::log(v0) - 0.5 * rhs_j * rhs_j / v0 + std::log(pi_[0]);
		for (std::size_t c = 1; c < nc; ++c) {
			double v1 = dj * ve_ + dj * dj * vb_ * gamma_[c];
			loglik[c] = -0.5 * std::log(v1) - 0.5 * rhs_j * rhs_j / v1 + std::log(pi_[c]);
		}
		double u = unif();
		double cum = 0.0;
		int cls = 0;
		for (std::size_t c = 0; c < nc; ++c) {
			double denom = 0.0;
			for (std::size_t k = 0; k < nc; ++k) {
				denom += std::exp(loglik[k] - loglik[c]);
			}
			cum += 1.0 / denom;
			if (u < cum) {
				cls = static_cast<int>(c);
				break;
			}
		}
		d_[j] = cls;
		if (data_.opts.rao_blackwell) {
			cond_mean_(j) = mixture_mean(loglik, rhs_j, dj);
		}
		if (set(j, cls == 0 ? 0.0 : slab_draw(rhs_j, dj, vb_ * gamma_[cls])) == 0.0) {
			d_[j] = 0;
		}
	}

	std::vector<double> class_counts() const {
		std::vector<double> counts(pi_.size(), 0.0);
		for (int cls : d_) {
			counts[cls] += 1.0;
		}
		return counts;
	}

	// The slab variance from the effects in the slab, each scaled by its
	// class (bayesN and bayesA have every effect in class 1, of scale 1).
	void update_marker_variance() {
		double ssb = 0.0;
		double dfb = 0.0;
		for (arma::uword j = 0; j < p_; ++j) {
			if (d_[j] != 0) {
				ssb += b_(j) * b_(j) / gamma_[d_[j]];
				dfb += 1.0;
			}
		}
		vb_ = (ssb + ssb_prior_ * nub_) / chisq(dfb + nub_);
	}

	// Legarra et al.: each effect's prior precision from its inverse-Gaussian
	// conditional, then lambda^2 from its gamma conditional, kept in [1, 1e6].
	void update_lasso_precisions() {
		const double eps = 1e-8;
		double rate = std::max(lambda2_, eps);
		double sum_tau = 0.0;
		for (arma::uword j = 0; j < p_; ++j) {
			double abs_b = std::max(std::abs(b_(j)), eps);
			double mu = std::max(std::sqrt(rate) / abs_b, eps);
			lambda_(j) = rinvgauss(mu, rate, engine());
			sum_tau += 1.0 / lambda_(j);
		}
		std::gamma_distribution<double> gamma(static_cast<double>(p_) + 1.0, 2.0 / sum_tau);
		lambda2_ = std::max(1.0, std::min(gamma(engine()), 1e6));
	}

	// y'y - 2 b'X'y + b'X'X b is the residual sum of squares only while the
	// statistics agree with the LD; where they disagree (an LD reference
	// unlike the study's) it can turn negative, which qgg does not guard
	// and which makes its chain NaN. The residual variance then keeps its
	// last value, as hibayes falls back on one too.
	void update_residual_variance() {
		arma::vec r = residuals();
		double sse = arma::median(data_.w) - arma::dot(b_, r + data_.wy);
		double scale = sse + sse_prior_ * nue_;
		if (scale > 0.0) {
			ve_ = scale / chisq(data_.n + nue_);
		}
	}
};

// hibayes' priors: bayesB as BayesBpi and Park and Casella's bayesL.
class HibayesFormChain : public AlphabetChainBase {
public:
	enum Method { B, L };

	HibayesFormChain(const AlphabetData& data, Method method, double pi, double h2,
	                 unsigned int seed)
	    : AlphabetChainBase(data, seed), method_(method), d_(p_, method == L ? 1 : 0) {
		init_priors(pi, h2);
	}

	void sweep() {
		if (method_ == B) {
			for (arma::uword j = 0; j < p_; ++j) {
				sweep_spike_slab(j);
			}
			std::vector<double> counts = {0.0, 0.0};
			for (int flag : d_) {
				counts[flag] += 1.0;
			}
			pi_ = dirichlet_plus_one(counts, engine());
		} else {
			for (arma::uword j = 0; j < p_; ++j) {
				sweep_lasso(j);
			}
			std::gamma_distribution<double> gamma(shape0_ + p_, 1.0 / (rate0_ + arma::sum(vargL_) / 2.0));
			lambda2_ = gamma(engine());
		}
		update_variances();
	}

	void keep(int n_kept) {
		keep_common(n_kept, d_, pi_[1]);
	}

private:
	Method method_;
	std::vector<int> d_;
	std::vector<double> pi_;
	double dfvara_ = 4.0;
	double dfvare_ = -2.0;
	double vara_ = 0.0;
	double s2vara_ = 0.0;
	double s2varg_ = 0.0;
	double yy_ = 0.0;
	double lambda2_ = 0.0;
	double shape0_ = 1.1;
	double rate0_ = 0.0;
	arma::vec vargL_;

	// SBayesD()'s defaults. Its genotype variances vx are R's diagonal,
	// (X'X)_jj / w_j, and its y'y, the mean of (X'X)_jj (b^2 + (n - 2) se^2),
	// is the mean of (X'X)_jj for standardized statistics.
	void init_priors(double pi, double h2) {
		double sumvx = 0.0;
		yy_ = 0.0;
		for (arma::uword j = 0; j < p_; ++j) {
			yy_ += ld_.diag(j);
			sumvx += ld_.diag(j) / data_.w(j);
		}
		yy_ /= static_cast<double>(p_);
		double vary = yy_ / (data_.n - 1.0);
		vara_ = ((dfvara_ - 2.0) / dfvara_) * vary * h2;
		ve_ = vary * (1.0 - h2);
		s2vara_ = vara_ * (dfvara_ - 2.0) / dfvara_;
		pi_ = (method_ == L) ? std::vector<double>{0.0, 1.0} : std::vector<double>{1.0 - pi, pi};
		double scale = (1.0 - pi_[0]) * sumvx;
		s2varg_ = s2vara_ / scale;
		double R2 = (dfvara_ - 2.0) / dfvara_;
		lambda2_ = 2.0 * (1.0 - R2) / R2 * sumvx;
		rate0_ = (shape0_ - 1.0) / lambda2_;
		vargL_ = arma::vec(p_, arma::fill::value(vara_ / scale));
	}

	void sweep_spike_slab(arma::uword j) {
		double xx = this->xx(j);
		double rhs_j = rhs(j);
		double gi = b_(j);
		double varg = (gi * gi + s2varg_ * dfvara_) / chisq(dfvara_ + 1.0);
		double logdetV = std::log(varg * xx / ve_ + 1.0);
		double uhat = rhs_j / (xx + ve_ / varg);
		double s1 = -0.5 * (logdetV - rhs_j * uhat / ve_) + std::log(pi_[1]);
		double accept = 1.0 / (1.0 + std::exp(s1 - std::log(pi_[0])));
		double v = xx + ve_ / varg;
		cond_mean_(j) = (1.0 - accept) * rhs_j / v;
		d_[j] = unif() < accept ? 0 : 1;
		if (d_[j] == 0) {
			set(j, 0.0);
			return;
		}
		if (set(j, normal(rhs_j / v, std::sqrt(ve_ / v))) == 0.0) {
			d_[j] = 0;
		}
	}

	// Park and Casella: b_j's prior variance is ve * tau_j^2, and 1 / tau_j^2
	// has an inverse-Gaussian conditional.
	void sweep_lasso(arma::uword j) {
		double xx = this->xx(j);
		double rhs_j = rhs(j);
		double v = xx + 1.0 / vargL_(j);
		cond_mean_(j) = rhs_j / v;
		double gi = normal(rhs_j / v, std::sqrt(ve_ / v));
		if (std::abs(gi) < 1e-6) {
			gi = 1e-6;
		}
		double lambda = std::sqrt(lambda2_);
		double vargi = 1.0 / rinvgauss(std::sqrt(ve_) * lambda / std::abs(gi), lambda2_, engine());
		if (vargi > 0.0) {
			vargL_(j) = vargi;
		}
		set(j, gi);
	}

	void update_variances() {
		arma::vec r = residuals();
		vara_ = (arma::dot(b_, data_.wy - r) + s2vara_ * dfvara_) / chisq(data_.n + dfvara_);
		ve_ = (yy_ - arma::dot(b_, data_.wy + r)) / chisq(data_.n + dfvare_);
		if (ve_ < 0.0) {
			ve_ = vara_ * 0.5;
		}
	}
};

// ==== LDpred2 ===================================================================

// The point-normal prior of LDpred2's Gibbs samplers in bigsnpr, read from a
// factor F (k x p) of R = F'F / (k - 1) through rss_ld.h's FactorLd, so a
// single-site update costs O(k) and nothing p x p is formed:
//   * LDpred2-grid holds the proportion p of causal variants and the
//     heritability h2 fixed; with `sparse`, an effect whose posterior
//     inclusion probability falls below p is set to zero;
//   * LDpred2-auto draws p from its Beta full conditional, with the counts
//     divided by the mean LD score, and takes h2 = b'R_s b after each sweep,
//     for R_s = c R + (1 - c) I at the LD shrinkage c.
// Either reads R_s, and with no sign jumps sets a drawn effect whose sign
// differs from the current one to zero; bigsnpr's grid has neither, which
// is c = 1 with sign jumps allowed. Both report what bigsnpr does:
// Rao-Blackwellized effects, the average over the kept sweeps of each
// variant's conditional posterior mean, with the inclusion probabilities and
// the fitted correlations R_s b. A chain whose squared effects sum past
// twice the squared marginal effects has diverged and stops.

// How a chain runs: LDpred2-auto's starting p and h2, or LDpred2-grid's fixed
// ones.
struct Ldpred2Settings {
	bool automatic = false;
	double p = 0.1;
	double h2 = 0.1;
	bool sparse = false;  // grid only
	double shrink = 1.0;  // the LD shrinkage c
	bool allow_jump_sign = true;
	double p_lower = 1e-5;  // auto only: the bounds on p
	double p_upper = 1.0;
	double mean_ld = 1.0;  // auto only: the mean LD score
};

// The chain's report besides the effects every chain reports: those are
// Rao-Blackwellized and the fitted correlations are R_s b, as bigsnpr has
// them.
struct Ldpred2Estimates : ChainEstimates {
	arma::vec postp;  // inclusion probabilities
	std::vector<double> path_p;  // auto: p after each sweep
	std::vector<double> path_h2;  // auto: h2 after each sweep
};

// The LD is FactorLd's R_s, its shrinkage s = 1 - c, so (R_s b)_j is its
// row(); the chain always Rao-Blackwellizes and tracks, as bigsnpr does.
class Ldpred2Chain : private ChainRandom {
public:
	Ldpred2Chain(const arma::mat& F, const arma::vec& bhat, const arma::vec& n,
	             const Ldpred2Settings& settings, unsigned int seed)
	    : ChainRandom(seed), bhat_(bhat), n_(n), s_(settings), m_(bhat.n_elem),
	      b_(m_, arma::fill::zeros), ld_(F, 1.0 - settings.shrink, b_),
	      cond_mean_(m_, arma::fill::zeros), cond_postp_(m_, arma::fill::zeros),
	      cond_corr_(m_, arma::fill::zeros), gap0_(divergence_gap(bhat)) {
		opts_.rao_blackwell = true;
		opts_.track = true;
		est_.beta.zeros(m_);
		est_.postp.zeros(m_);
		est_.corr.zeros(m_);
		if (s_.automatic) {
			s_.h2 = std::max(s_.h2, kMinH2);
			s_.p = std::min(std::max(s_.p, s_.p_lower), s_.p_upper);
		}
	}

	// One pass over the variants; LDpred2-auto then draws p and sets h2.
	void sweep() {
		const double inv_odd_p = (1.0 - s_.p) / s_.p;
		const double sigma2 = s_.h2 / (static_cast<double>(m_) * s_.p);
		arma::uword causal = 0;
		for (arma::uword j = 0; j < m_; ++j) {
			if (update(j, inv_odd_p, sigma2)) {
				++causal;
			}
		}
		if (s_.automatic) {
			draw_parameters(causal);
		}
	}

	// A kept sweep adds its conditional means to the posterior means. Every
	// effect was just redrawn, so b'b is bigsnpr's gap.
	void keep(int n_kept) {
		if (check_divergence(est_, opts_, b_, gap0_)) {
			return;
		}
		keep_effects(est_, opts_, b_, cond_mean_, n_kept);
		est_.postp += cond_postp_ / n_kept;
		est_.corr += cond_corr_ / n_kept;
	}

	void report(int itr, int n_iter) const {
		Rprintf("Iteration %4d of %d: p %g, h2 %g\n", itr, n_iter, s_.p, s_.h2);
	}

	bool stopped() const { return est_.diverged; }

	const Ldpred2Estimates& estimates() const { return est_; }

private:
	static constexpr double kMinH2 = 1e-3;

	// Variant j's update; true when it draws a non-zero effect. The fitted
	// correlation (R_s b)_j net of b_j's own term is the other effects' share
	// of the statistic.
	bool update(arma::uword j, double inv_odd_p, double sigma2) {
		const double fitted = ld_.row(j, b_);
		const double prev = b_(j);
		const double res = bhat_(j) - (fitted - ld_.diag(j) * prev);
		const double C1 = sigma2 * n_(j);
		const double C2 = 1.0 / (1.0 + 1.0 / C1);
		const double C3 = C2 * res;
		const double C4 = C2 / n_(j);
		const double postp =
		    1.0 / (1.0 + inv_odd_p * std::sqrt(1.0 + C1) * std::exp(-C3 * C3 / C4 / 2.0));
		cond_mean_(j) = C3 * postp;
		cond_postp_(j) = postp;
		cond_corr_(j) = fitted;
		double next = 0.0;
		bool drawn = false;
		if (s_.sparse && postp < s_.p) {
			cond_mean_(j) = 0.0;
		} else if (postp > unif()) {
			const double draw = normal(C3, std::sqrt(C4));
			if (s_.allow_jump_sign || draw * prev >= 0.0) {
				next = draw;
				drawn = true;
			}
		}
		const double diff = next - prev;
		if (diff != 0.0) {
			cur_h2_ += diff * (2.0 * fitted + diff * ld_.diag(j));
			ld_.update(j, diff);
			b_(j) = next;
		}
		return drawn;
	}

	// LDpred2-auto's p from its Beta full conditional, within its bounds, and
	// h2 = b'R_s b, floored as bigsnpr floors it.
	void draw_parameters(arma::uword causal) {
		const double k = static_cast<double>(causal);
		const double m = static_cast<double>(m_);
		beta_distribution beta(1.0 + k / s_.mean_ld, 1.0 + (m - k) / s_.mean_ld);
		s_.p = std::min(std::max(beta(engine()), s_.p_lower), s_.p_upper);
		s_.h2 = std::max(cur_h2_, kMinH2);
		est_.path_p.push_back(s_.p);
		est_.path_h2.push_back(s_.h2);
	}

	const arma::vec& bhat_;
	const arma::vec& n_;
	Ldpred2Settings s_;
	arma::uword m_;
	arma::vec b_;
	FactorLd ld_;
	ChainOptions opts_;
	arma::vec cond_mean_;
	arma::vec cond_postp_;
	arma::vec cond_corr_;
	double gap0_;
	double cur_h2_ = 0.0;
	Ldpred2Estimates est_;
};
