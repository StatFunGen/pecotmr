// rss_ld.h -- the LD a coordinate-descent RSS solve reads.
//
// The solvers in penalized_rss.cpp minimise
// beta'R_s beta - 2 beta'z + penalty(beta) for R_s = (1 - s) R + s I. They
// only ever ask R_s four things: a diagonal entry, one entry of R_s beta for
// the current beta, the effect of changing one coordinate, and (for the L0
// swaps) a whole column or product. Each class below answers those from one
// way of holding R_s, so a solver templated on it runs unchanged on either.
#pragma once

#include <cpp11armadillo.hpp>

// R held whole, with R_s = (1 - s) R + s I applied on the fly and Rbeta =
// R beta kept as a running sum (the genotype-free analogue of tracking
// yhat = X beta), so a shrinkage grid reuses the one matrix rather than
// forming R_s for each s. With s = 0 every answer is R's own, bit for bit.
class DenseLd {
public:
	DenseLd(const arma::mat& R, double s, const arma::vec& beta)
	    : R_(R), s_(s), c_(1.0 - s), Rbeta_(R * beta) {}

	double diag(arma::uword j) const { return c_ * R_(j, j) + s_; }

	// (R_s beta)_j for the current beta.
	double row(arma::uword j, const arma::vec& beta) const {
		return c_ * Rbeta_(j) + s_ * beta(j);
	}

	// Coordinate j of beta changed by `del`.
	void update(arma::uword j, double del) { Rbeta_ += del * R_.col(j); }

	// R_s beta for the current beta, and column j of R_s.
	arma::vec product(const arma::vec& beta) const { return c_ * Rbeta_ + s_ * beta; }
	arma::vec column(arma::uword j) const {
		arma::vec col = c_ * R_.col(j);
		col(j) += s_;
		return col;
	}

	// beta'R_s beta, recomputed from R.
	double quad(const arma::vec& beta) const {
		return c_ * arma::as_scalar(beta.t() * R_ * beta) + s_ * arma::dot(beta, beta);
	}

private:
	const arma::mat& R_;
	double s_;
	double c_;
	arma::vec Rbeta_;
};

// R = F'F / (m - 1) held as its factor F (m x p), with R_s = (1 - s) R + s I.
// It keeps u = F beta instead of R_s beta, so with c = (1 - s) / (m - 1),
// (R_s beta)_j = c f_j'u + s beta_j: a coordinate step costs O(m) rather than
// O(p), and neither R nor R_s is ever formed.
class FactorLd {
public:
	FactorLd(const arma::mat& F, double s, const arma::vec& beta)
	    : F_(F), s_(s), c_((1.0 - s) / (static_cast<double>(F.n_rows) - 1.0)), u_(F * beta),
	      diag_(c_ * arma::sum(arma::square(F), 0).t() + s) {}

	double diag(arma::uword j) const { return diag_(j); }

	double row(arma::uword j, const arma::vec& beta) const {
		return c_ * arma::dot(F_.col(j), u_) + s_ * beta(j);
	}

	void update(arma::uword j, double del) { u_ += del * F_.col(j); }

	arma::vec product(const arma::vec& beta) const { return c_ * (F_.t() * u_) + s_ * beta; }

	arma::vec column(arma::uword j) const {
		arma::vec col = c_ * (F_.t() * F_.col(j));
		col(j) += s_;
		return col;
	}

	double quad(const arma::vec& beta) const {
		arma::vec Fb = F_ * beta;
		return c_ * arma::dot(Fb, Fb) + s_ * arma::dot(beta, beta);
	}

private:
	const arma::mat& F_;
	double s_;
	double c_;
	arma::vec u_;
	arma::vec diag_;
};
