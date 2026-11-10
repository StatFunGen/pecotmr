// rss_ld.h -- the LD a coordinate-descent RSS solve reads.
//
// The solvers in penalized_rss.cpp minimise
// beta'R_s beta - 2 beta'z + penalty(beta) for R_s = (1 - s) R + s I. They
// only ever ask R_s four things: a diagonal entry, one entry of R_s beta for
// the current beta, the effect of changing one coordinate, and (for the L0
// swaps) a whole column or product. Each class below answers those from one
// way of holding R_s, so a solver templated on it runs unchanged on either.
#pragma once

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <vector>

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

// ---- Studies' variants in their union ----------------------------------------

// The solvers across studies -- penalized_rss.cpp's coupled path and the
// PRS-CS chain of rss_mcmc.h -- match a variant across studies by its place in
// the union of the studies' variants. Each study's places as R passes them,
// 1-based, made 0-based: one per value of the study's (`what` names them),
// each within the union and none twice in a study.
inline std::vector<arma::uvec> union_index(const cpp11::list& idxList,
                                           const std::vector<arma::vec>& values,
                                           int nUnion, const char* what) {
	if (idxList.size() != static_cast<R_xlen_t>(values.size())) {
		cpp11::stop("There are %d union indices for %d studies.", (int)idxList.size(),
		            (int)values.size());
	}
	std::vector<arma::uvec> out;
	for (arma::uword k = 0; k < values.size(); ++k) {
		cpp11::integers ix(idxList[static_cast<R_xlen_t>(k)]);
		if (static_cast<arma::uword>(ix.size()) != values[k].n_elem) {
			cpp11::stop("Study %d has %d union indices but %d %s.", (int)k + 1,
			            (int)ix.size(), (int)values[k].n_elem, what);
		}
		std::vector<unsigned char> seen(nUnion, 0);
		arma::uvec u(ix.size());
		for (R_xlen_t i = 0; i < ix.size(); ++i) {
			int v = ix[i];
			if (v < 1 || v > nUnion || seen[v - 1]) {
				cpp11::stop("Study %d's union indices must be distinct, in 1..%d.",
				            (int)k + 1, nUnion);
			}
			seen[v - 1] = 1;
			u(i) = static_cast<arma::uword>(v - 1);
		}
		out.push_back(u);
	}
	return out;
}
