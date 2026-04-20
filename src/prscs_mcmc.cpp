/**
 * @file prs_cs_wrapper.cpp
 * @brief cpp11 wrapper for the prs_cs function.
 */

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include "prscs_mcmc.h"

using namespace cpp11;
using namespace arma;

/**
 * @brief cpp11 wrapper for the prs_cs function.
 *
 * @param a Shape parameter for the prior distribution of psi.
 * @param b Scale parameter for the prior distribution of psi.
 * @param phi Global shrinkage parameter. If NULL, it will be estimated automatically.
 * @param bhat Vector of effect sizes.
 * @param maf Vector of minor allele frequencies. If NULL, it is assumed to be a vector of zeros.
 * @param n Sample size.
 * @param ld_blk List of LD blocks.
 * @param n_iter Number of MCMC iterations.
 * @param n_burnin Number of burn-in iterations.
 * @param thin Thinning interval.
 * @param verbose Whether to print verbose output.
 * @param seed Random seed. If NULL, no seed is set.
 * @return A list containing the posterior estimates.
 */
[[cpp11::register]]
cpp11::writable::list prs_cs_rcpp(double a, double b, sexp phi,
                       doubles bhat, sexp maf,
                       int n, list ld_blk,
                       int n_iter, int n_burnin, int thin,
                       bool verbose, sexp seed) {
	// Convert cpp11 types to C++ types
	std::vector<double> bhat_vec = cpp11::as_cpp<std::vector<double>>(bhat);
	std::vector<double> maf_vec;
	if (maf != R_NilValue) {
		maf_vec = cpp11::as_cpp<std::vector<double>>(maf);
	} else {
		maf_vec = std::vector<double>(bhat_vec.size(), 0.0); // Populate with zeros if maf is NULL
	}

	std::vector<mat> ld_blk_vec;
	for (int i = 0; i < ld_blk.size(); ++i) {
		ld_blk_vec.push_back(as_Mat(doubles_matrix<>(ld_blk[i])));
	}

	// Use stack variable to avoid heap allocation and memory leak risk.
	double phi_val = 0.0;
	double* phi_ptr = nullptr;
	if (phi != R_NilValue) {
		phi_val = cpp11::as_cpp<double>(phi);
		phi_ptr = &phi_val;
	}

	unsigned int seed_val = 0;
	if (seed != R_NilValue) {
		seed_val = cpp11::as_cpp<unsigned int>(seed);
	} else {
		seed_val = std::random_device{}();
	}

	std::map<std::string, vec> output = prs_cs_mcmc(a, b, phi_ptr, bhat_vec, maf_vec, n, ld_blk_vec,
	                                                      n_iter, n_burnin, thin, verbose, seed_val);

	// Convert the output to a list
	using namespace cpp11::literals;
	writable::list result({
		"beta_est"_nm = as_doubles(output["beta_est"]),
		"psi_est"_nm = as_doubles(output["psi_est"]),
		"sigma_est"_nm = cpp11::as_sexp(output["sigma_est"](0)),
		"phi_est"_nm = cpp11::as_sexp(output["phi_est"](0))
	});

	return result;
}
