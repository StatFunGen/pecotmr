#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <random>
#include <unordered_map>
#include "sdpr_mcmc.h"

using namespace cpp11;
using namespace arma;

// cpp11 interface function
[[cpp11::register]]
cpp11::writable::list sdpr_rcpp(
	const doubles& bhat_r,
	const list& LD,
	int n,
	sexp per_variant_sample_size = R_NilValue,
	sexp array = R_NilValue,
	double a = 0.1,
	double c = 1.0,
	int M = 1000,
	double a0k = 0.5,
	double b0k = 0.5,
	int iter = 1000,
	int burn = 200,
	int thin = 5,
	int n_threads = 1,
	int opt_llk = 1,
	bool verbose = true,
	sexp seed = R_NilValue
	) {
	// Convert inputs to C++ types
	std::vector<double> bhat = cpp11::as_cpp<std::vector<double>>(bhat_r);

	// Convert list to std::vector<arma::mat>
	std::vector<mat> ref_ld_mat;
	for (int i = 0; i < LD.size(); i++) {
		ref_ld_mat.push_back(as_Mat(doubles_matrix<>(LD[i])));
	}

	// Initialize per_variant_sample_size and array if NULL
	std::vector<double> sz;
	std::vector<int> arr;
	if (per_variant_sample_size != R_NilValue) {
		sz = cpp11::as_cpp<std::vector<double>>(per_variant_sample_size);
	} else {
		sz = std::vector<double>(bhat.size(), n);
	}
	if (array != R_NilValue) {
		arr = cpp11::as_cpp<std::vector<int>>(array);
	} else {
		arr = std::vector<int>(bhat.size(), 1);
	}

	// Resolve seed
	unsigned int seed_val = 0;
	if (seed != R_NilValue) {
		seed_val = cpp11::as_cpp<unsigned int>(seed);
	} else {
		seed_val = std::random_device{}();
	}

	// Create mcmc_data object
	mcmc_data data(bhat, ref_ld_mat, sz, arr);

	// Call the mcmc function
	std::unordered_map<std::string, vec> results = mcmc(
		data, n, a, c, M, a0k, b0k, iter, burn, thin, n_threads, opt_llk, verbose, seed_val
		);

	// Convert results to list
	using namespace cpp11::literals;
	writable::list output({
		"beta_est"_nm = as_doubles(results["beta"]),
		"h2"_nm = as_doubles(results["h2"])
	});

	return output;
}
