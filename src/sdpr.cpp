#include <RcppArmadillo.h>
#include <random>
#include <unordered_map>
#include "sdpr_mcmc.h"

// Rcpp interface function
// [[Rcpp::export]]
Rcpp::List sdpr_rcpp(
	const std::vector<double>&          bhat,
	const Rcpp::List&                   LD,
	int                                 n,
	Rcpp::Nullable<Rcpp::NumericVector> per_variant_sample_size = R_NilValue,
	Rcpp::Nullable<Rcpp::IntegerVector> array = R_NilValue,
	double                              a = 0.1,
	double                              c = 1.0,
	size_t                              M = 1000,
	double                              a0k = 0.5,
	double                              b0k = 0.5,
	int                                 iter = 1000,
	int                                 burn = 200,
	int                                 thin = 5,
	unsigned                            n_threads = 1,
	int                                 opt_llk = 1,
	bool                                verbose = true,
	Rcpp::Nullable<unsigned int>        seed = R_NilValue
	) {
	// Convert Rcpp::List to std::vector<arma::mat>
	std::vector<arma::mat> ref_ld_mat;
	for (int i = 0; i < LD.size(); i++) {
		ref_ld_mat.push_back(Rcpp::as<arma::mat>(LD[i]));
	}

	// Initialize per_variant_sample_size and array if NULL
	std::vector<double> sz;
	std::vector<int> arr;
	if (per_variant_sample_size.isNotNull()) {
		sz = Rcpp::as<std::vector<double> >(per_variant_sample_size);
	} else {
		sz = std::vector<double>(bhat.size(), n);
	}
	if (array.isNotNull()) {
		arr = Rcpp::as<std::vector<int> >(array);
	} else {
		arr = std::vector<int>(bhat.size(), 1);
	}

	// Resolve seed
	unsigned int seed_val = 0;
	if (seed.isNotNull()) {
		seed_val = Rcpp::as<unsigned int>(seed);
	} else {
		seed_val = std::random_device{}();
	}

	// Create mcmc_data object
	mcmc_data data(bhat, ref_ld_mat, sz, arr);

	// Call the mcmc function
	std::unordered_map<std::string, arma::vec> results = mcmc(
		data, n, a, c, M, a0k, b0k, iter, burn, thin, n_threads, opt_llk, verbose, seed_val
		);

	// Convert results to Rcpp::List
	Rcpp::List output = Rcpp::List::create(
		Rcpp::Named("beta_est") = results["beta"],
		Rcpp::Named("h2") = results["h2"]
		);

	return output;
}
