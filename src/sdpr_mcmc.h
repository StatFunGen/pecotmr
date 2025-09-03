#include <iostream>
#include <vector>
#include <current/armadillo>
#include <cmath>
#include <random>
#include <unordered_map>

typedef struct {
	std::vector<arma::mat> A;
	std::vector<arma::mat> B;
	std::vector<arma::mat> L;
	std::vector<arma::vec> beta_mrg;
	std::vector<arma::vec> calc_b_tmp;
	std::vector<double> num;
	std::vector<double> denom;
} ldmat_data;

class beta_distribution {
public:
// Constructor that takes alpha and beta parameters
beta_distribution(double alpha, double beta)
	: a_distribution(alpha, 1.0), b_distribution(beta, 1.0) {
}

// Generates a beta-distributed random number
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

/**
 * @class mcmc_data
 * @brief Represents the data used in the MCMC algorithm.
 *
 * The mcmc_data class holds the necessary data for the MCMC algorithm, including
 * the marginal beta values, LD matrix, sample sizes, and genotyping array information.
 *
 * @details
 * The mcmc_data class is designed to store and manage the data required for the MCMC
 * algorithm. It includes the following member variables:
 * - beta_mrg: A vector of marginal beta values for each SNP.
 * - ref_ld_mat: A vector of LD matrices, where each matrix corresponds to a subset of SNPs.
 * - sz: A vector of sample sizes for each SNP.
 * - array: A vector of genotyping array information for each SNP.
 *
 * The class also computes the boundary indices for each LD matrix, which represent the
 * start and end positions of the SNPs in the full set of SNPs.
 *
 * @param beta_mrg A vector of marginal beta values for each SNP.
 * @param ref_ld_mat A vector of LD matrices, where each matrix corresponds to a subset of SNPs.
 * @param sz A vector of sample sizes for each SNP.
 * @param array A vector of genotyping array information for each SNP.
 *
 * @note The `array` parameter is used to store information about the genotyping array for each SNP.
 * It is typically an integer value that represents the specific array used for genotyping.
 * For example, 1 may represent the Affymetrix array, while 2 may represent the Illumina array.
 *
 * @note The `sz` parameter represents the sample size for each SNP. It is used to scale the
 * elements of the LD matrix based on the sample sizes of the SNPs.
 *
 * Example usage:
 * @code
 * std::vector<double> beta_mrg = {0.1, 0.2, 0.3, 0.4, 0.5};
 * std::vector<arma::mat> ref_ld_mat = {arma::mat(2, 2, arma::fill::eye), arma::mat(3, 3, arma::fill::eye)};
 * std::vector<double> sz = {1000, 1000, 1500, 1500, 1500};
 * std::vector<int> array = {1, 1, 2, 2, 2};
 *
 * mcmc_data data(beta_mrg, ref_ld_mat, sz, array);
 * @endcode
 */
class mcmc_data {
public:
std::vector<double> beta_mrg;
std::vector<std::pair<size_t, size_t> > boundary;
std::vector<arma::mat> ref_ld_mat;
std::vector<double> sz;
std::vector<int> array;

/**
 * @brief Constructs an mcmc_data object.
 *
 * @param beta_mrg A vector of marginal beta values for each SNP.
 * @param ref_ld_mat A vector of LD matrices, where each matrix corresponds to a subset of SNPs.
 * @param sz A vector of sample sizes for each SNP.
 * @param array A vector of genotyping array information for each SNP.
 */
mcmc_data(const std::vector<double>&    beta_mrg,
          const std::vector<arma::mat>& ref_ld_mat,
          const std::vector<double>&    sz,
          const std::vector<int>&       array)
	: beta_mrg(beta_mrg), ref_ld_mat(ref_ld_mat), sz(sz), array(array) {
	compute_boundary();
}

private:
void compute_boundary() {
	size_t start = 0;
	for (const auto& mat : ref_ld_mat) {
		size_t end = start + mat.n_rows;
		boundary.emplace_back(start, end);
		start = end;
	}
}
};

class MCMC_state {
public:
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
MCMC_state(size_t num_snp, size_t max_cluster, \
           double a0, double b0, double sz) {
	a0k = a0; b0k = b0; N = sz;
	// Changed May 20 2021
	// Now N (sz) is absorbed into A, B; so set to 1.
	N = 1.0;

	n_snp = num_snp;
	M = max_cluster;
	alpha = 1;
	eta = 1;
	beta = arma::vec(num_snp, arma::fill::zeros);
	b = arma::vec(num_snp, arma::fill::zeros);
	p.assign(max_cluster, 1.0/max_cluster);
	log_p.assign(max_cluster, 0);
	for (size_t i=0; i<max_cluster; i++) {
		log_p[i] = logf(p[i] + 1e-40);
	}

	cluster_var.assign(max_cluster, 0.0);
	suff_stats.assign(max_cluster, 0);
	sumsq.assign(max_cluster, 0.0);
	V.assign(max_cluster, 0.0);
	cls_assgn.assign(num_snp, 0);
	std::random_device rd;
	r.seed(rd());
	std::uniform_int_distribution<int> dist(0, M-1);
	for (size_t i=0; i<num_snp; i++) {
		cls_assgn[i] = dist(r);
	}
}

void sample_sigma2();
void calc_b(size_t j, const mcmc_data &dat, const ldmat_data &ldmat_dat);
void sample_assignment(size_t j, const mcmc_data &dat, \
                       const ldmat_data &ldmat_dat);
void update_suffstats();
void sample_V();
void update_p();
void sample_alpha();
void sample_beta(size_t j, const mcmc_data &dat, \
                 ldmat_data &ldmat_dat);
void compute_h2(const mcmc_data &dat);
void sample_eta(const ldmat_data &ldmat_dat);

private:
double a0k;
double b0k;
size_t M, n_snp;
std::mt19937 r;
};

class MCMC_samples {
public:
arma::vec beta;
double h2;
MCMC_samples(size_t num_snps) {
	beta = arma::vec(num_snps, arma::fill::zeros);
	h2 = 0;
}
};

/**
 * @brief Perform Markov Chain Monte Carlo (MCMC) for estimating effect sizes.
 *
 * @param data The `mcmc_data` object containing the necessary data for the MCMC algorithm.
 * @param sz The sample size of the GWAS.
 * @param a Factor to shrink the reference LD matrix. Default is 0.1.
 * @param c Factor to correct for the deflation. Default is 1.
 * @param M Max number of variance components. Default is 1000.
 * @param a0k Hyperparameter for inverse gamma distribution. Default is 0.5.
 * @param b0k Hyperparameter for inverse gamma distribution. Default is 0.5.
 * @param iter Number of iterations for MCMC. Default is 1000.
 * @param burn Number of burn-in iterations for MCMC. Default is 200.
 * @param thin Thinning interval for MCMC. Default is 5.
 * @param n_threads Number of threads to use. Default is 1.
 * @param opt_llk Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs)
 *                and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort).
 *                Default is 1.
 * @param verbose Whether to print verbose output. Default is true.
 *
 * @return An `std::unordered_map` containing the estimated effect sizes (beta) and heritability (h2).
 *
 * @details
 * This function performs Markov Chain Monte Carlo (MCMC) to estimate effect sizes and heritability
 * based on the provided `mcmc_data` object and other parameters.
 *
 * The `mcmc_data` object should contain the following:
 * - `beta_mrg`: A vector of marginal beta values for each SNP.
 * - `ref_ld_mat`: A vector of LD matrices, where each matrix corresponds to a subset of SNPs.
 * - `sz`: A vector of sample sizes for each SNP.
 * - `array`: A vector of genotyping array information for each SNP.
 *
 * The function returns an `std::unordered_map` with the following key-value pairs:
 * - "beta": An `arma::vec` containing the estimated effect sizes for each SNP.
 * - "h2": An `arma::vec` containing the estimated heritability.
 *
 * Additional parameters:
 * - `sz`: The sample size of the GWAS.
 * - `a`: Factor to shrink the reference LD matrix. Default is 0.1.
 * - `c`: Factor to correct for the deflation. Default is 1.
 * - `M`: Max number of variance components. Default is 1000.
 * - `a0k`, `b0k`: Hyperparameters for the inverse gamma distribution. Default is 0.5 for both.
 * - `iter`: Number of iterations for MCMC. Default is 1000.
 * - `burn`: Number of burn-in iterations for MCMC. Default is 200.
 * - `thin`: Thinning interval for MCMC. Default is 5.
 * - `n_threads`: Number of threads to use. Default is 1.
 * - `opt_llk`: Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs)
 *              and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort).
 *              Default is 1.
 * - `verbose`: Whether to print verbose output. Default is true.
 *
 * @note The `mcmc` function assumes the existence of the `Function_pool` class and its member functions,
 *       as well as the `mcmc.h` header file with the necessary class and struct definitions.
 */
std::unordered_map<std::string, arma::vec> mcmc(
	mcmc_data& data,
	unsigned         sz,
	double           a,
	double           c,
	size_t           M,
	double           a0k,
	double           b0k,
	int              iter,
	int              burn,
	int              thin,
	unsigned         n_threads,
	int              opt_llk,
	bool             verbose
	);
