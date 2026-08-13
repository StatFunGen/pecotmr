#ifndef QTL_ENRICHMENT_HPP
#define QTL_ENRICHMENT_HPP
#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <random>
#include <omp.h>
#include <cmath>
#include <cstdio>

using namespace cpp11;
using namespace arma;

class SuSiEFit {
public:
std::vector<std::string> variable_names;
mat alpha;
std::vector<double> prior_variance;

SuSiEFit(SEXP r_susie_fit) {
	list susie_fit(r_susie_fit);

	doubles pip_vec(susie_fit["pip"]);
	// Get names from the pip vector
	cpp11::strings pip_names(pip_vec.attr("names"));
	variable_names.reserve(pip_names.size());
	for (int i = 0; i < pip_names.size(); ++i) {
		variable_names.push_back(std::string(pip_names[i]));
	}

	alpha = as_Mat(doubles_matrix<>(susie_fit["alpha"]));
	prior_variance = cpp11::as_cpp<std::vector<double>>(susie_fit["prior_variance"]);

	if (alpha.n_rows != prior_variance.size()) {
		cpp11::stop("The number of rows in alpha must match the length of prior_variance.");
	}

	// Check if all elements in prior_variance are not greater than 0
	if (std::all_of(prior_variance.begin(), prior_variance.end(), [](double x) {
			return x <= 0;
		})) {
		cpp11::stop("At least one element in prior_variance must be greater than 0.");
	}

	// Filter out rows with prior_variance = 0
	std::vector<uword> valid_rows;
	for (size_t i = 0; i < prior_variance.size(); ++i) {
		if (prior_variance[i] > 0) {
			valid_rows.push_back(i);
		}
	}
	alpha = alpha.rows(uvec(valid_rows));
	prior_variance.erase(std::remove(prior_variance.begin(), prior_variance.end(), 0), prior_variance.end());

	// Add a check to make sure each row of alpha sums to 1
	for (uword i = 0; i < alpha.n_rows; ++i) {
		double row_sum = sum(alpha.row(i));
		if (std::abs(row_sum - 1.0) > 1e-6) {
			cpp11::stop("Row %d of single effect PIP matrix (alpha) does not sum to 1. It is: %g",
			            (int)(i + 1), row_sum);
		}
	}
}

std::vector<std::string> impute_qtn(std::mt19937 &gen) const {
	std::vector<std::string> qtn_names;

	for (uword i = 0; i < alpha.n_rows; ++i) {
		arma::rowvec row_i = alpha.row(i);
		std::vector<double> alpha_row(row_i.begin(), row_i.end());
		std::discrete_distribution<> dist(alpha_row.begin(), alpha_row.end());
		int random_index = dist(gen);
		qtn_names.push_back(variable_names[random_index]);
	}

	return qtn_names;
}
};

std::vector<size_t> filter_outlier_indices(
    const std::vector<double>& estimates,
    const std::vector<double>& variances,
    double prior_variance,
    bool bessel_correction = true,
    double threshold = 3.0)
{
    size_t n = estimates.size();
    double mean = 0.0;
    double sd = 0.0;
    std::vector<double> shrinkage_ests;

    // Calculate shrinkage estimates and mean
    for(size_t i = 0; i < n; ++i) {
        double shrinkage = (estimates[i] * prior_variance) /
                          (prior_variance + variances[i]);
        shrinkage_ests.push_back(shrinkage);
        mean += shrinkage;
    }
    mean /= n;

    // Calculate standard deviation
    for(size_t i = 0; i < n; ++i) {
        sd += pow(shrinkage_ests[i] - mean, 2);
    }
    double denom = bessel_correction ? (n - 1) : n;
    sd = sqrt(sd / denom);

    // Return indices of non-outlier elements
    std::vector<size_t> kept;
    for(size_t i = 0; i < n; ++i) {
        if(fabs(shrinkage_ests[i] - mean) <= threshold * sd) {
            kept.push_back(i);
        }
    }

    return kept;
}

std::vector<double> run_EM(
	const std::vector<double> &gwas_pip,
	const std::vector<int> &   annotation_vector,
	double                     pi_gwas,
	double                     pi_qtl,
	double                     total_snp,
	int                        max_iter = 1000,
	double                     a1_tol = 0.01)
{
	double a0 = log(pi_gwas / (1 - pi_gwas));
	double a1 = 0;
	double var0 = 0;
	double var1 = 0;
	double r0, r1;
	r0 = r1 = exp(a0);
	double r_null = pi_gwas / (1 - pi_gwas);
	int iter = 0;

	while (true) {
		iter++;
		// E step
		double pseudo_count = 1.0;
		double e0g0 = pseudo_count * (1 - pi_gwas) * (1 - pi_qtl);
		double e0g1 = pseudo_count * (1 - pi_qtl) * pi_gwas;
		double e1g0 = pseudo_count * (1 - pi_gwas) * pi_qtl;
		double e1g1 = pseudo_count * pi_gwas * pi_qtl;

		for (size_t i = 0; i < gwas_pip.size(); i++) {
			double val = gwas_pip[i];
			if (val == 1)
				val = 1 - 1e-8;
			// posterior ratio
			val = val / (1 - val);
			// val/r_null is marginal likelihood/bayes factor
			if (annotation_vector[i] == 0) {
				val = r0 * (val / r_null);
				// updated posterior with current prior given eqtl = 0
				val = val / (1 + val);
				e0g1 += val;
				e0g0 += 1 - val;
			}

			if (annotation_vector[i] == 1) {
				val = r1 * (val / r_null);
				// updated posterior with current prior given eqtl = 1
				val = val / (1 + val);
				e1g1 += val;
				e1g0 += 1 - val;
			}
		}

		e0g0 += total_snp - (e0g0 + e0g1 + e1g0 + e1g1);

		double a1_new = log(e1g1 * e0g0 / (e1g0 * e0g1));

		if (fabs(a1_new - a1) < a1_tol || iter >= max_iter) {
			a1 = a1_new;
			a0 = log(e0g1 / e0g0);
			var1 = (1.0 / e0g0 + 1.0 / e1g0 + 1.0 / e1g1 + 1.0 / e0g1);
			var0 = (1.0 / e0g1 + 1.0 / e0g0);
			break;
		}

		a1 = a1_new;
		a0 = log(e0g1 / e0g0);
		r0 = exp(a0);
		r1 = exp(a0 + a1);

		if (iter % 100 == 0) {
			Rprintf("EM Iteration %d: a0 = %g, a1 = %g\n", iter, a0, a1);
		}
	}
	if (iter == max_iter) {
		Rprintf("WARNING: EM algorithm did not converge after %d iterations!\n", iter);
	}

	std::vector<double> av;
	av.push_back(a0);
	av.push_back(a1);
	av.push_back(var0);
	av.push_back(var1);

	return av;
}

std::map<std::string, double> qtl_enrichment_workhorse(
	const std::vector<SuSiEFit> &   qtl_susie_fits,
	const std::vector<double> &     gwas_pip,
	const std::vector<std::string> &gwas_variable_names,
	double                          pi_gwas,
	double                          pi_qtl,
	int                             ImpN,
	double                          shrinkage_lambda,
	bool                            double_shrinkage = false,
	bool                            bessel_correction = true,
	int                             num_threads = 4)
{

	std::vector<double> a0_vec(ImpN, 0.0);
	std::vector<double> v0_vec(ImpN, 0.0);
	std::vector<double> a1_vec(ImpN, 0.0);
	std::vector<double> v1_vec(ImpN, 0.0);

	std::map<std::string, int> gwas_variant_index;

	for (size_t i = 0; i < gwas_variable_names.size(); ++i) {
		gwas_variant_index[gwas_variable_names[i]] = i;
	}

	// pi_gwas = sum(gwas_pip) / total_snp
	double total_snp = std::accumulate(gwas_pip.begin(), gwas_pip.end(), 0.0) / pi_gwas;

	Rprintf("Fine-mapped GWAS and QTL data loaded successfully for enrichment analysis!\n");

	#pragma omp parallel for num_threads(num_threads)
	for (int k = 0; k < ImpN; k++) {
		// Initialize the RNG for this thread
		std::random_device rd;
		std::mt19937 gen(rd());

		// Use QTL to annotate GWAS variants
		std::vector<int> annotation_vector(gwas_pip.size(), 0);
		int missing_qtl_count = 0;
		int total_qtl_count = 0;

		for (size_t i = 0; i < qtl_susie_fits.size(); i++) {
			std::vector<std::string> variants = qtl_susie_fits[i].impute_qtn(gen);
			for (const auto &variant : variants) {
				auto it = gwas_variant_index.find(variant);
				if (it != gwas_variant_index.end()) {
					annotation_vector[it->second] = 1;
				} else {
					++missing_qtl_count;
				}
			}
			total_qtl_count += variants.size();
		}
		double missing_variant_proportion = static_cast<double>(missing_qtl_count) / total_qtl_count;
		std::vector<double> rst = run_EM(gwas_pip, annotation_vector, pi_gwas, pi_qtl, total_snp);

	#pragma omp critical
		{
			a0_vec[k] = rst[0];
			a1_vec[k] = rst[1];
			v0_vec[k] = rst[2];
			v1_vec[k] = rst[3];
			Rprintf("Proportion of xQTL missing from GWAS variants: %g in MI round %d\n",
			        missing_variant_proportion, k);
		}
	}

	Rprintf("EM updates completed!\n");

	// Apply outlier filtering if shrinkage is specified.
	// Use index-based filtering to keep all four vectors aligned.
	double pv = (shrinkage_lambda == 0) ? -1.0 : 1.0 / shrinkage_lambda;
	std::vector<size_t> kept_indices;
	if (pv > 0) {
		kept_indices = filter_outlier_indices(a1_vec, v1_vec, pv, bessel_correction);
		if (kept_indices.size() < static_cast<size_t>(ImpN)) {
			Rprintf("Outlier filtering removed %d MI round(s)\n",
			        ImpN - static_cast<int>(kept_indices.size()));
		}
	} else {
		kept_indices.resize(ImpN);
		std::iota(kept_indices.begin(), kept_indices.end(), 0);
	}

	size_t m = kept_indices.size();

	// MI combining using only surviving rounds
	double a0_est = 0;
	double a1_est = 0;
	double var0 = 0;
	double var1 = 0;
	for (size_t i = 0; i < m; i++) {
		size_t k = kept_indices[i];
		a0_est += a0_vec[k];
		a1_est += a1_vec[k];
		var0 += v0_vec[k];
		var1 += v1_vec[k];
	}
	a0_est /= m;
	a1_est /= m;

	double bv0 = 0;
	double bv1 = 0;
	for (size_t i = 0; i < m; i++) {
		size_t k = kept_indices[i];
		bv0 += pow(a0_vec[k] - a0_est, 2.0);
		bv1 += pow(a1_vec[k] - a1_est, 2.0);
	}
	bv0 /= (m - 1);
	bv1 /= (m - 1);
	var0 /= m;
	var1 /= m;

	double sd0 = sqrt(var0 + bv0 * (m + 1.0) / m);
	double sd1 = sqrt(var1 + bv1 * (m + 1.0) / m);

	double a1_est_ns = a1_est;
	double sd1_ns = sd1;

	// Apply shrinkage
	if (pv > 0) {
		if (double_shrinkage) {
			// Double shrinkage (matches upstream fastenloc):
			// 1. Shrink each MI estimate individually
			double a1_shrink_est = 0;
			double var1_shrink = 0;
			std::vector<double> a1_shrink_vec;
			for (size_t i = 0; i < m; i++) {
				size_t k = kept_indices[i];
				double post_a1 = (a1_vec[k] * pv) / (pv + v1_vec[k]);
				double post_var = 1.0 / (1.0 / pv + 1.0 / v1_vec[k]);
				a1_shrink_vec.push_back(post_a1);
				a1_shrink_est += post_a1;
				var1_shrink += post_var;
			}
			a1_shrink_est /= m;
			var1_shrink /= m;

			// Between-imputation variance of shrunk estimates
			double bv1_shrink = 0;
			for (size_t i = 0; i < m; i++) {
				bv1_shrink += pow(a1_shrink_vec[i] - a1_shrink_est, 2.0);
			}
			bv1_shrink /= (m - 1);
			double sd1_shrink = sqrt(var1_shrink + bv1_shrink * (m + 1.0) / m);

			// 2. Shrink the combined estimate again
			a1_est = (a1_shrink_est * pv) / (pv + sd1_shrink * sd1_shrink);
			sd1 = sqrt(1.0 / (1.0 / pv + 1.0 / (sd1_shrink * sd1_shrink)));
		} else {
			// Single shrinkage: shrink the MI-combined estimate once
			a1_est = (a1_est_ns * pv) / (pv + sd1_ns * sd1_ns);
			sd1 = sqrt(1.0 / (1.0 / pv + 1.0 / (sd1_ns * sd1_ns)));
		}
	}

	a0_est = log(pi_gwas / (1 + pi_qtl * exp(a1_est) - pi_qtl - pi_gwas));

	double pi1_ne = exp(a0_est) / (1 + exp(a0_est));
	double pi1_e = exp(a0_est + a1_est) / (1 + exp(a0_est + a1_est));

	double p1 = (1 - pi_qtl) * pi1_ne;
	double p2 = pi_qtl / (1 + exp(a0_est + a1_est));
	double p12 = pi_qtl * pi1_e;

	// Create the map to store output
	std::map<std::string, double> output_map;
	output_map["Intercept"] = a0_est;
	output_map["sd (intercept)"] = sd0;
	output_map["Enrichment (no shrinkage)"] = a1_est_ns;
	output_map["Enrichment (w/ shrinkage)"] = a1_est;
	output_map["sd (no shrinkage)"] = sd1_ns;
	output_map["sd (w/ shrinkage)"] = sd1;
	output_map["Alternative (coloc) p1"] = p1;
	output_map["Alternative (coloc) p2"] = p2;
	output_map["Alternative (coloc) p12"] = p12;
	output_map["Effective MI rounds"] = static_cast<double>(m);

	return output_map;
}

#endif // QTL_ENRICHMENT_HPP
