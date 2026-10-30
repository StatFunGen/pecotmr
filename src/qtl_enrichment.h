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
#include <limits>
#include <numeric>
#include <algorithm>

using namespace cpp11;
using namespace arma;

// One QTL region's fit as the R wrapper hands it over: alpha (effects x
// variables), the effects' prior variances, and each variable's integer code,
// which is the 0-based index of the same variant in the GWAS PIP vector, or a
// negative id for a variant the GWAS vector does not carry. Coding the
// variables in R keeps variant-name strings out of the kernel. Only effects
// with prior variance > 0 are kept; a fit with none of them is a null fit,
// which samples nothing and draws nothing from the random stream. The R
// wrapper validates the fit first; the checks here only guard the kernel.
class SuSiEFit {
public:
std::vector<int> codes;
mat alpha;

SuSiEFit(SEXP r_susie_fit) {
	list susie_fit(r_susie_fit);

	alpha = as_Mat(doubles_matrix<>(susie_fit["alpha"]));
	std::vector<double> prior_variance = cpp11::as_cpp<std::vector<double>>(susie_fit["prior_variance"]);
	codes = cpp11::as_cpp<std::vector<int>>(susie_fit["codes"]);

	if (alpha.n_rows != prior_variance.size()) {
		cpp11::stop("The number of rows in alpha must match the length of prior_variance.");
	}
	if (alpha.n_cols != codes.size()) {
		cpp11::stop("The number of columns in alpha must match the number of variant codes.");
	}

	std::vector<uword> valid_rows;
	for (size_t i = 0; i < prior_variance.size(); ++i) {
		if (prior_variance[i] > 0) {
			valid_rows.push_back(i);
		}
	}
	if (valid_rows.empty()) {
		alpha.set_size(0, alpha.n_cols);
	} else if (valid_rows.size() < alpha.n_rows) {
		alpha = alpha.rows(uvec(valid_rows));
	}

	for (uword i = 0; i < alpha.n_rows; ++i) {
		double row_sum = sum(alpha.row(i));
		if (std::abs(row_sum - 1.0) > 1e-6) {
			cpp11::stop("Row %d of single effect PIP matrix (alpha) does not sum to 1. It is: %g",
			            (int)(i + 1), row_sum);
		}
	}
}

// One sampled variable per kept effect, returned as its code.
std::vector<int> impute_qtn(std::mt19937 &gen) const {
	std::vector<int> qtn_codes;

	for (uword i = 0; i < alpha.n_rows; ++i) {
		arma::rowvec row_i = alpha.row(i);
		std::vector<double> alpha_row(row_i.begin(), row_i.end());
		std::discrete_distribution<> dist(alpha_row.begin(), alpha_row.end());
		int random_index = dist(gen);
		qtn_codes.push_back(codes[random_index]);
	}

	return qtn_codes;
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

    // A filter that removes EVERY round is not a filter: it means the
    // dispersion estimate was degenerate rather than that every round was an
    // outlier. With a single MI round and the Bessel correction the sd is
    // 0/0, so no |x - mean| <= 3 * NaN comparison holds and every round is
    // dropped -- which then divides by zero throughout the MI combination and
    // reports NaN for every estimate. Keep them all instead.
    if(kept.empty()) {
        kept.resize(n);
        std::iota(kept.begin(), kept.end(), 0);
    }

    return kept;
}

// One imputation round's EM estimates. Runs on OpenMP worker threads, so it
// neither prints nor raises; whether it converged, and whether the variant
// total left the (not QTL, not causal) cell positive, are returned instead.
struct EMResult {
	double a0;
	double a1;
	double var0;
	double var1;
	bool   converged;
	bool   valid;
	int    iterations;
};

// The GWAS variants the EM reads, prepared once for every round: the
// variants of non-zero PIP, each with its prior odds term
// (pip / (1 - pip)) / r_null, computed exactly as the EM loop used to compute
// it; and, for every GWAS variant, its place among them (-1 for a PIP of 0).
// A PIP of 0 has posterior 0 under either prior, so such a variant adds
// exactly one count to its annotation's (not causal) cell whatever the
// priors: the EM adds those counts once per iteration instead of visiting
// each of them.
struct GwasOdds {
	std::vector<double> odds;
	std::vector<int>    place;
	double              n_zero;
};

GwasOdds prepare_gwas_odds(const std::vector<double> &gwas_pip, double pi_gwas)
{
	double r_null = pi_gwas / (1 - pi_gwas);
	GwasOdds g;
	g.place.assign(gwas_pip.size(), -1);
	for (size_t i = 0; i < gwas_pip.size(); i++) {
		double val = gwas_pip[i];
		if (val == 0) {
			continue;
		}
		if (val == 1)
			val = 1 - 1e-8;
		// posterior ratio over the null prior odds: the marginal likelihood
		val = val / (1 - val);
		g.place[i] = static_cast<int>(g.odds.size());
		g.odds.push_back(val / r_null);
	}
	g.n_zero = static_cast<double>(gwas_pip.size() - g.odds.size());
	return g;
}

// `annotated` flags which of the non-zero-PIP variants (GwasOdds::odds) this
// round's QTL sample annotates, and `zero_annotated` how many PIP-0 variants
// it annotates. `qtl_only_sampled` counts this round's distinct sampled QTL
// variants that the GWAS PIP vector does not carry. As in fastenloc, where
// they are part of the SNP universe with GWAS PIP 0, each is a (QTL, not
// causal) variant.
EMResult run_EM(
	const GwasOdds &           gwas,
	const std::vector<char> &  annotated,
	double                     zero_annotated,
	double                     pi_gwas,
	double                     pi_qtl,
	double                     total_snp,
	double                     qtl_only_sampled,
	int                        max_iter = 1000,
	double                     a1_tol = 0.01)
{
	double a0 = log(pi_gwas / (1 - pi_gwas));
	double a1 = 0;
	double var0 = 0;
	double var1 = 0;
	double r0, r1;
	r0 = r1 = exp(a0);
	int iter = 0;
	const size_t m = gwas.odds.size();
	const double zero_unannotated = gwas.n_zero - zero_annotated;

	while (true) {
		iter++;
		// E step
		double pseudo_count = 1.0;
		double e0g0 = pseudo_count * (1 - pi_gwas) * (1 - pi_qtl);
		double e0g1 = pseudo_count * (1 - pi_qtl) * pi_gwas;
		double e1g0 = pseudo_count * (1 - pi_gwas) * pi_qtl;
		double e1g1 = pseudo_count * pi_gwas * pi_qtl;

		for (size_t j = 0; j < m; j++) {
			if (!annotated[j]) {
				double val = r0 * gwas.odds[j];
				// updated posterior with current prior given eqtl = 0
				val = val / (1 + val);
				e0g1 += val;
				e0g0 += 1 - val;
			} else {
				double val = r1 * gwas.odds[j];
				// updated posterior with current prior given eqtl = 1
				val = val / (1 + val);
				e1g1 += val;
				e1g0 += 1 - val;
			}
		}
		// The PIP-0 variants: one (not causal) count each.
		e0g0 += zero_unannotated;
		e1g0 += zero_annotated;

		// A GWAS PIP of 0 has posterior 0 under either prior, so each sampled
		// QTL-only variant adds exactly one (QTL, not causal) count.
		e1g0 += qtl_only_sampled;

		e0g0 += total_snp - (e0g0 + e0g1 + e1g0 + e1g1);
		if (!(e0g0 > 0)) {
			double nan = std::numeric_limits<double>::quiet_NaN();
			return EMResult{nan, nan, nan, nan, false, false, iter};
		}

		double a1_new = log(e1g1 * e0g0 / (e1g0 * e0g1));

		if (fabs(a1_new - a1) < a1_tol || iter >= max_iter) {
			bool converged = fabs(a1_new - a1) < a1_tol;
			a1 = a1_new;
			a0 = log(e0g1 / e0g0);
			var1 = (1.0 / e0g0 + 1.0 / e1g0 + 1.0 / e1g1 + 1.0 / e0g1);
			var0 = (1.0 / e0g1 + 1.0 / e0g0);
			return EMResult{a0, a1, var0, var1, converged, true, iter};
		}

		a1 = a1_new;
		a0 = log(e0g1 / e0g0);
		r0 = exp(a0);
		r1 = exp(a0 + a1);
	}
}

// The estimates, plus per-round diagnostics the R wrapper reports: the
// proportion of sampled QTL variants absent from the GWAS PIP vector, whether
// the round's EM converged, and how many rounds the variant total was too
// small for.
struct EnrichmentOutput {
	std::map<std::string, double> estimates;
	std::vector<double>           missing_proportion;
	std::vector<int>              em_converged;
	std::vector<int>              em_iterations;
	int                           invalid_rounds;
};

EnrichmentOutput qtl_enrichment_workhorse(
	const std::vector<SuSiEFit> &qtl_susie_fits,
	const std::vector<double> &  gwas_pip,
	double                       pi_gwas,
	double                       pi_qtl,
	double                       total_snp,
	int                          ImpN,
	double                       shrinkage_lambda,
	bool                         double_shrinkage = false,
	bool                         bessel_correction = true,
	int                          num_threads = 4,
	unsigned int                 base_seed = 0)
{

	std::vector<double> a0_vec(ImpN, 0.0);
	std::vector<double> v0_vec(ImpN, 0.0);
	std::vector<double> a1_vec(ImpN, 0.0);
	std::vector<double> v1_vec(ImpN, 0.0);
	std::vector<double> missing_proportion(ImpN, 0.0);
	std::vector<int>    em_converged(ImpN, 0);
	std::vector<int>    em_valid(ImpN, 0);
	std::vector<int>    em_iterations(ImpN, 0);
	const GwasOdds      gwas = prepare_gwas_odds(gwas_pip, pi_gwas);

	// QTL-only variants carry the codes -1 .. -n_qtl_only.
	int n_qtl_only = 0;
	for (const auto &fit : qtl_susie_fits) {
		for (int code : fit.codes) {
			n_qtl_only = std::max(n_qtl_only, -code);
		}
	}

	// Each round writes only its own slot of the vectors above, so the loop
	// needs no critical section, and nothing in it prints or raises.
	#pragma omp parallel for num_threads(num_threads)
	for (int k = 0; k < ImpN; k++) {
		// Per-round RNG seeded deterministically from the base seed plus the
		// round index, so results are reproducible for a given base_seed
		// regardless of which thread runs this round.
		std::mt19937 gen(base_seed + static_cast<unsigned int>(k));

		// Use QTL to annotate GWAS variants; a sampled variant the GWAS vector
		// does not carry is counted once, however many fits sampled it.
		std::vector<char> annotated(gwas.odds.size(), 0);
		std::vector<char> zero_seen;
		std::vector<char> qtl_only_seen(n_qtl_only, 0);
		int qtl_only_sampled = 0;
		int zero_annotated = 0;
		int missing_qtl_count = 0;
		int total_qtl_count = 0;

		for (size_t i = 0; i < qtl_susie_fits.size(); i++) {
			std::vector<int> qtn_codes = qtl_susie_fits[i].impute_qtn(gen);
			for (int code : qtn_codes) {
				if (code >= 0) {
					int at = gwas.place[code];
					if (at >= 0) {
						annotated[at] = 1;
					} else {
						// A PIP-0 variant, counted once however many fits
						// sampled it.
						if (zero_seen.empty()) {
							zero_seen.assign(gwas_pip.size(), 0);
						}
						if (!zero_seen[code]) {
							zero_seen[code] = 1;
							++zero_annotated;
						}
					}
				} else {
					++missing_qtl_count;
					char &seen = qtl_only_seen[-code - 1];
					if (!seen) {
						seen = 1;
						++qtl_only_sampled;
					}
				}
			}
			total_qtl_count += qtn_codes.size();
		}
		missing_proportion[k] = total_qtl_count > 0
		                            ? static_cast<double>(missing_qtl_count) / total_qtl_count
		                            : std::numeric_limits<double>::quiet_NaN();
		EMResult rst = run_EM(gwas, annotated, static_cast<double>(zero_annotated), pi_gwas, pi_qtl, total_snp, static_cast<double>(qtl_only_sampled));
		a0_vec[k] = rst.a0;
		a1_vec[k] = rst.a1;
		v0_vec[k] = rst.var0;
		v1_vec[k] = rst.var1;
		em_converged[k] = rst.converged ? 1 : 0;
		em_valid[k] = rst.valid ? 1 : 0;
		em_iterations[k] = rst.iterations;
	}
	int invalid_rounds = ImpN - std::accumulate(em_valid.begin(), em_valid.end(), 0);

	// Apply outlier filtering if shrinkage is specified.
	// Use index-based filtering to keep all four vectors aligned.
	double pv = (shrinkage_lambda == 0) ? -1.0 : 1.0 / shrinkage_lambda;
	std::vector<size_t> kept_indices;
	if (pv > 0) {
		kept_indices = filter_outlier_indices(a1_vec, v1_vec, pv, bessel_correction);
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
	// Rubin's rules need two imputations to estimate a between-round variance;
	// with one surviving round only the within-round variance is measurable,
	// and dividing by (m - 1) would make every downstream estimate NaN.
	if (m > 1) {
		bv0 /= (m - 1);
		bv1 /= (m - 1);
	} else {
		bv0 = 0;
		bv1 = 0;
	}
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
			bv1_shrink = (m > 1) ? bv1_shrink / (m - 1) : 0.0;
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

	// The intercept matches piGwas through the unshrunk enrichment, as in
	// fastenloc; the coloc priors below then pair it with the shrunk one.
	a0_est = log(pi_gwas / (1 + pi_qtl * exp(a1_est_ns) - pi_qtl - pi_gwas));

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

	return EnrichmentOutput{output_map, missing_proportion, em_converged, em_iterations, invalid_rounds};
}

#endif // QTL_ENRICHMENT_HPP
