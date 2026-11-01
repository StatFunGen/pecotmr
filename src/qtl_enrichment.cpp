#include "qtl_enrichment.h"

// `rGwasPip` holds the GWAS PIP values; `rQtlSusieFit` is a list of fits, each
// with `alpha`, `prior_variance` and integer `codes` that index its variables
// into `rGwasPip` (negative for variants it does not carry). `totalVariants`
// is the variant count the EM remainder is taken over.
[[cpp11::register]]
cpp11::writable::list qtlEnrichmentCpp(
	SEXP rGwasPip, SEXP rQtlSusieFit,
	double piGwas, double piQtl, double totalVariants,
	int ImpN = 25, double shrinkageLambda = 1.0,
	bool doubleShrinkage = false,
	bool besselCorrection = true,
	int numThreads = 1,
	cpp11::sexp seed = R_NilValue)
{
	std::vector<double> gwas_pip = cpp11::as_cpp<std::vector<double>>(doubles(rGwasPip));
	const int n_gwas = static_cast<int>(gwas_pip.size());

	// Convert rQtlSusieFit to C++ type, refusing a code that points past the
	// GWAS vector, since a sampled variant's code indexes it directly.
	list susie_fit_list(rQtlSusieFit);
	std::vector<SuSiEFit> susie_fits;
	susie_fits.reserve(susie_fit_list.size());
	for (int i = 0; i < susie_fit_list.size(); ++i) {
		susie_fits.emplace_back(susie_fit_list[i]);
		for (int code : susie_fits.back().codes) {
			if (code >= n_gwas) {
				cpp11::stop("A variant code points past the end of the GWAS PIP vector.");
			}
		}
	}

	// Resolve the base seed: use the user-provided value, otherwise a
	// nondeterministic seed from std::random_device.
	unsigned int base_seed = 0;
	if (seed != R_NilValue) {
		base_seed = cpp11::as_cpp<unsigned int>(seed);
	} else {
		base_seed = std::random_device{}();
	}

	EnrichmentOutput output = qtl_enrichment_workhorse(susie_fits, gwas_pip, piGwas, piQtl, totalVariants, ImpN, shrinkageLambda, doubleShrinkage, besselCorrection, numThreads, base_seed);

	// The estimates (std::map order), then the per-round diagnostics.
	writable::list output_list;
	writable::strings names;
	for (auto const& element : output.estimates) {
		output_list.push_back(cpp11::as_sexp(element.second));
		names.push_back(element.first);
	}
	writable::doubles missing(ImpN);
	writable::logicals converged(ImpN);
	writable::integers iterations(ImpN);
	for (int k = 0; k < ImpN; ++k) {
		missing[k] = output.missing_proportion[k];
		converged[k] = output.em_converged[k] != 0;
		iterations[k] = output.em_iterations[k];
	}
	output_list.push_back(missing);
	names.push_back("Missing proportion per MI round");
	output_list.push_back(converged);
	names.push_back("EM converged per MI round");
	output_list.push_back(iterations);
	names.push_back("EM iterations per MI round");
	output_list.push_back(cpp11::as_sexp(output.invalid_rounds));
	names.push_back("Invalid MI rounds");
	output_list.attr("names") = names;

	return output_list;
}
