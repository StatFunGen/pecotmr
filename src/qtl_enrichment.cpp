#include "qtl_enrichment.hpp"

[[cpp11::register]]
cpp11::writable::list qtl_enrichment_rcpp(
	SEXP r_gwas_pip, SEXP r_qtl_susie_fit,
	double pi_gwas = 0, double pi_qtl = 0,
	int ImpN = 25, double shrinkage_lambda = 1.0,
	bool double_shrinkage = false,
	bool bessel_correction = true,
	int num_threads = 1)
{
	// Convert r_gwas_pip to C++ type
	doubles gwas_pip_vec(r_gwas_pip);
	std::vector<double> gwas_pip = cpp11::as_cpp<std::vector<double>>(gwas_pip_vec);
	cpp11::strings pip_names(gwas_pip_vec.attr("names"));
	std::vector<std::string> gwas_pip_names;
	gwas_pip_names.reserve(pip_names.size());
	for (int i = 0; i < pip_names.size(); ++i) {
		gwas_pip_names.push_back(std::string(pip_names[i]));
	}

	// Convert r_qtl_susie_fit to C++ type
	list susie_fit_list(r_qtl_susie_fit);
	std::vector<SuSiEFit> susie_fits;

	for (int i = 0; i < susie_fit_list.size(); ++i) {
		SuSiEFit susie_fit(susie_fit_list[i]);
		susie_fits.push_back(susie_fit);
	}

	std::map<std::string, double> output = qtl_enrichment_workhorse(susie_fits, gwas_pip, gwas_pip_names, pi_gwas, pi_qtl, ImpN, shrinkage_lambda, double_shrinkage, bessel_correction, num_threads);

	// Convert std::map to list
	using namespace cpp11::literals;
	writable::list output_list;
	writable::strings names;
	for (auto const& element : output) {
		output_list.push_back(cpp11::as_sexp(element.second));
		names.push_back(element.first);
	}
	output_list.attr("names") = names;

	return output_list;
}
