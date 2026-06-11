#include "qtl_enrichment.hpp"

[[cpp11::register]]
cpp11::writable::list qtlEnrichmentRcpp(
	SEXP rGwasPip, SEXP rQtlSusieFit,
	double piGwas = 0, double piQtl = 0,
	int ImpN = 25, double shrinkageLambda = 1.0,
	bool doubleShrinkage = false,
	bool besselCorrection = true,
	int numThreads = 1)
{
	// Convert rGwasPip to C++ type
	doubles gwas_pip_vec(rGwasPip);
	std::vector<double> gwas_pip = cpp11::as_cpp<std::vector<double>>(gwas_pip_vec);
	cpp11::strings pip_names(gwas_pip_vec.attr("names"));
	std::vector<std::string> gwas_pip_names;
	gwas_pip_names.reserve(pip_names.size());
	for (int i = 0; i < pip_names.size(); ++i) {
		gwas_pip_names.push_back(std::string(pip_names[i]));
	}

	// Convert rQtlSusieFit to C++ type
	list susie_fit_list(rQtlSusieFit);
	std::vector<SuSiEFit> susie_fits;

	for (int i = 0; i < susie_fit_list.size(); ++i) {
		SuSiEFit susie_fit(susie_fit_list[i]);
		susie_fits.push_back(susie_fit);
	}

	std::map<std::string, double> output = qtl_enrichment_workhorse(susie_fits, gwas_pip, gwas_pip_names, piGwas, piQtl, ImpN, shrinkageLambda, doubleShrinkage, besselCorrection, numThreads);

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
