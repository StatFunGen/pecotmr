#include "qtl_enrichment.hpp"

// [[Rcpp::export]]
Rcpp::List qtl_enrichment_rcpp(
	SEXP r_gwas_pip, SEXP r_qtl_susie_fit,
	double pi_gwas = 0, double pi_qtl = 0,
	int ImpN = 25, double shrinkage_lambda = 1.0,
	int num_threads = 1)
{
	// Convert r_gwas_pip to C++ type
	Rcpp::NumericVector gwas_pip_vec = Rcpp::as<Rcpp::NumericVector>(r_gwas_pip);
	std::vector<double> gwas_pip = Rcpp::as<std::vector<double> >(gwas_pip_vec);
	std::vector<std::string> gwas_pip_names = Rcpp::as<std::vector<std::string> >(gwas_pip_vec.names());

	// Convert r_qtl_susie_fit to C++ type
	Rcpp::List susie_fit_list(r_qtl_susie_fit);
	std::vector<SuSiEFit> susie_fits;

	for (int i = 0; i < susie_fit_list.size(); ++i) {
		SuSiEFit susie_fit(Rcpp::wrap(susie_fit_list[i]));
		susie_fits.push_back(susie_fit);
	}

	std::map<std::string, double> output = qtl_enrichment_workhorse(susie_fits, gwas_pip, gwas_pip_names, pi_gwas, pi_qtl, ImpN, shrinkage_lambda, num_threads);

	// Convert std::map to Rcpp::List
	Rcpp::List output_list;
	for (auto const& element : output) {
		output_list[element.first] = element.second;
	}

	return output_list;
}
