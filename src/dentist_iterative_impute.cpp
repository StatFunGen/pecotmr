// This function implements, to our knowledge, the methods decribed in the DENTIST paper
// https://github.com/Yves-CHEN/DENTIST/tree/master#Citations
// Some codes are adapted and rewritten from https://github.com/Yves-CHEN/DENTIST/tree/master
// to fit the Rcpp implementation.
// The code reflects our understanding and interpretation of DENTIST method which may difer in details
// from the author's original proposal, although in various tests we find that our implementation and the
// original results are mostly identical
#include <RcppArmadillo.h>
#include <omp.h> // Required for parallel processing
#include <algorithm>
#include <random>
#include <vector>
#include <gsl/gsl_cdf.h>
#include <unordered_set>
#include <ctime>
#include <cstdlib>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

std::vector<size_t> generateSetOfNumbers(size_t size, unsigned int seed) {
	std::vector<size_t> indexes(size);
	std::iota(indexes.begin(), indexes.end(), 0);
	std::mt19937 gen(seed);
	std::shuffle(indexes.begin(), indexes.end(), gen);
	return indexes;
}

// Get a quantile value
double getQuantile(const std::vector<double>& dat, double whichQuantile) {
	std::vector<double> sortedData = dat;
	std::sort(sortedData.begin(), sortedData.end());
	size_t pos = ceil(sortedData.size() * whichQuantile) - 1;
	return sortedData.at(pos);
}

// Get a quantile value based on grouping
double getQuantile2(const std::vector<double>& dat, const std::vector<size_t>& grouping, double whichQuantile, bool invert_grouping = false) {
	std::vector<double> filteredData;
	for (size_t i = 0; i < dat.size(); ++i) {
		// Apply grouping filter with optional inversion
		if ((grouping[i] == 1) != invert_grouping) {
			filteredData.push_back(dat[i]);
		}
	}
	if (filteredData.size() < 50) return 0;
	return getQuantile(filteredData, whichQuantile);
}

// Get a quantile value based on grouping
double getQuantile2_chen_et_al(const std::vector<double> &dat, std::vector<size_t> grouping, double whichQuantile)
{
	size_t sum = std::accumulate(grouping.begin(), grouping.end(), 0);

	if (sum < 50)
	{
		return 0;
	}

	std::vector<double> diff2;
	for (size_t i = 0; i < dat.size(); i++)
	{
		if (grouping[i] == 1)
			diff2.push_back(dat[i]);
	}
	return getQuantile(diff2, whichQuantile);
}

// Calculate minus log p-value of chi-squared statistic
double minusLogPvalueChisq2(double stat) {
	double p = 1.0 - gsl_cdf_chisq_P(stat, 1.0);
	return -log10(p);
}

// Perform one iteration of the algorithm, assuming LD_mat is an arma::mat
void oneIteration(const arma::mat& LD_mat, const std::vector<size_t>& idx, const std::vector<size_t>& idx2,
                  const arma::vec& zScore, arma::vec& imputedZ, arma::vec& rsqList, arma::vec& zScore_e,
                  size_t nSample, float probSVD, int ncpus, bool verbose) {
	if (verbose) {
		Rcpp::Rcout << "LD_mat dimensions: " << LD_mat.n_rows << " x " << LD_mat.n_cols << std::endl;
		Rcpp::Rcout << "idx size: " << idx.size() << std::endl;
		Rcpp::Rcout << "idx2 size: " << idx2.size() << std::endl;
		Rcpp::Rcout << "zScore size: " << zScore.size() << std::endl;
		Rcpp::Rcout << "imputedZ size: " << imputedZ.size() << std::endl;
		Rcpp::Rcout << "rsqList size: " << rsqList.size() << std::endl;
		Rcpp::Rcout << "zScore_e size: " << zScore_e.size() << std::endl;
	}

	// DEBUG: Print thread information in oneIteration
	time_t rawtime2;
	struct tm * timeinfo2;
	char buffer2[80];
	time(&rawtime2);
	timeinfo2 = localtime(&rawtime2);
	strftime(buffer2, sizeof(buffer2), "%Y-%m-%d %H:%M:%S", timeinfo2);
	
	const char* omp_env2 = getenv("OMP_NUM_THREADS");
	Rcpp::Rcout << "[DEBUG oneIteration " << buffer2 << "] Thread settings:" << std::endl;
	Rcpp::Rcout << "  OMP_NUM_THREADS env var: " << (omp_env2 ? omp_env2 : "not set") << std::endl;
	Rcpp::Rcout << "  ncpus parameter: " << ncpus << std::endl;
	Rcpp::Rcout << "  omp_get_max_threads(): " << omp_get_max_threads() << std::endl;
	
	int nProcessors = omp_get_max_threads();
	if (ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);
	
	Rcpp::Rcout << "  Setting omp_set_num_threads(" << nProcessors << ")" << std::endl;
	Rcpp::Rcout << "  After setting, omp_get_num_threads(): " << omp_get_num_threads() << std::endl;

	size_t K = std::min(static_cast<size_t>(idx.size()), nSample) * probSVD;

	arma::vec zScore_eigen(idx.size());
	arma::mat LD_it(idx2.size(), idx.size());
	arma::mat VV(idx.size(), idx.size());

	// Check dimensions before filling LD_it and VV matrices
	if (idx2.size() > LD_mat.n_rows || idx.size() > LD_mat.n_cols) {
		Rcpp::stop("Inconsistent dimensions between LD_mat and idx2/idx in oneIteration()");
	}

	// Check if any index in idx or idx2 is greater than or equal to zscore size
	for (size_t i = 0; i < idx.size(); ++i) {
		if (idx[i] >= zScore.size()) {
			Rcpp::stop("Invalid index in idx: " + std::to_string(idx[i]));
		}
	}
	for (size_t i = 0; i < idx2.size(); ++i) {
		if (idx2[i] >= zScore.size()) {
			Rcpp::stop("Invalid index in idx2: " + std::to_string(idx2[i]));
		}
	}

	if (verbose) {
		Rcpp::Rcout << "Filling LD_it and VV matrices" << std::endl;
	}

	// Fill LD_it and VV matrices using direct indexing
#pragma omp parallel for collapse(2)
	for (size_t i = 0; i < idx2.size(); i++) {
		for (size_t k = 0; k < idx.size(); k++) {
			// if (verbose) {
			//	Rcpp::Rcout << "Filling LD_it: i = " << i << ", k = " << k << ", idx2[i] = " << idx2[i] << ", idx[k] = " << idx[k] << std::endl;
			//}
			LD_it(i, k) = LD_mat.at(idx2[i] * LD_mat.n_cols + idx[k]);
		}
	}

#pragma omp parallel for
	for (size_t i = 0; i < idx.size(); i++) {
		//if (verbose) {
		//	Rcpp::Rcout << "Filling VV: i = " << i << ", idx[i] = " << idx[i] << std::endl;
		//}
		zScore_eigen(i) = zScore[idx[i]];
		for (size_t j = 0; j < idx.size(); j++) {
			//if (verbose) {
			//	Rcpp::Rcout << "Filling VV: i = " << i << ", j = " << j << ", idx[i] = " << idx[i] << ", idx[j] = " << idx[j] << std::endl;
			//}
			VV(i, j) = LD_mat.at(idx[i] * LD_mat.n_rows + idx[j]);
		}
	}
	if (verbose) {
		Rcpp::Rcout << "Performing eigen decomposition" << std::endl;
	}

	// Eigen decomposition
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, VV);

	int nRank = eigvec.n_rows;
	int nZeros = arma::sum(eigval < 0.0001);
	nRank -= nZeros;
	K = std::min(K, static_cast<size_t>(nRank));

	if (verbose) {
		Rcpp::Rcout << "Rank: " << nRank << ", Zeros: " << nZeros << ", K: " << K << std::endl;
	}

	if (K <= 1) {
		Rcpp::stop("Rank of eigen matrix <= 1");
	}
	arma::mat ui = arma::eye<arma::mat>(eigvec.n_rows, K);
	arma::mat wi = arma::eye<arma::mat>(K, K);
	for (size_t m = 0; m < K; ++m) {
		int j = eigvec.n_rows - m - 1;
		ui.col(m) = eigvec.col(j);
		wi(m, m) = 1.0 / eigval(j);
	}

	if (verbose) {
		Rcpp::Rcout << "Calculating imputed Z scores and R squared values" << std::endl;
	}

	// Calculate imputed Z scores and R squared values
	arma::mat beta = LD_it * ui * wi;
	arma::vec zScore_eigen_imp = beta * (ui.t() * zScore_eigen);
	arma::mat product = beta * (ui.t() * LD_it.t());
	arma::vec rsq_eigen = product.diag();

#pragma omp parallel for
	for (size_t i = 0; i < idx2.size(); ++i) {
		imputedZ[idx2[i]] = zScore_eigen_imp(i);
		rsqList[idx2[i]] = rsq_eigen(i);
		rsqList[idx2[i]] = std::min(rsq_eigen(i), 1.0); // Ensure rsq does not exceed 1
		if (rsq_eigen(i) >= 1) {
			// Handle the case where rsq_eigen is unexpectedly high
			Rcpp::warning("Adjusted rsq_eigen value exceeding 1: " + std::to_string(rsq_eigen(i)));
		}
		size_t j = idx2[i];
		zScore_e[j] = (zScore[j] - imputedZ[j]) / std::sqrt(LD_mat(j, j) - rsqList[j]);
	}

}

/**
 * @brief Executes DENTIST algorithm for quality control in GWAS summary data: the iterative imputation function.
 *
 * DENTIST (Detecting Errors iN analyses of summary staTISTics) identifies and removes problematic variants
 * in GWAS summary data by comparing observed GWAS statistics to predicted values based on linkage disequilibrium (LD)
 * information from a reference panel. It helps detect genotyping/imputation errors, allelic errors, and heterogeneity
 * between GWAS and LD reference samples, improving the reliability of subsequent analyses.
 *
 * @param LD_mat The linkage disequilibrium (LD) matrix from a reference panel, as an arma::mat object.
 * @param nSample The sample size used in the GWAS whose summary statistics are being analyzed.
 * @param zScore A vector of Z-scores from GWAS summary statistics.
 * @param pValueThreshold Threshold for the p-value below which variants are considered for quality control.
 * @param propSVD Proportion of singular value decomposition (SVD) components retained in the analysis.
 * @param gcControl A boolean flag to apply genetic control corrections.
 * @param nIter The number of iterations to run the DENTIST algorithm.
 * @param gPvalueThreshold P-value threshold for grouping variants into significant and null categories.
 * @param ncpus The number of CPU cores to use for parallel processing.
 * @param seed Seed for random number generation, affecting the selection of variants for analysis.
 * @param verbose A boolean flag to enable verbose output for debugging.
 *
 * @return A List object containing:
 * - original_z: A vector of original Z-scores for each marker.
 * - imputed_z: A vector of imputed Z-scores for each marker.
 * - z_diff: A vector of outlier test z-scores
 * - rsq: A vector of R-squared values for each marker, indicating goodness of fit.
 * - iter_to_correct: An integer vector indicating the iteration in which each marker passed the quality control.
 *
 * @note The function is designed for use in Rcpp and requires Armadillo for matrix operations and OpenMP for parallel processing.
 */

// [[Rcpp::export]]
List dentist_iterative_impute(const arma::mat& LD_mat, size_t nSample, const arma::vec& zScore,
                              double pValueThreshold, float propSVD, bool gcControl, int nIter,
                              double gPvalueThreshold, int ncpus, int seed, bool correct_chen_et_al_bug,
                              bool verbose = false) {
	if (verbose) {
		Rcpp::Rcout << "LD_mat dimensions: " << LD_mat.n_rows << " x " << LD_mat.n_cols << std::endl;
		Rcpp::Rcout << "nSample: " << nSample << std::endl;
		Rcpp::Rcout << "zScore size: " << zScore.size() << std::endl;
		Rcpp::Rcout << "pValueThreshold: " << pValueThreshold << std::endl;
		Rcpp::Rcout << "propSVD: " << propSVD << std::endl;
		Rcpp::Rcout << "gcControl: " << gcControl << std::endl;
		Rcpp::Rcout << "nIter: " << nIter << std::endl;
		Rcpp::Rcout << "gPvalueThreshold: " << gPvalueThreshold << std::endl;
		Rcpp::Rcout << "ncpus: " << ncpus << std::endl;
		Rcpp::Rcout << "seed: " << seed << std::endl;
		Rcpp::Rcout << "correct_chen_et_al_bug: " << correct_chen_et_al_bug << std::endl;
	}

	// Set number of threads for parallel processing
	// DEBUG: Print thread information
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
	
	const char* omp_env = getenv("OMP_NUM_THREADS");
	Rcpp::Rcout << "[DEBUG dentist_iterative_impute " << buffer << "] Thread settings:" << std::endl;
	Rcpp::Rcout << "  OMP_NUM_THREADS env var: " << (omp_env ? omp_env : "not set") << std::endl;
	Rcpp::Rcout << "  ncpus parameter: " << ncpus << std::endl;
	Rcpp::Rcout << "  omp_get_max_threads(): " << omp_get_max_threads() << std::endl;
	
	int nProcessors = omp_get_max_threads();
	if (ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);
	
	Rcpp::Rcout << "  Setting omp_set_num_threads(" << nProcessors << ")" << std::endl;
	Rcpp::Rcout << "  After setting, omp_get_num_threads(): " << omp_get_num_threads() << std::endl;

	size_t markerSize = zScore.size();
	std::vector<size_t> randOrder = generateSetOfNumbers(markerSize, seed);
	std::vector<size_t> idx, idx2;
	idx.reserve(markerSize / 2);
	idx2.reserve(markerSize / 2);
	std::vector<size_t> fullIdx(randOrder.begin(), randOrder.end());

	// Determining indices for partitioning
	for (size_t i = 0; i < markerSize; ++i) {
		if (randOrder[i] > markerSize / 2) idx.push_back(i);
		else idx2.push_back(i);
	}

	if (verbose) {
		Rcpp::Rcout << "Indices partitioned" << std::endl;
	}

	std::vector<size_t> groupingGWAS(markerSize, 0);
	for (size_t i = 0; i < markerSize; ++i) {
		if (minusLogPvalueChisq2(std::pow(zScore(i), 2)) > -log10(gPvalueThreshold)) {
			groupingGWAS[i] = 1;
		}
	}

	if (verbose) {
		Rcpp::Rcout << "Grouping GWAS finished" << std::endl;
	}

	arma::vec imputedZ = arma::zeros<arma::vec>(markerSize);
	arma::vec rsq = arma::zeros<arma::vec>(markerSize);
	arma::vec zScore_e = arma::zeros<arma::vec>(markerSize);
	arma::ivec iterID = arma::zeros<arma::ivec>(markerSize);

	std::vector<double> diff(idx2.size());
	std::vector<size_t> grouping_tmp(idx2.size());

	for (int t = 0; t < nIter; ++t) {
		// Perform iteration with current subsets
		if (verbose) {
			Rcpp::Rcout << "\nIteration " << t << std::endl;
			Rcpp::Rcout << "Performing  with current subsets" << std::endl;
			Rcpp::Rcout << "Performing oneIteration()" << std::endl;
		}

		oneIteration(LD_mat, idx, idx2, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus, verbose);

		diff.resize(idx2.size());
		grouping_tmp.resize(idx2.size());

		// Assess differences and grouping for thresholding
		for (size_t i = 0; i < idx2.size(); ++i) {
			diff[i] = std::abs(zScore_e[idx2[i]]);
			grouping_tmp[i] = groupingGWAS[idx2[i]];
		}

		if (verbose) {
			Rcpp::Rcout << "Assessing differences and grouping for thresholding" << std::endl;
		}

		double threshold = getQuantile(diff, 0.995);
		double threshold1, threshold0;
		/*
		        In the original DENTIST method, whenever you call !grouping_tmp, it is going to change the original value of grouping_tmp as well.
		        For example, if grouping_tmp is (0,0,1,1,1), and you run:
		        double threshold0 = getQuantile2 <double> (diff,!grouping_tmp , (99.5/100.0)) ;
		        then your grouping_tmp will become (1,1,0,0,0) even you are just calling it in the function.
		        https://github.com/Yves-CHEN/DENTIST/blob/2fefddb1bbee19896a30bf56229603561ea1dba8/main/inversion.cpp#L647
		        https://github.com/Yves-CHEN/DENTIST/blob/2fefddb1bbee19896a30bf56229603561ea1dba8/main/inversion.cpp#L675
		        Thus if we correct the original DENTIST code, i.e., correct_chen_et_al_bug = TRUE,
		                we go through our function, getQuantile2, which doesn't have this issue
		                else, i.e., correct_chen_et_al_bug = TRUE, it goes through the original function getQuantile2_chen_et_al
		 */
		if (correct_chen_et_al_bug) {
			threshold1 = getQuantile2(diff, grouping_tmp, 0.995, false);
			threshold0 = getQuantile2(diff, grouping_tmp, 0.995, true);
		} else {
			threshold1 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
			std::transform(grouping_tmp.begin(), grouping_tmp.end(), grouping_tmp.begin(), [](size_t val) {
				return 1 - val;
			});
			threshold0 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
		}

		if (threshold1 == 0) {
			threshold1 = threshold;
			threshold0 = threshold;
		}
		if (correct_chen_et_al_bug || nIter - 2 >= 0) {
			/*In the original DENTIST method, if t=0 (first iteration) and nIter is 1,
			   t is defined as a size_t (unassigned integer)
			   https://github.com/Yves-CHEN/DENTIST/blob/2fefddb1bbee19896a30bf56229603561ea1dba8/main/inversion.cpp#L628
			   and it will treat t (which is 0) no larger than nIter-2 (which is -1) which is wrong
			   Thus if we correct the original DENTIST code, i.e., correct_chen_et_al_bug = TRUE, or when nIter - 2 >=0,
			   it will compare t and nIter as we expect.
			   and if we want to keep the original DENTIST code, i.e., correct_chen_et_al_bug = FALSE, then it will skip this if condition for t > nIter - 2
			 */
			if (t > nIter - 2) {
				threshold0 = threshold;
				threshold1 = threshold;
			}
		}

		if (verbose) {
			Rcpp::Rcout << "Thresholds calculated: " << threshold << ", " << threshold1 << ", " << threshold0 << std::endl;
			Rcpp::Rcout << "Applying threshold-based filtering for QC" << std::endl;
		}

		// Apply threshold-based filtering for QC
		std::vector<size_t> idx2_QCed;
		for (size_t i = 0; i < diff.size(); ++i) {
			if ((grouping_tmp[i] == 1 && diff[i] <= threshold1) ||
			    (grouping_tmp[i] == 0 && diff[i] <= threshold0)) {
				idx2_QCed.push_back(idx2[i]);
			}
		}

		// Perform another iteration with updated sets of indices (idx and idx2_QCed)
		if (verbose) {
			Rcpp::Rcout << "Performing oneIteration() with updated sets of indices" << std::endl;
		}

		oneIteration(LD_mat, idx2_QCed, idx, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus, verbose);

		if (verbose) {
			Rcpp::Rcout << "Recalculating differences and groupings after the iteration" << std::endl;
		}

		// Recalculate differences and groupings after the iteration
		diff.resize(fullIdx.size());
		grouping_tmp.resize(fullIdx.size());

		for (size_t i = 0; i < fullIdx.size(); ++i) {
			diff[i] = std::abs(zScore_e[fullIdx[i]]);
			grouping_tmp[i] = groupingGWAS[fullIdx[i]];
		}

		if (verbose) {
			Rcpp::Rcout << "Re-determining thresholds based on the recalculated differences and groupings" << std::endl;
		}

		// Re-determine thresholds based on the recalculated differences and groupings
		threshold = getQuantile(diff, 0.995);
		if (correct_chen_et_al_bug) {
			threshold1 = getQuantile2(diff, grouping_tmp, 0.995, false);
			threshold0 = getQuantile2(diff, grouping_tmp, 0.995, true);
		} else {
			threshold1 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
			std::transform(grouping_tmp.begin(), grouping_tmp.end(), grouping_tmp.begin(), [](size_t val) {
				return 1 - val;
			});
			threshold0 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
		}


		if (correct_chen_et_al_bug || nIter - 2 >= 0) {
			if (t > nIter - 2) {
				threshold0 = threshold;
				threshold1 = threshold;
			}
		}

		if (threshold1 == 0) {
			threshold1 = threshold;
			threshold0 = threshold;
		}

		if (verbose) {
			Rcpp::Rcout << "Adjusting for genetic control and inflation factor if necessary" << std::endl;
		}

		// Adjust for genetic control and inflation factor if necessary
		std::vector<double> chisq(fullIdx.size());
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			chisq[i] = std::pow(zScore_e[fullIdx[i]], 2);
		}

		if (chisq.size() <= 2) {
			Rcpp::stop("chisq.size() must be greater than 2.");
		}

		// Calculate the median chi-squared value as the inflation factor
		std::nth_element(chisq.begin(), chisq.begin() + chisq.size() / 2, chisq.end());
		double medianChisq = chisq[chisq.size() / 2];
		double inflationFactor = medianChisq / 0.46;

		std::vector<size_t> fullIdx_tmp;
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			if (gcControl) {
				// When gcControl is true, check if the variant passes the adjusted threshold
				if (!(diff[i] > threshold && minusLogPvalueChisq2(chisq[i] / inflationFactor) > -log10(pValueThreshold))) {
					fullIdx_tmp.push_back(fullIdx[i]);
				}
			} else {
				// When gcControl is false, simply check if the variant passes the basic threshold
				if (minusLogPvalueChisq2(chisq[i]) < -log10(pValueThreshold)) {
					if ((groupingGWAS[fullIdx[i]] == 1 && diff[i] <= threshold1) ||
					    (groupingGWAS[fullIdx[i]] == 0 && diff[i] <= threshold0)) {
						fullIdx_tmp.push_back(fullIdx[i]);
						iterID[fullIdx[i]]++;
					}
				}
			}
		}

		// Update the indices for the next iteration based on filtering criteria
		fullIdx = fullIdx_tmp;
		if (verbose) {
			Rcpp::Rcout << "fullIdx size: " << fullIdx.size() << std::endl;
		}
		randOrder = generateSetOfNumbers(fullIdx.size(), seed + t * seed);
		idx.clear();
		idx2.clear();
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			if (randOrder[i] > fullIdx.size() / 2) idx.push_back(fullIdx[i]);
			else idx2.push_back(fullIdx[i]);
		}
	}

	return List::create(Named("original_z") = zScore,
	                    Named("imputed_z") = imputedZ,
	                    Named("rsq") = rsq,
	                    Named("z_diff") = zScore_e,
	                    Named("iter_to_correct") = iterID);
}
