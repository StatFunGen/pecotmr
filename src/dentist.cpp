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
#include <numeric> // For std::iota
#include <gsl/gsl_cdf.h>
#include <fstream> // Include the header for file operations
#include <iostream>
#include <string>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;


// Assuming sort_indexes is defined as provided
std::vector<size_t> sort_indexes(const std::vector<int>& v, unsigned int theSize) {
	std::vector<size_t> idx(theSize);
	std::iota(idx.begin(), idx.end(), 0);
	std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
		return v[i1] < v[i2];
	});
	return idx;
}

// Improved generateSetOfNumbers function using C++11 random
std::vector<size_t> generateSetOfNumbers(int SIZE, int seed) {
	std::vector<int> numbers(SIZE, 0);
	std::mt19937 rng(seed); // Mersenne Twister: Good quality random number generator
	std::uniform_int_distribution<int> dist(0, INT_MAX);

	// Generate the first random number
	numbers[0] = dist(rng);
	for (int index = 1; index < SIZE; index++) {
		int tempNum;
		do {
			tempNum = dist(rng); // Generate a new random number
			// Check for uniqueness in the current list of generated numbers
			bool isUnique = true;
			for (int index2 = 0; index2 < index; index2++) {
				if (tempNum == numbers[index2]) {
					isUnique = false;
					break;
				}
			}
			// If the number is not unique, force the loop to try again
			if (!isUnique) tempNum = -1;
		} while (tempNum == -1);
		// Assign the unique number to the list
		numbers[index] = tempNum;
	}

	// Sort the indices of 'numbers' based on their values
	return sort_indexes(numbers, SIZE);
}

// Get a quantile value
double getQuantile(const std::vector<double>& dat, double whichQuantile) {
	std::vector<double> sortedData = dat;
	std::sort(sortedData.begin(), sortedData.end());
	size_t pos = ceil(sortedData.size() * whichQuantile) - 1;
	return sortedData.at(pos);
}

// Get a quantile value based on grouping
double getQuantile2(const std::vector<double>& dat, const std::vector<uint>& grouping, double whichQuantile, bool invert_grouping = false) {
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

// Calculate minus log p-value of chi-squared statistic
double minusLogPvalueChisq2(double stat) {
	double p = 1.0 - gsl_cdf_chisq_P(stat, 1.0);
	return -log10(p);
}

// Perform one iteration of the algorithm, assuming LDmat is an arma::mat
void oneIteration(const arma::mat& LDmat, const std::vector<uint>& idx, const std::vector<uint>& idx2,
                  const arma::vec& zScore, arma::vec& imputedZ, arma::vec& rsqList, arma::vec& zScore_e,
                  uint nSample, float probSVD, int ncpus, int t, int rrr) {
	int nProcessors = omp_get_max_threads();
	if(ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);

	uint K = std::min(static_cast<uint>(idx.size()), nSample) * probSVD;

	arma::mat LD_it(idx2.size(), idx.size());
	arma::vec zScore_eigen(idx.size());
	arma::mat VV(idx.size(), idx.size());

	// Fill LD_it and VV matrices using direct indexing
    #pragma omp parallel for collapse(2)
	for (size_t i = 0; i < idx2.size(); i++) {
		for (size_t k = 0; k < idx.size(); k++) {
			LD_it(i, k) = LDmat(idx2[i], idx[k]);
		}
	}
    

    /* 
    ============= Generate the output LD_it ============
    add by Rui on 20240308 (start)
    */
    // Add this code snippet within the oneIteration function after populating LD_it matrix
    // arguments:
    // t:   the index of interaction
    // rrr: the first (1) or second (2) time that the oneIteraction function is called within each iteraction
    std::string outputFileName = "/home/rd2972/private_data/20240300_rss_qc_imputation/DENTIST/per_iteration/dentist_Rcpp/DENTIST_Rcpp_output_" + std::to_string(t) +  std::to_string(rrr) + ".txt";
    std::ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFileName << std::endl;
        return;
    }

    // Write LD_it matrix to the output file for this iteration
    outputFile << LD_it << std::endl;

    // Close the output file
    outputFile.close();
    
    /* 
    ============= Generate the output LD_it ============
    add by Rui on 20240308 (end)
    */
    
    #pragma omp parallel for
	for (size_t i = 0; i < idx.size(); i++) {
		zScore_eigen(i) = zScore[idx[i]];
		for (size_t j = 0; j < idx.size(); j++) {
			VV(i, j) = LDmat(idx[i], idx[j]);
		}
	}

	// Eigen decomposition
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, VV);

	int nRank = eigvec.n_rows;
	int nZeros = arma::sum(eigval < 0.0001);
	nRank -= nZeros;
	K = std::min(K, static_cast<uint>(nRank));
	if (K <= 1) {
		Rcpp::stop("Rank of eigen matrix <= 1");
	}
	arma::mat ui = arma::eye(eigvec.n_rows, K);
	arma::mat wi = arma::eye(K, K);
	for (uint m = 0; m < K; ++m) {
		int j = eigvec.n_rows - m - 1;
		ui.col(m) = eigvec.col(j);
		wi(m, m) = 1.0 / eigval(j);
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
		uint j = idx2[i];
		zScore_e[j] = (zScore[j] - imputedZ[j]) / std::sqrt(LDmat(j, j) - rsqList[j]);
	}
}

/**
 * @brief Executes DENTIST algorithm for quality control in GWAS summary data.
 *
 * DENTIST (Detecting Errors iN analyses of summary staTISTics) identifies and removes problematic variants
 * in GWAS summary data by comparing observed GWAS statistics to predicted values based on linkage disequilibrium (LD)
 * information from a reference panel. It helps detect genotyping/imputation errors, allelic errors, and heterogeneity
 * between GWAS and LD reference samples, improving the reliability of subsequent analyses.
 *
 * @param LDmat The linkage disequilibrium (LD) matrix from a reference panel, as an arma::mat object.
 * @param nSample The sample size used in the GWAS whose summary statistics are being analyzed.
 * @param zScore A vector of Z-scores from GWAS summary statistics.
 * @param pValueThreshold Threshold for the p-value below which variants are considered for quality control.
 * @param propSVD Proportion of singular value decomposition (SVD) components retained in the analysis.
 * @param gcControl A boolean flag to apply genetic control corrections.
 * @param nIter The number of iterations to run the DENTIST algorithm.
 * @param gPvalueThreshold P-value threshold for grouping variants into significant and null categories.
 * @param ncpus The number of CPU cores to use for parallel processing.
 * @param seed Seed for random number generation, affecting the selection of variants for analysis.
 *
 * @return A List object containing:
 * - imputed_z: A vector of imputed Z-scores for each marker.
 * - rsq: A vector of R-squared values for each marker, indicating goodness of fit.
 * - corrected_z: A vector of adjusted Z-scores after error detection.
 * - iter_to_correct: An integer vector indicating the iteration in which each marker passed the quality control.
 * - is_problematic: A binary vector indicating whether each marker is considered problematic (1) or not (0).
 *
 * @note The function is designed for use in Rcpp and requires Armadillo for matrix operations and OpenMP for parallel processing.
 */

// [[Rcpp::export]]
List dentist_rcpp(const arma::mat& LDmat, uint nSample, const arma::vec& zScore,
                  double pValueThreshold, float propSVD, bool gcControl, int nIter,
                  double gPvalueThreshold, int ncpus, int seed) {
    
	// Set number of threads for parallel processing
	int nProcessors = omp_get_max_threads();
	if(ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);

	uint markerSize = zScore.size();
	// Initialization based on the seed input
	std::vector<size_t> randOrder = generateSetOfNumbers(markerSize, seed);
	std::vector<uint> idx, idx2, fullIdx(randOrder.begin(), randOrder.end());

	// Determining indices for partitioning
	for (uint i = 0; i < markerSize; ++i) {
		if (randOrder[i] > markerSize / 2) idx.push_back(i);
		else idx2.push_back(i);
	}

	std::vector<uint> groupingGWAS(markerSize, 0);
	for (uint i = 0; i < markerSize; ++i) {
		if (minusLogPvalueChisq2(std::pow(zScore(i),2)) > -log10(gPvalueThreshold)) {
			groupingGWAS[i] = 1;
		}
	}

	arma::vec imputedZ = arma::zeros<arma::vec>(markerSize);
	arma::vec rsq = arma::zeros<arma::vec>(markerSize);
	arma::vec zScore_e = arma::zeros<arma::vec>(markerSize);
	arma::ivec iterID = arma::zeros<arma::ivec>(markerSize);

	for (int t = 0; t < nIter; ++t) {
		std::vector<uint> idx2_QCed;
		std::vector<double> diff;
		std::vector<uint> grouping_tmp;

		// Perform iteration with current subsets
        // define rrr as the first time we run the oneIteration function during No.t iteraction, and initialize as 1
        int rrr = 1;
		oneIteration(LDmat, idx, idx2, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus, t, rrr);

		// Assess differences and grouping for thresholding
		diff.resize(idx2.size());
		grouping_tmp.resize(idx2.size());
		for (size_t i = 0; i < idx2.size(); ++i) {
			diff[i] = std::abs(zScore_e[idx2[i]]);
			grouping_tmp[i] = groupingGWAS[idx2[i]];
		}

		double threshold = getQuantile(diff, 0.995);
		double threshold1 = getQuantile2(diff, grouping_tmp, 0.995, false);
		double threshold0 = getQuantile2(diff, grouping_tmp, 0.995, true);

		if(threshold1 == 0) {
			threshold1 = threshold;
			threshold0 = threshold;
		}

		if(t > nIter - 2 ) {
			threshold0 = threshold;
			threshold1 = threshold;
		}

		// Apply threshold-based filtering for QC
		for (size_t i = 0; i < diff.size(); ++i) {
			if ((grouping_tmp[i] == 1 && diff[i] <= threshold1) ||
			    (grouping_tmp[i] == 0 && diff[i] <= threshold0)) {
				idx2_QCed.push_back(idx2[i]);
			}
		}

		// Perform another iteration with updated sets of indices (idx and idx2_QCed)
        // this is the second time we run the oneIteration function during No.t iteraction, so we switch the value of rrr from 1 to 2
        rrr = 2;
		oneIteration(LDmat, idx, idx2_QCed, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus, t, rrr);

		// Recalculate differences and groupings after the iteration
		diff.resize(fullIdx.size());
		grouping_tmp.resize(fullIdx.size());
		for(size_t i = 0; i < fullIdx.size(); ++i) {
			diff[i] = std::abs(zScore_e[fullIdx[i]]);
			grouping_tmp[i] = groupingGWAS[fullIdx[i]];
		}

		// Re-determine thresholds based on the recalculated differences and groupings
		threshold = getQuantile(diff, 0.995);
		threshold1 = getQuantile2(diff, grouping_tmp, 0.995, false);
		threshold0 = getQuantile2(diff, grouping_tmp, 0.995, true);

		if(threshold1 == 0) {
			threshold1 = threshold;
			threshold0 = threshold;
		}

		if(t > nIter - 2 ) {
			threshold0 = threshold;
			threshold1 = threshold;
		}

		// Adjust for genetic control and inflation factor if necessary
		std::vector<double> chisq(fullIdx.size());
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			chisq[i] = std::pow(zScore_e[fullIdx[i]], 2);
		}

		// Calculate the median chi-squared value as the inflation factor
		std::nth_element(chisq.begin(), chisq.begin() + chisq.size() / 2, chisq.end());
		double medianChisq = chisq[chisq.size() / 2];
		double inflationFactor = medianChisq / 0.46;

		std::vector<uint> fullIdx_tmp;
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			double currentDiffSquared = std::pow(diff[i], 2);

			if (gcControl) {
				// When gcControl is true, check if the variant passes the adjusted threshold
				if (!(diff[i] > threshold && minusLogPvalueChisq2(currentDiffSquared / inflationFactor) > -log10(pValueThreshold))) {
					fullIdx_tmp.push_back(fullIdx[i]);
				}
			} else {
				// When gcControl is false, simply check if the variant passes the basic threshold
				if (minusLogPvalueChisq2(currentDiffSquared) < -log10(pValueThreshold)) {
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
		randOrder = generateSetOfNumbers(fullIdx.size(), seed + t * seed); // Update seed for randomness
		idx.clear();
		idx2.clear();
		for (size_t i : fullIdx) {
			if (randOrder[i] > fullIdx.size() / 2) idx.push_back(i);
			else idx2.push_back(i);
		}
	}
	// Prepare and return results
	return List::create(Named("imputed_z") = imputedZ,
	                    Named("rsq") = rsq,
	                    Named("corrected_z") = zScore_e,
	                    Named("iter_to_correct") = iterID,
	                    Named("is_problematic") = wrap(groupingGWAS));
}