// cpp11 implementation of the DENTIST method (Chen et al.)
// https://github.com/Yves-CHEN/DENTIST/tree/master#Citations
//
// This code reproduces the original DENTIST binary's algorithm and logic,
// using Armadillo (LAPACK) for eigendecomposition and GCTA-style LD computation.
#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <omp.h>
#include <algorithm>
#include <random>
#include <vector>
#include <unordered_set>

using namespace cpp11;
using namespace arma;

// DENTIST RNG: uses srand/rand + sort-by-random-values to produce a permutation.
// This faithfully reproduces the original DENTIST C++ binary's random partitioning.
//
// A simpler modern alternative would be:
//   std::vector<size_t> indexes(size);
//   std::iota(indexes.begin(), indexes.end(), 0);
//   std::mt19937 gen(seed);
//   std::shuffle(indexes.begin(), indexes.end(), gen);
//   return indexes;
// However, using the original RNG ensures exact reproducibility with the binary.
std::vector<size_t> generateSetOfNumbers(size_t size, unsigned int seed) {
	std::vector<int> numbers(size, 0);
	srand(seed);
	numbers[0] = rand();
	for (size_t index = 1; index < size; index++) {
		int tempNum;
		do {
			tempNum = rand();
			for (size_t index2 = 0; index2 < size; index2++)
				if (tempNum == numbers[index2]) tempNum = -1;
		} while (tempNum == -1);
		numbers[index] = tempNum;
	}
	// sort_indexes: return indices that would sort the vector
	std::vector<size_t> idx(size);
	std::iota(idx.begin(), idx.end(), 0);
	std::sort(idx.begin(), idx.end(), [&numbers](size_t i1, size_t i2) {
		return numbers[i1] < numbers[i2];
	});
	return idx;
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
// Use R::pchisq with lower.tail=FALSE to compute upper tail directly,
// avoiding catastrophic cancellation from 1.0 - CDF for large stats.
// This matches the original DENTIST binary's use of Boost complement().
double minusLogPvalueChisq2(double stat) {
	double p = Rf_pchisq(stat, 1.0, 0, 0);  // lower.tail=FALSE, log.p=FALSE
	return -log10(p);
}

// Perform one iteration of the DENTIST algorithm using Armadillo's eig_sym
// (LAPACK dsyevd) for eigendecomposition. Both Eigen and Armadillo return
// eigenvalues in ascending order, so the logic is identical.
void oneIteration(const mat& LD_mat, const std::vector<size_t>& idx, const std::vector<size_t>& idx2,
                  const vec& zScore, vec& imputedZ, vec& rsqList, vec& zScore_e,
                  size_t nSample, float probSVD, int ncpus, bool verbose) {
	if (verbose) {
		Rprintf("LD_mat: %lux%lu idx: %lu idx2: %lu\n",
		        (unsigned long)LD_mat.n_rows, (unsigned long)LD_mat.n_cols,
		        (unsigned long)idx.size(), (unsigned long)idx2.size());
	}

	int nProcessors = omp_get_max_threads();
	if (ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);

	size_t K = std::min(static_cast<size_t>(idx.size()), nSample) * probSVD;

	// Validate dimensions
	if (idx2.size() > LD_mat.n_rows || idx.size() > LD_mat.n_cols)
		cpp11::stop("Inconsistent dimensions between LD_mat and idx2/idx in oneIteration()");
	for (size_t i = 0; i < idx.size(); ++i)
		if (idx[i] >= zScore.size())
			cpp11::stop("Invalid index in idx: %d", (int)idx[i]);
	for (size_t i = 0; i < idx2.size(); ++i)
		if (idx2[i] >= zScore.size())
			cpp11::stop("Invalid index in idx2: %d", (int)idx2[i]);

	// Convert to arma::uvec for idiomatic submatrix extraction
	uvec aidx(idx.size()), aidx2(idx2.size());
	for (size_t i = 0; i < idx.size(); i++) aidx(i) = idx[i];
	for (size_t i = 0; i < idx2.size(); i++) aidx2(i) = idx2[i];

	// Extract submatrices and z-score subset
	mat LD_it = LD_mat(aidx2, aidx);
	mat VV = LD_mat(aidx, aidx);
	vec zScore_sub = zScore(aidx);

	// Eigendecomposition (ascending order, same as Eigen's SelfAdjointEigenSolver)
	vec eigval;
	mat eigvec;
	eig_sym(eigval, eigvec, VV);

	// Determine effective rank (eigenvalues >= 0.0001)
	int n_eig = eigval.n_elem;
	int nZeros = 0;
	for (int j = 0; j < n_eig; j++)
		if (eigval(j) < 0.0001) nZeros++;
	int nRank = n_eig - nZeros;
	if (K > static_cast<size_t>(nRank)) K = nRank;

	if (verbose)
		Rprintf("Rank: %d, Zeros: %d, K: %lu\n", nRank, nZeros, (unsigned long)K);

	if (K <= 1) {
		cpp11::warning("Rank of eigen matrix <= 1, skipping imputation for this partition");
		for (size_t i = 0; i < idx2.size(); ++i) {
			imputedZ[idx2[i]] = 0.0;
			rsqList[idx2[i]] = 0.0;
			zScore_e[idx2[i]] = 0.0;
		}
		return;
	}

	// Build ui (top K eigenvectors, largest first) and wi (inverse eigenvalues)
	mat ui(n_eig, K);
	vec wi(K);
	for (size_t m = 0; m < K; m++) {
		int j = n_eig - m - 1;  // from largest eigenvalue
		ui.col(m) = eigvec.col(j);
		wi(m) = 1.0 / eigval(j);
	}

	// Imputation: beta = LD_it * ui * diag(wi), then imputed_z = beta * ui' * z
	mat beta = LD_it * (ui.each_row() % wi.t());
	vec zScore_imp = beta * (ui.t() * zScore_sub);
	vec rsq_vec = diagvec(beta * (ui.t() * LD_it.t()));

	// Store results
	for (size_t i = 0; i < idx2.size(); ++i) {
		imputedZ[idx2[i]] = zScore_imp(i);
		rsqList[idx2[i]] = rsq_vec(i);
		if (rsq_vec(i) >= 1) {
			rsqList[idx2[i]] = std::min(rsq_vec(i), 1.0);
			cpp11::warning("Adjusted rsq value exceeding 1: %g", rsq_vec(i));
		}
		size_t j = idx2[i];
		double denom_sq = LD_mat(j, j) - rsqList[j];
		if (denom_sq < 1e-8) denom_sq = 1e-8;
		zScore_e[j] = (zScore[j] - imputedZ[j]) / std::sqrt(denom_sq);
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
 * @param LD_mat_r The linkage disequilibrium (LD) matrix from a reference panel.
 * @param nSample The sample size used in the GWAS whose summary statistics are being analyzed.
 * @param zScore_r A vector of Z-scores from GWAS summary statistics.
 * @param pValueThreshold Threshold for the p-value below which variants are considered for quality control.
 * @param propSVD Proportion of singular value decomposition (SVD) components retained in the analysis.
 * @param gcControl A boolean flag to apply genetic control corrections.
 * @param nIter The number of iterations to run the DENTIST algorithm.
 * @param gPvalueThreshold P-value threshold for grouping variants into significant and null categories.
 * @param ncpus The number of CPU cores to use for parallel processing.
 * @param correct_chen_et_al_bug Whether to correct the original DENTIST bug.
 * @param verbose A boolean flag to enable verbose output for debugging.
 *
 * @return A List object containing:
 * - original_z: A vector of original Z-scores for each marker.
 * - imputed_z: A vector of imputed Z-scores for each marker.
 * - z_diff: A vector of outlier test z-scores
 * - rsq: A vector of R-squared values for each marker, indicating goodness of fit.
 * - iter_to_correct: An integer vector indicating the iteration in which each marker passed the quality control.
 */

[[cpp11::register]]
cpp11::writable::list dentist_iterative_impute(const doubles_matrix<>& LD_mat_r, int nSample, const doubles& zScore_r,
                              double pValueThreshold, double propSVD, bool gcControl, int nIter,
                              double gPvalueThreshold, int ncpus, bool correct_chen_et_al_bug,
                              bool verbose) {
	mat LD_mat = as_Mat(LD_mat_r);
	vec zScore = as_Col(zScore_r);

	if (verbose) {
		Rprintf("LD_mat dimensions: %lu x %lu\n", (unsigned long)LD_mat.n_rows, (unsigned long)LD_mat.n_cols);
		Rprintf("nSample: %d\n", nSample);
		Rprintf("zScore size: %lu\n", (unsigned long)zScore.size());
		Rprintf("pValueThreshold: %g\n", pValueThreshold);
		Rprintf("propSVD: %g\n", propSVD);
		Rprintf("gcControl: %d\n", gcControl);
		Rprintf("nIter: %d\n", nIter);
		Rprintf("gPvalueThreshold: %g\n", gPvalueThreshold);
		Rprintf("ncpus: %d\n", ncpus);
		Rprintf("correct_chen_et_al_bug: %d\n", correct_chen_et_al_bug);
	}

	// Set number of threads for parallel processing
	int nProcessors = omp_get_max_threads();
	if (ncpus < nProcessors) nProcessors = ncpus;
	omp_set_num_threads(nProcessors);

	size_t markerSize = zScore.size();
	// Original DENTIST hardcodes seed=10 for initial partitioning
	std::vector<size_t> randOrder = generateSetOfNumbers(markerSize, 10);
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
		Rprintf("Indices partitioned\n");
	}

	std::vector<size_t> groupingGWAS(markerSize, 0);
	for (size_t i = 0; i < markerSize; ++i) {
		if (minusLogPvalueChisq2(std::pow(zScore(i), 2)) > -log10(gPvalueThreshold)) {
			groupingGWAS[i] = 1;
		}
	}

	if (verbose) {
		Rprintf("Grouping GWAS finished\n");
	}

	vec imputedZ = zeros<vec>(markerSize);
	vec rsq = zeros<vec>(markerSize);
	vec zScore_e = zeros<vec>(markerSize);
	Col<int> iterID = zeros<Col<int>>(markerSize);

	std::vector<double> diff(idx2.size());
	std::vector<size_t> grouping_tmp(idx2.size());

	for (int t = 0; t < nIter; ++t) {
		// Perform iteration with current subsets
		if (verbose) {
			Rprintf("\n=== Iteration %d ===\n", t);
			Rprintf("idx.size()=%lu idx2.size()=%lu fullIdx.size()=%lu\n",
			        (unsigned long)idx.size(), (unsigned long)idx2.size(), (unsigned long)fullIdx.size());
			Rprintf("Performing oneIteration()\n");
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
			Rprintf("Assessing differences and grouping for thresholding\n");
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
			   and if we want to keep the original DENTIST code, i.e., correct_chen_et_al_bug = TRUE, then it will skip this if condition for t > nIter - 2
			 */
			if (t > nIter - 2) {
				threshold0 = threshold;
				threshold1 = threshold;
			}
		}

		if (verbose) {
			Rprintf("Thresholds calculated: %g, %g, %g\n", threshold, threshold1, threshold0);
			Rprintf("Applying threshold-based filtering for QC\n");
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
			Rprintf("idx2_QCed.size()=%lu\n", (unsigned long)idx2_QCed.size());
			Rprintf("Performing oneIteration() with updated sets of indices\n");
		}

		oneIteration(LD_mat, idx2_QCed, idx, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, ncpus, verbose);

		if (verbose) {
			Rprintf("Recalculating differences and groupings after the iteration\n");
		}

		// Recalculate differences and groupings after the iteration
		diff.resize(fullIdx.size());
		grouping_tmp.resize(fullIdx.size());

		for (size_t i = 0; i < fullIdx.size(); ++i) {
			diff[i] = std::abs(zScore_e[fullIdx[i]]);
			grouping_tmp[i] = groupingGWAS[fullIdx[i]];
		}

		if (verbose) {
			Rprintf("Re-determining thresholds based on the recalculated differences and groupings\n");
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
			Rprintf("Phase2 thresholds: %g, %g, %g\n", threshold, threshold1, threshold0);
			Rprintf("Adjusting for genetic control and inflation factor if necessary\n");
		}

		// Adjust for genetic control and inflation factor if necessary
		std::vector<double> chisq(fullIdx.size());
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			chisq[i] = std::pow(zScore_e[fullIdx[i]], 2);
		}

		// Original DENTIST does not check chisq.size(); it just computes the median.
		// We only need to guard against empty vectors to avoid undefined behavior.
		if (chisq.empty()) {
			if (verbose) {
				Rprintf("chisq is empty, breaking out of iteration loop.\n");
			}
			break;
		}

		// Calculate the median chi-squared value as the inflation factor
		std::nth_element(chisq.begin(), chisq.begin() + chisq.size() / 2, chisq.end());
		double medianChisq = chisq[chisq.size() / 2];
		double inflationFactor = medianChisq / 0.46;

		std::vector<size_t> fullIdx_tmp;
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			// Use diff[i]*diff[i] instead of chisq[i] because nth_element above
			// scrambled the chisq array. The binary uses diff[i]*diff[i] here.
			double chisq_i = diff[i] * diff[i];
			if (gcControl) {
				// When gcControl is true, check if the variant passes the adjusted threshold
				if (!(diff[i] > threshold && minusLogPvalueChisq2(chisq_i / inflationFactor) > -log10(pValueThreshold))) {
					fullIdx_tmp.push_back(fullIdx[i]);
				}
			} else {
				// When gcControl is false, simply check if the variant passes the basic threshold
				if (minusLogPvalueChisq2(chisq_i) < -log10(pValueThreshold)) {
					// In original DENTIST, grouping_tmp[i] is used here, which has been
					// inverted by the !operator. When correct_chen_et_al_bug=FALSE we must
					// use grouping_tmp to match original behavior; when TRUE we use the
					// un-mutated groupingGWAS.
					size_t grp = correct_chen_et_al_bug ? groupingGWAS[fullIdx[i]] : grouping_tmp[i];
					if ((grp == 1 && diff[i] <= threshold1) ||
					    (grp == 0 && diff[i] <= threshold0)) {
						fullIdx_tmp.push_back(fullIdx[i]);
						iterID[fullIdx[i]]++;
					}
				}
			}
		}

		// Update the indices for the next iteration based on filtering criteria
		fullIdx = fullIdx_tmp;
		if (verbose) {
			Rprintf("Iter %d: fullIdx=%lu threshold=%g threshold1=%g threshold0=%g\n",
			        t, (unsigned long)fullIdx.size(), threshold, threshold1, threshold0);
		}
		// Early exit if all variants were filtered out
		if (fullIdx.empty()) {
			if (verbose) {
				Rprintf("All variants filtered out at iteration %d, stopping early.\n", t);
			}
			break;
		}
		// Original DENTIST uses seed = 20000 + t*20000 for subsequent iterations
		randOrder = generateSetOfNumbers(fullIdx.size(), 20000 + t * 20000);
		idx.clear();
		idx2.clear();
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			if (randOrder[i] > fullIdx.size() / 2) idx.push_back(fullIdx[i]);
			else idx2.push_back(fullIdx[i]);
		}
	}

	using namespace cpp11::literals;
	writable::list result({
		"original_z"_nm = as_doubles(zScore),
		"imputed_z"_nm = as_doubles(imputedZ),
		"rsq"_nm = as_doubles(rsq),
		"z_diff"_nm = as_doubles(zScore_e),
		"iter_to_correct"_nm = as_integers(iterID)
	});

	return result;
}
