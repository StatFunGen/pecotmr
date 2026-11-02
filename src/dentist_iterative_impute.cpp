// cpp11 implementation of the DENTIST method (Chen et al.)
// https://github.com/Yves-CHEN/DENTIST/tree/master#Citations
//
// This code reproduces the original DENTIST binary's algorithm and logic,
// using Armadillo (LAPACK) for eigendecomposition and GCTA-style LD computation.
//
// Derived from DENTIST (https://github.com/Yves-CHEN/DENTIST), (C) 2018
// Wenhan Chen, Zhihong Zhu and Jian Yang, distributed under the GNU Lesser
// General Public License, version 3. See inst/COPYRIGHTS.
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

// Store one iteration's imputation for the variants in idx2: the imputed z,
// its R^2 (capped at 1, with the exceeding values collected for R to report
// once), and the outlier statistic (z - imputed) / sqrt(R_jj - rsq_j).
// `ldDiag` holds R_jj for each of idx2.
void storeIteration(const std::vector<size_t>& idx2, const vec& zScore,
                    const vec& zScore_imp, const vec& rsq_vec, const vec& ldDiag,
                    vec& imputedZ, vec& rsqList, vec& zScore_e,
                    std::vector<double>& rsqExceed) {
	for (size_t i = 0; i < idx2.size(); ++i) {
		imputedZ[idx2[i]] = zScore_imp(i);
		rsqList[idx2[i]] = rsq_vec(i);
		if (rsq_vec(i) >= 1) {
			rsqList[idx2[i]] = std::min(rsq_vec(i), 1.0);
			// Collect the exceeding value; returned to R (rsq_exceed) which
			// summarizes them into a single warning, instead of emitting a
			// per-occurrence warning from here.
			rsqExceed.push_back(rsq_vec(i));
		}
		size_t j = idx2[i];
		double denom_sq = ldDiag(i) - rsqList[j];
		if (denom_sq < 1e-8) denom_sq = 1e-8;
		zScore_e[j] = (zScore[j] - imputedZ[j]) / std::sqrt(denom_sq);
	}
}

// Perform one iteration of the DENTIST algorithm using Armadillo's eig_sym
// (LAPACK dsyevd) for eigendecomposition. Both Eigen and Armadillo return
// eigenvalues in ascending order, so the logic is identical.
void oneIteration(const mat& LD_mat, const std::vector<size_t>& idx, const std::vector<size_t>& idx2,
                  const vec& zScore, vec& imputedZ, vec& rsqList, vec& zScore_e,
                  size_t nSample, float probSVD, int numThreads, bool verbose,
                  std::vector<double>& rsqExceed) {
	if (verbose) {
		Rprintf("LD_mat: %lux%lu idx: %lu idx2: %lu\n",
		        (unsigned long)LD_mat.n_rows, (unsigned long)LD_mat.n_cols,
		        (unsigned long)idx.size(), (unsigned long)idx2.size());
	}

	int nProcessors = omp_get_max_threads();
	if (numThreads < nProcessors) nProcessors = numThreads;
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
	vec ldDiag(idx2.size());
	for (size_t i = 0; i < idx2.size(); ++i) ldDiag(i) = LD_mat(idx2[i], idx2[i]);
	storeIteration(idx2, zScore, zScore_imp, rsq_vec, ldDiag, imputedZ, rsqList,
	               zScore_e, rsqExceed);
}

// One iteration from a factor F of the LD (m x p, zero-mean columns,
// crossprod(F) / (m - 1) = R) instead of R itself. With the thin SVD
// F[, idx] = U D V', R[idx, idx] = V diag(d^2 / (m - 1)) V', so the dense
// route's top-K eigenpairs are the SVD's leading K components, and
//   R[idx2, idx] V_K diag(1 / lambda_K) = F[, idx2]' U_K diag(1 / d_K)
//   rsq = diag(beta V_K' R[idx, idx2]) = rowSums((F[, idx2]' U_K)^2) / (m - 1)
// -- the dense route's imputation exactly, without forming any of R. The
// eigenvalues beyond the SVD's are zero, below the dense route's 1e-4 rank
// cut, so both count the same rank.
void oneIterationFactor(const mat& F, const std::vector<size_t>& idx, const std::vector<size_t>& idx2,
                        const vec& zScore, vec& imputedZ, vec& rsqList, vec& zScore_e,
                        size_t nSample, float probSVD, int numThreads, bool verbose,
                        std::vector<double>& rsqExceed) {
	if (verbose) {
		Rprintf("factor: %lux%lu idx: %lu idx2: %lu\n",
		        (unsigned long)F.n_rows, (unsigned long)F.n_cols,
		        (unsigned long)idx.size(), (unsigned long)idx2.size());
	}

	int nProcessors = omp_get_max_threads();
	if (numThreads < nProcessors) nProcessors = numThreads;
	omp_set_num_threads(nProcessors);

	size_t K = std::min(static_cast<size_t>(idx.size()), nSample) * probSVD;

	if (idx2.size() > F.n_cols || idx.size() > F.n_cols)
		cpp11::stop("Inconsistent dimensions between the LD factor and idx2/idx in oneIterationFactor()");
	for (size_t i = 0; i < idx.size(); ++i)
		if (idx[i] >= zScore.size())
			cpp11::stop("Invalid index in idx: %d", (int)idx[i]);
	for (size_t i = 0; i < idx2.size(); ++i)
		if (idx2[i] >= zScore.size())
			cpp11::stop("Invalid index in idx2: %d", (int)idx2[i]);

	uvec aidx(idx.size()), aidx2(idx2.size());
	for (size_t i = 0; i < idx.size(); i++) aidx(i) = idx[i];
	for (size_t i = 0; i < idx2.size(); i++) aidx2(i) = idx2[i];

	double nm1 = static_cast<double>(F.n_rows) - 1.0;
	mat U, V;
	vec d;
	svd_econ(U, d, V, F.cols(aidx));
	vec eigval = square(d) / nm1;  // descending

	int nRank = 0;
	for (uword j = 0; j < eigval.n_elem; j++)
		if (eigval(j) >= 0.0001) nRank++;
	if (K > static_cast<size_t>(nRank)) K = nRank;

	if (verbose)
		Rprintf("Rank: %d, K: %lu\n", nRank, (unsigned long)K);

	if (K <= 1) {
		cpp11::warning("Rank of eigen matrix <= 1, skipping imputation for this partition");
		for (size_t i = 0; i < idx2.size(); ++i) {
			imputedZ[idx2[i]] = 0.0;
			rsqList[idx2[i]] = 0.0;
			zScore_e[idx2[i]] = 0.0;
		}
		return;
	}

	mat F2 = F.cols(aidx2);
	mat G = F2.t() * U.cols(0, K - 1);                // F[, idx2]' U_K
	mat beta = G.each_row() / d.subvec(0, K - 1).t();  // ... diag(1 / d_K)
	vec zScore_imp = beta * (V.cols(0, K - 1).t() * zScore(aidx));
	vec rsq_vec = sum(square(G), 1) / nm1;
	vec ldDiag = sum(square(F2), 0).t() / nm1;
	storeIteration(idx2, zScore, zScore_imp, rsq_vec, ldDiag, imputedZ, rsqList,
	               zScore_e, rsqExceed);
}

/**
 * @brief Executes DENTIST algorithm for quality control in GWAS summary data: the iterative imputation function.
 *
 * DENTIST (Detecting Errors iN analyses of summary staTISTics) identifies and removes problematic variants
 * in GWAS summary data by comparing observed GWAS statistics to predicted values based on linkage disequilibrium (LD)
 * information from a reference panel. It helps detect genotyping/imputation errors, allelic errors, and heterogeneity
 * between GWAS and LD reference samples, improving the reliability of subsequent analyses.
 *
 * The registered entry points take the reference LD either as a matrix (dentistIterativeImpute, ldMatR) or as
 * a factor F with R = F'F / (nrow(F) - 1) (dentistIterativeImputeFactor, factorR); both run dentistRun.
 *
 * @param one The per-iteration imputation: oneIteration on R, or oneIterationFactor on F.
 * @param nSample The LD reference panel's real sample size (NOT the GWAS size): it caps the SVD truncation
 * rank, K = min(|idx|, nSample) * propSVD.
 * @param zScoreR A vector of Z-scores from GWAS summary statistics.
 * @param pValueThreshold Threshold for the p-value below which variants are considered for quality control.
 * @param propSVD Proportion of singular value decomposition (SVD) components retained in the analysis.
 * @param gcControl A boolean flag to apply genetic control corrections.
 * @param nIter The number of iterations to run the DENTIST algorithm.
 * @param gPvalueThreshold P-value threshold for grouping variants into significant and null categories.
 * @param numThreads The number of CPU cores to use for parallel processing.
 * @param correctChenEtAlBug Whether to correct the original DENTIST bug.
 * @param verbose A boolean flag to enable verbose output for debugging.
 *
 * @return A List object containing:
 * - original_z: A vector of original Z-scores for each marker.
 * - imputed_z: A vector of imputed Z-scores for each marker.
 * - z_diff: A vector of outlier test z-scores
 * - rsq: A vector of R-squared values for each marker, indicating goodness of fit.
 * - iter_to_correct: An integer vector indicating the iteration in which each marker passed the quality control.
 */

// The DENTIST iterations around a per-iteration imputation `one`, which is
// the dense route (oneIteration on R) or the factor route (oneIterationFactor
// on a factor of R); everything else is shared.
template <class OneIteration>
cpp11::writable::list dentistRun(const OneIteration& one, const vec& zScore, int nSample,
                                 double pValueThreshold, double propSVD, bool gcControl, int nIter,
                                 double gPvalueThreshold, int numThreads, bool correctChenEtAlBug,
                                 bool verbose, cpp11::sexp seed) {
	if (verbose) {
		Rprintf("nSample: %d\n", nSample);
		Rprintf("zScore size: %lu\n", (unsigned long)zScore.size());
		Rprintf("pValueThreshold: %g\n", pValueThreshold);
		Rprintf("propSVD: %g\n", propSVD);
		Rprintf("gcControl: %d\n", gcControl);
		Rprintf("nIter: %d\n", nIter);
		Rprintf("gPvalueThreshold: %g\n", gPvalueThreshold);
		Rprintf("numThreads: %d\n", numThreads);
		Rprintf("correctChenEtAlBug: %d\n", correctChenEtAlBug);
	}

	// Set number of threads for parallel processing
	int nProcessors = omp_get_max_threads();
	if (numThreads < nProcessors) nProcessors = numThreads;
	omp_set_num_threads(nProcessors);

	size_t markerSize = zScore.size();
	// Resolve the RNG seed. NULL preserves the original DENTIST hard-coded seeds
	// (10 for the initial partition, 20000 + t*20000 per iteration) for exact
	// fidelity with the reference binary; a provided seed overrides both.
	bool userSeed = (seed != R_NilValue);
	unsigned int base_seed = userSeed ? cpp11::as_cpp<unsigned int>(seed) : 0u;
	unsigned int initSeed = userSeed ? base_seed : 10u;
	std::vector<size_t> randOrder = generateSetOfNumbers(markerSize, initSeed);
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
	// Rsq values that exceeded 1 (capped at 1.0); returned to R as rsq_exceed.
	std::vector<double> rsqExceed;

	for (int t = 0; t < nIter; ++t) {
		// Perform iteration with current subsets
		if (verbose) {
			Rprintf("\n=== Iteration %d ===\n", t);
			Rprintf("idx.size()=%lu idx2.size()=%lu fullIdx.size()=%lu\n",
			        (unsigned long)idx.size(), (unsigned long)idx2.size(), (unsigned long)fullIdx.size());
			Rprintf("Performing oneIteration()\n");
		}

		one(idx, idx2, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, numThreads, verbose, rsqExceed);

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
		        Thus if we correct the original DENTIST code, i.e., correctChenEtAlBug = TRUE,
		                we go through our function, getQuantile2, which doesn't have this issue
		                else, i.e., correctChenEtAlBug = TRUE, it goes through the original function getQuantile2_chen_et_al
		 */
		if (correctChenEtAlBug) {
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
		if (correctChenEtAlBug || nIter - 2 >= 0) {
			/*In the original DENTIST method, if t=0 (first iteration) and nIter is 1,
			   t is defined as a size_t (unassigned integer)
			   https://github.com/Yves-CHEN/DENTIST/blob/2fefddb1bbee19896a30bf56229603561ea1dba8/main/inversion.cpp#L628
			   and it will treat t (which is 0) no larger than nIter-2 (which is -1) which is wrong
			   Thus if we correct the original DENTIST code, i.e., correctChenEtAlBug = TRUE, or when nIter - 2 >=0,
			   it will compare t and nIter as we expect.
			   and if we want to keep the original DENTIST code, i.e., correctChenEtAlBug = TRUE, then it will skip this if condition for t > nIter - 2
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

		one(idx2_QCed, idx, zScore, imputedZ, rsq, zScore_e, nSample, propSVD, numThreads, verbose, rsqExceed);

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
		if (correctChenEtAlBug) {
			threshold1 = getQuantile2(diff, grouping_tmp, 0.995, false);
			threshold0 = getQuantile2(diff, grouping_tmp, 0.995, true);
		} else {
			threshold1 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
			std::transform(grouping_tmp.begin(), grouping_tmp.end(), grouping_tmp.begin(), [](size_t val) {
				return 1 - val;
			});
			threshold0 = getQuantile2_chen_et_al(diff, grouping_tmp, 0.995);
		}


		if (correctChenEtAlBug || nIter - 2 >= 0) {
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
					// inverted by the !operator. When correctChenEtAlBug=FALSE we must
					// use grouping_tmp to match original behavior; when TRUE we use the
					// un-mutated groupingGWAS.
					size_t grp = correctChenEtAlBug ? groupingGWAS[fullIdx[i]] : grouping_tmp[i];
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
		// NULL seed keeps the original per-iteration seed 20000 + t*20000; a
		// user seed derives a distinct per-iteration seed from base_seed.
		unsigned int iterSeed = userSeed
			? base_seed + static_cast<unsigned int>(t + 1) * 20000u
			: static_cast<unsigned int>(20000 + t * 20000);
		randOrder = generateSetOfNumbers(fullIdx.size(), iterSeed);
		idx.clear();
		idx2.clear();
		for (size_t i = 0; i < fullIdx.size(); ++i) {
			if (randOrder[i] > fullIdx.size() / 2) idx.push_back(fullIdx[i]);
			else idx2.push_back(fullIdx[i]);
		}
	}

	using namespace cpp11::literals;
	writable::list result({
		"originalZ"_nm = as_doubles(zScore),
		"imputedZ"_nm = as_doubles(imputedZ),
		"rsq"_nm = as_doubles(rsq),
		"zDiff"_nm = as_doubles(zScore_e),
		"iterToCorrect"_nm = as_integers(iterID),
		"rsqExceed"_nm = as_doubles(vec(rsqExceed))
	});

	return result;
}

[[cpp11::register]]
cpp11::writable::list dentistIterativeImpute(const doubles_matrix<>& ldMatR, int nSample, const doubles& zScoreR,
                              double pValueThreshold, double propSVD, bool gcControl, int nIter,
                              double gPvalueThreshold, int numThreads, bool correctChenEtAlBug,
                              bool verbose, cpp11::sexp seed = R_NilValue) {
	mat LD_mat = as_Mat(ldMatR);
	vec zScore = as_Col(zScoreR);
	if (verbose)
		Rprintf("LD_mat dimensions: %lu x %lu\n", (unsigned long)LD_mat.n_rows, (unsigned long)LD_mat.n_cols);
	auto one = [&LD_mat](const std::vector<size_t>& a, const std::vector<size_t>& b,
	                     const vec& z, vec& imp, vec& rsq, vec& ze, size_t n, float prop,
	                     int threads, bool verb, std::vector<double>& exceed) {
		oneIteration(LD_mat, a, b, z, imp, rsq, ze, n, prop, threads, verb, exceed);
	};
	return dentistRun(one, zScore, nSample, pValueThreshold, propSVD, gcControl, nIter,
	                  gPvalueThreshold, numThreads, correctChenEtAlBug, verbose, seed);
}

// DENTIST from a factor of the LD (see oneIterationFactor): the same
// iterations as dentistIterativeImpute(), never forming the LD matrix.
[[cpp11::register]]
cpp11::writable::list dentistIterativeImputeFactor(const doubles_matrix<>& factorR, int nSample, const doubles& zScoreR,
                              double pValueThreshold, double propSVD, bool gcControl, int nIter,
                              double gPvalueThreshold, int numThreads, bool correctChenEtAlBug,
                              bool verbose, cpp11::sexp seed = R_NilValue) {
	mat F = as_Mat(factorR);
	vec zScore = as_Col(zScoreR);
	if (F.n_cols != zScore.n_elem)
		cpp11::stop("The LD factor has %d columns but there are %d z-scores.",
		            (int)F.n_cols, (int)zScore.n_elem);
	if (verbose)
		Rprintf("LD factor dimensions: %lu x %lu\n", (unsigned long)F.n_rows, (unsigned long)F.n_cols);
	auto one = [&F](const std::vector<size_t>& a, const std::vector<size_t>& b,
	                const vec& z, vec& imp, vec& rsq, vec& ze, size_t n, float prop,
	                int threads, bool verb, std::vector<double>& exceed) {
		oneIterationFactor(F, a, b, z, imp, rsq, ze, n, prop, threads, verb, exceed);
	};
	return dentistRun(one, zScore, nSample, pValueThreshold, propSVD, gcControl, nIter,
	                  gPvalueThreshold, numThreads, correctChenEtAlBug, verbose, seed);
}
