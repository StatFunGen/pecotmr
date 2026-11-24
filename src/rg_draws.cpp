// rg_draws.cpp -- the per-draw statistics of estimateLocalRg()'s simulated
// p-values and intervals for a pair given the rest of its group.
//
// R draws every random number, in the order it always has, and hands the
// draws here whole; each function reads one draw at a time, as the R loops
// over draws did:
//   partial  .rgPartialConditionalP()'s draws of the components of
//            condition g and the rest (C) under LAVA's null: each draw's
//            conditional p-value of no partial covariance of g and f, NA
//            where the rest's covariance in the draw cannot be solved;
//   mixed    .rgMixedDrawChunk()'s draws of a group's H, or
//            .rgConditionalDrawQuantiles()'s Wishart draws of its Omega on
//            one reference, one column of its cells per draw: each draw's
//            partial covariance and rg of the pair given the rest
//            (.rgPartial(), .rgCorrelation()) and the two conditions' r2
//            given every other (.rgExplained()).
// A solve is R's solve(): the LU decomposition of LAPACK's dgesv, refused
// where it is exactly singular or its reciprocal condition number in the
// 1-norm is below machine epsilon. A positive definite matrix is one whose
// smallest eigenvalue, of its lower triangle as R's eigen() reads it, is
// above sqrt(epsilon) times the largest in size (.rgPositiveDefinite()).

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <cmath>
#include <limits>
#include <vector>

#define R_NO_REMAP_RMATH
#include <Rmath.h>

using namespace cpp11;
using namespace cpp11::literals;
using namespace arma;

namespace {

const double kEps = std::numeric_limits<double>::epsilon();
const double kNa = NA_REAL;

// R's 1-based indices as 0-based.
uvec zeroBased(const integers& x) {
  uvec out(x.size());
  for (R_xlen_t k = 0; k < x.size(); ++k) out[k] = x[k] - 1;
  return out;
}

// x solving A x = b as R's solve() does (La_solve()): dgetrf and dgetrs on a
// copy of A, then dgecon's reciprocal condition number in the 1-norm from
// A's norm; false where R's solve() stops, U having a zero pivot or the
// condition number being below epsilon.
bool rSolve(const mat& A, const mat& b, mat& x) {
  blas_int n = static_cast<blas_int>(A.n_rows);
  blas_int nrhs = static_cast<blas_int>(b.n_cols);
  blas_int info = 0;
  mat lu = A;
  x = b;
  podarray<blas_int> pivots(A.n_rows);
  lapack::getrf(&n, &n, lu.memptr(), &n, pivots.memptr(), &info);
  if (info != 0) return false;
  char trans = 'N';
  lapack::getrs(&trans, &n, &nrhs, lu.memptr(), &n, pivots.memptr(),
                x.memptr(), &n, &info);
  if (info != 0) return false;
  char one = '1';
  podarray<double> work(4 * A.n_rows);
  podarray<blas_int> iwork(A.n_rows);
  // dlange reads A without writing it.
  double anorm = lapack::lange(&one, &n, &n, const_cast<double*>(A.memptr()),
                               &n, work.memptr());
  double rcond = 0.0;
  lapack::gecon(&one, &n, lu.memptr(), &n, &anorm, &rcond, work.memptr(),
                iwork.memptr(), &info);
  return rcond >= kEps;
}

// A solve R's code makes without a handler stops as R's solve() would.
mat solveOrStop(const mat& A, const mat& b) {
  mat x;
  if (!rSolve(A, b, x)) {
    cpp11::stop("system is computationally singular");
  }
  return x;
}

// Whether a symmetric matrix is positive definite, to its scale: its
// eigenvalues read from its lower triangle, the smallest above
// sqrt(epsilon) times the largest in size (at least 1e-300).
bool positiveDefinite(const mat& x) {
  vec values;
  if (!eig_sym(values, symmatl(x))) {
    cpp11::stop("eigen() failed on a draw's covariance");
  }
  double largest = std::max(arma::max(arma::abs(values)), 1e-300);
  return values.min() > std::sqrt(kEps) * largest;
}

// The genetic correlation of a covariance and two heritabilities, NA unless
// both are positive (.rgCorrelation()).
double correlation(double covariance, double h1, double h2) {
  if (std::isnan(h1) || std::isnan(h2) || std::isnan(covariance) ||
      h1 <= 0.0 || h2 <= 0.0) {
    return kNa;
  }
  return covariance / std::sqrt(h1 * h2);
}

struct Partial {
  double cov, var1, var2;
};

// The covariance of conditions i and j and their two variances given the
// conditions `rest`, NA where rest's covariance is not positive definite:
// H_ij,ij - H_rest,ij' H_rest^-1 H_rest,ij (.rgPartial()).
Partial partialOf(const mat& H, uword i, uword j, const uvec& rest) {
  if (rest.n_elem == 0) return {H(i, j), H(i, i), H(j, j)};
  const mat Hrr = H.submat(rest, rest);
  if (!positiveDefinite(Hrr)) return {kNa, kNa, kNa};
  const uvec ij = {i, j};
  const mat across = H.submat(rest, ij);
  const mat left = H.submat(ij, ij) - across.t() * solveOrStop(Hrr, across);
  return {left(0, 1), left(0, 0), left(1, 1)};
}

// The share of condition k's variance the others explain together, H_ok'
// H_oo^-1 H_ok / H_kk; NA where H_kk is not positive or H_oo is not
// positive definite (.rgExplained()).
double explained(const mat& H, uword k) {
  uvec others(H.n_rows - 1);
  for (uword l = 0, at = 0; l < H.n_rows; ++l) {
    if (l != k) others[at++] = l;
  }
  const mat Hoo = H.submat(others, others);
  if (H(k, k) <= 0.0 || !positiveDefinite(Hoo)) return kNa;
  const uvec at = {k};
  const mat across = H.submat(others, at);
  return accu(across % solveOrStop(Hoo, across)) / H(k, k);
}

}  // namespace

// Each draw's conditional p-value of no partial covariance of conditions g
// and f given the rest, g = C[1] and rest = C[-1] (.rgPartialConditionalP()).
// Draw k's components of C are D_C = top_C + noise_k, P rows, with S = D_C'
// D_C + gram_k their cross-products over all K rows and D_M = D_C' top; its
// Omega_C = S / K - Sigma_CC gives a = Omega_rest^-1 Omega_rest,g, NA where
// that cannot be solved, and w = (1, -a). f's components given C's are
// normal, mean M_f + (D_C - M_C) B and variance tau^2 each, so the partial
// covariance w' D_f / K - (Sigma_gf - a' Sigma_rest,f) is normal with mean
// (w' D_M,f + w' (S - D_M,C) B) / K - (Sigma_gf - a' Sigma_rest,f) and
// standard deviation sqrt(tau^2 w' S w) / K, and the draw's p-value is its
// two-sided tail beyond the observed partial covariance.
[[cpp11::register]]
cpp11::writable::doubles rgPartialDrawPCpp(const doubles_matrix<>& topR,
                                           const doubles_matrix<>& noiseR,
                                           const doubles& gramR,
                                           const doubles_matrix<>& sigmaR,
                                           double K, const integers& CR,
                                           int f, const doubles& BR,
                                           double tau2, double observedR) {
  const mat top = as_Mat(topR), noise = as_Mat(noiseR);
  const mat Sigma = as_Mat(sigmaR);
  const uvec C = zeroBased(CR);
  const uword P = top.n_cols, q = C.n_elem, fi = f - 1, g = C[0];
  const uvec rest = C.tail(q - 1);
  const uword nDraws = noise.n_rows / P;
  const cube gram(REAL(static_cast<SEXP>(gramR)), q, q, nDraws);
  const vec B(REAL(static_cast<SEXP>(BR)), q);
  const mat topC = top.cols(C), sigmaCC = Sigma.submat(C, C);
  const vec sigmaRestF = Sigma.submat(rest, uvec{fi});
  const double observed = std::abs(observedR);
  vec out(nDraws);
  for (uword k = 0; k < nDraws; ++k) {
    const mat DC = topC + noise.rows(k * P, k * P + P - 1);
    const mat S = DC.t() * DC + gram.slice(k);
    const mat DM = DC.t() * top;
    const mat omegaC = S / K - sigmaCC;
    mat a;
    if (!rSolve(omegaC.submat(1, 1, q - 1, q - 1),
                omegaC.submat(1, 0, q - 1, 0), a)) {
      out[k] = kNa;
      continue;
    }
    const vec w = join_cols(vec{1.0}, -vec(a));
    const double shared = Sigma(g, fi) - accu(vec(a) % sigmaRestF);
    const double mu =
        (accu(w % DM.col(fi)) + accu(w % ((S - DM.cols(C)) * B))) / K - shared;
    const double sd = std::sqrt(tau2 * accu(w % (S * w))) / K;
    out[k] = Rf_pnorm5(-observed, mu, sd, 1, 0) +
             Rf_pnorm5(observed, mu, sd, 0, 0);
  }
  return as_doubles(out);
}

// Each draw's pair statistics from its cells, one column of the P x P H
// (column-major) per draw: the partial covariance and rg of conditions i
// and j given `rest` and the r2 of each given every other condition
// (.rgMixedDrawChunk(), .rgConditionalDrawQuantiles()).
[[cpp11::register]]
cpp11::writable::list rgMixedDrawPairsCpp(const doubles_matrix<>& cellsR,
                                          int P, int i, int j,
                                          const integers& restR) {
  const mat cells = as_Mat(cellsR);
  const uvec rest = zeroBased(restR);
  const uword n = cells.n_cols, ii = i - 1, jj = j - 1;
  vec cov(n), rg(n), r2i(n), r2j(n);
  for (uword d = 0; d < n; ++d) {
    const mat H = reshape(cells.col(d), P, P);
    const Partial partial = partialOf(H, ii, jj, rest);
    cov[d] = partial.cov;
    rg[d] = correlation(partial.cov, partial.var1, partial.var2);
    r2i[d] = explained(H, ii);
    r2j[d] = explained(H, jj);
  }
  return writable::list({"cov"_nm = as_doubles(cov), "rg"_nm = as_doubles(rg),
                         "r2i"_nm = as_doubles(r2i),
                         "r2j"_nm = as_doubles(r2j)});
}
