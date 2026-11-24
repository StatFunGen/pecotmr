// eigen_likelihood.cpp -- the Gaussian likelihoods of rotated statistics that
// criterion = "likelihood" fits, and their maximum.
//
// Every model reads statistics rotated onto an LD reference's eigenvectors
// as Gaussian terms, each given those before it:
//   variance   one condition's rotated statistics b ~ N(0, v), v = offset +
//              G theta floored: the heritabilities (one per annotation) and
//              the intercept, either held by moving it into the offset;
//   tied       b ~ N(0, s v), v = offset + g h2 floored, the noise tied to
//              1 - h2 and s a scale on the phenotype's variance where fitted;
//   nullSpace  the tied model over all n samples: the window's components
//              a_k with variance d (1 + h2 c_k) and the mean square of the n -
//              r directions outside them with d (1 - h2), weighted by their
//              number;
//   cross      a pair's second statistics given the first's at their
//              covariance v12 = offset + G theta: mean v12 / v11 b1, variance
//              v22 - v12^2 / v11 floored;
//   reduced    a pair on two references reduced to independent pairs (a_i,
//              b_i) of correlation g s_i: b_i given a_i, mean g s_i a_i and
//              variance 1 - (g s_i)^2 floored;
//   joint      a group on one reference: each condition's statistics given
//              those before it, through the LDL' decomposition of each
//              direction's covariance, with each conditional variance
//              floored; the genetic correlations from canonical partial
//              correlations and the pairs' shared noise fitted after them or
//              held. On several references (`mixed`), the block LDL'
//              decomposition of every condition's statistics at once, each
//              conditional covariance's eigenvalues floored.
// Each gives its log-likelihood and, where asked, its score and expected
// information: for terms e ~ N(0, D), each weighted by w,
//   score = sum(w (e^2 / D - 1) dD / (2 D) - e de / D),
//   info  = sum(w dD dD' / (2 D^2) + de de' / D),
// with dD and de the variances' and residuals' derivatives (a variance at
// its floor has none), and on several references the same for each block,
// r ~ N(0, V): -dr' V^-1 r + (r' V^-1 dV V^-1 r - tr(V^-1 dV)) / 2 and dr_j'
// V^-1 dr_l + tr(V^-1 dV_j V^-1 dV_l) / 2. A group's derivatives are carried
// through the decomposition; only its P x P covariances' are differenced.
//
// The maximum is dlib's box-constrained minimizer, find_min_box_constrained(),
// of the negative log-likelihood along a projected Newton direction
// (Bertsekas 1982) with the expected information as its Hessian, so
// projected Fisher scoring: a parameter within eps of a bound its gradient
// pushes past steps onto the bound, eps the smaller of 1e-3 of its range and
// the size of the diagonally scaled projected gradient step; the rest take
// the information's Newton step on their block; and the direction is cut
// back to the box, so that dlib's line search reads the slope of the path it
// walks. Each iteration is one of dlib's, its line search from the whole
// step. A run stops where the next step's predicted rise is below 1e-14 of
// the log-likelihood's size (1e-10 where only the value is read: a null's or
// a profile's) or an iteration does not raise it, and it has not converged
// where it reaches its limit first.
//
// Bertsekas DP (1982). Projected Newton methods for optimization problems
// with simple constraints. SIAM J Control Optim, 20:221-246.
//
// The variance likelihood follows HDL's llfun() and the cross likelihood
// HDL.L()'s step 2 (https://github.com/zhenin/HDL, Copyright (C) the HDL
// authors, GPL-3), and the nullSpace likelihood FastHer's fasther.R
// (https://github.com/gulsvi/FastHer, Copyright (c) 2025 Gulnara R.
// Svishcheva, MIT). Their license notices are in inst/COPYRIGHTS.

#include "dlib_compat.h"

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace cpp11;
using namespace cpp11::literals;
using namespace arma;

namespace {

typedef dlib::matrix<double, 0, 1> column;

const double kInf = std::numeric_limits<double>::infinity();

// ---- Reading the data -------------------------------------------------------

// Element `name` of `data`, R's NULL where it has none.
SEXP element(const list& data, const char* name) {
  if (Rf_isNull(data.names())) return R_NilValue;
  strings names(data.names());
  for (R_xlen_t i = 0; i < names.size(); ++i) {
    if (std::string(names[i]) == name) return data[i];
  }
  return R_NilValue;
}

// A numeric vector or matrix, whatever R's storage mode.
vec numbers(SEXP x) {
  sexp real(safe[Rf_coerceVector](x, REALSXP));
  return vec(REAL(real), Rf_xlength(real));
}

mat numberMatrix(SEXP x) {
  sexp real(safe[Rf_coerceVector](x, REALSXP));
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (Rf_isNull(dim)) return mat(REAL(real), Rf_xlength(real), 1);
  return mat(REAL(real), INTEGER(dim)[0], INTEGER(dim)[1]);
}

SEXP needed(const list& data, const char* name) {
  SEXP x = element(data, name);
  if (Rf_isNull(x)) cpp11::stop("The likelihood's data has no '%s'.", name);
  return x;
}

vec neededVec(const list& data, const char* name) {
  return numbers(needed(data, name));
}

mat neededMat(const list& data, const char* name) {
  return numberMatrix(needed(data, name));
}

double neededNumber(const list& data, const char* name) {
  return neededVec(data, name)[0];
}

bool flag(const list& data, const char* name) {
  SEXP x = element(data, name);
  return !Rf_isNull(x) && numbers(x)[0] != 0.0;
}

// A term-wise vector given as one value or one per term.
vec spread(const vec& x, uword n) {
  return x.n_elem == 1 ? vec(n, fill::value(x[0])) : x;
}

vec indicator(const uvec& b) { return conv_to<vec>::from(b); }

// x kept within [lower, upper], element by element.
vec boxed(const vec& x, const vec& lower, const vec& upper) {
  return arma::min(arma::max(x, lower), upper);
}

bool same(const vec& a, const vec& b) {
  return a.n_elem == b.n_elem && std::equal(a.begin(), a.end(), b.begin());
}

// ---- The terms' log-likelihood, score and information -----------------------

struct Evaluation {
  double value;
  vec score;
  mat information;
};

// The log-likelihood of terms e ~ N(0, D), each standing for w alike.
double termsValue(const vec& D, const vec& e2, const vec& w) {
  return -accu(w % (log(D) + e2 / D)) / 2.0;
}

// Their score and expected information through their variances, dD the
// variances' derivatives (terms x parameters).
void varianceScoring(Evaluation& out, const vec& D, const vec& e2,
                     const vec& w, const mat& dD) {
  out.score = dD.t() * (w % (e2 / D - 1.0) / (2.0 * D));
  mat scaled = dD.each_col() % (sqrt(w) / D);
  out.information = scaled.t() * scaled / 2.0;
}

// What the terms' means add where they move, de the residuals' derivatives.
void meanScoring(Evaluation& out, const vec& D, const vec& e, const mat& de) {
  out.score -= de.t() * (e / D);
  mat scaled = de.each_col() / sqrt(D);
  out.information += scaled.t() * scaled;
}

// An evaluation's score and information in the parameters `free` alone.
void restrict(Evaluation& out, const uvec& free) {
  out.score = vec(out.score.elem(free));
  out.information = mat(out.information.submat(free, free));
}

// ---- The models -------------------------------------------------------------

class Likelihood {
 public:
  virtual ~Likelihood() = default;
  // The log-likelihood at every parameter theta and, where `scoring` asks,
  // its score and expected information in the parameters `free`.
  virtual Evaluation at(const vec& theta, bool scoring,
                        const uvec& free) const = 0;
};

class VarianceLikelihood : public Likelihood {
 public:
  explicit VarianceLikelihood(const list& data)
      : G_(neededMat(data, "G")),
        b2_(neededVec(data, "b2")),
        offset_(spread(neededVec(data, "offset"), b2_.n_elem)),
        floor_(neededNumber(data, "floor")) {}

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const vec raw = offset_ + G_ * theta;
    const vec v = clamp(raw, floor_, kInf);
    const vec w(v.n_elem, fill::ones);
    Evaluation out{termsValue(v, b2_, w), vec(), mat()};
    if (scoring) {
      varianceScoring(out, v, b2_, w, G_.each_col() % indicator(raw > floor_));
      restrict(out, free);
    }
    return out;
  }

 private:
  mat G_;
  vec b2_, offset_;
  double floor_;
};

class TiedLikelihood : public Likelihood {
 public:
  explicit TiedLikelihood(const list& data)
      : offset_(neededVec(data, "offset")),
        g_(neededVec(data, "g")),
        b2_(neededVec(data, "b2")),
        floor_(neededNumber(data, "floor")),
        scaled_(flag(data, "scaled")) {}

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const vec raw = offset_ + g_ * theta[0];
    const vec v = clamp(raw, floor_, kInf);
    const double s = scaled_ ? theta[1] : 1.0;
    const vec D = s * v;
    const vec w(D.n_elem, fill::ones);
    Evaluation out{termsValue(D, b2_, w), vec(), mat()};
    if (scoring) {
      mat dD(D.n_elem, scaled_ ? 2 : 1);
      dD.col(0) = s * g_ % indicator(raw > floor_);
      if (scaled_) dD.col(1) = v;
      varianceScoring(out, D, b2_, w, dD);
      restrict(out, free);
    }
    return out;
  }

 private:
  vec offset_, g_, b2_;
  double floor_;
  bool scaled_;
};

class NullSpaceLikelihood : public Likelihood {
 public:
  explicit NullSpaceLikelihood(const list& data)
      : a_(neededVec(data, "a")),
        ck_(neededVec(data, "ck")),
        n_(neededNumber(data, "n")),
        scaled_(flag(data, "scaled")) {
    const double rest = n_ - a_.n_elem;
    e2_ = join_cols(a_, vec{(n_ - accu(a_)) / rest});
    w_ = join_cols(vec(a_.n_elem, fill::ones), vec{rest});
  }

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const double h2 = theta[0];
    const double d = scaled_ ? theta[1] : 1.0;
    const vec shape = join_cols(1.0 + h2 * ck_, vec{1.0 - h2});
    const vec D = d * shape;
    Evaluation out{termsValue(D, e2_, w_), vec(), mat()};
    if (scoring) {
      mat dD(D.n_elem, scaled_ ? 2 : 1);
      dD.col(0) = d * join_cols(ck_, vec{-1.0});
      if (scaled_) dD.col(1) = shape;
      varianceScoring(out, D, e2_, w_, dD);
      restrict(out, free);
    }
    return out;
  }

 private:
  vec a_, ck_, e2_, w_;
  double n_;
  bool scaled_;
};

class CrossLikelihood : public Likelihood {
 public:
  explicit CrossLikelihood(const list& data)
      : G_(neededMat(data, "G")),
        b1_(neededVec(data, "b1")),
        b2_(neededVec(data, "b2")),
        v11_(neededVec(data, "v11")),
        v22_(neededVec(data, "v22")),
        offset_(spread(neededVec(data, "offset"), b1_.n_elem)),
        floor_(neededNumber(data, "floor")) {}

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const vec v12 = offset_ + G_ * theta;
    const vec e = b2_ - v12 / v11_ % b1_;
    const vec raw = v22_ - square(v12) / v11_;
    const vec cv = clamp(raw, floor_, kInf);
    const vec w(cv.n_elem, fill::ones);
    Evaluation out{termsValue(cv, square(e), w), vec(), mat()};
    if (scoring) {
      const vec slope = 2.0 * v12 / v11_ % indicator(raw > floor_);
      varianceScoring(out, cv, square(e), w, -(G_.each_col() % slope));
      meanScoring(out, cv, e, -(G_.each_col() % (b1_ / v11_)));
      restrict(out, free);
    }
    return out;
  }

 private:
  mat G_;
  vec b1_, b2_, v11_, v22_, offset_;
  double floor_;
};

class ReducedLikelihood : public Likelihood {
 public:
  explicit ReducedLikelihood(const list& data)
      : s_(neededVec(data, "s")),
        a_(neededVec(data, "a")),
        b_(neededVec(data, "b")),
        floor_(neededNumber(data, "floor")) {}

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const double g = theta[0];
    const vec gs = g * s_;
    const vec e = b_ - gs % a_;
    const vec raw = 1.0 - square(gs);
    const vec cv = clamp(raw, floor_, kInf);
    const vec w(cv.n_elem, fill::ones);
    Evaluation out{termsValue(cv, square(e), w), vec(), mat()};
    if (scoring) {
      const vec dD = -2.0 * g * square(s_) % indicator(raw > floor_);
      varianceScoring(out, cv, square(e), w, mat(dD));
      meanScoring(out, cv, e, mat(-s_ % a_));
      restrict(out, free);
    }
    return out;
  }

 private:
  vec s_, a_, b_;
  double floor_;
};

// ---- A group's covariances --------------------------------------------------

// The lower-triangular factor of a correlation matrix from its canonical
// partial correlations, row by row below the diagonal.
mat cpcFactor(const vec& z, uword P) {
  mat F(P, P, fill::zeros);
  for (uword i = 0; i < P; ++i) {
    double length = 1.0;
    for (uword j = 0; j < i; ++j) {
      const double zij = z[i * (i - 1) / 2 + j];
      F(i, j) = zij * length;
      length *= std::sqrt(1.0 - zij * zij);
    }
    F(i, i) = length;
  }
  return F;
}

// The symmetric matrix of the pairs' values, row by row below the diagonal,
// 0 on it.
mat pairMatrix(const vec& values, uword P) {
  mat M(P, P, fill::zeros);
  for (uword i = 1; i < P; ++i) {
    for (uword j = 0; j < i; ++j) {
      M(i, j) = values[i * (i - 1) / 2 + j];
      M(j, i) = M(i, j);
    }
  }
  return M;
}

struct Covariances {
  mat G, C;
};

// What every group reads: its genetic covariances G at the partial
// correlations theta[1:nPairs] and the noise each pair shares C on the
// statistics' scale, fitted after them or held (less the covariance for tied
// noise), with their derivatives differenced a little either way within the
// parameters' bounds.
class GroupLikelihood : public Likelihood {
 public:
  explicit GroupLikelihood(const list& data) {
    const list joint(needed(data, "joint"));
    P_ = static_cast<uword>(neededNumber(joint, "P"));
    nPairs_ = P_ * (P_ - 1) / 2;
    crossScale_ = neededMat(joint, "crossScale");
    h2_ = neededVec(joint, "h2");
    tied_ = neededMat(joint, "tied");
    floor_ = neededNumber(joint, "floor");
    SEXP held = element(data, "held");
    held_ = Rf_isNull(held) ? mat(P_, P_, fill::zeros) : numberMatrix(held);
    SEXP lower = element(data, "lower");
    SEXP upper = element(data, "upper");
    lower_ = Rf_isNull(lower) ? vec() : numbers(lower);
    upper_ = Rf_isNull(upper) ? vec() : numbers(upper);
  }

 protected:
  Covariances covariances(const vec& theta) const {
    const mat F = cpcFactor(theta.head(nPairs_), P_);
    const mat G = (F * F.t()) % sqrt(h2_ * h2_.t());
    const mat shared =
        theta.n_elem > nPairs_
            ? pairMatrix(theta.tail(theta.n_elem - nPairs_), P_)
            : mat(held_ - tied_ % G);
    return {G, shared % crossScale_};
  }

  void covarianceSteps(const vec& theta, const uvec& free, cube& dG,
                       cube& dC) const {
    dG.set_size(P_, P_, free.n_elem);
    dC.set_size(P_, P_, free.n_elem);
    for (uword t = 0; t < free.n_elem; ++t) {
      const uword j = free[t];
      const double h = 1e-5 * std::max(1.0, std::abs(theta[j]));
      const double top = j < upper_.n_elem ? upper_[j] : kInf;
      const double bottom = j < lower_.n_elem ? lower_[j] : -kInf;
      const double up = std::min(theta[j] + h, top);
      const double down = std::max(theta[j] - h, bottom);
      vec hi = theta, lo = theta;
      hi[j] = up;
      lo[j] = down;
      const Covariances a = covariances(hi), b = covariances(lo);
      dG.slice(t) = (a.G - b.G) / (up - down);
      dC.slice(t) = (a.C - b.C) / (up - down);
    }
  }

  uword P_, nPairs_;
  mat crossScale_, tied_, held_;
  vec h2_, lower_, upper_;
  double floor_;
};

// A group on one reference: every block diagonal, so the decomposition runs
// direction by direction, a loading a vector over the directions and its
// derivatives a matrix over them and the parameters.
class JointLikelihood : public GroupLikelihood {
 public:
  explicit JointLikelihood(const list& data) : GroupLikelihood(data) {
    const list joint(needed(data, "joint"));
    u_ = neededMat(joint, "u");
    v_ = neededMat(joint, "v");
    lam_ = neededVec(joint, "lam");
    lam2M_ = neededVec(joint, "lam2M");
  }

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const uword K = u_.n_rows, P = P_, q = scoring ? free.n_elem : 0;
    const Covariances cov = covariances(theta);
    cube dG, dC;
    if (scoring) covarianceSteps(theta, free, dG, dC);
    std::vector<std::vector<vec>> L(P);
    std::vector<std::vector<mat>> dL(P);
    std::vector<vec> D(P), e(P);
    std::vector<mat> dD(P), de(P);
    double ll = 0.0;
    for (uword k = 0; k < P; ++k) {
      L[k].resize(k);
      dL[k].resize(k);
      for (uword j = 0; j < k; ++j) {
        const vec sigma = cov.G(k, j) * lam2M_ + cov.C(k, j) * lam_;
        vec acc(K, fill::zeros);
        for (uword m = 0; m < j; ++m) acc += L[k][m] % L[j][m] % D[m];
        L[k][j] = (sigma - acc) / D[j];
        if (scoring) {
          mat num = lam2M_ * vec(dG.tube(k, j)).t() +
                    lam_ * vec(dC.tube(k, j)).t();
          for (uword m = 0; m < j; ++m) {
            num -= dL[k][m].each_col() % (L[j][m] % D[m]) +
                   dL[j][m].each_col() % (L[k][m] % D[m]) +
                   dD[m].each_col() % (L[k][m] % L[j][m]);
          }
          num -= dD[j].each_col() % L[k][j];
          dL[k][j] = num.each_col() / D[j];
        }
      }
      vec explained(K, fill::zeros), predicted(K, fill::zeros);
      for (uword j = 0; j < k; ++j) {
        explained += square(L[k][j]) % D[j];
        predicted += L[k][j] % e[j];
      }
      const vec raw = v_.col(k) - explained;
      D[k] = clamp(raw, floor_, kInf);
      e[k] = u_.col(k) - predicted;
      ll -= accu(log(D[k]) + square(e[k]) / D[k]) / 2.0;
      if (scoring) {
        mat dDk(K, q, fill::zeros), dek(K, q, fill::zeros);
        for (uword j = 0; j < k; ++j) {
          dDk += dL[k][j].each_col() % (2.0 * L[k][j] % D[j]) +
                 dD[j].each_col() % square(L[k][j]);
          dek += dL[k][j].each_col() % e[j] + de[j].each_col() % L[k][j];
        }
        dD[k] = -(dDk.each_col() % indicator(raw > floor_));
        de[k] = -dek;
      }
    }
    Evaluation out{ll, vec(), mat()};
    if (!scoring) return out;
    const vec w(K, fill::ones);
    out.score = vec(q, fill::zeros);
    out.information = mat(q, q, fill::zeros);
    for (uword k = 0; k < P; ++k) {
      Evaluation term{0.0, vec(), mat()};
      varianceScoring(term, D[k], square(e[k]), w, dD[k]);
      meanScoring(term, D[k], e[k], de[k]);
      out.score += term.score;
      out.information += term.information;
    }
    return out;
  }

 private:
  mat u_, v_;
  vec lam_, lam2M_;
};

// A conditional covariance block as the decomposition reads it: itself where
// its Cholesky pivots clear the floor, else its eigenvalues floored; with its
// inverse and log determinant, and its eigenpairs where floored.
struct Block {
  mat matrix, inverse, vectors;
  vec values;
  double logDet;
  bool floored;
};

Block flooredBlock(const mat& S, double floor) {
  Block out;
  mat U;
  if (chol(U, S) && std::pow(min(U.diag()), 2) >= floor) {
    const mat Uinv = inv(trimatu(U));
    out.matrix = S;
    out.inverse = Uinv * Uinv.t();
    out.logDet = 2.0 * accu(log(U.diag()));
    out.floored = false;
    return out;
  }
  eig_sym(out.values, out.vectors, S);
  const vec f = clamp(out.values, floor, kInf);
  out.matrix = out.vectors * diagmat(f) * out.vectors.t();
  out.inverse = out.vectors * diagmat(1.0 / f) * out.vectors.t();
  out.logDet = accu(log(f));
  out.floored = true;
  return out;
}

// The floored block's derivative along dS: dS itself where nothing is
// floored, else Q (F o Q'dS Q) Q', F the divided differences of max(lambda,
// floor) over the eigenvalues (its slope where two meet).
mat flooredStep(const Block& block, const mat& dS, double floor) {
  if (!block.floored) return dS;
  const vec& l = block.values;
  const uword n = l.n_elem;
  mat F(n, n);
  for (uword i = 0; i < n; ++i) {
    for (uword j = 0; j < n; ++j) {
      const double gap = l[i] - l[j];
      if (std::abs(gap) > 1e-12 * std::max(1.0, std::abs(l[i]))) {
        F(i, j) = (std::max(l[i], floor) - std::max(l[j], floor)) / gap;
      } else {
        F(i, j) = l[i] > floor ? 1.0 : 0.0;
      }
    }
  }
  const mat& Q = block.vectors;
  return Q * (F % (Q.t() * dS * Q)) * Q.t();
}

// A group on several references: each condition's statistics on its
// reference's components, covarying with another's as G_kl B_ab + C_kl N_ab
// (a and b their references), and its own variances on its diagonal block.
class MixedLikelihood : public GroupLikelihood {
 public:
  explicit MixedLikelihood(const list& data) : GroupLikelihood(data) {
    const list joint(needed(data, "joint"));
    const list u(needed(joint, "u")), v(needed(joint, "v"));
    for (R_xlen_t k = 0; k < u.size(); ++k) {
      u_.push_back(numbers(u[k]));
      v_.push_back(numbers(v[k]));
    }
    const vec basisOf = neededVec(joint, "basisOf");
    basisOf_ = conv_to<uvec>::from(basisOf - 1.0);
    const list B(needed(joint, "B")), N(needed(joint, "N"));
    nBases_ = static_cast<uword>(std::lround(std::sqrt(B.size())));
    for (R_xlen_t i = 0; i < B.size(); ++i) {
      B_.push_back(numberMatrix(B[i]));
      N_.push_back(numberMatrix(N[i]));
    }
  }

  Evaluation at(const vec& theta, bool scoring,
                const uvec& free) const override {
    const uword P = P_;
    const Covariances cov = covariances(theta);
    std::vector<std::vector<mat>> L(P), LD(P);
    std::vector<Block> D(P);
    std::vector<vec> e(P);
    double ll = 0.0;
    for (uword k = 0; k < P; ++k) {
      L[k].resize(k);
      LD[k].resize(k);
      for (uword j = 0; j < k; ++j) {
        mat num = covariance(cov.G, cov.C, k, j);
        for (uword m = 0; m < j; ++m) {
          num -= LD[k][m] * L[j][m].t();
        }
        L[k][j] = num * D[j].inverse;
        LD[k][j] = L[k][j] * D[j].matrix;
      }
      mat S = diagmat(v_[k]);
      vec predicted(u_[k].n_elem, fill::zeros);
      for (uword j = 0; j < k; ++j) {
        S -= LD[k][j] * L[k][j].t();
        predicted += L[k][j] * e[j];
      }
      D[k] = flooredBlock((S + S.t()) / 2.0, floor_);
      e[k] = u_[k] - predicted;
      ll -= (D[k].logDet + dot(e[k], D[k].inverse * e[k])) / 2.0;
    }
    Evaluation out{ll, vec(), mat()};
    if (scoring) score(out, theta, free, L, LD, D, e);
    return out;
  }

 private:
  mat covariance(const mat& G, const mat& C, uword k, uword l) const {
    const uword at = basisOf_[k] + nBases_ * basisOf_[l];
    return G(k, l) * B_[at] + C(k, l) * N_[at];
  }

  // The score and information from each parameter's derivatives of the
  // loadings, conditional covariances and innovations, carried through the
  // decomposition one parameter at a time.
  void score(Evaluation& out, const vec& theta, const uvec& free,
             const std::vector<std::vector<mat>>& L,
             const std::vector<std::vector<mat>>& LD,
             const std::vector<Block>& D, const std::vector<vec>& e) const {
    const uword P = P_, q = free.n_elem;
    cube dG, dC;
    covarianceSteps(theta, free, dG, dC);
    std::vector<std::vector<mat>> dD(q, std::vector<mat>(P));
    std::vector<std::vector<vec>> de(q, std::vector<vec>(P));
    for (uword t = 0; t < q; ++t) {
      std::vector<std::vector<mat>> dL(P);
      for (uword k = 0; k < P; ++k) {
        dL[k].resize(k);
        for (uword j = 0; j < k; ++j) {
          mat num = covariance(dG.slice(t), dC.slice(t), k, j) -
                    L[k][j] * dD[t][j];
          for (uword m = 0; m < j; ++m) {
            num -= dL[k][m] * LD[j][m].t() +
                   L[k][m] * dD[t][m] * L[j][m].t() +
                   LD[k][m] * dL[j][m].t();
          }
          dL[k][j] = num * D[j].inverse;
        }
        mat dS(u_[k].n_elem, u_[k].n_elem, fill::zeros);
        vec dek(u_[k].n_elem, fill::zeros);
        for (uword j = 0; j < k; ++j) {
          const mat A = dL[k][j] * LD[k][j].t();
          dS -= A + A.t() + L[k][j] * dD[t][j] * L[k][j].t();
          dek -= dL[k][j] * e[j] + L[k][j] * de[t][j];
        }
        dD[t][k] = flooredStep(D[k], (dS + dS.t()) / 2.0, floor_);
        de[t][k] = dek;
      }
    }
    out.score = vec(q, fill::zeros);
    out.information = mat(q, q, fill::zeros);
    for (uword k = 0; k < P; ++k) {
      const vec a = D[k].inverse * e[k];
      std::vector<mat> W(q);
      std::vector<vec> Ide(q);
      for (uword t = 0; t < q; ++t) {
        W[t] = D[k].inverse * dD[t][k];
        Ide[t] = D[k].inverse * de[t][k];
      }
      for (uword t = 0; t < q; ++t) {
        out.score[t] += -dot(de[t][k], a) +
                        (dot(e[k], W[t] * a) - trace(W[t])) / 2.0;
        for (uword l = 0; l <= t; ++l) {
          const double entry =
              dot(de[t][k], Ide[l]) +
              accu(W[t] % W[l].t()) / 2.0;
          out.information(t, l) += entry;
          if (l != t) out.information(l, t) += entry;
        }
      }
    }
  }

  std::vector<vec> u_, v_;
  uvec basisOf_;
  uword nBases_;
  std::vector<mat> B_, N_;
};

std::unique_ptr<Likelihood> likelihoodOf(const std::string& model,
                                         const list& data) {
  if (model == "variance") {
    return std::unique_ptr<Likelihood>(new VarianceLikelihood(data));
  }
  if (model == "tied") {
    return std::unique_ptr<Likelihood>(new TiedLikelihood(data));
  }
  if (model == "nullSpace") {
    return std::unique_ptr<Likelihood>(new NullSpaceLikelihood(data));
  }
  if (model == "cross") {
    return std::unique_ptr<Likelihood>(new CrossLikelihood(data));
  }
  if (model == "reduced") {
    return std::unique_ptr<Likelihood>(new ReducedLikelihood(data));
  }
  if (model == "joint") {
    if (flag(list(needed(data, "joint")), "mixed")) {
      return std::unique_ptr<Likelihood>(new MixedLikelihood(data));
    }
    return std::unique_ptr<Likelihood>(new JointLikelihood(data));
  }
  cpp11::stop("Unknown likelihood: '%s'.", model.c_str());
}

// ---- The maximum ------------------------------------------------------------

vec toVec(const column& c) {
  vec out(c.size());
  for (long i = 0; i < c.size(); ++i) out[i] = c(i);
  return out;
}

column toColumn(const vec& v) {
  column out(v.n_elem);
  for (uword i = 0; i < v.n_elem; ++i) out(i) = v[i];
  return out;
}

// What a run reads of a likelihood at its free parameters x, the rest held
// at `full`: the value, and the scored evaluation, the last one kept since
// dlib asks for the gradient and then the Hessian at one point; and the best
// value it has seen, a line search that gives up being able to end below it.
class Problem {
 public:
  Problem(const Likelihood& model, const vec& full, const uvec& free)
      : model_(model), full_(full), free_(free) {}

  double value(const vec& x) {
    const double v = hasScored_ && same(x, scoredAt_)
                         ? scored_.value
                         : model_.at(thetaOf(x), false, free_).value;
    const double finite = std::isfinite(v) ? v : -kInf;
    if (finite > bestValue_) {
      bestValue_ = finite;
      best_ = x;
    }
    return finite;
  }

  const Evaluation& scored(const vec& x) {
    if (!hasScored_ || !same(x, scoredAt_)) {
      scored_ = model_.at(thetaOf(x), true, free_);
      scoredAt_ = x;
      hasScored_ = true;
    }
    return scored_;
  }

  double bestValue() const { return bestValue_; }
  const vec& best() const { return best_; }

 private:
  vec thetaOf(const vec& x) const {
    vec theta = full_;
    theta.elem(free_) = x;
    return theta;
  }

  const Likelihood& model_;
  vec full_;
  uvec free_;
  Evaluation scored_{0.0, vec(), mat()};
  vec scoredAt_, best_;
  bool hasScored_ = false;
  double bestValue_ = -kInf;
};

// The projected Newton step (above) at theta, where `at` is the scored
// evaluation, and its predicted rise in the log-likelihood: the information's
// quadratic model's on the free block, and the score's slope times the move
// for a parameter stepped onto its bound.
struct Direction {
  vec step;
  double rise;
};

Direction projectedNewton(const Evaluation& at, const vec& theta,
                          const vec& lower, const vec& upper) {
  const vec grad = -at.score;
  const vec curvature = at.information.diag();
  vec scaled = grad;
  for (uword i = 0; i < grad.n_elem; ++i) {
    if (curvature[i] > 0.0) scaled[i] /= curvature[i];
  }
  const vec projected = boxed(theta - scaled, lower, upper) - theta;
  const double size = projected.is_finite() ? max(abs(projected)) : kInf;
  Direction out{vec(theta.n_elem, fill::zeros), 0.0};
  std::vector<uword> freeIdx;
  for (uword i = 0; i < theta.n_elem; ++i) {
    const double eps = std::min(1e-3 * (upper[i] - lower[i]), size);
    if (theta[i] <= lower[i] + eps && grad[i] > 0) {
      out.step[i] = lower[i] - theta[i];
    } else if (theta[i] >= upper[i] - eps && grad[i] < 0) {
      out.step[i] = upper[i] - theta[i];
    } else {
      freeIdx.push_back(i);
    }
  }
  out.rise = -dot(grad, out.step);
  const uvec free = conv_to<uvec>::from(freeIdx);
  if (free.n_elem > 0) {
    mat H = at.information.submat(free, free);
    const double ridge = mean(abs(H.diag()));
    vec solved;
    if (std::isfinite(ridge) && ridge > 0.0) {
      H.diag() += 1e-12 * ridge;
      const vec g = grad.elem(free);
      if (solve(solved, H, g, solve_opts::no_approx)) {
        out.step.elem(free) = -solved;
        out.rise += dot(g, solved) / 2.0;
      }
    }
  }
  return out;
}

// The projected Newton step as dlib's search strategy reads it, cut back to
// the box where that still descends.
class ProjectedNewton {
 public:
  ProjectedNewton(Problem* problem, const vec& lower, const vec& upper)
      : problem_(problem), lower_(lower), upper_(upper) {}
  double get_wolfe_rho() const { return 0.01; }
  double get_wolfe_sigma() const { return 0.9; }
  unsigned long get_max_line_search_iterations() const { return 100; }

  template <typename T>
  const column get_next_direction(const T& x, const double, const T&) {
    const vec theta = toVec(x);
    const Evaluation& at = problem_->scored(theta);
    const vec d = projectedNewton(at, theta, lower_, upper_).step;
    const vec inside = boxed(theta + d, lower_, upper_) - theta;
    return toColumn(dot(-at.score, inside) < 0.0 ? inside : d);
  }

 private:
  Problem* problem_;
  vec lower_, upper_;
};

// dlib's stop strategy for a single iteration.
class OneIteration {
 public:
  template <typename T>
  bool should_continue_search(const T&, const double, const T&) {
    return calls_++ == 0;
  }

 private:
  int calls_ = 0;
};

struct Fit {
  vec par;
  double value;
  bool converged;
  mat information;
};

// The maximum from `start`, one of dlib's iterations at a time so that each
// line search starts from the whole step (dlib's own loop starts each from
// the last one's step length, which at a floor's kink can leave a run taking
// steps too small to read). It has converged where the next step's predicted
// rise is no more than `tolerance` of the log-likelihood's size or an
// iteration does not raise it, and not where `maxit` iterations come first
// or the likelihood is not finite.
Fit maximize(Problem& problem, const vec& start, const vec& lower,
             const vec& upper, unsigned long maxit, double tolerance) {
  column x = toColumn(boxed(start, lower, upper));
  bool converged = false, failed = false;
  try {
    for (unsigned long it = 0; it < maxit && !converged; ++it) {
      const double before = problem.value(toVec(x));
      dlib::find_min_box_constrained(
          ProjectedNewton(&problem, lower, upper), OneIteration(),
          [&problem](const column& c) { return -problem.value(toVec(c)); },
          [&problem](const column& c) {
            return toColumn(-problem.scored(toVec(c)).score);
          },
          x, toColumn(lower), toColumn(upper));
      const vec theta = toVec(x);
      const double after = problem.value(theta);
      const Direction next =
          projectedNewton(problem.scored(theta), theta, lower, upper);
      converged = !(after > before) ||
                  next.rise <= tolerance * std::max(1.0, std::abs(after));
    }
  } catch (const dlib::error&) { failed = true; }  // a non-finite likelihood
  vec par = toVec(x);
  double value = problem.value(par);
  if (problem.bestValue() > value) {
    par = problem.best();
    value = problem.bestValue();
  }
  return {par, value, converged && !failed, problem.scored(par).information};
}

uvec everyParameter(uword n) {
  return n == 0 ? uvec() : regspace<uvec>(0, n - 1);
}

}  // namespace

// The log-likelihood `model` gives `data` at theta and, where `scoring`
// asks, its score and expected information.
[[cpp11::register]]
cpp11::writable::list eigenLikelihoodCpp(const std::string& model,
                                         const list& data,
                                         const doubles& thetaR,
                                         bool scoring) {
  const std::unique_ptr<Likelihood> likelihood = likelihoodOf(model, data);
  const vec theta = numbers(thetaR);
  const Evaluation at =
      likelihood->at(theta, scoring, everyParameter(theta.n_elem));
  if (!scoring) {
    writable::list out({"value"_nm = as_sexp(at.value)});
    return out;
  }
  writable::list out({"value"_nm = as_sexp(at.value),
                      "score"_nm = as_doubles(at.score),
                      "information"_nm = as_doubles_matrix(at.information)});
  return out;
}

// The maximum of the log-likelihood `model` gives `data` from each of
// `starts` within [lower, upper], in the parameters `free` (every one where
// empty) with the rest held at `full`'s: one list(par, value, converged,
// information) per start.
[[cpp11::register]]
cpp11::writable::list eigenLikelihoodFitCpp(const std::string& model,
                                            const list& data,
                                            const list& starts,
                                            const doubles& lowerR,
                                            const doubles& upperR,
                                            const doubles& fullR,
                                            const integers& freeR,
                                            int maxit,
                                            bool valueOnly) {
  const std::unique_ptr<Likelihood> likelihood = likelihoodOf(model, data);
  const vec lower = numbers(lowerR), upper = numbers(upperR);
  uvec free = everyParameter(lower.n_elem);
  vec full(lower.n_elem, fill::zeros);
  if (freeR.size() > 0) {
    free = conv_to<uvec>::from(numbers(freeR) - 1.0);
    full = numbers(fullR);
  }
  const double tolerance = valueOnly ? 1e-10 : 1e-14;
  writable::list fits(starts.size());
  for (R_xlen_t i = 0; i < starts.size(); ++i) {
    Problem problem(*likelihood, full, free);
    const Fit fit = maximize(problem, numbers(starts[i]), lower, upper,
                             static_cast<unsigned long>(maxit), tolerance);
    fits[i] = writable::list(
        {"par"_nm = as_doubles(fit.par), "value"_nm = as_sexp(fit.value),
         "converged"_nm = as_sexp(fit.converged),
         "information"_nm = as_doubles_matrix(fit.information)});
  }
  return fits;
}
