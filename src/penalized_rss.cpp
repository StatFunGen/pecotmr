// penalized_rss.cpp -- coordinate descent for penalized regression on the RSS
// objective, along a path of penalty values.
//
// The package's one RSS coordinate-descent solver, for the penalties
//   LASSO (with glmnet's elastic-net mixing), MCP, SCAD, L0, L0L1, L0L2.
// It minimises
//   beta'R_s beta - 2 beta'z + penalty(beta),   R_s = (1 - s) R + s I,
// with z the marginal correlations and s = `shrink` (0 for an LD matrix
// already shrunk). On standardized data this is exactly twice ncvreg's
// objective, (1/2n)||y - X beta||^2 + P(beta) with (1/n)||x_j||^2 = 1, and
// twice L0Learn's, (1/2)||r||^2 + lambda0 ||beta||_0 + lambda1 ||beta||_1 +
// lambda2 ||beta||_2^2 with unit-norm columns and response, so a penalty
// value means here what it means there.
//
// LASSO, MCP and SCAD mix in a ridge term the way glmnet and ncvreg do:
// l1 = lambda * alpha, l2 = lambda * (1 - alpha), penalty = 2 P(beta; l1) +
// l2 beta'beta, and their path is over lambda. The L0 family's path is over
// lambda0 with lambda1 and lambda2 fixed: penalty = 2 lambda0 ||beta||_0 +
// 2 lambda1 ||beta||_1 + 2 lambda2 ||beta||_2^2.
//
// Two strategies walk the path. Both are warm-started down it, and both
// record the values at or above pathMax -- the smallest at which every
// coefficient is zero -- without solving:
//   strongRule      ncvreg's cdfit_gaussian: the sequential strong rule
//                   screens each step, coordinate descent cycles in index
//                   order over every coordinate ever active, and scans of the
//                   strong set and then of the rest add each coordinate that
//                   leaves zero. The L0 family has no strong rule, so there
//                   every coordinate is strong.
//   gradientRanked  L0Learn's (CD.h, CDL0.h, CDL012.h): each step cycles
//                   coordinates in decreasing |gradient| at the previous
//                   solution, screened to the support and the next
//                   screenSize, sweeps only the support once it has held for
//                   stableSweeps sweeps, and checks every zero coordinate at
//                   convergence. Unlike L0Learn, the support is always swept
//                   and the ranking is always the current gradient.
// A solve stops when no coefficient moved by `thr` or more in a sweep,
// ncvreg's test on this scale; L0Learn's relative change in the objective is
// not used, so `thr` means the same under both. The L0 family's default
// path, when the caller gives none, is L0Learn's adaptive lambda0 path under
// gradientRanked, solved at the midpoints of its values (off_edge()), and a
// geometric one from pathMax under strongRule, and L0Learn's swaps follow
// either.
//
// The coordinate update follows lassosum's elnet()
// (https://github.com/tshmak/lassosum, Copyright (c) 2016 tshmak, MIT). The
// MCP and SCAD operators and the strongRule path are ported from ncvreg
// (https://github.com/pbreheny/ncvreg, Copyright (c) Patrick Breheny, GPL-3),
// and the L0 thresholding, the gradientRanked path, the adaptive lambda0
// path and the swaps from L0Learn (https://github.com/hazimehh/L0Learn,
// Copyright (c) 2021 Hussein Hazimeh, MIT). Their license notices are in
// inst/COPYRIGHTS.
//
// CoupledPathSolver solves K conditions' paths together under a penalty on
// the difference between their effects, PROSPER's cross-ancestry coupling:
//   sum_k [beta_k'R_k,s beta_k - 2 beta_k'z_k + penalty(beta_k)]
//     + sum_j sum_{k < k' carrying j} c_kk' (beta_jk - beta_jk')^2,
// each condition against its own LD, with c_kk' the strength of the pair
// (zero for a pair not coupled). A coordinate then sees the partial
// correlation z_j + sum_k' c_kk' beta_jk' over its partners and the curvature
// (R_s)_jj + sum_k' c_kk' over them, so every thresholding operator
// stays exact; the conditions step through each point together, each
// refreshing its coupling from the others' current coefficients before every
// sweep and scan.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "rss_ld.h"

using namespace cpp11;
using namespace arma;

// ---- Penalties and strategies --------------------------------------------

enum Penalty { PEN_LASSO = 0, PEN_MCP = 1, PEN_SCAD = 2, PEN_L0 = 3 };

enum Strategy { STRONG_RULE = 0, GRADIENT_RANKED = 1 };

static Penalty parse_penalty(const std::string& s) {
  if (s == "lasso") return PEN_LASSO;
  if (s == "MCP")   return PEN_MCP;
  if (s == "SCAD")  return PEN_SCAD;
  if (s == "L0" || s == "L0L1" || s == "L0L2") return PEN_L0;
  cpp11::stop("Unknown penalty: '%s'. Use lasso, MCP, SCAD, L0, L0L1, or L0L2.", s.c_str());
  return PEN_LASSO; // unreachable
}

static Strategy parse_strategy(const std::string& s) {
  if (s == "strongRule")     return STRONG_RULE;
  if (s == "gradientRanked") return GRADIENT_RANKED;
  cpp11::stop("Unknown strategy: '%s'. Use strongRule or gradientRanked.", s.c_str());
  return STRONG_RULE; // unreachable
}

// The penalty weights at one point of the path: l1 and l2 multiply the L1
// and L2 terms, lambda0 the L0 term.
struct PointPenalty {
  double l1;
  double l2;
  double lambda0;
};

// ---- Thresholding operators ----------------------------------------------
// v = (R_s)_jj, the curvature of coordinate j's quadratic.
//
// Each operator is the exact minimizer of the one-coordinate problem
// v b^2 - 2 z b + 2 P(b; l1) + l2 b^2, whose curvature is a = v + l2. ncvreg
// writes the curvature v (1 + l2), which is the same only when v = 1: a
// correlation matrix's diagonal, but not a diagonal the shrinkage or a
// genotype factor leaves slightly off 1.

static inline double lasso_thresh(double z, double l1, double l2, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  return std::copysign(az - l1, z) / (v + l2);
}

// MCP (minimax concave penalty) — Breheny & Huang 2011
// gamma > 1 controls concavity; gamma -> inf recovers LASSO
static inline double mcp_thresh(double z, double l1, double l2,
                                double gamma, double v) {
  double az = std::abs(z);
  double a = v + l2;
  if (az <= l1) return 0.0;
  if (az <= gamma * l1 * a)
    return std::copysign(az - l1, z) / (a - 1.0 / gamma);
  return z / a;
}

// SCAD (smoothly clipped absolute deviation) — Fan & Li 2001
// gamma > 2 controls transition; default 3.7
static inline double scad_thresh(double z, double l1, double l2,
                                 double gamma, double v) {
  double az = std::abs(z);
  double a = v + l2;
  if (az <= l1) return 0.0;
  if (az <= l1 * (1.0 + a))
    return std::copysign(az - l1, z) / a;
  if (az <= gamma * l1 * a)
    return std::copysign(az - gamma * l1 / (gamma - 1.0), z) /
           (a - 1.0 / (gamma - 1.0));
  return z / a;
}

// L0 (+L1+L2) — Hazimeh & Mazumder 2020
// Soft-threshold for L1, then hard-threshold for L0. A coefficient exactly at
// its threshold stays zero, as L0Learn's does and as the lasso's does at l1.
// lambda0 = L0 penalty weight, l1 = L1 weight, l2 = L2 weight
static inline double l0_thresh(double z, double lambda0,
                               double l1, double l2, double v) {
  double az = std::abs(z);
  if (az <= l1) return 0.0;
  double denom = v + 2.0 * l2;
  double reg = (az - l1) / denom;
  double thr = std::sqrt(2.0 * lambda0 / denom);
  if (reg <= thr) return 0.0;
  return std::copysign(reg, z);
}

static inline double apply_threshold(Penalty pen, double t, double v,
                                     const PointPenalty& pp, double gamma) {
  switch (pen) {
    case PEN_LASSO: return lasso_thresh(t, pp.l1, pp.l2, v);
    case PEN_MCP:   return mcp_thresh(t, pp.l1, pp.l2, gamma, v);
    case PEN_SCAD:  return scad_thresh(t, pp.l1, pp.l2, gamma, v);
    case PEN_L0:    break;
  }
  return l0_thresh(t, pp.lambda0, pp.l1, pp.l2, v);
}

// ---- Full penalty value (for fbeta computation) --------------------------

static double mcp_penalty_single(double bj, double lambda, double gamma) {
  double ab = std::abs(bj);
  if (ab <= lambda * gamma)
    return lambda * ab - ab * ab / (2.0 * gamma);
  return 0.5 * lambda * lambda * gamma;
}

static double scad_penalty_single(double bj, double lambda, double gamma) {
  double ab = std::abs(bj);
  if (ab <= lambda)
    return lambda * ab;
  if (ab <= gamma * lambda)
    return (2.0 * gamma * lambda * ab - ab * ab - lambda * lambda) /
           (2.0 * (gamma - 1.0));
  return lambda * lambda * (gamma + 1.0) / 2.0;
}

// The penalty the thresholding operators minimize, read from the same
// weights the descent used: 2 P(beta; l1) + l2 beta'beta for LASSO, MCP and
// SCAD, and 2 lambda0 ||beta||_0 + 2 l1 ||beta||_1 + 2 l2 ||beta||_2^2 for
// the L0 family.
static double compute_penalty(Penalty pen, const vec& beta,
                              const PointPenalty& pp, double gamma) {
  double val = 0.0;
  uword p = beta.n_elem;
  switch (pen) {
    case PEN_LASSO:
      val = 2.0 * pp.l1 * sum(abs(beta)) + pp.l2 * dot(beta, beta);
      break;
    case PEN_MCP:
      for (uword j = 0; j < p; j++)
        val += 2.0 * mcp_penalty_single(beta(j), pp.l1, gamma);
      val += pp.l2 * dot(beta, beta);
      break;
    case PEN_SCAD:
      for (uword j = 0; j < p; j++)
        val += 2.0 * scad_penalty_single(beta(j), pp.l1, gamma);
      val += pp.l2 * dot(beta, beta);
      break;
    case PEN_L0: {
      double nnz = static_cast<double>(accu(beta != 0.0));
      val = 2.0 * pp.lambda0 * nnz + 2.0 * pp.l1 * sum(abs(beta)) +
            2.0 * pp.l2 * dot(beta, beta);
      break;
    }
  }
  return val;
}

// ---- Settings ------------------------------------------------------------

// The constants the path works to (R's .rssPathControl()).
struct PathControl {
  int points;        // length of a default geometric path
  double ratio;      // its last value over its first
  int screenSize;    // gradientRanked: coordinates swept outside the support
  int stableSweeps;  // gradientRanked: sweeps the support holds before only
                     // it is swept
  int maxSupport;    // adaptive lambda0 path: stop at this many non-zeros
  int maxPoints;     // adaptive lambda0 path: at most this many points
  bool offEdge;      // adaptive lambda0 path: solve it at the midpoints of
                     // its values (off_edge())
};

static PathControl parse_control(const list& control) {
  PathControl c;
  c.points       = as_cpp<int>(control["points"]);
  c.ratio        = as_cpp<double>(control["ratio"]);
  c.screenSize   = as_cpp<int>(control["screenSize"]);
  c.stableSweeps = as_cpp<int>(control["stableSweeps"]);
  c.maxSupport   = as_cpp<int>(control["maxSupport"]);
  c.maxPoints    = as_cpp<int>(control["maxPoints"]);
  c.offEdge      = as_cpp<bool>(control["offEdge"]);
  return c;
}

struct PathSettings {
  Penalty pen;
  Strategy strategy;
  double gamma;
  double alpha;    // LASSO, MCP, SCAD
  double lambda1;  // L0 family
  double lambda2;  // L0 family
  double thr;
  int maxiter;
  int maxSwaps;
  PathControl control;
};

// The path's points, in the order solved.
struct PathResult {
  std::vector<vec> beta;
  std::vector<double> value;
  std::vector<int> conv;
  std::vector<double> loss;
  std::vector<double> fbeta;
};

// Whether the swaps keep the columns of R_s they read: a factor's column
// costs O(mp) to form, a dense one only the O(p) of copying it.
template <class Ld> struct KeepsColumns { static const bool value = false; };
template <> struct KeepsColumns<FactorLd> { static const bool value = true; };

// ---- The path solver -----------------------------------------------------
// `ld` answers R_s's entries for the current beta (rss_ld.h): DenseLd tracks
// R * beta as a running sum, FactorLd tracks F * beta. It is built once for
// the whole path, its running product carried from point to point with beta.
// `shift_` and `extra_` carry the coupling to other studies' coefficients
// (CoupledPathSolver): zero for a study fitted alone, where adding them
// changes no value.

template <class Ld> class CoupledPathSolver;

template <class Ld>
class PathSolver {
  friend class CoupledPathSolver<Ld>;

public:
  PathSolver(const PathSettings& settings, Ld& ld, const vec& z)
      : st_(settings), ld_(ld), z_(z), p_(z.n_elem), v_(z.n_elem),
        beta_(z.n_elem, fill::zeros), absGrad_(abs(z)),
        shift_(z.n_elem, fill::zeros), extra_(z.n_elem, fill::zeros),
        ever_(z.n_elem, 0), strong_(z.n_elem, 0), supportChanged_(false),
        pathMax_(0.0) {
    for (uword j = 0; j < p_; j++) v_(j) = ld_.diag(j);
  }

  // The caller's path, or for the L0 family with none the default one.
  PathResult run(const vec& path) {
    pathMax_ = path_max();
    if (path.n_elem > 0) {
      double prev = pathMax_;
      for (uword k = 0; k < path.n_elem; k++) {
        solve_point(path(k), prev);
        prev = path(k);
      }
    } else if (st_.strategy == GRADIENT_RANKED) {
      adaptive_path();
    } else {
      geometric_path();
    }
    return result_;
  }

private:
  enum ScanSet { IN_STRONG, OUTSIDE_STRONG };

  const PathSettings& st_;
  Ld& ld_;
  const vec& z_;
  uword p_;
  vec v_;
  vec beta_;
  vec absGrad_;  // |z - R_s beta| at the last recorded point (or swap)
  vec shift_;    // coupling: the partners' coefficients, each times its c
  vec extra_;    // coupling: the partners' strengths c, summed
  std::vector<unsigned char> ever_;    // strongRule: ever active
  std::vector<unsigned char> strong_;  // strongRule: ever strong
  std::vector<uword> everList_;        // ever_ in index order
  std::map<uword, vec> columns_;       // the swaps' columns (FactorLd)
  bool supportChanged_;
  double pathMax_;
  PathResult result_;

  PointPenalty point_penalty(double value) const {
    if (st_.pen == PEN_L0) return PointPenalty{st_.lambda1, st_.lambda2, value};
    return PointPenalty{value * st_.alpha, value * (1.0 - st_.alpha), 0.0};
  }

  // The largest path value at which a zero coefficient with gradient |g|
  // leaves zero: |g| / alpha for LASSO, MCP and SCAD (infinite for a ridge),
  // and (|g| - l1)^2 / (2 (v + 2 l2)) for the L0 family, L0Learn's
  // (Xrmax - lambda1)^2 / (2 (1 + 2 lambda2)) at v = 1.
  double entry_value(double absg, double v) const {
    if (st_.pen == PEN_L0) {
      double excess = std::max(absg - st_.lambda1, 0.0);
      return excess * excess / (2.0 * (v + 2.0 * st_.lambda2));
    }
    if (st_.alpha <= 0.0) return std::numeric_limits<double>::infinity();
    return absg / st_.alpha;
  }

  double path_max() const {
    double best = 0.0;
    for (uword j = 0; j < p_; j++)
      best = std::max(best, entry_value(std::abs(z_(j)), v_(j) + extra_(j)));
    return best;
  }

  // The largest entry value among the zero coefficients, from the gradient
  // at the last recorded point.
  double largest_entry() const {
    double best = 0.0;
    for (uword j = 0; j < p_; j++)
      if (beta_(j) == 0.0)
        best = std::max(best, entry_value(absGrad_(j), v_(j) + extra_(j)));
    return best;
  }

  uword support_size() const { return accu(beta_ != 0.0); }

  // One point: recorded as zero at or above pathMax while beta is still zero
  // (ncvreg's closed-form first point), solved otherwise; `prev` is the
  // point before it, or pathMax for the first.
  void solve_point(double value, double prev) {
    PointPenalty pp = point_penalty(value);
    if (support_size() == 0 && value >= pathMax_ * (1.0 - 1e-12)) {
      finish_point(value, pp, 1);
      return;
    }
    int conv = solve(pp, value, prev);
    if (st_.pen == PEN_L0 && st_.maxSwaps > 0)
      conv = std::min(conv, swap_phase(pp, value, prev));
    finish_point(value, pp, conv);
  }

  int solve(const PointPenalty& pp, double value, double prev) {
    if (st_.strategy == STRONG_RULE) return solve_strong(pp, value, prev);
    return solve_ranked(pp);
  }

  // Record the point. The gradient z - R_s beta computed here feeds the next
  // point's screening, ranking and adaptive lambda0, and its loss
  // beta'R_s beta - 2 z'beta.
  void finish_point(double value, const PointPenalty& pp, int conv) {
    vec Rsb = ld_.product(beta_);
    absGrad_ = abs(gradient(Rsb));
    double loss = dot(beta_, Rsb) - 2.0 * dot(z_, beta_);
    result_.beta.push_back(beta_);
    result_.value.push_back(value);
    result_.conv.push_back(conv);
    result_.loss.push_back(loss);
    result_.fbeta.push_back(loss + compute_penalty(st_.pen, beta_, pp, st_.gamma));
  }

  // The gradient's negative half, z - R_s beta less the coupling's pull,
  // from the product R_s beta.
  vec gradient(const vec& Rsb) const {
    return z_ + shift_ - Rsb - extra_ % beta_;
  }

  // One coordinate's exact update; returns how far it moved.
  double update_coordinate(uword j, const PointPenalty& pp) {
    double bj = beta_(j);
    // Partial correlation: z_j - sum_{k != j} (R_s)_jk beta_k, plus the
    // coupling to the partners' coefficients
    double t = z_(j) + shift_(j) - ld_.row(j, beta_) + v_(j) * bj;
    double b = apply_threshold(st_.pen, t, v_(j) + extra_(j), pp, st_.gamma);
    if (b == bj) return 0.0;
    beta_(j) = b;
    ld_.update(j, b - bj);
    if ((bj == 0.0) != (b == 0.0)) supportChanged_ = true;
    return std::abs(b - bj);
  }

  // One sweep over `set`, in its order; returns the largest move.
  double sweep(const std::vector<uword>& set, const PointPenalty& pp) {
    double dlx = 0.0;
    for (uword j : set) dlx = std::max(dlx, update_coordinate(j, pp));
    cpp11::check_user_interrupt();
    return dlx;
  }

  // ---- strongRule (ncvreg's cdfit_gaussian) ----

  // ncvreg's sequential strong rule at `value` after `prev`: a coordinate
  // whose |gradient| at the previous solution is above the bound joins the
  // strong set. Minus infinity (every coordinate strong) where the rule does
  // not apply or is not finite.
  double strong_cutoff(double value, double prev) const {
    double cut = -std::numeric_limits<double>::infinity();
    switch (st_.pen) {
      case PEN_LASSO: cut = 2.0 * value - prev; break;
      case PEN_MCP:
        cut = value + st_.gamma / (st_.gamma - 1.0) * (value - prev);
        break;
      case PEN_SCAD:
        cut = value + st_.gamma / (st_.gamma - 2.0) * (value - prev);
        break;
      case PEN_L0: break;
    }
    cut *= st_.alpha;
    if (!std::isfinite(cut)) return -std::numeric_limits<double>::infinity();
    return cut;
  }

  // Descent over the ever-active set until it converges, a scan of the
  // strong set for coordinates that leave zero, and once that finds none a
  // scan of the rest; any found joins the ever-active set and the descent
  // resumes. 1 when a scan of the rest finds none within `maxiter` sweeps.
  int solve_strong(const PointPenalty& pp, double value, double prev) {
    double cut = strong_cutoff(value, prev);
    for (uword j = 0; j < p_; j++)
      if (!strong_[j] && absGrad_(j) > cut) strong_[j] = 1;
    int budget = st_.maxiter;
    while (true) {
      do {
        if (!descend(everList_, pp, budget)) return 0;
      } while (scan(pp, IN_STRONG));
      if (!scan(pp, OUTSIDE_STRONG)) return 1;
    }
  }

  // Coordinate descent over `set` until no coefficient moves by thr or more
  // in a sweep; false if the `budget` of sweeps runs out first.
  bool descend(const std::vector<uword>& set, const PointPenalty& pp,
               int& budget) {
    while (budget > 0) {
      budget--;
      if (sweep(set, pp) < st_.thr) return true;
    }
    return false;
  }

  // ncvreg's scans: each coordinate of the set not yet active (so at zero)
  // updated in index order, joining the active and strong sets if it leaves
  // zero. Whether any did.
  bool scan(const PointPenalty& pp, ScanSet which) {
    bool added = false;
    for (uword j = 0; j < p_; j++) {
      if (ever_[j]) continue;
      if ((which == IN_STRONG) != (strong_[j] != 0)) continue;
      if (update_coordinate(j, pp) > 0.0) {
        ever_[j] = 1;
        strong_[j] = 1;
        added = true;
      }
    }
    if (added) rebuild_ever_list();
    return added;
  }

  void rebuild_ever_list() {
    everList_.clear();
    for (uword j = 0; j < p_; j++)
      if (ever_[j]) everList_.push_back(j);
  }

  // ---- gradientRanked (L0Learn's CD) ----

  // The support and the screenSize coordinates outside it with the largest
  // |gradient| at the previous solution, in decreasing |gradient| (ties by
  // index).
  std::vector<uword> ranked_set() const {
    std::vector<uword> order(p_);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [this](uword a, uword b) {
      return absGrad_(a) > absGrad_(b);
    });
    std::vector<uword> set;
    int screened = 0;
    for (uword j : order) {
      if (beta_(j) != 0.0) {
        set.push_back(j);
      } else if (screened < st_.control.screenSize) {
        set.push_back(j);
        screened++;
      }
    }
    return set;
  }

  // Sweeps of the ranked set; once its support has held for stableSweeps
  // sweeps only the support is swept (L0Learn's RestrictSupport). At
  // convergence every zero coordinate is checked (L0Learn's CWMinCheck);
  // those that leave zero join the set and the sweeps resume. 1 when a check
  // finds none within `maxiter` sweeps.
  int solve_ranked(const PointPenalty& pp) {
    std::vector<uword> set = ranked_set();
    std::vector<unsigned char> inSet(p_, 0);
    for (uword j : set) inSet[j] = 1;
    int stable = 0;
    for (int budget = st_.maxiter; budget > 0; budget--) {
      supportChanged_ = false;
      double dlx = sweep(set, pp);
      stable = supportChanged_ ? 0 : stable + 1;
      if (stable == st_.control.stableSweeps) restrict_to_support(set, inSet);
      if (dlx >= st_.thr) continue;
      std::vector<uword> added = check_zeros(pp);
      if (added.empty()) return 1;
      for (uword j : added) {
        if (inSet[j]) continue;
        set.push_back(j);
        inSet[j] = 1;
      }
      stable = 0;
    }
    return 0;
  }

  void restrict_to_support(std::vector<uword>& set,
                           std::vector<unsigned char>& inSet) const {
    std::vector<uword> kept;
    for (uword j : set) {
      if (beta_(j) != 0.0) {
        kept.push_back(j);
      } else {
        inSet[j] = 0;
      }
    }
    set.swap(kept);
  }

  // Each zero coordinate updated in index order; those that left zero.
  std::vector<uword> check_zeros(const PointPenalty& pp) {
    std::vector<uword> added;
    for (uword j = 0; j < p_; j++) {
      if (beta_(j) != 0.0) continue;
      if (update_coordinate(j, pp) > 0.0) added.push_back(j);
    }
    return added;
  }

  // ---- Paths for the L0 family without one ----

  // L0Learn's lambda0 path (Grid1D.cpp): from pathMax, each next value is
  // 0.99 of the largest entry value outside the support, 0.97 of the last
  // value when that is not smaller, and after a point that left the support
  // unchanged the smaller of 0.8 of the last value and 0.97 of the entry
  // value. It stops at maxSupport non-zeros, maxPoints points, or when no
  // coordinate can enter.
  void adaptive_path() {
    double value = pathMax_;
    solve_point(value, value);
    for (int k = 1; k < st_.control.maxPoints; k++) {
      if (support_size() >= static_cast<uword>(st_.control.maxSupport)) break;
      double entry = largest_entry();
      if (!(entry > 0.0)) break;
      double next;
      if (same_support_as_previous()) {
        next = std::min(0.8 * value, 0.97 * entry);
      } else {
        next = 0.99 * entry >= value ? 0.97 * value : 0.99 * entry;
      }
      solve_point(next, value);
      value = next;
    }
  }

  // Whether the last two recorded points share a support.
  bool same_support_as_previous() const {
    std::size_t k = result_.beta.size();
    if (k < 2) return false;
    const vec& last = result_.beta[k - 1];
    const vec& before = result_.beta[k - 2];
    for (uword j = 0; j < p_; j++)
      if ((last(j) != 0.0) != (before(j) != 0.0)) return false;
    return true;
  }

  // `points` values from pathMax down to `ratio` of it, evenly on a log
  // scale.
  void geometric_path() {
    int points = st_.control.points;
    double prev = pathMax_;
    for (int k = 0; k < points; k++) {
      double frac = points > 1 ? static_cast<double>(k) / (points - 1) : 0.0;
      double value = pathMax_ * std::pow(st_.control.ratio, frac);
      solve_point(value, prev);
      prev = value;
    }
  }

  // ---- L0Learn's swaps (CDL012Swaps.cpp) ----

  // Rounds of swaps, each followed by a fresh solve under the strategy, until
  // a round finds none or maxSwaps rounds have run. The smallest convergence
  // flag of those solves.
  int swap_phase(const PointPenalty& pp, double value, double prev) {
    int conv = 1;
    for (int round = 0; round < st_.maxSwaps; round++) {
      if (!swap_once(pp)) break;
      conv = std::min(conv, solve(pp, value, prev));
    }
    return conv;
  }

  // One round: the first support member, in index order, whose best
  // replacement outside the support has a larger partial correlation than
  // it can be held at, (v + 2 l2) |beta_i| + l1, is zeroed and the
  // replacement set to its own L0 update. False when no member qualifies.
  bool swap_once(const PointPenalty& pp) {
    vec g = gradient(ld_.product(beta_));
    for (uword i = 0; i < p_; i++) {
      double bi = beta_(i);
      if (bi == 0.0) continue;
      vec col = column(i);
      double best = -1.0;
      uword jstar = p_;
      for (uword j = 0; j < p_; j++) {
        if (beta_(j) != 0.0) continue;
        double r = std::abs(g(j) + col(j) * bi);
        if (r > best) {
          best = r;
          jstar = j;
        }
      }
      if (jstar == p_) continue;
      double vstar = v_(jstar) + extra_(jstar);
      if (best <= (vstar + 2.0 * pp.l2) * std::abs(bi) + pp.l1) continue;
      beta_(i) = 0.0;
      ld_.update(i, -bi);
      double t = z_(jstar) + shift_(jstar) - ld_.row(jstar, beta_);
      double b = l0_thresh(t, pp.lambda0, pp.l1, pp.l2, vstar);
      if (b != 0.0) {
        beta_(jstar) = b;
        ld_.update(jstar, b);
        ever_[jstar] = 1;
        strong_[jstar] = 1;
        rebuild_ever_list();
      }
      absGrad_ = abs(gradient(ld_.product(beta_)));
      return true;
    }
    return false;
  }

  vec column(uword i) {
    if (!KeepsColumns<Ld>::value) return ld_.column(i);
    auto it = columns_.find(i);
    if (it == columns_.end()) it = columns_.emplace(i, ld_.column(i)).first;
    return it->second;
  }
};

// ---- Entry-point helpers -------------------------------------------------

static PathSettings parse_settings(const std::string& penaltyStr,
                                   const std::string& strategyStr,
                                   double gamma, double alpha, double lambda1,
                                   double lambda2, double thr, int maxiter,
                                   int maxSwaps, const list& control) {
  PathSettings st;
  st.pen = parse_penalty(penaltyStr);
  st.strategy = parse_strategy(strategyStr);
  st.gamma = gamma;
  st.alpha = alpha;
  st.lambda1 = lambda1;
  st.lambda2 = lambda2;
  st.thr = thr;
  st.maxiter = maxiter;
  st.maxSwaps = maxSwaps;
  st.control = parse_control(control);
  return st;
}

// Whether the path to solve is the L0 family's adaptive default, which is
// solved off its edges (off_edge()).
static bool adaptive_default(const PathSettings& st, const vec& path) {
  return path.n_elem == 0 && st.pen == PEN_L0 &&
         st.strategy == GRADIENT_RANKED && st.control.offEdge;
}

// The adaptive lambda0 path's values sit just below the full data's entry
// thresholds, each where one more coordinate leaves zero. A held-out refit at
// such a value often does not take the support the full data did, so
// selection would score a sparser model than it returns, and a coordinate
// that hurts would pass. The path is solved instead at the geometric
// midpoints of consecutive values, the first kept, which on real genotypes
// recovered all the accuracy the edges lost.
static vec off_edge(const std::vector<double>& values) {
  vec out(values.size());
  for (uword k = 0; k < values.size(); k++)
    out(k) = k == 0 ? values[0] : std::sqrt(values[k - 1] * values[k]);
  return out;
}

// A caller's path must be non-increasing (R's .rssSolvePath() sorts it), and
// only the L0 family has a default to fall back on when it is empty.
static void check_path(const vec& path, Penalty pen) {
  if (path.n_elem == 0 && pen != PEN_L0)
    cpp11::stop("A LASSO, MCP or SCAD path needs at least one lambda.");
  for (uword k = 1; k < path.n_elem; k++)
    if (path(k) > path(k - 1))
      cpp11::stop("The penalty path must be non-increasing.");
}

// The path's results, in the shape both entry points return.
static cpp11::writable::list path_result(const PathResult& res, uword p) {
  using namespace cpp11::literals;
  uword npoints = res.value.size();
  mat beta(p, npoints);
  for (uword k = 0; k < npoints; k++) beta.col(k) = res.beta[k];
  vec value(res.value);
  Col<int> conv(npoints);
  for (uword k = 0; k < npoints; k++) conv(k) = res.conv[k];
  vec loss(res.loss);
  vec fbeta(res.fbeta);
  writable::list result({
    "beta"_nm   = as_doubles_matrix(beta),
    "lambda"_nm = as_doubles(value),
    "conv"_nm   = as_integers(conv),
    "loss"_nm   = as_doubles(loss),
    "fbeta"_nm  = as_doubles(fbeta)
  });
  return result;
}

// ---- The coupled path solver ---------------------------------------------
// K studies' paths, solved together at each point of one shared path. A
// variant is the same variant in every study carrying it, matched by its
// index into the studies' union (`unionIdx`), in one allele frame. Each study
// keeps its own PathSolver -- its LD, screening sets and running product --
// and the studies step through every stage of a point in turn: a sweep, a
// scan, a check of the zeros, a swap. Before each, a study's coupling is
// refreshed from the others' current coefficients (Gauss-Seidel across the
// studies), so the coupling is felt at every pass rather than once per solve.

template <class Ld>
class CoupledPathSolver {
public:
  CoupledPathSolver(const PathSettings& settings,
                    std::vector<std::unique_ptr<Ld>>& lds,
                    const std::vector<vec>& zs,
                    const std::vector<uvec>& unionIdx, uword nUnion,
                    const mat& strength, bool grouped)
      : st_(settings), unionIdx_(unionIdx), members_(nUnion),
        strength_(strength), grouped_(grouped),
        everGroup_(nUnion, 0), strongGroup_(nUnion, 0) {
    for (uword k = 0; k < zs.size(); k++) {
      solvers_.emplace_back(new PathSolver<Ld>(settings, *lds[k], zs[k]));
      for (uword i = 0; i < unionIdx[k].n_elem; i++)
        members_[unionIdx[k](i)].push_back(std::make_pair(k, i));
    }
    // Under the group penalty a variant's coupling is solved inside its
    // group's update (update_group()), so the conditions carry none apart.
    if (grouped_) return;
    for (uword k = 0; k < solvers_.size(); k++) {
      vec& extra = solvers_[k]->extra_;
      for (uword i = 0; i < extra.n_elem; i++) {
        double sum = 0.0;
        for (const auto& m : members_[unionIdx_[k](i)])
          if (m.first != k) sum += strength_(k, m.first);
        extra(i) = sum;
      }
    }
  }

  // The caller's path, or with none a geometric one from the joint pathMax.
  void run(const vec& path) {
    pathMax_ = 0.0;
    if (grouped_) {
      for (uword j = 0; j < members_.size(); j++)
        pathMax_ = std::max(pathMax_, group_entry(j, true));
    } else {
      for (auto& s : solvers_) pathMax_ = std::max(pathMax_, s->path_max());
    }
    for (auto& s : solvers_) s->pathMax_ = pathMax_;
    vec points = path.n_elem > 0 ? path : geometric_path();
    double prev = pathMax_;
    for (uword k = 0; k < points.n_elem; k++) {
      solve_point(points(k), prev);
      prev = points(k);
    }
  }

  cpp11::writable::list result() const {
    using namespace cpp11::literals;
    writable::list studies(solvers_.size());
    for (uword k = 0; k < solvers_.size(); k++)
      studies[static_cast<R_xlen_t>(k)] =
          path_result(solvers_[k]->result_, solvers_[k]->p_);
    writable::list out({
      "studies"_nm  = studies,
      "lambda"_nm   = as_doubles(vec(solvers_[0]->result_.value)),
      "coupling"_nm = as_doubles(vec(coupling_))
    });
    return out;
  }

private:
  typedef typename PathSolver<Ld>::ScanSet ScanSet;

  const PathSettings& st_;
  const std::vector<uvec>& unionIdx_;
  std::vector<std::vector<std::pair<uword, uword>>> members_;
  const mat& strength_;  // c_kk', one per pair of conditions
  bool grouped_;         // each variant's coefficients one group
  std::vector<unsigned char> everGroup_;    // grouped: ever active
  std::vector<unsigned char> strongGroup_;  // grouped: ever strong
  double pathMax_ = 0.0;
  std::vector<std::unique_ptr<PathSolver<Ld>>> solvers_;
  std::vector<double> coupling_;  // the coupling term at each point

  // `points` values from the joint pathMax down to `ratio` of it, evenly on
  // a log scale, as PathSolver::geometric_path() takes them.
  vec geometric_path() const {
    int points = st_.control.points;
    vec out(points);
    for (int k = 0; k < points; k++) {
      double frac = points > 1 ? static_cast<double>(k) / (points - 1) : 0.0;
      out(k) = pathMax_ * std::pow(st_.control.ratio, frac);
    }
    return out;
  }

  // Condition k's coupling: the sum of each of its coefficients' partners in
  // the other conditions, each times the pair's strength. None under the
  // group penalty, whose updates solve the coupling themselves.
  void refresh_shift(uword k) {
    if (grouped_) return;
    vec& shift = solvers_[k]->shift_;
    for (uword i = 0; i < shift.n_elem; i++) {
      double sum = 0.0;
      for (const auto& m : members_[unionIdx_[k](i)])
        if (m.first != k)
          sum += strength_(k, m.first) * solvers_[m.first]->beta_(m.second);
      shift(i) = sum;
    }
  }

  // The sum over variants of c_kk' (beta_jk - beta_jk')^2 over every pair of
  // conditions carrying the variant.
  double coupling_term() const {
    double total = 0.0;
    for (const auto& carriers : members_) {
      for (uword a = 0; a < carriers.size(); a++) {
        for (uword b = a + 1; b < carriers.size(); b++) {
          double d = solvers_[carriers[a].first]->beta_(carriers[a].second) -
                     solvers_[carriers[b].first]->beta_(carriers[b].second);
          total += strength_(carriers[a].first, carriers[b].first) * d * d;
        }
      }
    }
    return total;
  }

  bool all_zero() const {
    for (const auto& s : solvers_)
      if (s->support_size() > 0) return false;
    return true;
  }

  // One point, as PathSolver::solve_point() takes it: recorded as zero at or
  // above the joint pathMax while every study is still zero.
  void solve_point(double value, double prev) {
    PointPenalty pp = solvers_[0]->point_penalty(value);
    if (all_zero() && value >= pathMax_ * (1.0 - 1e-12)) {
      finish(value, pp, 1);
      return;
    }
    int conv = solve(pp, value, prev);
    if (st_.pen == PEN_L0 && st_.maxSwaps > 0)
      conv = std::min(conv, grouped_ ? group_swap_phase(pp, value, prev)
                                     : swap_phase(pp, value, prev));
    finish(value, pp, conv);
  }

  int solve(const PointPenalty& pp, double value, double prev) {
    if (grouped_) return solve_grouped(pp, value, prev);
    if (st_.strategy == STRONG_RULE) return solve_strong(pp, value, prev);
    return solve_ranked(pp);
  }

  // ---- The group penalty across the conditions (grouped_) ----
  // Each variant's coefficients in the conditions carrying it are one group
  // under 2 P(||beta_j||; l1 sqrt(K_j)), lambda0 once per group for the L0
  // family, with the coupling between them, so a variant is in every
  // condition's weights or in none. The groups are solved one at a time,
  // each exactly (update_group()), under ncvreg's strong rule on the groups'
  // gradient norms; the L0 family swaps whole groups (group_swap_once()).

  // L0Learn's swaps over groups: rounds, each followed by a fresh solve,
  // until a round finds none or maxSwaps rounds have run. The smallest
  // convergence flag of those solves.
  int group_swap_phase(const PointPenalty& pp, double value, double prev) {
    int conv = 1;
    for (int round = 0; round < st_.maxSwaps; round++) {
      if (!group_swap_once(pp)) break;
      conv = std::min(conv, solve_grouped(pp, value, prev));
    }
    return conv;
  }

  // One round of swap_once() over groups: the first active group, in index
  // order, whose best replacement among the zero groups has a larger
  // gradient norm, once the active group is set to zero, than the group is
  // held at, (v + 2 l2) ||beta_i|| + l1 sqrt(K_i), and whose swap -- the
  // group zeroed, the replacement set to its own update -- lowers the
  // objective. A swap that does not is undone. False when no group swaps.
  bool group_swap_once(const PointPenalty& pp) {
    std::vector<vec> g;
    for (auto& s : solvers_) g.push_back(s->gradient(s->ld_.product(s->beta_)));
    for (uword i = 0; i < members_.size(); i++) {
      double held = group_norm(i);
      if (held == 0.0) continue;
      std::vector<vec> cols = group_columns(i);
      double best = -1.0;
      uword jstar = members_.size();
      for (uword j = 0; j < members_.size(); j++) {
        if (members_[j].empty() || group_norm(j) != 0.0) continue;
        double r = swapped_gradient(j, i, g, cols);
        if (r > best) {
          best = r;
          jstar = j;
        }
      }
      if (jstar == members_.size()) continue;
      double curvature = 1.0 + 2.0 * pp.l2;
      if (best <= curvature * held + pp.l1 * group_weight(i)) continue;
      double before = group_objective(pp);
      vec keptI = group_coefs(i);
      zero_group(i);
      update_group(jstar, pp);
      if (group_objective(pp) >= before) {
        zero_group(jstar);
        set_group(i, keptI);
        continue;
      }
      everGroup_[jstar] = 1;
      strongGroup_[jstar] = 1;
      for (auto& s : solvers_)
        s->absGrad_ = abs(s->gradient(s->ld_.product(s->beta_)));
      return true;
    }
    return false;
  }

  // The objective at the current coefficients: every condition's loss
  // beta'R_s beta - 2 z'beta, the coupling, and the group L0 family's
  // penalty, 2 lambda0 and 2 l1 sqrt(K_j) ||beta_j|| per group in the model
  // and 2 l2 ||beta||^2.
  double group_objective(const PointPenalty& pp) const {
    double total = coupling_term();
    for (const auto& s : solvers_) {
      vec Rsb = s->ld_.product(s->beta_);
      total += dot(s->beta_, Rsb) - 2.0 * dot(s->z_, s->beta_) +
               2.0 * pp.l2 * dot(s->beta_, s->beta_);
    }
    for (uword j = 0; j < members_.size(); j++) {
      double r = group_norm(j);
      if (r > 0.0)
        total += 2.0 * pp.lambda0 + 2.0 * pp.l1 * group_weight(j) * r;
    }
    return total;
  }

  // Group j's coefficients, one per condition carrying it.
  vec group_coefs(uword j) const {
    vec b(members_[j].size());
    for (uword m = 0; m < members_[j].size(); m++)
      b(m) = solvers_[members_[j][m].first]->beta_(members_[j][m].second);
    return b;
  }

  // Group j's coefficients set to `b`.
  void set_group(uword j, const vec& b) {
    for (uword m = 0; m < members_[j].size(); m++) {
      PathSolver<Ld>& s = *solvers_[members_[j][m].first];
      uword i = members_[j][m].second;
      double del = b(m) - s.beta_(i);
      if (del == 0.0) continue;
      s.beta_(i) = b(m);
      s.ld_.update(i, del);
    }
  }

  // The norm of group j's coefficients.
  double group_norm(uword j) const {
    double total = 0.0;
    for (const auto& m : members_[j]) {
      double b = solvers_[m.first]->beta_(m.second);
      total += b * b;
    }
    return std::sqrt(total);
  }

  // Group i's columns of R_s, one per condition carrying it (empty for one
  // that does not).
  std::vector<vec> group_columns(uword i) {
    std::vector<vec> cols(solvers_.size());
    for (const auto& m : members_[i])
      cols[m.first] = solvers_[m.first]->column(m.second);
    return cols;
  }

  // The gradient norm of zero group j were group i set to zero: each
  // condition's gradient at j plus i's pull on it there.
  double swapped_gradient(uword j, uword i, const std::vector<vec>& g,
                          const std::vector<vec>& cols) const {
    double total = 0.0;
    for (const auto& m : members_[j]) {
      double gj = g[m.first](m.second);
      for (const auto& n : members_[i])
        if (n.first == m.first)
          gj += cols[m.first](m.second) * solvers_[n.first]->beta_(n.second);
      total += gj * gj;
    }
    return std::sqrt(total);
  }

  // Group j's coefficients set to zero.
  void zero_group(uword j) {
    for (const auto& m : members_[j]) {
      PathSolver<Ld>& s = *solvers_[m.first];
      double b = s.beta_(m.second);
      if (b == 0.0) continue;
      s.beta_(m.second) = 0.0;
      s.ld_.update(m.second, -b);
    }
  }

  // A group's weight: the square root of its size.
  double group_weight(uword j) const {
    return std::sqrt(static_cast<double>(members_[j].size()));
  }

  // The norm of group j's gradient: at zero (`atZero`, its z-scores), else
  // at the last recorded point.
  double group_gradient(uword j, bool atZero) const {
    double total = 0.0;
    for (const auto& m : members_[j]) {
      const PathSolver<Ld>& s = *solvers_[m.first];
      double g = atZero ? s.z_(m.second) : s.absGrad_(m.second);
      total += g * g;
    }
    return std::sqrt(total);
  }

  // The largest path value at which zero group j leaves zero, as
  // PathSolver::entry_value() reads it for one coefficient: its gradient
  // norm over its weight, or for the L0 family the norm's excess over
  // lambda1's weight, squared, over twice the curvature.
  double group_entry(uword j, bool atZero) const {
    if (members_[j].empty()) return 0.0;
    double w = group_weight(j);
    double g = group_gradient(j, atZero);
    if (st_.pen == PEN_L0) {
      double excess = std::max(g - st_.lambda1 * w, 0.0);
      return excess * excess / (2.0 * (1.0 + 2.0 * st_.lambda2));
    }
    if (st_.alpha <= 0.0) return std::numeric_limits<double>::infinity();
    return g / (st_.alpha * w);
  }

  // The strong rule over the groups, then descent over the ever-active
  // groups until no coefficient moves by thr, a scan of the strong groups
  // not yet active, and once that adds none a scan of the rest. 1 when a
  // scan of the rest adds none within `maxiter` sweeps.
  int solve_grouped(const PointPenalty& pp, double value, double prev) {
    double cut = solvers_[0]->strong_cutoff(value, prev);
    for (uword j = 0; j < members_.size(); j++)
      if (!strongGroup_[j] && group_gradient(j, false) > cut * group_weight(j))
        strongGroup_[j] = 1;
    int budget = st_.maxiter;
    while (true) {
      while (true) {
        if (budget-- <= 0) return 0;
        double dlx = 0.0;
        for (uword j = 0; j < members_.size(); j++)
          if (everGroup_[j]) dlx = std::max(dlx, update_group(j, pp));
        cpp11::check_user_interrupt();
        if (dlx < st_.thr && !scan_groups(pp, true)) break;
      }
      if (!scan_groups(pp, false)) return 1;
    }
  }

  // Each group not yet active, among the strong ones (`strong`) or the rest,
  // updated in index order, joining the active and strong groups if it
  // leaves zero. Whether any did.
  bool scan_groups(const PointPenalty& pp, bool strong) {
    bool added = false;
    for (uword j = 0; j < members_.size(); j++) {
      if (everGroup_[j] || (strongGroup_[j] != 0) != strong) continue;
      if (update_group(j, pp) > 0.0) {
        everGroup_[j] = 1;
        strongGroup_[j] = 1;
        added = true;
      }
    }
    return added;
  }

  // Group j's exact update: the minimizer over its coefficients b of
  //   b'Hb - 2 t'b + 2 P(||b||),
  // H the conditions' curvatures plus the coupling between them and t their
  // partial correlations, found by majorize-minimize steps (H below h I,
  // each step the one-coefficient operator on a norm), a single step when
  // H is h I. Returns how far its coefficients moved.
  double update_group(uword j, const PointPenalty& pp) {
    const auto& mem = members_[j];
    uword K = mem.size();
    if (K == 0) return 0.0;
    double ridge = st_.pen == PEN_L0 ? 2.0 * pp.l2 : pp.l2;
    vec t(K), b0(K), a(K);
    mat C(K, K, fill::zeros);
    for (uword m = 0; m < K; m++) {
      PathSolver<Ld>& s = *solvers_[mem[m].first];
      uword i = mem[m].second;
      b0(m) = s.beta_(i);
      t(m) = s.z_(i) - s.ld_.row(i, s.beta_) + s.v_(i) * b0(m);
      a(m) = s.v_(i) + ridge;
      for (uword n = 0; n < K; n++)
        if (n != m) C(m, n) = strength_(mem[m].first, mem[n].first);
    }
    vec b = group_minimizer(t, a, C, b0, pp, group_weight(j));
    double move = 0.0;
    for (uword m = 0; m < K; m++) {
      if (b(m) == b0(m)) continue;
      PathSolver<Ld>& s = *solvers_[mem[m].first];
      uword i = mem[m].second;
      s.beta_(i) = b(m);
      s.ld_.update(i, b(m) - b0(m));
      if ((b0(m) == 0.0) != (b(m) == 0.0)) s.supportChanged_ = true;
      move = std::max(move, std::abs(b(m) - b0(m)));
    }
    return move;
  }

  // The majorize-minimize steps of update_group(), from b0: with H =
  // diag(a + row sums of C) - C and h = max(a + 2 row sums of C) above its
  // largest eigenvalue (Gershgorin), each step minimizes h b'b - 2 u'b +
  // 2 P(||b||) for u = t + (h I - H) b, whose solution is u's direction at
  // the one-coefficient operator's value on ||u||. The steps stop when no
  // coefficient moves by thr / 100, or at once when H is h I.
  vec group_minimizer(const vec& t, const vec& a, const mat& C, const vec& b0,
                      const PointPenalty& pp, double w) const {
    vec rows = sum(C, 1);
    mat H = diagmat(a + rows) - C;
    double h = max(a + 2.0 * rows);
    bool exact = accu(C) == 0.0 && max(a) - min(a) <= 1e-12 * h;
    PointPenalty scaled{pp.l1 * w, 0.0, pp.lambda0};
    vec b = b0;
    for (int it = 0; it < st_.maxiter; it++) {
      vec u = t + h * b - H * b;
      double nu = norm(u);
      double r = nu > 0.0
                     ? std::abs(apply_threshold(st_.pen, nu, h, scaled,
                                                st_.gamma))
                     : 0.0;
      vec next = nu > 0.0 ? vec((r / nu) * u) : vec(b.n_elem, fill::zeros);
      double move = max(abs(next - b));
      b = next;
      if (exact || move < st_.thr * 1e-2) break;
    }
    return b;
  }

  void finish(double value, const PointPenalty& pp, int conv) {
    for (uword k = 0; k < solvers_.size(); k++) {
      refresh_shift(k);
      solvers_[k]->finish_point(value, pp, conv);
    }
    coupling_.push_back(coupling_term());
  }

  // PathSolver::solve_strong() across the studies: each one's strong set
  // marked from its own gradient, then descent over every study's
  // ever-active set until no coefficient of any moves by thr, and the scans
  // of every study.
  int solve_strong(const PointPenalty& pp, double value, double prev) {
    for (auto& s : solvers_) {
      double cut = s->strong_cutoff(value, prev);
      for (uword j = 0; j < s->p_; j++)
        if (!s->strong_[j] && s->absGrad_(j) > cut) s->strong_[j] = 1;
    }
    int budget = st_.maxiter;
    while (true) {
      do {
        if (!descend(pp, budget)) return 0;
      } while (scan_all(pp, PathSolver<Ld>::IN_STRONG));
      if (!scan_all(pp, PathSolver<Ld>::OUTSIDE_STRONG)) return 1;
    }
  }

  bool descend(const PointPenalty& pp, int& budget) {
    while (budget > 0) {
      budget--;
      double dlx = 0.0;
      for (uword k = 0; k < solvers_.size(); k++) {
        refresh_shift(k);
        dlx = std::max(dlx, solvers_[k]->sweep(solvers_[k]->everList_, pp));
      }
      if (dlx < st_.thr) return true;
    }
    return false;
  }

  bool scan_all(const PointPenalty& pp, ScanSet which) {
    bool added = false;
    for (uword k = 0; k < solvers_.size(); k++) {
      refresh_shift(k);
      if (solvers_[k]->scan(pp, which)) added = true;
    }
    return added;
  }

  // PathSolver::solve_ranked() across the studies: each study sweeps its own
  // ranked set, restricted to its support once that has held, and at
  // convergence every study checks its zeros.
  int solve_ranked(const PointPenalty& pp) {
    uword K = solvers_.size();
    std::vector<std::vector<uword>> sets(K);
    std::vector<std::vector<unsigned char>> inSet(K);
    std::vector<int> stable(K, 0);
    for (uword k = 0; k < K; k++) {
      sets[k] = solvers_[k]->ranked_set();
      inSet[k].assign(solvers_[k]->p_, 0);
      for (uword j : sets[k]) inSet[k][j] = 1;
    }
    for (int budget = st_.maxiter; budget > 0; budget--) {
      double dlx = 0.0;
      for (uword k = 0; k < K; k++) {
        PathSolver<Ld>& s = *solvers_[k];
        refresh_shift(k);
        s.supportChanged_ = false;
        dlx = std::max(dlx, s.sweep(sets[k], pp));
        stable[k] = s.supportChanged_ ? 0 : stable[k] + 1;
        if (stable[k] == st_.control.stableSweeps)
          s.restrict_to_support(sets[k], inSet[k]);
      }
      if (dlx >= st_.thr) continue;
      if (!add_zeros(pp, sets, inSet, stable)) return 1;
    }
    return 0;
  }

  // Every study's zero coordinates checked; those that left zero join its
  // set. Whether any did.
  bool add_zeros(const PointPenalty& pp, std::vector<std::vector<uword>>& sets,
                 std::vector<std::vector<unsigned char>>& inSet,
                 std::vector<int>& stable) {
    bool any = false;
    for (uword k = 0; k < solvers_.size(); k++) {
      refresh_shift(k);
      std::vector<uword> added = solvers_[k]->check_zeros(pp);
      if (added.empty()) continue;
      any = true;
      stable[k] = 0;
      for (uword j : added) {
        if (inSet[k][j]) continue;
        sets[k].push_back(j);
        inSet[k][j] = 1;
      }
    }
    return any;
  }

  // L0Learn's swap rounds across the studies: each round gives every study
  // one swap, then a fresh solve of all of them.
  int swap_phase(const PointPenalty& pp, double value, double prev) {
    int conv = 1;
    for (int round = 0; round < st_.maxSwaps; round++) {
      bool swapped = false;
      for (uword k = 0; k < solvers_.size(); k++) {
        refresh_shift(k);
        if (solvers_[k]->swap_once(pp)) swapped = true;
      }
      if (!swapped) break;
      conv = std::min(conv, solve(pp, value, prev));
    }
    return conv;
  }
};

// ---- Registered entry points ---------------------------------------------

// The path against the LD matrix R, fitted as R_s = (1 - s) R + s I with s =
// `shrink`, which is never formed (DenseLd). `pathR` is the path, in lambda
// for LASSO, MCP and SCAD and in lambda0 for the L0 family, non-increasing,
// and for the L0 family may be empty for the default path. `control` holds
// the path constants (R's .rssPathControl()).
[[cpp11::register]]
cpp11::writable::list penalizedRssCpp(
    const doubles& zR,
    const doubles_matrix<>& ldR,
    double shrink,
    const doubles& pathR,
    const std::string& penaltyStr,
    const std::string& strategyStr,
    double gamma,
    double alpha,
    double lambda1,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps,
    const list& control) {

  PathSettings st = parse_settings(penaltyStr, strategyStr, gamma, alpha,
                                   lambda1, lambda2, thr, maxiter, maxSwaps,
                                   control);
  vec z = as_Col(zR);
  mat R = as_Mat(ldR);
  if (R.n_rows != z.n_elem || R.n_cols != z.n_elem)
    cpp11::stop("The LD matrix is %d x %d but there are %d z-scores.",
                (int)R.n_rows, (int)R.n_cols, (int)z.n_elem);
  vec path = as_Col(pathR);
  check_path(path, st.pen);
  if (adaptive_default(st, path)) {
    DenseLd scout(R, shrink, vec(z.n_elem, fill::zeros));
    PathSolver<DenseLd> first(st, scout, z);
    path = off_edge(first.run(path).value);
  }
  DenseLd ld(R, shrink, vec(z.n_elem, fill::zeros));
  PathSolver<DenseLd> solver(st, ld, z);
  return path_result(solver.run(path), z.n_elem);
}

// penalizedRssCpp() with the LD given by a factor F of R = F'F / (m - 1),
// shrunk to R_s = (1 - s) R + s I: the same path, updates and swaps, without
// forming R or R_s.
[[cpp11::register]]
cpp11::writable::list penalizedRssFactorCpp(
    const doubles& zR,
    const doubles_matrix<>& factorR,
    double shrink,
    const doubles& pathR,
    const std::string& penaltyStr,
    const std::string& strategyStr,
    double gamma,
    double alpha,
    double lambda1,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps,
    const list& control) {

  PathSettings st = parse_settings(penaltyStr, strategyStr, gamma, alpha,
                                   lambda1, lambda2, thr, maxiter, maxSwaps,
                                   control);
  vec z = as_Col(zR);
  mat F = as_Mat(factorR);
  if (F.n_cols != z.n_elem)
    cpp11::stop("The LD factor has %d columns but there are %d z-scores.",
                (int)F.n_cols, (int)z.n_elem);
  vec path = as_Col(pathR);
  check_path(path, st.pen);
  if (adaptive_default(st, path)) {
    FactorLd scout(F, shrink, vec(z.n_elem, fill::zeros));
    PathSolver<FactorLd> first(st, scout, z);
    path = off_edge(first.run(path).value);
  }
  FactorLd ld(F, shrink, vec(z.n_elem, fill::zeros));
  PathSolver<FactorLd> solver(st, ld, z);
  return path_result(solver.run(path), z.n_elem);
}

// ---- Coupled entry points --------------------------------------------------

// The studies' z vectors.
static std::vector<vec> coupled_z(const list& zList) {
  std::vector<vec> zs;
  for (R_xlen_t k = 0; k < zList.size(); k++)
    zs.push_back(as_Col(doubles(zList[k])));
  return zs;
}

// The studies' LD matrices or factors, each checked against its z-scores:
// a matrix square over them, a factor with a column for each.
static std::vector<mat> coupled_mats(const list& matList,
                                     const std::vector<vec>& zs,
                                     bool factor) {
  if (matList.size() != static_cast<R_xlen_t>(zs.size()))
    cpp11::stop("There are %d LD sources for %d studies.",
                (int)matList.size(), (int)zs.size());
  std::vector<mat> out;
  for (uword k = 0; k < zs.size(); k++) {
    mat M = as_Mat(doubles_matrix<>(matList[static_cast<R_xlen_t>(k)]));
    bool fits = factor ? M.n_cols == zs[k].n_elem
                       : M.n_rows == zs[k].n_elem && M.n_cols == zs[k].n_elem;
    if (!fits)
      cpp11::stop("Study %d's LD is %d x %d but it has %d z-scores.",
                  (int)k + 1, (int)M.n_rows, (int)M.n_cols,
                  (int)zs[k].n_elem);
    out.push_back(M);
  }
  return out;
}

// The pairs' strengths, K x K for K conditions: symmetric, finite and not
// negative. The diagonal is not read.
static mat coupled_strength(const doubles_matrix<>& strengthR, uword K) {
  mat W = as_Mat(strengthR);
  if (W.n_rows != K || W.n_cols != K)
    cpp11::stop("The coupling strengths are %d x %d for %d conditions.",
                (int)W.n_rows, (int)W.n_cols, (int)K);
  if (!W.is_finite() || W.min() < 0.0 || !W.is_symmetric())
    cpp11::stop("The coupling strengths must be symmetric, finite and not "
                "negative.");
  return W;
}

// The coupled path on `mats` held as Ld. The matrices outlive the solver,
// which reads them by reference.
template <class Ld>
static cpp11::writable::list run_coupled(const PathSettings& st,
                                         const std::vector<mat>& mats,
                                         const std::vector<vec>& zs,
                                         const std::vector<uvec>& idx,
                                         int nUnion, const mat& strength,
                                         bool grouped, double shrink,
                                         const vec& path) {
  std::vector<std::unique_ptr<Ld>> lds;
  for (uword k = 0; k < mats.size(); k++)
    lds.emplace_back(new Ld(mats[k], shrink, vec(zs[k].n_elem, fill::zeros)));
  CoupledPathSolver<Ld> solver(st, lds, zs, idx, nUnion, strength, grouped);
  solver.run(path);
  return solver.result();
}

// The path of K conditions coupled by `strengthR` (CoupledPathSolver), the
// K x K strengths c_kk' of their pairs (zero for a pair not coupled), each
// against its own LD matrix (DenseLd), shrunk to (1 - s) R + s I with s =
// `shrink`. `unionIdxList` gives each study's variants' 1-based indices into
// the union of the studies' variants, of size `nUnion`, in one allele frame.
// The path is shared, and must be given for LASSO, MCP and SCAD; the L0
// family's default is geometric from the joint pathMax. Returns each study's
// path (`studies`, as penalizedRssCpp() returns one), the path values and the
// coupling term at each.
[[cpp11::register]]
cpp11::writable::list penalizedRssCoupledCpp(
    const list& zList,
    const list& ldList,
    const list& unionIdxList,
    int nUnion,
    const doubles_matrix<>& strengthR,
    bool grouped,
    double shrink,
    const doubles& pathR,
    const std::string& penaltyStr,
    const std::string& strategyStr,
    double gamma,
    double alpha,
    double lambda1,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps,
    const list& control) {

  PathSettings st = parse_settings(penaltyStr, strategyStr, gamma, alpha,
                                   lambda1, lambda2, thr, maxiter, maxSwaps,
                                   control);
  std::vector<vec> zs = coupled_z(zList);
  std::vector<uvec> idx = union_index(unionIdxList, zs, nUnion, "z-scores");
  std::vector<mat> Rs = coupled_mats(ldList, zs, false);
  vec path = as_Col(pathR);
  check_path(path, st.pen);
  mat strength = coupled_strength(strengthR, zs.size());
  return run_coupled<DenseLd>(st, Rs, zs, idx, nUnion, strength, grouped,
                              shrink, path);
}

// penalizedRssCoupledCpp() with each study's LD given by a factor F of
// R = F'F / (m - 1) (FactorLd).
[[cpp11::register]]
cpp11::writable::list penalizedRssCoupledFactorCpp(
    const list& zList,
    const list& factorList,
    const list& unionIdxList,
    int nUnion,
    const doubles_matrix<>& strengthR,
    bool grouped,
    double shrink,
    const doubles& pathR,
    const std::string& penaltyStr,
    const std::string& strategyStr,
    double gamma,
    double alpha,
    double lambda1,
    double lambda2,
    double thr,
    int maxiter,
    int maxSwaps,
    const list& control) {

  PathSettings st = parse_settings(penaltyStr, strategyStr, gamma, alpha,
                                   lambda1, lambda2, thr, maxiter, maxSwaps,
                                   control);
  std::vector<vec> zs = coupled_z(zList);
  std::vector<uvec> idx = union_index(unionIdxList, zs, nUnion, "z-scores");
  std::vector<mat> Fs = coupled_mats(factorList, zs, true);
  vec path = as_Col(pathR);
  check_path(path, st.pen);
  mat strength = coupled_strength(strengthR, zs.size());
  return run_coupled<FactorLd>(st, Fs, zs, idx, nUnion, strength, grouped,
                               shrink, path);
}
