# OTTERS Lassosum LD-Quadratic Selector Fix

## Summary

This PR fixes the OTTERS lassosum regression by removing the old `min(fbeta)`
selector from the default OTTERS path and replacing it with the LD-quadratic
selector

```text
score(beta) = (c^T beta) / sqrt(beta^T R beta)
```

where:

- `c` is the aligned summary-statistics correlation vector
- `R` is the supplied LD correlation matrix
- `beta` is one candidate on the lassosum `(s, lambda)` path

The important conclusion from the follow-up diagnostics is that genotype is not
fundamentally required for lassosum selection once the selector is written in
this LD form. The remaining sketch issue is about how `R` is constructed or
standardized, not about the selector formula itself.

## What Failed

Old OTTERS did not select lassosum models by `min(fbeta)`. It fit the beta path
and then used lassosum pseudovalidation to choose the final `(s, lambda)`.

The refactor changed that selector to `min(fbeta)`, and the OTTERS wrapper also
double-scaled the lassosum input before it reached the low-level solver.

Those two changes were enough to move the selected model to a very different
part of the same grid.

## Example 206

Fixture `chr1_206088859_208088859__ENSG00000123843` isolates the selector bug.

- old saved vs old direct published `lassosum`: Pearson `1.0`, `0` opposite-sign variants
- corrected-scaling + `min(fbeta)`: Pearson about `0.360`, `1309` opposite-sign variants

On this fixture:

- published lassosum selected `s = 0.2`, `lambda = 1e-4`
- `min(fbeta)` selected `s = 1`, `lambda = 1e-4`

This is not a grid-definition problem. Both candidates are already on the old
grid. The regression comes from changing the selection rule over the same
candidate path.

## Mathematical Rationale

Old pseudovalidation can be written as:

```text
scaled_beta = beta / sd
pred        = X * scaled_beta
score       = (c^T beta) / sqrt(Var(pred))
```

After centering and standardizing the columns of `X` by the same per-variant
scale, this becomes:

```text
score(beta) = (c^T beta) / sqrt(beta^T R beta)
```

So the selector can be evaluated directly from summary-statistics correlation
and LD, without using genotype explicitly. That is the selector implemented in
this PR.

## Validation

We ran the LD-quadratic score directly against genotype-backed pseudovalidation
on the same candidate matrix.

### PLINK1-derived LD

The LD-quadratic score matches genotype pseudovalidation essentially exactly.

- Fixture `161`:
  - genotype best: `soft_lambda=0.041050213`
  - LD-quadratic PLINK1 best: `soft_lambda=0.041050213`
  - Pearson `0.9999999`
  - same best candidate `TRUE`
- Fixture `206`:
  - genotype best: `soft_lambda=0.029906976`
  - LD-quadratic PLINK1 best: `soft_lambda=0.029906976`
  - Pearson `1.0000000`
  - same best candidate `TRUE`

This validates the selector formula itself.

### Sketch-derived LD

The same LD-quadratic selector applied to the current sketch-derived `R` gives:

- Fixture `161`:
  - sketch-LD quadratic best: `soft_lambda=0.021788613`
  - Pearson vs genotype pseudovalidation `0.9971044`
  - same best candidate `FALSE`
- Fixture `206`:
  - sketch-LD quadratic best: `soft_lambda=0.029906976`
  - Pearson vs genotype pseudovalidation `0.9999503`
  - same best candidate `TRUE`

That result is consistent with the sketch sample-matrix diagnostic: the same
winner shift appears there too. So the remaining mismatch is not between
"sample-matrix pseudovalidation" and "quadratic LD scoring". It is between the
current sketch representation and the PLINK1/genotype-backed standardized LD
path.

## What This PR Changes

## `R/regularized_regression.R`

- fixes the OTTERS lassosum scaling contract so correlation input is only
  converted once before the low-level solver
- removes genotype-format-specific selector dispatch from
  `lassosum_rss_weights()`
- makes the default selector `ld_quadratic`
- keeps `min(fbeta)` only as an explicit debug option
- preserves first-max tie behavior for equal selector scores

## `R/otters.R`

- passes correlation-scale statistics into lassosum explicitly via
  `stat$cor` / `stat$z`
- removes the temporary genotype-source / variant-metadata plumbing that was
  only needed for the earlier compatibility patch

## Scope

This PR only changes the OTTERS lassosum selection path.

- It does not change PRS-CS or SDPR behavior.
- It does not claim the current sketch-derived LD is already identical to the
  PLINK1/genotype-backed path.
- It does claim that the selector should be implemented in LD form, and that
  `min(fbeta)` was the wrong default for OTTERS parity.

## Evidence Paths

- `temp_reference/otters_regression/lassosum_oldR_direct_206/`
- `temp_reference/otters_regression/lassosum_forensics_206/`
- `temp_reference/otters_regression/pseudovalidation_sketch_vs_bfile/`
