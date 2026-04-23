# OTTERS Lassosum Selector Fix

## Summary

This change makes OTTERS lassosum selection source-aware again.

- PLINK1 `bed/bim/fam` inputs use old genotype-backed `lassosum` pseudovalidation.
- PLINK2 sketch `pgen/pvar/psam` inputs use sketch-genotype pseudovalidation on restored `U`.
- Inputs with no genotype-backed source still fall back to `min(fbeta)`, with an explicit warning that the fallback is not old-OTTERS-compatible.
- The OTTERS wrapper also stops double-scaling lassosum input before it reaches the low-level solver.

## Example 206

Fixture `chr1_206088859_208088859__ENSG00000123843` isolates the regression.

- old saved vs old direct published `lassosum`: Pearson `1.0`, `0` opposite-sign variants
- corrected-scaling + `min(fbeta)`: Pearson about `0.360`, `1309` opposite-sign variants

Old published lassosum selected `s=0.2, lambda=1e-4`. Corrected-scaling plus `min(fbeta)` still selected `s=1, lambda=1e-4`.

Implication:

- this is not a grid-definition problem
- both candidates are already on the old grid
- the regression comes from changing the selection rule over the same grid

For this fixture, old pseudovalidation is a degenerate tie case, so old OTTERS selects the first maximum. Replacing that with `min(fbeta)` moves the selected model to the opposite end of the same grid.

## Example 161

Fixture `chr1_161023868_163023868__ENSG00000162745` is the informative selector case.

On the common matched set of `4298` variants:

- exact PLINK1 bfile pseudovalidation best candidate: `soft_lambda=0.041050213`
- sketch-genotype pseudovalidation best candidate: `soft_lambda=0.021788613`
- score agreement remains high:
  - Pearson `0.9971`
  - Spearman `0.9412`
  - max absolute difference `0.0474`

When the sketch path is forced to use `sd_bfile`, the winning candidate returns to `soft_lambda=0.041050213`.

Implication:

- this is not a variant-overlap problem
- this is a predictive normalization problem inside pseudovalidation
- restored sketch `U` recovers the score surface closely
- sketch-derived `sd` is not identical to old `lassosum:::sd.bfile()`

So PLINK1 is the exact old-parity path, while PLINK2 sketch is the right source-aware selector for sketch inputs but remains an approximation relative to old PLINK1 behavior.

## What Is Fixed

## `R/regularized_regression.R`

- `lassosum_rss_weights()` no longer defaults every OTTERS run to `min(fbeta)`.
- It now selects by source:
  - PLINK1: published `lassosum` pseudovalidation
  - PLINK2 sketch: restored-`U` pseudovalidation
  - no genotype-backed source: warning + `min(fbeta)` fallback
- It preserves old first-max tie behavior for pseudovalidation selectors.
- It passes lassosum input in correlation units only once before the low-level solver converts to `z / sqrt(n)`.

What actually failed:

- old OTTERS-compatible selection was replaced with `min(fbeta)`
- the OTTERS wrapper also double-scaled lassosum input

Why it worked earlier:

- old OTTERS used published `lassosum.pipeline(..., test.bfile = bfile)` followed by `pseudovalidate(out)`

Why it failed here:

- the refactor changed the selector and the scaling contract
- fixture `206` exposed the selector regression clearly because it is a tied pseudovalidation table where first-max and `min(fbeta)` pick different candidates

Classification:

- selector replacement: real behavioral regression
- double scaling: real implementation bug
- fixture `206` tie case: edge-case trigger that exposed the regression cleanly

Compatibility:

- exact old behavior is restored for PLINK1 genotype-backed OTTERS runs
- PLINK2 sketch inputs now use a source-aware selector instead of the generic fallback
- pure precomputed-LD inputs remain on the fallback path and now say so explicitly

## `R/otters.R`

- `otters_weights()` now forwards aligned `chrom/pos/A1/A2` variant metadata to lassosum when available.
- This restores the variant-alignment information needed for source-aware selection.

Classification:

- wrapper or interface drift

## Tests

Focused validation for this change:

- `tests/testthat/test_rr_lassosum.R`
- `tests/testthat/test_rr_dispatch.R`
- `tests/testthat/test_otters.R`

Local evidence paths used during debugging:

- `temp_reference/otters_regression/lassosum_oldR_direct_206/`
- `temp_reference/otters_regression/lassosum_forensics_206/`
- `temp_reference/otters_regression/pseudovalidation_sketch_vs_bfile/`

## Narrow Compatibility Note

This change only alters the OTTERS lassosum selection path.

- It restores exact old behavior when a legacy PLINK1 genotype source is available.
- It makes PLINK2 sketch inputs use a sketch-genotype selector instead of the old generic fallback.
- It does not change PRS-CS or SDPR behavior.
