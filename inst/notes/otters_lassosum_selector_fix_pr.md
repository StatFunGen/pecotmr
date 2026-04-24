# OTTERS Lassosum Selector Fix

## Summary

This change fixes the immediate OTTERS lassosum regression, but the updated diagnostics also sharpen the longer-term conclusion: genotype is not fundamentally required for lassosum selection once the selector is written in LD form.

- The regression was caused by replacing old pseudovalidation-based model selection with `min(fbeta)`.
- The compatibility patch in this PR restores old behavior conservatively:
  - PLINK1 `bed/bim/fam` inputs use old genotype-backed `lassosum` pseudovalidation.
  - PLINK2 sketch `pgen/pvar/psam` inputs use sketch-side pseudovalidation on restored `U`.
  - inputs with no genotype-backed source still fall back to `min(fbeta)`, with an explicit warning that the fallback is not old-OTTERS-compatible.
- The OTTERS wrapper also stops double-scaling lassosum input before it reaches the low-level solver.
- New diagnostics after the patch show that the selector itself can be expressed as an LD quadratic form:
  `score(beta) = (c^T beta) / sqrt(beta^T R beta)`.
  On PLINK1-derived LD, that LD-only score reproduces old genotype-backed pseudovalidation almost exactly.

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

So PLINK1 is the exact old-parity oracle in this PR, but it is no longer the conceptual end-state. The real conclusion is:

- genotype is not fundamentally needed for lassosum selection
- the LD quadratic selector is the mathematically correct formulation
- the remaining mismatch is specific to how the sketch-side `R` / scaling is constructed, not to the selector formula itself

## LD-Only Follow-up

After the compatibility patch, we ran a direct genotype-vs-LD pseudovalidation diagnostic using the same candidate matrix and the same aligned `cor` vector.

For each candidate `beta_k`, old pseudovalidation computes:

```text
scaled_beta = beta / sd
pred        = X * scaled_beta
score       = (c^T beta) / sqrt(Var(pred))
```

After centering and standardizing the matrix columns by the same `sd`, this becomes:

```text
score(beta) = (c^T beta) / sqrt(beta^T R beta)
```

where `R` is the standardized LD correlation matrix.

### PLINK1-derived LD

The LD quadratic score matches genotype-backed pseudovalidation essentially exactly.

- Fixture `161`:
  - genotype best: `soft_lambda=0.041050213`, score `0.3922803981`
  - LD-only PLINK1 best: `soft_lambda=0.041050213`, score `0.3922601179`
  - Pearson `0.9999999`, same best candidate `TRUE`
- Fixture `206`:
  - genotype best: `soft_lambda=0.029906976`, score `0.3895932657`
  - LD-only PLINK1 best: `soft_lambda=0.029906976`, score `0.3897370312`
  - Pearson `1.0000000`, same best candidate `TRUE`

This validates the LD-only quadratic selector. Genotype was only one way the old code evaluated the score; it is not inherently required once the standardized `R` is available.

### Sketch-derived LD

The same LD-only formulation,
`score(beta) = (c^T beta) / sqrt(beta^T R_sketch beta)`,
does **not** yet reproduce the genotype selector on the informative `161` fixture.

- Fixture `161`:
  - genotype best: `soft_lambda=0.041050213`
  - LD-only sketch best: `soft_lambda=0.021788613`
  - Pearson `0.9971044`, same best candidate `FALSE`
- Fixture `206`:
  - genotype best: `soft_lambda=0.029906976`
  - LD-only sketch best: `soft_lambda=0.029906976`
  - Pearson `0.9999503`, same best candidate `TRUE`

So the unresolved problem is no longer the selector formula. It is the sketch-side construction of a standardized `R` that matches the PLINK1/genotype-backed path.

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

Important clarification from the new diagnostics:

- this PR should not be read as claiming genotype is necessary for lassosum
- it should be read as a compatibility patch that restores the old selector while the cleaner LD-only selector is validated
- the LD-only selector is now validated on PLINK1-derived LD; the remaining work is to make the sketch-derived LD path produce the same standardized `R`

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
- It makes PLINK2 sketch inputs use a sketch-side selector instead of the old generic fallback.
- It does not change PRS-CS or SDPR behavior.

The new diagnostics imply that a future cleanup should replace the current genotype/sketch split with one format-generic LD-only selector once the sketch-side `R` construction is brought into line with the PLINK1 path.
