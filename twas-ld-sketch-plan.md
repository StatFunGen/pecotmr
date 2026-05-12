# TWAS Pipeline Refactoring: LD Sketch + SVD

## Context

The current TWAS pipeline loads precomputed LD correlation matrices R from independent LD blocks (or computes R from genotype files), then computes `z_twas = (wᵀz) / √(wᵀRw)`. This requires aligning gene windows with LD block boundaries, concatenating blocks, deduplicating boundary variants, and extracting per-gene submatrices — the dominant source of complexity in twas.R.

The refactoring replaces R with an SVD of "LD sketch" genotypes — random projections of the reference panel that preserve LD structure. Each gene independently defines its own window, loads the sketch for that window, and computes the TWAS statistic via:

```
z_twas = (wᵀz) / ‖Λ^(1/2) Vᵀw‖     where Λ = D²/(n_sketch - 1)
```

This eliminates the LD block alignment problem entirely.

---

## Implementation Steps

### Step 1: Thin wrapper `load_ld_sketch()` in R/LD.R

Add after `load_LD_from_genotype()` (~line 480). This is a thin wrapper around the existing `load_LD_matrix()` with `return_genotype = TRUE` — it reuses all existing genotype loading, variant ID normalization, ref_panel construction (allele frequencies, variance), and duplicate removal. The only new logic is standardizing X and computing SVD.

```r
#' @export
load_ld_sketch <- function(ld_meta_file_path, region, n_sample = NULL)
```

Logic:
1. Call `load_LD_matrix(ld_meta_file_path, region, return_genotype = TRUE, n_sample = n_sample)` — this returns the genotype matrix X, ref_panel (with allele_freq), LD_variants, etc.
2. Standardize X using HWE-based formula with `p = ref_panel$allele_freq`:
   - Center: `X_c = sweep(X, 2, 2*p)` (subtract theoretical mean `2p`)
   - Scale: `X_std = sweep(X_c, 2, sqrt(2*p*(1-p)), "/")` (divide by theoretical SD)
   - Drop zero-variance columns where `p == 0` or `p == 1` (and corresponding ref_panel rows / variant_ids)
3. `safe_svd(X_std, tol = 0)` — full SVD, no truncation
4. Return: `list(V = svd$v, D = svd$d, n_sketch = nrow(X), ref_panel = result$ref_panel, variant_ids = result$LD_variants)`

Key: `n_sketch` = nrow(X_sketch) used for Λ normalization. `n_sample` = original panel size, passed through to `load_LD_matrix()` for variance computation.

### Step 2: Modify `twas_z()` in R/twas.R (lines 653-667)

New signature:
```r
twas_z <- function(weights, z, R = NULL, X = NULL, V = NULL, D = NULL, n_sketch = NULL)
```

Add SVD branch before existing R/X branches:
- If V, D, n_sketch provided: `Lambda = D^2 / (n_sketch - 1)`, `denom = sum(Lambda * (crossprod(V, weights))^2)`, `zscore = stat / sqrt(denom)`
- Existing R and X paths remain unchanged

### Step 3: Modify `twas_analysis()` in R/twas.R (lines 747-773)

New signature:
```r
twas_analysis <- function(weights_matrix, gwas_sumstats_db, LD_matrix = NULL,
                          extract_variants_objs, V = NULL, D = NULL,
                          n_sketch = NULL, ld_variant_ids = NULL)
```

SVD path: subset V rows via `match(valid_variants_objs, ld_variant_ids)`, pass V_subset/D/n_sketch to `twas_z()`. Existing LD_matrix path remains.

### Step 4: Refactor `harmonize_twas()` in R/twas.R (lines 30-226)

**Remove:**
- `group_contexts_by_region()` inner function (lines 33-91)
- Region-wide `load_LD_matrix()` call (lines 106-115)
- Per-gene LD submatrix extraction (lines 211-220)
- Context-grouping loop structure

**New per-gene loop:**
```
for (molecular_id in molecular_ids):
  1. Collect variant positions across all contexts → build gene window region
  2. load_ld_sketch(ld_meta_file_path, gene_region, n_sample) — internally calls load_LD_matrix(return_genotype=TRUE), standardizes, SVDs
  3. Use sketch$variant_ids as reference for allele QC
  
  for (study in gwas_studies):
    4. harmonize_gwas(gwas_file, gene_region, sketch$variant_ids, ...)
    
    for (context in contexts):
      5. match_ref_panel(weights, sketch$variant_ids, ...)
      6. adjust_susie_weights() if needed
      7. Scale weights by sqrt(ref_panel$variance)
  
  8. Store sketch SVD: mol_res[["svd_V"]], [["svd_D"]], [["n_sketch"]], [["ld_variant_ids"]]
```

Return value changes: `mol_res[["LD"]]` (p×p matrix) → `mol_res[["svd_V"]]`, `[["svd_D"]]`, `[["n_sketch"]]`, `[["ld_variant_ids"]]`

### Step 5: Update `twas_pipeline()` in R/twas.R

- Lines 529-532: Pass SVD components to `twas_analysis()` instead of `[["LD"]]`
- Line 569: Clear SVD components instead of LD matrix after use

### Step 6: Tests

New tests in a new file `tests/testthat/test_twas_sketch.R`:

1. **Mathematical equivalence test**: Generate X, compute R=cor(X), compute SVD(X_std). Run `twas_z()` both ways, assert identical z-scores within floating-point tolerance.
2. **`load_ld_sketch()` unit test**: Mock `load_LD_matrix()` with known X. Verify SVD standardization, zero-variance column removal, ref_panel passthrough.
3. **`twas_analysis()` SVD path**: Verify partial variant overlap works correctly.
4. **End-to-end mock test**: Full pipeline with SVD path produces correct results.

### Step 7: NAMESPACE and documentation

Run `devtools::document()` after adding roxygen to `load_ld_sketch()` and updating param docs on modified functions.

---

## Files to modify

| File | Changes |
|------|---------|
| `R/LD.R` | Add `load_ld_sketch()` (thin wrapper around existing `load_LD_matrix(return_genotype=TRUE)` + SVD) |
| `R/twas.R` | Modify `twas_z()`, `twas_analysis()`, `harmonize_twas()`, `twas_pipeline()`. Remove `group_contexts_by_region()` |
| `tests/testthat/test_twas_sketch.R` | New test file |
| `NAMESPACE` | Auto-generated via roxygen |

## Files NOT modified

| File | Reason |
|------|--------|
| `R/LD.R` (load_LD_matrix, load_LD_from_genotype, load_LD_from_blocks) | Reused as-is; `load_ld_sketch()` calls `load_LD_matrix()` |
| `R/allele_qc.R` | match_ref_panel unchanged |
| `R/file_utils.R` | load_genotype_region unchanged |
| `R/misc.R` | safe_svd, compute_LD unchanged |
| `R/susie_wrapper.R` | adjust_susie_weights unchanged |
| `R/ctwas_wrapper.R` | cTWAS deferred |
| xqtl-protocol repo | Not changed in this PR |

## Existing functions reused (not reimplemented)

- `load_LD_matrix()` — R/LD.R:265 — main LD loader; `load_ld_sketch()` calls this with `return_genotype=TRUE`
- `load_LD_from_genotype()` — R/LD.R:408 — called by `load_LD_matrix()` for genotype sources; handles genotype loading, variant ID normalization, ref_panel with allele frequencies/variance, .afreq sidecar
- `safe_svd()` — R/misc.R:156 — SVD with tolerance filtering
- `match_ref_panel()` — R/allele_qc.R:26 — allele QC (reference = sketch variants)
- `ensure_chr_match()` — R/misc.R:480 — chr prefix harmonization

## Verification

1. Run `pixi run Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test_twas_sketch.R")'`
2. Run existing tests to verify no regressions: `pixi run Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test_twas_method_fallback.R")'`
3. Mathematical equivalence: the key test generates a genotype matrix, computes TWAS z-scores via both R and SVD paths, asserts they match within 1e-10.
