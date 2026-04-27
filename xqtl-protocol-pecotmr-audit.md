# xqtl-protocol pecotmr Usage Audit

## 1. API Mismatches (Calls Expected to Fail)

### A. `load_quantile_twas_weights` does not exist in pecotmr

**File:** `code/pecotmr_integration/twas_ctwas.ipynb` (quantile_twas cell)

```r
twas_weights_results[[gene_db]] <- load_quantile_twas_weights(
    weight_db_files = weight_dbs, tau_values = tau_values,
    between_cluster = 0.8, num_intervals = 3)
```

This function is not defined anywhere in pecotmr. It was likely planned as part of the quantile TWAS feature but never implemented (or was removed). This cell will fail outright.

### B. `twas_pipeline` called with nonexistent `quantile_twas` parameter

**File:** `code/pecotmr_integration/twas_ctwas.ipynb` (quantile_twas cell)

```r
twas_results_db <- twas_pipeline(..., quantile_twas = TRUE, ...)
```

`twas_pipeline()` at `R/twas.R:290` has no `quantile_twas` parameter and no `...` in its signature. This will error with "unused argument."

### C. `rss_analysis_pipeline` called with renamed parameter `stochastic_ld_sample`

**File:** `code/mnm_analysis/mnm_methods/rss_analysis.ipynb` (univariate_rss cell)

```r
rss_analysis_pipeline(..., stochastic_ld_sample = ${stochastic_ld_sample}, ...)
```

This parameter was renamed to `sketch_samples` in the current pecotmr API (`R/univariate_pipeline.R:214`). The function has no `...`, so this will error with "unused argument."

### D. `load_multitrait_tensorqtl_sumstat` called with wrong parameter name

**File:** `code/multivariate_genome/MASH/mash_preprocessing.ipynb` (random_null_tensorqtl_1 cell)

```r
pecotmr::load_multitrait_tensorqtl_sumstat(
    phenotype_path = phenotype_path, ..., na_remove = T/F)
```

Two problems:
- First parameter is named `sumstats_paths` in pecotmr (`R/mash_wrapper.R:141`), not `phenotype_path`
- The parameter `na_remove` was renamed to `nan_remove` (`R/mash_wrapper.R:143`)

Both will cause "unused argument" errors since there's no `...`.

### E. `mash_ran_null_sample` - typo and removed parameters

**File:** `code/multivariate_genome/MASH/mash_preprocessing.ipynb` (random_null_tensorqtl_1 cell)

```r
pecotmr::mash_ran_null_sample(dat, n_random, n_null,
    expected_ncondition, exclude_condition, z_only = TRUE, seed = ...)
```

Three problems:
- Function name is `mash_rand_null_sample` (with a "d") -- `R/mash_wrapper.R:568`
- `expected_ncondition` parameter no longer exists
- `z_only` parameter no longer exists

The current signature is `mash_rand_null_sample(dat, n_random, n_null, exclude_condition, seed = NULL)`.

### F. `get_ctwas_meta_data` is deprecated

**File:** `code/pecotmr_integration/twas_ctwas.ipynb` (ctwas_1 and ctwas_3 cells)

Used extensively. Still works but emits deprecation warnings and will eventually be removed. The replacement is `ld_loader()` per `R/ctwas_wrapper.R:58`.

---

## 2. Safety/Sanity Checks That Could Move to pecotmr

### A. Weight file pre-validation

**File:** twas_ctwas.ipynb (twas cell, quantile_twas cell)

Before calling `load_twas_weights()`, xqtl-protocol:
- Checks `file.size(file) > 200` (non-trivial file)
- Wraps `readRDS(file)` in `tryCatch` to filter corrupt files
- Validates nested structure (`twas_variant_names` key exists)
- Filters out NULL/empty results

**Recommendation:** `load_twas_weights()` should do this validation internally -- skip files that are too small, corrupt, or structurally invalid, rather than requiring every caller to implement the same filter.

### B. NA/Inf z-score filtering after TWAS

**File:** twas_ctwas.ipynb (ctwas_1 cell)

```r
z_gene[[study]] <- z_gene[[study]][
    !is.na(z_gene[[study]]$z) & !is.infinite(z_gene[[study]]$z) &
    z_gene[[study]]$id %in% names(weight_list[[study]]),]
```

**Recommendation:** `twas_pipeline()` or the TWAS z-score computation itself should guarantee clean output. Downstream consumers shouldn't need to re-filter.

### C. Duplicate LD variant removal

**File:** twas_ctwas.ipynb (ctwas_1 cell)

```r
dup_idx <- which(duplicated(LD_list$LD_variants))
if (length(dup_idx) >= 1) LD_list$LD_matrix <- LD_list$LD_matrix[-dup_idx, -dup_idx]
```

**Recommendation:** `load_LD_matrix()` should handle this internally. Duplicate variants in the LD matrix are a data integrity issue that the loader should resolve before returning.

### D. GWAS sample size validation

**File:** twas_ctwas.ipynb (ctwas_1 cell)

```r
if(length(z_snp[['sample_size']][[study]]!=1) | z_snp[['sample_size']][[study]] <= 0) {
    stop("Please check sample size provided for ", study, " at --gwas_meta_data. ")
}
```

**Recommendation:** Could be validated inside `harmonize_gwas()` or the GWAS metadata loading step in pecotmr.

### E. chr prefix normalization

Scattered across multiple locations:
- ctwas_1: `ifelse(grepl("^chr", snp_map$id), snp_map$id, paste0("chr", snp_map$id))`
- mnm_postprocessing: `if(any(grepl("chr", qtl_all_var))) add_chr_prefix(gwas_all_var) else gsub("chr", "", gwas_all_var)`
- mash_preprocessing: retry with `gsub("chr", "", region)` on failure

**Recommendation:** pecotmr already has `normalize_variant_id()` and internal `strip_chr_prefix()`, but variant ID harmonization at the "chr" level should be consistently handled in all loading functions rather than requiring callers to do it.

### F. Genomic region overlap detection

**File:** SuSiE_enloc.ipynb (susie_coloc cell)

Manual region parsing and overlap checking:
```r
split_region <- unlist(strsplit(region, "_"))
block_chrom <- as.numeric(split_region[1] %>% gsub("chr","",.))
block_start <- ...
if (gene_region$chrom == block_chrom &&
    (gene_region$start <= block_end | gene_region$end >= block_start))
```

**Recommendation:** pecotmr has `parse_region()` and `region_to_df()` but lacks a simple `regions_overlap(a, b)` utility. This pattern is repeated enough to justify one.

---

## 3. Generalizable Pipeline Logic Worth Moving to pecotmr

### A. TWAS method selection with fallback (high value)

**File:** twas_ctwas.ipynb -- `update_twas_method()` function

This ~40-line function handles a real problem: the "best" TWAS method (by cross-validation) sometimes produces NA/Inf results for a specific GWAS. It falls back to the next-best method by rsq. This logic is not xQTL-specific -- it applies to any TWAS analysis.

**What it does:** For each gene-context-GWAS group, if the selected method yielded invalid z/p-values, pick the best alternative method that has valid results and meets the rsq threshold.

### B. TWAS-to-cTWAS region assembly orchestration (high value)

**File:** twas_ctwas.ipynb (ctwas_1 cell)

The entire workflow of:
1. Loading per-region TWAS results
2. Trimming variants via `trim_ctwas_variants()`
3. Getting chromosome-wide LD variant info
4. Harmonizing GWAS via `harmonize_gwas()`
5. Re-computing TWAS z-scores when variants are trimmed (calling `twas_analysis()` with fresh LD)
6. Assembling into cTWAS region data via `assemble_region_data()`

This is ~200 lines of orchestration that any TWAS-to-cTWAS pipeline would need. It currently depends on a few ctwas package functions but the overall flow is generalizable.

### C. cTWAS fine-mapping with LD diagnosis and recovery (high value)

**File:** twas_ctwas.ipynb (ctwas_3 cell)

The workflow of:
1. Screen regions -> fine-map -> diagnose LD mismatch -> identify problematic genes -> re-fine-map without LD -> merge boundary regions

This is a robust, production-tested recipe for dealing with real-world LD mismatches. It's not specific to xQTL data at all.

### D. GWAS metadata loading and per-study LD caching (medium value)

**Files:** rss_analysis.ipynb, twas_ctwas.ipynb

The pattern of:
- Reading a GWAS metadata TSV with study_id, chrom, file_path, sample_size columns
- Mapping studies to per-region file paths
- Caching LD matrices by study to avoid re-loading

This is boilerplate that every multi-study RSS analysis repeats. The Python `load_regional_rss_data()` function in rss_analysis.ipynb does something similar -- it could inform an R equivalent.

### E. MASH data batching and merging (medium value)

**File:** mash_preprocessing.ipynb (susie_to_mash_1, susie_to_mash_2 cells)

The pipeline of:
1. Processing regions in chunks (`per_chunk = 100`)
2. Extracting strong/random/null z-score matrices per region
3. Renaming rownames with region IDs for uniqueness
4. Merging across regions with `merge_mash_data()`
5. Filtering invalid entries with `filter_invalid_summary_stat()`
6. Computing `ZtZ = t(Z) %*% Z / n`

The batching/merging/filtering logic is generalizable. pecotmr already has `merge_mash_data()` and `filter_invalid_summary_stat()`, but the end-to-end orchestration (batch -> merge -> filter -> ZtZ) could be a single pipeline function.

### F. QTL-GWAS overlap analysis pipeline (medium value, Python)

**File:** mnm_postprocessing.ipynb (overlap_qtl_gwas cells, Python)

Loads QTL and GWAS metadata, groups by chromosome, checks region overlap, intersects variant lists with chr-prefix harmonization. This is applicable to any pairwise colocalization setup, not just xQTL.

### G. Variant feature engineering (lower value, Python)

**File:** gems_pipeline.py

Parsing "chr:pos:ref:alt" variant IDs into structured fields and classifying as SNP/indel/insertion/deletion. Simple but broadly useful. pecotmr has `parse_variant_id()` already but the SNP/indel classification could be an addition.

---

## Summary

| Category | Count | Severity |
|----------|-------|----------|
| Calls that will fail | 5 (A-E) | Blocking |
| Deprecated but still working | 1 (F) | Warning |
| Sanity checks to absorb | 6 (A-F) | Robustness |
| Generalizable pipeline logic | 7 (A-G) | Feature opportunities |

The most impactful items are the **quantile TWAS breakage** (the entire quantile_twas workflow is dead -- both `load_quantile_twas_weights` and the `quantile_twas` param to `twas_pipeline` don't exist), the **`stochastic_ld_sample` rename**, and the **MASH parameter name changes**. For generalizable logic, the **TWAS method fallback** and **cTWAS assembly/diagnosis pipelines** are the highest-value candidates to move into pecotmr.
