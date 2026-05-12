# TWAS Pipeline Analysis

## End-to-end data flow

```
Weight RDS files ──→ load_twas_weights() ──┐
GWAS TSV (tabix) ──→ harmonize_gwas()  ──┤──→ harmonize_twas() ──→ twas_analysis() ──→ results
LD metadata TSV  ──→ load_LD_matrix()  ──┘         │                      │
                                                    │                      ↓
                                              match_ref_panel()      z = (wᵀz) / √(wᵀRw)
                                              (allele QC engine)
```

## The pipeline has 3 layers

### Layer 1: `twas_pipeline()` — orchestrator (R/twas.R:326-632)

The main entry point. It:
1. Optionally filters molecular events via `event_filters`
2. Calls `harmonize_twas()` to load + QC everything
3. For each gene: picks the best model by CV R² (via `pick_best_model()`)
4. For each gene-context-study triple: calls `twas_analysis()` to compute z-scores
5. Optionally runs MR analysis if SuSiE credible sets exist and p-value is small enough
6. Merges CV metrics with TWAS results into a single table
7. Applies `apply_method_fallback()` for NA/Inf z-scores
8. Optionally formats output for cTWAS via `format_twas_data()`

### Layer 2: `harmonize_twas()` — data harmonization (R/twas.R:30-226)

The most complex function. It:
1. **Groups contexts by genomic position** — contexts whose weight variants are within 5kb of each other get clustered so they share a single LD query region
2. **Loads LD once** for the combined region via `load_LD_matrix()`
3. For each GWAS study: calls `harmonize_gwas()` which loads data via tabix, standardizes columns, then calls `match_ref_panel()` to align alleles to the LD reference
4. For each context: calls `match_ref_panel()` again to align the weight matrix to LD, then optionally `adjust_susie_weights()` to recalculate SuSiE weights for the variant subset
5. Scales weights by `sqrt(variance)` where variance comes from the LD ref panel
6. Extracts a per-gene LD submatrix for the intersection of GWAS + weight variants

**The three-way intersection is the core issue** — every variant must be present in all three sources (weights, GWAS, LD) to participate, and alleles must be harmonized across all three.

### Layer 3: `twas_analysis()` / `twas_z()` — the actual computation (R/twas.R:747-773, 653-667)

This is tiny and straightforward:
- Subset everything to shared variants
- For each method: `z_twas = (wᵀ * z_gwas) / sqrt(wᵀ * R * w)`
- Return z-score and chi-squared p-value

## Where the complexity lives

| Area | What happens | Why it's complex |
|------|-------------|------------------|
| **Data loading** | `load_twas_weights()` merges multivariate (mnm_rs) + univariate weights from multiple RDS files | Deeply nested list structures, context name cleaning, weight alignment across methods |
| **Allele QC** | `match_ref_panel()` handles exact match, sign flip, strand flip, INDEL matching | 6+ boolean flags, inner join on (chrom, pos) then allele matching logic |
| **Context grouping** | `group_contexts_by_region()` clusters contexts by variant position overlap | Hierarchical clustering + IRanges interval merging |
| **SuSiE adjustment** | `adjust_susie_weights()` recalculates weights from log Bayes factors when variants are dropped | Re-derives alpha from LBF, recomputes posterior means |
| **Model selection** | Best model chosen by CV R², with fallback if z-score is NA/Inf | Two-pass: first in `pick_best_model()`, then `apply_method_fallback()` |

## The xqtl-protocol adds another layer on top

The `twas_ctwas.ipynb` notebook wraps `twas_pipeline()` with:
1. **Weight validation** — file size checks, tryCatch on readRDS (now ported into pecotmr)
2. **Batch loading** — `batch_load_twas_weights()` splits genes by memory
3. **`update_twas_method()`** — a SECOND method fallback pass after `twas_pipeline()` returns (the in-pecotmr `apply_method_fallback()` was ported from this)
4. **Context name merging** — `merge_context_names()` strips region suffixes from weight names
5. **cTWAS assembly** — ~200 lines loading LD again, re-running `harmonize_gwas()`, calling `trim_ctwas_variants()`, `assemble_region_data()`, then `est_param()`, `screen_regions()`, `finemap_regions()`

## Key observations

1. **`match_ref_panel()` is called 3+ times per gene** — once for GWAS vs LD, once for weights vs LD, and once inside `adjust_susie_weights()`. Each call re-parses variant IDs and re-joins.

2. **LD is loaded once but variant subsetting happens repeatedly** — `harmonize_twas()` loads the full region LD, then extracts per-gene submatrices. But in the cTWAS cell, LD is loaded *again* for the same chromosome.

3. **The data structure is deeply nested** — `twas_weights_data[[molecular_id]]$weights[[context]]` requires `get_nested_element()` and `find_data()` utilities just to navigate. The same molecular_id/context/study keys appear at multiple levels.

4. **Model selection happens in two places** — `pick_best_model()` inside `twas_pipeline()` selects the best method before TWAS, then `apply_method_fallback()` fixes bad selections after TWAS. xqtl-protocol had a third pass via `update_twas_method()`.

5. **Weight scaling is tightly coupled to LD loading** — `scaled_weights = weights * sqrt(variance)` depends on allele frequencies and sample size from the LD ref panel, computed during `load_LD_matrix()`.

6. **The actual TWAS math is ~10 lines** — the other ~600 lines in twas.R are data loading, harmonization, QC, and result formatting.
