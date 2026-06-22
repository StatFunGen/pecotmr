context("Collection-level accessors")

# ===========================================================================
# Shared fixture
# ===========================================================================

.ca_makeTopLoci <- function(n = 3, withCs = TRUE) {
  tl <- data.frame(
    variant_id     = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    chrom          = rep("1", n),
    pos            = as.integer(100 * seq_len(n)),
    A1             = rep("G", n),
    A2             = rep("A", n),
    N              = rep(1000, n),
    MAF            = rep(0.1, n),
    marginal_beta  = rep(0.1, n),
    marginal_se    = rep(0.05, n),
    marginal_z     = rep(2.0, n),
    marginal_p     = rep(0.05, n),
    pip            = seq(0.9, by = -0.1, length.out = n),
    posterior_mean = rep(0.05, n),
    posterior_sd   = rep(0.02, n),
    stringsAsFactors = FALSE)
  if (withCs) tl$cs_95 <- paste0("susie_", c(1L, 1L, 0L)[seq_len(n)])
  tl
}

.ca_makeFmEntry <- function(n = 3) {
  FineMappingEntry(
    variantIds = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    susieFit = list(payload = sprintf("fit_n=%d", n)),
    topLoci    = .ca_makeTopLoci(n))
}

# ===========================================================================
# QtlFineMappingResult collection accessors
# ===========================================================================

