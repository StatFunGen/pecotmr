context("Collection-level accessors")

# ===========================================================================
# Shared fixture
# ===========================================================================

.ca_makeTopLoci <- function(n = 3, withCs = TRUE) {
  tl <- data.frame(
    variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    pip        = seq(0.9, by = -0.1, length.out = n),
    stringsAsFactors = FALSE)
  if (withCs) tl$cs <- c(1L, 1L, 0L)[seq_len(n)]
  tl
}

.ca_makeFmEntry <- function(n = 3) {
  FineMappingEntry(
    variantIds = paste0("chr1:", 100 * seq_len(n), ":A:G"),
    trimmedFit = list(payload = sprintf("fit_n=%d", n)),
    topLoci    = .ca_makeTopLoci(n))
}

# ===========================================================================
# QtlFineMappingResult collection accessors
# ===========================================================================

