# Shared S-LDSC fixture generators for test_sldscWrapper.R (helper unit tests)
# and test_sldscPostprocessingPipeline.R (integration tests).
#
# Fixture convention:
#   - 2 chromosomes (1, 2), 50 SNPs each -> 100 total
#   - 2 target annotations: "annot_A" (binary), "annot_B" (continuous)
#   - baseline annotations (baselineLD_0 ..) in joint run
#   - 10 jackknife blocks
#   - Polyfun appends "_0" to target annotation names in .results

# Create a single .annot.gz file for one chromosome.
# Real polyfun .annot.gz files have CHR, SNP, BP, CM + annotation columns only
# (no MAF/A1/A2 -- those come from the .frq / PLINK files).
.make_annot_gz <- function(dir, chrom, nSnps = 50) {
  df <- data.frame(
    CHR = chrom,
    SNP = paste0("rs", (chrom - 1L) * 100L + seq_len(nSnps)),
    BP  = seq_len(nSnps) * 1000L,
    CM  = seq_len(nSnps) * 0.01,
    annot_A = sample(c(0L, 1L), nSnps, replace = TRUE),
    annot_B = rnorm(nSnps, 2, 0.5),
    stringsAsFactors = FALSE
  )
  path <- file.path(dir, sprintf("target.%d.annot.gz", chrom))
  gz <- gzfile(path, "wb")
  vroom::vroom_write(df, gz, delim = "\t")
  close(gz)
  invisible(df)
}

# Create a PLINK .frq file for one chromosome
.make_frq <- function(dir, chrom, plinkName ="ref_chr", nSnps = 50) {
  df <- data.frame(
    CHR = chrom,
    SNP = paste0("rs", (chrom - 1L) * 100L + seq_len(nSnps)),
    A1  = "A",
    A2  = "G",
    MAF = runif(nSnps, 0.01, 0.49),
    NCHROBS = 200L,
    stringsAsFactors = FALSE
  )
  path <- file.path(dir, sprintf("%s%d.frq", plinkName, chrom))
  vroom::vroom_write(df, path, delim = "\t")
  invisible(df)
}

# Create the three polyfun output files (.results, .log, .part_delete)
# for a single-target run. Real polyfun output includes baseline categories
# even in single-target mode, so we add 2 dummy baseline categories.
.make_polyfun_single <- function(dir, prefix, target_name, nBlocks = 10,
                                 h2g = 0.3, tau = 1e-7, enrichment = 2.5,
                                 n_baseline = 2) {
  target_cat <- paste0(target_name, "_0")
  baseline_cats <- paste0("baselineLD_", seq_len(n_baseline) - 1L)
  all_cats <- c(target_cat, baseline_cats)
  n_cats <- length(all_cats)

  taus_all <- c(tau, rep(1e-8, n_baseline))
  enrichments_all <- c(enrichment, rep(1.0, n_baseline))

  results <- data.frame(
    Category                = all_cats,
    Coefficient             = taus_all,
    Coefficient_std_error   = abs(taus_all) * 0.3,
    Enrichment              = enrichments_all,
    Enrichment_std_error    = enrichments_all * 0.2,
    Enrichment_p            = rep(0.01, n_cats),
    `Prop._h2`              = c(0.15, rep(0.425, n_baseline)),
    `Prop._SNPs`            = c(0.06, rep(0.47, n_baseline)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  vroom::vroom_write(results, paste0(prefix, ".results"), delim = "\t")

  writeLines(c(
    "Analysis started at 2024-01-01",
    sprintf("Total Observed scale h2: %g (0.05)", h2g),
    "Analysis finished"
  ), paste0(prefix, ".log"))

  blocks <- matrix(rnorm(nBlocks * n_cats,
                         mean = rep(taus_all, each = nBlocks),
                         sd = abs(rep(taus_all, each = nBlocks)) * 0.5),
                   nrow = nBlocks, ncol = n_cats)
  colnames(blocks) <- all_cats
  vroom::vroom_write(as.data.frame(blocks), paste0(prefix, ".part_delete"), delim = "\t")
  invisible(NULL)
}

# Create polyfun output files for a joint run (target + baseline annotations)
.make_polyfun_joint <- function(dir, prefix, target_names,
                                n_baseline = 3, nBlocks = 10, h2g = 0.3) {
  target_cats <- paste0(target_names, "_0")
  baseline_cats <- paste0("baselineLD_", seq_len(n_baseline) - 1L)
  all_cats <- c(target_cats, baseline_cats)
  n_cats <- length(all_cats)

  taus <- c(rep(1e-7, length(target_cats)), rep(1e-8, n_baseline))
  enrichments <- c(rep(2.0, length(target_cats)), rep(1.0, n_baseline))

  results <- data.frame(
    Category                = all_cats,
    Coefficient             = taus,
    Coefficient_std_error   = abs(taus) * 0.3,
    Enrichment              = enrichments,
    Enrichment_std_error    = enrichments * 0.2,
    Enrichment_p            = rep(0.05, n_cats),
    `Prop._h2`              = rep(1 / n_cats, n_cats),
    `Prop._SNPs`            = rep(1 / n_cats, n_cats),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  vroom::vroom_write(results, paste0(prefix, ".results"), delim = "\t")

  writeLines(c(
    "Analysis started at 2024-01-01",
    sprintf("Total Observed scale h2: %g (0.05)", h2g),
    "Analysis finished"
  ), paste0(prefix, ".log"))

  blocks <- matrix(rnorm(nBlocks * n_cats, mean = rep(taus, each = nBlocks),
                         sd = abs(rep(taus, each = nBlocks)) * 0.5),
                   nrow = nBlocks, ncol = n_cats)
  colnames(blocks) <- all_cats
  vroom::vroom_write(as.data.frame(blocks), paste0(prefix, ".part_delete"), delim = "\t")
  invisible(NULL)
}

# Build a complete fixture directory for the full pipeline
.make_sldsc_fixtures <- function(envir = parent.frame()) {
  base_dir <- withr::local_tempdir(.local_envir = envir)

  anno_dir <- file.path(base_dir, "annot")
  frq_dir  <- file.path(base_dir, "frq")
  out_dir  <- file.path(base_dir, "output")
  dir.create(anno_dir)
  dir.create(frq_dir)
  dir.create(out_dir)

  plink_name <- "ref_chr"

  # Annotation + freq files for 2 chromosomes
  for (chr in 1:2) {
    .make_annot_gz(anno_dir, chr)
    .make_frq(frq_dir, chr, plinkName =plink_name)
  }

  targets <- c("annot_A", "annot_B")

  # Single-target runs: 2 targets x 2 traits
  for (trait in c("traitX", "traitY")) {
    for (i in seq_along(targets)) {
      pref <- file.path(out_dir, sprintf("%s_single_%s", trait, targets[i]))
      .make_polyfun_single(out_dir, pref, targets[i], h2g = 0.3 + (i - 1) * 0.05)
    }
  }

  # Joint runs: 1 per trait
  for (trait in c("traitX", "traitY")) {
    pref <- file.path(out_dir, sprintf("%s_joint", trait))
    .make_polyfun_joint(out_dir, pref, targets, h2g = 0.3)
  }

  list(
    base_dir   = base_dir,
    anno_dir   = anno_dir,
    frq_dir    = frq_dir,
    out_dir    = out_dir,
    plinkName =plink_name,
    targets    = targets,
    trait_names = c("traitX", "traitY")
  )
}

# =============================================================================
# In-memory builders for SldscData / compute / pipeline tests (no file I/O)
# =============================================================================

# Build one readSldscTrait-shaped run list (the in-memory shape the pipeline
# and standardizeSldscTrait consume).
.sldscMkRun <- function(cats, h2g = 0.3, nBlocks = 10L,
                        tau = 1e-7, enrichment = 2.0, enrichmentP = 0.01) {
  n <- length(cats)
  list(
    categories   = cats,
    tau          = setNames(rep(tau, n), cats),
    tauSe        = setNames(rep(abs(tau) * 0.3, n), cats),
    enrichment   = setNames(rep(enrichment, n), cats),
    enrichmentSe = setNames(rep(enrichment * 0.2, n), cats),
    enrichmentP  = setNames(rep(enrichmentP, n), cats),
    propH2       = setNames(rep(0.2, n), cats),
    propSnps     = setNames(rep(0.1, n), cats),
    h2g          = h2g,
    tauBlocks    = matrix(rep(tau, nBlocks * n), nBlocks, n,
                          dimnames = list(NULL, cats)),
    nBlocks      = nBlocks)
}

# Build a small, valid in-memory SldscData: 2 traits, 2 target annotations
# (annot_A binary, annot_B continuous), single + optional joint runs.
.sldscMkData <- function(withJoint = TRUE, withFrq = TRUE,
                         traitNames = c("traitX", "traitY")) {
  annot <- data.frame(
    CHR     = c(1, 1, 1, 2, 2, 2),
    SNP     = paste0("rs", 1:6),
    annot_A = c(1, 0, 1, 0, 1, 0),
    annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2.0),
    stringsAsFactors = FALSE)
  frq <- if (withFrq)
    data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
               MAF = rep(0.2, 6), stringsAsFactors = FALSE)
  else NULL
  mkTrait <- function() {
    tr <- list(single = list(.sldscMkRun(c("annot_A_0", "baselineLD_0")),
                             .sldscMkRun(c("annot_B_0", "baselineLD_0"))))
    tr$joint <- if (withJoint)
      .sldscMkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")) else NULL
    tr
  }
  traits <- setNames(lapply(traitNames, function(.) mkTrait()), traitNames)
  SldscData(annot = annot, frq = frq, traits = traits)
}
