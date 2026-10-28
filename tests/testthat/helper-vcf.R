# Small VCF and SeqArray GDS files for the VCF / GDS readers' tests
# (R/GenotypeDosageSeed.R, R/genotypeIo.R): written, bgzipped and
# tabix-indexed at test time, so every record's expected dosage is known.

# A bgzipped, tabix-indexed VCF of `records` (a data frame of CHROM, POS, ID,
# REF, ALT and one FORMAT column per sample, holding the FORMAT values as
# text). `format` names the FORMAT keys the records carry, `ds` the Number
# FORMAT/DS is declared with (NULL declares no DS), and `gt = FALSE` leaves
# FORMAT/GT undeclared. Returns the .vcf.gz path.
.writeTestVcf <- function(
    records,
    samples,
    format = "GT",
    ds = NULL,
    gt = TRUE,
    contigs = unique(records$CHROM),
    dir = tempfile("vcf_")
) {
    dir.create(dir, showWarnings = FALSE)
    header <- c(
        "##fileformat=VCFv4.2",
        paste0("##contig=<ID=", contigs, ">"),
        if (gt) {
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        },
        if (!is.null(ds)) {
            paste0(
                "##FORMAT=<ID=DS,Number=",
                ds,
                ',Type=Float,Description="Dosage">'
            )
        },
        paste(
            c(
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                samples
            ),
            collapse = "\t"
        )
    )
    fields <- c(
        list(
            records$CHROM,
            records$POS,
            records$ID,
            records$REF,
            records$ALT,
            ".",
            "PASS",
            ".",
            paste(format, collapse = ":")
        ),
        unname(as.list(records[samples]))
    )
    plain <- file.path(dir, "test.vcf")
    writeLines(c(header, do.call(paste, c(fields, sep = "\t"))), plain)
    gz <- Rsamtools::bgzip(plain, paste0(plain, ".gz"), overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "vcf")
    gz
}

# The records the readers' tests share: a SNP; a deletion whose REF spans the
# SNP after it; two records at one position, the halves of a split
# multi-allelic site, both named rs4; and one record with two ALT alleles.
# Every way a record can overlap another, which is what a range read must
# see through. GT and DS (one value per ALT allele) agree, apart from s1's
# first record, whose DS carries an imputed fraction.
.vcfOverlapRecords <- function() {
    data.frame(
        CHROM = "22",
        POS = c(100L, 200L, 202L, 300L, 300L, 400L),
        ID = c(".", "rs2", "rs3", "rs4", "rs4", "rs6"),
        REF = c("A", "ACGTA", "G", "C", "C", "A"),
        ALT = c("G", "A", "T", "T", "G", "C,G"),
        s1 = c("0/0:0.1", "0/1:1", "1|1:2", "0/1:1", "0/0:0", "1/2:1,1"),
        s2 = c("0/1:0.9", "0/0:0", "0|1:1", "0/0:0", "0/1:1", "0/2:0,1"),
        s3 = c("1/1:1.8", "./.:.", "0|0:0", "0/1:1", "1/1:2", "0|1:1,0"),
        stringsAsFactors = FALSE
    )
}

# The dosage of each of .vcfOverlapRecords()'s seven variants (rows: the
# two-ALT record splits in two) from its GT calls and from its DS values.
.vcfOverlapGt <- function() {
    rbind(
        c(0, 1, 2),
        c(1, 0, NA),
        c(2, 1, 0),
        c(1, 0, 1),
        c(0, 1, 2),
        c(1, 0, 1),
        c(1, 1, 0)
    )
}

# Keep only the GT subfield of each sample column.
.vcfGtOnly <- function(records, samples = c("s1", "s2", "s3")) {
    gt <- lapply(records[samples], sub, pattern = ":.*$", replacement = "")
    `[<-`(records, samples, value = gt)
}

# A SeqArray GDS converted from the VCF at `vcf`, importing FORMAT/DS when
# `ds` is TRUE.
.writeTestSeqArray <- function(vcf, ds = FALSE) {
    out <- tempfile(fileext = ".gds")
    SeqArray::seqVCF2GDS(
        vcf,
        out,
        fmt.import = if (ds) "DS" else character(),
        verbose = FALSE
    )
    out
}
