#!/usr/bin/env Rscript
# =============================================================================
# Prepare PLINK 2.0 logistic GWAS results for METAL
# - Project: Genetic determinants of Lewy body pathology
# - Version: R; data.table, dplyr
# - Last Updated: 16-SEP-2026
#
# Overview
# - Keeps additive-model rows without a PLINK error code, biallelic SNVs with a
#   consistent effect allele and non-missing statistics; writes a METAL input file.
# - Run once per cohort (e.g. GP2, NIHNBB), then run 01_METAL_meta-analysis.txt.
#
# Usage
#   Rscript 00_prepare_METAL_input.R <plink2 .glm.logistic(.hybrid) file> <output .tab>
#   e.g. Rscript 00_prepare_METAL_input.R LB_sumstats.txt GP2_METAL.tab
#   Output columns: SNP (CHROM:POS) effect_allele noneffect_allele beta se Pvalue N EAF
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript 00_prepare_METAL_input.R <plink2 glm file> <output .tab>", call. = FALSE)
in_file <- args[1]; out_file <- args[2]
if (!file.exists(in_file)) stop("file not found: ", in_file, call. = FALSE)

# Read PLINK2 glm results
GP2 <- fread(in_file)

# Clean column name
setnames(GP2, "#CHROM", "CHROM", skip_absent = TRUE)
needed <- c("CHROM", "POS", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "LOG(OR)_SE", "P", "A1_FREQ", "ERRCODE")
missing <- setdiff(needed, names(GP2))
if (length(missing)) stop("input is missing column(s): ", paste(missing, collapse = ", "),
                          " (run plink2 --glm with cols=+a1freq,+err)", call. = FALSE)
n0 <- nrow(GP2)

# Keep only additive model results
GP2a <- GP2 %>%
  filter(TEST == "ADD")

# Keep rows with no PLINK error code (ERRCODE is NA, "", or "." when there is no error)
GP2a <- GP2a %>%
  filter(is.na(ERRCODE) | ERRCODE == "" | ERRCODE == ".")

# Check A1 consistency and define effect / non-effect allele
GP2a <- GP2a %>%
  mutate(
    effect_allele = case_when(
      A1 == ALT ~ ALT,
      A1 == REF ~ REF,
      TRUE ~ NA_character_
    ),
    noneffect_allele = case_when(
      A1 == ALT ~ REF,
      A1 == REF ~ ALT,
      TRUE ~ NA_character_
    )
  )

# Remove allele inconsistencies
GP2a <- GP2a %>%
  filter(!is.na(effect_allele), !is.na(noneffect_allele))

# Remove indels / non-SNVs
GP2a <- GP2a %>%
  filter(
    nchar(REF) == 1,
    nchar(ALT) == 1,
    nchar(effect_allele) == 1,
    nchar(noneffect_allele) == 1
  )

# Remove rows with missing or impossible values
GP2a <- GP2a %>%
  filter(
    !is.na(OR),
    !is.na(`LOG(OR)_SE`),
    !is.na(P),
    !is.na(OBS_CT),
    !is.na(A1_FREQ),
    OR > 0
  )

# Create beta, effect-allele frequency and SNP ID
GP2a <- GP2a %>%
  mutate(
    beta = log(OR),
    EAF = A1_FREQ,
    SNP = paste(CHROM, POS, sep = ":")
  )

# Keep METAL-ready columns
metal <- GP2a %>%
  select(
    SNP,
    effect_allele,
    noneffect_allele,
    beta,
    se = `LOG(OR)_SE`,
    Pvalue = P,
    N = OBS_CT,
    EAF
  )

n_dup <- sum(duplicated(metal$SNP))
if (n_dup)
  warning(n_dup, " SNP IDs (CHROM:POS) occur more than once (multi-allelic sites); ",
          "METAL uses marker names only, so check these before meta-analysis", call. = FALSE)

# Write output
fwrite(metal, out_file, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("%s: %d rows read, %d written to %s\n", in_file, n0, nrow(metal), out_file))
