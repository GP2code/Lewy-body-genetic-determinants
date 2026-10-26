#!/usr/bin/env bash
# =============================================================================
# Logistic-regression GWAS (PLINK 2.0)
# - Project: Genetic determinants of Lewy body pathology
# - Version: bash; PLINK v2.0
# - Last Updated: 16-SEP-2026
#
# Overview
# - Per-chromosome logistic GWAS with Firth fallback, adjusted for age, sex and PC1-5
#   (covariates variance-standardised), minor allele count >= 20.
# - Run once per cohort (GP2, Mayo, NIHNBB) and per comparison:
#     GWAS #1: autopsy-confirmed LBD cases vs controls
#     GWAS #2: cortical vs subcortical LBD cases
#   Only the PHENO column of the covariate file differs between comparisons.
#
# Usage
#   bash 00_logistic_GWAS_plink2.sh <workdir>
#     <workdir>/chr{1..22}.pgen/.pvar/.psam  genotypes
#     <workdir>/COVAR.txt                    FID IID PHENO AGE sex PC1..PC5 (PHENO binary, PLINK coding)
#   Output: <workdir>/chr{1..22}.GWAS_RESULTS.PHENO.glm.logistic.hybrid
#
# =============================================================================
set -euo pipefail

WORKDIR="$1"
PLINK="plink2"      
THREADS=32
MEMORY_MB=120000

for chr in {1..22}; do
  echo "Running chromosome ${chr}..."

  "$PLINK" \
    --pfile "${WORKDIR}/chr${chr}" \
    --pheno "${WORKDIR}/COVAR.txt" \
    --pheno-name PHENO \
    --covar "${WORKDIR}/COVAR.txt" \
    --covar-name AGE,sex,PC1,PC2,PC3,PC4,PC5 \
    --covar-variance-standardize \
    --glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
    --mac 20 \
    --memory "$MEMORY_MB" \
    --threads "$THREADS" \
    --out "${WORKDIR}/chr${chr}.GWAS_RESULTS"
done
