# =============================================================================
# Ordinal GWAS of Lewy body pathology stage (OrdinalGWAS.jl)
# - Project: Genetic determinants of Lewy body pathology
# - Version: Julia; OrdinalGWAS.jl, GLM.jl
# - Last Updated: 16-SEP-2026
#
# Overview
# - Proportional-odds (logit link) GWAS comparing brainstem-, limbic- and neocortical-
#   predominant LBD cases, with likelihood-ratio tests. One chromosome per call.
# - Response `stage` is 1/2/3 = brainstem / limbic / neocortical.
#
# Usage
#   julia 01_ordinal_GWAS.jl <chr> <dir> <tag> <model> <bed_pattern> <formula>
#   e.g. julia 01_ordinal_GWAS.jl 1 /path/to/data LBD stage_model "chr{}" "stage ~ AGE + sex + PC1 + PC2 + PC3 + PC4 + PC5"
#     <dir>/<bed_pattern with {} = chr>.bed/.bim/.fam  PLINK 1 genotypes
#     <dir>/ALIGNED/<tag>_<model>_COVAR_aligned.csv  covariates aligned to the genotype sample order
#     <dir>/ALIGNED/<tag>_<model>_rowinds.txt        row indices (one per line) of the samples to analyse
#   Output: <dir>/ALIGNED/<tag>_<model>_chr<chr>.pval.txt
# =============================================================================

using OrdinalGWAS, GLM

if length(ARGS) != 6
    println(stderr, "usage: julia 01_ordinal_GWAS.jl <chr> <dir> <tag> <model> <bed_pattern> <formula>")
    exit(1)
end

c, DIR, TAG, MODEL, BEDPAT, FORM = ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5], ARGS[6]
ALIGNED = joinpath(DIR, "ALIGNED")
idx = parse.(Int, readlines(joinpath(ALIGNED, "$(TAG)_$(MODEL)_rowinds.txt")))
f   = eval(Meta.parse("@formula(" * FORM * ")"))

println("chr$c  model=$MODEL  n=", length(idx), "  ", FORM); flush(stdout)
ordinalgwas(f, joinpath(ALIGNED, "$(TAG)_$(MODEL)_COVAR_aligned.csv"),
            joinpath(DIR, replace(BEDPAT, "{}" => c));
            covrowinds = idx, geneticrowinds = idx,
            test = :lrt, link = LogitLink(),
            pvalfile = joinpath(ALIGNED, "$(TAG)_$(MODEL)_chr$(c).pval.txt"))
println("chr$c DONE"); flush(stdout)
