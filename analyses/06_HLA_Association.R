#!/usr/bin/env Rscript
# =============================================================================
# HLA allele association analysis (GP2 WGS HLA calls)
# - Project: Genetic determinants of Lewy body pathology
# - Version: R; data.table, fastglm, parallel
# - Last Updated: 16-SEP-2026
#
# Usage
#   Rscript 06_HLA_Association.R --samples samples.txt --manifest manifest.csv \
#     --hla-calls hla_calls.csv --hla-stats hla_call_stats.csv --pcs pcs.eigenvec \
#     [--related pcs.related] --out-dir outputs [--n-pcs 5] [--min-reads-mean 90]
#     [--min-carrier-freq 0.05] [--min-control-age 21] [--kinship 0.0884] [--cores N]
#   --samples    FID and IID columns, no header
#   --manifest   GP2ID, GP2sampleID, PD_PHENO (1 case, 0 control), qc_sex (1/2), AGE_ANALYSIS, wgs
#   --hla-calls  GP2sampleID and one dosage column (0/1/2) per allele named HLA-<allele>
#   --hla-stats  allele (without "HLA-") and reads_mean
#   --pcs        PLINK .eigenvec (IID, PC1, PC2, ...); --related: IID1, IID2, KINSHIP
#   Outputs: with_age/hla_results_with_age.csv, with_age/aoo_results.csv,
#            without_age/hla_results_without_age.csv (note column explains rows without a p-value)
#
# Overview
# For every HLA allele passing QC, fits
#   1. case-control logistic regression, with and without age as a covariate,
#      under five genotype codings (see MODEL_CODINGS below)
#   2. age-at-onset linear regression in cases (additive coding)
# using fastglm, adjusting for sex and the first --n-pcs principal components.
#
# Sample QC 
#   in --samples and the manifest, has HLA calls and PCs, passed WGS QC,
#   case/control status and sex known. 
#
# Allele QC
#   reads_mean >= --min-reads-mean, carrier frequency > --min-carrier-freq,
#   and at least one copy in both cases and controls. Individuals with a
#   missing call are excluded for that allele only.
#
# Interpretation: beta/OR refer to allele_term = 1 in each coding, so for the inverse
# codings (recessive_100, dominant_110) "Risk" means the allele is protective.
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fastglm)
  library(parallel)
})

MODEL_CODINGS <- list(          # genotype dosage 0 / 1 / 2  ->  allele_term
  additive      = c(0, 1, 2),
  recessive_001 = c(0, 0, 1),
  recessive_100 = c(1, 0, 0),
  dominant_011  = c(0, 1, 1),
  dominant_110  = c(1, 1, 0)
)
CLASSICAL_GENES <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")


# ---- command line -----------------------------------------------------------
OPTIONS <- list(          # name = default (NULL = required, NA = optional)
  samples = NULL, manifest = NULL, hla_calls = NULL, hla_stats = NULL, pcs = NULL,
  related = NA, out_dir = NULL,
  n_pcs = 5, min_reads_mean = 90, min_carrier_freq = 0.05,
  min_control_age = 21, kinship = 0.0884,
  cores = max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
)

usage <- function() {
  lines <- readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  rule <- grep("^# =+$", lines)
  cat(sub("^# ?", "", lines[(rule[1] + 1):(rule[2] - 1)]), sep = "\n")
}

parse_args <- function(args) {
  if (any(args %in% c("-h", "--help"))) { usage(); quit(status = 0) }
  opt <- OPTIONS
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, " (see --help)", call. = FALSE)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[1])
    if (!key %in% names(OPTIONS)) stop("unknown option: --", kv[1], " (see --help)", call. = FALSE)
    if (length(kv) > 1) {
      val <- paste(kv[-1], collapse = "=")
    } else {
      if (i == length(args)) stop("option --", kv[1], " needs a value", call. = FALSE)
      i <- i + 1; val <- args[i]
    }
    opt[[key]] <- if (is.numeric(OPTIONS[[key]])) as.numeric(val) else val
    if (is.numeric(OPTIONS[[key]]) && is.na(opt[[key]]))
      stop("option --", kv[1], " must be a number, got: ", val, call. = FALSE)
    i <- i + 1
  }
  missing <- names(opt)[vapply(opt, is.null, logical(1))]
  if (length(missing))
    stop("missing required option(s): ", paste0("--", gsub("_", "-", missing), collapse = ", "),
         " (see --help)", call. = FALSE)
  for (f in c("samples", "manifest", "hla_calls", "hla_stats", "pcs"))
    if (!file.exists(opt[[f]])) stop("file not found (--", gsub("_", "-", f), "): ", opt[[f]], call. = FALSE)
  if (!is.na(opt$related) && !file.exists(opt$related))
    stop("file not found (--related): ", opt$related, call. = FALSE)
  if (opt$n_pcs < 0 || opt$n_pcs != round(opt$n_pcs)) stop("--n-pcs must be a whole number >= 0", call. = FALSE)
  opt$n_pcs <- as.integer(opt$n_pcs)
  opt$cores <- if (.Platform$OS.type == "windows") 1L else max(1L, as.integer(opt$cores))
  opt
}


# ---- inputs -----------------------------------------------------------------
# fread with ID columns kept as character; `#` is stripped from header names
fread_ids <- function(path, id_cols) {
  header <- sub("^#", "", names(fread(path, nrows = 0)))
  dt <- fread(path, colClasses = list(character = which(header %in% id_cols)))
  setnames(dt, header)
  dt
}

require_cols <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(what, " is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
}

load_inputs <- function(opt) {
  samples <- fread(opt$samples, header = FALSE, colClasses = "character")
  if (ncol(samples) < 2) stop("--samples needs two columns (FID IID)", call. = FALSE)
  samples <- unique(data.table(IID = samples[[2]]))

  manifest <- fread_ids(opt$manifest, c("GP2ID", "GP2sampleID"))
  require_cols(manifest, c("GP2ID", "GP2sampleID", "PD_PHENO", "qc_sex", "AGE_ANALYSIS", "wgs"), "--manifest")
  as_code <- function(v, col, allowed) {
    num <- suppressWarnings(as.numeric(as.character(v)))
    bad <- unique(v[!is.na(v) & (is.na(num) | !num %in% allowed)])
    if (length(bad))
      stop("manifest ", col, " must be ", paste(allowed, collapse = "/"), " or NA; found: ",
           paste(head(bad, 5), collapse = ", "), call. = FALSE)
    num
  }
  manifest <- manifest[, .(IID = GP2ID, GP2sampleID,
                           case_control = as_code(PD_PHENO, "PD_PHENO", c(0, 1)),
                           sex = as_code(qc_sex, "qc_sex", c(1, 2)),
                           AGE_ANALYSIS = suppressWarnings(as.numeric(as.character(AGE_ANALYSIS))),
                           wgs)]

  hla <- fread_ids(opt$hla_calls, "GP2sampleID")
  require_cols(hla, "GP2sampleID", "--hla-calls")
  allele_cols <- grep("^HLA-", names(hla), value = TRUE)
  if (!length(allele_cols)) stop("--hla-calls has no columns named HLA-*", call. = FALSE)
  hla <- hla[, c("GP2sampleID", allele_cols), with = FALSE]
  bad <- allele_cols[vapply(allele_cols, function(a) {
    v <- hla[[a]]; any(!is.na(v) & !v %in% 0:2)
  }, logical(1))]
  if (length(bad))
    stop("--hla-calls dosages must be 0, 1, 2 or NA; check column(s): ",
         paste(head(bad, 5), collapse = ", "), call. = FALSE)

  stats <- fread(opt$hla_stats)
  require_cols(stats, c("allele", "reads_mean"), "--hla-stats")

  pcs <- fread_ids(opt$pcs, "IID")
  pc_cols <- if (opt$n_pcs > 0) paste0("PC", seq_len(opt$n_pcs)) else character(0)
  require_cols(pcs, c("IID", pc_cols), "--pcs")
  pcs <- pcs[, c("IID", pc_cols), with = FALSE]

  related <- NULL
  if (!is.na(opt$related)) {
    related <- fread_ids(opt$related, c("IID1", "IID2"))
    require_cols(related, c("IID1", "IID2", "KINSHIP"), "--related")
  }

  list(samples = samples, manifest = manifest, hla = hla, allele_cols = allele_cols,
       stats = stats, pcs = pcs, pc_cols = pc_cols, related = related)
}


# ---- sample QC --------------------------------------------------------------
report <- function(label, dat) {
  cat(sprintf("  %-38s %6d  (cases %d, controls %d)\n", label, nrow(dat),
              sum(dat$case_control == 1, na.rm = TRUE), sum(dat$case_control == 0, na.rm = TRUE)))
}

prune_related <- function(dat, related, kinship) {
  if (is.null(related)) {
    cat("  no --related file: relatedness filter skipped\n")
    return(dat)
  }
  pairs <- related[KINSHIP >= kinship & IID1 %in% dat$IID & IID2 %in% dat$IID, .(IID1, IID2)]
  cases <- dat[case_control == 1, IID]
  exclude <- character(0)
  while (nrow(pairs) > 0) {
    ids <- unique(c(pairs$IID1, pairs$IID2))
    counts <- table(c(pairs$IID1, pairs$IID2))
    candidates <- ids[!ids %in% cases]
    if (!length(candidates)) candidates <- ids           # only case-case pairs left
    to_remove <- names(which.max(counts[names(counts) %in% candidates]))
    exclude <- c(exclude, to_remove)
    pairs <- pairs[IID1 != to_remove & IID2 != to_remove]
  }
  cat(sprintf("  related pairs (KINSHIP >= %g) resolved by removing %d individual(s)\n",
              kinship, length(exclude)))
  dat[!IID %in% exclude]
}

build_cohort <- function(inp, opt, with_age) {
  dat <- merge(inp$samples, inp$manifest, by = "IID");  report("in --samples and manifest", dat)
  dat <- merge(dat, inp$hla, by = "GP2sampleID");       report("with HLA calls", dat)
  dat <- merge(dat, inp$pcs, by = "IID");               report("with PCs", dat)
  dup <- dat[duplicated(IID), unique(IID)]
  if (length(dup))
    stop(length(dup), " individual(s) appear more than once after merging manifest, HLA calls ",
         "and PCs (e.g. ", dup[1], "); de-duplicate the inputs first", call. = FALSE)
  dat <- dat[wgs == 1];                                 report("passed WGS QC", dat)
  dat <- dat[!is.na(case_control)];                     report("case or control", dat)
  bad <- setdiff(unique(dat$case_control), c(0, 1))
  if (length(bad)) stop("PD_PHENO must be 1 (case), 0 (control) or NA; found ",
                        paste(bad, collapse = ", "), call. = FALSE)
  dat <- dat[!is.na(sex)];                              report("sex known", dat)
  if (with_age) {
    dat <- dat[!is.na(AGE_ANALYSIS)];                   report("age known", dat)
    dat <- dat[case_control == 1 | AGE_ANALYSIS >= opt$min_control_age]
    report(sprintf("controls aged >= %g", opt$min_control_age), dat)
  }
  if (nrow(dat) && anyNA(dat[, inp$pc_cols, with = FALSE])) {
    dat <- dat[complete.cases(dat[, inp$pc_cols, with = FALSE])]; report("PCs complete", dat)
  }
  dat <- prune_related(dat, inp$related, opt$kinship);  report("after relatedness filter", dat)
  cat(sprintf("  male %d, female %d\n", sum(dat$sex == 1), sum(dat$sex == 2)))
  if (sum(dat$case_control == 1) == 0 || sum(dat$case_control == 0) == 0)
    stop("no cases or no controls left after sample QC", call. = FALSE)
  dat
}


# ---- allele QC --------------------------------------------------------------
select_alleles <- function(dat, inp, opt) {
  good_reads <- paste0("HLA-", inp$stats[reads_mean >= opt$min_reads_mean, allele])
  alleles <- intersect(inp$allele_cols, good_reads)
  cat(sprintf("  alleles with reads_mean >= %g: %d\n", opt$min_reads_mean, length(alleles)))
  if (!length(alleles)) return(character(0))

  geno <- as.matrix(dat[, alleles, with = FALSE])
  cf <- colMeans(geno > 0, na.rm = TRUE)
  alleles <- alleles[which(cf > opt$min_carrier_freq)]
  cat(sprintf("  alleles with carrier frequency > %g: %d\n", opt$min_carrier_freq, length(alleles)))

  is_case <- dat$case_control == 1
  in_both <- vapply(alleles, function(a) {
    d <- dat[[a]]
    sum(d[is_case], na.rm = TRUE) >= 1 && sum(d[!is_case], na.rm = TRUE) >= 1
  }, logical(1))
  alleles <- alleles[in_both]
  cat(sprintf("  alleles carried in both cases and controls (tested): %d\n", length(alleles)))
  alleles
}


# ---- models -----------------------------------------------------------------
covariate_matrix <- function(dat, covars) {
  mm <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))),
                     model.frame(~ ., dat[, covars, with = FALSE], na.action = na.fail))
  colnames(mm) <- make.names(colnames(mm))
  mm
}

genotype_counts <- function(d) c(n_0 = sum(d == 0), n_1 = sum(d == 1), n_2 = sum(d == 2))

fit_term <- function(mm, y, family) {
  fit <- try(fastglm(x = mm, y = y, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  coefs <- summary(fit)$coefficients
  if (!"allele_term" %in% rownames(coefs)) return(NULL)
  coefs["allele_term", ]
}

join_notes <- function(...) {
  x <- c(...)
  x <- x[!is.na(x)]
  if (length(x)) paste(x, collapse = "; ") else NA_character_
}

# One case-control model for allele `a` under coding `model_type`
case_control_model <- function(a, model_type, dat, covar_mm, covars_label) {
  dosage <- dat[[a]]
  ok <- !is.na(dosage)
  d <- dosage[ok]
  y <- dat$case_control[ok]
  x <- MODEL_CODINGS[[model_type]][d + 1]
  cc <- genotype_counts(d[y == 1]); ct <- genotype_counts(d[y == 0])
  n_cases <- sum(y == 1); n_controls <- sum(y == 0)

  row <- data.table(
    allele = a, model_type = model_type,
    coding = paste(MODEL_CODINGS[[model_type]], collapse = "/"),
    beta = NA_real_, se = NA_real_, p = NA_real_, OR = NA_real_,
    n_cases = n_cases, n_controls = n_controls,
    case_n_0 = cc[["n_0"]], case_n_1 = cc[["n_1"]], case_n_2 = cc[["n_2"]],
    control_n_0 = ct[["n_0"]], control_n_1 = ct[["n_1"]], control_n_2 = ct[["n_2"]],
    case_af = (cc[["n_1"]] + 2 * cc[["n_2"]]) / (2 * n_cases),
    control_af = (ct[["n_1"]] + 2 * ct[["n_2"]]) / (2 * n_controls),
    direction = NA_character_, covars = covars_label,
    keep = as.integer(any(startsWith(a, CLASSICAL_GENES))),
    note = if (any(!ok)) paste0("n_missing_calls=", sum(!ok)) else NA_character_
  )

  if (length(unique(x)) < 2) {
    row$note <- join_notes("no_variation", row$note)
    return(row)
  }
  if (model_type != "additive" && any(table(factor(x, 0:1), factor(y, 0:1)) == 0))
    row$note <- join_notes("possible_separation", row$note)

  est <- fit_term(cbind(covar_mm[ok, , drop = FALSE], allele_term = x), y, binomial())
  if (is.null(est)) {
    row$note <- join_notes("model_failed", row$note)
    return(row)
  }
  row[, `:=`(beta = est[[1]], se = est[[2]], p = est[[4]], OR = exp(est[[1]]),
             direction = if (is.na(est[[1]])) NA_character_ else if (est[[1]] > 0) "Risk" else "Protective")]
  row
}

# Age-at-onset model (cases only, additive coding) for allele `a`
age_at_onset_model <- function(a, cases, covar_mm, covars_label) {
  dosage <- cases[[a]]
  ok <- !is.na(dosage)
  d <- dosage[ok]
  cc <- genotype_counts(d)
  n_cases <- length(d)

  row <- data.table(
    allele = a, model_type = "age_of_onset", coding = "0/1/2",
    beta = NA_real_, se = NA_real_, p = NA_real_,
    n_cases = n_cases,
    case_n_0 = cc[["n_0"]], case_n_1 = cc[["n_1"]], case_n_2 = cc[["n_2"]],
    n_carriers = cc[["n_1"]] + cc[["n_2"]],
    carrier_freq = (cc[["n_1"]] + cc[["n_2"]]) / n_cases,
    case_af = (cc[["n_1"]] + 2 * cc[["n_2"]]) / (2 * n_cases),
    direction = NA_character_, covars = covars_label,
    keep = as.integer(any(startsWith(a, CLASSICAL_GENES))),
    note = if (any(!ok)) paste0("n_missing_calls=", sum(!ok)) else NA_character_
  )

  if (length(unique(d)) < 2) {
    row$note <- join_notes("no_variation", row$note)
    return(row)
  }
  est <- fit_term(cbind(covar_mm[ok, , drop = FALSE], allele_term = d), cases$AGE_ANALYSIS[ok], gaussian())
  if (is.null(est)) {
    row$note <- join_notes("model_failed", row$note)
    return(row)
  }
  row[, `:=`(beta = est[[1]], se = est[[2]], p = est[[4]],
             direction = if (is.na(est[[1]])) NA_character_ else if (est[[1]] > 0) "Later onset" else "Earlier onset")]
  row
}

parallel_rows <- function(X, FUN, cores) {
  out <- mclapply(X, FUN, mc.cores = cores)
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("worker error: ", as.character(out[[which(failed)[1]]]), call. = FALSE)
  rbindlist(out, fill = TRUE)
}

add_fdr <- function(res, by = NULL) {
  res[, fdr_all := p.adjust(p, method = "fdr"), by = by]
  res[, fdr_keep := NA_real_]
  res[keep == 1, fdr_keep := p.adjust(p, method = "fdr"), by = by]
  res
}

run_case_control <- function(dat, alleles, covars, cores) {
  covar_mm <- covariate_matrix(dat, covars)
  label <- paste(covars, collapse = " + ")
  res <- rbindlist(lapply(names(MODEL_CODINGS), function(mt)
    parallel_rows(alleles, function(a) case_control_model(a, mt, dat, covar_mm, label), cores)))
  add_fdr(res, by = "model_type")
}

run_age_at_onset <- function(dat, alleles, covars, cores) {
  cases <- dat[case_control == 1]
  covar_mm <- covariate_matrix(cases, covars)
  label <- paste(covars, collapse = " + ")
  add_fdr(parallel_rows(alleles, function(a) age_at_onset_model(a, cases, covar_mm, label), cores))
}

write_results <- function(res, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(res, path)
  n_bad <- sum(is.na(res$p))
  cat(sprintf("  wrote %s (%d rows%s)\n", path, nrow(res),
              if (n_bad) sprintf(", %d without a p-value; see the note column", n_bad) else ""))
}


# ---- main -------------------------------------------------------------------
main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  inp <- load_inputs(opt)
  cat(sprintf("n_pcs = %d, cores = %d\n", opt$n_pcs, opt$cores))

  cat("\n== With age ==\n")
  dat <- build_cohort(inp, opt, with_age = TRUE)
  alleles <- select_alleles(dat, inp, opt)
  if (length(alleles)) {
    write_results(run_case_control(dat, alleles, c("AGE_ANALYSIS", "sex", inp$pc_cols), opt$cores),
                  file.path(opt$out_dir, "with_age", "hla_results_with_age.csv"))
    write_results(run_age_at_onset(dat, alleles, c("sex", inp$pc_cols), opt$cores),
                  file.path(opt$out_dir, "with_age", "aoo_results.csv"))
  }

  cat("\n== Without age ==\n")
  dat <- build_cohort(inp, opt, with_age = FALSE)
  alleles <- select_alleles(dat, inp, opt)
  if (length(alleles))
    write_results(run_case_control(dat, alleles, c("sex", inp$pc_cols), opt$cores),
                  file.path(opt$out_dir, "without_age", "hla_results_without_age.csv"))
}

main()
