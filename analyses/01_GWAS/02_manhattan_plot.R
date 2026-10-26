#!/usr/bin/env Rscript
# =============================================================================
# Manhattan plot of GWAS / meta-analysis summary statistics
# - Project: Genetic determinants of Lewy body pathology
# - Version: R; data.table, ggplot2, ggrepel
# - Last Updated: 16-SEP-2026
#
# Overview
# - Genome-wide Manhattan plot (autosomes) with genome-wide (5e-8) and suggestive
#   (1e-5) lines; the top variant within each labelled locus window is highlighted.
#
# Usage
#   Rscript 02_manhattan_plot.R <summary_stats.tab> [manhattan.png]
#   Input: tab-separated, header line, 18 columns in this order:
#     CHROM BP Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr P Direction
#     HetISq HetChiSq HetDf HetPVal MAF_variability SNP   (header names are ignored)
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 02_manhattan_plot.R <summary_stats.tab> [manhattan.png]", call. = FALSE)
in_file  <- args[1]
out_file <- if (length(args) >= 2) args[2] else "manhattan.png"
if (!file.exists(in_file)) stop("file not found: ", in_file, call. = FALSE)

col_names <- c("CHROM", "BP", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq",
               "MaxFreq", "Effect", "StdErr", "P", "Direction", "HetISq",
               "HetChiSq", "HetDf", "HetPVal", "MAF_variability", "SNP")

# Loci to label (hg38): the top variant inside each window is highlighted
LOCI <- data.table(
  gene  = c("GBA1",        "TMEM175",   "SNCA",       "HLA",        "APOE"),
  CHROM = c(1L,            4L,          4L,           6L,           19L),
  start = c(154.2e6,       0,           88.7e6,       28.5e6,       43.9e6),
  end   = c(156.2e6,       2.0e6,       90.8e6,       33.5e6,       45.9e6)
)

gwas <- fread(in_file, header = FALSE, skip = 1, col.names = col_names)

# Autosomes only
gwas[, CHROM := suppressWarnings(as.integer(sub("^chr", "", as.character(CHROM), ignore.case = TRUE)))]
gwas <- gwas[!is.na(P) & CHROM %in% 1:22]
if (!nrow(gwas)) stop("no autosomal variants with a P-value in ", in_file, call. = FALSE)
gwas[, POS := as.numeric(BP)]
gwas[, logP := -log10(P)]
setorder(gwas, CHROM, POS)

# Cumulative position
chr_sizes <- gwas[, .(chr_len = max(POS, na.rm = TRUE)), by = CHROM]
chr_sizes[, tot := cumsum(chr_len) - chr_len]
gwas <- merge(gwas, chr_sizes, by = "CHROM")
gwas[, BPcum := POS + tot]

# Axis centres
axisdf <- gwas[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHROM]

# Top variant per labelled locus
top_hits_labeled <- rbindlist(lapply(seq_len(nrow(LOCI)), function(i) {
  l <- LOCI[i]
  hit <- gwas[CHROM == l$CHROM & POS >= l$start & POS <= l$end][which.min(P)]
  if (nrow(hit)) hit[, gene := l$gene] else NULL
}))

# Blue-grey palette for alternating chromosomes
chr_colors <- c("0" = "#8DA0CB", "1" = "#4D4D4D")

p <- ggplot(gwas, aes(x = BPcum, y = logP)) +
  geom_point(aes(color = as.factor(CHROM %% 2)), size = 0.6, alpha = 0.8) +
  scale_color_manual(values = chr_colors) +
  geom_point(data = top_hits_labeled, aes(x = BPcum, y = logP), color = "#D62728", size = 2) +
  geom_text_repel(
    data = top_hits_labeled, aes(x = BPcum, y = logP, label = gene),
    color = "#D62728", fontface = "italic", size = 3.5,
    nudge_y = 1.5, segment.color = "grey40", segment.size = 0.3
  ) +
  scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "grey60", linewidth = 0.3) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(color = "black"),
    text = element_text(family = "sans")
  )

ggsave(out_file, plot = p, width = 10, height = 5, dpi = 300)
cat("Saved", out_file, "\n")
