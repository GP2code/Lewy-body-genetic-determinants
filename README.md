# Genetic determinants of Lewy body pathology

`GP2 ❤️ Open Science 😍`

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23653763.svg)](https://doi.org/10.5281/zenodo.23653763)

**Last updated:** September 2026

## Summary

This is the online repository for the manuscript titled ***"Genetic determinants of Lewy body pathology"***. To investigate the genetic determinants of neuropathologically defined Lewy body disease (LBD), we:
- Performed a genome-wide association study (GWAS) of a large cohort of autopsy-confirmed LBD cases, independent of clinical phenotype.
- Compared genetic associations in neuropathologically defined LBD with those reported for clinically diagnosed Parkinson's disease (PD).
- Identified genetic determinants associated with cortical versus predominantly subcortical Lewy body pathology.

## Data statement
Data used in the preparation of this article were obtained from the Global Parkinson's Genetics Program (GP2; https://gp2.org), Mayo Brain Bank and NIH NeuroBioBank.

All GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson's disease, and are available via application on the website (https://amp-pd.org/register-for-amp-pd). For up-to-date information on GP2 data acquisition, access, and policies, visit https://gp2.org/. Tier 1 data can be accessed by completing a form on the Accelerating Medicines Partnership in Parkinson's Disease (AMP®-PD) website (https://amp-pd.org/register-for-amp-pd). Tier 2 data access requires approval and a Data Use Agreement signed by your institution.

In this analysis we used Tier 2 GP2 Release 11 data ([10.5281/zenodo.17753486](https://doi.org/10.5281/zenodo.17753486)) and GP2 Release 12 whole-genome sequencing data ([10.5281/zenodo.20932193](https://doi.org/10.5281/zenodo.20932193)).

### Helpful Links
- [GP2 Website](https://gp2.org/)
  - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
  - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

## Citation
*(pending publication)*

## Repository Orientation
- The `analyses/` directory includes all analyses discussed in the manuscript.

<pre>
THIS_REPO/
├── analyses/
│   ├── 01_GWAS/
│   │   ├── 00_logistic_GWAS_plink2.sh
│   │   ├── 01_ordinal_GWAS.jl
│   │   └── 02_manhattan_plot.R
│   ├── 02_Meta-analysis/
│   │   ├── 00_prepare_METAL_input.R
│   │   └── 01_METAL_meta-analysis.txt
│   ├── 03_Conditional_analysis.ipynb
│   ├── 04_Beta-beta.py
│   ├── 05_LDSC.py
│   └── 06_HLA_Association.R
└── README.md
</pre>

## Analysis Notebooks
### Languages: bash, Julia, R, Python
| Directory | Notebooks / scripts | Description |
|-----------|---------------------|-------------|
| `analyses/01_GWAS/` | `00_logistic_GWAS_plink2.sh` | Per-cohort logistic GWAS (PLINK 2.0) of autopsy-confirmed LBD vs controls, and cortical vs subcortical LBD |
| `analyses/01_GWAS/` | `01_ordinal_GWAS.jl` | Ordinal GWAS of brainstem-, limbic- and neocortical-predominant LBD (OrdinalGWAS.jl) |
| `analyses/01_GWAS/` | `02_manhattan_plot.R` | Manhattan plot of the summary statistics |
| `analyses/02_Meta-analysis/` | `00_prepare_METAL_input.R` | QC of per-cohort PLINK 2.0 results and METAL input files |
| `analyses/02_Meta-analysis/` | `01_METAL_meta-analysis.txt` | Standard-error weighted meta-analysis across cohorts (METAL) |
| `analyses/` | `03_Conditional_analysis.ipynb` | GCTA-COJO conditional and joint analyses at APOE and SNCA, with GP2 EUR controls as the LD reference |
| `analyses/` | `04_Beta-beta.py` | Effect-size comparison of clinical PD GWAS signals with the LBD GWAS |
| `analyses/` | `05_LDSC.py` | Genetic correlation of LBD with PD, Alzheimer's disease and clinical LBD (LDSC) |
| `analyses/` | `06_HLA_Association.R` | HLA allele association and age-at-onset analyses using GP2 WGS HLA calls |

## Software
| **Software** | **Version(s)** | **Resource URL** | **RRID** | **Notes** |
|--------------|----------------|------------------|----------|-----------|
| PLINK | v.2.0 | http://www.nitrc.org/projects/plink | RRID:SCR_001757 | Used for GWAS, genotype QC and LD calculation |
| OrdinalGWAS.jl | | https://github.com/OpenMendel/OrdinalGWAS.jl | NA | Used for ordinal GWAS |
| Julia Programming Language | | https://julialang.org/ | RRID:SCR_021666 | Used to run OrdinalGWAS.jl |
| METAL | | http://csg.sph.umich.edu/abecasis/Metal/ | RRID:SCR_002013 | Used for meta-analysis |
| GCTA | v.1.94.1 | https://yanglab.westlake.edu.cn/software/gcta/ | NA | Used for COJO conditional and joint analyses |
| LDSC | v.1.0.1 | https://github.com/bulik/ldsc | RRID:SCR_022801 | Used for genetic correlation |
| Python Programming Language | 3.10, 3.11 | http://www.python.org/ | RRID:SCR_008394 | pandas (RRID:SCR_018214); NumPy (RRID:SCR_008633); SciPy (RRID:SCR_008058); Matplotlib (RRID:SCR_008624); adjustText. Used for data wrangling, analyses and plotting |
| Jupyter Notebook | | https://jupyter.org/ | RRID:SCR_018315 | Used to run analysis notebooks |
| R Project for Statistical Computing | | http://www.r-project.org/ | RRID:SCR_001905 | data.table; dplyr (RRID:SCR_016708); ggplot2 (RRID:SCR_014601); ggrepel (RRID:SCR_017393); fastglm. Used for data wrangling, regression and plotting |
