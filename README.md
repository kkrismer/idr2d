# IDR2D: Irreproducible Discovery Rate for Genomic Interactions

[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fnar%2Fgkaa030-blue.svg)](https://doi.org/10.1093/nar/gkaa030) [![BioC](https://img.shields.io/badge/BioC-1.16.0-brightgreen.svg)](https://doi.org/doi:10.18129/B9.bioc.idr2d) [![platforms](https://www.bioconductor.org/shields/availability/release/idr2d.svg)](https://www.bioconductor.org/packages/release/bioc/html/idr2d.html#archives) [![Coverage Status](https://coveralls.io/repos/github/kkrismer/idr2d/badge.svg?branch=master)](https://coveralls.io/github/kkrismer/idr2d?branch=master)

https://idr2d.mit.edu

Chromatin interaction data from protocols such as ChIA-PET and HiChIP provide valuable insights into genome organization and gene regulation, but can include spurious interactions that do not reflect underlying genome biology. We introduce a generalization of the Irreproducible Discovery Rate (IDR) method called IDR2D that identifies replicable interactions shared by experiments. IDR2D provides a principled set of interactions and eliminates artifacts from single experiments.

## Installation

The *idr2d* package is part of Bioconductor since release 3.10. To install it on your system, enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("idr2d")
```

Alternatively, the development version can be installed directly from this repository:

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("kkrismer/idr2d")
```

R 3.6 (or higher) and Bioconductor 3.10 (or higher) is required in both cases. Additionally, the 64-bit version of Python 3.5 (or higher) and the Python package [hic-straw](https://pypi.org/project/hic-straw/) are required for Hi-C analysis from Juicer *.hic* files. 

## Usage

There are two vignettes available on Bioconductor, focusing on [*idr2d* and ChIA-PET data](https://www.bioconductor.org/packages/release/bioc/vignettes/idr2d/inst/doc/idr2d.html) and [*idr2d* and ChIP-seq data](https://www.bioconductor.org/packages/release/bioc/vignettes/idr2d/inst/doc/idr1d.html).

The [reference manual](https://bioc.ism.ac.jp/packages/devel/bioc/manuals/idr2d/man/idr2d.pdf) might also be helpful if you know what you are looking for.

### Example code for ChiP-seq, ChIA-PET and Hi-C experiments

Analyzing results from replicate **ChIP-seq** experiments
(stored in tab-delimited files *chip-seq-rep1.txt* and *chip-seq-rep2.txt*):
```
library(idr2d)

rep1_df <- read.table("chip-seq-rep1.txt", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
rep2_df <- read.table("chip-seq-rep2.txt", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)

idr_results <- estimate_idr1d(rep1_df, rep2_df, 
                              value_transformation = "identity")
summary(idr_results)

rep1_idr_df <- idr_results$rep1_df
draw_idr_distribution_histogram(rep1_idr_df)
draw_rank_idr_scatterplot(rep1_idr_df)
draw_value_idr_scatterplot(rep1_idr_df)
```

Analyzing results from replicate **ChIA-PET** experiments
(stored in tab-delimited files *chia-pet-rep1.txt* and *chia-pet-rep2.txt*):
```
library(idr2d)

rep1_df <- read.table("chia-pet-rep1.txt", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
rep2_df <- read.table("chia-pet-rep2.txt", header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)

idr_results <- estimate_idr2d(rep1_df, rep2_df, 
                              value_transformation = "identity")
summary(idr_results)

rep1_idr_df <- idr_results$rep1_df
draw_idr_distribution_histogram(rep1_idr_df)
draw_rank_idr_scatterplot(rep1_idr_df)
draw_value_idr_scatterplot(rep1_idr_df)
```

Analyzing chromosome 1 results in 1 Mbp resolution from replicate **Hi-C** experiments
(stored in Juicer .hic files *hic-rep1.hic* and *hic-rep2.hic*):
```
library(idr2d)

rep1_df <- parse_juicer_matrix("hic-rep1.hic", resolution = 1e+06, chromosome = "chr1")
rep2_df <- parse_juicer_matrix("hic-rep2.hic", resolution = 1e+06, chromosome = "chr1")

idr_results_df <- estimate_idr2d_hic(rep1_df, rep2_df)
summary(idr_results_df)

draw_idr_distribution_histogram(idr_results_df)
draw_rank_idr_scatterplot(idr_results_df)
draw_value_idr_scatterplot(idr_results_df)
draw_hic_contact_map(idr_results_df, idr_cutoff = 0.05, chromosome = "chr1")
```

Analyzing chromosome 1 results in 1 Mbp resolution from replicate **Hi-C** experiments
(stored in ICE normalized HiC-Pro .matrix and .bed files *rep1_1000000_iced.matrix*, *rep1_1000000_abs.bed* and *rep2_1000000_iced.matrix*, *rep2_1000000_abs.bed*):
```
library(idr2d)

rep1_df <- parse_hic_pro_matrix("rep1_1000000_iced.matrix", "rep1_1000000_abs.bed", chromosome = "chr1")
rep2_df <- parse_hic_pro_matrix("rep2_1000000_iced.matrix", "rep2_1000000_abs.bed", chromosome = "chr1")

idr_results_df <- estimate_idr2d_hic(rep1_df, rep2_df)
summary(idr_results_df)

draw_idr_distribution_histogram(idr_results_df)
draw_rank_idr_scatterplot(idr_results_df)
draw_value_idr_scatterplot(idr_results_df)
draw_hic_contact_map(idr_results, idr_cutoff = 0.05, chromosome = "chr1")
```

## Build status

| Platform | Status |
|------|------|
| Travis CI | [![Travis build status](https://travis-ci.com/kkrismer/idr2d.svg?branch=master)](https://travis-ci.com/kkrismer/idr2d) |
| Bioconductor 3.18 (release) | [![BioC release](https://bioconductor.org/shields/build/release/bioc/idr2d.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/idr2d/) |
| Bioconductor 3.19 (devel) | [![BioC devel](https://bioconductor.org/shields/build/devel/bioc/idr2d.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/idr2d/) |

## Citation

If you use IDR2D in your research, please cite:

**IDR2D identifies reproducible genomic interactions**  
Konstantin Krismer, Yuchun Guo, and David K. Gifford  
Nucleic Acids Research, Volume 48, Issue 6, 06 April 2020, Page e31; DOI: https://doi.org/10.1093/nar/gkaa030

## Funding

The development of this method was supported by National Institutes of Health (NIH) grants 1R01HG008363 and 1R01NS078097, and the MIT Presidential Fellowship.
