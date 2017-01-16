# breedTools

> R package encapsulating methods presented in *Estimation of genome-wide and
> locus-specific breed composition in pigs* (doi: 10.2527/tas2016.0003)
> published in Translational Animal Science (TAS).

breedTools can be used to estimate breed composition from 
DNA microarray data and has been applied to the estimation of pig breed
composition and purity.

The methods in this package require a set of reference data (genotypes).

## Table of Contents

1. [Installation](#installation)
2. [Usage](#Usage)

## Installation

**NOTE**: Several functions in `breedTools` rely on another package, 
[`snpTools`](https://github.com/funkhou9/snpTools).
First ensure [`devtools`](https://github.com/hadley/devtools) is installed,
then install `snpTools` with:

```r
devtools::install_github("funkhou9/snpTools")
```

then install `breedTools` with:

```r
devtools::install_github("funkhou9/breedTools")
```

attaching the `breedTools` package with:

```r
library(breedTools)
```

will allow access to reference pig datasets.

## Usage

Key functions include:

* `solve_composition()` to estimate genome-wide breed composition (GWBC) for a 
set of animals, given a reference panel.

* `solve_KBP()` to estimate local, *KIT*-based breed probabilities (KBP) for a 
set of animals, based on looking at haplotypes within the *KIT* region on
chromosome 8.

* `screen_purity()` combines GWBC and KBP methods and provides a full report
for a set of animals.

* `build_KBP()` is used to build a new reference panel for KBP estimation.

* `allele_freq()` is used to build a new reference panel for GWBC estimation.

Note other functions in `snpTools` such as `snpTools::filter_geno()`, which
filters out SNPs based on MAF.

Example GWBC and KBP reference panels `GWBC_ref_A`, `GWBC_ref_B`, `KBP_ref_A`,
and `KBP_ref_B` can be loaded with, for example:

```r
data("GWBC_ref_A")
```

Note that `GWBC_ref_A` and `KBP_ref_A` are used to generate the results in the
TAS paper.