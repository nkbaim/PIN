# PIN

### Installation of development version from GitHub

```sh
$ library(devtools)
$ install_github("shaowenguang/PIN")
```

### Installation of release version from CRAN

Needs to wait for the official release...

```sh
$ install.packages("PIN")
```

### Description

This document describes a computational method, namely PIN (proteome integrity number), to assess a quantitative measure of protein stability directly from bottom-up proteomic datasets.

By quantifying the relative abundance of semi-tryptic peptides (i.e. the likely products of protein degradation), an individual protein integrity score (iPIS) was first calculated for each measured protein. A per-sample PIN was then derived to summarize a proteome-wide measurement indicating the level of protein degradation.

This method has been used in a clinical cohort of prostate tissue samples, in which degradation is usually inevitable during the sample preparation.




### Usage



First: to load the package.
```sh
$ library(PIN)
```

Second: to extract the peptide list, together with their quantities, from OpenSWATH outputs.
```sh
$ peptides <- generate_peptide_table(search_results="./openswath_search_results.tsv", sample_annotation="./sample_annotation_table", sptxt="./spectral_library.sptxt")
```

Third: to perform PIN analysis.
```sh
$ perform_PIN_analysis(peptides)
```
#### Outputs
| File Name | Description |
| ------ | ------ |
| PIN.tsv | PINs and their p-values, for each sample |
| iPIS.tsv | individual Protein Integrity Score, for each individual proteins |
| lib_peptide.tsv | a temp file that records the assay library information |





### Todos

 -  
 -  

License
----
