# Variant_specific_inflation
 Computing and studying variant-specific inflation expected due to stratified variance
 
 Please see variant_infl_tutorial.md for tutorial (or you can also download the pdf file). 
 
Useful functions are provided in the file variant_specific_inflation_functions.R.


# Content of the repository
- A tutorial about computing variant-specific inflation factors, includes multiple types of files called variant_infl_tutorial.html, .md, .pdf, .Rmd. The code used to generate the tutorial is in the .Rmd file. 
- Functions to compute variant-specific inflation factors (a single function), and to make a multi-panel QQ-plot, each panel corresponding to a different set of genetic variants. These functions are found in the file variant_specific_inflation_functions.R
- A script that was used to generate simulated data: Generate_simulated_data.R
- Simulated genotype and phenotype data: Genotypes.gds, Phenotypes.RData.


# Required software
This repository provides an R function for computing variant-specific inflation factors, and an example workflow on a mock dataset, in the tutorial. The example workflow uses the R/Bionconductor package GENESIS (used for fitting "null models" and for genetic association analyses), and the example genetic data is saved in a gds file. All R package versions are provided in the tutorial file. 

## Workflows and pipelines can be developed in many other ways
We note that workflows can be developed in which allele frequencies are computed using other softwares, other file formats for genetic data, and other ways to fit null models and perform genetic association testing. The function computing variant-specific inflation factor is implemented in R, so if other software is used compute allele frequencies, etc, the data should be read by R. 

# Installation guide
Please see the the tutorial for installation and usage. Currently, we do not provide a software but rather a source code in R. The function for computing variant-specific inflation factors for a given variant does not use any package (but packages are required to prepare input data to the function). The plotting function uses the packages ggplot2, data.table, and tidyr. 

One needs to follow standard R commands to install packages used, e.g. install.packages(). All packages used, by the plotting function and by the example workflow, are professionally maintained. 

# Running the tutorial
To run the tutorial one needs to install and load the R package rmarkdown, and use the command render("variant_infl_turorial.Rmd"). It takes less than a minute to run it. It runs over the example data and generate figures. Users can also run the R function from the tutorial line-by-line to look at intermediate outputs and to fit adapt the workflow to their needs. 






