# enterics-seroepi
Enteropathogen antibody dynamics and force of infection among children in low-resource settings


## Description

This repository includes data and computational notebooks to support the paper:

Arnold BF, Martin DL, Juma J, Mkocha H, Ochieng JB, et al. Enteropathogen antibody dynamics and force of infection among children in low-resource settings. _eLife_  2019; 8:e45594 doi: 10.7554/eLife.45594. https://elifesciences.org/articles/45594

The repository includes all files needed to replicate the analyses, and the repository is cross-referenced with the Open Science Framework (https://osf.io/r4av7).  For questions, please write Ben Arnold at UCSF (bfarnold@gmail.com). 

## Details

All notebooks are written in R markdown, and are located in the `R` subdirectory (.Rmd with corresponding, compiled .html file).  You should be able to replicate all analyses by cloning this directory and creating two new subdirectories called `figs` and `output` to store results (those files excluded to save space).  Notebooks reference relevant analyses from the paper in their names and titles.  

The R script `enterics-seroepi-run-all.R` runs the notebooks in order and will re-create the entire analysis. Before running the notbooks, you will need to create directories entitled `output` and `figs` one level up, alongside the `data` directory, to store saved results (binary files other than datasets not pushed to GitHub).

De-identifed data from the three countries is provided in the `data` subdirectory. Each dataset is provided in .csv and .rds format. The files with `2.rds` suffixes are datasets that include a few derived variables (e.g., seropositivity cutoffs and indicators) created by the analysis notebooks.

