# enterics-seroepi
Seroepidemiology of enteropathogens among children in low-resource settings


## Description

This repository includes data and computational notebooks to support the paper entitled _Enteropathogen seroepidemiology among children in low-resource settings_ (in review). All files required to replicate the analyses are included herein, and the repository is cross-referenced with the Open Science Framework (https://osf.io/r4av7).  For questions, please write Ben Arnold at UC Berkeley (benarnold@berkeley.edu). 

## Details

All notebooks are written in R markdown, and are located in the `R` subdirectory (.Rmd with corresponding, compiled .html file).  You should be able to replicate all analyses by cloning this directory and creating two new subdirectories called `figs` and `output` to store results (those files excluded to save space).  Notebooks reference relevant analyses from the paper in their names and titles.  

De-identifed data from the three countries is provided in the `data` subdirectory. Each dataset is provided in .csv and .rds format. The files with `2.rds` suffixes are datasets that include a few derived variables (e.g., seropositivity cutoffs and indicators) created by the analysis notebooks.

