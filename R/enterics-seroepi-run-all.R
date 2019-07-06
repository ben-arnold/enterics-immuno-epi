
#-----------------------------------
# This script runs all computational
# notebooks used in the article:
#
# Enteropathogen antibody dynamics
# and force of infection among children
# in low-resource settings
#-----------------------------------
library(here)
here::here()

#-----------------------------------
# Raw study data processing
# to create public data files
#
# (not run for public replication)
#-----------------------------------
rmarkdown::render(here::here("R/dm/haiti-enteric-ab-data-format.Rmd"),
                  output_file = here::here("R/dm/haiti-enteric-ab-data-format.html"))
rmarkdown::render(here::here("R/dm/asembo-enteric-ab-data-format.Rmd"),
                  output_file = here("R/dm/asembo-enteric-ab-data-format.html"))
rmarkdown::render(here::here("R/dm/kongwa-enteric-ab-data-format.Rmd"),
                  output_file = here::here("R/dm/kongwa-enteric-ab-data-format.html"))

#-----------------------------------
# Figure 1 - Antibody distributions
# (plus seropositivity indicators 
#  and creation of analysis2 files)
#-----------------------------------
rmarkdown::render(here::here("R/Fig1-haiti-ab-distributions.Rmd"),
                  output_file = here::here("R/Fig1-haiti-ab-distributions.html"))

rmarkdown::render(here::here("R/Fig1sup1sup2-kongwa-ab-distributions.Rmd"),
                  output_file = here::here("R/Fig1sup1sup2-kongwa-ab-distributions.html"))

rmarkdown::render(here::here("R/Fig1sup3-asembo-ab-distributions.Rmd"),
                  output_file = here::here("R/Fig1sup3-asembo-ab-distributions.html"))

#-----------------------------------
# Table 1 - summary of datasets
#-----------------------------------
rmarkdown::render(here::here("R/Table1-enteric-ab-list.Rmd"),
                             output_file = here::here("R/Table1-enteric-ab-list.html"))


#-----------------------------------
# Figure 2 - joint antibody distributions
# + SI File 3
#-----------------------------------
rmarkdown::render(here::here("R/SI-File3-Fig2-ab-joint-distributions.Rmd"),
                  output_file = here::here("R/SI-File3-Fig2-ab-joint-distributions.html"))

#-----------------------------------
# Figure 3 - age curves
# + SI File 6
#-----------------------------------
rmarkdown::render(here::here("R/SI-File6-Fig3-agecurve-analyses.Rmd"),
                  output_file = here::here("R/SI-File6-Fig3-agecurve-analyses.html"))

#-----------------------------------
# Figure 4 - antibody trajectories
#-----------------------------------
rmarkdown::render(here::here("R/Fig4-Fig4sup2-haiti-ab-profiles.Rmd"),
                  output_file = here::here("R/Fig4-Fig4sup2-haiti-ab-profiles.html"))

rmarkdown::render(here::here("R/Fig4sup1-Fig4sup3-asembo-ab-profiles.Rmd"),
                  output_file = here::here("R/Fig4sup1-Fig4sup3-asembo-ab-profiles.html"))


#-----------------------------------
# Table 2 - Haiti incidence analysis
#-----------------------------------
rmarkdown::render(here::here("R/Table2-haiti-incidence-analysis.Rmd"),
                  output_file = here::here("R/Table2-haiti-incidence-analysis.html"))


#-----------------------------------
# Figure 6 - Kenya incidence analysis,
# + SI File 7 cross-sectional estimators
#-----------------------------------
rmarkdown::render(here::here("R/Fig6-asembo-incidence-analysis.Rmd"),
                  output_file = here::here("R/Fig6-asembo-incidence-analysis.html"))

rmarkdown::render(here::here("R/SI-File7-Fig6-asembo-cross-sectional-FOI-estimation.Rmd"),
                  output_file = here::here("R/SI-File7-Fig6-asembo-cross-sectional-FOI-estimation.html"))


#-----------------------------------
# Figure 5 - seroprevalence vs FOI
#-----------------------------------
rmarkdown::render(here::here("R/Fig5-foi-v-prev-asembo-haiti.Rmd"),
                  output_file = here::here("R/Fig5-foi-v-prev-asembo-haiti.html"))


#-----------------------------------
# Figure 7 - sampling interval simulation
# + SI File 8
#-----------------------------------
rmarkdown::render(here::here("R/SI-File8-sampling-interval-simulation.Rmd"),
                  output_file = here::here("R/SI-File8-sampling-interval-simulation.html"))

rmarkdown::render(here::here("R/Fig7-sampling-interval-adjusted-foi.Rmd"),
                  output_file = here::here("R/Fig7-sampling-interval-adjusted-foi.html"))

#-----------------------------------
# ADDITIONAL SI FILES BELOW THAT
# DO NOT GENERATE PRIMARY MANUSCRIPT
# FIGURES OR TABLES
#-----------------------------------

#-----------------------------------
# SI File 1
# intervention, bead, and season effects
#-----------------------------------
rmarkdown::render(here::here("R/SI-File1-intervention-bead-season-effects.Rmd"),
                  output_file = here::here("R/SI-File1-intervention-bead-season-effects.html"))

#-----------------------------------
# SI File 2
# seropositivity cutoff agreement
#-----------------------------------
rmarkdown::render(here::here("R/SI-File2-cutoff-agreement.Rmd"),
                  output_file = here::here("R/SI-File2-cutoff-agreement.html"))

#-----------------------------------
# SI File 4 
# serology vs stool infections in Kenya
#-----------------------------------
rmarkdown::render(here::here("R/SI-File4-serology-by-infection-status.Rmd"),
                  output_file = here::here("R/SI-File4-serology-by-infection-status.html"))

#-----------------------------------
# SI File 5
# fold-change in IgG sensitivity analyses
#-----------------------------------
rmarkdown::render(here::here("R/SI-File5-IgG-fold-change-sensitivity-analyses.Rmd"),
                  output_file = here::here("R/SI-File5-IgG-fold-change-sensitivity-analyses.html"))


