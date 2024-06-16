## Chicken colour discrimination depends on background colour

This repository contains files relevant to 'Chicken colour discrimination depends on background colour', Journal of Experimental Biology DOI: 10.1242/jeb.209429 (https://doi.org/10.1242/jeb.209429). This article is written by Peter Olsson, Robin D. Johnsson, James J. Foster, John D. Kirwan, Olle Lind, and Almut Kelber.

The source code is primarily R markdown files run using the R language (with the package 'knitr'). HTML outputs are provided for these files in results.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11836380.svg)](https://doi.org/10.5281/zenodo.11836380)


### File manifest

The **src** directory contains the R code used for the statistical analysis of the behavioural data, within R markdown (.Rmd) files. Code for the Psychometric MLE experiments (method i) is within the files 'Psychometric_MLE_green.Rmd' and 'Psychometric_MLE_orange.Rmd'. Code for the Psychometric Bayesian model (method ii) is contained in the file 'Psychometric_Bayesian.Rmd'. Code for the Non-psychometric MLE model (method iii) is contained in the file 'Non-psychometric_MLE.Rmd'. The Stan model code for the main model in method ii is contained in 'Psychometric_Bayesian_model.Stan'.

The response data derived from the behavioural experiments are also included in the src directory. Data for the Psychometric MLE experiments (method i) is within the files 'green colours.csv' and 'orange colours.csv'. Data for the Psychometric Bayesian model (method ii) is contained in the file 'colour_discriminate_long_format.txt' and data for the Non-psychometric MLE model (method iii) is contained in the file 'colour_discriminate_short_format.txt'.

The **results** directory contains HTML outputs of the R markdown files in the src dir, with corresponding names. It also contains the model object of the Psychometric Bayesian main model as 'Psychometric_Bayesian_model_object.rds'. The document 'Colour discrimination on different colour backgrounds in chicken.pdf', within results, contains supplementary figures.

The **data** directory contains only the file 'Olsson_et_al_stimuli_ and_background_spectra.xlsx': Spreadsheets of light measurements used in these experiments.
