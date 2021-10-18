# Causal Modelling of Heavy-Tailed Variables and Confounders with Application to River Flow

Repository for the homonym paper, by O. C. Pasche, V. Chavez-Demoulin and A. C. Davison ([arXiv:2110.06686](https://arxiv.org/abs/2110.06686)).


## Abstract

Confounding variables are a recurrent challenge for causal discovery and inference. In many situations, complex causal mechanisms only manifest themselves in extreme events, or take simpler forms in the extremes.  Stimulated by data on extreme river flows and precipitation, we introduce a new causal discovery methodology for heavy-tailed variables that allows the use of a known potential confounder as a covariate and allows its effect to be almost entirely removed when the variables have comparable tails, and also decreases it sufficiently to enable correct causal inference when the confounder has a heavier tail.
We introduce a new parametric estimator for the existing causal tail coefficient and a permutation test. Simulations show that the methods work well and the ideas are applied to the motivating dataset.


## Repository contents

The `R_functions/` folder contains the necessary R functions to compute the causal tail coefficient parametric estimators, partial confidence intervals, modified peaks over threshold GPD fitting allowing box or linear constraints, and other helpers.

The `Simulations/` folder contains all the code and results related to the simulations for the non-parametric and parametric estimators.
- `CTCs_simulations.R` is the script used for the simulations of all causal tail coefficient estimators on the the three distributions and different tails (to change on lines 29-31). The outputs are avaliable in .RData format in the subfolder `CTCs_Simulations/`.
- `plot_CTCs_simulations_onelab_ALL.R` is used to plot the results for these simulations, for one or multiple distribution choices. The figures are avaliable in the `Results/CTCs_Simulations/` subfolder.
- `CTCs_Permutation_test_power.R` is the script used for the simulations of the permutation test for causality. The outputs are avaliable in .RData format in the subfolder `CTCs_Simulations/Permutation_tests/`.
- `plot_test_power_merged.R` is used to plot the results for the permutation test power simulations. The figures are avaliable in the `Results/Permutation_tests/` subfolder.
- `bootstrap_CTC_boxplots.R` is used to create the boxplots illustrating the bootstrap distribution of the difference statistic and its limitations.
- `old/LGPD_scale_constraints.R` (deprecated) contains tests for the behaviour of the GPD fit with box and linear constraints. Some outputs are avaliable in the `Results/constrained_lgpdctc/` subfolder.
- `old/LGPD_CTC_ParametersAnalysis.R` (deprecated) outputs the maximum liklihood estimated parameters for all versions of the LGPD causal tail coefficient estimators (post-fit correction, linear constraints and box-constraints in the student case). Some outputs are avaliable in the `Results/lgpdctc_ParametersAnalysis/` subfolder.


The `Switzerland/` folder contains the code and results related to the exploratory and causal analysis of the river discharge data.
- `Switzerland_exploration.ipynb` contain the explorataory data analysis for both discharge and precipitation data (basic statistics, dates range and missing values, creation of the interactive maps, time series plots, ...)
- `Stations_webscraping.ipynb` contains the code used to web-scrape additional information about the sations and to convert coordinates, from the official Swiss websites and online services (links in the code).
- `Swiss_SumerDischarges_CTC_causal.R` and `Swiss_SumerDischarges_CTC_indep.R` contain the functional code for the causal analysis of causal and non-causal station pairs, respectively.
- `Stations_runner.R` runs the causal analysis functions for a list of given station pairs.
- `Swiss_SumerDischarges_estims_plots.R` creates some summary plots and shape parameter estimates analyses. 
- The `Results/` subfolders contain some outputs from the python and R codes, such as the interactive maps in .html format (need to be downloaded and opened in a web-browser), and results discussed in the paper.


By Olivier Colin PASCHE\
EPFL Lausanne (CH), Spring 2020,\
Research Center for Statistics, University of Geneva (CH), 2021.\
Supported by the Swiss National Science Foundation.


