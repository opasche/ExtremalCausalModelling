# Causal Modelling of Heavy-Tailed Variables and Confounders with Application to River Flow

Repository for the homonym paper, by O. C. Pasche, V. Chavez-Demoulin and A. C. Davison, published in *Extremes*.
* Publication available in open access at: [https://doi.org/10.1007/s10687-022-00456-4](https://doi.org/10.1007/s10687-022-00456-4).
* Arxiv: [arXiv:2110.06686](https://arxiv.org/abs/2110.06686).


## Abstract

Confounding variables are a recurrent challenge for causal discovery and inference. In many situations, complex causal mechanisms only manifest themselves in extreme events, or take simpler forms in the extremes.  Stimulated by data on extreme river flows and precipitation, we introduce a new causal discovery methodology for heavy-tailed variables that allows the effect of a known potential confounder to be almost entirely removed when the variables have comparable tails, and also decreases it sufficiently to enable correct causal inference when the confounder has a heavier tail.  We also introduce a new parametric estimator for the existing causal tail coefficient and a permutation test. Simulations show that the methods work well and the ideas are applied to the motivating dataset.


## Repository contents

The `./R_functions/` folder contains the necessary R functions to compute the causal tail coefficient parametric estimators, partial confidence intervals, modified peaks over threshold GPD fitting allowing box or linear constraints, and other helpers.

The `./Simulations/` folder contains all the code and results related to the simulations for the non-parametric and parametric estimators.
- `CTCs_simulations.R` is the script used for the simulations of all causal tail coefficient estimators on the the three distributions and different tails (to change on lines 29-31). The outputs are avaliable in .RData format in the subfolder `CTCs_Simulations/`.
- `plot_CTCs_simulations_onelab_ALL.R` is used to plot the results for these simulations, for one or multiple distribution choices. The figures are avaliable in the `Results/CTCs_Simulations/` subfolder.
- `CTCs_Permutation_test_power.R` is the script used for the simulations of the permutation test for causality. The outputs are avaliable in .RData format in the subfolder `CTCs_Simulations/Permutation_tests/`.
- `plot_test_power_merged.R` is used to plot the results for the permutation test power simulations. The figures are avaliable in the `Results/Permutation_tests/` subfolder.
- `bootstrap_CTC_boxplots.R` is used to create the boxplots illustrating the bootstrap distribution of the difference statistic and its limitations.
- `old/LGPD_scale_constraints.R` (deprecated) contains tests for the behaviour of the GPD fit with box and linear constraints. Some outputs are avaliable in the `Results/constrained_lgpdctc/` subfolder.
- `old/LGPD_CTC_ParametersAnalysis.R` (deprecated) outputs the maximum liklihood estimated parameters for all versions of the LGPD causal tail coefficient estimators (post-fit correction, linear constraints and box-constraints in the student case). Some outputs are avaliable in the `Results/lgpdctc_ParametersAnalysis/` subfolder.


The `./Switzerland/` folder contains the code and results related to the exploratory and causal analysis of the river discharge data.
- `Switzerland_exploration.ipynb` contain the explorataory data analysis for both discharge and precipitation data (basic statistics, dates range and missing values, creation of the interactive maps, time series plots, ...)
- `Stations_webscraping.ipynb` contains the code used to web-scrape additional information about the sations and to convert coordinates, from the official Swiss websites and online services (links in the code).
- `Swiss_SumerDischarges_CTCs_routine.R` contains the functional code for the causal analysis of station pairs, with (optionnal) esthetical wrappers for causal and non-causal station pair types, respectively.
- `Stations_runner.R` runs the causal analysis functions for a list of given station pairs.
- `Swiss_SumerDischarges_estims_plots.R` creates some summary plots and shape parameter estimates analyses. 
- The `Results/` subfolders contain some outputs from the python and R codes, such as the interactive maps in .html format (need to be downloaded and opened in a web-browser), and results discussed in the paper.
- The `Data/` subfolder should contain the Swiss river discharge dataset (`river_dat.Rdata`) and Swiss precipitation dataset (`precip.csv`). We unfortunately do not have rights to share them directly here, but they are available to academics upon request, from the [Swiss Federal Office for the Environment (FOEN)](https://www.hydrodaten.admin.ch/) and [MeteoSwiss](https://gate.meteoswiss.ch/idaweb), respectively. Some additional processed data is then output in the `data_wrangled/` subfolder by, in particular, `Stations_webscraping.ipynb` and `RData_to_csv.R`.


By Olivier C. PASCHE\
EPFL Lausanne (CH), Spring 2020,\
Research Center for Statistics, University of Geneva (CH), 2021.\
Supported by the Swiss National Science Foundation.


