# Metabolic signatures of youth exposure to mixtures of per- and polyfluoroalkyl substances

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Demo](#demo)


# Overview  

Human exposure to perfluoroalkyl substances (PFAS) is ubiquitous and has been associated with several cardiometabolic diseases. However, the metabolic pathways linking PFAS exposure and human disease are unclear. Here, we examined associations of PFAS exposure and metabolic pathways in two independent cohorts of children and young adults. We analyzed untargeted metabolomics using Bayesian hierarchical regression modeling approach (BHRMA) with a g-prior and g-computation for modeling exposure-mixtures to estimate the impact of exposure to a mixture of six ubiquitous PFAS. This repository contains all code used in the analysis of this project. It also contains a reproducible example of the BHRMA, described in detail in the manuscript.


# Repo contents  

- [0_project_setup](./0_project_setup): Code for setting up the project, including loading packages and data. 
- [1_descriptive_statistics](./1_descriptive_statistics): Code for generating descriptive statistics.  
- [2_mixtures_MWAS](./2_mixtures_MWAS): Code for performing a metabolome wide association study (MWAS) in both the SOLAR and CHS cohort, using a Bayesian regression with g-estimation and g-computation for modeling exposure-mixtures to estimate the impact of exposure to a mixture of six PFAS on each individual metabolite feature.  
- [3_pathway_analysis](./3_pathway_analysis): Code for performing a pathway analysis and analyzing the results from the functional analysis module from MetaboAnalyst 5.0.
- [4_results_figures](./4_results_figures): Code for creating the final figures for the manuscript.  
- [5_reproducible_example](./5_reproducible_example): Code for running an example Bayesian regression with g-estimation and g-computation model using simulated data. 


# System requirements

## Hardware Requirements

The reproducible example of the BHRMA will take approximately 4 minutes to run using a standard computer with at least 8 GB of ram and a CPU with 4+ cores and 3.3+ GHz/core. 

The full mixtures MWAS analysis, which included running the BHRMA individually across 23,166 metabolites in two individual cohorts, was performed on the University of Southern California Center for Advanced Research Computing's High Performance Computing Cluster, using 128 cores running in parallel (https://www.carc.usc.edu/). 


## Software

### OS 

This code was written with a computer running a Windows 11 operating system. This code was tested with a computer running macOS 11. 

### Programs

In order to run the BHRMA, users must download JAGS (just another gibbs sampler), which is available from sourceforge (https://sourceforge.net/projects/mcmc-jags/). Typical install time for JAGS on a normal desktop computer is less than 5 minutes. In addition, users must install r packages `rjags`,  `R2jags`, and `tidyverse` (details below).


### R Package dependencies

Users should install the following packages from CRAN prior to running the BHRMA. From an `R` terminal: 

```
install.packages(c('tidyverse', 'rjags', 'R2jags'))
```


Additional libraries used throughout this project can be installed from an `R` terminal using the following code:  

```
# Packages for analysis, parralelization, and other tasks.
install.packages(c('purrr', 'broom', 'pbdMPI', 'here', 'fs',
                   'devtools', 'janitor', 'tidylog'))
                   
# Packages for data visualization                   
install.packages(c('kableExtra', 'gplots', 'cowplot', 'correlation', 'GGally', 'RColorBrewer', 'dendextend','ggrepel','colorspace','ggExtra'))

# Github package containing functions for summarizing and formatting numeric and categorical variables
devtools::install_github("JAGoodrich/jag2")
```

Typical install time for all `R` package dependencies on a normal desktop computer is less than 10 minutes. 



# Demo  

Detailed instructions for running the BHRMA model using simulated data is provided in (./5_reproducible_example/BHRMA_g_analysis_example.R). On a normal desktop computer, the expected run time for a single model using the simulated data is approximately 3-5 minutes.  

The output of this function will be a data frame named "fit" with 6 columns and 15 rows. Column headers are: 
```
c("var.names", "Mean", "SD", "X2.5.", "X97.5.", "p.val")  
```   

Where the mean, SD, X2.5, and X97.5 are the mean, SD, lower, and upper 95% confidence intervals for the posterior distribution of the specified model parameter. 

Model parameters are provided in rows. Specificially, each PFAS in the mixture has two rows in the data frame (provided in the var.names column): a row ending in ".beta", which is the effect estimate for each individual PFAS towards the overall mixture effect, and a row ending in ".gamma", which is the posterior inclusion probability (PIP) for the individual PFAS. Additionally, the output contains a row for:    
- psi: the overall mixture effect for the specified counterfactual profile    
- eta.high: the estimated value of the outcome at the upper limit of the counterfactual profile   
- eta.low: the predicted value of the outcome at the lower limit of the counterfactual profile   
