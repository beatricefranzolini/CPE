# CPE

R codes used for the paper "Conditional partial exchangeability: a probabilistic framework for multi-view clustering" by Franzolini, De Iorio, Eriksson.

**Author contact**: [Beatrice Franzolini](https://beatricefranzolini.github.io)

## Overview 

This repository contains the R code used for simulations and data analysis in the manuscript and supplementary materials. 
The codes include the implementation of posterior MCMC for telescopic hierarchical Dirichlet processes (t-HDP) for different kernel choices and poly-tree structures.

## Contents 

### Main scripts 

`Simulation_n1.R` to run to get results of the simulation study for Scenario n.1 (results in Section 7.1 of the main paper and Section S5 of the Supplement) (The code can be sourced in full or executed line by line.)

`Simulation_n2.R` to run to get results of the simulation study for Scenario n.2 (results in Section 7.1 of the main paper and Section S5 of the Supplement) (The code can be sourced in full or executed line by line.)

`Simulation_nA.R` to run to get results of the simulation study for Scenario n.A (results in Section S5 of the Supplement) (The code can be sourced in full or executed line by line.)

`Simulation_nB.R` to run to get results of the simulation study for Scenario n.B (results in Section S5 of the Supplement) (The code can be sourced in full or executed line by line.)

`Simulation_nC.R` to run to get results of the simulation study for Scenario n.3 (The code can be sourced in full or executed line by line.)

`Comparison_VariantwithDP.R` to run to get the results corresponding to the DP-t-HDP model (and simulated data Scenario A. B. and 1) presented in Section S6 of the Supplement (The code can be sourced in full or executed line by line.)

`Enriched_NNIG_uni.R` to run to get the results corresponding to the enriched Dirichlet process with 10 layers normal kernel and Normal-Inverse-Gamma base measure (and simulated data Scenario n.1) (presented in Section 7.1 of paper and Section S5 and S6 of the Supplement) (The code can be sourced in full or executed line by line.)

`main_metabolomics.R` to run to get results on real data as presented in Section 7.2 of the main paper and Section S7 of the Supplement (real data and t-HDP model)
WARNING: real data are not available, not all Figures are reproducible without real data!

### MCMC functions (These functions are called by main scripts.)

`HDP.R` MCMC for a standard HDP (as in Teh et. al, 2006) Normal prior

`telescopic_HDP_NN_uni.R` MCMC for t-HDP Markovian Normal prior

`telescopic_HDP_NN_uni_3L1P.R` MCMC for t-HDP triangular Normal prior

`telescopic_HDP_NN_uni_VariantwithDP` MCMC for DP-t-HDP Markovian Normal prior

`telescopic_HDP_NNIG_uni.R` MCMC for t-HDP Markovian Normal-Inverse-Gamma prior

`telescopic_HDP_NNIG_uni_VariantwithDP` MCMC for DP-t-HDP Markovian Normal-Inverse-Gamma prior

`telescopic_HDP_NNIW_multi.R` MCMC for t-HDP Markovian Normal-Inverse-Wishart prior

`telescopic_HDP_NNIW_multi_3L1P.R` MCMC for t-HDP triangular Normal-Inverse-Wishart prior

`telescopic_HDP_NNIX_multi.R` MCMC for t-HDP Markovian Normal-Inverse-Chisquared prior

`telescopic_HDP_NNIX_multi_3L1P.R` MCMC for t-HDP triangular Normal-Inverse-Chisquared prior

### Utility Scripts (These scripts are called by main script.)

`Conditional_t-HDP.R` code sourcing MCMC algorithms for the t-HDP models

`functions_to_overwrite_to_extract_cluster_configurations_fromLSBP.R` used in  `Simulation_n1.R` to extract cluster configurations from LSBP-based models.

### Output

`output` directory containing csv files with MCMC chains for all models and data + a pdf with the codes used to map csv to each pair of data and model

### Others

`Simulation_triangular.R` to run to perform a simulation not used in the paper

`toy_example_dimensionality.R` to run to get Figure 1 with toy-example (Section 1 of the main paper) (The code can be sourced in full or executed line by line.)

#### Questions or bugs
For bug reporting purposes, e-mail [Beatrice Franzolini](https://beatricefranzolini.github.io) (franzolini@pm.me).




