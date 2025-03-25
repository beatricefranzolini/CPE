# CPE

R codes used for the paper "Conditional partial exchangeability: a probabilistic framework for multi-view clustering" by Franzolini, De Iorio, Eriksson.

**Authors**: [Beatrice Franzolini](https://beatricefranzolini.github.io)

#### Overview 

The repository contains the following:

`Comparison_VariantwithDP.R` -> 	to run to get the results in Section S6 of the Supplement (simulation studies A, B, and, 1 and DP-t-HDP model)

`Conditional_t-HDP.R` ->		source code to conditional MCMC algorithms 

`Enriched_NNIG_uni.R` ->		MCMC for enriched Dirichlet process with 10 layers normal kernel and Normal-Inverse-Gamma base measure

`functions_to_overwrite_to_extract_cluster_configurations_fromLSBP.R` -> used by Simulation_n3.R

`HDP.R' -> 			MCMC code for a standard HDP (as in Teh et. al)  

`main_metabolomics.R` -> 		to run to get results on real data as presented in Section 7.2 of the main paper and Section S7 of the Supplement (real data and t-HDP model)

`Simulation_n1.R` -> 		to run to get results of the simulation study for Scenario n.A (results in Section S5 of the Supplement)

`Simulation_n2.R` -> 		to run to get results of the simulation study for Scenario n.B (results in Section S5 of the Supplement)

`Simulation_n3.R` -> 		to run to get results of the simulation study for Scenario n.1 (results in Section 7.1 of the main paper and Section S5 of the Supplement)

`Simulation_n4.R` -> 		to run to get results of the simulation study for Scenario n.2 (results in Section 7.1 of the main paper and Section S5 of the Supplement)

`Simulation_n5_missspec.R` -> 	to run to get results of the simulation study for Scenario n.3

`Simulation_triangular.R` -> 	simulation not used in the paper

`telescopic_HDP_NN_uni.R` -> 		MCMC for t-HDP Markovian Normal prior

`telescopic_HDP_NN_uni_3L1P.R` -> 	MCMC for t-HDP triangular Normal prior

`telescopic_HDP_NN_uni_VariantwithDP` ->	MCMC for DP-t-HDP Markovian Normal prior

`telescopic_HDP_NNIG_uni.R` -> 		MCMC for t-HDP Markovian Normal-Inverse-Gamma prior

`telescopic_HDP_NNIG_uni_VariantwithDP` ->MCMC for DP-t-HDP Markovian Normal-Inverse-Gamma prior

`telescopic_HDP_NNIW_multi.R` -> 		MCMC for t-HDP Markovian Normal-Inverse-Wishart prior

`telescopic_HDP_NNIW_multi_3L1P.R` -> 	MCMC for t-HDP triangular Normal-Inverse-Wishart prior

`telescopic_HDP_NNIX_multi.R` -> 		MCMC for t-HDP Markovian Normal-Inverse-Chisquared prior

`telescopic_HDP_NNIX_multi_3L1P.R` -> 	MCMC for t-HDP triangular Normal-Inverse-Chisquared prior

`toy_example_dimensionality.R` -> to run to get Figure 1 with toy-example (Section 1 of the main paper)


#### Questions or bugs
For bug reporting purposes, e-mail [Beatrice Franzolini](https://beatricefranzolini.github.io) (franzolini@pm.me).




