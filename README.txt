# CPE
Codes used for the paper "Conditional partial exchangeability: a probabilistic framework for multi-view clustering" by Beatrice Franzolini, Maria De Iorio, and Johan Eriksson
####
Conditional_t-HDP.R ->sorce code to conditional MCMC algorithms
Enriched_NNIG_uni.R ->mcmc for enriched Dirichlet process with 10 layers normal kernel and normal inverse gamma base measure
functions_to_overwrite_to_extract_cluster_configurations_fromLSBP.R -> used by Simulation_n3
HDP.R -> MCMC code for a standard HDP 
main_metabolomics.R -> to run to get results on real data
Simulation_n1.R -> to run to get results of the simulation study for Scenario n.1
Simulation_n2.R -> to run to get results of the simulation study for Scenario n.2
Simulation_n3.R -> to run to get results of the simulation study for Scenario n.3
Simulation_n4.R -> to run to get results of the simulation study for Scenario n.4
Simulation_n5_missspec.R -> to run to get results of the simulation study for Scenario n.5
Simulation_triangular.R -> simulation not used in the paper
telescopic_HDP_NNIG_uni.R -> MCMC for t-HDP markovian Normal-inverse-Gamma prior
telescopic_HDP_NNIW_multi.R -> MCMC for t-HDP markovian Normal-inverse-Wishart prior
telescopic_HDP_NNIW_multi_3L1P.R -> MCMC for t-HDP triangular Normal-inverse-Wishart prior
telescopic_HDP_NNIX_multi.R -> MCMC for t-HDP markovian Normal-inverse-Chisquared prior
telescopic_HDP_NNIX_multi_3L1P.R -> MCMC for t-HDP triangular Normal-inverse-Chisquared prior
telescopic_HDP_NN_uni.R -> MCMC for t-HDP markovian Normal prior
telescopic_HDP_NN_uni_3L1P.R -> MCMC for t-HDP triangular Normal prior
toy_example_dimensionality.R -> to run to get Figure 1 with toy-example
