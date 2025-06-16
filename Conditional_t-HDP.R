#t-HDP _ conditional algorithm 

#A. Markovian dependence across L layers
#A.1 univariate layers
source("telescopic_HDP_NN_uni.R")   #used in simulation n.A and n.B
source("telescopic_HDP_NNIG_uni.R") #used in simulation n.1 and n.3
#A.2 multivariate layers
source("telescopic_HDP_NNIX_multi.R") #used in simulation n.2
source("telescopic_HDP_NNIW_multi.R") #not used

#B. Triangular structure with 3 layers
#B.1 univariate layers
source("telescopic_HDP_NN_uni_3L1P.R") #not used
#B.2 multivariate layers
source("telescopic_HDP_NNIW_multi_3L1P.R") #not used
source("telescopic_HDP_NNIX_multi_3L1P.R") #used in metabolomics application

#C. Variant with DP at first layer 
source("telescopic_HDP_NN_uni_VariantwithDP.R")  #used in simulation n.A and n.B
source("telescopic_HDP_NNIG_uni_VariantwithDP.R")#used in simulation n.1
