# The following code simulate data (scenario 1) 
# and estimate k-means, the t-HDP and the LSBP
# (note that estimate of Enriched and DP-t-HDP are obtained running separate
# codes, named specifically "Enriched_NNIG_uni.R" and 
# "Comparison_VariantwithDP.R")
# results obtained from this code are presented in Table 1 of Section 7.1
# and Figure S5.5
# and Section S5 and S6 of the Supplement
# of the paper.

rm(list = ls())

library(salso)  # version 0.3.35          to comupte psm and point estimate
library(fossil) # version 0.4.0           to compute rand indexes
library(coda)   # version 0.19-4.1        to compute the effective sample size
library(factoextra)# version              to compute gap stat for kmeans

library(cowplot)  # version 1.1.3         to plot
library(ggpubr)   # version 0.6.0         to plot
library(dplyr)    # version 1.1.4         to plot
library(ggplot2)  # version 3.5.0         to plot
library(GGally)   # version 2.2.1         to plot

library(progress)      # version 1.2.3    to draw the progress bar
library(mvtnorm)       # version 1.2-4    for multivariate normal density
library(LaplacesDemon) # version 16.1.6   for inverse Wishart


library(rstudioapi) # version 0.15.0
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

devtools::install_github("tommasorigon/LSBP")
library(LSBP)
library(Formula) 

#source mcmc codes for telescopic clustering models:
source("Conditional_t-HDP.R")

#if true the code performs 100k iterations of the mcmc, 
#if false the code upload the outcome of the mcmc from a csv
run_MCMC = FALSE 


##SCENARIO N.1 ##########################################################################
#simulate data
set.seed(1)
n = 200
L = 10
true_layer = matrix(NA, nrow = n, ncol = L)
data = true_layer
true_layer[,1] = c(rep(1,floor(n/2)), rep(2,n - floor(n/2)))
data[true_layer[,1]==1, 1] = rnorm(sum(true_layer[,1]==1), 0, 1)
data[true_layer[,1]==2, 1] = rnorm(sum(true_layer[,1]==2), 4, 1)

for(l in 2:L){
  true_layer[,l] = true_layer[,l-1]
  moving = sample(1:n, floor(n*0.05))
  true_layer[moving,l] = true_layer[moving,l]%%2 + 1
  data[true_layer[,l]==1, l] = rnorm(sum(true_layer[,l]==1), 0, 1)
  data[true_layer[,l]==2, l] = rnorm(sum(true_layer[,l]==2), 4, 1)
  print(adj.rand.index(true_layer[,l], true_layer[,l-1]))
}

data = scale(data, scale = FALSE)
#ggpairs(as.data.frame(data))+ theme_bw()
data_plot = as.data.frame(data[,seq(1,10,2)])
colnames(data_plot) = c("Layer 1", "Layer 3", "Layer 5", "Layer 7","Layer 9")
print(ggpairs(data_plot, upper = list(continuous = "density"), 
        aes(color = as.factor(true_layer[,1]), 
            shape = as.factor(true_layer[,1])))+ theme_classic())
print("Figure S5.4 just printed")

#k-means 
#fviz_nbclust(as.matrix(data[,l]), kmeans, method = "gap", k.max = 5)
print("Table 1 and Table S6.1 entries for K-MEANS")
set.seed(0)
m = matrix(NA, nrow = L, ncol = n) #clustering configuration
for(l in 1:L){
  m[l, ] = kmeans(as.matrix(data[,l]), centers = 2)$cluster
}
for(l in 1:L){
  print(paste("Rand Index k-means Layer n.", l ,": ",
              rand.index(m[l, ],true_layer[,l])))
}
for(l in 1:L){
  if(sum(m[l,]%in%c(1,2))!=200){print(paste("label switching problem at layer",
                                            l))}
  mist = sum(m[l, ]!=true_layer[,l])
  print(paste("Mistakes Layer n.", l ,": ", 
              min(mist, n - mist )))
}

#estimate telescopic - conditional ##################################
totiter = 100000; burnin = 50000
if(run_MCMC){
  set.seed(1)
  m_saved =  telescopic_HDP_NNIG_uni(data, totiter = totiter)
  
  write.table(m_saved[,1,], file = paste("021_layer_1.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,2,], file = paste("021_layer_2.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,3,], file = paste("021_layer_3.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,4,], file = paste("021_layer_4.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,5,], file = paste("021_layer_5.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,6,], file = paste("021_layer_6.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,7,], file = paste("021_layer_7.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,8,], file = paste("021_layer_8.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,9,], file = paste("021_layer_9.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,10,], file = paste("021_layer_10.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
}else{
  pb <- progress_bar$new(
    format = " MCMC output [:bar] :percent Estimated completion time: :eta",
    total = L, clear = FALSE, width= 100)
  n_tot = 200
  m_saved = array(NA,c(totiter, L, n))
  for(l in 1:L){
    temp = as.matrix(read.table(paste("output/021_layer_",l,".csv", sep="")))
    m_saved[,l,] = temp
    pb$tick()
  }
}

# psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
# psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
# psm3 = psm(m_saved[(totiter-burnin+1):totiter,3,])
# psm4 = psm(m_saved[(totiter-burnin+1):totiter,4,])
# psm5 = psm(m_saved[(totiter-burnin+1):totiter,5,])
# psm6 = psm(m_saved[(totiter-burnin+1):totiter,6,])
# psm7 = psm(m_saved[(totiter-burnin+1):totiter,7,])
# psm8 = psm(m_saved[(totiter-burnin+1):totiter,8,])
# psm9 = psm(m_saved[(totiter-burnin+1):totiter,9,])
# psm10 = psm(m_saved[(totiter-burnin+1):totiter,10,])

#compute minVI
print("computing point estimate")
VIminlayer1 = salso(m_saved[(totiter-burnin+1):totiter,1,])
VIminlayer2 = salso(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer3 = salso(m_saved[(totiter-burnin+1):totiter,3,])
VIminlayer4 = salso(m_saved[(totiter-burnin+1):totiter,4,])
VIminlayer5 = salso(m_saved[(totiter-burnin+1):totiter,5,])
VIminlayer6 = salso(m_saved[(totiter-burnin+1):totiter,6,])
VIminlayer7 = salso(m_saved[(totiter-burnin+1):totiter,7,])
VIminlayer8 = salso(m_saved[(totiter-burnin+1):totiter,8,])
VIminlayer9 = salso(m_saved[(totiter-burnin+1):totiter,9,])
VIminlayer10 = salso(m_saved[(totiter-burnin+1):totiter,10,])

# write.table(VIminlayer1, file = paste("layer_1_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer2, file = paste("layer_2_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer3, file = paste("layer_3_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer4, file = paste("layer_4_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer5, file = paste("layer_5_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer6, file = paste("layer_6_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer7, file = paste("layer_7_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer8, file = paste("layer_8_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer9, file = paste("layer_9_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer10, file = paste("layer_10_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)

# effective sample size and mixing
par(mfrow=c(2,1))
print("Simulation study n.1, Table S8.1 Layer 1,5,10:")
pb <- progress_bar$new(
  format = " rand index per iter l=50[:bar] :percent Estimated completion time: :eta",
  total = (totiter - burnin + 1)*3, clear = FALSE, width= 100)
for (l in c(1,5,10)){
  ri = NULL
  for (iter in burnin:totiter){
    ri[iter-burnin+1] = c(ri, rand.index(m_saved[iter,l,], true_layer[,l]))
    pb$tick()
  }
  # plot(ri, type = "l", xlab="MCMC iterations")
  # abline(v = burnin, col="red", lwd=3, lty=2)
  # title(main=paste("Rand Index Layer ",l," - Trace plot \n Conditional algorithm"),
  #                cex.lab=0.75)
  print(effectiveSize(ri/(totiter-burnin) ))
}

#produce Figure S5.5 in the supplement:
par(mfrow=c(3,1))
print("computing rand index")
pb <- progress_bar$new(
  format = " Computing rand index per iteration [:bar] :percent Estimated completion time: :eta",
  total = (totiter-burnin)*3, clear = FALSE, width= 100)
for (l in c(1,5,10)){
  ri = NULL
  for (iter in (burnin+1):totiter){
    ri = c(ri, rand.index(m_saved[iter,l,], true_layer[,l]))
    pb$tick()
  }
  hist(ri, xlab=paste(" "), ylab=paste(" "), xlim=c(0,1), col=l,
       main=paste("Posterior distribution Rand Index Layer n.",l),
       cex.lab=0.75)
}
print("Figure S5.5 just printed")


print("Table 1 and Table S6.1 entries for t-HDP")
print(paste("Rand Index Layer n.", 1 ,": ", 
            rand.index(VIminlayer1,true_layer[,1])))
print(paste("Rand Index Layer n.", 2 ,": ", 
            rand.index(VIminlayer2,true_layer[,2])))
print(paste("Rand Index Layer n.", 3 ,": ", 
            rand.index(VIminlayer3,true_layer[,3])))
print(paste("Rand Index Layer n.", 4 ,": ", 
            rand.index(VIminlayer4,true_layer[,4])))
print(paste("Rand Index Layer n.", 5 ,": ", 
            rand.index(VIminlayer5,true_layer[,5])))
print(paste("Rand Index Layer n.", 6 ,": ", 
            rand.index(VIminlayer6,true_layer[,6])))
print(paste("Rand Index Layer n.", 7 ,": ", 
            rand.index(VIminlayer7,true_layer[,7])))
print(paste("Rand Index Layer n.", 8 ,": ", 
            rand.index(VIminlayer8,true_layer[,8])))
print(paste("Rand Index Layer n.", 9 ,": ", 
            rand.index(VIminlayer9,true_layer[,9])))
print(paste("Rand Index Layer n.", 10 ,": ", 
            rand.index(VIminlayer10,true_layer[,10])))


print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3!=true_layer[,3]) ))
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4!=true_layer[,4]) ))
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5!=true_layer[,5]) ))
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6!=true_layer[,6]) ))
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7!=true_layer[,7]) ))
VIminlayer8[VIminlayer8==1] = 99; VIminlayer8[VIminlayer8==2] = 1
VIminlayer8[VIminlayer8==99] = 2
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8!=true_layer[,8]) )) #label switching
VIminlayer9[VIminlayer9==1] = 99; VIminlayer9[VIminlayer9==2] = 1
VIminlayer9[VIminlayer9==99] = 2
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9!=true_layer[,9]) )) #label switching
VIminlayer10[VIminlayer10==1] = 99; VIminlayer10[VIminlayer10==2] = 1
VIminlayer10[VIminlayer10==99] = 2
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10!=true_layer[,10]) )) #label switching



################################################################################
#logit stick breaking 

#modify functions to get cluster allocations: 
source("functions_to_overwrite_to_extract_cluster_configurations_fromLSBP.R")


environment(my_prior_LSBP) <- asNamespace('LSBP')
assignInNamespace("prior_LSBP", my_prior_LSBP, ns = "LSBP")

environment(my_LSBP_Gibbs_univ) <- asNamespace('LSBP')
assignInNamespace("LSBP_Gibbs_univ", my_LSBP_Gibbs_univ, ns = "LSBP")

environment(my_LSBP_Gibbs_multi) <- asNamespace('LSBP')
assignInNamespace("LSBP_Gibbs_multi", my_LSBP_Gibbs_multi, ns = "LSBP")

environment(LSBP_Gibbs) <- asNamespace('LSBP')
assignInNamespace("LSBP_Gibbs", my_LSBP_Gibbs, ns = "LSBP")



data_flat = data.frame( X = (c(rep(1,n), rep(2,n), rep(3,n),
                         rep(4,n), rep(5,n), rep(6,n),
                         rep(7,n), rep(8,n), rep(8,n),
                         rep(10,n))), Y = as.numeric(data))

if(run_MCMC){
  set.seed(0)
  R <- 10000 #iterations
  burnin <- 2000
  fit_gibbs <- LSBP_Gibbs(Y ~ X|X, data=data_flat, H=5, 
                          control=control_Gibbs(R=R,burn_in=burnin,
                                                method_init="random") )
  clust_LSBP = fit_gibbs$param$G
  write.table(clust_LSBP, file = paste("061_layer_1to10.csv",sep=""), 
              row.names = FALSE,
              col.names= FALSE)
}else{
  clust_LSBP = as.matrix(read.table(paste("output/061_layer_1to10.csv", sep="")))
}

# psm1 = psm(fit_gibbs$param$G[,1:200 ])
# psm2 = psm(fit_gibbs$param$G[,201:400 ])
# psm3 = psm(fit_gibbs$param$G[,401:600 ])
# psm4 = psm(fit_gibbs$param$G[,601:800 ])
# psm5 = psm(fit_gibbs$param$G[,801:1000 ])
# psm6 = psm(fit_gibbs$param$G[,1001:1200 ])
# psm7 = psm(fit_gibbs$param$G[,1201:1400 ])
# psm8 = psm(fit_gibbs$param$G[,1401:1600 ])
# psm9 = psm(fit_gibbs$param$G[,1601:1800 ])
# psm10 = psm(fit_gibbs$param$G[,1801:2000 ])


VIminlayer1_LSBP = salso(clust_LSBP[,1:200], loss = VI())
VIminlayer2_LSBP = salso(clust_LSBP[,201:400], loss = VI())
VIminlayer3_LSBP = salso(clust_LSBP[,401:600], loss = VI())
VIminlayer4_LSBP = salso(clust_LSBP[,601:800], loss = VI())
VIminlayer5_LSBP = salso(clust_LSBP[,801:1000], loss = VI())
VIminlayer6_LSBP = salso(clust_LSBP[,1001:1200], loss = VI())
VIminlayer7_LSBP = salso(clust_LSBP[,1201:1400], loss = VI())
VIminlayer8_LSBP = salso(clust_LSBP[,1401:1600], loss = VI())
VIminlayer9_LSBP = salso(clust_LSBP[,1601:1800], loss = VI())
VIminlayer10_LSBP = salso(clust_LSBP[,1801:2000], loss = VI())

# VIminlayer1_LSBP = minVI(psm1, fit_gibbs$param$G[,1:200], method = "all")$cl[1,]
# VIminlayer2_LSBP = minVI(psm2, fit_gibbs$param$G[,201:400], method = "all")$cl[1,]
# VIminlayer3_LSBP = minVI(psm3, fit_gibbs$param$G[,401:600], method = "all")$cl[1,]
# VIminlayer4_LSBP = minVI(psm4, fit_gibbs$param$G[,601:800], method = "all")$cl[1,]
# VIminlayer5_LSBP = minVI(psm5, fit_gibbs$param$G[,801:1000], method = "all")$cl[1,]
# VIminlayer6_LSBP = minVI(psm6, fit_gibbs$param$G[,1001:1200], method = "all")$cl[1,]
# VIminlayer7_LSBP = minVI(psm7, fit_gibbs$param$G[,1201:1400], method = "all")$cl[1,]
# VIminlayer8_LSBP = minVI(psm8, fit_gibbs$param$G[,1401:1600], method = "all")$cl[1,]
# VIminlayer9_LSBP = minVI(psm9, fit_gibbs$param$G[,1601:1800], method = "all")$cl[1,]
# VIminlayer10_LSBP = minVI(psm10, fit_gibbs$param$G[,1801:2000], method = "all")$cl[1,]


print("Values Table 1 and Table S6.1 - LSBP")
print(paste("Rand Index Layer n.", 1 ,": ", 
            rand.index(VIminlayer1_LSBP,true_layer[,1])))
print(paste("Rand Index Layer n.", 2 ,": ", 
            rand.index(VIminlayer2_LSBP,true_layer[,2])))
print(paste("Rand Index Layer n.", 3 ,": ", 
            rand.index(VIminlayer3_LSBP,true_layer[,3])))
print(paste("Rand Index Layer n.", 4 ,": ", 
            rand.index(VIminlayer4_LSBP,true_layer[,4])))
print(paste("Rand Index Layer n.", 5 ,": ", 
            rand.index(VIminlayer5_LSBP,true_layer[,5])))
print(paste("Rand Index Layer n.", 6 ,": ", 
            rand.index(VIminlayer6_LSBP,true_layer[,6])))
print(paste("Rand Index Layer n.", 7 ,": ", 
            rand.index(VIminlayer7_LSBP,true_layer[,7])))
print(paste("Rand Index Layer n.", 8 ,": ", 
            rand.index(VIminlayer8_LSBP,true_layer[,8])))
print(paste("Rand Index Layer n.", 9 ,": ", 
            rand.index(VIminlayer9_LSBP,true_layer[,9])))
print(paste("Rand Index Layer n.", 10 ,": ", 
            rand.index(VIminlayer10_LSBP,true_layer[,10])))

print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1_LSBP!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2_LSBP!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3_LSBP!=true_layer[,3]) ))
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4_LSBP!=true_layer[,4]) ))
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5_LSBP!=true_layer[,5]) ))
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6_LSBP!=true_layer[,6]) ))
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7_LSBP!=true_layer[,7]) ))
VIminlayer8_LSBP[VIminlayer8_LSBP==1] = 99; VIminlayer8_LSBP[VIminlayer8_LSBP==2] = 1
VIminlayer8_LSBP[VIminlayer8_LSBP==99] = 2
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8_LSBP!=true_layer[,8]) ))
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9_LSBP!=true_layer[,9]) )) 
VIminlayer10_LSBP[VIminlayer10_LSBP==1] = 99; VIminlayer10_LSBP[VIminlayer10_LSBP==2] = 1
VIminlayer10_LSBP[VIminlayer10_LSBP==99] = 2
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10_LSBP!=true_layer[,10]) )) #label switching

