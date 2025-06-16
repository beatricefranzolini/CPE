# Simulation scenario C 
# The following code simulate data (scenario C) 
# and estimate k-mean and the t-HDP
# results obtained from this code are presented 
# in Section S.5 of the Supplementary
# Figure S5.6 and S5.7 

#Warning 2 :the code work need a version of fastmap >= 1.2.0 
#if version of fastmap < 1.2.0, please reinstall it, uploading and then rerun
#the following code 

rm(list = ls())

library(salso)  # version 0.3.35          to comupte psm 
library(fossil) # version 0.4.0           to compute rand indexes
library(coda)   # version 0.19-4.1        to compute the effective sample size
library(factoextra)# version 1.0.7        to compute gap stat for kmeans

library(cowplot)  # version 1.1.3         to plot
library(ggpubr)   # version 0.6.0         to plot
library(dplyr)    # version 1.1.4         to plot
library(ggplot2)  # version 3.5.0         to plot
library(GGally)   # version 2.2.1         to plot

library(progress)      # version 1.2.3    to draw the progress bar
library(mvtnorm)       # version 1.2-4    for multivariate normal density
library(LaplacesDemon) # version 16.1.6   for inverse Wishart
library(pheatmap)      # version 1.0.12

devtools::install_github("sarawade/mcclust.ext") #package for point estimates
library(mcclust.ext)

library(rstudioapi) # version 0.15.0
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source mcmc codes for telescopic clustering models:
source("Conditional_t-HDP.R")

#if true the code performs 100k iterations of the mcmc, 
#if false the code upload the outcome of the mcmc from a csv
run_MCMC = FALSE 


##SCENARIO N.4 ##########################################################################
#simulate data
set.seed(1)
sigma1 = MCMCpack::riwish(10, 0.1 *diag(2))
sigma2 = MCMCpack::riwish(10, 0.1* diag(2))
sigma3 = MCMCpack::riwish(10, 0.1 *diag(2))
n = 200
L = 2 #bivariate
true_layer = matrix(NA, nrow = n, ncol = L)
data = matrix(NA, nrow = n, ncol = L*2)
true_layer[,1] = c(rep(1,floor(n/3)), rep(2,floor(n/3)), rep(3,n - 2*floor(n/3)))
data[true_layer[,1]==1, c(1,2)] = mvtnorm::rmvt(sum(true_layer[,1]==1), sigma = sigma1)
data[true_layer[,1]==2, c(1,2)] = mvtnorm::rmvt(sum(true_layer[,1]==2), sigma = sigma2) +4
data[true_layer[,1]==3, c(1,2)] = mvtnorm::rmvt(sum(true_layer[,1]==3), sigma = sigma3) +8

l = 2
true_layer[,l] = true_layer[,l-1]
moving = sample(1:(2*floor(n/3)), floor((2*floor(n/3))*0.05), replace = FALSE)
true_layer[moving,l] = true_layer[moving,l]%%2 + 1
data[true_layer[,l]==1, c(3,4)] = mvtnorm::rmvt(sum(true_layer[,l]==1), sigma = sigma1)
data[true_layer[,l]==2, c(3,4)] = mvtnorm::rmvt(sum(true_layer[,l]==2), sigma = sigma2) +4
data[true_layer[,l]==3, c(3,4)] = mvtnorm::rmvt(sum(true_layer[,l]==3), sigma = sigma3) +8
print(adj.rand.index(true_layer[,l], true_layer[,l-1]))

data_plot = as.data.frame(data)
colnames(data_plot) = c("Layer 1 - v1", "Layer 1 - v2", "Layer 2 - v1", "Layer 2 - v2")
print(ggpairs(data_plot, upper = list(continuous = "density"),
              aes(color = as.factor(true_layer[,1]), 
              shape = as.factor(true_layer[,1])))+ theme_classic())
print("Figure S5.6 just printed")

data = scale(data)

#k-means 
fviz_nbclust(as.matrix(data[,1:2]), kmeans, method = "gap", k.max = 5)
fviz_nbclust(as.matrix(data[,3:4]), kmeans, method = "gap", k.max = 5)

set.seed(0)
m = matrix(NA, nrow = L, ncol = n) #clustering configuration
m[1, ] = stats::kmeans(as.matrix(data[,1:2]), centers = 3)$cluster
m[2, ] = stats::kmeans(as.matrix(data[,3:4]), centers = 1)$cluster
for(l in 1:L){
  print(paste("Rand Index k-means Layer n.", l ,": ",
              rand.index(m[l, ],true_layer[,l])))
}
for(l in 1:L){
  if(sum(m[l,]%in%c(1,2,3))!=200){print(paste("label switching problem at layer", l))}
  mist = sum(m[l, ]!=true_layer[,l])
  print(paste("Mistakes Layer n.", l ,": ", 
              min(mist, n - mist )))
}


#estimate telescopic - conditional ##################################
colayer = c(1,1,2,2)

if(run_MCMC){
  totiter = 100000
  burnin = 50000
  set.seed(1)
  m_saved = telescopic_HDP_NNIX_multi(data, colayer, alpharandom = TRUE, totiter = totiter)
  write.table(m_saved[,1,], file = paste("02C_layer_1.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,2,], file = paste("02C_layer_2.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
}else{
  print("MCMC output")
  totiter = 100000; burnin = 50000
  m_saved = array(NA,c(totiter, L, n))
  pb <- progress_bar$new(
    format = " MCMC output [:bar] :percent Estimated completion time: :eta",
    total = L, clear = FALSE, width= 100)
  for(l in 1:L){
    temp = as.matrix(read.table(paste("output/02C_layer_",l,".csv", sep="")))
    m_saved[,l,] = temp
    pb$tick()
  }
}
print("computing point estimate")
psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer1 = minVI(psm1, m_saved[(totiter-burnin+1):totiter,1,], method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, m_saved[(totiter-burnin+1):totiter,2,], method = "all")$cl[1,]
# write.table(VIminlayer1, file = paste("layer_1_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer2, file = paste("layer_2_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)

print(rand.index(VIminlayer1, true_layer[,1]))
print(rand.index(VIminlayer2, true_layer[,2]))

#effective sample size, mixing and uncertainty 
ri = matrix(NA, nrow = totiter-burnin + 1, ncol = L)
pb <- progress_bar$new(
  format = " rand index per iter l=1[:bar] :percent Estimated completion time: :eta",
  total = (totiter - burnin + 1)*2, clear = FALSE, width= 100)
for (l in 1:L){
  for (iter in burnin:totiter){
    ri[iter-burnin+1, l] = rand.index(m_saved[iter,l,], true_layer[,l])
    pb$tick()
  }
  #print(effectiveSize(ri[(burnin+1):totiter, l])/(totiter-burnin) )
}

# par(mfrow=c(2,1))
# for (l in 1:L){
#   plot(ri[,l], type = "l", xlab="MCMC iterations")
#   abline(v = burnin, col="red", lwd=3, lty=2)
#   title(main=paste("Rand Index Layer ",l," - Trace plot \n Conditional algorithm"),
#       cex.lab=0.75)
# }

par(mfrow=c(2,1))
for (l in 1:L){
  hist(ri[,l], xlab=paste(" "), ylab=paste(" "), xlim=c(0,1), col=l,
       main=paste("Posterior distribution Rand Index Layer n.",l),
       cex.lab=0.75)
  abline(v = rand.index(m[l, ],true_layer[,l]), col="red", lwd=3, lty=2)
}
print("Figure S5.7 just printed")
