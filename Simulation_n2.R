# Simulation scenario 2 with 100 layers
# The following code simulate data (scenario 2) 
# and estimate k-mean and the t-HDP
# results obtained from this code are presented in Section 7.1 of the paper 
# Figure 4 and Figure 5 and 
# Table S8.1 in the Supplementary

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


##SCENARIO N.2 ##########################################################################
#simulate data
print("SIMULATION N.2")
set.seed(1)
n = 200
L = 100
true_layer = matrix(NA, nrow = n, ncol = L)
data = true_layer
true_layer[,1] = c(rep(1,floor(n/2)), rep(2,n - floor(n/2)))
data[true_layer[,1]==1, 1] = rnorm(sum(true_layer[,1]==1), 0, 1)
data[true_layer[,1]==2, 1] = rnorm(sum(true_layer[,1]==2), 3, 1)

for(l in 2:L){
  true_layer[,l] = true_layer[,l-1]
  moving = sample(1:n, floor(n*0.02))
  true_layer[moving,l] = true_layer[moving,l]%%2 + 1
  data[true_layer[,l]==1, l] = rnorm(sum(true_layer[,l]==1), 0, 1)
  data[true_layer[,l]==2, l] = rnorm(sum(true_layer[,l]==2), 3, 1)
  #print(adj.rand.index(true_layer[,l], true_layer[,l-1]))
}

print("compute rand index between any pair of layers")
pb <- progress_bar$new(
  format = " RI [:bar] :percent Estimated completion time: :eta",
  total = L*L, clear = FALSE, width= 100)
RI_matrix = matrix(0, nrow = L, ncol = L)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix[l,ll] = rand.index(true_layer[,l], true_layer[,ll])
    pb$tick()
  }
}
RI_matrix = t(RI_matrix) + RI_matrix

#kmeans 
num_clust_kmeans = rep(NA, L)
b <- progress_bar$new(
  format = " k-means num of clusters [:bar] :percent Estimated completion time: :eta",
  total = L, clear = FALSE, width= 100)
for(l in 1:L){
  temp = fviz_nbclust(as.matrix(data[,l]), kmeans, method = "gap", k.max = 5)
  num_clust_kmeans[l] = seq(1,5)[temp$data$gap==max(temp$data$gap)]
  #print(num_clust_kmeans[l])
  pb$tick()
}
set.seed(0)
m = matrix(NA, nrow = L, ncol = n) #clustering configuration
performance_kmeans = rep(NA, L)
pb <- progress_bar$new(
  format = " k-means solution [:bar] :percent Estimated completion time: :eta",
  total = L, clear = FALSE, width= 100)
for(l in 1:L){
  m[l, ] = stats::kmeans(as.matrix(data[,l]), centers = num_clust_kmeans[l])$cluster
  performance_kmeans[l] = rand.index(m[l, ],true_layer[,l])
  pb$tick()
}

if(run_MCMC){
  totiter = 70000; burnin = 20000
  set.seed(1)
  telescopic_HDP_NNIG_uni(data, alpharandom= TRUE, totiter = totiter)
  # for(l in 1:L){
  #   write.table(m_saved[,l,], file = paste("022_layer_",l,".csv", sep=""), 
  #             row.names = FALSE,
  #             col.names= FALSE)
  # }
}else{
  print("MCMC output")
  totiter = 73328; burnin = 20000
  m_saved = array(NA,c(totiter, L, n))
  pb <- progress_bar$new(
    format = " MCMC output [:bar] :percent Estimated completion time: :eta",
    total = L, clear = FALSE, width= 100)
  for(l in 1:L){
    temp = as.matrix(read.table(paste("output/022_layer_",l,".csv", sep="")))
    m_saved[,l,] = temp
    pb$tick()
  }
}

# psm = array(NA,c(L, n, n)) 
# for(l in 1:L){
#   psm[l,,] = psm(m_saved[(totiter-burnin+1):totiter,l,])
#   print(l)
# }

VIminlayer = array(NA,c(L, n)) 
pb <- progress_bar$new(
  format = " point estimate [:bar] :percent Estimated completion time: :eta",
  total = L, clear = FALSE, width= 100)
for(l in 1:L){
  #VIminlayer[l,] = minVI(psm[l,,], m_saved[(totiter-burnin+1):totiter,l,], method = "all")$cl[1,]
  VIminlayer[l,] = salso(m_saved[(totiter-burnin+1):totiter,l,])
  pb$tick()
}

performance = rep(NA, L)
pb <- progress_bar$new(
  format = " rand index per layer [:bar] :percent Estimated completion time: :eta",
  total = L, clear = FALSE, width= 100)
for(l in 1:L){
  performance[l] = rand.index(VIminlayer[l,],true_layer[,l])
  pb$tick()
}


print("Simulation study n.2, Table S8.1 Layer 1,50,100:")
performance1 = rep(NA, totiter-burnin)
pb <- progress_bar$new(
  format = " rand index per iter l=1[:bar] :percent Estimated completion time: :eta",
  total = totiter - burnin + 1, clear = FALSE, width= 100)
for(iter in burnin:totiter){
  performance1[iter-burnin+1] = rand.index(true_layer[,1],m_saved[iter,1,])
  pb$tick()
}

print(effectiveSize(performance1)/(totiter-burnin)) 

performance50 = rep(NA, totiter-burnin)
pb <- progress_bar$new(
  format = " rand index per iter l=50[:bar] :percent Estimated completion time: :eta",
  total = totiter - burnin + 1, clear = FALSE, width= 100)
for(iter in burnin:totiter){
  performance50[iter-burnin+1] = rand.index(true_layer[,50],m_saved[iter,50,])
  pb$tick()
}

print(effectiveSize(performance50)/(totiter-burnin)) 

performance100 = rep(NA, totiter-burnin)
pb <- progress_bar$new(
  format = " rand index per iter l=100[:bar] :percent Estimated completion time: :eta",
  total = totiter - burnin + 1, clear = FALSE, width= 100)
for(iter in burnin:totiter){
  performance100[iter-burnin+1] = rand.index(true_layer[,100],m_saved[iter,100,])
  pb$tick()
}

print(effectiveSize(performance100)/(totiter-burnin) )

results_plot = data.frame(Layer = c(rep(1:100,2)),
                          RandIndex = c(performance, performance_kmeans), 
                          Model = factor(c(rep("t-HDP",100), rep("k-means", 100))))

print(ggplot(results_plot, aes(x = Layer, y = RandIndex)) + 
  geom_line(aes(color = Model, linetype = Model), linewidth=0.8) + 
  scale_color_manual(values = c("blue", "red"))+
  theme_classic() + ylim(0,1) + ylab("Rand Index"))
print("Figure 4 just printed")

print(pheatmap(RI_matrix, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA))
print("Figure 5 - panel A just printed")

RI_matrix_est = matrix(0, nrow = L, ncol = L)
pb <- progress_bar$new(
  format = " RI t-HDP[:bar] :percent Estimated completion time: :eta",
  total = L*L, clear = FALSE, width= 100)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix_est [l,ll] = rand.index(VIminlayer[l,], VIminlayer[ll,])
    pb$tick()
  }
}
RI_matrix_est  = t(RI_matrix_est ) + RI_matrix_est 
print(pheatmap(RI_matrix_est, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA))
print("")
print("Figure 5 - panel B just printed")

RI_matrix_est_kmeans = matrix(0, nrow = L, ncol = L)
pb <- progress_bar$new(
  format = " RI k-means[:bar] :percent Estimated completion time: :eta",
  total = L*L, clear = FALSE, width= 100)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix_est_kmeans [l,ll] = rand.index(m[l, ], m[ll, ])
    pb$tick()
  }
}

RI_matrix_est_kmeans  = t(RI_matrix_est_kmeans ) + RI_matrix_est_kmeans 
print(pheatmap(RI_matrix_est_kmeans, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA))
print("")
print("Figure 5 - panel C just printed")