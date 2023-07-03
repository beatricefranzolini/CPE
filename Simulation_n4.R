#Simulation scenario 5 with 100 layers (true clustering)

library(mcclust.ext) #to minimize VI
library(T4cluster) #to compute psm with psm() function
library(fossil) #to compute rand indexes
library(coda) #to compute the effective sample size
library(factoextra)#for gap stat and fvz_nbclust

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(ggplot2)
library(GGally) #to plot multivariate continous data
library(pheatmap) #to plot

source("Conditional_t-HDP.R")

##SCENARIO N.5 ##########################################################################
#simulate data
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
  print(adj.rand.index(true_layer[,l], true_layer[,l-1]))
}

RI_matrix = matrix(0, nrow = L, ncol = L)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix[l,ll] = rand.index(true_layer[,l], true_layer[,ll])
  }
}
RI_matrix = t(RI_matrix) + RI_matrix

#kmeans 
num_clust_kmeans = rep(NA, L)
for(l in 1:L){
  temp = fviz_nbclust(as.matrix(data[,l]), kmeans, method = "gap", k.max = 5)
  num_clust_kmeans[l] = seq(1,5)[temp$data$gap==max(temp$data$gap)]
  print(num_clust_kmeans[l])
}
set.seed(0)
m = matrix(NA, nrow = L, ncol = n) #clustering configuration
performance_kmeans = rep(NA, L)
for(l in 1:L){
  m[l, ] = stats::kmeans(as.matrix(data[,l]), centers = num_clust_kmeans[l])$cluster
  performance_kmeans[l] = rand.index(m[l, ],true_layer[,l])
}

#totiter= 73328
totiter = 70000; burnin = 20000
set.seed(1)
telescopic_HDP_NNIG_uni(data, alpharandom= TRUE, totiter = totiter)
#load results:
m_saved = array(NA,c(totiter, L, n))
for(l in 1:L){
  temp = as.matrix(read.table(paste("Simulation_num4/layer",l,".csv", sep="")))
  m_saved[,l,] = temp
  print(l)
}


psm = array(NA,c(L, n, n)) 
for(l in 1:L){
  psm[l,,] = psm(m_saved[(totiter-burnin+1):totiter,l,])
  print(l)
}

VIminlayer = array(NA,c(L, n)) 
for(l in 1:L){
  VIminlayer[l,] = minVI(psm[l,,], m_saved[(totiter-burnin+1):totiter,l,], method = "all")$cl[1,]
  print(l)
}

performance = rep(NA, L)
for(l in 1:L){
  performance[l] = rand.index(VIminlayer[l,],true_layer[,l])
  print(l)
}

performance1 = rep(NA, totiter-burnin)
for(iter in burnin:totiter){
  performance1[iter-burnin+1] = rand.index(true_layer[,1],m_saved[iter,1,])
  print(iter)
}

effectiveSize(performance1)/(totiter-burnin) 

performance50 = rep(NA, totiter-burnin)
for(iter in burnin:totiter){
  performance50[iter-burnin+1] = rand.index(true_layer[,50],m_saved[iter,50,])
  print(iter)
}

effectiveSize(performance50)/(totiter-burnin) 

performance100 = rep(NA, totiter-burnin)
for(iter in burnin:totiter){
  performance100[iter-burnin+1] = rand.index(true_layer[,100],m_saved[iter,100,])
  print(iter)
}

effectiveSize(performance100)/(totiter-burnin) 

results_plot = data.frame(Layer = c(rep(1:100,2)),
                          RandIndex = c(performance, performance_kmeans), 
                          Model = factor(c(rep("t-HDP",100), rep("k-means", 100))))

ggplot(results_plot, aes(x = Layer, y = RandIndex)) + 
  geom_line(aes(color = Model, linetype = Model), size=0.8) + 
  scale_color_manual(values = c("blue", "red"))+
  theme_classic() + ylim(0,1) + ylab("Rand Index")

pheatmap(RI_matrix, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA)

RI_matrix_est = matrix(0, nrow = L, ncol = L)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix_est [l,ll] = rand.index(VIminlayer[l,], VIminlayer[ll,])
  }
}
RI_matrix_est  = t(RI_matrix_est ) + RI_matrix_est 
pheatmap(RI_matrix_est, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA)

RI_matrix_est_kmeans = matrix(0, nrow = L, ncol = L)
for(l in 1:L){
  for(ll in l:L){
    RI_matrix_est_kmeans [l,ll] = rand.index(m[l, ], m[ll, ])
  }
}

RI_matrix_est_kmeans  = t(RI_matrix_est_kmeans ) + RI_matrix_est_kmeans 
pheatmap(RI_matrix_est_kmeans, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F,border_color=NA)
