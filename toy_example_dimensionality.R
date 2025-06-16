#example of dimensionality influence 
# REPRODUCE FIGURE 1
rm(list = ls())
library(tidyverse)   # data manipulation # version 2.0.0
library(cluster)     # clustering algorithms # version 2.1.4
library(factoextra)  # clustering algorithms & visualization # version 1.0.7       
library(mcclust)     # for point estimate DP # version 1.0.1
library(mcclust.ext) # for point estimate DP # version 1.0 
library(progress)   # version 1.2.3    to draw the progress bar

set.seed(1)
data_toy_X1 <- c(rnorm(100, 0, 1), rnorm(100, 2, 1))
data_toy_Y1 <- c(rnorm(50, 2, 1), rnorm(50, 0, 1), rnorm(50, 2 , 1), rnorm(50, 0, 1))
data_toy_Y2 <- c(rnorm(50, 2, 1), rnorm(50, 0, 1), rnorm(50, 2 , 1), rnorm(50, 0, 1))
data_toy_Y3 <- c(rnorm(50, 2, 1), rnorm(50, 0, 1), rnorm(50, 2 , 1), rnorm(50, 0, 1))
#true clusters are:
#X: 1:100(0) 2:100(4)
#Y: 1:50(4) 2:50(0) 1:50(4) 2:50(0)
#X+Y: 1:50(0,4) 2:50(0,0) 3:50(4,4) 4:50(4,0)

data = data.frame(matrix(c( 
  data_toy_X1 - mean(data_toy_X1),
  data_toy_Y1 - mean(data_toy_Y1),
  data_toy_Y2 - mean(data_toy_Y2),
  data_toy_Y3 - mean(data_toy_Y3)), ncol = 4))


#k-means
set.seed(1)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(data, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:6

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(data, kmeans, method = "silhouette")

fviz_nbclust(data, kmeans, method = "gap", k.max = 8)

kmeans_clusters = kmeans(data, 2)

data.pca <- prcomp(data[,2:4])
data["Dim1_PCA_feature2"] = get_pca_ind(data.pca)$coord[,1]
print(fviz_cluster(kmeans_clusters, data = data, 
             choose.vars = c("X1", "Dim1_PCA_feature2"), geom="point") + theme_classic()+
  labs(x="Feature 1", y = "First PC of feature 2") + ggtitle(" ") + 
  theme(text = element_text(size = 20)))
print("Figure 1 - panel B just printed")

#Dirichlet process
#compute the marginal likelihood of a norm-norm model (marg lik of a cluster) multivariate
margnn_multi <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  J = dim(y)[2]; n = dim(y)[1]
  sigma = sqrt(sigma2); s = sqrt(s2)
  p = 0
  for (jj in (1:J)){
    Y2 = sum(y[,jj]**2); Y = sum(y[,jj])
    p = p - n / 2 * log(2 * pi) - (n - 1) * log(sigma) - 1/2 * log(n * s2 + sigma2) - 
      Y2 / (2 * sigma2) - m**2 / (2 * s2) + (s * Y / sigma + sigma * m / s ) ** 2 / (2 * (n * s2 + sigma2))
  }
  return(p)
} 
mcmc_DPM_norm_norm_multi<- function(y, hyper=c(0.1, 0, 1, 1), c = FALSE, totiter = 20000){
  pb <- progress_bar$new(
    format = " MCMC output [:bar] :percent Estimated completion time: :eta",
    total = totiter, clear = FALSE, width= 100)
  n = dim(y)[1]
  alpha = hyper[1]
  m = hyper[2]
  s2 = hyper[3]
  sigma2 = hyper[4]
  if(!c[1]){c = kmeans(y, centers = floor(n/10), nstart = 20)$cluster}
  print(paste("initialization to", length(unique(c)), "clusters"))
  temp = 0; c_saved = NULL
  for (iter in (1:totiter)){
    for (i in (1:n)){
      nc = table(c)
      nc[as.character(c[i])] = nc[as.character(c[i])] - 1; c[i] = NA; c_unique = unique(c[!is.na(c)])
      prob = NULL; j = 1
      for (cc in c_unique){
        prob[j] = margnn_multi(rbind(y[t(c == cc & !is.na(c)),,drop = FALSE], y[i,,drop = FALSE]), m, s2, sigma2, log = TRUE ) - 
          margnn_multi(y[t(c == cc & !is.na(c)),,drop = FALSE], m, s2, sigma2, log = TRUE ) +
          log(nc[as.character(cc)])
        j = j + 1
      } 
      prob[j] = margnn_multi(y[i,,drop = FALSE], m, s2, sigma2, log = TRUE) + log(alpha)
      prob = prob - max(prob)
      prob = exp(prob) / sum( exp(prob) )
      c[i] = sample(c(c_unique, setdiff(1:n,c_unique)[1]), 1, prob = prob)
    }
    #print(c(iter,length(unique(c))))
    c_saved = rbind(c_saved,c)
    pb$tick()
  }
  return(c_saved)
}


data = as.matrix(data[,1:4])

set.seed(1)
DP_clusters_chain = mcmc_DPM_norm_norm_multi(data, totiter = 1000)
psm = comp.psm(DP_clusters_chain)
DP_clusters = minVI(psm, DP_clusters_chain, method= "all")$cl[1,]
#label_switiching
#DP_clusters[DP_clusters==1] = 0 
#DP_clusters[DP_clusters==2] = 1 
#DP_clusters[DP_clusters==0] = 2 

data = data.frame(data)
data.pca <- prcomp(data[,2:4])
data["Dim1_PCA_feature2"] = get_pca_ind(data.pca)$coord[,1]
DP_clusters = list(data = data, cluster = DP_clusters)

print(fviz_cluster(DP_clusters, data = data, 
             choose.vars = c("X1", "Dim1_PCA_feature2"), geom="point") + theme_classic()+
  labs(x="Feature 1", y = "First PC of feature 2") + ggtitle(" ") + 
  theme(text = element_text(size = 20)))
print("Figure 1 - panel C just printed")


#true clusters 
true_clusters = c(rep(1,50),rep(2,50),rep(3,50),rep(4,50))
true_clusters = list(data = data, cluster = true_clusters)

print(fviz_cluster(true_clusters, data = data, 
             choose.vars = c("X1", "Dim1_PCA_feature2"), geom="point") + theme_classic()+
  labs(x="Feature 1", y = "First PC of feature 2") + ggtitle(" ") + 
  theme(text = element_text(size = 20)))
print("Figure 1 - panel A just printed")