#Simulation scenario 4 with 3 layers (true clustering)

library(mcclust.ext) #to minimize VI
library(T4cluster) #to compute psm with psm() function
library(fossil) #to compute rand indexes
library(coda) #to compute the effective sample size

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(ggplot2)
library(GGally) #to plot multivariate continous data

source("Conditional_t-HDP.R")

##SCENARIO N.4 ##########################################################################
#simulate data
set.seed(1)
n = 200
L = 3 #bivariate
true_layer = matrix(NA, nrow = n, ncol = L*2)
data = true_layer
true_layer[,1] = c(rep(1,floor(n/2)), rep(2,n - floor(n/2)))
data[true_layer[,1]==1, 1] = rnorm(sum(true_layer[,1]==1), 0, 1)
data[true_layer[,1]==1, 2] = rnorm(sum(true_layer[,1]==1), 0, 1)
data[true_layer[,1]==2, 1] = rnorm(sum(true_layer[,1]==2), 4, 1)
data[true_layer[,1]==2, 2] = rnorm(sum(true_layer[,1]==2), 4, 1)

l = 2
true_layer[,l] = true_layer[,l-1]
moving = sample(1:n, floor(n*0.05))
true_layer[moving,l] = true_layer[moving,l]%%2 + 1
data[true_layer[,l]==1, 3] = rnorm(sum(true_layer[,l]==1), 0, 1)
data[true_layer[,l]==1, 4] = rnorm(sum(true_layer[,l]==1), 0, 1)
data[true_layer[,l]==2, 3] = rnorm(sum(true_layer[,l]==2), 4, 1)
data[true_layer[,l]==2, 4] = rnorm(sum(true_layer[,l]==2), 4, 1)
print(rand.index(true_layer[,l], true_layer[,l-1]))

l = 3
true_layer[,l] = true_layer[,l-2]
moving = sample(1:n, floor(n*0.05))
true_layer[moving,l] = true_layer[moving,l]%%2 + 1
data[true_layer[,l]==1, 5] = rnorm(sum(true_layer[,l]==1), 0, 1)
data[true_layer[,l]==1, 6] = rnorm(sum(true_layer[,l]==1), 0, 1)
data[true_layer[,l]==2, 5] = rnorm(sum(true_layer[,l]==2), 4, 1)
data[true_layer[,l]==2, 6] = rnorm(sum(true_layer[,l]==2), 4, 1)
print(rand.index(true_layer[,l], true_layer[,l-2]))


data = scale(data)
#ggpairs(as.data.frame(data))+ theme_bw()
#ggpairs(as.data.frame(data[,1:4]))+ theme_bw()


source("Conditional_t-HDP.R")
#estimate telescopic - conditional ##################################
colayer = c(1,1,2,2,3,3)

totiter = 100000
set.seed(1)
m_saved = telescopic_HDP_NormWish_cond_3L1P(data, colayer, totiter = totiter)
#telescopic_HDP_normnorm_cond_3L1P(data, colayer, totiter = totiter)
