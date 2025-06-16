# The following code simulate data (scenario 1) 
# and estimate the Enriched Dirichlet process
# (note that estimate of k-means, the t-HDP, the LSBP and the DP-t-HDP 
# are obtained running separate
# codes, named specifically "Simulation_n1.R" and 
# "Comparison_VariantwithDP.R")
# results obtained from this code are presented in Table 1 of Section 7.1
# and Section S5 and S6 of the Supplement of the paper.
# Table 1 of main and Table S6.1 of Supplement *E-DP's columns

rm(list = ls())

library(salso)  # version 0.3.35          to comupte psm and point estimate
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

library(rstudioapi) # version 0.15.0
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source mcmc codes for telescopic clustering models:
source("Conditional_t-HDP.R")

#few useful functions used in the MCMC: 
#has to work for one X and many param
kernel_eval = function(x, hyper1, hyper2, log = TRUE){
  return(dnorm(x, hyper1, sqrt(hyper2), log = log))
}

posterior_sampler = function(x, k0 = 0.1, a = 0.1, b = 0.1){
  n_temp = length(x)
  if(n_temp>0){
    xbar = mean(x)
    sum_squared = sum((x - xbar)^2)
    sigma = 1 / rgamma(1, a + n_temp / 2, b + sum_squared/2 +
                         xbar^2/2 * n*k0 / (k0+n_temp))
    mu = rnorm(1, mean = n_temp / (n_temp + k0) * xbar,
               sd = sqrt(sigma / (n_temp+k0)))
    return( list(mu = mu, sigma = sigma) )
  }else{
    sigma = 1 / rgamma(1, a , b )
    mu = rnorm(1, mean = 0,
               sd = sqrt(sigma / k0))
    return( list(mu = mu, sigma = sigma) )
  }
}

#if true the code performe 100k iterations of the mcmc, 
#if false the code upload the outcome of the mcmc from a csv
run_MCMC = FALSE 

#1. simulate data ##############################################################
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
  #print(adj.rand.index(true_layer[,l], true_layer[,l-1]))
}

data = scale(data, scale = FALSE)
#ggpairs(as.data.frame(data))+ theme_bw()
# data_plot = as.data.frame(data[,seq(1,10,2)])
# colnames(data_plot) = c("Layer 1", "Layer 3", "Layer 5", "Layer 7","Layer 9")
# ggpairs(data_plot, upper = list(continuous = "density"), aes(color = as.factor(true_layer[,1]), 
#                                                              shape = as.factor(true_layer[,1])))+ theme_classic()


if(run_MCMC){
  #the data in a matrix X, which is nxL
  X = data 
  n_tot = dim(X)[1]
  
  H = 2 #number of components
  for(l in 1:L){
    assign(paste0("prob", l), rep(1,H^(l-1)))
    assign(paste0("atomtheta", l), rnorm(H^(l-1)))
    assign(paste0("atomsigma", l), rgamma(H^(l-1), 0.1, 0.1))
  }
  
  
  m = matrix(0, nrow = L, ncol = n_tot) #clustering configuration
  
  m[1, ] = stats::kmeans(as.matrix(X[,1]),min(2,H))$cluster
  for(l in 2:L){
    temp = 0
    for (mm in unique(m[l-1,])){
      n_m = sum(m[l-1,]==mm) 
      if(n_m>floor(n_m/2)+1){
        m[l, m[l-1,]==mm ] = temp + 
                  stats::kmeans(as.matrix(X[m[l-1,]==mm ,l]),
                                min(2, floor(n_m/2)+1, H))$cluster
      }else{
        m[l, m[l-1,]==mm ] = max(m[l, ]) + seq(1,n_m)
      }
      temp = temp + H
    }
  }
  
  for(l in 1:L){
    assign(paste0("nbar", l), rep(0,H^(l-1)))
  }
  
  set.seed(0)
  totiter = 10000; burnin = 2000
  alpha = rep(1, L)
  k0 = 0.1; a_sigma = 0.1; b_sigma = 0.1 
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = totiter, clear = FALSE, width= 100)
  
  for (iter in 1:totiter){
    pb$tick()
    for (l in 1:L){
      print(c(iter,l))
      #update cluster frequencies
      assign((paste0("nbar", l)), table(factor(m[l,], levels = 1:(H^(l)))) )
      howmanyDP = H^(l-1); init = 0 
      probtemp = NULL; thetatemp = NULL; sigmatemp = NULL
      for (dp in 1:howmanyDP){
        probtemp2 = NULL
        nbar = get(paste0("nbar", l))
        nl_sum = c(rev(cumsum(rev(nbar[(init+1):(init+H)])))[2:H],0)
        b = rep(0, H); 
          for (hh in 1:H){
            b[hh] = rbeta(1, 1 + nbar[init + hh], alpha[l] + nl_sum[hh])
            if(hh>1){
              probtemp2 = c(probtemp2, log(b[hh]) + sum(log(1 - b[1:(hh-1)])))
            }else{
              probtemp2 = c(probtemp2, log(b[hh]) )
            }
            out = posterior_sampler(X[m[l,]==(init + hh),l], k0, a_sigma, b_sigma)
            thetatemp = c(thetatemp, out$mu)
            sigmatemp = c(sigmatemp, out$sigma)
          }
        probtemp2 = exp(probtemp2) / sum(exp(probtemp2))
        probtemp = c(probtemp, probtemp2)
        init = init + H
      }
      assign((paste0("prob", l)), probtemp )
      assign((paste0("atomtheta", l)), thetatemp )
      assign((paste0("atomsigma", l)), sigmatemp )
    }
  
      
    for (l in 1:L){
      print(c(iter,l))
      for (i in 1:n_tot){
        if(l>1){past = m[l-1, i]}else{past = 1}
          fut = rep(0,H)
        if(l<L){
          for (hh in 1:H){
            present = (past-1)*H + hh; tempfut = 0
            for(ll in (l+1):L){
              tempthe = get(paste0("atomtheta", 
                              ll))[(H^(ll-l)*(present-1)+1):(H^(ll-l)*present)]
              tempsig = get(paste0("atomsigma", 
                              ll))[(H^(ll-l)*(present-1)+1):(H^(ll-l)*present)]
              tempprob = get(paste0("prob", 
                                ll))[(H^(ll-l)*(present-1)+1):(H^(ll-l)*present)]
              tempfut = tempfut + log(sum(kernel_eval(X[i,ll], tempthe,
                                                 tempsig, log = FALSE)*
                                       tempprob))
            }
          fut[hh] = tempfut
          }
        }
        tempthe = get(paste0("atomtheta", 
                              l))[(H*(past-1)+1):(H*past)]
        tempsig = get(paste0("atomsigma", 
                              l))[(H*(past-1)+1):(H*past)]
        tempprob = get(paste0("prob", 
                      l))[(H*(past-1)+1):(H*past)]
        prob = log(tempprob) + kernel_eval(X[i,l], tempthe, tempsig) + fut
        if(max(prob)==-Inf){prob[]=1}
        prob = prob - max(prob)
        prob = exp(prob) 
        if(sum(exp(fut))!=0){ #if fut is zero it cannot move
          m[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
        }
      }
      print(c(iter,table(m[l,])))
    }
  
    for (layer in (1:L)){
      write.table(t(m[layer,]), file = paste("031_layer_",layer,".csv",sep=""),
                  append = TRUE, row.names = FALSE,
                  col.names= FALSE)
    }
  }
}else{
  totiter = 10000; burnin = 2000; n_tot = 200
  pb <- progress_bar$new(
    format = " MCMC output [:bar] :percent Estimated completion time: :eta",
    total = L, clear = FALSE, width= 100)
  #load results:
  m_saved = array(NA,c(totiter, L, n_tot))
  for(l in 1:L){
    temp = as.matrix(read.table(paste("output/031_layer_",l,".csv", sep="")))
    m_saved[,l,] = temp
    pb$tick()
  }
}

psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
psm3 = psm(m_saved[(totiter-burnin+1):totiter,3,])
psm4 = psm(m_saved[(totiter-burnin+1):totiter,4,])
psm5 = psm(m_saved[(totiter-burnin+1):totiter,5,])
psm6 = psm(m_saved[(totiter-burnin+1):totiter,6,])
psm7 = psm(m_saved[(totiter-burnin+1):totiter,7,])
psm8 = psm(m_saved[(totiter-burnin+1):totiter,8,])
psm9 = psm(m_saved[(totiter-burnin+1):totiter,9,])
psm10 = psm(m_saved[(totiter-burnin+1):totiter,10,])

#compute minVI
VIminlayer1 = salso(m_saved[(totiter-burnin+1):totiter,1,], loss = VI())
VIminlayer2 = salso(m_saved[(totiter-burnin+1):totiter,2,], loss = VI())
VIminlayer3 = salso(m_saved[(totiter-burnin+1):totiter,3,], loss = VI())
VIminlayer4 = salso(m_saved[(totiter-burnin+1):totiter,4,], loss = VI())
VIminlayer5 = salso(m_saved[(totiter-burnin+1):totiter,5,], loss = VI())
VIminlayer6 = salso(m_saved[(totiter-burnin+1):totiter,6,], loss = VI())
VIminlayer7 = salso(m_saved[(totiter-burnin+1):totiter,7,], loss = VI())
VIminlayer8 = salso(m_saved[(totiter-burnin+1):totiter,8,], loss = VI())
VIminlayer9 = salso(m_saved[(totiter-burnin+1):totiter,9,], loss = VI())
VIminlayer10 = salso(m_saved[(totiter-burnin+1):totiter,10,], loss = VI())

# write.table(VIminlayer1, file = paste("layer_1.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer2, file = paste("layer_2.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer3, file = paste("layer_3.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer4, file = paste("layer_4.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer5, file = paste("layer_5.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer6, file = paste("layer_6.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer7, file = paste("layer_7.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer8, file = paste("layer_8.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer9, file = paste("layer_9.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer10, file = paste("layer_10.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)

print("Values fifth column table 1/table S6.1")
print(paste("Rand Index Layer n.", 1 ,": ", rand.index(VIminlayer1,true_layer[,1])))
print(paste("Rand Index Layer n.", 2 ,": ", rand.index(VIminlayer2,true_layer[,2])))
print(paste("Rand Index Layer n.", 3 ,": ", rand.index(VIminlayer3,true_layer[,3])))
print(paste("Rand Index Layer n.", 4 ,": ", rand.index(VIminlayer4,true_layer[,4])))
print(paste("Rand Index Layer n.", 5 ,": ", rand.index(VIminlayer5,true_layer[,5])))
print(paste("Rand Index Layer n.", 6 ,": ", rand.index(VIminlayer6,true_layer[,6])))
print(paste("Rand Index Layer n.", 7 ,": ", rand.index(VIminlayer7,true_layer[,7])))
print(paste("Rand Index Layer n.", 8 ,": ", rand.index(VIminlayer8,true_layer[,8])))
print(paste("Rand Index Layer n.", 9 ,": ", rand.index(VIminlayer9,true_layer[,9])))
print(paste("Rand Index Layer n.", 10 ,": ", rand.index(VIminlayer10,true_layer[,10])))

temp = 0
for(l in 1:L){
  temp = temp + rand.index(get(paste0("VIminlayer",l)),true_layer[,l])
}
temp / L

print("Values last column table 1/ tenth column of table S6.1")
print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3!=true_layer[,3]) ))
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4!=true_layer[,4]) ))
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5!=true_layer[,5]) ))
#table(VIminlayer6); 
#correct label switching
VIminlayer6[VIminlayer6==3]=10; VIminlayer6[VIminlayer6==2]=3
VIminlayer6[VIminlayer6==10]=2;
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6!=true_layer[,6]) ))
#table(VIminlayer7); 
#correct label switching
VIminlayer7[VIminlayer7==3]=10; VIminlayer7[VIminlayer7==2]=3
VIminlayer7[VIminlayer7==10]=2;
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7!=true_layer[,7]) ))
#table(VIminlayer8); 
#correct label switching
VIminlayer8[VIminlayer8==3]=10; VIminlayer8[VIminlayer8==2]=3
VIminlayer8[VIminlayer8==10]=2;
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8!=true_layer[,8]) )) 
#table(VIminlayer9); 
#correct label switching
VIminlayer9[VIminlayer9==4]=10; VIminlayer9[VIminlayer9==2]=4
VIminlayer9[VIminlayer9==10]=2;
VIminlayer9[VIminlayer9==4]=10; VIminlayer9[VIminlayer9==1]=4
VIminlayer9[VIminlayer9==10]=1;
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9!=true_layer[,9]) )) 
#table(VIminlayer10); 
#correct label switching
VIminlayer10[VIminlayer10==4]=10; VIminlayer10[VIminlayer10==2]=4
VIminlayer10[VIminlayer10==10]=2;
VIminlayer10[VIminlayer10==4]=10; VIminlayer10[VIminlayer10==1]=4
VIminlayer10[VIminlayer10==10]=1;
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10!=true_layer[,10]) ))

