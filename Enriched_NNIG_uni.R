rm(list=ls())
#Simulation scenario 3 with 10 layers (Enriched)
library(progress)#to draw the progress bar
library(mcclust.ext) #to minimize VI
library(T4cluster) #to compute psm with psm() function
library(fossil) #to compute rand indexes
library(coda) #to compute the effective sample size
library(factoextra)#for gap stat

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(ggplot2)
library(GGally) #to plot multivariate continous data

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
  print(adj.rand.index(true_layer[,l], true_layer[,l-1]))
}

data = scale(data, scale = FALSE)
#ggpairs(as.data.frame(data))+ theme_bw()
data_plot = as.data.frame(data[,seq(1,10,2)])
colnames(data_plot) = c("Layer 1", "Layer 3", "Layer 5", "Layer 7","Layer 9")
ggpairs(data_plot, upper = list(continuous = "density"), aes(color = as.factor(true_layer[,1]), 
                                                             shape = as.factor(true_layer[,1])))+ theme_classic()

H = 2 #number of components

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

#the data in a matrix X, which is nxL
X = data 
n_tot = dim(X)[1]

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
    write.table(t(m[layer,]), file = paste("Enriched",layer,".csv",sep=""), append = TRUE, row.names = FALSE,
                col.names= FALSE)
  }
}


#load results:
m_saved = array(NA,c(totiter, L, n_tot))
for(l in 1:L){
  temp = as.matrix(read.table(paste("Enriched",l,".csv", sep="")))
  m_saved[,l,] = temp
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

VIminlayer1 = minVI(psm1, m_saved[(totiter-burnin+1):totiter,1,], method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, m_saved[(totiter-burnin+1):totiter,2,], method = "all")$cl[1,]
VIminlayer3 = minVI(psm3, m_saved[(totiter-burnin+1):totiter,3,], method = "all")$cl[1,]
VIminlayer4 = minVI(psm4, m_saved[(totiter-burnin+1):totiter,4,], method = "all")$cl[1,]
VIminlayer5 = minVI(psm5, m_saved[(totiter-burnin+1):totiter,5,], method = "all")$cl[1,]
VIminlayer6 = minVI(psm6, m_saved[(totiter-burnin+1):totiter,6,], method = "all")$cl[1,]
VIminlayer7 = minVI(psm7, m_saved[(totiter-burnin+1):totiter,7,], method = "all")$cl[1,]
VIminlayer8 = minVI(psm8, m_saved[(totiter-burnin+1):totiter,8,], method = "all")$cl[1,]
VIminlayer9 = minVI(psm9, m_saved[(totiter-burnin+1):totiter,9,], method = "all")$cl[1,]
VIminlayer10 = minVI(psm10, m_saved[(totiter-burnin+1):totiter,10,], method = "all")$cl[1,]

write.table(VIminlayer1, file = paste("layer_1.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer3, file = paste("layer_3.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer4, file = paste("layer_4.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer5, file = paste("layer_5.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer6, file = paste("layer_6.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer7, file = paste("layer_7.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer8, file = paste("layer_8.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer9, file = paste("layer_9.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer10, file = paste("layer_10.csv",sep=""), row.names = FALSE,
            col.names= FALSE)

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

print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3!=true_layer[,3]) ))
table(VIminlayer4); 
VIminlayer4[VIminlayer4==9]=1;VIminlayer4[VIminlayer4==14]=2
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4!=true_layer[,4]) ))
table(VIminlayer5); 
VIminlayer5[VIminlayer5==18]=1;VIminlayer5[VIminlayer5==28]=2
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5!=true_layer[,5]) ))
table(VIminlayer6); 
VIminlayer6[VIminlayer6==36]=1;VIminlayer6[VIminlayer6==56]=2
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6!=true_layer[,6]) ))
table(VIminlayer7); 
VIminlayer7[VIminlayer7==71]=1;VIminlayer7[VIminlayer7==112]=2
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7!=true_layer[,7]) ))
table(VIminlayer8); 
VIminlayer8[VIminlayer8==142]=1;VIminlayer8[VIminlayer8==224]=2
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8!=true_layer[,8]) )) 
table(VIminlayer9); 
VIminlayer9[VIminlayer9==283]=1;VIminlayer9[VIminlayer9==447]=2
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9!=true_layer[,9]) )) #label switching
table(VIminlayer10); VIminlayer10 = VIminlayer10+100
VIminlayer10[VIminlayer10==102]=1;VIminlayer10[VIminlayer10==104]=2
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10!=true_layer[,10]) )) #label switching

