#Simulation scenario 3 with 10 layers (true clustering)

#devtools::install_github(“sarawade/mcclust.ext”)
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

source("Conditional_t-HDP.R")

##SCENARIO N.3 ##########################################################################
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
ggpairs(data_plot, upper = list(continuous = "density"), aes(color = as.factor(true_layer[,1]), 
                                                                             shape = as.factor(true_layer[,1])))+ theme_classic()

#k-means 
#fviz_nbclust(as.matrix(data[,l]), kmeans, method = "gap", k.max = 5)
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
  if(sum(m[l,]%in%c(1,2))!=200){print(paste("label switching problem at layer", l))}
  mist = sum(m[l, ]!=true_layer[,l])
  print(paste("Mistakes Layer n.", l ,": ", 
              min(mist, n - mist )))
}

#estimate telescopic - conditional ##################################
totiter = 100000; burnin = 50000
set.seed(1)
m_saved =  telescopic_HDP_NNIG_uni(data, totiter = totiter)

write.table(m_saved[,1,], file = paste("layer_1_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,2,], file = paste("layer_2_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,3,], file = paste("layer_3_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,4,], file = paste("layer_4_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,5,], file = paste("layer_5_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,6,], file = paste("layer_6_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,7,], file = paste("layer_7_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,8,], file = paste("layer_8_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,9,], file = paste("layer_9_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,10,], file = paste("layer_10_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
#reload results:
#m_saved = array(NA,c(totiter, L, n_tot))
#for(l in 1:L){
#  temp = as.matrix(read.table(paste("layer_",l,"_cond.csv", sep="")))
#  m_saved[,l,] = temp
#}
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

write.table(VIminlayer1, file = paste("layer_1_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer3, file = paste("layer_3_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer4, file = paste("layer_4_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer5, file = paste("layer_5_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer6, file = paste("layer_6_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer7, file = paste("layer_7_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer8, file = paste("layer_8_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer9, file = paste("layer_9_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer10, file = paste("layer_10_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)

#effective sample size and mixing
par(mfrow=c(2,1))
for (l in 1:L){
  ri = NULL
  for (iter in 1:totiter){
    ri = c(ri, rand.index(m_saved[iter,l,], true_layer[,l]))
  }
  plot(ri, type = "l", xlab="MCMC iterations")
  abline(v = burnin, col="red", lwd=3, lty=2)
  title(main=paste("Rand Index Layer ",l," - Trace plot \n Conditional algorithm"),
                 cex.lab=0.75)
  if(l%in%c(1,5,10)){print(effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) )}
}

par(mfrow=c(3,1))
for (l in c(1,5,10)){
  ri = NULL
  for (iter in (burnin+1):totiter){
    ri = c(ri, rand.index(m_saved[iter,l,], true_layer[,l]))
  }
  hist(ri, xlab=paste(" "), ylab=paste(" "), xlim=c(0,1), col=l,
       main=paste("Posterior distribution Rand Index Layer n.",l),
       cex.lab=0.75)
}



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


print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3!=true_layer[,3]) ))
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4!=true_layer[,4]) ))
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5!=true_layer[,5]) ))
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6!=true_layer[,6]) ))
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7!=true_layer[,7]) ))
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8!=true_layer[,8]) )) #label switching
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9!=true_layer[,9]) )) #label switching
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10!=true_layer[,10]) )) #label switching




#logit stick breaking 
#devtools::install_github("tommasorigon/LSBP")
library(LSBP)

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


set.seed(0)
R <- 10000 #iterations
burnin <- 2000
fit_gibbs <- LSBP_Gibbs(Y ~ X|X, data=data_flat, H=5, 
                        control=control_Gibbs(R=R,burn_in=burnin,method_init="random") )

psm1 = psm(fit_gibbs$param$G[,1:200 ])
psm2 = psm(fit_gibbs$param$G[,201:400 ])
psm3 = psm(fit_gibbs$param$G[,401:600 ])
psm4 = psm(fit_gibbs$param$G[,601:800 ])
psm5 = psm(fit_gibbs$param$G[,801:1000 ])
psm6 = psm(fit_gibbs$param$G[,1001:1200 ])
psm7 = psm(fit_gibbs$param$G[,1201:1400 ])
psm8 = psm(fit_gibbs$param$G[,1401:1600 ])
psm9 = psm(fit_gibbs$param$G[,1601:1800 ])
psm10 = psm(fit_gibbs$param$G[,1801:2000 ])

VIminlayer1_LSBP = minVI(psm1, fit_gibbs$param$G[,1:200], method = "all")$cl[1,]
VIminlayer2_LSBP = minVI(psm2, fit_gibbs$param$G[,201:400], method = "all")$cl[1,]
VIminlayer3_LSBP = minVI(psm3, fit_gibbs$param$G[,401:600], method = "all")$cl[1,]
VIminlayer4_LSBP = minVI(psm4, fit_gibbs$param$G[,601:800], method = "all")$cl[1,]
VIminlayer5_LSBP = minVI(psm5, fit_gibbs$param$G[,801:1000], method = "all")$cl[1,]
VIminlayer6_LSBP = minVI(psm6, fit_gibbs$param$G[,1001:1200], method = "all")$cl[1,]
VIminlayer7_LSBP = minVI(psm7, fit_gibbs$param$G[,1201:1400], method = "all")$cl[1,]
VIminlayer8_LSBP = minVI(psm8, fit_gibbs$param$G[,1401:1600], method = "all")$cl[1,]
VIminlayer9_LSBP = minVI(psm9, fit_gibbs$param$G[,1601:1800], method = "all")$cl[1,]
VIminlayer10_LSBP = minVI(psm10, fit_gibbs$param$G[,1801:2000], method = "all")$cl[1,]

print(paste("Rand Index Layer n.", 1 ,": ", rand.index(VIminlayer1_LSBP,true_layer[,1])))
print(paste("Rand Index Layer n.", 2 ,": ", rand.index(VIminlayer2_LSBP,true_layer[,2])))
print(paste("Rand Index Layer n.", 3 ,": ", rand.index(VIminlayer3_LSBP,true_layer[,3])))
print(paste("Rand Index Layer n.", 4 ,": ", rand.index(VIminlayer4_LSBP,true_layer[,4])))
print(paste("Rand Index Layer n.", 5 ,": ", rand.index(VIminlayer5_LSBP,true_layer[,5])))
print(paste("Rand Index Layer n.", 6 ,": ", rand.index(VIminlayer6_LSBP,true_layer[,6])))
print(paste("Rand Index Layer n.", 7 ,": ", rand.index(VIminlayer7_LSBP,true_layer[,7])))
print(paste("Rand Index Layer n.", 8 ,": ", rand.index(VIminlayer8_LSBP,true_layer[,8])))
print(paste("Rand Index Layer n.", 9 ,": ", rand.index(VIminlayer9_LSBP,true_layer[,9])))
print(paste("Rand Index Layer n.", 10 ,": ", rand.index(VIminlayer10_LSBP,true_layer[,10])))

print(paste("Mistakes Layer n.", 1 ,": ", sum(VIminlayer1_LSBP!=true_layer[,1]) ))
print(paste("Mistakes Layer n.", 2 ,": ", sum(VIminlayer2_LSBP!=true_layer[,2]) ))
print(paste("Mistakes Layer n.", 3 ,": ", sum(VIminlayer3_LSBP!=true_layer[,3]) ))
print(paste("Mistakes Layer n.", 4 ,": ", sum(VIminlayer4_LSBP!=true_layer[,4]) ))#label switching
print(paste("Mistakes Layer n.", 5 ,": ", sum(VIminlayer5_LSBP!=true_layer[,5]) ))
print(paste("Mistakes Layer n.", 6 ,": ", sum(VIminlayer6_LSBP!=true_layer[,6]) ))
print(paste("Mistakes Layer n.", 7 ,": ", sum(VIminlayer7_LSBP!=true_layer[,7]) ))
#VIminlayer8_LSBP[VIminlayer8_LSBP==5]=2
print(paste("Mistakes Layer n.", 8 ,": ", sum(VIminlayer8_LSBP!=true_layer[,8]) )) #label switching
print(paste("Mistakes Layer n.", 9 ,": ", sum(VIminlayer9_LSBP!=true_layer[,9]) )) 
print(paste("Mistakes Layer n.", 10 ,": ", sum(VIminlayer10_LSBP!=true_layer[,10]) )) #label switching
