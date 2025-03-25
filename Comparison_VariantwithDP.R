#Simulation comparisons with DP at first layer #Section S6 of the supplement
rm(list = ls())
#library(progress)#to draw the progress bar

#devtools::install_github("sarawade/mcclust.ext")
#library(mcclust.ext) #to minimize VI
library(salso) #version 0.3.35                  to comupte psm and point estimate
#library(T4cluster) #to compute psm with psm() function
library(fossil) #version 0.4.0                            to compute rand indexes
library(coda) #version 0.19-4.1              to compute the effective sample size

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(ggplot2)  #to plot
library(GGally) #to plot multivariate continous data


library(rstudioapi) # version 0.15.0
#set working directory to Source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Conditional_t-HDP.R")

##SCENARIO N.1 ##########################################################################
#simulate data
set.seed(1)
data_toy_layer1 <- c(rnorm(100, 0, 1), rnorm(100, 4, 1))
data_toy_layer2 <- c(rnorm(100, 4, 1), rnorm(100, 0, 1))

true_layer1 = c(rep(1,100), rep(2,100))
true_layer2 = c(rep(1,100), rep(2,100))
print(rand.index(true_layer1, true_layer2))

true_chr1 = factor(true_layer1, levels = c("1", "2", "Cluster A", "Cluster B"))
true_chr1[true_chr1 == 1] = "Cluster A"
true_chr1[true_chr1 == 2] = "Cluster B"
true_chr1 = factor(true_chr1, levels = c("Cluster A", "Cluster B"))

true_chr2 = factor(true_layer2, levels = c("1", "2", "Cluster C", "Cluster D"))
true_chr2[true_chr2 == 1] = "Cluster C"
true_chr2[true_chr2 == 2] = "Cluster D"
true_chr2 = factor(true_chr2, levels = c("Cluster C", "Cluster D"))

data_plot = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                       "ClustersL1" = true_chr1,
                       "ClustersL2" = true_chr2)
adj.rand.index(true_layer1,true_layer2)
#plot truth
pmain <- ggplot(data_plot, aes(x = Layer1, y = Layer2, 
                               color = ClustersL1, shape = ClustersL2))+
  geom_point(size = 3 )+theme_minimal(base_size = 18)+
  xlab("Layer n.1") + ylab("Layer n.2") +
  ggpubr::color_palette("jco")  +
  labs(colour="True \nclusters at \nlayer n.1", 
       shape="True \nclusters at \nlayer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot, aes(x = Layer1, fill = ClustersL1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot, aes(x = Layer2, fill = ClustersL2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

data1 = as.matrix(data_toy_layer1 - mean(data_toy_layer1))
data2 = as.matrix(data_toy_layer2 - mean(data_toy_layer2))
data = matrix(c( data_toy_layer1 - mean(data_toy_layer1),
                 data_toy_layer2 - mean(data_toy_layer2)), ncol = 2)

#estimate telescopic - conditional ##################################
totiter = 100000; burnin = 50000
set.seed(0)
m_saved =  telescopic_HDP_NN_uni_VariantwithDP(data, totiter = totiter)
write.table(m_saved[,1,], file = paste("layer_1_cond_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,2,], file = paste("layer_2_cond_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)

psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer1 = salso(m_saved[(totiter-burnin+1):totiter,1,], loss = VI())
VIminlayer2 = salso(m_saved[(totiter-burnin+1):totiter,2,], loss = VI())
write.table(VIminlayer1, file = paste("layer_1_cond_pointEst_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2_cond_pointEst_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
#take care of label switching to detect mistakes (if needed):
#temp = VIminlayer2
#VIminlayer2[temp==1] = 2
#VIminlayer2[temp==2] = 1
mistakes = (VIminlayer1!=true_layer1 | VIminlayer2!=true_layer2)

est_chr1 = factor(VIminlayer1, levels = c("1", "2", "Cluster A", "Cluster B"))
est_chr1[est_chr1 == 1] = "Cluster A"
est_chr1[est_chr1 == 2] = "Cluster B"
est_chr1 = factor(est_chr1, levels = c("Cluster A", "Cluster B"))

est_chr2 = factor(VIminlayer2, levels = c("1", "2", "Cluster C", "Cluster D"))
est_chr2[est_chr2 == 1] = "Cluster C"
est_chr2[est_chr2 == 2] = "Cluster D"
est_chr2 = factor(est_chr2, levels = c("Cluster C", "Cluster D"))

data_plot_tHDP_cond = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                                 "Cluster1" = est_chr1,
                                 "Cluster2" = est_chr2,
                                 "mistakes" = mistakes)
#plot
pmain <- ggplot(data_plot_tHDP_cond, aes(x = Layer1, y = Layer2,color = Cluster1, shape = Cluster2))+
  geom_point(size = 3 )+ geom_point(data = data_plot_tHDP_cond %>% 
                                      filter(mistakes == TRUE),
                                    pch=21, 
                                    size=7, stroke = 2,
                                    colour="red")+
  theme_minimal(base_size = 18)+
  xlab("Layer n.1") + ylab("Layer n.2") +
  ggpubr::color_palette("jco")  +
  labs(colour="Est. \nclusters at \nlayer n.1", 
       shape="Est. \nclusters at \nlayer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot_tHDP_cond, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_tHDP_cond, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

#effective sample size - layer 1
ri= NULL
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,1,], true_layer1))
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 1 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)

effectiveSize(ri)/totiter 
effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) 

#effective sample size - layer 2
ri= NULL
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,2,], true_layer2))
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 2 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)

effectiveSize(ri)/totiter 
effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) 

##SCENARIO N.2 ##########################################################################
#simulate data
set.seed(1)
data_toy_layer1 <- c(rnorm(100, 0, 1), rnorm(100, 4, 1))
data_toy_layer2 <- c(rnorm(50, 4, 1), rnorm(50, 0, 1),rnorm(50, 4 , 1), rnorm(50, 0, 1))

true_layer1 = c(rep(1,100), rep(2,100))
true_layer2 = c(rep(1,50), rep(2,50),rep(1,50), rep(2,50))
print(rand.index(true_layer1, true_layer2))

true_chr1 = factor(true_layer1, levels = c("1", "2", "Cluster A", "Cluster B"))
true_chr1[true_chr1 == 1] = "Cluster A"
true_chr1[true_chr1 == 2] = "Cluster B"
true_chr1 = factor(true_chr1, levels = c("Cluster A", "Cluster B"))

true_chr2 = factor(true_layer2, levels = c("1", "2", "Cluster C", "Cluster D"))
true_chr2[true_chr2 == 1] = "Cluster C"
true_chr2[true_chr2 == 2] = "Cluster D"
true_chr2 = factor(true_chr2, levels = c("Cluster C", "Cluster D"))

data_plot = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                       "Cluster1" = true_chr1,
                       "Cluster2" = true_chr2)

print(adj.rand.index(true_layer1, true_layer2))


#plot truth
pmain <- ggplot(data_plot, aes(x = Layer1, y = Layer2, color = Cluster1, shape = Cluster2))+
  geom_point(size = 3 )+theme_minimal(base_size = 18)+
  xlab("Layer n.1") + ylab("Layer n.2") +
  ggpubr::color_palette("default") +
  labs(colour="True \nclusters at \nlayer n.1", 
       shape="True \nclusters at \nlayer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("default")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("default")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

data1 = as.matrix(data_toy_layer1 - mean(data_toy_layer1))
data2 = as.matrix(data_toy_layer2 - mean(data_toy_layer2))
data = matrix(c( data_toy_layer1 - mean(data_toy_layer1),
                 data_toy_layer2 - mean(data_toy_layer2)), ncol = 2)

#estimate telescopic - conditional ##################################
totiter = 100000; burnin = 50000
set.seed(0)
m_saved =  telescopic_HDP_NN_uni_VariantwithDP(data, totiter = totiter)
write.table(m_saved[,1,], file = paste("layer_1_cond_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,2,], file = paste("layer_2_cond_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)

psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer1 = salso(m_saved[(totiter-burnin+1):totiter,1,], loss = VI())
VIminlayer2 = salso(m_saved[(totiter-burnin+1):totiter,2,], loss = VI())
write.table(VIminlayer1, file = paste("layer_1_cond_pointEst_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2_cond_pointEst_DP.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
#take care of label switching to detect mistakes (if needed):
#temp = VIminlayer1
#VIminlayer1[temp==1] = 2
#VIminlayer1[temp==2] = 1
mistakes = (VIminlayer1!=true_layer1 | VIminlayer2!=true_layer2)

est_chr1 = factor(VIminlayer1, levels = c("1", "2", "Cluster A", "Cluster B"))
est_chr1[est_chr1 == 1] = "Cluster A"
est_chr1[est_chr1 == 2] = "Cluster B"
est_chr1 = factor(est_chr1, levels = c("Cluster A", "Cluster B"))

est_chr2 = factor(VIminlayer2, levels = c("1", "2", "Cluster C", "Cluster D"))
est_chr2[est_chr2 == 1] = "Cluster C"
est_chr2[est_chr2 == 2] = "Cluster D"
est_chr2 = factor(est_chr2, levels = c("Cluster C", "Cluster D"))

data_plot_tHDP_cond = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                                 "Cluster1" = est_chr1,
                                 "Cluster2" = est_chr2,
                                 "mistakes" = mistakes)
#plot
pmain <- ggplot(data_plot_tHDP_cond, aes(x = Layer1, y = Layer2,color = Cluster1, shape = Cluster2))+
  geom_point(size = 3 )+ geom_point(data = data_plot_tHDP_cond %>% 
                                      filter(mistakes == TRUE),
                                    pch=21, 
                                    size=7, stroke = 2,
                                    colour="red")+
  theme_minimal(base_size = 18)+
  xlab("Layer n.1") + ylab("Layer n.2") +
  ggpubr::color_palette("default")  +
  labs(colour="Est.\nclusters at \nlayer n.1", 
       shape="Est. \nclusters at \nlayer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot_tHDP_cond, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("default")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_tHDP_cond, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("default")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

#effective sample size - layer 1
ri= NULL
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,1,], true_layer1))
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 1 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)

effectiveSize(ri)/totiter 
effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) 

#effective sample size - layer 2
ri= NULL
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,2,], true_layer2))
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 2 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)

effectiveSize(ri)/totiter 
effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) 


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

#estimate telescopic - conditional ##################################
totiter = 100000; burnin = 50000
set.seed(2)
m_saved =  telescopic_HDP_NNIG_uni_VariantwithDP(data, totiter = totiter)

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




