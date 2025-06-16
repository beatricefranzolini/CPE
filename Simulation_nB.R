# Simulation scenario B with two layers 
# The following code simulate data (scenario B) 
# and estimate the t-HDP, the independent and the constant clustering model
# results obtained from this code are presented in Section S5 of the Supplement.
# Figure S5.1(panel B), Figure S5.3 and Figure S6.2(panel A) 
# Figure S8.2, Second row Table S8.1

#Warning 2 :the code work need a version of fastmap >= 1.2.0 
#if version of fastmap < 1.2.0, please reinstall it, uploading and then rerun
#the following code 

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

devtools::install_github("sarawade/mcclust.ext") #package for point estimates
library(mcclust.ext)

library(rstudioapi) # version 0.15.0
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#source mcmc codes for telescopic clustering models:
source("Conditional_t-HDP.R")
source("HDP.R")

#if true the code performs 100k iterations of the mcmc, 
#if false the code upload the outcome of the mcmc from a csv
run_MCMC = FALSE 

##SCENARIO N.B ##########################################################################

print("SIMULATION N.B")
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

#print(adj.rand.index(true_layer1, true_layer2))


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
print(ggdraw(p2))
print("Figure S5.1 - panel B just printed")

data1 = as.matrix(data_toy_layer1 - mean(data_toy_layer1))
data2 = as.matrix(data_toy_layer2 - mean(data_toy_layer2))
data = matrix(c( data_toy_layer1 - mean(data_toy_layer1),
                 data_toy_layer2 - mean(data_toy_layer2)), ncol = 2)

#estimate telescopic -  ##################################
print("SIMULATION N.B - t-HDP")
if(run_MCMC){
  totiter = 100000; burnin = 50000
  set.seed(0)
  m_saved =  telescopic_HDP_NN_uni(data, totiter = totiter)
  write.table(m_saved[,1,], file = paste("02B_layer_1.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
  write.table(m_saved[,2,], file = paste("02B_layer_2.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
}else{
  L = 2
  n_tot = 200
  totiter = 100000; burnin = 50000
  m_saved = array(NA,c(totiter, L, n_tot))
  print("MCMC")
  for(l in 1:L){
    temp = as.matrix(read.table(paste("output/02B_layer_",l,".csv", sep="")))
    m_saved[,l,] = temp
  }
}
print("Computing point estimate")
psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer1 = minVI(psm1, m_saved[(totiter-burnin+1):totiter,1,], method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, m_saved[(totiter-burnin+1):totiter,2,], method = "all")$cl[1,]
# write.table(VIminlayer1, file = paste("layer_1_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer2, file = paste("layer_2_cond_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
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
print(ggdraw(p2))
print("Figure S5.3 - panel C (same as Figure S6.2 - panel A) just printed")

#effective sample size - layer 1
ri= NULL
pb <- progress_bar$new(
  format = " Computing rand index per iteration [:bar] :percent Estimated completion time: :eta",
  total = totiter, clear = FALSE, width= 100)
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,1,], true_layer1))
  pb$tick()
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 1 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)
print("Figure S8.2 - panel A just printed")

#effectiveSize(ri)/totiter 
print("Simulation study n.B, Table S8.1 Layer 1:")
print(effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) )

#effective sample size - layer 2
ri= NULL
pb <- progress_bar$new(
  format = " Computing rand index per iteration [:bar] :percent Estimated completion time: :eta",
  total = totiter, clear = FALSE, width= 100)
for (iter in 1:totiter){
  ri = c(ri, rand.index(m_saved[iter,2,], true_layer2))
  pb$tick()
}
plot(ri, type = "l", xlab="MCMC iterations")
title(main="Rand Index layer 2 - Trace plot \n Conditional algorithm",
      cex.lab=0.75)
abline(v = burnin, col="red", lwd=3, lty=2)
print("Figure S8.2 - panel B just printed")

#effectiveSize(ri)/totiter 
print("Simulation study n.B, Table S8.1 Layer 2:")
print(effectiveSize(ri[(totiter-burnin+1):totiter])/(totiter-burnin) )


#estimate independent layers ##################################
print("SIMULATION N.B - Independent")
if(run_MCMC){
  group1 = c(rep(1, 200)); totiter = 1000; burnin = 500
  set.seed(1)
  c_layer1 = HDP_normal_normal(data1, group1, totiter = totiter)
  set.seed(2)
  c_layer2 = HDP_normal_normal(data2, group1, totiter = totiter)
  write.table(c_layer1, file = paste("04B_layer_1.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
  write.table(c_layer2, file = paste("04B_layer_2.csv",sep=""), row.names = FALSE,
              col.names= FALSE)
}else{
  L = 2
  n_tot = 200
  totiter = 1000; burnin = 500
  c_layer1 = as.matrix(read.table(paste("output/04B_layer_1.csv", sep="")))
  c_layer2 = as.matrix(read.table(paste("output/04B_layer_2.csv", sep="")))
}
#VI estimation of the clusters:
psm1 = psm(c_layer1[(totiter - burnin + 1):totiter,])
psm2 = psm(c_layer2[(totiter - burnin + 1):totiter,])
VIminlayer1 = minVI(psm1, c_layer1, method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, c_layer2, method = "all")$cl[1,]
# write.table(VIminlayer1, file = paste("layer_1_indep_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
# write.table(VIminlayer2, file = paste("layer_2_indep_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
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

data_plot_indep = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                             "Cluster1" = est_chr1,
                             "Cluster2" = est_chr2,
                             "mistakes" = mistakes)
#plot indep
pmain <- ggplot(data_plot_indep, aes(x = Layer1, y = Layer2,color = Cluster1, shape = Cluster2))+
  geom_point(size = 3 )+ geom_point(data = data_plot_indep %>% 
                                      filter(mistakes == TRUE),
                                    pch=21, 
                                    size=7, stroke = 2,
                                    colour="red")+
  theme_minimal(base_size = 18)+
  xlab("Layer n.1") + ylab("Layer n.2") +
  ggpubr::color_palette("default")  +
  labs(colour="Est. \nclusters at \nlayer n.1", 
       shape="Est. \nclusters at \nlayer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot_indep, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("default")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_indep, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("default")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
print(ggdraw(p2))
print("Figure S5.3 - panel A just printed")

#estimate constant layers#############################################################
print("SIMULATION N.B - Constant")
if(run_MCMC){
  group1 = c(rep(1, 200)); totiter = 1000; burnin = 500
  set.seed(1)
  c_layer_constant = HDP_normal_normal(data, group1, totiter = totiter)
  write.table(c_layer_constant, file = paste("05B_layer_1and2.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
}else{
  L = 2
  n_tot = 200
  totiter = 1000; burnin = 500
  c_layer_constant = as.matrix(read.table(paste("output/05B_layer_1and2.csv", sep="")))
}
psm = psm(c_layer_constant[(totiter - burnin + 1):totiter,])

VIminlayer = salso(c_layer_constant[(totiter - burnin + 1):totiter,], maxNClusters = 2)
VIminlayer1 = VIminlayer
VIminlayer2 = VIminlayer
# write.table(VIminlayer1, file = paste("layer_constant_pointEst.csv",sep=""), row.names = FALSE,
#             col.names= FALSE)
#take care of label switching to detect mistakes:
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

data_plot_const = data.frame("Layer1" = data_toy_layer1, "Layer2" = data_toy_layer2, 
                             "Cluster1" = est_chr1,
                             "Cluster2" = est_chr2,
                             "mistakes" = mistakes)
#plot indep
pmain <- ggplot(data_plot_const, aes(x = Layer1, y = Layer2,color = Cluster1, shape = Cluster2))+
  geom_point(size = 3 )+ geom_point(data = data_plot_const %>% 
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
  geom_density(data = data_plot_const, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("default")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_const, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("default")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
print(ggdraw(p2))
print("Figure S5.3 - panel B just printed")
