#Simulation scenario 1 with two layers (true clustering)

#library(progress)#to draw the progress bar

library(mcclust.ext) #to minimize VI
library(T4cluster) #to compute psm with psm() function
library(fossil) #to compute rand indexes
library(coda) #to compute the effective sample size

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(ggplot2)

source("Conditional_t-HDP.R")
source("HDP.R")

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
set.seed(1)
m_saved =  telescopic_HDP_NN_uni(data, totiter = totiter)
write.table(m_saved[,1,], file = paste("layer_1_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(m_saved[,2,], file = paste("layer_2_cond.csv",sep=""), row.names = FALSE,
            col.names= FALSE)

psm1 = psm(m_saved[(totiter-burnin+1):totiter,1,])
psm2 = psm(m_saved[(totiter-burnin+1):totiter,2,])
VIminlayer1 = minVI(psm1, m_saved[(totiter-burnin+1):totiter,1,], method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, m_saved[(totiter-burnin+1):totiter,2,], method = "all")$cl[1,]
write.table(VIminlayer1, file = paste("layer_1_cond_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2_cond_pointEst.csv",sep=""), row.names = FALSE,
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

#estimate independent layers ##################################
group1 = c(rep(1, 200)); totiter = 10000; burnin = 5000
set.seed(1)
c_layer1 = HDP_normal_normal(data1, group1, totiter = totiter)
set.seed(2)
c_layer2 = HDP_normal_normal(data2, group1, totiter = totiter)
#VI estimation of the clusters:
psm1 = psm(c_layer1[(totiter - burnin + 1):totiter,])
psm2 = psm(c_layer2[(totiter - burnin + 1):totiter,])
VIminlayer1 = minVI(psm1, c_layer1, method = "all")$cl[1,]
VIminlayer2 = minVI(psm2, c_layer2, method = "all")$cl[1,]
write.table(VIminlayer1, file = paste("layer_1_indep_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
write.table(VIminlayer2, file = paste("layer_2_indep_pointEst.csv",sep=""), row.names = FALSE,
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
  ggpubr::color_palette("jco")  +
  labs(colour="Est. \nclusters at \nlayer n.1", 
       shape="Est. \nclusters at \nlayer n.2") 
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot_indep, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_indep, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

#estimate constant layers#############################################################
group1 = c(rep(1, 200)); totiter = 10000; burnin = 5000
set.seed(3)
c_layer_constant = HDP_normal_normal(data, group1, totiter = totiter)

psm = psm(c_layer_constant[(totiter - burnin + 1):totiter,])

VIminlayer = minVI(psm, c_layer_constant[(totiter - burnin + 1):totiter,], method = "all")$cl[1,]   
VIminlayer1 = VIminlayer
VIminlayer2 = VIminlayer
write.table(VIminlayer1, file = paste("layer_constant_pointEst.csv",sep=""), row.names = FALSE,
            col.names= FALSE)
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
  ggpubr::color_palette("jco")  +
  labs(colour="Est. clusters at \n layer n.1", 
       shape="Est. clusters at \n layer n.2")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = data_plot_const, aes(x = Layer1, fill = Cluster1),
               alpha = 0.7, size = 0.2)+
  ggpubr::fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = data_plot_const, aes(x = Layer2, fill = Cluster2),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  ggpubr::fill_palette("jco")
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
