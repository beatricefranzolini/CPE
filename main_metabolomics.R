
rm(list = ls())
library(purrr) #for k-means exploratory
library(cluster) #for k-means exploratory
library(rstatix) #for t-test exploratory

library(progress)#to draw the progress bar

library(mcclust.ext) #to minimize Binder and VI
library(T4cluster) #to compute psm with psm() function
library(fossil) #to compute rand indexes
library(coda) #to compute the effective sample size

library(cowplot) #to plot
library(ggpubr) #to plot 
library(dplyr) #used in plot
library(pheatmap) #to plot

library(MASS)

library(CholWishart)
library(GGally) #to plot multivariate continuous data
library(zoo) #for linear interpolation (missing data)
library("factoextra")#to extact PCA

library(tidyverse) #for test-t
library(zoo) #linear interpolation

source("Conditional_t-HDP.R")

BMI_rawdata = read.table("metabolomics/data/BMI_zscores_20200608.csv", sep=",", header = TRUE)

bmi =  BMI_rawdata[,c(1,13:22)]


MUM_rawdata = read.table("metabolomics/data/FormA229_glycemia_demographics_2019-09-04.csv",
                            sep=";", header = TRUE)
mum = data.frame("SubjectID" = MUM_rawdata$SubjectID, 
                   "ogtt" = MUM_rawdata$ogtt_fasting_pw26,
                 "BMIpp" = MUM_rawdata$ppBMI)

METABOL_rawdata = read.table("metabolomics/data/CHME_LAB_META_20200714_yr8.csv",
                             sep=",", header = TRUE)

metabolites = data.frame("SubjectID" = METABOL_rawdata$SubjectID, 
                         "mb1" = METABOL_rawdata$CME26CLDL,
                         "mb2" = METABOL_rawdata$CME26HDL,
                         "mb3" = METABOL_rawdata$CME26TTG,
                         "mb4" = METABOL_rawdata$CME26PG,
                         "mb5" = METABOL_rawdata$CME26TOCH,
                         "mb6" = METABOL_rawdata$CME26SPH,
                         "mb7" = METABOL_rawdata$CME26APOA,
                         "mb8" = METABOL_rawdata$CME26APOB,
                         "mb9" = METABOL_rawdata$CME26OME3,
                         "mb10" = METABOL_rawdata$CME26OME6,
                         "mb11" = METABOL_rawdata$CME26PUFP,
                         "mb12" = METABOL_rawdata$CME26MUFP,
                         "mb13" = METABOL_rawdata$CME26SFAP,
                         "mb14" = METABOL_rawdata$CME26LAP,
                         "mb15" = METABOL_rawdata$CME26DHAP,
                         "mb16" = METABOL_rawdata$CME26ALA,
                         "mb17" = as.numeric(METABOL_rawdata$CME26GLN),
                         "mb18" = METABOL_rawdata$CME26GLY,
                         "mb19" = METABOL_rawdata$CME26HIS,
                         "mb20" = METABOL_rawdata$CME26ILE,
                         "mb21" = METABOL_rawdata$CME26LEU,
                         "mb22" = as.numeric(METABOL_rawdata$CME26VAL),
                         "mb23" = METABOL_rawdata$CME26PHE,
                         "mb24" = METABOL_rawdata$CME26TYR,
                         "mb25" = METABOL_rawdata$CME26GLU,
                         "mb26" = METABOL_rawdata$CME26LAC,
                         "mb27" = as.numeric(METABOL_rawdata$CME26PRY),
                         "mb28" = METABOL_rawdata$CME26CITR,
                         "mb29" = as.numeric(METABOL_rawdata$CME26BHB),
                         "mb30" = METABOL_rawdata$CME26OAC,
                         "mb31" = METABOL_rawdata$CME26ACA,
                         "mb32" = METABOL_rawdata$CME26AC,
                         "mb33" = as.numeric(METABOL_rawdata$CME26CRN),
                         "mb34" = METABOL_rawdata$CME26ALB,
                         "mb35" = METABOL_rawdata$CME26GLYA)

data_withNA = merge(merge(bmi, mum, by = "SubjectID"), metabolites, by = "SubjectID")
sum(rowSums(is.na(data_withNA))>=1)/dim(data_withNA)[1]
dim(data_withNA)
set.seed(0)
data = data_withNA[rowSums(is.na(data_withNA))<floor(0.2*(dim(data_withNA)[2]-1)),-1]
data = na.approx(data)
data = data[rowSums(is.na(data))==0,]

data = scale(data)

data.pca <- prcomp(data[,13:47])
summary(data.pca)
temp = get_pca_ind(data.pca)
eigs <- data.pca$sdev^2
plot(eigs / sum(eigs), type = "l")
summary(lm(data[,9]~temp$coord))
#6 PCA
data_to_clust = data[,1:18]
data_to_clust[,13:18]=temp$coord[,1:6]
corrplot::corrplot(cor(data_to_clust))

colayer = c(rep(1,10),rep(2,2),rep(3,6))
  
#run MCMC 
totiter = 100000
set.seed(1)
m_saved = telescopic_HDP_NNIX_multi_3L1P(data_to_clust, colayer, 
                                                  totiter = totiter)

#point estimates
#l=1
layer1_chain = read.csv("layer1.csv", sep =" ", header = FALSE)
layer1_chain = as.matrix(layer1_chain)
dim(layer1_chain)

set.seed(1)
layer1 = layer1_chain[seq(50005,100000,5),]
psm1 = psm(layer1)
hist(psm1)
tempBI1 = minbinder(psm1, layer1, method = "all")$cl
tempIV1 = minVI(psm1, layer1, method = "all")$cl

clust_layer1 = tempIV1[1,]
write.table(clust_layer1, file = paste("layer1_point.csv"), row.names = FALSE,
            col.names= FALSE)

#l=2
layer2_chain = read.csv("layer2.csv", sep =" ", header = FALSE)
layer2_chain = as.matrix(layer2_chain)
dim(layer2_chain)

set.seed(1)
layer2 = layer2_chain[seq(50005,100000,5),]
psm2 = psm(layer2)
hist(psm2)
tempBI2 = minbinder(psm2, layer2, method = "all")$cl
tempIV2 = minVI(psm2, layer2, method = "all")$cl

clust_layer2 = tempBI2[3,]
write.table(clust_layer2, file = paste("layer2_point.csv"), row.names = FALSE,
            col.names= FALSE)

#l=3
layer3_chain = read.csv("layer3.csv", sep =" ", header = FALSE)
layer3_chain = as.matrix(layer3_chain)
dim(layer3_chain)

set.seed(1)
layer3 = layer3_chain[seq(50005,100000,5),]
psm3 = psm(layer3)
hist(psm3)
tempIV3 = minVI(psm3, layer3, method = "all")$cl
tempBI3 = minbinder(psm3, layer3, method = "all")$cl

clust_layer3 = tempBI3[3,]
write.table(clust_layer3, file = paste("layer3_point.csv"), row.names = FALSE,
            col.names= FALSE)

#clusters' frequencies
prop.table(table(clust_layer1))
prop.table(table(clust_layer2))
prop.table(table(clust_layer3))

prop.table(table(clust_layer1,clust_layer2), 1)
prop.table(table(clust_layer1,clust_layer3), 1)


prop.table(table(clust_layer1,clust_layer2), 2)

#FIGURES 
data =  t(t(data)*attr(data,"scaled:scale")) +  rep(attr(data,"scaled:center"),
                                                    each = dim(data)[1]) 

#scatter plots
ggpairs(as.data.frame(data[,seq(2,10,2)]), aes(color = as.factor(clust_layer1),  # Color by group (cat. variable)
                                       alpha = 0.5))+ theme_classic()

ggpairs(as.data.frame(data[,11:12]), aes(color = as.factor(clust_layer2),  # Color by group (cat. variable)
                                         alpha = 0.5))+ theme_classic()

set.seed(1)
ggpairs(as.data.frame(data[,sample(seq(13,47),7)]), aes(color = as.factor(clust_layer3),  # Color by group (cat. variable)
                                         alpha = 0.5))+ theme_classic()

#Figures #LAYER 1
Z_times = c(3, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9)
quantile10 = aggregate(data[,1:10]~factor(clust_layer1), 
                 FUN = 'quantile', probs=c(0.1) )
quantile90 = aggregate(data[,1:10]~factor(clust_layer1), 
                       FUN = 'quantile', probs=c(0.90) )
meanZBMI = aggregate(data[,1:10]~factor(clust_layer1), FUN=mean) 
for(clust in unique(clust_layer1)){
  ZBMI = data.frame(x = Z_times, 
                    meanZBMI = as.numeric(meanZBMI[clust, 2:11]), 
            p10 = as.numeric(quantile10[clust, 2:11]),
            p90 = as.numeric(quantile90[clust, 2:11]))
  print(ggplot(ZBMI, aes(x = x, y = meanZBMI)) + geom_line(aes(y = meanZBMI),linetype = "dashed") + 
    geom_line(aes(y = p10),linetype = "dashed") + 
    geom_line(aes(y = p90),linetype = "dashed") +
    geom_point( ) + 
    theme_classic()  + ylim(-4, 6) + theme(text = element_text(size = 20))+
    labs(x="Years", y = "Z-BMI") + 
    annotate("rect", xmin = 3, xmax = 9, ymin = -2, ymax = 1,
             alpha = .15, fill = "green") + 
    annotate("rect", xmin = 3, xmax = 9, ymin = 1, ymax = 2,
             alpha = .15, fill = "orange") + 
    geom_text(x = 6, y = 1.5, label="Overweight", size=7) +
    annotate("rect", xmin = 3, xmax = 9, ymin = -2, ymax = -3,
             alpha = .15, fill = "orange") + 
    geom_text(x = 6, y = -2.5, label="Thinness", size=7) +
    annotate("rect", xmin = 3, xmax = 9, ymin = -3, ymax = -4,
             alpha = .15, fill = "red") + 
    geom_text(x = 6, y = -3.5, label="Severe thinness", size=7) +
    annotate("rect", xmin = 3, xmax = 9, ymin = 2, ymax = 6,
             alpha = .15, fill = "red") + 
    geom_text(x = 6, y = 2.5, label="Obesity", size=7) +
    geom_ribbon(aes(x = x,
                    ymin = p10,
                    ymax = p90),
                fill = "gray",
                alpha = 0.5) + 
    geom_text(x = 6, y = 5, 
      label=paste(round(prop.table(table(clust_layer1))[clust],2)*100,"%"),
      size=7, color="red"))
}

#small figures

for(clust in unique(clust_layer1)){
  ZBMI = data.frame(x = Z_times, 
                    meanZBMI = as.numeric(meanZBMI[clust, 2:11]), 
                    p10 = as.numeric(quantile10[clust, 2:11]),
                    p90 = as.numeric(quantile90[clust, 2:11]))
  print(ggplot(ZBMI, aes(x = x, y = meanZBMI)) + geom_line(aes(y = meanZBMI),linetype = "dashed") + 
          geom_line(aes(y = p10),linetype = "dashed") + 
          geom_line(aes(y = p90),linetype = "dashed") +
          geom_point( ) + 
          theme_classic()  + ylim(-4, 6) + theme(text = element_text(size = 20))+
          labs(x="Years", y = "Z-BMI") + 
          annotate("rect", xmin = 3, xmax = 9, ymin = -2, ymax = 1,
                   alpha = .15, fill = "green") + 
          annotate("rect", xmin = 3, xmax = 9, ymin = 1, ymax = 2,
                   alpha = .15, fill = "orange") +
          annotate("rect", xmin = 3, xmax = 9, ymin = -2, ymax = -3,
                   alpha = .15, fill = "orange") +
          annotate("rect", xmin = 3, xmax = 9, ymin = -3, ymax = -4,
                   alpha = .15, fill = "red")  +
          annotate("rect", xmin = 3, xmax = 9, ymin = 2, ymax = 6,
                   alpha = .15, fill = "red")  +
          geom_ribbon(aes(x = x,
                          ymin = p10,
                          ymax = p90),
                      fill = "gray",
                      alpha = 0.5))
}

#Figures #LAYER 2

MUM = data.frame(Cluster = factor(clust_layer2), 
                  ogtt = data[,11], 
                  ppBMI = data[,12])

ggplot(MUM, aes(x = Cluster, y = ogtt, fill=Cluster)) +
  geom_boxplot() + theme_classic() + theme(text = element_text(size = 20)) + 
  geom_text(x = 1.3, y = 6.8, label="1%", size=7, color="red") +
  geom_text(x = 2.3, y = 4.6, label="71%", size=7, color="red") +
  geom_text(x = 3.3, y = 5.1, label="28%", size=7, color="red")

ggplot(MUM, aes(x = Cluster, y = ppBMI, fill=Cluster)) +
  geom_boxplot() + theme_classic() + theme(text = element_text(size = 20)) + 
  geom_text(x = 1.3, y = 29.7, label="1%", size=7, color="red") +
  geom_text(x = 2.3, y = 23.8, label="71%", size=7, color="red") +
  geom_text(x = 3.3, y = 31.4, label="28%", size=7, color="red")

#small figures #check colors
var1_cluster1 = ggplot() + 
  geom_boxplot(aes(y = MUM$ogtt[clust_layer2==1]), fill="#F8766D" ) + 
  ylim(c(3, 8.5)) +
  labs(x=" ", y = " ") + theme_classic()

var1_cluster2 = ggplot() + 
  geom_boxplot(aes(y = MUM$ogtt[clust_layer2==2]), fill="#00BA38" )+
  ylim(c(3, 8)) +
  labs(x=" ", y = " ") + theme_classic()

var1_cluster3 = ggplot() + 
  geom_boxplot(aes(y = MUM$ogtt[clust_layer2==3]), fill="#619CFF" )+
  ylim(c(3, 8)) +
  labs(x=" ", y = " ") + theme_classic()


var2_cluster1 = ggplot() + 
  geom_boxplot(aes(y = MUM$ppBMI[clust_layer2==1]), fill="#F8766D" ) + 
  ylim(c(15, 42)) +
  labs(x=" ", y = " ") + theme_classic()

var2_cluster2 = ggplot() + 
  geom_boxplot(aes(y = MUM$ppBMI[clust_layer2==2]), fill="#00BA38" )+
  ylim(c(15, 42)) +
  labs(x=" ", y = " ") + theme_classic()

var2_cluster3 = ggplot() + 
  geom_boxplot(aes(y = MUM$ppBMI[clust_layer2==3]), fill="#619CFF" )+
  ylim(c(15, 42)) +
  labs(x=" ", y = " ") + theme_classic()


plot_grid(var1_cluster1, var2_cluster1)
plot_grid(var1_cluster2, var2_cluster2)
plot_grid(var1_cluster3, var2_cluster3)


#Figures #LAYER 3
#small figures
temp = aggregate(scale(data[,13:47])~factor(clust_layer3), FUN=mean)
mat1 = matrix(as.numeric(temp[1,2:36]), nrow = 7)
pheatmap(mat1, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F)
mat2 = matrix(as.numeric(temp[2,2:36]), nrow = 7)
pheatmap(mat2, display_numbers=F, show_colnames=F, 
         cluster_rows=F, cluster_cols=F, color = hcl.colors(50, "BluYl"),
         legend = F)

data_temp = data.frame(cbind(clust_layer3,data[,13:47])) 
mydata = data_temp %>% as_tibble()
mydata.long = mydata %>%
  pivot_longer(-clust_layer3, names_to = "variables", values_to = "value")
stat.test = mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ clust_layer3) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

clust_layer3 = read.table("metabolomics/results/layer3_point.csv")
sum_met_clust = data.frame(cbind(data[,13:47],factor(clust_layer3$V1)))
names_vector <- c("Clinical LDL Cholesterol", "HDL Cholesterol", "Triglycerides", "Phosphoglycerides", "Cholines Phosphoglycerides", "Sphingomyelins", "APO A1", "APO B", "Omega 3", "Omega 6", "Poly-Unsaturated FA (PUFA)", "Mono-Unsaturated FA (MUFA)", "Saturated FA (SFA)", "Linoleic acid", "Docosahexaenoic acid (DHA)", "Alanine", "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Valine", "Phenylalanine", "Tyrosine", "Glucose", "Lactate", "Pyruvate", "Citrate", "beta-Hydroxybutyric acid (bOHbutyrate)", "Acetate", "Acetoacetate", "Acetone", "Creatinine", "Albumin", "Glycoprotein acetyls")
colnames(sum_met_clust) = c(names_vector, "cluster")

tab_mean = t(sum_met_clust %>% group_by(cluster) %>%
  summarise(across(everything(), mean)))
tab_sd = t(sum_met_clust %>% group_by(cluster) %>%
               summarise(across(everything(), sd)))
tab_median = t(sum_met_clust %>% group_by(cluster) %>%
             summarise(across(everything(), median)))
tab_IQR = t(sum_met_clust %>% group_by(cluster) %>%
                 summarise(across(everything(), IQR)))
tab_kruskal = NULL
for (met in 1:35){
  tab_kruskal = c(tab_kruskal, kruskal.test(sum_met_clust[,met] ~ sum_met_clust$cluster)$p.value)
}

sum_met = cbind(tab_mean[-1,], tab_IQR[-1,], tab_kruskal)

write.table(round(sum_met,4), file = paste("table_metabolites.csv"))

group_by(sum_met_clust, cluster) %>%
  summarise(
    count = n(),
    mean = mean(weight, na.rm = TRUE),
    sd = sd(weight, na.rm = TRUE),
    median = median(weight, na.rm = TRUE),
    IQR = IQR(weight, na.rm = TRUE)
  )
library(dplyr)
