# Clustering Analysis 

#Package

library(tidyverse)
library(factoextra)
library(compareGroups)
library(FactoMineR)
library(factoextra)
library(survival)
library(survminer)
library(ggplot2)

#data with lipids only 

#Subset data with HFpEF patients 
HFpEF<-subset (Became_short_complete, filtre_lipido==1 & groupe_final==2)
names(HFpEF)

# Subset data with only lipids  = d
d<-HFpEF[,407:649]

names(d)
view(d)

#Hierchical clustering analysis 

#Compenent principal analysis 
res.pca <- PCA((d), scale=TRUE, ncp = 242, graph = TRUE)

# Percentage variance per dimension 
res.pca$eig

#Hierchical clustering with the opitmal number of cluster (using the inertia gain)
res.hcpc <- HCPC(res.pca, graph = FALSE)

#Dendidogramm with the ineria gain
res.PCA.hcpc <- HCPC(res=res.pca, nb.clust=-1, consol=FALSE, min=3, max=10, graph=TRUE, metric="euclidean", method="ward", iter.max=10, nb.par=5)

#Visualization of clusters on PCA
fviz_cluster(res.hcpc, repel = TRUE, show.clust.cent = TRUE,palette = "jco", ggtheme = theme_minimal(), 
             main = "Factor map", ggrepel.max.overlaps = Inf)

#"Clust" row added to the data with clinical variables 
x<-cbind(HFpEF, res.hcpc$data.clust) 
names(x) 

#Survival analysis (Kaplan Meier Curves)
x$TimeCombinedYears<-x$TimeCombinedMonths/12
x$TimeDeath_Nat_Hist_Year

# Primary endpoint (All-causes mortality and HF hospitalization)
sfit<-survfit(Surv(x$TimeCombinedYears, Combined)~clust, data=x)
sfit

# Survival curve 1
ggsurvplot(sfit, conf.int =FALSE, pval=TRUE, risk.table = TRUE,
           legend.labs=c("Cluster 1", "Cluster 2", "Cluster 3"), legend.title="HFpEF",
           palette=c("blue","grey","#FF9933"), 
           title="Kaplan-Meier Curve for mortality from all-causes and for HF hospitalization", 
           risk.table.height=.20, break.time.by=5,xlim=c(0,5))

# Secondary endpoint (All-causes mortality)
sfit<-survfit(Surv(x$TimeDeath_Nat_Hist_Year, Death)~clust, data=x)
sfit

# Survival curve 2
ggsurvplot(sfit, conf.int =FALSE, pval=TRUE, risk.table = TRUE,
           legend.labs=c("Cluster 1", "Cluster 2", "Cluster 3"), legend.title="HFpEF",
           palette=c("blue","grey","#FF9933"), 
           title="Kaplan-Meier Curve for mortality from all-causes", 
           risk.table.height=.20, break.time.by=5,xlim=c(0,5))
