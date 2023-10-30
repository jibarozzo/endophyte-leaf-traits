library(ggplot2)
library(gridExtra)
library(cowplot)
library(vegan)
library(lme4)
library(MuMIn)
library(nlme)
library(plyr)
library(tidyverse)
library(ggpubr)
#=====================#
#Variance Partitioning
#=====================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#1. FUNGI#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
fungi<-read.csv(file="FinalData_Soil2_ITS_rarefied.csv")
#select summer
#fungi.summer<-fungi[fungi$season=="SUMMER",]
#remove NA
fungi <- na.omit(fungi)
spp.list.fungi<-names(fungi)[grep("SH",names(fungi))] ## find and list all variables
#spp.list
fung.mat<-fungi[,spp.list.fungi]  ## create species matrix

#calculate distance matrix
fung.dist<-vegdist(fung.mat, method="bray")

#transform and calculate some variables that Bri used
fungi$CNag<-log(fungi$AGC/fungi$AGN)
fungi$CNbg<-log(fungi$BGC/fungi$BGN)
fungi$BGN.AGN<-log(fungi$BGN/fungi$AGN)
fungi$BGC.AGC<-log(fungi$BGC/fungi$AGC)

#transform variables
fungi$AGBsqrt<-1/sqrt(fungi$AGBiomass)
fungi$logBGTotal<-log(fungi$BGTotal)
fungi$logunder30N<-log(fungi$under30N)
fungi$logdiam<-log(fungi$diam)
fungi$lognoLeaves<-log(fungi$noLeaves)

fungi$loglengthLeaf<-log(fungi$lengthLeaf)
fungi$logover30N<-log(fungi$over30N)
fungi$logHeightOver30<-log(fungi$HeightOver30)
fungi$logRtSt<-log(fungi$RtSt)

fungi$logseedweight<-log(fungi$seedweight)
fungi$logHeightUnder30<-log(fungi$HeightUnder30)
fungi$CNag<-log(fungi$AGC/fungi$AGN)
fungi$CNbg<-log(fungi$BGC/fungi$BGN)


#specify fungi the position is not numeric, but as factor
fungi$Position<-as.factor(fungi$Position)
fungi$provenance<-as.factor(fungi$provenance)

#partition by trait, provenance
fungi.varpart<-varpart(fung.dist,~logBGTotal+logdiam+AGBsqrt+lognoLeaves+RtSt+
                         logseedweight+CNag+CNbg,
                       ~Total.relevant.PAHs,
                          ~provenance,
                          ~season,
                          sqrt.dist = FALSE, chisquare = TRUE,data=fungi)
fungi.varpart


#plot in pdf
#pdf("Fungi_all_varpart.pdf",width=10.5,height=5) 
#par(mfrow=c(1,2), mar=c(1,2,1,2), bty='n')
plot(fungi.varpart,Xnames=c("Traits","Residual oil" ,"Provenance","Season"))
title(main="Fungi")

#rda.result <- rda(fung.dist ~ X1 + Condition(X2) +
#                    Condition(as.matrix(mite.pcnm)))
#anova(rda.result, step=200, perm.max=200)

#========================================================================================================#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#1.1 Summer FUNGI#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
fungi<-read.csv(file="FinalData_Soil2_ITS_rarefied.csv")
#select summer
fungi.summer<-fungi[fungi$season=="SUMMER",]
#remove NA
fungi.summer <- na.omit(fungi.summer)
spp.list.fungsum<-names(fungi.summer)[grep("SH",names(fungi.summer))] ## find and list all variables
#spp.list
fungsum.mat<-fungi.summer[,spp.list.fungsum]  ## create species matrix

#calculate distance matrix
dist.fungisum<-vegdist(fungsum.mat, method="bray")

#transform variables
fungi.summer$AGBsqrt<-1/sqrt(fungi.summer$AGBiomass)
fungi.summer$logBGTotal<-log(fungi.summer$BGTotal)
fungi.summer$logunder30N<-log(fungi.summer$under30N)
fungi.summer$logdiam<-log(fungi.summer$diam)
fungi.summer$lognoLeaves<-log(fungi.summer$noLeaves)

fungi.summer$loglengthLeaf<-log(fungi.summer$lengthLeaf)
fungi.summer$logover30N<-log(fungi.summer$over30N)
fungi.summer$logHeightOver30<-log(fungi.summer$HeightOver30)
fungi.summer$logRtSt<-log(fungi.summer$RtSt)

fungi.summer$logseedweight<-log(fungi.summer$seedweight)
fungi.summer$logHeightUnder30<-log(fungi.summer$HeightUnder30)
fungi.summer$CNag<-log(fungi.summer$AGC/fungi.summer$AGN)
fungi.summer$CNbg<-log(fungi.summer$BGC/fungi.summer$BGN)


#specify that the position is not numeric, but as factor
fungi.summer$Position<-as.factor(fungi.summer$Position)
fungi.summer$provenance<-as.factor(fungi.summer$provenance)
#partition by trait, provenance
#combine Salinity and Water Level as the abiotic properties
#HostDensity, forest cover and Woody volume as the biotic properties
fungi.sum.varpart<-varpart(dist.fungisum,~logBGTotal+logdiam+AGBsqrt+lognoLeaves+RtSt+
                             logseedweight+CNag+CNbg,
                           ~provenance,
                           ~Total.relevant.PAHs,
                           sqrt.dist = FALSE, chisquare = TRUE,data=fungi.summer)
fungi.sum.varpart

#plot in pdf
#pdf("Fungi_Summer_VariancePartition.pdf",width=10.5,height=5) 
#par(mfrow=c(1,2), mar=c(1,2,1,2), bty='n')
plot(fungi.sum.varpart,Xnames=c("Traits","Provenance"))
#title(main="Fungi_summer")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#1. BACTERIA #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
bact<-read.csv(file="Final_Soil2_16S_rarefied.csv")
#select summer

#remove NA
bact <- na.omit(bact)
#write.table(bact, file="SubsetFinal_soil_16s_rarefied.csv", sep=',')

#read in the distance matrix without samples with NA to match that number of samples in bact
bac.dist<-read.csv(file="unifrac_weighted_subset_no_NA.csv", header=FALSE)
#dist<-as.dist(dist, diag=FALSE, upper=FALSE)

#dist<-read.table(file="weighted_unifrac_Soil2_16S_100.txt")
#bac.dist <- as.dist(as(dist, "matrix"))
#bac.dist<-as.matrix(dist, diag=TRUE, upper=TRUE)

#transform and calculate some variables that Bri used
bact$CNag<-log(bact$AGC/bact$AGN)
bact$CNbg<-log(bact$BGC/bact$BGN)
bact$BGN.AGN<-log(bact$BGN/bact$AGN)
bact$BGC.AGC<-log(bact$BGC/bact$AGC)

#transform variables
bact$AGBsqrt<-1/sqrt(bact$AGBiomass)
bact$logBGTotal<-log(bact$BGTotal)
bact$logunder30N<-log(bact$under30N)
bact$logdiam<-log(bact$diam)
bact$lognoLeaves<-log(bact$noLeaves)

bact$loglengthLeaf<-log(bact$lengthLeaf)
bact$logover30N<-log(bact$over30N)
bact$logHeightOver30<-log(bact$HeightOver30)
bact$logRtSt<-log(bact$RtSt)

bact$logseedweight<-log(bact$seedweight)
bact$logHeightUnder30<-log(bact$HeightUnder30)


#specify fungi the position is not numeric, but as factor
bact$Position<-as.factor(bact$Position)
bact$provenance<-as.factor(bact$provenance)

#partition by trait, provenance
bact.varpart<-varpart(bac.dist,~logBGTotal+logdiam+AGBsqrt+lognoLeaves+RtSt+
                        logseedweight+CNag+CNbg,
                      ~Total.relevant.PAHs,
                      ~provenance,
                      ~season,
                       sqrt.dist = FALSE, chisquare = TRUE, data=bact)
bact.varpart

#plot in pdf
#pdf("Fungi_all_varpart.pdf",width=10.5,height=5) 
#par(mfrow=c(1,2), mar=c(1,2,1,2), bty='n')
plot(bact.varpart,Xnames=c("Traits","Provenance","Season"))
title(main="Bacteria")

#=