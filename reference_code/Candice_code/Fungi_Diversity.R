library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(vegan)
library(VennDiagram)
library(picante)
library(nlme)
library(plyr)
library(pls)
library(MuMIn)

#+++++++++++++++++++++++++++++++++++++++++++
#FUNGI
#+++++++++++++++++++++++++++++++++++++++++++
#********************************************
## Diversity and Richness, with rarefaction
#rarefaction done
#********************************************
#read data
seq<-read.csv(file="FinalData_Soil2_ITS_unrarefied.csv")

#specify columns containing the OTUs
sp.cols <- grep("SH", names(seq))
#rarefy to 2500, add add to a column
sp.col.rare<-rrarefy(seq[,sp.cols],1000)

#create new dataframe with sampleinfo merged data
samp<-seq[,1:64]
#merge the rarefied matrix with the data info
merge<-cbind(samp,sp.col.rare)
#write. as table
write.table(merge, file="Soil2_ITS_Rarefied_data.csv", sep=',', row.names=FALSE)

#=========================================#
#add a column, and read in the new table
rare.d<-read.csv(file="Soil2_ITS_Rarefied_data.csv")
#specify columns containing the OTUs
sp.cols <- grep("SH", names(rare.d))
#add a column
rare.d$rich.vegan <-specnumber(rare.d[,sp.cols])
#take log of richness
rare.d$logrich<-log(rare.d$rich.vegan)
#shannon diversity
rare.d$shan<-diversity(rare.d[,sp.cols], index="shannon")
#simpson
rare.d$simpson<-diversity(rare.d[,sp.cols], index="simpson")
#log of shannon diversity
rare.d$logshan<-log(rare.d$shan)


write.table(rare.d, file="FinalData_Soil2_ITS_rarefied.csv", sep=',', row.names=FALSE)

#calculate chao1, then manually add to the column in the final data
chao<-estimateR(rare.d[,sp.cols])
write.table(chao, file="Fungi_chao.csv", sep=',')

#++++++++++++++++++++++++++#
#t-test by season#
#==========================#
fungi<-read.csv(file="Final_Soil2_ITS_rarefied.csv")
t.test(shan~season,data=fungi)
permTS(shan~season, alternative = "two.sided", method="exact.mc",B=999, exact=TRUE, data=fungi)
t.test(simpson~season,data=fungi)
t.test(S.chao1~season,data=fungi)

#=============================================#
#Test for influencing alpha diversity of fungi
#=============================================#

fundiv<-read.csv(file="FinalData_Soil2_ITS_rarefied_noNAs.csv")
