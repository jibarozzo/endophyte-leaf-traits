library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(vegan)
library(VennDiagram)
library(picante)
library(nlme)
library(plyr)
library(phyloseq)
library(DESeq2)
library(perm)
#+++++++++++++++++++++++++++++++++++++++++++
#BACTERIA
#+++++++++++++++++++++++++++++++++++++++++++
#********************************************
## Diversity and Richness, with rarefaction
#rarefaction done
#********************************************
#read data, below is old file with no plant trait data
seq.b<-read.csv(file="Soil2_16S_unrarefied.csv")

#below is data with plant trait but not working as of 1-29-2019
#seq.b<-read.csv(file="FinalData_Soil2_16S_unrarefied_2.csv")

#specify columns containing the OTUs which starts with X and New
sp.cols.b <- grep("X|New", names(seq.b))
#find out which number to rarefy by smallest sum
#rowSums(seq.b[,30:8146])

if (!identical(all.equal(sp.cols.b, round(sp.cols.b)), TRUE)) 
  stop("function is meaningful only for integers (counts)")
sp.cols.b <- as.matrix(sp.cols.b)
   sp.cols.b <- round(sp.cols.b)
if (ncol(sp.cols.b) == 1)
  sp.cols.b <- t(sp.cols.b)
if (length(sample) > 1 && length(sample) != nrow(sp.cols.b))


#colSums(seq.b[,30:8146])
#rarefy to 2500, add to a column
sp.col.rare.b<-rrarefy(seq.b[,sp.cols.b],15000)

#create new dataframe with sampleinfo merged data
#samp.b<-seq.b[,1:64]
samp.b<-seq.b[,1:29]
#merge the rarefied matrix with the data info
merge.b<-cbind(samp.b,sp.col.rare.b)
#write. as table
write.table(merge.b, file="Soil2_16S_Rarefied_data_notrait.csv", sep=',', row.names=FALSE)

#=========================================#
#add a column, and read in the new table
rare.bac<-read.csv(file="Soil2_16S_Rarefied_data_notrait.csv")
#specify columns containing the OTUs
sp.cols.bac <- grep("X|New", names(rare.bac))
#add a column
rare.bac$rich.vegan <-specnumber(rare.bac[,sp.cols.bac])
#take log of richness
rare.bac$logrich<-log(rare.bac$rich.vegan)
#shannon diversity
rare.bac$shan<-diversity(rare.bac[,sp.cols.bac], index="shannon")
#simpson
rare.bac$simpson<-diversity(rare.bac[,sp.cols.bac], index="simpson")
#log of shannon diversity
rare.bac$logshan<-log(rare.bac$shan)

write.table(rare.bac, file="Soil2_16S_Rarefied_data_notrait.csv", sep=',', row.names=FALSE)

#calculate chao1, then manually add to the column in the final data
chao.b<-estimateR(rare.bac[,sp.cols.bac])
write.table(chao.b, file="Bacteria_chao.csv", sep=',')

#++++++++++++++++++++++++++#
#t-test by season#
#==========================#
bac<-read.csv(file="Final_Soil2_16S_rarefied.csv")

t.test(shan~season,data=bac)
permTS(shan~season, alternative = "two.sided", method="exact.mc",B=999, data=bac)
t.test(simpson~season,data=bac)
t.test(S.chao1~season,data=bac)
