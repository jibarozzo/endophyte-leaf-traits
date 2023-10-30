library(reshape2)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(vegan)
library(VennDiagram)
library(picante)
library(ape)
library(nlme)
library(plyr)

##################################################
#Since sample are labeled differently in the OTU table than in the 
#environmental data, I'll have to merge
#All OTU tables used were the ones filtered for < 100 sequences 
#i.e. samples with <100 seq were excluded from the dataset
##################################################

#FUNGI
#+++++++++++++++++++++++++++++++++++++++++++
#OTU table
otu<-read.csv(file="Soil2_ITS.csv")
otu.t<-t(otu)
write.table(otu.t, file="Soil2_ITS_transposed.csv", sep=',', col.names=FALSE)

#merge
otu1<-read.csv(file="Soil2_ITS_transposed.csv")
sampleinfo<-read.csv(file="S2_sample_map_BJ.csv")
data.t<-merge(sampleinfo, otu1, by="SampleID")
write.table(data.t, file="Soil2_ITS_unrarefied.csv", sep=',', row.names=FALSE)

#add plant trait data from Bri
fun.d<-read.csv(file="Soil2_ITS_unrarefied.csv")
trait<-read.csv(file="PlantTrait_Data.csv")
merged.fun<-merge(trait, fun.d, by="SampleID")
write.table(merged.fun, file="FinalData_Soil2_ITS_unrarefied.csv", sep=',', row.names=FALSE)

#+++++++++++++++++++++++++++++++++++++++++++
#BACTERIA
#+++++++++++++++++++++++++++++++++++++++++++
#OTU table
botu<-read.csv(file="Soil2_16S.csv")
botu.t<-t(botu)
write.table(botu.t, file="Soil2_16S_transposed_2.csv", sep=',', col.names=FALSE)

#merge
botu1<-read.csv(file="Soil2_16S_transposed_2.csv")
sampleinfo<-read.csv(file="S2_sample_map_BJ.csv")
data.t<-merge(sampleinfo, botu1, by="SampleID")
write.table(data.t, file="Soil2_16S_unrarefied_2.csv", sep=',', row.names=FALSE)

#add plant trait data from Bri
bac.d<-read.csv(file="Soil2_16S_unrarefied_2.csv")
trait<-read.csv(file="PlantTrait_Data.csv")
merged.bac<-merge(trait, bac.d, by="SampleID")
write.table(merged.bac, file="FinalData_Soil2_16S_unrarefied_2.csv", sep=',', row.names=FALSE)
