#HEATMAP
library(ggplot2)

setwd("J:/Working/UTK_PostdocProject/Bri's_Paper/MicrobiomeData/StatisticalAnalyses/BLAST_heatmap")

#+++++++++++++++++++++++++++++++++++++++++++#
#FUNGI
#+++++++++++++++++++++++++++++++++++++++++++#
#OTU table
otu<-read.csv(file="S2_ITS_assigned_class.csv")
otu.t<-t(otu)
write.table(otu.t, file="S2_ITS_assigned_class_transposed.csv", sep=',', col.names=FALSE)

#merge
otu1<-read.csv(file="S2_ITS_assigned_class_transposed.csv")
sampleinfo<-read.csv(file="S2_sample_map_BJ.csv")
data.t<-merge(sampleinfo, otu1, by="SampleID")
write.table(data.t, file="S2_ITS_assigned_class_unrarefied.csv", sep=',', row.names=FALSE)

#add plant trait data from Bri
fun.d<-read.csv(file="S2_ITS_assigned_class_unrarefied.csv")
trait<-read.csv(file="PlantTrait_Data.csv")
merged.fun<-merge(trait, fun.d, by="SampleID")
write.table(merged.fun, file="FinalS2_ITS_assigned_class_unrarefied.csv", sep=',', row.names=FALSE)

#Aggregate class abundances by provenance
its.class<-read.csv(file="S2_ITS_assigned_class_unrarefied.csv")
agg<-aggregate(.~provenance, data=its.class, sum)
write.table(agg, file="S2_ITS_aggregate_class.csv", sep=',', col.names=NA)

#manually edit the
agg.f<-read.csv(file="S2_ITS_aggregate_class.csv")
agg.f.t<-t(agg.f)
write.table(agg.f.t, file="S2_ITS_aggregate_class_transposed.csv", sep=',', col.names=NA)

#aggregate by class across all samples
class.d<-read.csv(file="S2_ITS_aggregate_class_transposed.csv")
class.agg<-aggregate(.~taxa, data=class.d, sum)
write.table(class.agg, file="Final_S2_ITS_aggregate_class_transposed.csv", sep=',', col.names=NA)

#PLOT HEATMAP
cls.prov<-read.csv(file="Final_S2_ITS_aggregate_class_transposed.csv")
levels(cls.prov)[levels(cls.prov$provenance) == "CNTRL"]<-"NC"

its.plot<-ggplot(cls.prov, aes(provenance, taxa)) +geom_tile(aes(fill = Percent), color = "white") +
    scale_fill_gradientn(colours = c("white","blue", "purple", "red"),
                     breaks=c(5,25,50,75,100),
                     na.value = "gray")+
  scale_x_discrete(labels=c("BJ","CL","CP","V", "NC"))+
  scale_y_discrete(limits = rev(levels(as.factor(cls.prov$taxa))))+
    theme(plot.background=element_blank(),panel.border=element_blank(),legend.text=element_text
          (size=5), legend.title = element_text(size=5),legend.key.size = unit(1.2, "lines"), 
          axis.text=element_text(size=10))+
  labs(x="",y="")
its.plot
#ggsave(its.plot,file="ITS_heatmap_2.tiff", width=120, height=70, units="mm", dpi=500)
ggsave(its.plot,file="ITS_heatmap_3.tiff", width=120, height=70, units="mm", dpi=500)

#+++++++++++++++++++++++++++++++++++++++++++#
#BACTERIA
#this is by phyla not but class, mislabeled
#+++++++++++++++++++++++++++++++++++++++++++#
#OTU table
otu.bac<-read.csv(file="S2_16S_assigned_class.csv")
otu.bac.t<-t(otu.bac)
otu.bac.t$SampleID<=otu.bac.t$taxonomy
write.table(otu.bac.t, file="S2_16S_assigned_class_transposed.csv", sep=',', col.names=FALSE)

#merge
otu1<-read.csv(file="S2_16S_assigned_class_transposed.csv")
sampleinfo<-read.csv(file="S2_sample_map_BJ.csv")
data.t<-merge(sampleinfo, otu1, by="SampleID")
write.table(data.t, file="S2_16S_assigned_class_unrarefied.csv", sep=',', row.names=FALSE)

#add plant trait data from Bri
fun.d<-read.csv(file="S2_16S_assigned_class_unrarefied.csv")
trait<-read.csv(file="PlantTrait_Data.csv")
merged.fun<-merge(trait, fun.d, by="SampleID")
write.table(merged.fun, file="FinalS2_16S_assigned_class_unrarefied.csv", sep=',', row.names=FALSE)

#Aggregate class abundances by provenance
its.class<-read.csv(file="FinalS2_16S_assigned_class_unrarefied.csv")
agg<-aggregate(.~provenance, data=its.class, sum)
write.table(agg, file="S2_16S_aggregate_class.csv", sep=',', col.names=NA)

#manually edit the
agg.f<-read.csv(file="S2_16S_aggregate_class.csv")
agg.f.t<-t(agg.f)
write.table(agg.f.t, file="S2_16S_aggregate_class_transposed.csv", sep=',', col.names=NA)

#===========================================================================#
#for some reason, control was not recognize during aggregating by provenance
cntrl<-read.csv(file="S2_16S_cntrl.csv")
cntrl.t<-t(cntrl)
write.table(cntrl.t, file="S2_16S_cntrl_transposed.csv", sep=',', col.names=NA)

cntrl.d<-read.csv(file="S2_16S_cntrl_transposed.csv")
cntrl.agg<-aggregate(.~taxa, data=cntrl.d, sum)
write.table(cntrl.agg, file="S2_16S_cntrl_aggregate_class_transposed.csv", sep=',', row.names=NA)
#===========================================================================#

#aggregate by class across all samples
class.d<-read.csv(file="S2_16S_aggregate_phyla_transposed.csv")
class.agg<-aggregate(.~taxa, data=class.d, sum)
write.table(class.agg, file="Final_S2_16S_aggregate_phyla_transposed.csv", sep=',', row.names=NA)

#================================#
#PLOT HEATMAP
#not below plots only the top 21 phyla for each provenance
#21 because the 20 and 21 are close in numners, it's arbitrary
#================================#

cls.prov<-read.csv(file="Final_S2_16S_top21plyla_eachprov.csv")

#change control to NC
levels(cls.prov)[levels(cls.prov$provenance) == "CNTRL"]<-"NC"

bac.plot<-ggplot(cls.prov, aes(provenance, taxa)) +geom_tile(aes(fill = Percent), color = "white") +
  scale_fill_gradientn(colours = c("white","blue", "purple", "red"),
                       breaks=c(5,25,50,75,100),
                       na.value = "gray")+
  scale_x_discrete(labels=c("BJ","CL","CP","V", "NC"))+
  scale_y_discrete(limits = rev(levels(as.factor(cls.prov$taxa))))+
  theme(plot.background=element_blank(),panel.border=element_blank(),legend.text=element_text
        (size=10), legend.title = element_text(size=10),legend.key.size = unit(1.5, "lines"), 
        axis.text=element_text(size=14))+
  labs(x="",y="")
bac.plot
ggsave(bac.plot,file="16S_heatmap_top21_a.tiff", width=140, height=170, units="mm", dpi=700)
