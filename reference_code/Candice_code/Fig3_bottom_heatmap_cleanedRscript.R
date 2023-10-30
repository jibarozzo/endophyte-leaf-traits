#Figures
library(ggplot2)

setwd("J:/Working/UTK_PostdocProject/Bri's_Paper/MicrobiomeData/StatisticalAnalyses/BLAST_heatmap")

#+++++++++++++++++++++++++++++++++++++++++++#
#FUNGI
#+++++++++++++++++++++++++++++++++++++++++++#

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
