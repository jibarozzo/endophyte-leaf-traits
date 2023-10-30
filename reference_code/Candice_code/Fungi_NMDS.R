#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#NMDS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#====#
#Fungi
#====#
its.d<-read.csv(file="Final_Soil2_ITS_rarefied.csv")
its<-aggregate(.~provenance+season, data=its.d, mean)
levels(its$provenance)[levels(its$provenance)=="CTRL"] <- "NC"

#its<-its.d[its.d$season == "Fall",]
#perform nmds
sp.cols <- grep("SH", names(its))
md.all.s <- metaMDS(its[,sp.cols], distance="bray")
#save file as .tiff in 300 dpi
tiff("Fungi_NMDS.tiff", width = 120, height = 120, units = 'mm', res = 600)
#plot
all.plot.s<-plot(md.all.s, type="n", display="sites")
#assign grouping for color and shape
grp.sea<-factor(its$season)
grp.pr<-factor(its$provenance)
cols=c('magenta1', 'darkgreen','darkorange','red','blue')
shps=c(19,17)
#ordihull(md.all.s,groups=grp.pr,draw="polygon",col="white",label=F)
points(scores(all.plot.s), pch = shps[grp.sea], col = cols[grp.pr], cex=1.5)
ordiarrows(md.all.s, groups=grp.pr,levels=2, replicates=1, code=1, lty=1)
legend("topright", legend=levels(grp.sea), bty = "n", col= c("black"), pch = shps)
legend("bottomleft", legend = levels(grp.pr), bty = "n", col = cols, pch=c(20))
dev.off()
