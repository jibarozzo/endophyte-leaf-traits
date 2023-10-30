library(Hmisc)

var<-read.csv(file="PlantTrait__correlation.csv")
na.omit(var)
#var<-var[complete.cases(var),]
var.corr<-rcorr(as.matrix(var), type="pearson")
var.corr
#extract only the correlation coefficient
r<-var.corr$r
write.table(r, file="trait_corr_coeff.csv", sep=',')

#too large to plot
pairs(var) 
#main="Simple Scatterplot Matrix") 


#clustering
#using hierarchical clustering using Ward.D
library(picante)
library(dendextend)

#correlation test among environmental variables
var<-read.csv(file="PlantTrait__correlation.csv")
trans<-t(var)

#cluster all variables
clus<-hclust(dist(trans), method="ward.D2")
plot(clus)

#cluster only the relevant variables, excluding chemical stuff
var<-read.csv(file="PlantTrait__correlation.csv")
var.2<-var[,3:23]
trans.2<-t(var.2)
clus<-hclust(dist(trans.2), method="ward.D2")
plot(clus)
