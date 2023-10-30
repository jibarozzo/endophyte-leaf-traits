library(Hmisc)
library(MASS)
library(ggplot2)
library(plyr)

#****************************************************#
#====================================================#
#LDA analysis with standardized estimates
#This was the final analysis included in the manuscript
#=====================================================#
#****************************************************#

#First scale the data with mean = 0 and sd = 1
trait<-read.csv(file="PlantTrait_Data_original.csv")

#trait.d2<-trait.d[!trait.d$Genotype=="CTRL",]
#remove rows with NAs
#trait.d[complete.cases(trait.d), ]

#transform and calculate some variables that Bri used
trait$CNag<-scale(trait$AGC/trait$AGN)
trait$CNbg<-scale(trait$BGC/trait$BGN)
trait$BGN.AGN<-scale(trait$BGN/trait$AGN)
trait$BGC.AGC<-scale(trait$BGC/trait$AGC)

trt.std<-trait[,4:40]
scaled.trait<-scale(trt.std)

#get column sample/genotype info
gent<-trait[,1:3]
#combine with scaled trait values
newd.trait<-cbind(gent,scaled.trait)

#write.table(newd.trait, file="Standardized_plant_traits_combined.csv", sep=',', row.names=FALSE)

#CNTRL was removed manually. This can also be done in R
#I also removed the provenance with NA data in relevant traits

#trait.newd2<-trait.newd[!trait.newd$Genotype=="CTRL",]


#-------------#
#-------------#
#LDA analysis
#no soil data
#-------------#
#-------------#

trait.newd<-read.csv(file="Standardized_plant_traits_combined_noCTRNL.csv")

trait.newd$Genotype<-as.factor(trait.newd$Genotype)

#first change the variable names
names(trait.newd)[names(trait.newd) =="BGTotal"]<-"BG"
names(trait.newd)[names(trait.newd) =="AGBiomass"]<-"AG"
names(trait.newd)[names(trait.newd) =="under30N"]<-"shoot_density"
names(trait.newd)[names(trait.newd) =="diam"]<-"shoot_diameter"
names(trait.newd)[names(trait.newd) =="noLeaves"]<-"leaf_number"
names(trait.newd)[names(trait.newd) =="over30N"]<-"tiller_density"
names(trait.newd)[names(trait.newd) =="lengthLeaf"]<-"leaf_length"
names(trait.newd)[names(trait.newd) =="HeightOver30"]<-"tiller_height"
names(trait.newd)[names(trait.newd) =="HeightUnder30"]<-"shoot_height"
names(trait.newd)[names(trait.newd) =="seedweight"]<-"seedmass"

ld2a <- lda(Genotype ~ BG+shoot_density+AG+shoot_diameter +leaf_number + RtSt+
              tiller_density+leaf_length +tiller_height +CNag+
              shoot_height+ seedmass + BGC.AGC+BGN.AGN+CNbg, 
            data = trait.newd, na.action=na.omit, method="moment")

ld2a
#plot(ld2a)
ld2b<-data.frame(varnames=rownames(coef(ld2a)), coef(ld2a))

#traitsLD123 <-as.data.frame(ld2a$scaling)
#traitsLD123

#first change the variable names
names(ld2a)[names(ld2a) =="seedmass"]<-paste(expression('seed mass'))
names(ld2a)[names(ld2a) =="CNag"]<-paste(expression('C:N'[AG]))
names(ld2a)[names(ld2a) =="shoot_height"]<-paste(expression('shoot height'))
#names(trait.newd)[names(trait.newd) =="BGN.AGN"]<-"N[BG:AG]"
#names(trait.newd)[names(trait.newd) =="BGC.AGC"]<-"C[BG:AG]"
#names(trait.newd)[names(trait.newd) =="CNbg"]<-"C:N[BG]"


lda.values <- predict(ld2a)
lda.values

#scatterplot LD1 vs LD2
#plot(lda.values$x[,1], lda.values$x[,2], bty='n')
#LD1 vs LD3
#plot(lda.values$x[,1], lda.values$x[,3], bty='n')


head(lda.values)
#only the x variables from predict
lda123<-as.data.frame(lda.values$x)
row.names(lda123)

#write as table the predicted values for each sample
LDtraits<-cbind(trait.newd,lda123)
LDtraits
#write.table(LDtraits, file="Standardized_PlantTraits_LDs_noNAs.csv", sep=',', row.names=FALSE)
#first change the variable names
#names(traitsLD123)[names(traitsLD123) =="seedmass"]<-paste(expression('seed mass'))
#names(trait.newd)[names(trait.newd) =="CNag"]<-paste(expression('C:N'[AG]))
#names(trait.newd)[names(trait.newd) =="shoot_height"]<-paste(expression('shoot height'))
#names(trait.newd)[names(trait.newd) =="seedmass"]<-paste(expression('seed mass'))

row.names(ld2a$scaling)
#labels<-data.frame(factor)
#labels
# change name into numbers
levels(labels$factor)[levels(labels$factor) =="BG"]<-"1"
levels(labels$factor)[levels(labels$factor) =="AG"]<-"2"
labels

#+++++++++++++++++++++++++++++++++++++++#
#Figure 1
#plot the data for Figure 1 in the paper
#R code adapted from Bri's
#+++++++++++++++++++++++++++++++++++++++#

#LD1 vs LD2
#add the arrow function
lda.arrows <- function(x, myscale =0.8 , tex = 0.70, choices = c(1,2), ...){
  ## adds `biplot` arrows to an lda using the discriminant function values
  heads <- coef(x)
  labels.1<-factor
  arrows(x0 = 0, y0 = 0, 
         x1 = myscale * heads[,choices[1]], 
         y1 = myscale * heads[,choices[2]], length=unit(0.02, "npc"),...)#angle=30
  text(myscale * heads[,choices], labels = labels.1, 
       cex = tex)
}

lda12=as.data.frame(lda.values$x)[,1:2]
cols=c('magenta1', 'darkgreen','darkorange','blue')
#tiff("Fig1_top_nosoil.tiff", res=400, width=150, height=110,unit="mm")
plot(lda12,type="n",xlim=c(-7,7),ylim=c(-8,4), bty='n',xlab("LD1 (71%)"), ylab("LD2 (24%)"))
points(x=lda12[,1],y=lda12[,2],col=cols[trait.newd$Genotype],pch=15,
       bg=trait.newd$Genotype)
lda.arrows(ld2a, col = 2, myscale = 0.9)
dev.off()

factor
ggbiplot(ld2a, choices = 1:2, scale = 1, pc.biplot = TRUE)
         
#LD1 vs LD3
lda13=as.data.frame(lda.values$x)[,1:3]
lda13$LD2=NULL
lda.arrows <- function(x, myscale =1 , tex = 0.70, choices = c(1,3), ...){
  ## adds `biplot` arrows to an lda using the discriminant function values
  heads <- coef(x)
  arrows(x0 = 0, y0 = 0, 
         x1 = myscale * heads[,choices[1]], 
         y1 = myscale * heads[,choices[2]], length=.07,...)#angle=30
  text(myscale * heads[,choices], labels = row.names(heads), 
       cex = tex)
}
cols=c('magenta1', 'darkgreen','darkorange','blue')
tiff("Fig1_bottom_b_nosoil.tiff", res=400, width=150, height=110,unit="mm")
plot(lda13,type="n",ylim=c(-3,4),xlim=c(-7, 8),bty='n',xlab("LD1 (71%)"), ylab("LD3 (4%)"))
points(x=lda13[,1],y=lda13[,2],col=cols[trait.newd$Genotype],pch=15,
       bg=trait.newd$Genotype)
lda.arrows(ld2a, col = 2, myscale = 1)
dev.off()


#+++++++++++++++++++++++#
#OR
#+++++++++++++++++++++++#

library(devtools)
install_github('fawda123/ggord') # Used to install ggord from github we need to run devtools to achieve this.
library(ggord)




#----------------------#
#check correlations
#----------------------#
#test for correlation between original and new discriminat function
#the following are based on Bri's R coding
#Some I get, some I totally don't

#traits are in columns 6-20 in the trait.newd file, without the soil
corrLDA=cor(trait.newd[6:20],lda123)
write.table(corrLDA, file="Correlation_coefficient_nosoil.csv", sep=',')

#cor(cbind(trait.newd$BGTotal,trait.newd$under30N,trait.newd$AGBiomass,trait.newd$diam,trait.newd$noLeaves,
#          trait.newd$RtSt,trait.newd$over30N,
#          trait.newd$lengthLeaf,trait.newd$HeightOver30,trait.newd$CNag,trait.newd$HeightUnder30,
#          trait.newd$soilC,trait.newd$soilN,
#          trait.newd$seedweight,trait.newd$BGN.AGN,trait.newd$BGC.AGC,trait.newd$CNbg, lda.values))

#determine which is significant
#there's got to be a better way to code this
#cor.test(trait.newd[,4],lda123[,1])
#cor.test(trait.newd[,5],lda123[,1])
cor.test(trait.newd[,6],lda123[,1])
cor.test(trait.newd[,7],lda123[,1])
cor.test(trait.newd[,8],lda123[,1])
cor.test(trait.newd[,9],lda123[,1])
cor.test(trait.newd[,10],lda123[,1])
cor.test(trait.newd[,11],lda123[,1])
cor.test(trait.newd[,12],lda123[,1])
cor.test(trait.newd[,13],lda123[,1])
cor.test(trait.newd[,14],lda123[,1])
cor.test(trait.newd[,15],lda123[,1])
cor.test(trait.newd[,16],lda123[,1])
cor.test(trait.newd[,17],lda123[,1])
cor.test(trait.newd[,18],lda123[,1])
cor.test(trait.newd[,19],lda123[,1])
cor.test(trait.newd[,20],lda123[,1])
#cor.test(trait.newd[,4],lda123[,2])
#cor.test(trait.newd[,5],lda123[,2])
cor.test(trait.newd[,6],lda123[,2])
cor.test(trait.newd[,7],lda123[,2])
cor.test(trait.newd[,8],lda123[,2])
cor.test(trait.newd[,9],lda123[,2])
cor.test(trait.newd[,10],lda123[,2])
cor.test(trait.newd[,11],lda123[,2])
cor.test(trait.newd[,12],lda123[,2])
cor.test(trait.newd[,13],lda123[,2])
cor.test(trait.newd[,14],lda123[,2])
cor.test(trait.newd[,15],lda123[,2])
cor.test(trait.newd[,16],lda123[,2])
cor.test(trait.newd[,17],lda123[,2])
cor.test(trait.newd[,18],lda123[,2])
cor.test(trait.newd[,19],lda123[,2])
cor.test(trait.newd[,20],lda123[,2])
#cor.test(trait.newd[,4],lda123[,3])
#cor.test(trait.newd[,5],lda123[,3])
cor.test(trait.newd[,6],lda123[,3])
cor.test(trait.newd[,7],lda123[,3])
cor.test(trait.newd[,8],lda123[,3])
cor.test(trait.newd[,9],lda123[,3])
cor.test(trait.newd[,10],lda123[,3])
cor.test(trait.newd[,11],lda123[,3])
cor.test(trait.newd[,12],lda123[,3])
cor.test(trait.newd[,13],lda123[,3])
cor.test(trait.newd[,14],lda123[,3])
cor.test(trait.newd[,15],lda123[,3])
cor.test(trait.newd[,16],lda123[,3])
cor.test(trait.newd[,17],lda123[,3])
cor.test(trait.newd[,18],lda123[,3])
cor.test(trait.newd[,19],lda123[,3])
cor.test(trait.newd[,20],lda123[,3])
#ldahist(data = lda.values$x[,1], g=trait.newd$Genotype)

##########################################
#correlation of soil with plant provenance
##########################################
trait<-read.csv(file="PlantTrait_Data_original.csv")

#remove rows with NAs
trait()
trait[complete.cases(trait), ]
trait$logsoilC<-log(trait$soilC)
trait$logsoilN<-log(trait$soilN)
as.factor(trait$Genotype)

#
soilC<-aov(logsoilC ~ Genotype, data = trait)
summary(soilC)

soilN<-aov(logsoilN ~ Genotype, data = trait)
summary(soilN)

#=========================================#
#Correlation soil c&N with traits
#=========================================#
trait<-read.csv(file="PlantTrait_Data_original.csv")
var<-trait[,5:27]
var.corr<-rcorr(as.matrix(var), type="pearson")
var.corr

#extract only the correlation coefficient
r<-var.corr$r
write.table(r, file="trait_and_soil_corr_coeff.csv", sep=',')

#extract only the P value
p<-var.corr$P
p.dist<-as.matrix(p, upper=FALSE)
write.table(p.dist, file="trait_and_soil_p-value_dist.csv", sep=',')

write.table(p, file="trait_and_soil_p-value_correlation.csv", sep=',')

#####################################
#LDA plant traits, untransformed
#####################################
library(MASS)
library(ggpubr)
#library(tidyverse)
#library(klaR)

trait<-read.csv(file="PlantTrait_Data_original.csv")

#transform and calculate some variables that Bri used
trait$CNag<-log(trait$AGC/trait$AGN)
trait$CNbg<-log(trait$BGC/trait$BGN)
trait$BGN.AGN<-log(trait$BGN/trait$AGN)
trait$BGC.AGC<-log(trait$BGC/trait$AGC)

#transform variables
trait$AGBsqrt<-1/sqrt(trait$AGBiomass)
trait$logAGB<-log(trait$AGBiomass)

#visualize normality
ggdensity(trait$logAGB)
ggdensity(trait$AGBsqrt)

#transform variables
trait$AGBsqrt<-1/sqrt(trait$AGBiomass)
trait$logAGB<-log(trait$AGBiomass)

trait$logBGTotal<-log(trait$BGTotal)
trait$logunder30N<-log(trait$under30N)
trait$logdiam<-log(trait$diam)
trait$lognoLeaves<-log(trait$noLeaves)

trait$loglengthLeaf<-log(trait$lengthLeaf)
trait$logover30N<-log(trait$over30N)
trait$logHeightOver30<-log(trait$HeightOver30)
trait$logRtSt<-log(trait$RtSt)

trait$logseedweight<-log(trait$seedweight)
trait$logHeightUnder30<-log(trait$HeightUnder30)

#trait$logBGTotal<-log(trait$BGTotal)
#trait$logunder30N<-log(trait$under30N)
#trait$logdiam<-log(trait$diam)
#trait$lognoLeaves<-log(trait$noLeaves)

#trait$loglengthLeaf<-log(trait$lengthLeaf)
#trait$logover30N<-log(trait$over30N)
#trait$logHeightOver30<-log(trait$HeightOver30)
#trait$logRtSt<-log(trait$RtSt)

#trait$logseedweight<-log(trait$seedweight)
#trait$logHeightUnder30<-log(trait$HeightUnder30)

#untransformed
lda <- lda(Genotype ~ BGTotal+under30N+AGBiomass+diam+noLeaves+RtSt+over30N+
             lengthLeaf+HeightOver30+CNag+HeightUnder30+soilC+soilN+
             seedweight+BGN.AGN+BGC.AGC+CNbg, data = trait, na.action=na.omit, method="moment")
lda


#save as different file name
write.table(trait, file="PlantTrait_Data_transformed_all.csv", sep=',', row.names=FALSE)

#============================#
#Using ade4 discrimin function
#============================#
library(ade4)
fsum.d<-read.csv(file="Fungi_summer_transformed_noNAs.csv")
#fsum.d<-fsum.d[complete.cases(fsum.d),]
fsum.d$provenance=as.factor(fsum.d$provenance)
pca1 <- dudi.pca(fsum.d[,19:28], scannf = FALSE)
dis1 <- discrimin(pca1, fsum.d$provenance, scannf = FALSE)
dis1
