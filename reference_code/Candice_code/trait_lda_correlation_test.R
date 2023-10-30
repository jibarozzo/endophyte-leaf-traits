library(Hmisc)
library(MASS)
library(ggplot2)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Using the LDA values, obtained from predict(ld2a) functions
#soil data were the standardized values, 0 mean, sd =1
#Used in Appendix S6
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

ld.soil<-read.csv(file="SoilData_LDA.csv")

#make dummy variable
contrasts(ld.soil$Genotype) = c(1, 0.5, 0.25,0.25)
with(ld.soil, Genotype)

#as.factor(ld.soil$Genotype)

soilC.lda<-lm(soilC ~ Genotype+LD1+LD2+LD3, 
              data = ld.soil)
summary(soilC.lda)


soilN.lda<-lm(soilN ~ Genotype+LD1+LD2+LD3, 
              data = ld.soil)
summary(soilN.lda)


##########################################
#correlation of soil with plant provenance
##########################################
trait<-read.csv(file="PlantTrait_Data_original_noNAs.csv")

trait$logsoilC<-log(trait$soilC)
trait$logsoilN<-log(trait$soilN)
as.factor(trait$Genotype)

#transform and calculate some variables that Bri used
trait$CNag<-log(trait$AGC/trait$AGN)
trait$CNbg<-log(trait$BGC/trait$BGN)
trait$BGN.AGN<-log(trait$BGN/trait$AGN)
trait$BGC.AGC<-log(trait$BGC/trait$AGC)
#
soilC<-aov(logsoilC ~ Genotype, data = trait)
summary(soilC)

soilN<-aov(logsoilN ~ Genotype, data = trait)
summary(soilN)


#provenance with traits standardized values
trait.new<-read.csv(file="Standardized_plant_traits_combined_noCTRNL.csv")

as.factor(trait.new$Genotype)

soilC<-lm(soilC ~ Genotype+BGTotal+under30N+AGBiomass+diam+noLeaves+RtSt+over30N+
             lengthLeaf+HeightOver30+CNag+HeightUnder30+seedweight+BGN.AGN+BGC.AGC+CNbg, 
          data = trait.new)
summary(soilC)

soilN<-lme(soilN ~ Genotype, data = trait)
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
