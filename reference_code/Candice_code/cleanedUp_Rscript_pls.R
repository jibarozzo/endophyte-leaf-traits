library(pls)
library(ggplot2)

#++++++++++++++++++++++++++++++++++#
#Partial Least regression Analysis
#Microbial diversity (alpha)
#++++++++++++++++++++++++++++++++++#
#skip below and proceed staight to standardized/scaled traits
fundiv<-read.csv(file="FinalData_Soil2_ITS_rarefied_noNAs.csv")

#first scale the PAHs values, plant traits have already been scaled
pah<-fundiv[,37]
fundiv$scaled.pah<-scale(pah)

#merged scaled plant traits with ITS data above

traits<-read.csv(file="Standardized_plant_traits_combined_noCTRNL.csv")
merged.d<-merge(fundiv, traits, by="SampleID", keep_y=TRUE)
write.table(merged.d, file="FungalITS_finalData_scaledVariables_pls.csv", sep=',', row.names=FALSE)


#=========================================#
#Conduct PLS using scaled variables
#=========================================#
#===========#
#1. FUNGI
#===========#

fundiv<-read.csv(file="FungalITS_finalData_scaledVariables_pls.csv")

#create dummy variable for season since this is regression
fundiv$Season<-NA
fundiv$Season[fundiv$season=="SUMMER"] <-0
fundiv$Season[fundiv$season=="WINTER"] <-1

#provenance
fundiv$Provenance<-NA
fundiv$Provenance[fundiv$provenance=="BJ"] <-0
fundiv$Provenance[fundiv$provenance=="CL"] <-1
fundiv$Provenance[fundiv$provenance=="CP"] <-2
fundiv$Provenance[fundiv$provenance=="V"] <-3

#rename some variables
names(fundiv)[names(fundiv) =="under30N"]<-"shoot_density"
names(fundiv)[names(fundiv) =="over30N"]<-"tiller_density"
names(fundiv)[names(fundiv) =="leaf_length"]<-"leaf_length"
names(fundiv)[names(fundiv) =="HeightOver30"]<-"tiller_height"
names(fundiv)[names(fundiv) =="HeightUnder30"]<-"shoot_height"

#perform plsr
fun.div.pls<-plsr(shan~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                    tiller_density+tiller_height+shoot_height+shoot_density+
                    seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                  data=fundiv, validation="CV", scale=TRUE)

# Find the number of dimensions with lowest cross validation error
cv<-RMSEP(fun.div.pls)
best.dims<-which.min(cv$val[estimate = "adjCV", , ]) - 1
best.dims

# Rerun the model, best.dims=3
new.fundiv.plsr<-plsr(shan~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                    tiller_density+tiller_height+shoot_height+shoot_density+
                    seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                    data=fundiv, ncomp = 3,validation="CV", scale=TRUE,jackknife = TRUE)
summary(new.fundiv.plsr)
jack.test(new.fundiv.plsr,ncomp=1:3, use.mean = TRUE)

#extract coefficient
coefficients<-coef(new.fundiv.plsr, ncomp=3)
coefficients

#extract explained variance attributed to components
compnames(new.fundiv.plsr, comps = 1:2, explvar = TRUE)

factor=(row.names(coefficients))
coef<-data.frame(factor)
coef

# change name into numbers
levels(coef$factor)[levels(coef$factor) =="BG"]<-"1"
levels(coef$factor)[levels(coef$factor) =="AG"]<-"2"
levels(coef$factor)[levels(coef$factor) =="shoot_height"]<-"3"
levels(coef$factor)[levels(coef$factor) =="shoot_density"]<-"4"
levels(coef$factor)[levels(coef$factor) =="tiller_density"]<-"5"
levels(coef$factor)[levels(coef$factor) =="tiller_height"]<-"6"
levels(coef$factor)[levels(coef$factor) =="leaf_length"]<-"7"
levels(coef$factor)[levels(coef$factor) =="leaf_num"]<-"8"
levels(coef$factor)[levels(coef$factor) =="RtSt"]<-"9"
levels(coef$factor)[levels(coef$factor) =="diameter"]<-"10"
levels(coef$factor)[levels(coef$factor) =="seed_mass"]<-"11"
levels(coef$factor)[levels(coef$factor) =="CNag"]<-"12"
levels(coef$factor)[levels(coef$factor) =="CNbg"]<-"13"
levels(coef$factor)[levels(coef$factor) =="BGN.AGN"]<-"14"
levels(coef$factor)[levels(coef$factor) =="BGC.AGC"]<-"15"
levels(coef$factor)[levels(coef$factor) =="soilC"]<-"16"
levels(coef$factor)[levels(coef$factor) =="soilN"]<-"17"
levels(coef$factor)[levels(coef$factor) =="PAH"]<-"18"
coef


#plot correlation loadings 
tiff("PLS_Shandiv_fungi_number_c.tiff", width = 100, height = 100, units = 'mm', res = 700)
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.65)
#mtext("D", outer=TRUE)
par(new=TRUE)
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=17, col="darkred")
dev.off()


#^^^^^^^^^^^^^^^^^^#
#1.2 Chao1 fungi
#^^^^^^^^^^^^^^^^^^#
#perform plsr
fun.chao.pls<-plsr(S.chao1~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                    tiller_density+tiller_height+shoot_height+shoot_density+
                    seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                  data=fundiv, validation="CV", scale=TRUE,jackknife = TRUE)

# Find the number of dimensions with lowest cross validation error
cv<-RMSEP(fun.chao.pls)
best.dims<-which.min(cv$val[estimate = "adjCV", , ]) - 1
best.dims

# Rerun the model, best.dims=3
new.funchao.plsr<-plsr(S.chao1~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                       tiller_density+tiller_height+shoot_height+shoot_density+
                       seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                     data=fundiv, ncomp = 3,validation="CV", scale=TRUE,jackknife = TRUE)
summary(new.funchao.plsr)

jack.test(new.funchao.plsr,ncomp=3, use.mean = TRUE)

#extract coefficient
coefficients.chaofun<-coef(new.funchao.plsr, ncomp=3)
coefficients.chaofun


#extract explained variance attributed to components
compnames(new.funchao.plsr, comps = 1:2, explvar = TRUE)

factor=(row.names(coefficients.chaofun))
coef<-data.frame(factor)
coef

# change name into numbers
levels(coef$factor)[levels(coef$factor) =="BG"]<-"1"
levels(coef$factor)[levels(coef$factor) =="AG"]<-"2"
levels(coef$factor)[levels(coef$factor) =="shoot_height"]<-"3"
levels(coef$factor)[levels(coef$factor) =="shoot_density"]<-"4"
levels(coef$factor)[levels(coef$factor) =="tiller_density"]<-"5"
levels(coef$factor)[levels(coef$factor) =="tiller_height"]<-"6"
levels(coef$factor)[levels(coef$factor) =="leaf_length"]<-"7"
levels(coef$factor)[levels(coef$factor) =="leaf_num"]<-"8"
levels(coef$factor)[levels(coef$factor) =="RtSt"]<-"9"
levels(coef$factor)[levels(coef$factor) =="diameter"]<-"10"
levels(coef$factor)[levels(coef$factor) =="seed_mass"]<-"11"
levels(coef$factor)[levels(coef$factor) =="CNag"]<-"12"
levels(coef$factor)[levels(coef$factor) =="CNbg"]<-"13"
levels(coef$factor)[levels(coef$factor) =="BGN.AGN"]<-"14"
levels(coef$factor)[levels(coef$factor) =="BGC.AGC"]<-"15"
levels(coef$factor)[levels(coef$factor) =="soilC"]<-"16"
levels(coef$factor)[levels(coef$factor) =="soilN"]<-"17"
levels(coef$factor)[levels(coef$factor) =="PAH"]<-"18"
coef


#plot correlation loadings 
tiff("PLS_chao1_fungi_final_b.tiff", width = 110, height = 110, units = 'mm', res = 300)
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.7)
par(new=TRUE)
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.7)
dev.off()

#=============#
#2. BACTERIA
#=============#

bactdiv<-read.csv(file="Final_Soil2_16S_rarefied_noNAS.csv")
#first scale the PAHs values, plant traits have already been scaled
pah<-bactdiv[,19]
bactdiv$scaled.pah<-scale(pah)

#merged scaled plant traits with ITS data above

traits<-read.csv(file="Standardized_plant_traits_combined_noCTRNL.csv")
merged.d<-merge(bactdiv, traits, by="SampleID")
write.table(merged.d, file="Bacteria_16S_finalData_scaledVariables_pls.csv", sep=',', row.names=FALSE)

#=========================================#
#Conduct PLS using scaled variables
#=========================================#
bactdiv<-read.csv(file="Bacteria_16S_finalData_scaledVariables_pls.csv")

#create dummy variable for season since this is regression
bactdiv$Season<-NA
bactdiv$Season[bactdiv$season=="SUMMER"] <-0
bactdiv$Season[bactdiv$season=="WINTER"] <-1

#provenance
bactdiv$Provenance<-NA
bactdiv$Provenance[bactdiv$provenance=="BJ"] <-0
bactdiv$Provenance[bactdiv$provenance=="CL"] <-1
bactdiv$Provenance[bactdiv$provenance=="CP"] <-2
bactdiv$Provenance[bactdiv$provenance=="V"] <-3


#rename some variables
names(bactdiv)[names(bactdiv) =="under30N"]<-"shoot_density"
names(bactdiv)[names(bactdiv) =="over30N"]<-"tiller_density"
names(bactdiv)[names(bactdiv) =="lengthLeaf"]<-"leaf_length"
names(bactdiv)[names(bactdiv) =="HeightOver30"]<-"tiller_height"
names(bactdiv)[names(bactdiv) =="HeightUnder30"]<-"shoot_height"

#perform plsr
bact.div.pls<-plsr(shan~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                    tiller_density+tiller_height+shoot_height+shoot_density+
                    seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                   data=bactdiv, validation="CV", scale=TRUE,jackknife = TRUE)
summary(bact.div.pls, what = "validation")
# Find the number of dimensions with lowest cross validation error
cv<-RMSEP(bact.div.pls)
best.dims<-which.min(cv$val[estimate = "adjCV", , ]) - 1
best.dims
# Rerun the model
new.bactdiv.plsr<-plsr(shan~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                        tiller_density+tiller_height+shoot_height+shoot_density+
                        seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                       data=bactdiv, ncomp = 2,validation="CV", scale=TRUE,jackknife = TRUE)
summary(new.bactdiv.plsr)
jack.test(new.bactdiv.plsr,ncomp=2, use.mean = TRUE)

#bar plot of coefficient
coefficients.bac<-coef(new.bactdiv.plsr)
coefficients.bac
factor=(row.names(coefficients.bac))
coef<-data.frame(factor)
coef
# change name into numbers
levels(coef$factor)[levels(coef$factor) =="BG"]<-"1"
levels(coef$factor)[levels(coef$factor) =="AG"]<-"2"
levels(coef$factor)[levels(coef$factor) =="shoot_height"]<-"3"
levels(coef$factor)[levels(coef$factor) =="shoot_density"]<-"4"
levels(coef$factor)[levels(coef$factor) =="tiller_density"]<-"5"
levels(coef$factor)[levels(coef$factor) =="tiller_height"]<-"6"
levels(coef$factor)[levels(coef$factor) =="leaf_length"]<-"7"
levels(coef$factor)[levels(coef$factor) =="leaf_num"]<-"8"
levels(coef$factor)[levels(coef$factor) =="RtSt"]<-"9"
levels(coef$factor)[levels(coef$factor) =="diameter"]<-"10"
levels(coef$factor)[levels(coef$factor) =="seed_mass"]<-"11"
levels(coef$factor)[levels(coef$factor) =="CNag"]<-"12"
levels(coef$factor)[levels(coef$factor) =="CNbg"]<-"13"
levels(coef$factor)[levels(coef$factor) =="BGN.AGN"]<-"14"
levels(coef$factor)[levels(coef$factor) =="BGC.AGC"]<-"15"
levels(coef$factor)[levels(coef$factor) =="soilC"]<-"16"
levels(coef$factor)[levels(coef$factor) =="soilN"]<-"17"
levels(coef$factor)[levels(coef$factor) =="PAH"]<-"18"
coef


#plot correlation loadings 
tiff("PLS_Shandiv_bacteria_final_numbers_c.tiff", width = 100, height = 100, units = 'mm', res = 700)
plot(new.bactdiv.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.70)
par(new=TRUE)
plot(new.bactdiv.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=17, col="darkred")
dev.off()

#based on scores
plot(scores(new.bactdiv.plsr))

#based on coef
#tiff("PLS_Shandiv_bacteria_coef.tiff", width = 170, height = 170, units = 'mm', res = 700)
#plot(new.bactdiv.plsr, plottype = "coef",ncomp=1:2,
#     labels=row.names(coefficients.bac), cex=0.9)
#dev.off()

#------------#
#2.2 Chao1
#------------#
bactdiv$chao1<-log(bactdiv$S.chao1)
bact.divchao.pls<-plsr(chao1~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                     tiller_density+tiller_height+shoot_height+shoot_density+
                     seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                   data=bactdiv, validation="CV", scale=TRUE,jackknife = TRUE)

cv<-RMSEP(bact.divchao.pls)
best.dims<-which.min(cv$val[estimate = "adjCV", , ]) - 1
best.dims

new.bactdiv.chao.plsr<-plsr(chao1~Season+Provenance+PAH+BG+diameter+AG+leaf_num+RtSt+leaf_length+
                         tiller_density+tiller_height+shoot_height+shoot_density+
                         seed_mass+CNag+CNbg+BGN.AGN+BGC.AGC+soilC+soilN,
                         data=bactdiv, ncomp = 2,validation="CV", scale=TRUE,jackknife = TRUE)
summary(new.bactdiv.chao.plsr)
jack.test(new.bactdiv.chao.plsr,ncomp=2, use.mean = TRUE)

#bar plot of coefficient
coefficients.bac.chao<-coef(new.bactdiv.chao.plsr)
coefficients.bac.chao
factor=(row.names(coefficients.bac.chao))
coef<-data.frame(factor)
coef
# change name into numbers
levels(coef$factor)[levels(coef$factor) =="BG"]<-"1"
levels(coef$factor)[levels(coef$factor) =="AG"]<-"2"
levels(coef$factor)[levels(coef$factor) =="shoot_height"]<-"3"
levels(coef$factor)[levels(coef$factor) =="shoot_density"]<-"4"
levels(coef$factor)[levels(coef$factor) =="tiller_density"]<-"5"
levels(coef$factor)[levels(coef$factor) =="tiller_height"]<-"6"
levels(coef$factor)[levels(coef$factor) =="leaf_length"]<-"7"
levels(coef$factor)[levels(coef$factor) =="leaf_num"]<-"8"
levels(coef$factor)[levels(coef$factor) =="RtSt"]<-"9"
levels(coef$factor)[levels(coef$factor) =="diameter"]<-"10"
levels(coef$factor)[levels(coef$factor) =="seed_mass"]<-"11"
levels(coef$factor)[levels(coef$factor) =="CNag"]<-"12"
levels(coef$factor)[levels(coef$factor) =="CNbg"]<-"13"
levels(coef$factor)[levels(coef$factor) =="BGN.AGN"]<-"14"
levels(coef$factor)[levels(coef$factor) =="BGC.AGC"]<-"15"
levels(coef$factor)[levels(coef$factor) =="soilC"]<-"16"
levels(coef$factor)[levels(coef$factor) =="soilN"]<-"17"
levels(coef$factor)[levels(coef$factor) =="PAH"]<-"18"
coef

#plot
tiff("PLS_chao1_bacteria_final_b.tiff", width = 110, height = 110, units = 'mm', res = 300)
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.7)
par(new=TRUE)
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.80)
dev.off()

#************************#
#Plot all four plsr plots
#Shannon and Chao
#************************#

tiff("FigS7_Combined_plsr_4_divplots.tiff", width = 140, height = 140, units = 'mm', res = 400)

#plot all four plsr
par(mfrow=c(2,2))
par(mar=c(4,4,2,1))

#1. Bacteria Shan
plot(new.bactdiv.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.70)
par(new=TRUE)
plot(new.bactdiv.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=17, col="darkred")

#2. Fungi
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.65)
par(new=TRUE)
plot(new.fundiv.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=17, col="darkred")

#3. Bacteria Chao
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.7)
par(new=TRUE)
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.80)

#4. Fungi Chao
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.7)
par(new=TRUE)
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.7)

dev.off()

#Combined chao
tiff("Combined_chao_plsr.tiff", width = 150, height = 90, units = 'mm', res = 300)
par(mfrow=c(1,2))
par(mar=c(2,3,1,4))
#1. Bacteria
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.65)
par(new=TRUE)
plot(new.bactdiv.chao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.80)

#Fungi
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=FALSE, labels=coef$factor, cex=0.65)
par(new=TRUE)
plot(new.funchao.plsr, plottype = "correlation", 
     ploty=TRUE, plotx=FALSE, pch=9, col="blue", cex=0.7)

dev.off()
