#### Betapart -- Beta diversity over time 
##Code for Boli adapted from Ch1 - Mareli's dissertation. 


library(betapart)

#making pse objects into otu tables
otue2.hei <- data.frame(phyloseq::otu_table(pse2HEI))
otue2.gus <- data.frame(phyloseq::otu_table(pse2GUS))
otue2.cal <- data.frame(phyloseq::otu_table(pse2CAL))
otue2.vir <- data.frame(phyloseq::otu_table(pse2VIR))
otue2.hyb <- data.frame(phyloseq::otu_table(pse2HYB))
otue2.cor <- data.frame(phyloseq::otu_table(pse2COR))

# Run this for each species 
### HEISCO 
pair2.HEI <- beta.pair.abund(otu2.hei, index.family = "bray") 
#output is 3 distance matrices (balance, gradient and total turnvoer). Output is a matrix of pairwise dissimilarity values among individuals for all time point combinations. 
#since I was only interested in dissimilarity values among individuals at T1-T2, T2-T3, and T1-T3, and not for all combinations, I manually pulled out those comparisons of interest.
pair2.HEI

#extracting each matrix from function result
HEI2.p.bal <- as.matrix(pair2.HEI$beta.bray.bal)
HEI2.p.gra <- as.matrix(pair2.HEI$beta.bray.gra)
HEI2.p.beta <- as.matrix(pair2.HEI$beta.bray)

#exporting dist matrices from pair.HEI to manually extract pairwise comparisons. 
write.csv(HEI2.p.bal, "Exp2_betapart/HEI2.p.bal_beta.pair.Jan2023.csv")
write.csv(HEI2.p.gra, "Exp2_betapart/HEI2.p.gra_beta.pair.Jan2023.csv")
write.csv(HEI2.p.beta, "Exp2_betapart/HEI2.p.beta_beta.pair.Jan2023.csv")

#reuploading the compiled table. 
HEI2_beta_combo <- read_excel("Exp2_betapart/HEI2.p.beta_beta.pair.Jan2023_SUMMARIZED_Jan2023.xlsx", sheet="HEI2_summary_combo")

### do for all species.... ###

###### rbinding to make beta_combo -- main turnover dataset
beta2_combo <- rbind(HEI2_beta_combo,CAL2_beta_combo, GUS2_beta_combo, VIR2_beta_combo, HYB2_beta_combo, CORD2_beta_combo)
beta2_combo$Lifespan_days <- as.character(beta2_combo$Lifespan_days)
beta2_combo$TimePoint <- as.character(beta2_combo$TimePoint)
View(beta2_combo)
str(beta2_combo)
saveRDS(beta2_combo, "./beta2_combo.rds")


#### Spearman's ro correlations #### Assessing how the trajectory of change (slope) changes over time in each species for each turnover fraction

#### gradient
turn.short <- cor.test(beta2.short$Days_since_Inoc, beta2.short$Gradient_Turnover, method="spearman",exact=FALSE)
turn.short #no, negative

turn.long <- cor.test(beta2.long$Days_since_Inoc, beta2.long$Gradient_Turnover, method="spearman",exact=FALSE)
turn.long #yes, negative

#Turnover by species
gra2.HEI <- cor.test(beta2_combo.HEI.4TP$Days_since_Inoc, beta2_combo.HEI.4TP$Gradient_Turnover, method="spearman", exact=FALSE)
gra2.HEI #no
gra2.CAL <- cor.test(beta2_combo.CAL.4TP$Days_since_Inoc, beta2_combo.CAL.4TP$Gradient_Turnover, method="spearman", exact=FALSE)
gra2.CAL #yes, negative

gra2.GUS <- cor.test(beta2_combo.GUS.4TP$Days_since_Inoc, beta2_combo.GUS.4TP$Gradient_Turnover, method="spearman", exact=FALSE)
gra.GUS #no

gra2.VIR <- cor.test(beta2_combo.VIR.4TP$Days_since_Inoc, beta2_combo.VIR.4TP$Gradient_Turnover, method="spearman", exact=FALSE)
gra2.VIR #yes, negative

gra2.COR<- cor.test(beta2_combo.COR$Days_since_Inoc, beta2_combo.COR$Gradient_Turnover, method="spearman", exact=FALSE)
gra2.COR #no

gra2.HYB<- cor.test(beta2_combo.HYB$Days_since_Inoc, beta2_combo.HYB$Gradient_Turnover, method="spearman", exact=FALSE)
gra2.HYB #no

#balancing  - same as gradient but with balancing turnover variable 

#total turnover - same but with total turnover variable


####### How does each turnover fraction differ among species? 
#boxplots and non parametric anova 
beta_all.complete$Species <- as.factor(beta_all.complete$Species)

gra.box <- ggplot(beta_all.complete, aes(x=Species, y=Gradient_Turnover)) + geom_boxplot(show.legend = TRUE) + geom_point(alpha=0.1, size=3,aes(color=Species), show.legend = TRUE, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.1))  + scale_x_discrete(limits=c("HEISCO", "GUSTSU", "CALOLO", "VIROSU", "HYBAPR", "CORDAL")) + theme_bw() + ylab("Gradient Turnover ") + xlab("Species")
gra.box
kruskal.test(Gradient_Turnover~Species, data=beta_all.complete)
pairwise.wilcox.test(beta_all.complete$Gradient_Turnover, beta_all.complete$Species, p.adjust.method = "bonferroni")


bal.box <- ggplot(beta_all.complete, aes(x=Species, y=Balancing_Turnover)) + geom_boxplot(show.legend = TRUE) + geom_point(alpha=0.1, size=3,aes(color=Species), show.legend = TRUE, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.1))  + scale_x_discrete(limits=c("HEISCO", "GUSTSU", "CALOLO", "VIROSU", "HYBAPR", "CORDAL")) + theme_bw() + ylab("Balancing Turnover ") + xlab("Species")
bal.box
kruskal.test(Balancing_Turnover~Species, data=beta_all.complete)
pairwise.wilcox.test(beta_all.complete$Balancing_Turnover, beta_all.complete$Species, p.adjust.method = "bonferroni")

turn.box <- ggplot(beta_all.complete, aes(x=Species, y=Turnover)) + geom_boxplot(show.legend = TRUE) + geom_point(alpha=0.1, size=5,aes(color=Species), show.legend = TRUE, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0.1)) + scale_x_discrete(limits=c("HEISCO", "GUSTSU", "CALOLO", "VIROSU", "HYBAPR", "CORDAL")) + theme_bw() + ylab("Turnover ") + xlab("Species")
turn.box
kruskal.test(Balancing_Turnover~Species, data=beta_all.complete)
pairwise.wilcox.test(beta_all.complete$Turnover, beta_all.complete$Species, p.adjust.method = "bonferroni")
