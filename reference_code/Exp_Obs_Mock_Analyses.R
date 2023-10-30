##########################################################################
###Script for evaluating Mock Communities in Panama Leaf Trait project.###
##########################################################################
##Simple regression analyses of observed vs expected OTUs in mock communities
##Bolívar Aponte Rolón
##June 18, 2020


rm(list=ls())

###Potentially used packages throughout the script
library("vegan")
library("ggplot2")
library("ggfortify")
library("devtools")
library("tidyverse")
library("ggthemes")
library("ggiraphExtra")
library(lme4)
library(lattice)
library(nlme)

###Set WD
setwd("G:/My Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Sample_Sequencing/Post_Sequencing/Sequence_Analyses")

####################################
###Cleaning the data a little bit###
####################################



#######################################################################################################################
###Cleaning data for for quality control. Check blanks and negatives and look for anomalies.### 
###Average the sample replicates (e.g. E1+E2+E3/3) to find out the average read. The average read threshold is >=11.### 
###Eliminate OTUs that are below this threshold.Also consider the MOckBlan and MockNegative###
####Cleaning using Betsy Arnold's approach###
#######################################################################################################################
data<- read.csv("Exp_Ob_MockCommunitySequences_BAR_09142020.csv")
str(data)
as_tibble(data)
###This ^ data set has been created in Excel with the calculation sheet from Nicole Colón. 
###All the EXPECTED calculations were made without susbtracting reads from MockNegR1. Now i will recalculate that here.
options(scipen = 999) #This command eliminates scientific notation.

#dataclean <- data %>%
  #rename(OTU_ID = X.OTU.ID, EVEN1_OBS = E1_OBS, EVEN2_OBS = E2_OBS, EVEN3_OBS = E3_OBS, 
  #       TIER1_OBS = T1_OBS, TIER2_OBS = T2_OBS, TIER3_OBS = T3_OBS,
  #       EVEN1_EXP = E1_EXP, EVEN2_EXP = E2_EXP, EVEN3_EXP = E3_EXP,
  #       TIER1_EXP = T1_EXP, TIER2_EXP = T2_EXP, TIER3_EXP = T3_EXP,
  #       ) %>%
  #select(1, 10:12, 2:9, 14:25, -NOTES, -Ratio) %>%
  # rowwise() %>%
  # mutate(Tier_OBS_avg = mean(c(TIER1_OBS, TIER2_OBS, TIER3_OBS)), 
  #       Even_OBS_avg = mean(c(EVEN1_OBS, EVEN2_OBS, EVEN3_OBS)), 
  #       Blank_Neg_SUM = sum(c(MockBlank_R1, MockNeg_R1)))

#Dataclean 2
dataclean2 <- data %>%
  rename(OTU_ID = X.OTU.ID, EVEN1_OBS = E1_OBS, EVEN2_OBS = E2_OBS, EVEN3_OBS = E3_OBS, 
         TIER1_OBS = T1_OBS, TIER2_OBS = T2_OBS, TIER3_OBS = T3_OBS,
         EVEN1_EXP = E1_EXP, EVEN2_EXP = E2_EXP, EVEN3_EXP = E3_EXP,
         TIER1_EXP = T1_EXP, TIER2_EXP = T2_EXP, TIER3_EXP = T3_EXP) %>%
  select(1, 10:12, 2:9, 13:25, -NOTES, -MockBlank_R1) %>%
  rowwise() %>%
  mutate(TIER1_OBS = TIER1_OBS - MockNeg_R1, 
         TIER2_OBS = TIER2_OBS - MockNeg_R1,
         TIER3_OBS = TIER3_OBS - MockNeg_R1,
         EVEN1_OBS = EVEN1_OBS - MockNeg_R1,
         EVEN2_OBS = EVEN2_OBS - MockNeg_R1,
         EVEN3_OBS = EVEN3_OBS - MockNeg_R1) %>%
  mutate(TIER1_OBS = replace(TIER1_OBS, which(TIER1_OBS < 0), 0), 
         TIER2_OBS = replace(TIER2_OBS, which(TIER2_OBS < 0), 0),
         TIER3_OBS = replace(TIER2_OBS, which(TIER2_OBS < 0), 0),
         EVEN1_OBS = replace(EVEN1_OBS, which(EVEN1_OBS < 0), 0),
         EVEN2_OBS = replace(EVEN2_OBS, which(EVEN2_OBS < 0), 0),
         EVEN3_OBS = replace(EVEN3_OBS, which(EVEN3_OBS < 0), 0)) %>%
  mutate(EVEN1_EXP = Ratio * (sum(EVEN1_OBS)) , 
         EVEN2_EXP = Ratio * (sum(EVEN2_OBS)), 
         EVEN3_EXP = Ratio * (sum(EVEN3_OBS)),
         TIER1_EXP = Ratio * (sum(TIER1_OBS)),
         TIER2_EXP = Ratio * (sum(TIER2_OBS)), 
         TIER3_EXP = Ratio * (sum(TIER3_OBS))) %>%
  mutate(Tier_OBS_avg = mean(c(TIER1_OBS, TIER2_OBS, TIER3_OBS)), 
         Even_OBS_avg = mean(c(EVEN1_OBS, EVEN2_OBS, EVEN3_OBS)))

#Calculating Relative Abundances again from the decontaminated reads!
dataclean3 <- dataclean2 %>%
  group_by(OTU_ID)%>%
  summarise(T1_OBS_PCT = sum(TIER1_OBS),
            T2_OBS_PCT = sum(TIER2_OBS),
            T3_OBS_PCT = sum(TIER3_OBS),
            E1_OBS_PCT = sum(EVEN1_OBS),
            E2_OBS_PCT = sum(EVEN2_OBS),
            E3_OBS_PCT = sum(EVEN3_OBS),
            ) %>%
  mutate(T1_OBS_PCT = T1_OBS_PCT/sum(T1_OBS_PCT)*100,
         T2_OBS_PCT = T2_OBS_PCT/sum(T2_OBS_PCT)*100,
         T3_OBS_PCT = T3_OBS_PCT/sum(T3_OBS_PCT)*100,
         E1_OBS_PCT = E1_OBS_PCT/sum(E1_OBS_PCT)*100,
         E2_OBS_PCT = E2_OBS_PCT/sum(E2_OBS_PCT)*100,
         E3_OBS_PCT = E3_OBS_PCT/sum(E3_OBS_PCT)*100)

###Joining###
dataclean2 <- dataclean2 %>%
  select(-c(19:24)) %>%
  left_join(dataclean3, by = "OTU_ID")

###RAW data check point###
write.csv(dataclean2, file = "G:\\My Drive\\VBL_users\\Grad_Students\\Bolivar\\Dissertation\\Leaf_Traits_Panama\\Data\\Sample_Sequencing\\Post_Sequencing\\Sequence_Analyses\\Exp_Ob_Mock_Analyses_CLEAN_RAW_09152020.csv", row.names = FALSE)

####Introducing NAs to replace with new label.
#Dataclean
#dataclean[dataclean == "#N/A"] <- "XXX" #Changing "#N/A" to NA. Thiss needs improvement
#dataclean[dataclean == "na"] <- "XXX"
#dataclean$Phylum[dataclean$Phylum  == "0"] <-  "XXX"

#as.tibble(dataclean)

#dataclean <- dataclean %>%
#  mutate(
#    Code = as.character(Code),
#    Code = ifelse(is.na(Code), "UNKNOWN", Code),
#    Code = as.factor(Code),
#    Phylum = as.character(Phylum),
#    Phylum = ifelse(is.na(Phylum), "UNKNOWN", Phylum),
#    Phylum = as.factor(Phylum),
#    Species = as.character(Species),
#    Species = ifelse(is.na(Species), "UNKNOWN", Species),
#    Species = as.factor(Species))%>%
#  select(-MockBlank_R1) %>%
#  arrange(desc(Even_OBS_avg))%>%
# filter(Even_OBS_avg > 11)

#dataclean <- dataclean %>%
#  arrange(OTU_ID) %>%
#  filter(!OTU_ID %in% c("OTU_14")) 
###It's OTU14. Given that it is so prevalent in the Mock_Neg I think we should flag it and remove it. 
###Therefore we have 27 OTU remaining. 
###This leaves you with 27 OTU. In percentage terms, out of the sum of the reads for the 28 OTU, this removes everything < 0.02%.


####Introducing NAs to replace with new label.
#Dataclean2
dataclean2[dataclean2 == "#N/A"] <- "XXX" #Changing "#N/A" to NA. Thiss needs improvement
dataclean2[dataclean2 == "na"] <- "XXX"
dataclean2$Phylum[dataclean2$Phylum  == "0"] <-  "XXX"

as.tibble(dataclean2)

dataclean2 <- dataclean2 %>%
  mutate(
    Code = as.character(Code),
    Code = ifelse(is.na(Code), "UNKNOWN", Code),
    Code = as.factor(Code),
    Phylum = as.character(Phylum),
    Phylum = ifelse(is.na(Phylum), "UNKNOWN", Phylum),
    Phylum = as.factor(Phylum),
    Species = as.character(Species),
    Species = ifelse(is.na(Species), "UNKNOWN", Species),
    Species = as.factor(Species))

dataclean2 <- dataclean2 %>%
  arrange(desc(Even_OBS_avg))%>%
  filter(Even_OBS_avg > 11)



###Checking if Dataclean and datclean2 are the same###
dataclean_check <- anti_join(dataclean, dataclean2, by = "OTU_ID")

###Saving###
write.csv(dataclean2, file = "G:\\My Drive\\VBL_users\\Grad_Students\\Bolivar\\Dissertation\\Leaf_Traits_Panama\\Data\\Sample_Sequencing\\Post_Sequencing\\Sequence_Analyses\\Exp_Ob_Mock_Analyses_CLEAN_28OTU_09262020.csv", row.names = FALSE)

########################################################################################################################################################
###Dataclean2 incorporates Betsy Anrold's approach but first it is decontaminated. I substract the MockNeg_R1 reads from the Observed reads then###
###I proceded to apply the read approach from Betsy. I compared dataclean (just Betsy's approach without first substracting MockNeg_R1) and dataclean2 ###
###and found that it has the same OTU's as compare using semi_join and anti_join. All rows match in both data sets. Of course, dataclean2 has different ###
###expected read due to decontamination. ***I kept dataclean code just in case. 09/27/2020 -BAR
########################################################################################################################################################


#############################################################################################################
###Cleaning for reshaping the data to longer format. This makes it easier to plot and do statistical analyses.
#############################################################################################################
newdata <- dataclean2 %>%
  select(-Code, -Phylum, -Species, -MockNeg_R1, -Tier_OBS_avg, -Even_OBS_avg, -Ratio) %>%
  rename(E1_OBS = EVEN1_OBS, E2_OBS = EVEN2_OBS, E3_OBS = EVEN3_OBS, 
       T1_OBS = TIER1_OBS, T2_OBS = TIER2_OBS, T3_OBS = TIER3_OBS,
       E1_EXP = EVEN1_EXP, E2_EXP = EVEN2_EXP, E3_EXP = EVEN3_EXP,
       T1_EXP = TIER1_EXP, T2_EXP = TIER2_EXP, T3_EXP = TIER3_EXP) %>%
  pivot_longer(!OTU_ID, names_to = "ET", values_to = "Value")%>%
  separate(ET, c("ET", "OBS_EXP", "Relative_Abundance")) %>%
  pivot_wider(names_from = OBS_EXP:Relative_Abundance, values_from = Value) %>%
  rename(Even_Tier = ET, Observed = OBS_NA, Expected = EXP_NA, Relative_Abundance = OBS_PCT)%>%
  mutate(log10Expected = log10(Expected), log10Observed = log10(Observed))
  
 
#Replacing -inf in log10Expected column
newdata <- do.call(data.frame,  lapply(newdata,                    
                          function(x) replace(x, is.infinite(x), NA))) %>% # Replace Inf in data by NA
  replace_na(list(log10Expected = 0))


str(newdata)

#Checking how many distinct OTUs are present
n_distinct(newdata$OTU_ID)

#Saving data frame 
write.csv(newdata, file = "G:\\My Drive\\VBL_users\\Grad_Students\\Bolivar\\Dissertation\\Leaf_Traits_Panama\\Data\\Sample_Sequencing\\Post_Sequencing\\Sequence_Analyses\\Exp_Ob_Mock_Analyses_CLEAN_LONG_09152020.csv", row.names = FALSE)


#Histogram and distribution of read numbers
hist(newdata$log10Expected)
hist(newdata$log10Observed)
is.na(newdata$log10Expected)

qqnorm(rnorm(n = length(newdata$log10Expected), 
             mean =mean(newdata$log10Expected, na.rm = TRUE), 
             sd = sd(newdata$log10Expected, na.rm = TRUE)))
shapiro.test(rnorm(n = length(newdata$log10Expected), 
             mean =mean(newdata$log10Expected, na.rm = TRUE), 
             sd = sd(newdata$log10Expected, na.rm = TRUE)))

qqnorm(rnorm(n = length(newdata$log10Observed), 
             mean =mean(newdata$log10Observed, na.rm = TRUE), 
             sd = sd(newdata$log10Observed, na.rm = TRUE)))
shapiro.test(rnorm(n = length(newdata$log10Observed), 
                   mean =mean(newdata$log10Observed, na.rm = TRUE), 
                   sd = sd(newdata$log10Observed, na.rm = TRUE)))


###Subsetting for TIER community
newdata.T <- newdata %>%
  filter(Even_Tier %in% c("T1","T2", "T3"))

###Simple linear regressions###
newdata.lm <- lm(log10Expected ~ log10Observed,  data = newdata.T)
summary(newdata.lm)
AIC(newdata.lm)

newdata.gls <- gls(log10Expected ~ log10Observed, method = "REML",  data = newdata.T)
summary(newdata.lm)


newdata.rdm <- lme(log10Expected ~ log10Observed, random = ~1 | Even_Tier, method = "REML", data = newdata.T)
summary(newdata.rdm)

newdata.rmdwt <- lme(log10Expected ~ log10Observed, weights = varIdent(form =~ 1| Even_Tier), random = ~1 | Even_Tier, method = "REML", data = newdata.T)

anova(newdata.rmdwt, type = "marginal")
AIC(newdata.lm, newdata.gls, newdata.rdm, newdata.rmdwt)

###From this comparison of linear models, it looks like a simple linear regression with no randome ffects is the best fit model. 
#It has the lowest AIC value (59). All other models have more observations, thus higher AIC values. 


###GGPLOT###
ggplot(data = newdata.T, aes(y = Expected, x = Observed, color = Even_Tier, shape = Even_Tier))+
  geom_point(position = "jitter", size = 3, alpha = 0.5)+
  geom_smooth(method= "lm", se=F, formula = y ~ x, aes(fill = Even_Tier)) +
  scale_x_log10() +
  scale_y_log10() +
  expand_limits(y = 0, x = 0) +
  geom_abline(slope = newdata.lm$coefficients[2], intercept = newdata.lm$coefficients[1], size = 1) +
  annotate("text", size = 3.5, x = 125, y=0.1, label="y = 1.19*x - 2.17, Adj.R^2 = 0.87") +
  annotate("text", size = 3.5 , x = 10000, y = 0.1, label = "OTUs = 27") +
  labs(title = "Aim3 Mock Community without MockBlank1 ", subtitle = "Bolívar Aponte Rolón Sept 28, 2020", 
       y = "log10(Expected reads)", x = "log10(Observed reads)")

#ALTERNATIVE MODE OF PLOTTING. THIS REQUIRES THE USE OF MUTATE() IN THE PIPELINE TO CREATE
#TWO NEW COLUMNS WITH LOG10 VALUES. 
#ggplot(data = newdata, aes(y = log10_Expected, x = log10_Observed, color = EVEN_TIER, shape = EVEN_TIER))+
#    geom_point(size = 3)+
#    geom_jitter()+
#    geom_smooth(method= "lm", se=F, formula = y ~ x, aes(fill = EVEN_TIER)) +
#    geom_abline(slope = newdata.lm$coefficients[2], intercept = newdata.lm$coefficients[1], size = 1)
#  annotate("text", x = 2, y=4, label="y = -0.03*x + 3.45, R2 =")
  

###############################################################################################################################################

