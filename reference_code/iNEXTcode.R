#####   Hill 0, 1, 2 - iNEXT    ##### 
#columns are samples, rows are otus (species x site matrix)
#raw abundances
#calculates hill numbers of raw data (Observed) + hill numbers extrapolated to full sampling coverage (Estimated)
#using iNEXT package 

#extract otu data from ps1
otu2.raw <- data.frame(phyloseq::otu_table(ps2))
rowSums(otu2.raw) 
colSums(otu2.raw)
otu2.raw.t <- t(otu2.raw) 
otu2.raw.t <- round(otu2.raw.t, 0) #raw data was rounded to nearest integer! 
rowSums(otu2.raw.t)
colSums(otu2.raw.t)
extDiv2.all <- iNEXT(otu2.raw.t, q=c(0,1,2), size=NULL, endpoint=NULL, knots=30, datatype = "abundance", se=TRUE)
extDiv2.all$DataInfo
extDiv2.all$AsyEst

plot(extDiv2.all, type=1) #type=1 is a sample-size based rarefaction/extrapolation curve
#plots all three qs, cool.
plot(extDiv.all, type=2) #sample completeness curve

#saving data as dataframe - in format given by iNEXT
extDiv2.all.df <- extDiv2.all$AsyEst #Estimator tab is our extrapolated diversity, Observed is hill diversity of raw data


#reshaping to long format - observed diversity
Observed2_div.reshaped <- dcast(extDiv2.all.df, Assemblage~Diversity, value.var="Observed")

#changing names to be more meaningful
names(Observed2_div.reshaped)[names(Observed2_div.reshaped) == "Assemblage"] <- "Sample_names"
names(Observed2_div.reshaped)[names(Observed2_div.reshaped) == "Shannon diversity"] <- "q1.O"
names(Observed2_div.reshaped)[names(Observed2_div.reshaped) == "Simpson diversity"] <- "q2.O"
names(Observed2_div.reshaped)[names(Observed2_div.reshaped) == "Species richness"] <- "q0.O"


#reshaping to long format - estimated diversity
Estimated2_div.reshaped <- dcast(extDiv2.all.df, Assemblage~Diversity, value.var="Estimator")
#changing names to be more meaningful
names(Estimated2_div.reshaped)[names(Estimated2_div.reshaped) == "Assemblage"] <- "Sample_names"
names(Estimated2_div.reshaped)[names(Estimated2_div.reshaped) == "Shannon diversity"] <- "q1.E"
names(Estimated2_div.reshaped)[names(Estimated2_div.reshaped) == "Simpson diversity"] <- "q2.E"
names(Estimated2_div.reshaped)[names(Estimated2_div.reshaped) == "Species richness"] <- "q0.E"


#adding both observed and estimated to sampledata1
sampledata2 <- full_join(sampledata2,Observed2_div.reshaped, by="Sample_names")
sampledata2 <- full_join(sampledata2,Estimated2_div.reshaped, by="Sample_names")
