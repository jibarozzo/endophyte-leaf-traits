#########################################################################################################

#Aim 2 - removing Eminus OTUs
#Created on: April 2, 2022 
#Created by: Mareli

#What it do?
#It removes all OTUs that are only found in E- plants, as they can be considered greenhouse only endophytes.Data already has CORDBI samples removed from OTU table and samples table

#libraries
library(readxl)
library(dplyr)
library(tibble)
library(phyloseq)
library(ggpubr)


#Description/Samples table
#matches AW2 OTU table rownames
#run this first because OTU rownames are set to descript rownames
descript_exp2 <- read_excel("Aim2_data/Samples_Otu_Tax_Exp2.xlsx", sheet="SamplesData_noBI")
row.names(descript_exp2) <- descript_exp2$Sample_names
View(descript_exp2)
descript_table2 <- sample_data(descript_exp2)
View(sample_data(descript_table2))
row.names(descript_table2) <- descript_exp2$Sample_names
View(sample_data(descript_table2))
#this gives you column names
names(descript_exp2)

# TAX TABLE 
tax2 <- read_excel("Aim2_data/Samples_Otu_Tax_Exp2.xlsx", sheet="Tax")
View(tax2)
row.names(tax2) <- tax2$OTU_ID
tax2 <- as.matrix(tax2)
class(tax2) 
tax_table2 <- tax_table(tax2)
View(tax_table(tax_table2))

#OTU - AW2
otu.AW2.2 <- read_excel("Aim2_data/Samples_Otu_Tax_Exp2.xlsx", sheet="Out_AW2_noBI")
otu.AW2.2 <- otu.AW2.2 %>% select(-Sample_names)
otu.AW2.2 <- as.matrix(otu.AW2.2)
class(otu.AW2.2) <- "numeric"
row.names(otu.AW2.2) <- descript_table2$Sample_names
otu_table.AW2.2 <- otu_table(otu.AW2.2, taxa_are_rows = FALSE)
View(otu_table(otu_table.AW2.2))

#make phyloseq object
ps2 <- phyloseq(otu_table.AW2.2, descript_table2, tax_table2)
ps2

#any OTUs with no reads?
ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2) # ALOT. 1139 OTUs in Experiment 2 only. 
ps2

# subset of ps2 -- only E+ samples and otus. Important -- dont remove OTUs with 0 reads here, because we want to know which ones have no reads in E+ but DO have reads in E-. 
ps2_e <- subset_samples(ps2, Eload == "E+")
ps2_e #1139 taxa

ps2_e_EminusRemoved <- prune_taxa(taxa_sums(ps2_e) > 0, ps2_e) #removed 406 OTUs, so thats how many I should expect to be removed when done manually. But im not doing it manually. 

#if OTU has 0 reads in Eplus, it will have reads in Eminus, so this prunes ps2_e and identifies all otus that have 0 reads in E+, which by default, are the OTUs that DO have reads in E-. 
ps2_e_EminuOnly <- prune_taxa(taxa_sums(ps2_e) == 0, ps2_e) 
ps2_e_EminuOnly

#gets names of OTUs and make them in a data frame with column name the same as in ps2
EminusOTUs2 <- taxa_names(ps2_e_EminuOnly)
View(EminusOTUs2)
EminusOTUsList2 <- as.data.frame(EminusOTUs2)
class(EminusOTUsList2)
View(EminusOTUsList2)
write.table(EminusOTUsList2, "Take3_Aim2_RemovedOTUsList.txt")


#Get vectors (names) of numbered OTUs from ps2
OTU_ID <- colnames(otu_table(ps2))
View(OTU_ID)


## returns all otu names that dont have a match in Eminus OTUs. 
notShared <- setdiff(OTU_ID, EminusOTUs2)
notShared

## Subset phyloseq object to only these OTUs
trimmedOTU <- subset(otu_table(ps2), select = colnames(otu_table(ps2)) %in% notShared)
View(trimmedOTU)

#exported into "Samples_Otu_Tax_Exp2.xlsx" as "OTU_AW2_trimmed"
write.csv(trimmedOTU, "Take3_Aim2_OTU_AW2_trimmed.csv")

#Outputs of this script -- OTU table without E minus and updated samples table without cordbi and with a column of OTU sample richness that matches the OTU table is in "Take3_Aim1.2Data_Final_Jan2022.xslx". This will be used in the final analysis. 


