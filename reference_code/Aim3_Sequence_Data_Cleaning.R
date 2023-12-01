#Script for cleaning and reshaping Sample and PCR blanks Sequence data for Leaf Traits Panama
#Bolívar Aponte Rolón
#February 16, 2022
#Last modified: November 30, 2023

rm(list = ls()) 
setwd("G:/My Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Sample_Sequencing/Post_Sequencing/Sequence_analyses")


###Potentially used packages throughout the script

library("ggplot2")
library("ggfortify")
library("devtools")
library("tidyverse")
library("ggthemes")
library("ggiraphExtra")

####################################################
###Cleaning and reshaping Blanks and PCR_CONTROLS###
####################################################

pcrdata <- read.csv("clean_data/post_biopipeline/All_PCR_R1_otuFINAL_95.csv") #PCR2NEGPLATE1 is missing. This is due to mislabeling in the original sequence files.
#There is a q.fastq.gz file that might be this sample.INSPECT THIS FILE. G:\My Drive\VBL_users\Grad_Students\Bolivar\Dissertation\Leaf_Traits_Panama\Data\Sample_Sequencing\Post_Sequencing\Original_Sample_gzip_files
#The issue stated in the above line has been resolved as of March 10, 2022. -BAR

pcrdata <- pcrdata %>%
  rename(OTU_ID = X.OTU.ID)

###Cleaning and reshaping Samples###
####################################

realsamp <- read.csv("clean_data/post_biopipeline/All_Samples_R1_otuFINAL95.csv")

realsamp <- realsamp %>%
  na.omit() %>%
  rename(OTU_ID = X.OTU.ID)

###############################################################################
###Full_join of the data frames to retain all values and row and transposing###
###############################################################################

#joined <- full_join(realsamp,pcrdata, by = "OTU_ID" ) %>%
# mutate(across(everything(), replace_na, replace = 0))

joined <- full_join(realsamp,pcrdata, by = "OTU_ID" ) %>%
  dplyr::mutate(across(c(2:168),~replace_na(.x, 0)))
str(joined)

joined.t <- as.data.frame(t(joined))
str(joined.t)

names(joined.t) <- joined.t %>% slice(1) %>% unlist()

joined.t <-  joined.t %>% 
  rownames_to_column("Sample_names") %>%
  slice(-1)

str(joined.t)

#####################################################
###Make subset tables from each sequence plate run###
#####################################################

################### Plate1 #####################

filter_table1<- read.csv("filters/Filter_Table_Plate1.csv")

names(filter_table1)
filter_table1 <- filter_table1 %>%
  rename(A = ?..A)
str(filter_table1)

filter_table1 <-filter_table1 %>% 
  pivot_longer(cols = 1:12, names_to = "samples", values_to = "values") %>%
  select(-samples) %>%
  slice(-5, -6)

plate1 <- semi_join(joined.t, filter_table1, by= c("Sample_names" = "values"))

plate1 <- plate1 %>%
  mutate(across(!Sample_names, as.numeric))
plate1 <- as_tibble(plate1)

filter_filter <- anti_join( filter_table1, plate1, by = c("values"="Sample_names"))#Corroborating if indeed all match. Supposed to be ZERO.
str(plate1)

#From factor dataframe to numeric

indx <- sapply(plate1, is.tibble)
plate1[indx] <- lapply(plate1[indx], function(x) as.numeric(as.character(x)))

################ Plate2 ####################

filter_table2<- read.csv("filters/Filter_Table_Plate2.csv")

names(filter_table2)

filter_table2 <- filter_table2 %>%
  rename(A = ?..A)

filter_table2 <- filter_table2 %>% 
  pivot_longer(cols = 1:10, names_to = "samples", values_to = "values") %>%
  select(-samples) %>%
  slice(-69, -70, -71, -77, -78, -79, -80)

plate2 <- semi_join(joined.t, filter_table2, by= c("Sample_names" = "values"))

plate2 <- plate2 %>%
  mutate(across(!Sample_names, as.numeric))
plate2 <- as_tibble(plate2)
filter_filter <- anti_join( filter_table2, plate2, by = c("values"="Sample_names")) #Corroborating if indeed all match. 
str(plate2)

#From factor dataframe to numeric

indx <- sapply(plate2, is.tibble)
plate2[indx] <- lapply(plate2[indx], function(x) as.numeric(as.character(x)))

###################################
###CC_fungi Contaminant removal.### 
###Adapted from MARIE and Shuzo Oita#####
###This code works by substracting the SUM of contaminants from each OTU cell. 
###If the value is negative it replaces it with 0. 

########MARELI MAGIC#########
##### Substracts the sum of the average of the ExtractionBlank avg and PCR negative control of Batch 001 from OTU sample reads, then removed control rows #####
######################################

################### Plate1 #####################

as.data.frame(plate1)
plate1 <- column_to_rownames(plate1, var = "Sample_names")

#add column with sums for each OTU

cont <- row.names(plate1)
cont <- cont[92:94] #change accordingly to your data - these are the negative controls
contamination <- c()
for(c in 1:ncol(plate1)){
  contamination[c]<- mean(plate1[rownames(plate1) %in% cont, c], na.rm = TRUE)
}
plate1 <-rbind(plate1, contamination)
row.names(plate1)[95] <- "contamination" #change the name of row 104

###subtract total contaminants from each OTU, if it is a negative number make it 0
cont2 <- c(cont, "contamination")
row <- which(!rownames(plate1) %in% cont2)
for(r in row){
  for(c in 1:ncol(plate1)){
    if(plate1[r,c] > plate1["contamination",c]) {
      new_reads <- plate1[r,c] - plate1["contamination",c]
      plate1[r,c] <- new_reads
    } else {plate1[r,c] <- 0}
  }
}


##remove controls from dataframe and makes a text file to paste into existing excel sheet

new_data_decontaminated001 <- plate1[!rownames(plate1) %in% cont2, ]
write.table(new_data_decontaminated001,"clean_data/otu_data/MinusNegsBlanks_001.txt", sep="\t") ### This became tab "Filtered-BlankNegbyCode" in "ExtractionBatch0001_PCR102NC.xlsx". 


################### Plate2 #####################

as.data.frame(plate2)
plate2 <- column_to_rownames(plate2, var = "Sample_names")

#add column with sums for each OTU

cont <- row.names(plate2)
cont <- cont[66:73] #change accordingly to your data - these are the negative controls
contamination <- c()
for(c in 1:ncol(plate2)){
  contamination[c]<- mean(plate2[rownames(plate2) %in% cont, c], na.rm = TRUE)
}
plate2 <-rbind(plate2, contamination)
row.names(plate2)[74] <- "contamination" #change the name of row 104

###subtract total contaminants from each OTU, if it is a negative number make it 0

cont2 <- c(cont, "contamination")
row <- which(!rownames(plate2) %in% cont2)
for(r in row){
  for(c in 1:ncol(plate2)){
    if(plate2[r,c] > plate2["contamination",c]) {
      new_reads <- plate2[r,c] - plate2["contamination",c]
      plate2[r,c] <- new_reads
    } else {plate2[r,c] <- 0}
  }
}



##remove controls from dataframe and makes a text file to paste into exisitng excel sheet

new_data_decontaminated002 <- plate2[!rownames(plate2) %in% cont2, ]
write.table(new_data_decontaminated002,"clean_data/otu_data/MinusNegsBlanks_002.txt", sep="\t")


###################################
####Saving decontaminated file ####
###################################
all_complete_decontaminated <- bind_rows(new_data_decontaminated001, new_data_decontaminated002)

write.csv(all_complete_decontaminated,"clean_data/otu_data/all_Decontaminated_byBlankandNegs_complete.csv")
trans <- t(all_complete_decontaminated)
trans <- as.data.frame(trans)
write.csv(trans, "clean_data/otu_data/all_toFilt_10_percent.csv") ## in the format to remove <10% OTUs -- see below



##########################################################################
#### Script to remove < 0.10% abundance per sample #### From Shuzo Oita###
##########################################################################
## csv file should have row = OTU, col = sample

otu.data <- read.csv('clean_data/otu_data/all_toFilt_10_percent.csv', row.names = 1,as.is = T)
otu.data2 <- apply(otu.data, 2, function(x) ifelse({100*x/sum(x)} < 0.10, 0, x))
OTU_table.cleaned <- subset((otu.data2), rowSums(otu.data2)>1)
OTU_table.cleaned <- as.data.frame(OTU_table.cleaned)
write.csv(OTU_table.cleaned, "all_cleaned_at_10_percent.csv")



###################
#Code for creating .csv files out of .txt files


#setwd("")
#filelist = list.files(pattern = ".txt")
#for (i in 1:length(filelist)){
#  input<-filelist[i]
#  output <- paste0(gsub("\\.txt$", "", input), ".csv")
#  print(paste("Processing the file:", input))
#  data = read.delim(input, header = TRUE)   
#  setwd("")
#  write.table(data, file=output, sep=",", col.names=TRUE, row.names=FALSE)
#  setwd("")
#}

#"G:/My Drive/VBL_users/Grad_Students/Bolivar/Dissertation/Leaf_Traits_Panama/Data/Sample_Sequencing/Post_Sequencing/Sequence_Analyses/Sample_Analyses/Panama_NGS_Sample_Analyses.R"
