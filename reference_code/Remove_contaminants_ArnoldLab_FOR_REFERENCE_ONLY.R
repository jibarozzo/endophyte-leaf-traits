rm(list = ls()) #Remove all objects

###CC_fungi Contaminant removal. Adapted from MARIE and Shuzo oita
###THis is a good code. As of now, it works by substracting the SUM of contaminants from each OTU cell. If the value is negative it replaces it with 0. 

otu.cf <-read.csv('All_PCR_R1_otuFINAL_95.csv', as.is = T, row.names = 1)

##need OTUs on top for numeric cells/contamination removal
otu.cf.t <- as.data.frame(t(otu.cf))

#add column with sums for each OTU
cont <- row.names(otu.cf.t)
cont <- cont[9:11] #change accordingly to your data - these are the negative controls
contamination <- c()
for(c in 1:ncol(otu.cf.t)){
  contamination[c]<- sum(otu.cf.t[rownames(otu.cf.t) %in% cont, c], na.rm = TRUE)
}
otu.cf.t <-rbind(otu.cf.t, contamination)
row.names(otu.cf.t)[12] <- "contamination" #change the name of row 12

##subtract total contaminants from each OTU, if it is a negative number make it 0

cont2 <- c(cont, "contamination")
row <- which(!rownames(otu.cf.t) %in% cont2)
for(r in row){
  for(c in 1:ncol(otu.cf.t)){
    if(otu.cf.t[r,c] > otu.cf.t["contamination",c]) {
      new_reads <- otu.cf.t[r,c] - otu.cf.t["contamination",c]
      otu.cf.t[r,c] <- new_reads
    } else {otu.cf.t[r,c] <- 0}
  }
}

##remove controls from dataframe
otu.def <- otu.cf.t[!rownames(otu.cf.t) %in% cont2, ]


# # just to check, extract max and min values by column
# max <- apply(otu.def, 2, max)
# min <- apply(otu.def, 2,  FUN = function(x) {min(x[x > 0])})
# #and make a single dataframe
# otu.def <- rbind(otu.def,max)
# otu.def <- rbind(otu.def,min)

##remove all OTUs that have fewer than 10 ? reads 
otu.10 <- otu.def[,colSums(otu.def) > 10]

#export

write.csv(otu.10, "Data/CC_fungi_10_reads_99_forward_0.25_220_noContaminants.csv")

