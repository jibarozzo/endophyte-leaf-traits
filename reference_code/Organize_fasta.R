### Output fasta file has read No. after sample name so we cannot do OTU clustering and making OTU tables ###
### this script removes the read number after the sample name ###
### By Shuzo Oita. 20190807. Do not share with anyone else without Shuzo's permission ###

### 1. Copy the folder "filtered" and rename like "filtered2"
### 2. Remove "all_for_filtered" file from the new copied folder "filtered2".
### 3. If your fasta file has a name with something else, please run the line 30 after you changed the number of character you would like to keep
### 4. Run the code below. While running, you can see the files changed/overwritten one by one.
### 5. If you fasta file include a non-character symbol such as -, please remove them using text edit
### 6. Change directory of terminal to the new copied folder
### 7. Run the script of cat *.fa > merge.fa (you can change the name of "merge.fa" on the right side of the code)
### 8. Check if merge.fa file is okay (no read No. but has sample No.), and go back to the original directory using cd ..
### 9. Run the code of "Making OTU table" as the HTS protocol says


# install packages if not #
library(stringr)
library(dplyr)
library (tools)

### Get the list of the files with ".fa", and run a loop ###
files  <- list.files() 
for (file.name in files) {
  if (regexpr('\\.fa$', file.name)  < 0) { # Does the file has '.fa' at the bottom？
    next                                 # If not, skip．
  }
  x <- read.csv(file.name, header=F,stringsAsFactors = F) ##Read data
  y <- file_path_sans_ext(file.name)##Get the file name without format
  z <- paste(">", y, sep = "")
  # y2 <- substr(y,14,nchar(y))　# if the sample name has something else, run this code and remove them
  x <- if_else(grepl(">", x$V1), z, x$V1)　##replace the row with > with file name
  write(x, file=file.name)
}




