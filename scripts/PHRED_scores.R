
# Script to get PHRED scores:

# load libs
library(Rsubread)
library(tidyverse)
library(tidyr)
library(dplyr)
library(data.table)
library(splitstackshape)


# set directories
mydir <- "~/Desktop/MG1655/raw_libs"

########################################
########################################

reads = grep(list.files(mydir,pattern="(R1.fastq.gz|R2.fastq.gz)"), # pattern to match against
             pattern='interleaved', # pattern to ignore
             invert=TRUE, value=TRUE)

reads <- reads[39:62]

# construct an empty dataframe to build on 
df_to_fill <- data.frame(
  read_position = numeric(),
  PHRED_means = numeric(),
  library = character(),
  read = numeric(),
  stringsAsFactors = FALSE
)

for (read in reads) {
  

  # get quality scores from each file (all the reads)
  quality_scores <- qualityScores(file.path(mydir,read), input_format="gzFASTQ", offset=33, nreads = -1)
  
  # subset to framents < 300 bp long 
  quality_scores_300 <- quality_scores[,1:280]
  
  # collect the PHRED means 
  PHRED_means <- colMeans(quality_scores_300, na.rm = FALSE, dims = 1)
  PHRED_means <- as.data.frame(PHRED_means)
  
  # pick up lib name from loop and save in df
  PHRED_means$lib <- paste0(as.character(read))
  
  # some parsing
  PHRED_means <- cSplit(PHRED_means, "lib", sep = "_")
  colnames(PHRED_means) <- c("PHRED_means","library","read")
  PHRED_means$read <- gsub(".fastq.gz","",PHRED_means$read)
  
  # add read position 
  read_position <- seq(1,280)
  PHRED_means <- cbind(read_position, PHRED_means)
  
  # concatenate each new dataframe created 
  df_to_fill <- rbind(df_to_fill, PHRED_means)
  
}

# save the results 
fwrite(x = df_to_fill, file= file.path(mydir,"PHRED_scores.csv"))
