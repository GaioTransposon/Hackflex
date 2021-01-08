
# Script to get PHRED scores:

# load libs
library(Rsubread)
library(tidyverse)
library(tidyr)
library(dplyr)
library(data.table)


# set directories
mydir <- "~/Desktop/MG1655/raw_libs"

########################################
########################################

reads = grep(list.files(mydir,pattern="(R1.fastq.gz|R2.fastq.gz)"), # pattern to match against
             pattern='interleaved', # pattern to ignore
             invert=TRUE, value=TRUE)


# construct an empty dataframe to build on 
df_to_fill <- data.frame(
  read_position = numeric(),
  PHRED_means = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (read in reads) {
  
  print(read)

  # get quality scores from each file (all the reads)
  quality_scores <- qualityScores(file.path(mydir,read), input_format="gzFASTQ", offset=33, nreads = -1)
  
  # subset to framents < 300 bp long 
  quality_scores_300 <- quality_scores[,1:280]
  
  # collect the PHRED means 
  PHRED_means <- colMeans(quality_scores_300, na.rm = FALSE, dims = 1)
  
  PHRED_means <- as.data.frame(PHRED_means)
  PHRED_means$library <- paste0(as.character(read))
  
  read_position <- seq(1,280)
  PHRED_means <- cbind(read_position, PHRED_means)
  
  df_to_fill <- rbind(df_to_fill, PHRED_means)
  
}

# save the results 
fwrite(x = df_to_fill, file= file.path(mydir,"PHRED_scores.csv"))
