
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

reads = grep(list.files(mydir,pattern="(.fastq.gz)"), # pattern to match against
             pattern='interleaved', # pattern to ignore
             invert=TRUE, value=TRUE)

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
  
  # collect the PHRED means 
  PHRED_means <- colMeans(quality_scores, na.rm = FALSE, dims = 1)
  PHRED_means <- as.data.frame(PHRED_means)
  
  # add read position
  PHRED_means$read_position <- seq(1,NROW(PHRED_means))
  
  # save name of library
  read <- gsub(".fastq.gz","",read)
  PHRED_means$library <- paste0(read)
  
  # save PHRED means
  fwrite(PHRED_means, file = file.path(mydir,paste0(read,".csv")))
}


