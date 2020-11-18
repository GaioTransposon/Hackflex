
# Script to get Fragment sizes plotted:

# load libs
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# set directories
middle_dir <- "~/Desktop/MG1655/hackflex_libs/ecoli/"
out_dir <- "~/Desktop/MG1655/hackflex_libs/ecoli/"

########################################
########################################

frag_files = list.files(middle_dir,pattern=".txt")

# construct an empty dataframe to build on 
df_to_fill <- data.frame(
  V1 = numeric(),
  order = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (frag_file in frag_files) {
  
  # read in file 
  frag_df <- read.table(file.path(middle_dir,frag_file), quote="\"", comment.char="")
  
  #subset to positive numbers only
  pos_df <- subset(frag_df, frag_df$V1 > 0)
  
  #add column with sequential numbers to each dataset
  new_df <- cbind( pos_df, order=seq(nrow(pos_df)) ) 
  
  new_df$library <- as.character(basename(frag_file))
  
  df_to_fill <- rbind(
    df_to_fill, 
    new_df
  )
  
}

head(df_to_fill)
tail(df_to_fill)


#replace library names 
df_to_fill$library <- str_replace_all(df_to_fill$library, "frag_sizes_nanop_vs_K1.txt", "Standard Flex EPM")
df_to_fill$library <- str_replace_all(df_to_fill$library, "frag_sizes_nanop_vs_K2.txt", "Standard Flex PS")
df_to_fill$library <- str_replace_all(df_to_fill$library, "frag_sizes_nanop_vs_K3.txt", "Flex EPM 1:50")
df_to_fill$library <- str_replace_all(df_to_fill$library, "frag_sizes_nanop_vs_K9.txt", "Hackflex")

colnames(df_to_fill) <- c("value","order","library")

head(df_to_fill)

#plot
pdf(paste0(out_dir,'fragment_size_new.pdf'))
libs <- c("Standard Flex EPM", 
          "Standard Flex PS", 
          "Flex EPM 1:50", 
          "Hackflex")
ds <- subset(df_to_fill, library %in% libs & ! is.na(value))
p  <- ggplot(ds, aes(value, colour=library, fill=library)) + scale_x_continuous(limits=c(0,1000)) + xlab("fragment size (bp)")
q  <- p + geom_density(alpha=0.55)
r <- q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
s <- r + theme(legend.title = element_blank()) 
t <- s + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14))
t
dev.off()


###############################################################
###############################################################