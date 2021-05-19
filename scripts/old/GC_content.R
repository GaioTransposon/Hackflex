
# Script to get GC content plotted: 

# load libs
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggpubr)
library(weights)


# set directories
middle_dir <- "~/Desktop/MG1655/main/"
out_dir <- "~/Desktop/MG1655/main/"

########################################
########################################

gc_files = list.files(middle_dir,pattern="GC_")
flagstat_files = list.files(middle_dir,pattern="flagstat")

myDF <- data.frame(gc_files,flagstat_files, stringsAsFactors = F)

datalist = list()
for (row in 1:nrow(myDF)) {
  
  gc <- myDF[row,1]
  flag <- myDF[row,2]
  
  id <- str_replace_all(gc, "GC_qc_reduced_trimmed2_trimmed_interleaved2_", "")
  
  # read in files
  gc_df <- read.table(file.path(middle_dir,gc), quote="\"", comment.char="", header = TRUE)
  flag_df <- read.delim(file.path(middle_dir,flag), header=FALSE)
  mapped_min_dups <- as.numeric(str_extract(flag_df[5,], "[^+]+"))
  
  gc_df$fractionOfReads <- (gc_df$fractionOfReads + 0.000000001)
  gc_df1 <- spread(gc_df, GCcontent, fractionOfReads)
  
  rownames(gc_df1) <- gc_df1$Sample
  gc_df1$Library <- NULL
  gc_df2 <- as.data.frame(t(gc_df1[,-1]))
  colnames(gc_df2) <- gc_df1$Sample
  
  gc_df2$diff <- gc_df2[,1]/gc_df2[,2]
  
  # mapped_min_dups is the tot. number of mapped reads (as per samtools flagstat) minus the eliminated duplicates
  gc_df2$reads <- gc_df2[,1]* mapped_min_dups
  gc_df3 <- cbind(gc_df2, gc_df[,3, drop=FALSE])
  
  rho <- wtd.cor(gc_df3$GCcontent,gc_df3$diff,weight=gc_df3$reads)
  rho <- rho[1,1]
  
  gc_df3$lib <- id 
  gc_df3$rho <- rho
  
  colnames(gc_df3) <- c("deduped_bam","Reference","diff","reads","GCcontent","lib","rho")
  
  datalist[[row]] <- gc_df3 # add it to your list
  
}

big_data = do.call(rbind, datalist)

#plot
pdf(paste0(out_dir,'GC_content.pdf'))
plot(as.data.frame(datalist[1])$GCcontent[as.data.frame(datalist[1])$diff!=1], as.data.frame(datalist[1])$diff[as.data.frame(datalist[1])$diff!=1], 
     type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(as.data.frame(datalist[2])$GCcontent[as.data.frame(datalist[2])$diff!=1], as.data.frame(datalist[2])$diff[as.data.frame(datalist[2])$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(as.data.frame(datalist[3])$GCcontent[as.data.frame(datalist[3])$diff!=1], as.data.frame(datalist[3])$diff[as.data.frame(datalist[3])$diff!=1], col="green", cex = 0.5, pch = 1) 
abline(lm(as.data.frame(datalist[1])$diff ~ as.data.frame(datalist[1])$GCcontent, weights=as.data.frame(datalist[1])$reads), col="blue")
abline(lm(as.data.frame(datalist[2])$diff ~ as.data.frame(datalist[2])$GCcontent, weights=as.data.frame(datalist[2])$reads), col="red")
abline(lm(as.data.frame(datalist[3])$diff ~ as.data.frame(datalist[3])$GCcontent, weights=as.data.frame(datalist[3])$reads), col="green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", "...")),
                  expression(paste("1:50 Flex           ", rho, " = ", "...")),
                  expression(paste("Hackflex            ", rho, " = ", "..."))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)

big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=lib))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  geom_smooth() +
  ylab("ratio observed/expected reads") +
  theme(legend.title = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)


big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=lib))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  stat_smooth(method="lm", se=FALSE) +
  ylab("ratio observed/expected reads") +
  theme(legend.title = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)

plot(as.data.frame(datalist[1])$GCcontent, as.data.frame(datalist[1])$diff, 
     type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(as.data.frame(datalist[2])$GCcontent, as.data.frame(datalist[2])$diff, col="red", cex = 0.5, pch = 2) +
  points(as.data.frame(datalist[3])$GCcontent, as.data.frame(datalist[3])$diff, col="green", cex = 0.5, pch = 1) 
abline(lm(as.data.frame(datalist[1])$diff ~ as.data.frame(datalist[1])$GCcontent, weights=as.data.frame(datalist[1])$reads), col="blue")
abline(lm(as.data.frame(datalist[2])$diff ~ as.data.frame(datalist[2])$GCcontent, weights=as.data.frame(datalist[2])$reads), col="red")
abline(lm(as.data.frame(datalist[3])$diff ~ as.data.frame(datalist[3])$GCcontent, weights=as.data.frame(datalist[3])$reads), col="green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", "...")),
                  expression(paste("1:50 Flex           ", rho, " = ", "...")),
                  expression(paste("Hackflex            ", rho, " = ", "..."))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()

