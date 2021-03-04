

# load libs
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)
library(weights)
library(readr)
library(splitstackshape)
library(stringr)
library(corrr)
library(data.table)
library(stringr)



# set directories and select samples: 

mydir <- "~/Desktop/MG1655/goal_barcode/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"


########################################
########################################


# barcode distribution: 

barcode_count <- read.table(file.path(mydir,"all_fastq_headers_clean.tsv"), quote="\"", comment.char="", header = FALSE)
head(barcode_count)




barcodes <- barcode_count %>%
  dplyr::select(V1,V3) 



head(barcodes)


# take the original barcode (this is the most frequent one)
barcodes1 <- barcodes %>%
  drop.levels() %>%
  dplyr::mutate(n=n()) %>%
  group_by(V1,V3) %>%
  dplyr::summarise(barcode_sum=sum(n)) %>%
  group_by(V1) %>%
  slice_max(order_by = barcode_sum,n=1) # change to n>1 if you want the top n most frequent barcode sequences

head(barcodes1)

barcodes1 <- cSplit(barcodes1, "V3", sep = ":")
barcodes1 <- cSplit(barcodes1, "V3_4", sep = "+")

barcodes2 <- barcodes1 %>%
  dplyr::select(V1, barcode_sum, V3_1, V3_4_1, V3_4_2) %>%
  dplyr::rename(library = V1,
                barcode_count = barcode_sum,
                read_direction = V3_1,
                i5 = V3_4_1,
                i7 = V3_4_2)

head(barcodes2)

barcodes3 <- barcodes2 %>%
  dplyr::mutate(G_i5 = str_count(i5, "G"),
                C_i5 = str_count(i5, "C"),
                G_i7 = str_count(i7, "G"),
                C_i7 = str_count(i7, "C"),
                GC_i5 = ((G_i5 + C_i5) / str_length(i5) * 100),
                GC_i7 = ((G_i7 + C_i7) / str_length(i7) * 100)) %>%
  dplyr::select(library, barcode_count,read_direction,GC_i5,GC_i7) %>%
  pivot_longer(cols=c("GC_i5","GC_i7"), names_to="oligo")
barcodes3$oligo <- gsub("GC_","", barcodes3$oligo)
head(barcodes3)
NROW(barcodes3)




# PLOT: 

pdf(paste0(mydir,"barcodes.pdf"))
# barcode distribution: 
ggplot(barcodes3, aes(x=barcode_count)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
# barcode counts vs barcode GC content
ggplot(barcodes3, aes(x=value,y=barcode_count, color=oligo))+
  geom_point()+
  stat_smooth(method="lm", se=FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(x="barcode GC content (%)",
       y="barcode count")
# barcodes of each library, their count, and their GC content
barcodes3 %>%
  mutate(library = fct_reorder(library, value)) %>%
  dplyr::filter(read_direction==1) %>%
  ggplot(., aes(x=library,y=barcode_count, color=value))+ 
  geom_point(alpha=0.8)+
  facet_grid(rows = vars(oligo)) +
  theme(axis.text.x=element_blank()) +
  scale_color_gradientn(colours = rainbow(5))
dev.off()






# need to get read count per lib! and plot this: 
# Plot number of reads assigned to barcode vs 96 barcode (rank sorted)
a <- barcodes3 %>% dplyr::filter(read_direction==1)
barcode_count_2 <- a$barcode_count
barcode_count_3 <- data.frame(sort(barcode_count_2))
#add column for wells
barcode_count_3$X2 <- 1:nrow(barcode_count_3) 
barplot(barcode_count_3$sort.barcode_count_2~barcode_count_3$X2, xlab="", ylab="", 
        xaxt='n', main="96-plex barcode counts") 
title(xlab="96 barcodes (rank-sorted)", line=1.0, cex.lab=1.5)
title(ylab="number of reads assigned to barcode", line=2.5, cex.lab=1.5)
grid(nx = NA, ny = 10, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)






# coeff of variation without the two outliers (percentage)
# subset the dataframe to exclude outliers with read count lower than 50
barcode_count_no_outliers <- subset(barcode_count, V2>50)

a = paste0("Barcode count summary:")
b = paste0("Min: ", summary(barcode_count$V2)[1])
c = paste0("1st Qu.: ", summary(barcode_count$V2)[2])
d = paste0("Median: ", summary(barcode_count$V2)[3])
e = paste0("Mean: ", summary(barcode_count$V2)[4])
f = paste0("3rd Qu.: ", summary(barcode_count$V2)[5])
g = paste0("Max: ", summary(barcode_count$V2)[6])
h = paste0("Standard deviation: ", sd(barcode_count$V2))
i = paste0("Sum: ", sum(barcode_count$V2))
l = paste0("Coeff. of variation (incl. outliers): ", sd(barcode_count$V2)/mean(barcode_count$V2))
m = paste0("Coeff. of variation (excl. outliers): ", sd(barcode_count_no_outliers$V2)/mean(barcode_count_no_outliers$V2))

pdf(paste0(mydir,'barcode_summary.pdf'),width=11,height=11)
plot(NA, xlim=c(0,11), ylim=c(0,11), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,11,a, pos=4)
text(1,10,b, pos=4)
text(1,9,c, pos=4)
text(1,8,d, pos=4)
text(1,7,e, pos=4)
text(1,6,f, pos=4)
text(1,5,g, pos=4)
text(1,4,h, pos=4)
text(1,3,i, pos=4)
text(1,2,l, pos=4)
text(1,1,m, pos=4)
dev.off()



########################################

# Barcode GC content:




#regress.lm = lm(barcode_count$X2 ~ barcodes_GC$X4, weights=barcode_count$X2)
#summary(regress.lm)
NROW(barcodes)
head(barcodes)
extreme_GC <- barcodes %>% dplyr::filter(GC_content1<20|GC_content1>80)
NROW(extreme_GC)/NROW(barcodes)*100
#it means that 10.75% of all the barcodes have more than 80% or less than 20% GC content

#pearson coefficient of correlation
barcodes_GC_pearson1 <- wtd.cor(barcodes$GC_content1,barcodes$barcode_count,weight=barcodes$barcode_count)
barcodes_GC_pearson2 <- wtd.cor(barcodes$GC_content2,barcodes$barcode_count,weight=barcodes$barcode_count)

a = paste0("Barcodes GC content:")
b = paste0("Min: ", min(barcodes_GC$X4))
c = paste0("Mean: ", mean(barcodes_GC$X4))
d = paste0("Median: ", median(barcodes_GC$X4))
e = paste0("Max: ", max(barcodes_GC$X4))
f = paste0("% of barcodes between 30-70% GC: ", nrow(new_frame)/96*100)
g = paste0("rho: ", barcodes_GC_pearson[1])
h = paste0("p-value: ", barcodes_GC_pearson[4])

pdf('barcodes_GC_summary.pdf',width=10,height=10)
plot(NA, xlim=c(0,10), ylim=c(0,10), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,8,a, pos=4)
text(1,7,b, pos=4)
text(1,6,c, pos=4)
text(1,5,d, pos=4)
text(1,4,e, pos=4)
text(1,3,f, pos=4)
text(1,2,g, pos=4)
text(1,1,h, pos=4)
dev.off()



################################################################################
################################################################################
################################################################################


