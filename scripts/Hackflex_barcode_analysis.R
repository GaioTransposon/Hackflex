

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
library(readxl)
library(reshape2)
library(textshape)


# set directories and select samples: 

barcode_libs <- "~/Desktop/MG1655/goal_barcode/"
barcode_dir <- "~/Desktop/MG1655/HF-barcode/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"


########################################
########################################


demux_clean <- read_csv(paste0(barcode_dir, "demux_clean.tsv"), col_names = FALSE)
demux_clean <- demux_clean[6:677,]
demux_clean$seq <- rep(seq(1,7), 96)

demux_clean <- as.data.frame(demux_clean)

counts <- demux_clean %>%
  dplyr::filter(seq=="5") # barcode count (all)

counts$X1 <- gsub("</BarcodeCount>","",counts$X1)
counts$X1 <- gsub("<BarcodeCount>","",counts$X1)

counts <- as.data.frame(as.numeric(counts$X1))
colnames(counts) <- "counts"

TOT_BARCODE_COUNTS <- sum(counts$counts)


NROW(R2_unpaired_headers)/TOT_BARCODE_COUNTS



a <- counts %>%
  dplyr::summarise(barcodes_combos=n(),
                   Min=min(counts),
                   Mean=mean(counts),
                   Median=median(counts),
                   Sd=sd(counts),
                   Max=max(counts),
                   Sum=sum(counts),
                   Coeff.variation=sd(counts)/mean(counts))
counts_transposed <- as.data.frame(t(a), make.names = TRUE)
counts_transposed$summary <- rownames(counts_transposed)
colnames(counts_transposed) <- c("value","summary")

counts_transposed <- counts_transposed %>%
  dplyr::select(summary, everything())



fwrite(x=counts_transposed, file=paste0(barcode_dir,"barcode_stats.csv"))


pdf(paste0(barcode_dir,"_barcode_distribution.pdf"))
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(counts$counts,
        main = NULL,
        xlab = NULL,
        ylab = NULL,
        axes = FALSE,
        col = "grey",
        border = "black",
        horizontal = TRUE,
        notch = TRUE
)
#par(mar=c(4, 3.1, 1.1, 2.1))
#bootm left top right
par(mar=c(5,6,4,2)+0.1)
#x and y labels font size with
opar=par(ps=14)
hist(counts$counts, xlim=c(0,200000),
     main = NULL, xlab = "barcode count", ylab = "Frequency")
opar
dev.off()




########################################
########################################


# barcodes GC content coverage bias : 
fastq <- read_delim(file.path(barcode_libs,"all_fastq_headers_clean.tsv"),":", 
           escape_double = FALSE, col_names = FALSE,
           trim_ws = TRUE)
barcodes <- fastq

barcodes <- barcodes %>% 
  dplyr::select(X1,X8) 

head(barcodes)
NROW(barcodes)

barcodes$count <- 1
z <- barcodes %>% 
  group_by(X1,X8) %>%
  dplyr::summarise(barcode_count=sum(count))

head(z)

# take the original barcode (this is the most frequent one)
z1 <- z %>%
  group_by(X1) %>%
  slice_max(order_by = barcode_count,n=1) # change to n>1 if you want the top n most frequent barcode sequences

head(z1)

z2 <- z1[,1:3] %>%
  distinct() %>%
  dplyr::mutate(i5 = sub("\\+.*", "", X8)) %>%
  dplyr::mutate(i7 = str_extract(X8, '\\b[^+]+$')) %>%
  dplyr::select(X1,X8,i5,i7,barcode_count) %>%
  dplyr::mutate(Gs = str_count(X8, "G"),
                Cs = str_count(X8, "C"),
                GC_content = ((Gs+Cs) / (str_length(X8)-1))*100) %>%
  dplyr::select(X1,i5,i7,barcode_count,GC_content)
z2 <- as.data.frame(z2)
head(z2)

true_i5 <- as.list(z2$i5)
true_i7 <- as.list(z2$i7)

extreme_GC <- z2 %>% dplyr::filter(GC_content<30|GC_content>70)
NROW(extreme_GC)/NROW(z2)*100
#it means that 13.54% of all the barcodes have more than 70% or less than 30% GC content

# Pearson coefficient of correlation
barcodes_GC_pearson <- wtd.cor(z2$GC_content,z2$barcode_count,weight=z2$barcode_count)
extreme_GC_pearson <- wtd.cor(extreme_GC$GC_content,extreme_GC$barcode_count,weight=extreme_GC$barcode_count)


########################################
########################################


# PLOT: 

pdf(paste0(barcode_dir,"barcodes_distribution&GC_content.pdf"))
# barcode distribution: 
z2 %>%
  dplyr::select(barcode_count) %>%
  distinct() %>%
  ggplot(., aes(x=barcode_count)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")
# barcode counts vs barcode GC content
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  ggplot(., aes(x=GC_content,y=barcode_count))+
  geom_point()+
  stat_smooth(method="lm", se=FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(x="barcode GC content (%)",
       y="barcode count")
# barcodes of each library, their count, and their GC content
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  dplyr::mutate(lib=as.character(seq(1:NROW(.)))) %>%
  ggplot(., aes(x=fct_reorder(lib, GC_content),y=barcode_count, color=GC_content))+ 
  geom_point(alpha=0.8)+
  theme(axis.text.x=element_blank()) +
  scale_color_gradientn(colours = rainbow(5))
# barcode counts - rank sorted 
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  dplyr::mutate(lib=as.character(seq(1:NROW(.)))) %>%
  ggplot(., aes(x=fct_reorder(lib, barcode_count), y=barcode_count))+
  geom_boxplot()+
  labs(x="Barcodes (rank-sorted)",
       y="number of reads assigned to barcode")
dev.off()



########################################
########################################



# barcode cross-contamination


# list of all sequencing libs from batch 3 
X20210212_HF_16S_V3 <- read_csv("Desktop/MG1655/HF-barcode/20210212_HF_16S_V3.csv", 
                                skip = 19)

x <- X20210212_HF_16S_V3 %>%
  dplyr::select(Sample_ID,index,index2)

# exclude all that contain the word "barcode" 
x <- dplyr::filter(x, !grepl("barcode",Sample_ID))
x$X10 <- paste0(x$index,"+",x$index2)


R1_undet_headers <- read_delim(paste0(barcode_dir,"R1_undet_headers"),":", 
                               escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE)
R2_undet_headers <- read_delim(paste0(barcode_dir,"R2_undet_headers"),":", 
                               escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE)

NROW(R1_undet_headers)==NROW(R2_undet_headers)

R1_undet0 <- R1_undet_headers %>%
  dplyr::select(X10) %>%
  dplyr::mutate(i5 = sub("\\+.*", "", X10)) %>%
  dplyr::mutate(i7 = str_extract(X10, '\\b[^+]+$')) %>%
  dplyr::select(X10,i5,i7)

# R2_unpaired <- R2_unpaired_headers %>%
#   dplyr::select(X10) %>%
#   dplyr::mutate(i5 = sub("\\+.*", "", X10)) %>%
#   dplyr::mutate(i7 = str_extract(X10, '\\b[^+]+$')) %>%
#   dplyr::select(i5,i7)


# keep the true i5 and i7 present in list: 
NROW(R1_undet0)
R1_undet1 <- R1_undet0 %>% 
  dplyr::filter(i5 %in% true_i5) %>%
  dplyr::filter(i7 %in% true_i7) 
NROW(R1_undet1)




# exclude all the barcode combos that ended up in undetermined, that actually derive from other libraries (not all libs were demuxed): 
R1_undet2 <- R1_undet1 %>%
  dplyr::filter(!X10 %in% x$X10)






u  <- as.data.frame(R1_undet2)
u
u$count <- 1

NROW(unique(u$i5))
NROW(unique(u$i7))

u_sum <- u %>%
  group_by(i5,i7) %>%
  dplyr::summarise(sum=sum(count))


NROW(unique(u_sum[,1:2]))

# it worked. 
u_sum <- as.data.frame(u_sum)
head(u_sum)
sum(u_sum$sum)


u_summ <- acast(u_sum, i5~i7, value.var="sum", fill = 0)


swap <- function(matrixRow,x,y){
  #x is diagonal index
  #y is max of the row
  indexY <- which(matrixRow == y)
  valX <- matrixRow[x]
  matrixRow[x] <- y
  matrixRow[indexY] <- valX
  return(matrixRow)
}


mat <- u_summ
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}


# total barcode cross-contamination
sum(u_summ)/NROW(fastq)*100



pdf(paste0(barcode_dir,"barcode_cross_contamination.pdf"))
mat <- u_summ
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}
mat %>%
  cluster_matrix() %>%
  tidy_matrix('i5', 'i7') %>%
  dplyr::mutate(
    i5 = factor(i5, levels = unique(i5)),
    i7 = factor(i7, levels = unique(i7))        
  ) %>%
  group_by(i5) %>%
  ggplot(aes(i5, i7, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red", 
                      na.value = "white") +
  geom_text(aes(label = round(value, 1)), size=0.8) +
  theme(axis.text.x=element_text(angle=90))+
  theme(
    axis.text.y = element_text(size = 5)   ,
    axis.text.x = element_text(
      size = 5, 
      hjust = 1, 
      vjust = 1, 
      angle = 45
    ),   
    legend.position = 'bottom',
    legend.key.height = grid::unit(.1, 'cm'),
    legend.key.width = grid::unit(.5, 'cm'),
    legend.text = element_text(angle = 90,size=5)
  ) 
dev.off()







