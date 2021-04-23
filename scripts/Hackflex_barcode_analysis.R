

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
barcode_source_data <- "~/Desktop/MG1655/goal_barcode/source_data/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"


########################################
########################################


demux_clean <- read_csv(paste0(barcode_source_data, "demux_clean.tsv"), col_names = FALSE)
demux_clean <- demux_clean[6:677,]
demux_clean$seq <- rep(seq(1,7), 96)

demux_clean <- as.data.frame(demux_clean)

counts <- demux_clean %>%
  dplyr::filter(seq=="5") # filter barcode counts (all)

counts$X1 <- gsub("</BarcodeCount>","",counts$X1)
counts$X1 <- gsub("<BarcodeCount>","",counts$X1)

counts <- as.data.frame(as.numeric(counts$X1))
colnames(counts) <- "counts"

TOT_BARCODE_COUNTS <- sum(counts$counts)
sort(counts$counts)
head(counts)
counts$counts<45000

a <- counts %>%
  dplyr::summarise(barcodes_combos=n(),
                   Min=min(counts),
                   Mean=mean(counts),
                   Median=median(counts),
                   Sd=sd(counts),
                   Max=max(counts),
                   Sum=sum(counts),
                   Coeff.variation=sd(counts)/mean(counts))
a_wo_outliers <- counts %>%
  dplyr::filter(counts>45000) %>%
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

counts_stats_transposed <- counts_transposed %>%
  dplyr::select(summary, everything())

fwrite(x=counts_stats_transposed, file=paste0(barcode_libs,"barcode_summary.csv"))


pdf(paste0(barcode_libs,"barcode_distribution.pdf"))
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(counts$counts, ylim=c(0,200000), 
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
     breaks = seq(from=1, to=200000, by=10000),
     main = NULL, xlab = "barcode count", ylab = "Frequency")
dev.off()

# there's none with more than 200,000 barcode counts 
which(counts$counts>200000)

########################################
########################################


# barcodes GC content coverage bias : 
fastq <- read_delim(file.path(barcode_source_data,"all_fastq_headers_clean.tsv"),":", 
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
  slice_max(order_by = barcode_count,n=1) %>% # change to n>1 if you want the top n most frequent barcode sequences
  group_by(X1) %>%
  dplyr::mutate(barcode_count=sum(barcode_count))

head(z1)

z2 <- z1[,2:3] %>%
  distinct() %>%
  dplyr::mutate(i5 = sub("\\+.*", "", X8)) %>%
  dplyr::mutate(i7 = str_extract(X8, '\\b[^+]+$')) %>%
  dplyr::select(X8,i5,i7,barcode_count) %>%
  dplyr::mutate(Gs = str_count(X8, "G"),
                Cs = str_count(X8, "C"),
                GC_content = ((Gs+Cs) / (str_length(X8)-1))*100) %>%
  dplyr::select(i5,i7,barcode_count,GC_content)
z2 <- as.data.frame(z2)
head(z2)


#contamination_rate <- sum(u_summ)*2/sum(z1$barcode_count)*100

true_i5 <- as.list(z2$i5)
true_i7 <- as.list(z2$i7)


###
# No barcode contained 3 or more identical bases in a row: 
stringr::str_count(z2$i7, "AAA+")
stringr::str_count(z2$i7, "TTT+")
stringr::str_count(z2$i7, "CCC+")
stringr::str_count(z2$i7, "GGG+")

stringr::str_count(z2$i5, "AAA+")
stringr::str_count(z2$i5, "TTT+")
stringr::str_count(z2$i5, "CCC+")
stringr::str_count(z2$i5, "GGG+")
###

distribution_text <- z2 %>%
  dplyr::summarise(mean=mean(barcode_count),
                   median=median(barcode_count),
                   sd=sd(barcode_count)) %>%
  dplyr::mutate(label=paste0("mean=",round(mean,2),
                             "\n median=",round(median,2),
                             "\n sd=",round(sd,2)))



extreme_GC <- z2 %>% dplyr::filter(GC_content<20|GC_content>80)

GC_text <- z2 %>%
  dplyr::summarise(mean=mean(GC_content),
                   median=median(GC_content),
                sd=sd(GC_content),
                `<20%`=NROW(which((GC_content<20)==TRUE))/NROW(GC_content)*100,
                `>80%`=NROW(which((GC_content>80)==TRUE))/NROW(GC_content)*100) %>%
  dplyr::mutate(label=paste0("mean=",round(mean,2),
                             "\n median=",round(median,2),
                             "\n sd=",round(sd,2),
                             "\n <20%=",paste0(round(`<20%`,2),"%"),
                             "\n >80%=",paste0(round(`>80%`,2)),"%"))


# Pearson coefficient of correlation
barcodes_GC_rho <- wtd.cor(z2$GC_content,z2$barcode_count,weight=z2$barcode_count)
pval <- barcodes_GC_rho[1,4]
rho <- barcodes_GC_rho[1,1]
#extreme_GC_pearson <- wtd.cor(extreme_GC$GC_content,extreme_GC$barcode_count,weight=extreme_GC$barcode_count)
barcodes_GC_rho <- as.data.frame(barcodes_GC_rho) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))

########################################
########################################


# PLOT: 

pdf(paste0(barcode_libs,"barcodes_distribution&GC_content.pdf"))
# barcode distribution: 
z2 %>%
  dplyr::select(barcode_count) %>%
  distinct() %>%
  ggplot(., aes(x=barcode_count)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  ggtitle("barcode distribution (only perfect barcode match)") +
  geom_text(
    data    = distribution_text,
    mapping = aes(x = Inf, y = Inf, label = label, hjust=1.0, vjust=1), #vjust=1 hjust=1.0
    size=3
  )
# barcode counts vs barcode GC content
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  ggplot(., aes(x=GC_content,y=barcode_count))+
  geom_point()+
  #geom_smooth(method="lm", se=FALSE) +
  #stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = "R", method = "pearson") +
  labs(x="barcode GC content (%)",
       y="barcode count") +
  ggtitle("GC content coverage bias (only perfect barcode match) \n with barcode counts as weights") +
  geom_text(
    data    = GC_text,
    mapping = aes(x = 70, y = Inf, label = label, hjust=1.0, vjust=1), #vjust=1 hjust=1.0
    size=3
  ) +
  geom_text(
    data    = barcodes_GC_rho,
    mapping = aes(x = 20, y = Inf, label = label, hjust=1.0, vjust=1), #vjust=1 hjust=1.0
    size=3
  )
# barcodes of each library, their count, and their GC content
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  dplyr::mutate(lib=as.character(seq(1:NROW(.)))) %>%
  ggplot(., aes(x=fct_reorder(lib, GC_content),y=barcode_count, color=GC_content))+ 
  geom_point(alpha=0.8)+
  theme(axis.text.x=element_blank()) +
  scale_color_gradientn(colours = rainbow(5)) +
  ggtitle("Barcode count vs barcode GC content (perfect barcode match only)")
# barcode counts - rank sorted 
z2 %>%
  dplyr::select(barcode_count,GC_content) %>%
  distinct() %>%
  dplyr::mutate(lib=as.character(seq(1:NROW(.)))) %>%
  ggplot(., aes(x=fct_reorder(lib, barcode_count), y=barcode_count))+
  geom_boxplot()+
  labs(x="Barcodes (rank-sorted)",
       y="number of reads assigned to barcode") +
  ggtitle("Rank-sorted barcodes and their count (perfect barcode match only)")
dev.off()



########################################
########################################



# barcode cross-contamination


# list of all sequencing libs from batch 3 
sequencing_run_barcodes <- read_xlsx(paste0(barcode_source_data,"sequencing_run_barcodes_20210212.xlsx"), 
                                skip = 19)

x <- sequencing_run_barcodes %>%
  dplyr::select(Sample_ID,index,index2)

# exclude all that contain the word "barcode" 
x <- dplyr::filter(x, !grepl("barcode",Sample_ID))
x$X10 <- paste0(x$index,"+",x$index2)


R1_undet_headers <- read_delim(paste0(barcode_source_data,"R1_undet_headers"),":", 
                               escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE)
# R2_undet_headers <- read_delim(paste0(barcode_source_data,"R2_undet_headers"),":", 
#                                escape_double = FALSE, col_names = FALSE,
#                                trim_ws = TRUE)
# 
# NROW(R1_undet_headers)==NROW(R2_undet_headers)

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
100-sum(u_summ)/NROW(fastq)*100


pdf(paste0(barcode_libs,"barcode_cross_contamination.pdf"))
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
  ggplot(aes(i5, i7, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "gold1", high = "red4",
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



################################################################################
################################################################################


# Analysis of Hackflex barcode libraries based on insert size, GC bias, and coverage 


################################################################################

# Read stats: 

stats <-read.table(file.path(barcode_libs,"reads_stats.tsv"))
head(stats)
stats$V5 <- NULL
stats$V9 <- NULL

stats$V1 <- gsub("interleaved_","",stats$V1)
stats$V1 <- gsub(".fastq","",stats$V1)

stats <- stats %>% 
  dplyr::mutate(library=as.factor(recode_fun(V1))) %>%
  drop.levels()

read_lengths <- stats %>%
  dplyr::select(library,V2,V6,V10)
colnames(read_lengths) <- c("library","before_cleaning","after_cleaning","after_resizing")
read_lengths$var <- "read_length"

n_reads <- stats %>%
  dplyr::select(library,V3,V7,V11)
colnames(n_reads) <- c("library","before_cleaning","after_cleaning","after_resizing")
n_reads$var <- "read_count"

n_bp <- stats %>%
  dplyr::select(library,V4,V8,V12)
colnames(n_bp) <- c("library","before_cleaning","after_cleaning","after_resizing")
n_bp$var <- "bp_count"


x <- cbind(read_lengths[,1:4],n_reads[,2:4], n_bp[,2:4])

to_paste <- c("", rep("avg read length",3), rep("# reads",3), rep("# bp",3))
colnames(x) <- paste(to_paste, colnames(x), sep = " ")

fwrite(x=x, file=paste0(barcode_libs,paste0(species,"_",goal,"_lib_processing_stats.csv")))

# coverage 
mean(x$`# bp after_cleaning`)/4600000
mean(x$`# bp after_resizing`)/4600000

################################################################################


# Insert size: 


IS_files <- grep(list.files(barcode_libs,pattern="^picardIS"), 
                 pattern='.txt', value=TRUE)


# construct an empty dataframe to build on 
df_to_fill_insert_size <- data.frame(
  insert_size = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (IS_file in IS_files) {
  
  # read in file 
  IS_df <- read_delim(file.path(barcode_libs,IS_file), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 10)
  
  IS_df <- IS_df %>%
    uncount(All_Reads.fr_count) 
  
  id <- sub(".*interleaved_", "", IS_file)
  id <- gsub(".dedup.txt","",id)
  IS_df$library=paste0(as.character(id))
  
  df_to_fill_insert_size <- rbind(
    df_to_fill_insert_size, 
    IS_df
  )
  
}


# expand rows based on Count, compute median insert sizes
med_IS <- df_to_fill_insert_size %>%
  group_by(library) %>%
  dplyr::summarize(median=round(median(insert_size),2),
                   mean=round(mean(insert_size),2),
                   sd=round(sd(insert_size),2),
                   var=round(var(insert_size),2))


IS_stats <- med_IS
IS_stats$mean_bins <- cut(IS_stats$mean, breaks = 3)
IS_stats1 <- IS_stats %>%
  group_by(mean_bins) %>%
  tally() %>%
  dplyr::mutate(perc_mean=n/sum(n)*100)
IS_stats1$IS_mean <- paste0(IS_stats1$mean_bins,": n=",IS_stats1$n," ",round(IS_stats1$perc_mean,2),"%")

fwrite(IS_stats1,file=paste0(barcode_libs,"HF-barcodes_IS_size_mean.csv"))


IS_stats$median_bins <- cut(IS_stats$median, breaks = 3)
IS_stats2 <- IS_stats %>%
  group_by(median_bins) %>%
  tally() %>%
  dplyr::mutate(perc_median=n/sum(n)*100)
IS_stats2$IS_median <- paste0(IS_stats2$median_bins,": n=",IS_stats2$n," ",round(IS_stats2$perc_median,2),"%")


fwrite(IS_stats2,file=paste0(barcode_libs,"HF-barcodes_IS_size_median.csv"))


########################################
########################################


# GC content - Picard


gc_files <- grep(list.files(barcode_libs,pattern="^picardGC_red"), 
                 pattern='.txt', value=TRUE)


GC_DF <- data.frame(
  GC = numeric(),
  WINDOWS = numeric(),
  MEAN_BASE_QUALITY = numeric(),
  NORMALIZED_COVERAGE = numeric(),
  library = character(),
  rho = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)

for (gc_file in gc_files) {
  
  gc <- gc_file
  
  id <- sub(".*interleaved_", "",gc)
  id <- gsub(".dedup.txt","",id)
  
  # read in files
  gc_df <- read_delim(file.path(barcode_libs,gc_file),
                      "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)
  
  # rename the library and select cols
  gc_df <- gc_df %>%
    dplyr::mutate(library=paste0(id)) %>%
    dplyr::select(GC,WINDOWS,MEAN_BASE_QUALITY,NORMALIZED_COVERAGE,library) %>%
    dplyr::filter(NORMALIZED_COVERAGE>0)
  
  head(gc_df)
  
  # get correlation between GC content and ratio, using fraction of reads as weight
  rho <- wtd.cor(gc_df$GC,
                 gc_df$NORMALIZED_COVERAGE,
                 weight=gc_df$WINDOWS)
  pval <- rho[1,4]
  rho <- rho[1,1]
  
  myDF <- gc_df %>%
    dplyr::mutate(rho=rho,
                  pval=pval) 
  
  head(myDF)
  
  GC_DF <- rbind(GC_DF, myDF)
  
}



GC_DF_text <- GC_DF %>% dplyr::select(library,rho,pval) %>%
  distinct() %>%
  dplyr::arrange(library) %>% # to get the factors order
  dplyr::mutate(label=paste0(library),
                label_rho=paste0("R=",
                                 round(rho,3)),
                label_pval=paste0("p=",
                                  round(pval,3)))


GC_stats <- GC_DF_text
GC_stats$rho_bins <- cut(GC_stats$rho, c(-1, -0.8, -0.6, -0.4, -0.2, 0, 
                                               0.2, 0.4, 0.6, 0.8, 1))
GC_stats <- GC_stats %>%
  group_by(rho_bins) %>%
  tally() %>%
  dplyr::mutate(perc=n/sum(n)*100)

GC_stats$R_GC <- paste0(GC_stats$rho_bins,": n=",GC_stats$n," ",round(GC_stats$perc,2),"%")

fwrite(GC_stats,file=paste0(barcode_libs,"HF-barcodes_GC_bias.csv"))


GC_DF_text$pval <- cut(GC_DF_text$pval, c(0, 0.01, 0.05, 1))

pdf(paste0(barcode_libs,"GC_plot.pdf"), width=3,height=4)
ggplot(GC_DF_text, aes(x=reorder(library,rho),y=rho)) +
  geom_col(width = 0.5)+
  coord_flip() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.title.y=element_text(size=10), 
        legend.title = element_text()) +
  #ggtitle("Weighted correlation of coverage with GC content for 96 Hackflex libraries") +
  labs(fill="p-value",
       x="HF-barcode library",
       y="Pearson R") #+
  #scale_fill_manual(values = c("red", 'orange1', 'darkmagenta', 'navyblue'))  
dev.off()





########################################
########################################
