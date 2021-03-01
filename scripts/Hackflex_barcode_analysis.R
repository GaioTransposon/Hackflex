

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


# Coverage: 


cov_files <- grep(list.files(mydir,pattern="^reduced"), 
                  pattern='.tsv', value=TRUE)

# Coverage:

# construct an empty dataframe to build on 
df_to_fill_coverage <- data.frame(
  position = numeric(),
  coverage = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (cov_file in cov_files) {
  
  print(cov_file)
  
  coverage <-read.table(file.path(mydir,cov_file))
  
  colnames(coverage) <- c("position","coverage")
  
  id <- sub(".*interleaved2_", "", cov_file)
  id <- sub(".dedup.tsv", "", id)
  id <- gsub("_R1","",id)
  coverage$library=paste0(as.character(id))
  
  df_to_fill_coverage <- rbind(df_to_fill_coverage,coverage)
  
}


# function to get numbers per library: 
get_me_stats <- function(DF) {
  out <- DF %>% 
    dplyr::summarise(min=min(coverage),
                     median=median(coverage),
                     mean=mean(coverage),
                     sd=sd(coverage),
                     max=max(coverage),
                     n_sites=n())
  return(out)
}



sink(file = paste0(mydir,"library_analysis.txt"), 
     append = FALSE, type = c("output"))
paste0("########### coverage ###########")
paste0("##################################")
df_to_fill_coverage %>% 
  group_by(library) %>% 
  get_me_stats()
sink()


p1 <- ggplot(df_to_fill_coverage, aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=0.5)+
  facet_grid(cols = vars(library)) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(0,100) +
  labs(x="coverage per site",
       y="Frequency")

p2 <- df_to_fill_coverage %>% 
  dplyr::filter(coverage <= 3) %>%
  ggplot(., aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=1)+
  facet_grid(cols = vars(library)) +
  theme(legend.position="none")+
  xlim(0,4) +
  labs(x="coverage per site",
       y="Frequency")

z <- df_to_fill_coverage %>%
  group_by(library) %>% 
  dplyr::mutate(mean=mean(coverage),
                sd=sd(coverage))

dat_text <- data.frame(
  label = paste0("mean=", #\t doesn't work
                 as.character(round(unique(z$mean),2)),
                 " sd=",
                 as.character(round(unique(z$sd),2))),
  library   = unique(z$library)
)

p <- ggplot(z, aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=0.5)+
  facet_grid(cols = vars(library)) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(0,100) +
  labs(x="coverage per site",
       y="Frequency")

p1 <- p + geom_text(
  data    = dat_text,
  mapping = aes(x = Inf, y = Inf, label = label, hjust=1.0, vjust=1), #vjust=1 hjust=1.0
  size=2
)


########################################

# Correlation between libraries based on their coverage :
df_to_fill_coverage_wide <- df_to_fill_coverage %>%
  dplyr::distinct(position, library, .keep_all = TRUE) %>% 
  pivot_wider(names_from=library, values_from=coverage) # values_fn = {summary_fun}

df_to_fill_coverage_wide$position <- NULL

########################################
########################################

# Fragment size: 

frag_files = list.files(mydir,pattern="frag_sizes")

# construct an empty dataframe to build on 
df_to_fill_fragm_size <- data.frame(
  V1 = numeric(),
  order = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (frag_file in frag_files) {
  
  # read in file 
  frag_df <- read.table(file.path(mydir,frag_file), quote="\"", comment.char="")
  
  # transform fragment size values to absolute values (because pos and neg stand for forward and reverse reads)
  frag_df$V1 <- abs(frag_df$V1)
  
  # add column with sequential numbers to each dataset
  new_df <- cbind( frag_df, order=seq(nrow(frag_df)) ) 
  
  id <- sub(".*interleaved2_", "", frag_file)
  id <- sub(".dedup.txt", "", id)
  id <- gsub("_R1","",id)
  new_df$library=paste0(as.character(id))
  
  df_to_fill_fragm_size <- rbind(
    df_to_fill_fragm_size, 
    new_df
  )
  
}

########################################
########################################


# GC content: 

gc_files = list.files(mydir,pattern="GC_qc")
RL_files = list.files(mydir,pattern="^RL")

myDF <- data.frame(gc_files,RL_files, stringsAsFactors = F)

datalist = list()
for (row in 1:nrow(myDF)) {
  
  gc <- myDF[row,1]
  rl <- myDF[row,2]
  
  id <- gsub("GC_qc_reduced_trimmed2_trimmed_interleaved2_", "",gc)
  id <- gsub("_R1.tsv","",id)
  
  # read in files
  gc_df <- read.table(file.path(mydir,gc), quote="\"", comment.char="", header = TRUE)
  rl_df <- read.delim(file.path(mydir,rl), header=TRUE)
  
  mapped <- as.numeric(sum(rl_df$Count))
  
  gc_df$fractionOfReads <- (gc_df$fractionOfReads + 0.000000001)
  gc_df1 <- spread(gc_df, GCcontent, fractionOfReads)
  
  rownames(gc_df1) <- gc_df1$Sample
  gc_df1$Library <- NULL
  gc_df2 <- as.data.frame(t(gc_df1[,-1]))
  colnames(gc_df2) <- gc_df1$Sample
  
  gc_df2$diff <- gc_df2[,1]/gc_df2[,2]
  
  # mapped 
  gc_df2$reads <- gc_df2[,1]* mapped
  gc_df3 <- cbind(gc_df2, gc_df[,3, drop=FALSE])
  head(gc_df2)
  rho <- wtd.cor(gc_df3$GCcontent,gc_df3$diff,weight=gc_df3$reads)
  rho <- rho[1,1]
  
  gc_df3$library <- id 
  gc_df3$rho <- rho
  
  colnames(gc_df3) <- c("deduped_bam","Reference","diff","reads","GCcontent","library","rho")
  
  datalist[[row]] <- gc_df3 # add it to your list
  
}

big_data = do.call(rbind, datalist)

smooth_GC <- big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=library))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  geom_smooth() +
  ylab("ratio observed/expected reads") +
  theme(legend.position="none") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)

straight_GC <- big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=library))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  stat_smooth(method="lm", se=FALSE) +
  ylab("ratio observed/expected reads") +
  theme(legend.title = element_blank(),
        legend.position="top") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)


# ########################################
# ########################################
# 
# # PHRED scores plotting
# 
# # open file contaning all the PHRED scores:
# PHRED_df <- read_csv(file.path(phred_dir,"PHRED_scores.csv"))
# 
# PHRED_df$lib <- as.character(paste0(PHRED_df$library,"_",PHRED_df$read))
# 
# 
# ########################################
# ########################################


# RL content: 

rl_files = list.files(mydir,pattern="^RL")

# construct an empty dataframe to build on 
rl_to_fill <- data.frame(
  library = character(),
  Readlength = numeric(),
  Count = numeric(),
  Fraction = numeric(),
  stringsAsFactors = FALSE
)

for (rl_file in rl_files) {
  
  # read in file
  rl_df <- read.table(file.path(mydir,rl_file), quote="\"", comment.char="", header = TRUE)
  
  rl_df <- rl_df %>% 
    dplyr::select(Sample,Readlength,Count, Fraction)
  
  rl_df$library <- str_replace_all(rl_df$Sample, "reduced_trimmed2_trimmed_interleaved2_", "")
  rl_df$library <- gsub("_R1.dedup","",rl_df$library)
  rl_df$Sample <- NULL
  
  rl_to_fill <- rbind(rl_to_fill, rl_df)
}

# expand rows based on count (to create a density plot, like the fragm size one)
rl_to_fill_expand <- setDT(expandRows(rl_to_fill, "Count"))


########################################
########################################


# ME data (Alfred output): 

me_files = list.files(mydir,pattern="^ME")

myDF <- data.frame(me_files, stringsAsFactors = F)


MElist = list()
for (row in 1:nrow(myDF)) {
  
  me_file <- myDF[row,1]
  
  # read in file
  me_df <- read.table(file.path(mydir,me_file), quote="\"", comment.char="", header = TRUE)
  
  me_df$Sample <- gsub("reduced_trimmed2_trimmed_interleaved2_", "",me_df$Sample)
  me_df$Sample <- gsub("_R1.dedup", "",me_df$Sample)
  
  me_df <- me_df %>%
    dplyr::mutate(library=Sample) %>%
    select(library, everything()) # bring at first position
  
  MElist[[row]] <- me_df # add it to your list
  
}


ME_data = do.call(rbind, MElist)

# store names
n <- ME_data$Sample
# transpose all but the first column (name)
ME_data <- as.data.frame(t(ME_data[,-2]))
colnames(ME_data) <- n
ME_data$Sample <- factor(row.names(ME_data))
rownames(ME_data) <- NULL

fwrite(x=ME_data, file=paste0(mydir,"Alfred_ME_data.csv"))


########################################
########################################


# barcode distribution: 

# number of mapped reads 

barcode_count <- read.table(file.path(mydir,"reads_after_cleaning.tsv"), quote="\"", comment.char="", header = FALSE)
head(barcode_count)

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

#barcode distribution (with boxplot)
pdf(paste0(mydir,"barcode_distribution.pdf"))
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(barcode_count$V2,
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
hist(barcode_count$V2, 
     col=rgb(0.2,0.8,0.5,0.5), 
     breaks=20, 
     main = NULL, xlab = "barcode count", ylab = "Frequency")
opar
dev.off()

##############

# Plot number of reads assigned to barcode vs 96 barcode (rank sorted)

#sort number of reads (ascending)
barcode_count_2 <- barcode_count$V2
barcode_count_3 <- data.frame(sort(barcode_count_2))

#add column for wells
barcode_count_3$X2 <- 1:nrow(barcode_count_3) 

pdf(paste0(mydir,"read_count_per_barcode.pdf"))
barplot(barcode_count_3$sort.barcode_count_2~barcode_count_3$X2, xlab="", ylab="", 
        xaxt='n', main="96-plex barcode counts") 
title(xlab="96 barcodes (rank-sorted)", line=1.0, cex.lab=1.5)
title(ylab="number of reads assigned to barcode", line=2.5, cex.lab=1.5)
grid(nx = NA, ny = 10, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
dev.off()


########################################

# Barcode GC content:

# Barcode GC content:

all_fastq_headers <- read_table2(paste0(mydir,"all_fastq_headers.tsv"), col_names = FALSE)

NROW(all_fastq_headers)
head(all_fastq_headers)

barcodes_seq <- all_fastq_headers %>%
  dplyr::select(X1,X3) 

NROW(barcodes_seq)
head(barcodes_seq)

# cleaning up (to reduce file size)
barcodes_seq$X1 <- gsub("fastq_files_headers_HF-barcode-","HF",barcodes_seq$X1)

barcodes <- as.data.frame(barcodes_seq)

head(barcodes)
NROW(barcodes)

# # take the top 5 unique barcodes to each lib, get the sd among those. larger the sd the better
# barcodes_sd <- barcodes %>%
#   drop.levels() %>%
#   dplyr::mutate(n=n()) %>%
#   group_by(X1,X3) %>%
#   dplyr::summarise(freq=sum(n)) %>%
#   group_by(X1) %>%
#   slice_max(order_by = freq,n=5) %>%
#   group_by(X1) %>%
#   dplyr::summarise(sd=sd(freq))
# hist(barcodes_sd$sd)

# take the original barcode
barcodes <- barcodes %>%
  drop.levels() %>%
  dplyr::mutate(n=n()) %>%
  group_by(X1,X3) %>%
  dplyr::summarise(freq=sum(n)) %>%
  group_by(X1) %>%
  slice_max(order_by = freq,n=1)

head(barcodes)

barcodes <- cSplit(barcodes, "X3", sep = ":")
barcodes <- cSplit(barcodes, "X3_4", sep = "+")

barcodes <- barcodes %>%
  dplyr::select(X1, freq, X3_1, X3_4_1, X3_4_2) %>%
  dplyr::rename(library = X1,
                barcode_count = freq,
                read_direction = X3_1,
                barcode_seq1 = X3_4_1,
                barcode_seq2 = X3_4_2)


barcodes <- barcodes %>%
  dplyr::mutate(G1 = str_count(barcode_seq1, "G"),
                C1 = str_count(barcode_seq1, "C"),
                G2 = str_count(barcode_seq2, "G"),
                C2 = str_count(barcode_seq2, "C"),
                GC_content1 = ((G1 + C1) / str_length(barcode_seq1) * 100),
                GC_content2 = ((G2 + C2) / str_length(barcode_seq2) * 100))

head(barcodes)

#plot
pdf("BarcodeGC_vs_barcode.pdf")
opar=par(ps=14)
plot(barcodes$GC_content1, barcodes$barcode_count, xlab="GC content", ylab="barcode count")
abline(lm(barcodes$barcode_count ~ barcodes$GC_content1, weights=barcodes$barcode_count))
opar 
plot(barcodes$GC_content2, barcodes$barcode_count, xlab="GC content", ylab="barcode count")
abline(lm(barcodes$barcode_count ~ barcodes$GC_content2, weights=barcodes$barcode_count))
opar
dev.off()

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


# plot
pdf(paste0(mydir,'out.pdf'))
# plot coverage 
ggarrange(p1,p2,                                    
          labels = c("A","B"), 
          nrow = 2             
)
# plot correlation between libs based on their coverage
df_to_fill_coverage_wide %>%
  corrr::correlate(method = "pearson") %>%
  # Re-arrange a correlation data frame 
  # to group highly correlated variables closer together.
  corrr::rearrange(method = "MDS", absolute = FALSE) %>%
  corrr::shave() %>% 
  corrr::rplot(shape = 19, colors = c("red", "green"), print_cor = TRUE) # dot color = corr; dot size = absolute value of corr
# plot fragment size 
ggplot(df_to_fill_fragm_size, aes(V1, colour=library, fill=library)) + 
  scale_x_continuous(limits=c(0,1000)) + 
  xlab("fragment size (bp)")+
  geom_density(alpha=0.01) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
# plot GC content bias
ggarrange(smooth_GC, straight_GC, 
          nrow = 2)
# # plot PHRED scores
# PHRED_df %>% 
#   ggplot(.,
#          aes(x=read_position, y=PHRED_means, colour=library, shape = read,
#              group=interaction(library, read))) + 
#   geom_point(alpha=0.8) + 
#   geom_line()+
#   xlab("read position (bp)") +
#   ylab("average PHRED score") +
#   theme_bw() +
#   scale_x_continuous(breaks = c(0,50,100,150,200,250,300), lim = c(0, 300)) +
#   scale_y_continuous(breaks = c(30,32,34,36,38), lim = c(29, 38))+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"),
#         legend.position="top",
#         legend.title = element_blank(),
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=18)) 
# plot read lengths
ggplot(rl_to_fill_expand, aes(Readlength, colour=library, fill=library)) + 
  scale_x_continuous(limits=c(100,400)) + 
  xlab("fragment size (bp)")+
  geom_density(alpha=0.01) +
  facet_grid(cols = vars(library)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=8),
        axis.title=element_text(size=14))
dev.off()


################################################################################
################################################################################
################################################################################


