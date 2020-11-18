
# set working directory for all input and output files
setwd("~/Desktop/MG1655/MG1655_final")

# Price per sample comparison:

library(readr)
library(dplyr)
library(weights)

price_per_sample_comparison <- read_csv("price_per_sample_comparison.csv")

pdf("price_per_sample.pdf")
library(ggplot2)
a <- ggplot(price_per_sample_comparison, aes(x=samples, y=price, color=method)) + geom_step() 
b <- a + scale_x_continuous(breaks = seq(0,4800,300)) + scale_y_continuous(breaks = seq(0,75,5))
c <- b + xlab("number of samples") + ylab("price per sample (A$)") 
d <- c + scale_color_manual(values=c("#7CAE00", "#CC79A7", "#00BFC4"))
e <- d + geom_point(aes(shape=method, color=method, size=method))+
  scale_shape_manual(values=c(1, 2, 0))+
  scale_color_manual(values=c("#7CAE00", "#CC79A7", "#00BFC4"))+
  scale_size_manual(values=c(2,2,2))
f <- e + theme_minimal() + theme(legend.title=element_blank())
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 16),
  legend.text = element_text(size=13))
g <- f + My_Theme
h <- g + theme(legend.position="top")
h
dev.off()


# price fold difference : 

price_fold_difference <- read_csv("price_fold_difference.csv")
library(ggthemes)

# save manually as pdf
library(ggthemes)
pdf("price_fold_difference.pdf")
v <- ggplot(price_fold_difference, aes(x=samples, y=fold)) + geom_bar(stat = "identity") +
  xlab("number of samples") + ylab("Fold difference in price") 
w <- v + theme(axis.text=element_text(size=11),
               axis.title=element_text(size=16))
z <- w + theme_minimal()
z
dev.off()
########################################

# Barcode count distribution

barcode_count <- read_csv("barcode_count.csv", col_names = FALSE)

summary(barcode_count$X2)
sd(barcode_count$X2)
sum(barcode_count$X2)
#coeff of variation 
sd(barcode_count$X2)/mean(barcode_count$X2)

# coeff of variation without the two outliers (percentage)
# subset the dataframe to exclude outliers with read count lower than 50
barcode_count_no_outliers <- subset(barcode_count, X2>50)

a = paste0("Barcode count summary:")
b = paste0("Min: ", summary(barcode_count$X2)[1])
c = paste0("1st Qu.: ", summary(barcode_count$X2)[2])
d = paste0("Median: ", summary(barcode_count$X2)[3])
e = paste0("Mean: ", summary(barcode_count$X2)[4])
f = paste0("3rd Qu.: ", summary(barcode_count$X2)[5])
g = paste0("Max: ", summary(barcode_count$X2)[6])
h = paste0("Standard deviation: ", sd(barcode_count$X2))
i = paste0("Sum: ", sum(barcode_count$X2))
l = paste0("Coeff. of variation (incl. outliers): ", sd(barcode_count$X2)/mean(barcode_count$X2))
m = paste0("Coeff. of variation (excl. outliers): ", sd(barcode_count_no_outliers$X2)/mean(barcode_count_no_outliers$X2))

pdf('barcode_summary.pdf',width=11,height=11)
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
pdf("barcode_distribution.pdf")
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(barcode_count$X2,
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
hist(barcode_count$X2, col=rgb(0.2,0.8,0.5,0.5), breaks=seq(0, 10000, 500), main = NULL, xlab = "barcode count", ylab = "Frequency")
opar
dev.off()


##############

# Barcode GC content:

barcodes_GC <- read_delim("barcodes_GC.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

#plot
pdf("BarcodeGC_vs_barcode.pdf")
opar=par(ps=14)
plot(barcodes_GC$X4, barcode_count$X2, xlab="GC content", ylab="barcode count")
abline(lm(barcode_count$X2 ~ barcodes_GC$X4, weights=barcode_count$X2))
opar 
dev.off()

#regress.lm = lm(barcode_count$X2 ~ barcodes_GC$X4, weights=barcode_count$X2)
#summary(regress.lm)

new_frame <- barcodes_GC%>% filter(X4>0.3&X4<0.7)
nrow(new_frame)/96*100
#it means that 79.2% of the barcodes (76/96) have between a 30% and 70% GC content

#pearson coefficient of correlation
barcodes_GC_pearson <- wtd.cor(barcodes_GC$X4,barcode_count$X2,weight=barcode_count$X2)

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

##############

# Plot number of reads assigned to barcode vs 96 barcode (rank sorted)

#sort number of reads (ascending)
barcode_count_2 <- barcode_count$X2
barcode_count_3 <- data.frame(sort(barcode_count_2))

#add column for wells
barcode_count_3$X2 <- 1:nrow(barcode_count_3) 

pdf("read_count_per_barcode.pdf")
barplot(barcode_count_3$sort.barcode_count_2~barcode_count_3$X2, xlab="", ylab="", ylim=c(0,10000), xaxt='n', main="96-plex barcode counts") 
title(xlab="96 barcodes (rank-sorted)", line=1.0, cex.lab=1.5)
title(ylab="number of reads assigned to barcode", line=2.5, cex.lab=1.5)
grid(nx = NA, ny = 10, col = "gray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
dev.off()


########################################

# PHRED scores

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# install.packages("Rsubread")
library(Rsubread)

#here already swapped L1 L2 (incl edited nreads)
L2f <- qualityScores("~/Desktop/MG1655/quality_scores/forward_E.coli_1.fastq", input_format="FASTQ", offset=33, nreads=662835)
L2r <- qualityScores("~/Desktop/MG1655/quality_scores/reverse_E.coli_1.fastq", input_format="FASTQ", offset=33, nreads=662835)
L1f <- qualityScores("~/Desktop/MG1655/quality_scores/forward_E.coli_2.fastq", input_format="FASTQ", offset=33, nreads=550557)
L1r <- qualityScores("~/Desktop/MG1655/quality_scores/reverse_E.coli_2.fastq", input_format="FASTQ", offset=33, nreads=550557)
Hf <- qualityScores("~/Desktop/MG1655/quality_scores/Undetermined_S0_L001_R1_001.fastq", input_format="FASTQ", offset=33, nreads=611234)
Hr <- qualityScores("~/Desktop/MG1655/quality_scores/Undetermined_S0_L001_R2_001.fastq", input_format="FASTQ", offset=33, nreads=611234)

# subset to framents < 150 bp long (as Standard Flex 
# and 1:50 Flex are 300 bp paired end MiSeq runs)
L1f_150 <- L1f[,1:150]
L1r_150 <- L1r[,1:150]
L2f_150 <- L2f[,1:150]
L2r_150 <- L2r[,1:150]
Hf_150 <- Hf[,1:150]
Hr_150 <- Hr[,1:150]

# collect the PHRED means 
a <- colMeans(L1f_150, na.rm = FALSE, dims = 1)
b <- colMeans(L1r_150, na.rm = FALSE, dims = 1)
c <- colMeans(L2f_150, na.rm = FALSE, dims = 1)
d <- colMeans(L2r_150, na.rm = FALSE, dims = 1)
e <- colMeans(Hf_150, na.rm = FALSE, dims = 1)
f <- colMeans(Hr_150, na.rm = FALSE, dims = 1)

# make dataframe containing all the PHRED means 
x <- seq(1,150)
df <- data.frame(x, a, b, c, d, e, f)

# make dataframe longer (1 row : 1 value)
library(tidyverse)
df_long <- df %>%
  gather(libraries, value, a:f)

# give sensible names 
library(dplyr)
df_long2 <- df_long %>% 
  mutate(libraries = recode(libraries, `a` = "Standard Flex R1", `b` = "Standard Flex R2", `c` = "1:50 Flex R1", `d` = "1:50 Flex R2", `e` = "Hackflex R1", `f`= "Hackflex R2" ))

pdf("PHRED_scores_raw_reads.pdf")
require(ggplot2)
d <- ggplot(df_long2, aes(x=x, y=value, group=libraries, colour=libraries ) ) + 
  geom_line(size=1) + 
  xlab("read position (bp)") +
  ylab("average PHRED score")
e <- d + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("firebrick 4", "firebrick 3", "#006600", "#009900", "dodgerblue 4", "dodgerblue 3")) +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(0,25,50,75,100,125,150), lim = c(0, 150)) +
  scale_y_continuous(breaks = c(30,32,34,36,38), lim = c(29, 38))
f <- e + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=18))
f
dev.off()

###############################################################

# Coverage:

#Hackflex

library(readr)
nanop_vs_hackflex <- read_delim("nanop_vs_hackflex.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

dec <-read.table("nanop_vs_hackflex.tsv")
trimdec <- dec[200:4466382,]
par(mar=c(1,1,1,1))
kkk <- hist(trimdec$V2[trimdec$V2<100], breaks=100)
kk5 <- hist(trimdec$V2[trimdec$V2<=4], breaks=4, xlim=c(0,4))

#######

#Standard Flex

library("readr")
nanop_vs_nextera1 <- read_delim("nanop_vs_nextera2.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

ilm <-read.table("nanop_vs_nextera2.tsv")
trimilm <- ilm[200:4466382,]
par(mar=c(1,1,1,1))
ooo <- hist(trimilm$V2[trimilm$V2<100], breaks=100)
oo5 <- hist(trimilm$V2[trimilm$V2<=4], breaks=4, xlim=c(0,4))

#######

#1:50 Flex

library(readr)
nanop_vs_nextera2 <- read_delim("nanop_vs_nextera1.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

nop <-read.table("nanop_vs_nextera1.tsv")
trimnop <- nop[200:4466382,]
par(mar=c(1,1,1,1))
ppp <- hist(trimnop$V2[trimnop$V2<100], breaks=100)
pp5 <- hist(trimnop$V2[trimnop$V2<=4], breaks=4, xlim=c(0,4))

#######

#plot
pdf("coverage.pdf")
par(mfrow=c(3,1))
plot(ooo, xlim=c(0,100), ylim=c(0,400000), main="Standard Flex", xlab = "coverage per site", cex.lab = 1.5)
plot(ppp, xlim=c(0,100), ylim=c(0,400000), main="1:50 Flex", xlab = "coverage per site", cex.lab = 1.5)
plot(kkk, xlim=c(0,100), ylim=c(0,400000), main="Hackflex", xlab = "coverage per site", cex.lab = 1.5)
par(fig=c(0.61, 0.99, 0.75, 1), new = T) 
plot(oo5, xlim=c(0,4), ylim=c(0,4000), xlab = NULL, ylab = NULL, main = NULL, cex.axis = 0.8)
par(fig=c(0.61, 0.99, 0.4, 0.65), new = T) 
plot(pp5, xlim=c(0,4), ylim=c(0,4000), xlab = NULL, ylab = NULL, main = NULL, cex.axis = 0.8)
par(fig=c(0.61, 0.99, 0.1, 0.35), new = T) 
plot(kk5, xlim=c(0,4), ylim=c(0,4000), xlab = NULL, ylab = NULL, main = NULL, cex.axis = 0.8)
dev.off()

summary_L1 <- summary(trimilm$V2) 
summary_L2 <- summary(trimnop$V2) 
summary_H <- summary(trimdec$V2) 

library(weights)
L1_L2_pearson <- wtd.cor(trimilm$V2, trimnop$V2)
H_L1_pearson <- wtd.cor(trimdec$V2, trimilm$V2)
H_L2_pearson <- wtd.cor(trimdec$V2, trimnop$V2)

wtd.cor(trimilm$V2, trimnop$V2)
wtd.cor(trimdec$V2, trimilm$V2)
wtd.cor(trimdec$V2, trimnop$V2)

# get the non-zeros pvalues this way
cor.test(trimilm$V2, trimnop$V2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)
cor.test(trimdec$V2, trimilm$V2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
cor.test(trimdec$V2, trimnop$V2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

v_1 = paste0("Coverage:") 
v_2 = paste0("Standard Flex:") 
v_3 = paste0("Min: ", summary_L1[1]) 
v_4 = paste0 ("1st Qu: ", summary_L1[2]) 
v_5 = paste0("Median: ", summary_L1[3])
v_6 = paste0("Mean: ", summary_L1[4])
v_7 = paste0("3rd Qu: ", summary_L1[5])
v_8 = paste0("Max: ", summary_L1[6])
v_9 = paste0("standard deviation: ", sd(trimilm$V2))
v_10 = paste0("coefficient of variation: ", sd(trimilm$V2)/mean(trimilm$V2))
v_11 = paste0("1:50 Flex:")
v_12 = paste0("Min: ", summary_L2[1]) 
v_13 = paste0("1st Qu: ", summary_L2[2]) 
v_14 = paste0("Median: ", summary_L2[3]) 
v_15 = paste0("Mean: ", summary_L2[4]) 
v_16 = paste0("3rd Qu: ", summary_L2[5]) 
v_17 = paste0("Max: ", summary_L2[6]) 
v_18 = paste0("standard deviation: ", sd(trimnop$V2)) 
v_19 = paste0("coefficient of variation: ", sd(trimnop$V2)/mean(trimnop$V2)) 
v_20 = paste0("Hackflex: ") 
v_21 = paste0("Min: ", summary_H[1]) 
v_22 = paste0("1st Qu: ", summary_H[2]) 
v_23 = paste0("Median: ", summary_H[3]) 
v_24 = paste0("Mean: ", summary_H[4]) 
v_25 = paste0("3rd Qu: ", summary_H[5]) 
v_26 = paste0("Max: ", summary_H[6]) 
v_27 = paste0("standard deviation: ", sd(trimdec$V2)) 
v_28 = paste0("coefficient of variation: ", sd(trimdec$V2)/mean(trimdec$V2)) 
v_29 = paste0("rho Standard Flex - 1:50 Flex: ", L1_L2_pearson[1]) 
v_30 = paste0("rho Hackflex - Standard Flex: ", H_L1_pearson[1]) 
v_31 = paste0("rho Hackflex - 1:50 Flex: ", H_L2_pearson[1]) 
v_32 = paste0("p-value Standard Flex - 1:50 Flex: ", L1_L2_pearson[4]) 
v_33 = paste0("p-value Hackflex - Standard Flex: ", H_L1_pearson[4]) 
v_34 = paste0("p-value Hackflex - 1:50 Flex: ", H_L2_pearson[4]) 

# coverage summary: 
all <- c(v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13, v_14, v_15, v_16, v_17, v_18, v_19, v_20, v_21, v_22, v_23, v_24, v_25, v_26, v_27, v_28, v_29, v_30, v_31, v_32, v_33, v_34)
write.table(all, file = 'coverage_summary.dsv', sep = ":", row.names = FALSE, col.names = FALSE)
#try <- read.table("coverage_summary.dsv", sep=":")
#try

Low_cov_Hackflex <- subset(trimdec, trimdec$V2<3)
Low_cov_Standard_Flex <- subset(trimilm, trimilm$V2<3)
Low_cov_1_50_Flex <- subset(trimnop, trimnop$V2<3)

#add column to describe Library before joining datasets
Low_cov_Standard_Flex["Library"]<- NA
Low_cov_1_50_Flex["Library"]<- NA
Low_cov_Hackflex["Library"]<- NA
Low_cov_Standard_Flex[is.na(Low_cov_Standard_Flex)] <- "Standard_Flex"
Low_cov_1_50_Flex[is.na(Low_cov_1_50_Flex)] <- "1:50_Flex"
Low_cov_Hackflex[is.na(Low_cov_Hackflex)] <- "Hackflex"

#concatenate dataframes
Low_cov_all <- rbind(Low_cov_Standard_Flex, Low_cov_1_50_Flex, Low_cov_Hackflex)

#sensible heaades names 
names(Low_cov_all) <- c("Reference genome position", "Reads", "Library")

#save low coverage areas as a .csv
write.csv(Low_cov_all, file = "Low_cov_all.csv",row.names=FALSE)

########################################

# Fragment size:

library(readr)
frag_sizes_nextera1 <- read.table("frag_sizes_nextera2.txt", quote="\"", comment.char="")
frag_sizes_nextera2 <- read.table("frag_sizes_nextera1.txt", quote="\"", comment.char="")
frag_sizes_hackflex <- read.table("frag_sizes_hackflex.txt", quote="\"", comment.char="")

#subset to positive numbers only
pos_L1 <- subset(frag_sizes_nextera2, frag_sizes_nextera2$V1 > 0)
pos_L2 <- subset(frag_sizes_nextera1, frag_sizes_nextera1$V1 > 0)
pos_H <- subset(frag_sizes_hackflex, frag_sizes_hackflex$V1 > 0)

#add column with sequential numbers to each dataset
library(dplyr)
newL1 <- cbind( pos_L1, order=seq(nrow(pos_L1)) ) 
newL2 <- cbind( pos_L2, order=seq(nrow(pos_L2)) ) 
newH <- cbind( pos_H, order=seq(nrow(pos_H)) ) 

#join datasets (missing rows will be filled with NA)
L1L2<- dplyr::full_join(newL1, newL2, by = "order")
L1L2H <- dplyr::full_join(L1L2, newH, by = "order")

#rename columns
newcolnames_L1L2H <- L1L2H %>% 
  rename(
    L1 = V1.x,
    L2 = V1.y,
    H = V1
  )

#convert to long format (1 variable per row)
library("reshape2")
test_data_long <- melt(newcolnames_L1L2H, id="order") 

#replace variable names with names for plot
library(stringr)
test_data_long$variable <- str_replace_all(test_data_long$variable, "L1", "Standard Flex")
test_data_long$variable <- str_replace_all(test_data_long$variable, "L2", "1:50 Flex")
test_data_long$variable <- str_replace_all(test_data_long$variable, "H", "Hackflex")

#plot
library(ggplot2)
pdf('fragment_size.pdf')
libs <- c("Standard Flex", "1:50 Flex", "Hackflex")
ds <- subset(test_data_long, variable %in% libs & ! is.na(value))
p  <- ggplot(ds, aes(value, colour=variable, fill=variable)) + scale_x_continuous(limits=c(0,1000)) + xlab("fragment size (bp)")
q  <- p + geom_density(alpha=0.55)
r <- q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
s <- r + theme(legend.title = element_blank()) 
t <- s + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=14))
t
dev.off()


###############################################################


# Number of reads from .bam : cleaned reads 

library(readr)
RL_nanop_vs_nextera1 <- read_delim("RL_nanop_vs_nextera2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_nanop_vs_nextera2 <- read_delim("RL_nanop_vs_nextera1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_nanop_vs_hackflex <- read_delim("RL_nanop_vs_hackflex.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_Saureus_vs_nextera1 <- read_delim("RL_Saureus_vs_S.aureus1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_Saureus_vs_nextera2 <- read_delim("RL_Saureus_vs_S.aureus2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_PAO1_vs_nextera1 <- read_delim("RL_Paeruginosa_vs_P.aeruginosa1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_PAO1_vs_nextera2 <- read_delim("RL_Paeruginosa_vs_P.aeruginosa2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

RL_L1 <- sum(RL_nanop_vs_nextera1$Count)
RL_L2 <- sum(RL_nanop_vs_nextera2$Count)
RL_H <- sum(RL_nanop_vs_hackflex$Count)
RL_M1 <- sum(RL_Saureus_vs_nextera1$Count)
RL_M2 <- sum(RL_Saureus_vs_nextera2$Count)
RL_N1 <- sum(RL_PAO1_vs_nextera1$Count)
RL_N2 <- sum(RL_PAO1_vs_nextera2$Count)


###############################################################

# GC content vs coverage : 

library(readr)
GC_nanop_vs_nextera1 <- read_delim("GC_nanop_vs_nextera2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_nanop_vs_nextera2 <- read_delim("GC_nanop_vs_nextera1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_nanop_vs_hackflex <- read_delim("GC_nanop_vs_hackflex.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

library(tidyverse)
library(tidyr)

# GC content vs coverage : Standard Flex, 1:50 Flex, Hackflex 
GC_nanop_vs_nextera1$fractionOfReads <- (GC_nanop_vs_nextera1$fractionOfReads + 0.000000001)
L1 <- spread(GC_nanop_vs_nextera1, GCcontent, fractionOfReads)
rownames(L1) <- L1$Sample
L1$Library <- NULL
L1_transpose <- as.data.frame(t(L1[,-1]))
colnames(L1_transpose) <- L1$Sample
L1_transpose$diff_L1 <- (L1_transpose$nanop_vs_nextera2 / L1_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
L1_transpose$reads_L1 <- (L1_transpose$nanop_vs_nextera2 * RL_L1)
L1_transpose<-cbind(L1_transpose, GC_nanop_vs_nextera2[,3, drop=FALSE])

GC_nanop_vs_nextera2$fractionOfReads <- (GC_nanop_vs_nextera2$fractionOfReads + 0.000000001)
L2 <- spread(GC_nanop_vs_nextera2, GCcontent, fractionOfReads)
rownames(L2) <- L2$Sample
L2$Library <- NULL
L2_transpose <- as.data.frame(t(L2[,-1]))
colnames(L2_transpose) <- L2$Sample
L2_transpose$diff_L2 <- (L2_transpose$nanop_vs_nextera1 / L2_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
L2_transpose$reads_L2 <- (L2_transpose$nanop_vs_nextera1 * RL_L2)
L2_transpose<-cbind(L2_transpose, GC_nanop_vs_nextera1[,3, drop=FALSE])

GC_nanop_vs_hackflex$fractionOfReads <- (GC_nanop_vs_hackflex$fractionOfReads + 0.000000001)
H <- spread(GC_nanop_vs_hackflex, GCcontent, fractionOfReads)
rownames(H) <- H$Sample
H$Library <- NULL
H_transpose <- as.data.frame(t(H[,-1]))
colnames(H_transpose) <- H$Sample
H_transpose$diff_H <- (H_transpose$nanop_vs_hackflex / H_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
H_transpose$reads_H <- (H_transpose$nanop_vs_hackflex * RL_H)
H_transpose<-cbind(H_transpose, GC_nanop_vs_hackflex[,3, drop=FALSE])

library(weights)
wtd.cor(L1_transpose$GCcontent,L1_transpose$diff,weight=L1_transpose$reads_L1)
wtd.cor(L2_transpose$GCcontent,L2_transpose$diff,weight=L2_transpose$reads_L2)
wtd.cor(H_transpose$GCcontent,H_transpose$diff,weight=H_transpose$reads_H)

#GC vs coverage plot
pdf("GC_vs_coverage.pdf")
plot(L1_transpose$GCcontent[L1_transpose$diff!=1], L1_transpose$diff[L1_transpose$diff!=1], xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(L2_transpose$GCcontent[L2_transpose$diff!=1], L2_transpose$diff[L2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(H_transpose$GCcontent[H_transpose$diff!=1], H_transpose$diff[H_transpose$diff!=1], col="green", cex = 0.5, pch = 1)
abline(lm(L1_transpose$diff ~ L1_transpose$GCcontent, weights=L1_transpose$reads_L1), col="blue")
abline(lm(L2_transpose$diff ~ L2_transpose$GCcontent, weights=L2_transpose$reads_L2), col="red")
abline(lm(H_transpose$diff ~ H_transpose$GCcontent, weights=H_transpose$reads_H), col="green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.9502)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.9587)),
                  expression(paste("Hackflex            ", rho, " = ", -0.7696))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()


# GC content vs coverage : Standard Flex and 1:50 Flex with S. aureus

library(readr)
GC_Saureus_vs_nextera1 <- read_delim("GC_Saureus_vs_S.aureus1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_Saureus_vs_nextera2 <- read_delim("GC_Saureus_vs_S.aureus2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

library(tidyverse)

GC_Saureus_vs_nextera1$fractionOfReads <- (GC_Saureus_vs_nextera1$fractionOfReads + 0.000000001)
M1 <- spread(GC_Saureus_vs_nextera1, GCcontent, fractionOfReads)
rownames(M1) <- M1$Sample
M1$Library <- NULL
M1_transpose <- as.data.frame(t(M1[,-1]))
colnames(M1_transpose) <- M1$Sample
M1_transpose$diff_M1 <- (M1_transpose$Saureus_vs_S.aureus1 / M1_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
M1_transpose$reads_M1 <- (M1_transpose$Saureus_vs_S.aureus1 * RL_M1)
M1_transpose<-cbind(M1_transpose, GC_Saureus_vs_nextera1[,3, drop=FALSE])

GC_Saureus_vs_nextera2$fractionOfReads <- (GC_Saureus_vs_nextera2$fractionOfReads + 0.000000001)
M2 <- spread(GC_Saureus_vs_nextera2, GCcontent, fractionOfReads)
rownames(M2) <- M2$Sample
M2$Library <- NULL
M2_transpose <- as.data.frame(t(M2[,-1]))
colnames(M2_transpose) <- M2$Sample
M2_transpose$diff_M2 <- (M2_transpose$Saureus_vs_S.aureus2 / M2_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
M2_transpose$reads_M2 <- (M2_transpose$Saureus_vs_S.aureus2 * RL_M2)
M2_transpose<-cbind(M2_transpose, GC_Saureus_vs_nextera2[,3, drop=FALSE])

library(weights)
wtd.cor(M1_transpose$GCcontent,M1_transpose$diff,weight=M1_transpose$reads_M1)
wtd.cor(M2_transpose$GCcontent,M2_transpose$diff,weight=M2_transpose$reads_M2)

#GC vs coverage plot
pdf("GC_vs_coverage_Saureus.pdf")
plot(M1_transpose$GCcontent[M1_transpose$diff!=1], M1_transpose$diff[M1_transpose$diff!=1], xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(M2_transpose$GCcontent[M2_transpose$diff!=1], M2_transpose$diff[M2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
abline(lm(M1_transpose$diff ~ M1_transpose$GCcontent, weights=M1_transpose$reads_M1), col="blue")
abline(lm(M2_transpose$diff ~ M2_transpose$GCcontent, weights=M2_transpose$reads_M2), col="red")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.9234)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.9398))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()


# GC content vs coverage : Standard Flex and 1:50 Flex with PAO1

library(readr)
GC_PAO_vs_nextera1 <- read_delim("GC_Paeruginosa_vs_P.aeruginosa1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_PAO_vs_nextera2 <- read_delim("GC_Paeruginosa_vs_P.aeruginosa2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

GC_PAO_vs_nextera1$fractionOfReads <- (GC_PAO_vs_nextera1$fractionOfReads + 0.000000001)
N1 <- spread(GC_PAO_vs_nextera1, GCcontent, fractionOfReads)
rownames(N1) <- N1$Sample
N1$Library <- NULL
N1_transpose <- as.data.frame(t(N1[,-1]))
colnames(N1_transpose) <- N1$Sample
N1_transpose$diff_N1 <- (N1_transpose$Paeruginosa_vs_P.aeruginosa1 / N1_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
N1_transpose$reads_N1 <- (N1_transpose$Paeruginosa_vs_P.aeruginosa1 * RL_N1)
N1_transpose<-cbind(N1_transpose, GC_PAO_vs_nextera1[,3, drop=FALSE])

GC_PAO_vs_nextera2$fractionOfReads <- (GC_PAO_vs_nextera2$fractionOfReads + 0.000000001)
N2 <- spread(GC_PAO_vs_nextera2, GCcontent, fractionOfReads)
rownames(N2) <- N2$Sample
N2$Library <- NULL
N2_transpose <- as.data.frame(t(N2[,-1]))
colnames(N2_transpose) <- N2$Sample
N2_transpose$diff_N2 <- (N2_transpose$Paeruginosa_vs_P.aeruginosa2 / N2_transpose$Reference)
# read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
N2_transpose$reads_N2 <- (N2_transpose$Paeruginosa_vs_P.aeruginosa2 * RL_N2)
N2_transpose<-cbind(N2_transpose, GC_PAO_vs_nextera2[,3, drop=FALSE])

library(weights)
wtd.cor(N1_transpose$GCcontent,N1_transpose$diff,weight=N1_transpose$reads_N1)
wtd.cor(N2_transpose$GCcontent,N2_transpose$diff,weight=N2_transpose$reads_N2)

#GC vs coverage plot
pdf("GC_vs_coverage_Paeruginosa.pdf")
plot(N1_transpose$GCcontent[N1_transpose$diff!=1], N1_transpose$diff[N1_transpose$diff!=1], xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(N2_transpose$GCcontent[N2_transpose$diff!=1], N2_transpose$diff[N2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
  abline(lm(N1_transpose$diff ~ N1_transpose$GCcontent, weights=N1_transpose$reads_N1), col="blue")
abline(lm(N2_transpose$diff ~ N2_transpose$GCcontent, weights=N2_transpose$reads_N2), col="red")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.9916)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.9903))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()



pdf("Standard_Flex_vs_1_50_Flex_3_species.pdf")
par(mfrow=c(3,1))
plot(M1_transpose$GCcontent[M1_transpose$diff!=1], M1_transpose$diff[M1_transpose$diff!=1], main="S. aureus", xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab="obs/exp reads", cex.lab=1.5) +
  points(M2_transpose$GCcontent[M2_transpose$diff!=1], M2_transpose$diff[M2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
  abline(lm(M1_transpose$diff ~ M1_transpose$GCcontent, weights=M1_transpose$reads_M1), col="blue")
legend("topright", 
       legend = c("Standard Flex", "1:50 Flex", "Hackflex"), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
abline(lm(M2_transpose$diff ~ M2_transpose$GCcontent, weights=M2_transpose$reads_M2), col="red")
plot(L1_transpose$GCcontent[L1_transpose$diff!=1], L1_transpose$diff[L1_transpose$diff!=1], main="E. coli", xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab="obs/exp reads", cex.lab=1.5) +
  points(L2_transpose$GCcontent[L2_transpose$diff!=1], L2_transpose$diff[L2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(H_transpose$GCcontent[H_transpose$diff!=1], H_transpose$diff[H_transpose$diff!=1], col="green", cex = 0.5, pch = 1)
abline(lm(L1_transpose$diff ~ L1_transpose$GCcontent, weights=L1_transpose$reads_L1), col="blue")
abline(lm(L2_transpose$diff ~ L2_transpose$GCcontent, weights=L2_transpose$reads_L2), col="red")
abline(lm(H_transpose$diff ~ H_transpose$GCcontent, weights=H_transpose$reads_H), col="green")
legend("topright", 
       legend = c("Standard Flex", "1:50 Flex", "Hackflex"), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
plot(N1_transpose$GCcontent[N1_transpose$diff!=1], N1_transpose$diff[N1_transpose$diff!=1], main="P. aeruginosa", xlim=c(0.1,0.9), ylim=c(0.3,2.5), type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab="obs/exp reads", cex.lab=1.5) +
  points(N2_transpose$GCcontent[N2_transpose$diff!=1], N2_transpose$diff[N2_transpose$diff!=1], col="red", cex = 0.5, pch = 2) +
  abline(lm(N1_transpose$diff ~ N1_transpose$GCcontent, weights=N1_transpose$reads_N1), col="blue")
abline(lm(N2_transpose$diff ~ N2_transpose$GCcontent, weights=N2_transpose$reads_N2), col="red")
legend("topright", 
       legend = c("Standard Flex", "1:50 Flex", "Hackflex"), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()


###############################################################

# ALFRED summary: 

setwd("~/Desktop/MG1655/MG1655_final")

library(readr)
ME_nanop_vs_nextera1 <- read_delim("ME_nanop_vs_nextera2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ME_nanop_vs_nextera2 <- read_delim("ME_nanop_vs_nextera1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
ME_nanop_vs_hackflex <- read_delim("ME_nanop_vs_hackflex.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

# join datasets
ME_all <- rbind(ME_nanop_vs_nextera1, ME_nanop_vs_nextera2)
ME_all <- rbind(ME_all, ME_nanop_vs_hackflex)

# transpose and keep headers
ME_all_transposed <- as.data.frame(t(ME_all[,-1]))
colnames(ME_all_transposed) <- ME_all$Sample

#remove row "Library"
ME_all_transposed <- ME_all_transposed[-c(1), ] 

#rename headers
names(ME_all_transposed) <- c("Standard Flex", "1:50 Flex", "Hackflex")

View(ME_all_transposed)

install.packages("formattable")
install.packages("htmltools")
install.packages("webshot")
library("formattable")
library("htmltools")
library("webshot")
webshot::install_phantomjs()

# function to export
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

FT <- formattable(ME_all_transposed)
export_formattable(FT,"ALFRED_summary.pdf")

