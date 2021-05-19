
setwd("~/Desktop/MG1655/hackflex_libs/")
########################################

# PHRED scores

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# install.packages("Rsubread")
library(Rsubread)

L1f <- qualityScores("K1_R1.fastq", input_format="FASTQ", offset=33, nreads=919146)
L1r <- qualityScores("K1_R2.fastq", input_format="FASTQ", offset=33, nreads=919146)
L2f <- qualityScores("K3_R1.fastq", input_format="FASTQ", offset=33, nreads=977948)
L2r <- qualityScores("K3_R2.fastq", input_format="FASTQ", offset=33, nreads=977948)
Hf <- qualityScores("K9_R1.fastq", input_format="FASTQ", offset=33, nreads=954935)
Hr <- qualityScores("K9_R2.fastq", input_format="FASTQ", offset=33, nreads=954935)

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
getwd()
pdf("PHRED_scores_raw_reads_new.pdf")
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
###############################################################

# Coverage:

#Hackflex

library(readr)
nanop_vs_hackflex <- read_delim("nanop_vs_K9.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

dec <-read.table("nanop_vs_K9.tsv")
trimdec <- dec[200:4466382,] # 4466382 is the E. coli genome length
par(mar=c(1,1,1,1))
kkk <- hist(trimdec$V2[trimdec$V2<100], breaks=100)
kk5 <- hist(trimdec$V2[trimdec$V2<=4], breaks=4, xlim=c(0,4))

#######

#Standard Flex

library("readr")
nanop_vs_nextera1 <- read_delim("nanop_vs_K1.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

ilm <-read.table("nanop_vs_K1.tsv")
trimilm <- ilm[200:4466382,] # 4466382 is the E. coli genome length
par(mar=c(1,1,1,1))
ooo <- hist(trimilm$V2[trimilm$V2<100], breaks=100)
oo5 <- hist(trimilm$V2[trimilm$V2<=4], breaks=4, xlim=c(0,4))

#######

#1:50 Flex

library(readr)
nanop_vs_nextera2 <- read_delim("nanop_vs_K3.tsv", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

nop <-read.table("nanop_vs_K3.tsv")
trimnop <- nop[200:4466382,] # 4466382 is the E. coli genome length
par(mar=c(1,1,1,1))
ppp <- hist(trimnop$V2[trimnop$V2<100], breaks=100)
pp5 <- hist(trimnop$V2[trimnop$V2<=4], breaks=4, xlim=c(0,4))

#######

#plot
pdf("coverage_new.pdf")
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
write.table(all, file = 'coverage_summary_new.csv', sep = ":", row.names = FALSE, col.names = FALSE)
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
write.csv(Low_cov_all, file = "Low_cov_all_new.csv",row.names=FALSE)




# Number of reads from .bam : cleaned reads 

library(readr)
RL_EC_vs_SF <- read_delim("RL_nanop_vs_K1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_EC_vs_SF150 <- read_delim("RL_nanop_vs_K3.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_EC_vs_HF <- read_delim("RL_nanop_vs_K9.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

RL_SA_vs_SF <- read_delim("~/Desktop/MG1655/MG1655_final/RL_Saureus_vs_S.aureus1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_SA_vs_SF150 <- read_delim("~/Desktop/MG1655/MG1655_final/RL_Saureus_vs_S.aureus2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_SA_vs_HF_rep1 <- read_delim("RL_SA_vs_K7.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_SA_vs_HF_rep2 <- read_delim("RL_SA_vs_K10.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

RL_PA_vs_SF <- read_delim("~/Desktop/MG1655/MG1655_final/RL_Paeruginosa_vs_P.aeruginosa1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_PA_vs_SF150 <- read_delim("~/Desktop/MG1655/MG1655_final/RL_Paeruginosa_vs_P.aeruginosa2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_PA_vs_HF_rep1 <- read_delim("RL_PAO_vs_K8.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
RL_PA_vs_HF_rep2 <- read_delim("RL_PAO_vs_K11.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


###############################################################
###############################################################

# GC content vs coverage : 

library(readr)
GC_EC_vs_SF <- read_delim("GC_nanop_vs_K1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_EC_vs_SF150 <- read_delim("GC_nanop_vs_K3.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_EC_vs_HF <- read_delim("GC_nanop_vs_K9.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

GC_SA_vs_SF <- read_delim("~/Desktop/MG1655/MG1655_final/GC_Saureus_vs_S.aureus1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_SA_vs_SF150 <- read_delim("~/Desktop/MG1655/MG1655_final/GC_Saureus_vs_S.aureus2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_SA_vs_HF_rep1 <- read_delim("GC_SA_vs_K7.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_SA_vs_HF_rep2 <- read_delim("GC_SA_vs_K10.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

GC_PA_vs_SF <- read_delim("~/Desktop/MG1655/MG1655_final/GC_Paeruginosa_vs_P.aeruginosa1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_PA_vs_SF150 <- read_delim("~/Desktop/MG1655/MG1655_final/GC_Paeruginosa_vs_P.aeruginosa2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_PA_vs_HF_rep1 <- read_delim("GC_PAO_vs_K8.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
GC_PA_vs_HF_rep2 <- read_delim("GC_PAO_vs_K11.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


library(tidyverse)
library(tidyr)

# GC content vs coverage : 
# takes as input: GC and RL outputs from alfred
# returns df to plot GC content diff with reference genome
myfun <- function(GC_df,RL_df) {
  
  GC_df <- as.data.frame(GC_df)
  GC_df$fractionOfReads <- (GC_df$fractionOfReads + 0.000000001)
  DF <- spread(GC_df, GCcontent, fractionOfReads)
  DF <- as.data.frame(DF)
  rownames(DF) <- DF$Sample
  DF$Library <- NULL
  DF_transpose <- as.data.frame(t(DF[,-1]))
  colnames(DF_transpose) <- DF$Sample
  DF_transpose$diff <- (DF_transpose[,1]/DF_transpose[,2])
  # read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
  DF_transpose$reads <- (DF_transpose[,1] * sum(RL_df$Count))
  DF_transpose<-cbind(DF_transpose, GC_df[,3, drop=FALSE])
  
  return(DF_transpose)
}

# # testing function: 
# GC_df <- as.data.frame(GC_EC_vs_SF)
# GC_df$fractionOfReads <- (GC_df$fractionOfReads + 0.000000001)
# DF <- spread(GC_df, GCcontent, fractionOfReads)
# DF <- as.data.frame(DF)
# rownames(DF) <- DF$Sample
# DF$Library <- NULL
# DF_transpose <- as.data.frame(t(DF[,-1]))
# #colnames(DF_transpose) <- DF$Sample
# DF_transpose$diff <- (DF_transpose[,1]-DF_transpose[,2])
# # read numbers here based on : tot. number of mapped reads (as per samtools flagstat) minus eliminated duplicates
# DF_transpose$reads <- (DF_transpose[,1] * sum(RL_EC_vs_SF$Count))
# DF_transpose<-cbind(DF_transpose, GC_df[,3, drop=FALSE])
# 

GC_EC_SF <- myfun(GC_EC_vs_SF,RL_EC_vs_SF)   
GC_EC_SF150 <- myfun(GC_EC_vs_SF150,RL_EC_vs_SF150)   
GC_EC_HF <- myfun(GC_EC_vs_SF,RL_EC_vs_HF)   

GC_SA_SF <- myfun(GC_SA_vs_SF,RL_SA_vs_SF)   
GC_SA_SF150 <- myfun(GC_SA_vs_SF150,RL_SA_vs_SF150)   
GC_SA_HF_rep1 <- myfun(GC_SA_vs_HF_rep1,RL_SA_vs_HF_rep1)   
GC_SA_HF_rep2 <- myfun(GC_SA_vs_HF_rep2,RL_SA_vs_HF_rep2)   

GC_PA_SF <- myfun(GC_PA_vs_SF,RL_PA_vs_SF)   
GC_PA_SF150 <- myfun(GC_PA_vs_SF150,RL_PA_vs_SF150)   
GC_PA_HF_rep1 <- myfun(GC_PA_vs_HF_rep1,RL_PA_vs_HF_rep1)   
GC_PA_HF_rep2 <- myfun(GC_PA_vs_HF_rep2,RL_PA_vs_HF_rep2)     


library(weights)
wtd.cor(GC_EC_SF$GCcontent,GC_EC_SF$diff,weight=GC_EC_SF$reads)[1,1]
wtd.cor(GC_EC_SF150$GCcontent,GC_EC_SF150$diff,weight=GC_EC_SF150$reads)[1,1]
wtd.cor(GC_EC_HF$GCcontent,GC_EC_HF$diff,weight=GC_EC_HF$reads)[1,1]

wtd.cor(GC_SA_SF$GCcontent,GC_SA_SF$diff,weight=GC_SA_SF$reads)[1,1]
wtd.cor(GC_SA_SF150$GCcontent,GC_SA_SF150$diff,weight=GC_SA_SF150$reads)[1,1]
wtd.cor(GC_SA_HF_rep1$GCcontent,GC_SA_HF_rep1$diff,weight=GC_SA_HF_rep1$reads)[1,1]
wtd.cor(GC_SA_HF_rep2$GCcontent,GC_SA_HF_rep2$diff,weight=GC_SA_HF_rep2$reads)[1,1]


wtd.cor(GC_PA_SF$GCcontent,GC_PA_SF$diff,weight=GC_PA_SF$reads)[1,1]
wtd.cor(GC_PA_SF150$GCcontent,GC_PA_SF150$diff,weight=GC_PA_SF150$reads)[1,1]
wtd.cor(GC_PA_HF_rep1$GCcontent,GC_PA_HF_rep1$diff,weight=GC_PA_HF_rep1$reads)[1,1]
wtd.cor(GC_PA_HF_rep2$GCcontent,GC_PA_HF_rep2$diff,weight=GC_PA_HF_rep2$reads)[1,1]



#GC vs coverage plot
pdf("GC_vs_coverage_new.pdf")
plot(GC_EC_SF$GCcontent[GC_EC_SF$diff!=1], GC_EC_SF$diff[GC_EC_SF$diff!=1], 
     xlim=c(0.1,0.9), 
     #ylim=c(0.3,2.5), 
     type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expected reads", cex.lab=1.5) +
  points(GC_EC_SF150$GCcontent[GC_EC_SF150$diff!=1], GC_EC_SF150$diff[GC_EC_SF150$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(GC_EC_HF$GCcontent[GC_EC_HF$diff!=1], GC_EC_HF$diff[GC_EC_HF$diff!=1], col="green", cex = 0.5, pch = 1)
abline(lm(GC_EC_SF$diff ~ GC_EC_SF$GCcontent, weights=GC_EC_SF$reads), col="blue")
abline(lm(GC_EC_SF150$diff ~ GC_EC_SF150$GCcontent, weights=GC_EC_SF150$reads), col="red")
abline(lm(GC_EC_HF$diff ~ GC_EC_HF$GCcontent, weights=GC_EC_HF$reads), col="green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.53)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.57)),
                  expression(paste("Hackflex            ", rho, " = ", -0.53))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
plot(GC_SA_SF$GCcontent[GC_SA_SF$diff!=1], GC_SA_SF$diff[GC_SA_SF$diff!=1], 
     xlim=c(0.1,0.9), 
     #ylim=c(0.5,1.5), 
     type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expSAted reads", cex.lab=1.5) +
  points(GC_SA_SF150$GCcontent[GC_SA_SF150$diff!=1], GC_SA_SF150$diff[GC_SA_SF150$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(GC_SA_HF_rep1$GCcontent[GC_SA_HF_rep1$diff!=1], GC_SA_HF_rep1$diff[GC_SA_HF_rep1$diff!=1], col="green", cex = 0.5, pch = 1) +
  points(GC_SA_HF_rep2$GCcontent[GC_SA_HF_rep2$diff!=1], GC_SA_HF_rep2$diff[GC_SA_HF_rep2$diff!=1], col="dark green", cex = 0.5, pch = 1)
abline(lm(GC_SA_SF$diff ~ GC_SA_SF$GCcontent, weights=GC_SA_SF$reads), col="blue")
abline(lm(GC_SA_SF150$diff ~ GC_SA_SF150$GCcontent, weights=GC_SA_SF150$reads), col="red")
abline(lm(GC_SA_HF_rep1$diff ~ GC_SA_HF_rep1$GCcontent, weights=GC_SA_HF_rep1$reads), col="green")
abline(lm(GC_SA_HF_rep2$diff ~ GC_SA_HF_rep2$GCcontent, weights=GC_SA_HF_rep2$reads), col="dark green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.49)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.47)),
                  expression(paste("Hackflex            ", rho, " = ", -0.49))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
plot(GC_PA_SF$GCcontent[GC_PA_SF$diff!=1], GC_PA_SF$diff[GC_PA_SF$diff!=1], 
     xlim=c(0.1,0.9), 
     #ylim=c(0.5,1.5), 
     type="p",col="blue",cex = 0.5, pch = 0, xlab="GC content of region", ylab=" ratio observed/expPAted reads", cex.lab=1.5) +
  points(GC_PA_SF150$GCcontent[GC_PA_SF150$diff!=1], GC_PA_SF150$diff[GC_PA_SF150$diff!=1], col="red", cex = 0.5, pch = 2) +
  points(GC_PA_HF_rep1$GCcontent[GC_PA_HF_rep1$diff!=1], GC_PA_HF_rep1$diff[GC_PA_HF_rep1$diff!=1], col="green", cex = 0.5, pch = 1) +
  points(GC_PA_HF_rep2$GCcontent[GC_PA_HF_rep2$diff!=1], GC_PA_HF_rep2$diff[GC_PA_HF_rep2$diff!=1], col="dark green", cex = 0.5, pch = 1)
abline(lm(GC_PA_SF$diff ~ GC_PA_SF$GCcontent, weights=GC_PA_SF$reads), col="blue")
abline(lm(GC_PA_SF150$diff ~ GC_PA_SF150$GCcontent, weights=GC_PA_SF150$reads), col="red")
abline(lm(GC_PA_HF_rep1$diff ~ GC_PA_HF_rep1$GCcontent, weights=GC_PA_HF_rep1$reads), col="green")
abline(lm(GC_PA_HF_rep2$diff ~ GC_PA_HF_rep2$GCcontent, weights=GC_PA_HF_rep2$reads), col="dark green")
legend("topright", 
       legend = c(expression(paste("Standard Flex   ", rho, " = ", -0.49)),
                  expression(paste("1:50 Flex           ", rho, " = ", -0.47)),
                  expression(paste("Hackflex            ", rho, " = ", -0.49))), 
       pch = c(0,2,1), 
       col = c("blue", "red", "green"), 
       bty = "n", 
       pt.cex = 1, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = 0.01)
dev.off()



