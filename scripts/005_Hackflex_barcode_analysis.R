

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
library(Rmisc)


# set directories and select samples: 

barcode_libs <- "~/Desktop/MG1655/goal_barcode/"
barcode_source_data <- "~/Desktop/MG1655/goal_barcode/source_data/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"


########################################
########################################

# barcodes demux


demux_v0 <- read_csv(paste0(barcode_source_data, "DemultiplexingStats_v0.xml"), col_names = FALSE)
demux_v1 <- read_csv(paste0(barcode_source_data, "DemultiplexingStats_v1.xml"), col_names = FALSE)

# function to return details of read counts
give_counts <- function(demux) {
  
  sample_name <- demux[demux$X1 %like% "Sample name", ] %>% distinct() 
  counts <- demux[demux$X1 %like% "<BarcodeCount>", ] %>%
    distinct()
  
  df <- cbind(sample_name,counts)
  colnames(df) <- c("sample_name","counts")
  
  df$counts <- gsub("</BarcodeCount>","",df$counts)
  df$counts <- gsub("<BarcodeCount>","",df$counts)
  df$counts <- as.numeric(df$counts)
  df$sample_name <- gsub("<Sample name=","",df$sample_name)
  df$sample_name <- gsub(">","",df$sample_name)
  
  # we multiply by 2 because each entry is repeated twice (R1 and R2) and we use distinct() above. 
  df <- df %>% dplyr::mutate(counts=counts*2)

  return(df)
  
}



counts_v0 <- give_counts(demux_v0)
counts_v1 <- give_counts(demux_v1)
View(counts_v0)

#####

# barcode contamination : 
undet <- counts_v0 %>% dplyr::filter(sample_name == "\"Undetermined\"")
all <- counts_v0 %>% dplyr::filter(sample_name == "\"all\"")
barcode_contamination_rate_v0 <- undet$counts/all$counts
sample_misassignment_rate_v0 <- barcode_contamination_rate_v0^2


undet <- counts_v1 %>% dplyr::filter(sample_name == "\"Undetermined\"")
all <- counts_v1 %>% dplyr::filter(sample_name == "\"all\"")
barcode_contamination_rate_v1 <- undet$counts/all$counts
sample_misassignment_rate_v1 <- barcode_contamination_rate_v1^2

#####

# keep only barcode libs

# v0
counts_v0 <- counts_v0 %>% dplyr::filter(!sample_name == "\"Undetermined\"" &
                                           !sample_name == "\"all\"")

# v1
to_remove <- rbind(counts_v1[98:126,1], # these are the rest of the libs on that sequencing run, and the "all" 
                   counts_v1[1,1]) # this is the undetermined
counts_v1 <- counts_v1 %>%
  dplyr::filter(str_detect(sample_name, "HF_barcode"))


#####


# function to return details of read counts
give_deets <- function(counts) {
  
  # set minimum: a fourth of the read count of the lower confidence interval (CI set at 0.999)  
  CIs <- CI(counts$counts,ci = 0.999)
  mymin <- CIs[3]/4
  
  mymin
  sort(counts$counts)

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
  counts_a <- counts_transposed %>%
    dplyr::select(summary, everything())
  counts_a$selection <- "all"
  
  b <- counts %>%
    dplyr::filter(counts>mymin) %>%
    dplyr::summarise(barcodes_combos=n(),
                     Min=min(counts),
                     Mean=mean(counts),
                     Median=median(counts),
                     Sd=sd(counts),
                     Max=max(counts),
                     Sum=sum(counts),
                     Coeff.variation=sd(counts)/mean(counts))
  counts_transposed <- as.data.frame(t(b), make.names = TRUE)
  counts_transposed$summary <- rownames(counts_transposed)
  colnames(counts_transposed) <- c("value","summary")
  counts_b <- counts_transposed %>%
    dplyr::select(summary, everything())
  counts_b$selection <- "without outliers"
  
  fin <- rbind(counts_a, counts_b)
  
  return(fin)
}


deets_v0 <- give_deets(counts_v0)
deets_v1 <- give_deets(counts_v1)

deets_v0 <- deets_v0 %>% dplyr::mutate(value=round(value,3))
deets_v1 <- deets_v1 %>% dplyr::mutate(value=round(value,3))

fwrite(x=deets_v0, file=paste0(barcode_libs,"barcodes_v0_summary.csv"))
fwrite(x=deets_v1, file=paste0(barcode_libs,"barcodes_v1_summary.csv"))


pdf(paste0(barcode_libs,"barcode_v0_distribution.pdf"))
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(counts_v0$counts, 
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
hist(counts_v0$counts, 
     breaks = seq(from=1, to=20000, by=500),
     main = NULL, xlab = "barcode count", ylab = "Frequency")
dev.off()


pdf(paste0(barcode_libs,"barcode_v1_distribution.pdf"))
# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram 
par(mar=c(0, 6, 1.1, 2))
boxplot(counts_v1$counts, ylim=c(0,400000), 
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
hist(counts_v1$counts, xlim=c(0,400000), 
     breaks = seq(from=1, to=400000, by=10000),
     main = NULL, xlab = "barcode count", ylab = "Frequency")
dev.off()



########################################
########################################

# GC content of oligos: 

complete_barcodes_v0 <- read_csv(file.path(barcode_source_data,"complete_barcodes_v0.csv"))
complete_barcodes_v1 <- read_csv(file.path(barcode_source_data,"complete_barcodes_v1.csv"))


give_GC_oligo <- function(complete_barcodes) {
  
  complete_barcodes <- cSplit(complete_barcodes, "BarcodeDesign_Well_Oligo", "_")
  complete_barcodes$oligo <- NULL
  colnames(complete_barcodes) <- c("Barcode","F","N","entire_barcode","barcode_design","plate_well","oligo")
  
  complete_barcodes <- complete_barcodes %>%
    distinct() %>%
    dplyr::mutate(Gs = str_count(Barcode, "G"),
                  Cs = str_count(Barcode, "C"),
                  GC_content = ((Gs+Cs) / (str_length(Barcode))*100))
  
  return(complete_barcodes)
}

give_GC_entire_barcode <- function(complete_barcodes) {
  
  complete_barcodes <- cSplit(complete_barcodes, "BarcodeDesign_Well_Oligo", "_")
  complete_barcodes$oligo <- NULL
  colnames(complete_barcodes) <- c("Barcode","F","N","entire_barcode","barcode_design","plate_well","oligo")
  
  complete_barcodes <- complete_barcodes %>%
    distinct() %>%
    dplyr::mutate(Gs = str_count(entire_barcode, "G"),
                  Cs = str_count(entire_barcode, "C"),
                  GC_content = ((Gs+Cs) / (str_length(entire_barcode))*100))
  
  return(complete_barcodes)
  
}


# GC content of oligos

v0_GC_oligo <- give_GC_oligo(complete_barcodes_v0)
summary(v0_GC_oligo$GC_content)
v0_GC_oligo %>% dplyr::filter(GC_content>80)
v0_GC_oligo %>% dplyr::filter(GC_content<10)


v1_GC_oligo <- give_GC_oligo(complete_barcodes_v1)
summary(v1_GC_oligo$GC_content)
v1_GC_oligo %>% dplyr::filter(GC_content>80)
v1_GC_oligo %>% dplyr::filter(GC_content<10)



# GC content of entire barcodes 

v0_GC_entire_barcode <- give_GC_entire_barcode(complete_barcodes_v0)
summary(v0_GC_entire_barcode$GC_content)
v0_GC_entire_barcode %>% dplyr::filter(GC_content>80)
v0_GC_entire_barcode %>% dplyr::filter(GC_content<10)


v1_GC_entire_barcode <- give_GC_entire_barcode(complete_barcodes_v1)
summary(v1_GC_entire_barcode$GC_content)
v1_GC_entire_barcode %>% dplyr::filter(GC_content>80)
v1_GC_entire_barcode %>% dplyr::filter(GC_content<10)



# GC content of entire barcodes - correlation with read counts: 

v0_df <- v0_GC_entire_barcode %>% 
  dplyr::select(plate_well,oligo, GC_content) %>%
  pivot_wider(names_from=c(oligo), values_from=GC_content)
GC_bias_v0 <- cbind(v0_df,counts_v0)
GC_v0_plot <- GC_bias_v0 %>%
  pivot_longer(cols=c(i5,i7)) %>%
  ggplot(., aes(x=value,y=counts,color=name))+
  geom_point(size=0.3)+
  geom_smooth(size=0.5)+
  labs(x="Barcode GC content (%)",
       y="Barcode counts",
       color="oligo")+
  theme_bw()


v1_df <- v1_GC_entire_barcode %>% 
  dplyr::select(plate_well,oligo, GC_content) %>%
  pivot_wider(names_from=c(oligo), values_from=GC_content)
GC_bias_v1 <- cbind(v1_df,counts_v1)
GC_v1_plot <- GC_bias_v1 %>%
  pivot_longer(cols=c(i5,i7)) %>%
  ggplot(., aes(x=value,y=counts,color=name))+
  geom_point(size=0.3)+
  geom_smooth(size=0.5)+
  labs(x="Barcode GC content (%)",
       y="Barcode counts",
       color="oligo")+
  theme_bw()

# Pearson coefficient of correlation - v0
i5_GC_rho_v0 <- as.data.frame(wtd.cor(GC_bias_v0$i5,GC_bias_v0$counts,weight=GC_bias_v0$counts)) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))
i7_GC_rho_v0 <- as.data.frame(wtd.cor(GC_bias_v0$i7,GC_bias_v0$counts,weight=GC_bias_v0$counts)) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))


# Pearson coefficient of correlation - v1
i5_GC_rho_v1 <- as.data.frame(wtd.cor(GC_bias_v1$i5,GC_bias_v1$counts,weight=GC_bias_v1$counts)) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))
i7_GC_rho_v1 <- as.data.frame(wtd.cor(GC_bias_v1$i7,GC_bias_v1$counts,weight=GC_bias_v1$counts)) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))


GC_plots <- ggarrange(GC_v0_plot,GC_v1_plot, common.legend = TRUE, labels=c("C","D"))


pdf(paste0(barcode_libs,"barcode_v0v1_GC_bias.pdf"), width=7,height=3)
GC_plots
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

read_lengths <- stats %>%
  dplyr::select(V1,V2,V6,V10)
colnames(read_lengths) <- c("library","before_cleaning","after_cleaning","after_resizing")
read_lengths$var <- "read_length"

n_reads <- stats %>%
  dplyr::select(V1,V3,V7,V11)
colnames(n_reads) <- c("library","before_cleaning","after_cleaning","after_resizing")
n_reads$var <- "read_count"

n_bp <- stats %>%
  dplyr::select(V1,V4,V8,V12)
colnames(n_bp) <- c("library","before_cleaning","after_cleaning","after_resizing")
n_bp$var <- "bp_count"


x <- cbind(read_lengths[,1:4],n_reads[,2:4], n_bp[,2:4])

to_paste <- c("", rep("avg read length",3), rep("# reads",3), rep("# bp",3))
colnames(x) <- paste(to_paste, colnames(x), sep = " ")

fwrite(x=x, file=paste0(barcode_libs,paste0("HF-barcode_v1_lib_processing_stats.csv")))

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

# for manuscript: 
mean(med_IS$mean)
min(med_IS$mean)
max(med_IS$mean)
min(med_IS$sd)
max(med_IS$sd)


hist(med_IS$mean)

IS_stats <- med_IS
IS_stats$mean_bins <- cut(IS_stats$mean, breaks = 3)
IS_stats1 <- IS_stats %>%
  group_by(mean_bins) %>%
  tally() %>%
  dplyr::mutate(perc_mean=n/sum(n)*100)
IS_stats1$IS_mean <- paste0(IS_stats1$mean_bins,": n=",IS_stats1$n," ",round(IS_stats1$perc_mean,2),"%")

fwrite(IS_stats1,file=paste0(barcode_libs,"HF-barcode_v1_IS_size_mean.csv"))


IS_stats$median_bins <- cut(IS_stats$median, breaks = 3)
IS_stats2 <- IS_stats %>%
  group_by(median_bins) %>%
  tally() %>%
  dplyr::mutate(perc_median=n/sum(n)*100)
IS_stats2$IS_median <- paste0(IS_stats2$median_bins,": n=",IS_stats2$n," ",round(IS_stats2$perc_median,2),"%")


fwrite(IS_stats2,file=paste0(barcode_libs,"HF-barcode_v1_IS_size_median.csv"))


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

fwrite(GC_stats,file=paste0(barcode_libs,"HF-barcode_v1_GC_bias.csv"))


GC_DF_text$pval <- cut(GC_DF_text$pval, c(0, 0.01, 0.05, 1))

pdf(paste0(barcode_libs,"HF-barcode_v1_GC_plot.pdf"), width=3,height=4)
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


# OligoAnalyzer: 


OligoAnalyzer_output_v0 <- read_excel("Hackflex/middle/OligoAnalyzer_output_v0.xlsx")
OligoAnalyzer_output_v1 <- read_excel("Hackflex/middle/OligoAnalyzer_output_v1.xlsx")

o_v0 <- OligoAnalyzer_output_v0
o_v1 <- OligoAnalyzer_output_v1

# replace name
names(o_v0) <- gsub(x = names(o_v0), pattern = "\\Δ", replacement = "delta")  
names(o_v1) <- gsub(x = names(o_v1), pattern = "\\Δ", replacement = "delta")  

# merge counts details to OligoAnalyzer output
c_v0 <- counts_v0 %>% dplyr::filter(!sample_name == "\"Undetermined\"" & !sample_name == "\"all\"")
o_v0 <- cbind(rbind(c_v0,c_v0),o_v0)
o_v0 <- cSplit(o_v0, "Name", "_")

# merge counts details to OligoAnalyzer output
c_v1 <- counts_v1 %>% dplyr::filter(!sample_name == "\"Undetermined\"" & !sample_name == "\"all\"")
o_v1 <- cbind(rbind(c_v1,c_v1),o_v1)
o_v1 <- cSplit(o_v1, "Name", "_")

fwrite(x=o_v0, file=paste0(barcode_libs,"OligoAnalyzer_output_v0.csv"))
fwrite(x=o_v1, file=paste0(barcode_libs,"OligoAnalyzer_output_v1.csv"))


#######################################

# MANUSCRIPT TEXT: 

# barcodes v0

# The estimated ΔG ranged between
summary(o_v0$deltaG)

# the average melting temperature was
mean(o_v0$Tm)
sd(o_v0$Tm)

# low read count libs: 
z <- o_v0 %>% dplyr::filter(counts<20)
NROW(z)

mean(z$deltaG)
sd(z$deltaG)

z %>% dplyr::select(Tm,Name_2,Name_3)
mean(z$Tm)
sd(z$Tm)



# barcodes v1

# The estimated ΔG ranged between
summary(o_v1$deltaG)

# the average melting temperature was
mean(o_v1$Tm)
sd(o_v1$Tm)

# low read count libs: 
z <- o_v1 %>% dplyr::filter(counts<50000)
NROW(z)

mean(z$deltaG)
sd(z$deltaG)

mean(z$Tm)
sd(z$Tm)

#######################################







# manuscript info - replace 
df <- o_v1 # or o_v1

summary(df$GcContent)
summary(df$Tm)
summary(df$Tm)

sd(df$Tm)
summary(df$deltaG)


# study hairpin formation of barcodes: 
# filter barcodes that produced the lowest read counts
# and with lowest deltaG
z <- df %>% 
  top_n(-6,counts) %>%
  dplyr::select(deltaG,
                Tm,
                counts,
                Name_1,
                Name_2,
                Name_3) %>%
  dplyr::arrange(counts)
z
summary(z$Tm)
summary(z$deltaG)



# DeltaG close to 10 are dangerous
as.data.frame(wtd.cor(df$deltaG,df$counts)) %>%
  dplyr::mutate(label=paste0("rho=",round(correlation,2),
                             "\n p-value=",round(p.value,2)))


df %>%
  ggplot(., aes(x=deltaG,y=counts, color=Name_3))+
  geom_point(size=0.3)+
  geom_smooth(size=0.5)+
  labs(x="deltaG",
       y="counts",
       color="oligo")+
  theme_bw() +
  geom_text(data = dplyr::filter(df, deltaG < -8),aes(label=Name_2), 
            size=3, check_overlap = TRUE, nudge_x = 0.1)

df %>%
  ggplot(., aes(x=Tm,y=counts, color=Name_3))+
  geom_point(size=0.3)+
  geom_smooth(size=0.5)+
  labs(x="Tm",
       y="counts",
       color="oligo")+
  theme_bw()

df %>%
  ggplot(., aes(x=Tm,y=counts, color=Name_3))+
  geom_point(size=0.3)+
  geom_smooth(size=0.5)+
  labs(x="Tm",
       y="counts",
       color="oligo")+
  theme_bw()


