

# load libs
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)
library(weights)
library(readr)
library(stringr)




# set directories
# mydir <- "~/Desktop/MG1655/goal_dilution/ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"

# mydir <- "~/Desktop/MG1655/goal_dilution/paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"

mydir <- "~/Desktop/MG1655/goal_dilution/saureus/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"

# mydir <- "~/Desktop/MG1655/goal_polymerase/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
# mydir <- "~/Desktop/MG1655/goal_KAPA/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
# mydir <- "~/Desktop/MG1655/goal_hackflex/ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
# mydir <- "~/Desktop/MG1655/goal_hackflex/paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
mydir <- "~/Desktop/MG1655/goal_hackflex/saureus/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"

########################################
########################################


# Coverage: 


cov_files <- grep(list.files(mydir,pattern="^reduced_trimmed2_trimmed_interleaved2_"), 
                  pattern='.tsv', value=TRUE)

# Coverage:

# construct an empty dataframe to build on 
df_to_fill <- data.frame(
  position = numeric(),
  coverage = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (cov_file in cov_files) {
  
  print(cov_file)
  
  coverage <-read.table(file.path(mydir,cov_file))
  
  colnames(coverage) <- c("position","coverage")
  
  id <- str_replace_all(cov_file, "reduced_trimmed2_trimmed_interleaved2_", "")
  coverage$library=paste0(as.character(id))
  
  df_to_fill <- rbind(df_to_fill,coverage)
  
}

p1 <- ggplot(df_to_fill, aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=1)+
  facet_grid(rows = vars(library)) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(0,100) +
  labs(x="coverage per site",
       y="Frequency")

p2 <- df_to_fill %>% 
  dplyr::filter(coverage <= 3) %>%
  ggplot(., aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=1)+
  facet_grid(rows = vars(library)) +
  theme(legend.position="none")+
  xlim(0,4) +
  labs(x="coverage per site",
       y="Frequency")

p3 <- df_to_fill %>% 
  dplyr::filter(coverage <= 3) %>%
  ggplot(., aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=1)+
  facet_grid(rows = vars(library), scales = "free_y") +
  theme(legend.position="none")+
  xlim(0,4) +
  labs(x="coverage per site",
       y="Frequency")




# plot
pdf(paste0(mydir,'coverage.pdf'))
ggarrange(p1,                                                 # First row with first plot
          ggarrange(p2, p3, ncol = 2, labels = c("B", "C")), # Second row with 2 plots
          nrow = 2, 
          labels = "A"                                        # Labels of the first plot
)
dev.off()


########################################
########################################


# Fragment size: 

frag_files = list.files(mydir,pattern="frag_sizes")

# construct an empty dataframe to build on 
df_to_fill <- data.frame(
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
  
  new_df$library <- as.character(basename(frag_file))
  
  df_to_fill <- rbind(
    df_to_fill, 
    new_df
  )
  
}

# reduce lib name
df_to_fill$library <- str_replace_all(df_to_fill$library, "frag_sizes_reduced_trimmed2_trimmed_interleaved2_", "")

# plot
pdf(paste0(mydir,'fragment_size.pdf'))
#ds <- subset(df_to_fill, ! is.na(value))
ds <- df_to_fill
ggplot(ds, aes(V1, colour=library, fill=library)) + 
  scale_x_continuous(limits=c(0,1000)) + 
  xlab("fragment size (bp)")+
  geom_density(alpha=0.55) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
dev.off()


########################################
########################################


# GC content: 

gc_files = list.files(mydir,pattern="GC_qc")
flagstat_files = list.files(mydir,pattern="flagstat")

myDF <- data.frame(gc_files,flagstat_files, stringsAsFactors = F)

datalist = list()
for (row in 1:nrow(myDF)) {
  
  gc <- myDF[row,1]
  flag <- myDF[row,2]
  
  id <- str_replace_all(gc, "GC_qc_reduced_trimmed2_trimmed_interleaved2_", "")
  
  # read in files
  gc_df <- read.table(file.path(mydir,gc), quote="\"", comment.char="", header = TRUE)
  flag_df <- read.delim(file.path(mydir,flag), header=FALSE)
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


smooth_GC <- big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=lib))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  geom_smooth() +
  ylab("ratio observed/expected reads") +
  theme(legend.position="none") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)

straight_GC <- big_data %>% 
  dplyr::filter(diff!=1) %>% 
  ggplot(.,aes(x=GCcontent,y=diff,color=lib))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  stat_smooth(method="lm", se=FALSE) +
  ylab("ratio observed/expected reads") +
  theme(legend.title = element_blank(),
        legend.position="top") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)


# plot
pdf(paste0(mydir,'GC_content.pdf'))
ggarrange(smooth_GC, straight_GC, 
          nrow = 2)
dev.off()


########################################
########################################

# PHRED scores plotting

# open file contaning all the PHRED scores:
PHRED_df <- read_csv(file.path(phred_dir,"PHRED_scores.csv"))
head(PHRED_df)

# subset to contain only libraries that are of interest in this script
libs_here <- substr(tools::file_path_sans_ext(unique(big_data$lib)), start = 1, stop = 3)

PHRED_df_sub <- PHRED_df %>% 
  dplyr::filter(str_detect(library, str_c(libs_here, collapse="|")))
  
  
# #plot
pdf(paste0(mydir,'PHRED_scores_raw_reads.pdf'))
ggplot(PHRED_df_sub,
       aes(x=read_position, y=PHRED_means, group=library, colour=library ) ) +
  geom_line(size=1, alpha=0.5) +
  xlab("read position (bp)") +
  ylab("average PHRED score") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="top",
        legend.title = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
scale_x_continuous(breaks = c(0,50,100,150,200,250,300), lim = c(0, 300)) +
scale_y_continuous(breaks = c(30,32,34,36,38), lim = c(29, 38))
dev.off()



