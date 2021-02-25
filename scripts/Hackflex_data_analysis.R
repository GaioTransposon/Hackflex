

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




# set directories
# mydir <- "~/Desktop/MG1655/goal_dilution/ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"

# mydir <- "~/Desktop/MG1655/goal_dilution/paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"

# mydir <- "~/Desktop/MG1655/goal_dilution/saureus/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"

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
# mydir <- "~/Desktop/MG1655/goal_hackflex/saureus/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
# mydir <- "~/Desktop/MG1655/goal_ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
# mydir <- "~/Desktop/MG1655/goal_saureus/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# 
mydir <- "~/Desktop/MG1655/goal_paeruginosa/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"

# insert here libs that you want to plot
my_subset <- c("Ec.SF_1","Ec.SF_1:50","Ec.HF","Ec.HF_55A","Ec.HF_06x","Ec.SF_1_rep","Ec.SF_1_PS","Ec.SF_1:50_rep") # all from E. coli
#my_subset <- c("Sa.SF_1","Sa.SF_1:50","Sa.HF","Sa.HF_55A","Sa.HF_06x") # all from S. aureus
#my_subset <- c("Pa.SF_1","Pa.SF_1:50","Pa.HF","Pa.HF_55A","Pa.HF_55A72E","Pa.HF_06x") # all from P. aeruginosa

########################################
########################################


recode_fun <- function(chars) {
  
  x <- recode(chars, J1 = "Ec.SF_1")
  x <- recode(x, J2 = "Ec.SF_1:50")
  x <- recode(x, K13 = "Ec.HF")
  x <- recode(x, K9 = "Ec.HF_55A")
  x <- recode(x, K14 = "Ec.HF_06x")
  x <- recode(x, K1 = "Ec.SF_1_rep")
  x <- recode(x, K2 = "Ec.SF_1_PS")
  x <- recode(x, K3 = "Ec.SF_1:50_rep")
  
  x <- recode(x, J5 = "Pa.SF_1")
  x <- recode(x, J6 = "Pa.SF_1:50")
  x <- recode(x, K8 = "Pa.HF")
  x <- recode(x, K11 = "Pa.HF_55A")
  x <- recode(x, K5 = "Pa.HF_55A72E")
  x <- recode(x, K15 = "Pa.HF_06x")
  
  x <- recode(x, J9 = "Sa.SF_1")
  x <- recode(x, J10 = "Sa.SF_1:50")
  x <- recode(x, K7 = "Sa.HF")
  x <- recode(x, K10 = "Sa.HF_55A")
  x <- recode(x, K16 = "Sa.HF_06x")
  
  return(x)
  
}

# reorder lib factor - function 
reorder_lib_fun <- function(df) {

  df$library  = factor(df$library, levels=c("Ec.SF_1",
                                         "Ec.SF_1:50",
                                         "Ec.HF",
                                         "Ec.HF_55A",
                                         "Ec.HF_06x",
                                         "Ec.SF_1_rep",
                                         "Ec.SF_1_PS",
                                         "Ec.SF_1:50_rep",
                                         "Pa.SF_1",
                                         "Pa.SF_1:50",
                                         "Pa.HF",
                                         "Pa.HF_55A",
                                         "Pa.HF_55A72E",
                                         "Pa.HF_06x",
                                         "Sa.SF_1",
                                         "Sa.SF_1:50","Sa.HF",
                                         "Sa.HF_55A",
                                         "Sa.HF_06x"))
  # drop unused levels
  df <- df %>% 
    drop.levels()
  
  return(df)
}

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
  head(coverage)
  colnames(coverage) <- c("position","coverage")
  
  id <- sub(".*interleaved2_", "", cov_file)
  id <- sub(".dedup.tsv", "", id)
  id <- gsub("_R1","",id)
  id <- recode_fun(id)
  coverage$library=paste0(as.character(id))
  
  df_to_fill_coverage <- rbind(df_to_fill_coverage,coverage)
  
}


# subset
df_to_fill_coverage <- df_to_fill_coverage %>% dplyr::filter(library %in% my_subset)

# re-order libs
df_to_fill_coverage <- reorder_lib_fun(df_to_fill_coverage)

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
  id <- recode_fun(id)
  new_df$library=paste0(as.character(id))
  
  df_to_fill_fragm_size <- rbind(
    df_to_fill_fragm_size, 
    new_df
  )
  
}

# subset
df_to_fill_fragm_size <- df_to_fill_fragm_size %>% dplyr::filter(library %in% my_subset)

# re-order libs
df_to_fill_fragm_size <- reorder_lib_fun(df_to_fill_fragm_size)

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
  id <- recode_fun(id)
  
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

# subset
big_data <- big_data %>% dplyr::filter(library %in% my_subset)

# re-order libs
big_data <- reorder_lib_fun(big_data)

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


########################################
########################################

# PHRED scores plotting

# open file contaning all the PHRED scores:
PHRED_df <- read_csv(file.path(phred_dir,"PHRED_scores.csv"))

# subset
PHRED_df <- PHRED_df %>% 
  dplyr::mutate(library=as.factor(recode_fun(library))) %>%
  dplyr::filter(library %in% my_subset) %>%
  drop.levels()

# re-order libs
PHRED_df <- reorder_lib_fun(PHRED_df)

PHRED_df$lib <- as.character(paste0(PHRED_df$library,"_",PHRED_df$read))


########################################
########################################


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
  rl_df$library <- recode_fun(rl_df$library)
  rl_df$Sample <- NULL
  
  rl_to_fill <- rbind(rl_to_fill, rl_df)
}

# subset
rl_to_fill <- rl_to_fill %>% dplyr::filter(library %in% my_subset)

# re-order libs
rl_to_fill <- reorder_lib_fun(rl_to_fill)

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
    dplyr::mutate(library=as.factor(recode_fun(Sample))) %>%
    select(library, everything()) # bring at first position
  
  MElist[[row]] <- me_df # add it to your list
  
}


ME_data = do.call(rbind, MElist)

# subset
ME_data <- ME_data %>% dplyr::filter(library %in% my_subset)

# re-order libs
ME_data <- reorder_lib_fun(ME_data)

# store names
n <- ME_data$Sample
# transpose all but the first column (name)
ME_data <- as.data.frame(t(ME_data[,-2]))
colnames(ME_data) <- n
ME_data$Sample <- factor(row.names(ME_data))
rownames(ME_data) <- NULL

fwrite(x=ME_data, file=paste0(mydir,"Alfred_ME_data.csv"))


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
  corrr::rplot(shape = 19, colors = c("red", "green")) 
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
# plot PHRED scores
PHRED_df %>% 
  ggplot(.,
         aes(x=read_position, y=PHRED_means, colour=library, shape = read,
             group=interaction(library, read))) + 
  geom_point(alpha=0.8) + 
  geom_line()+
  xlab("read position (bp)") +
  ylab("average PHRED score") +
  theme_bw() +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300), lim = c(0, 300)) +
  scale_y_continuous(breaks = c(30,32,34,36,38), lim = c(29, 38))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="top",
        legend.title = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18)) 
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


