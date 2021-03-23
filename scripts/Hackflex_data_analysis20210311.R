

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
library(gridExtra)
library(grid)
library(stringr)
library(corrplot)
library(kableExtra)


# args = commandArgs(trailingOnly=TRUE)
# mydir <- args[1]
# phred_dir <- args[2]
# my_subset <- args[3]

# set directories and select samples: 

# mydir <- "~/Desktop/MG1655/goal_ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Ec.SF_1.B1",
#                "Ec.SF_1:50.B1",
#                "Ec.HF.B3",
#                "Ec.HF_55A.B2",
#                "Ec.SF_1.B2",
#                "Ec.SF_1_PS.B2",
#                "Ec.SF_1:50.B2") # all from E. coli

# mydir <- "~/Desktop/MG1655/goal_paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Pa.SF_1.B1",
#                "Pa.SF_1:50.B1",
#                "Pa.HF.B2",
#                "Pa.HF_55A.B2",
#                "Pa.HF_55A72E.B2") # all from P. aeruginosa
# 
# mydir <- "~/Desktop/MG1655/goal_saureus/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Sa.SF_1.B1",
#                "Sa.SF_1:50.B1",
#                "Sa.HF.B2",
#                "Sa.HF_55A.B2") # all from S. aureus

########################################

# mydir <- "~/Desktop/MG1655/goal_size_selection/ecoli/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Ec.HF.B3",
#                "Ec.HF_06x.B3")
# 
# mydir <- "~/Desktop/MG1655/goal_size_selection/paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Pa.HF.B2",
#                "Pa.HF_06x.B3")

mydir <- "~/Desktop/MG1655/goal_size_selection/saureus/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"
my_subset <- c("Sa.HF.B2",
               "Sa.HF_06x.B3")

########################################

# mydir <- "~/Desktop/MG1655/goal_barcode/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- paste0("HF-barcode-", seq(1,96))

########################################
########################################


recode_fun <- function(chars) {
  
  x <- recode(chars, J1 = "Ec.SF_1.B1")
  x <- recode(x, J2 = "Ec.SF_1:50.B1")
  x <- recode(x, K13 = "Ec.HF.B3")
  x <- recode(x, K9 = "Ec.HF_55A.B2")
  x <- recode(x, K14 = "Ec.HF_06x.B3")
  x <- recode(x, K1 = "Ec.SF_1.B2")
  x <- recode(x, K2 = "Ec.SF_1_PS.B2")
  x <- recode(x, K3 = "Ec.SF_1:50.B2")
  
  x <- recode(x, J5 = "Pa.SF_1.B1")
  x <- recode(x, J6 = "Pa.SF_1:50.B1")
  x <- recode(x, K8 = "Pa.HF.B2")
  x <- recode(x, K11 = "Pa.HF_55A.B2")
  x <- recode(x, K5 = "Pa.HF_55A72E.B2")
  x <- recode(x, K15 = "Pa.HF_06x.B3")
  
  x <- recode(x, J9 = "Sa.SF_1.B1")
  x <- recode(x, J10 = "Sa.SF_1:50.B1")
  x <- recode(x, K7 = "Sa.HF.B2")
  x <- recode(x, K10 = "Sa.HF_55A.B2")
  x <- recode(x, K16 = "Sa.HF_06x.B3")
  
  return(x)
  
}

# reorder lib factor - function 
reorder_lib_fun <- function(df) {
  
  df$library  = factor(df$library, levels=c("Ec.SF_1.B1",
                                            "Ec.SF_1:50.B1",
                                            "Ec.HF.B3",
                                            "Ec.HF_55A.B2",
                                            "Ec.HF_06x.B3",
                                            "Ec.SF_1.B2",
                                            "Ec.SF_1_PS.B2",
                                            "Ec.SF_1:50.B2",
                                            "Pa.SF_1.B1",
                                            "Pa.SF_1:50.B1",
                                            "Pa.HF.B2",
                                            "Pa.HF_55A.B2",
                                            "Pa.HF_55A72E.B2",
                                            "Pa.HF_06x.B3",
                                            "Sa.SF_1.B1",
                                            "Sa.SF_1:50.B1",
                                            "Sa.HF.B2",
                                            "Sa.HF_55A.B2",
                                            "Sa.HF_06x.B3"))
  # drop unused levels
  df <- df %>% 
    drop.levels()
  
  return(df)
}

########################################

# grab first two chars, return species (used in ktable)
species <- stringr::str_extract(my_subset[1], "^.{2}")
if (species=="Ec") {  sp <- "Escherichia coli" }
if (species=="Pa") {  sp <- "Pseudomonas aeruginosa" }
if (species=="Sa") {  sp <- "Staphylococcus aureus" }
if (species=="HF") {  sp <- "Hackflex libraries (n=96)" }



########################################
########################################


# Read stats: 

stats <-read.table(file.path(mydir,"reads_stats.tsv"))
head(stats)
stats$V5 <- NULL
stats$V9 <- NULL

stats$V1 <- gsub("interleaved_","",stats$V1)
stats$V1 <- gsub(".fastq","",stats$V1)

stats <- stats %>% 
  dplyr::mutate(library=as.factor(recode_fun(V1))) %>%
  dplyr::filter(library %in% my_subset) %>%
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


########################################
########################################



# Coverage: 


cov_files <- grep(list.files(mydir,pattern="^CO"), 
                  pattern='.tsv', value=TRUE)


# Coverage:

# construct an empty dataframe to build on 
df_to_fill_coverage <- data.frame(
  Sample = character(),
  library = character(),
  Coverage = numeric(),
  Quantile = numeric(),
  stringsAsFactors = FALSE
)



for (cov_file in cov_files) {
  
  print(cov_file)
  
  coverage <-read_table2(file.path(mydir,cov_file))

  # expand rows and select cols
  coverage <- coverage %>% 
    uncount(Count) %>%
    dplyr::select(Sample,Coverage,Quantile)
  
  # clean lib names
  id <- sub(".*interleaved_", "", cov_file)
  id <- sub(".dedup.tsv.tsv", "", id)
  id <- recode_fun(id)
  coverage$library=paste0(as.character(id))
  
  df_to_fill_coverage <- rbind(df_to_fill_coverage,coverage)
  
}

head(coverage)

# subset
df_to_fill_coverage <- df_to_fill_coverage %>% dplyr::filter(library %in% my_subset)


# re-order libs
df_to_fill_coverage$library <- as.factor(df_to_fill_coverage$library)
df_to_fill_coverage <- reorder_lib_fun(df_to_fill_coverage)

# function to get numbers per library: 
get_me_stats <- function(DF) {
  out <- DF %>% 
    dplyr::summarise(min=min(Coverage),
                     median=median(Coverage),
                     mean=mean(Coverage),
                     sd=sd(Coverage),
                     max=max(Coverage),
                     n_sites=n())
  return(out)
}




text <- df_to_fill_coverage %>%
  group_by(library) %>% 
  dplyr::mutate(mean=mean(Coverage),
                sd=sd(Coverage),
                zeros=NROW(which(Coverage==0))) %>%
  dplyr::select(library,mean,sd,zeros) %>%
  distinct() %>%
  dplyr::mutate(label=paste0("mean=",round(mean,2),
                             "\n sd=",round(sd,2),
                             "\n zero=", zeros))

p1 <- ggplot(df_to_fill_coverage, aes(x=Coverage, color=library)) +
  geom_density(alpha=0.01) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(0,100) +
  labs(x="coverage per site",
       y="Frequency") +
  theme_bw()+
  theme(legend.position="none")


p2 <- df_to_fill_coverage %>%
  dplyr::filter(Coverage < 4) %>%
  ggplot(., aes(x=Coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=0.5)+
  facet_grid(rows = vars(library)) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(-0.5,5.5) +
  labs(x="coverage per site",
       y="Frequency") +
  theme_bw()+
  theme(legend.position="none",
        strip.text.y = element_text(size = 7, 
                                    colour = "black", 
                                    angle = 0)) +
  geom_text(
    data    = text,
    mapping = aes(x = Inf, y = Inf, label = label, hjust=1.0, vjust=1), #vjust=1 hjust=1.0
    size=3
  )


# final cov plot
cov_plot <- ggarrange(p1, p2, 
                      ncol=2, nrow=1, common.legend = TRUE,
                      labels=c("A","B"))



########################################


bed_files <- grep(list.files(mydir), 
                  pattern='.bedgraph', value=TRUE)


# Coverage:

# construct an empty dataframe to build on 
df_to_fill_bed <- data.frame(
  contig = character(),
  position = numeric(),
  coverage = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)



for (bed_file in bed_files) {
  
  bed <-read_delim(file.path(mydir,bed_file), 
                   "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  
  # clean lib names
  id <- sub(".*interleaved_", "", bed_file)
  id <- sub(".bedgraph", "", id)
  id <- recode_fun(id)
  bed$library=paste0(as.character(id))
  
  colnames(bed) <- c("contig","position","coverage","library")
  
  df_to_fill_bed <- rbind(df_to_fill_bed,bed)
  
}


# subset
df_to_fill_bed <- df_to_fill_bed %>%
  dplyr::filter(library %in% my_subset) 

# correlation by coverage:
bed_corr <- df_to_fill_bed
head(bed_corr)

bed_corr <- bed_corr %>%
  pivot_wider(names_from=library, values_from=coverage, values_fill = 0)

unique(bed_corr$contig)
tail(bed_corr)

bed_corr_libs <- bed_corr[,3:ncol(bed_corr)]


# compute correlation between libs based on their coverage (all contigs)
res <- cor.mtest(bed_corr_libs, conf.level = .95)
M <- cor(bed_corr_libs, use = "pairwise.complete.obs", method = "pearson")

M



# keep n largest contigs
df_to_fill_bed$contig <- as.character(df_to_fill_bed$contig)
keep <- unique(df_to_fill_bed$contig)[1:4]
# now I will use these IDs to plot interesting stuff
df_to_fill_bed <- subset(df_to_fill_bed, (contig %in% keep))


mycolors<- c('#999999', # grey
             '#E69F00', # orange
             '#56B4E9', # blue
             "#FF3300", # red
             "#333399", # purple
             "#33CC33", # green
             "#FFFF00") # yellow

myshapes<- c(0, 1, 8, 9, 15, 17, 19)

lowest_cov_plot <- df_to_fill_bed %>%
  dplyr::filter(coverage < 10) %>%
  ggplot(., aes(x=position,y=coverage,group=library)) +
  geom_point(aes(shape=library, color=library), size=1)+
  scale_shape_manual(values=myshapes)+
  scale_color_manual(values=mycolors)+
  #facet_grid(rows = vars(contig)) +
  #facet_grid(rows = contig, scales="free_x") +
  facet_grid(library~contig) +
  theme_bw()+
  theme(axis.text.x=element_text(size=6, angle=90),
        strip.text.y = element_text(size = 8, 
                                    colour = "black", 
                                    angle = 0),
        legend.position="right")+
  labs(x="Genomic position (bp)",
       y="Coverage", 
       title = "Lowest coverage regions (4 largest contigs)")

########################################
########################################

# Insert size: 


IS_files <- grep(list.files(mydir,pattern="^picard"), 
                 pattern='.txt', value=TRUE)


# construct an empty dataframe to build on 
df_to_fill_insert_size <- data.frame(
  insert_size = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (IS_file in IS_files) {
  
  # read in file 
  IS_df <- read_delim(file.path(mydir,IS_file), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 10)
  
  IS_df <- IS_df %>%
    uncount(All_Reads.fr_count) 
  
  id <- sub(".*interleaved_", "", IS_file)
  id <- gsub(".dedup.txt","",id)
  id <- recode_fun(id)
  IS_df$library=paste0(as.character(id))
  
  df_to_fill_insert_size <- rbind(
    df_to_fill_insert_size, 
    IS_df
  )
  
}

# subset
df_to_fill_insert_size <- df_to_fill_insert_size %>% dplyr::filter(library %in% my_subset)

# re-order libs
df_to_fill_insert_size$library <- as.factor(df_to_fill_insert_size$library)
df_to_fill_insert_size <- reorder_lib_fun(df_to_fill_insert_size)


# expand rows based on Count, compute median insert sizes (to add to plot)
med_IS <- df_to_fill_insert_size %>%
  group_by(library) %>%
  dplyr::summarize(median=round(median(insert_size),2),
                   mean=round(mean(insert_size),2))

# expand rows based on Count, select cols, and plot
insert_size_plot_facets <- df_to_fill_insert_size %>% 
  dplyr::select(library,insert_size) %>%
  ggplot(., aes(insert_size, colour=library, fill=library)) + 
  scale_x_continuous(limits=c(0,1000)) + 
  facet_grid(rows = vars(library), scales = "free") +
  xlab("insert size (bp)")+
  geom_density(alpha=0.01) +
  geom_vline(data = med_IS, aes(xintercept = median, 
                                color = library), size=0.5)+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=8),
        axis.title=element_text(size=14), 
        strip.text.y = element_text(size = 8, 
                                    colour = "black", 
                                    angle = 0)) +
  geom_text(data=med_IS, 
            aes(x=900, y=Inf, label=paste0("median=",median,"\n","mean=",mean), color=library), size=2, hjust=1, vjust=1.5)

insert_size_plot_together <- df_to_fill_insert_size %>% 
  dplyr::select(library,insert_size) %>%
  ggplot(., aes(insert_size, colour=library, fill=library)) + 
  scale_x_continuous(limits=c(0,1000)) + 
  xlab("insert size (bp)")+
  geom_density(alpha=0.01) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=8),
        axis.title=element_text(size=14))

# final plot
insert_size_plot <- ggarrange(insert_size_plot_together, insert_size_plot_facets, 
                              ncol=2, nrow=1, common.legend = TRUE,
                              labels=c("A","B"))

########################################
########################################


# GC content: 

gc_files = list.files(mydir,pattern="GC_qc")

GC_DF <- data.frame(
  REF_GCcontent = numeric(),
  LIB_Sample = character(),
  diff = numeric(),
  rho = numeric(),
  stringsAsFactors = FALSE
)

for (gc_file in gc_files) {
  
  gc <- gc_file
  
  id <- sub(".*interleaved_", "",gc)
  id <- gsub(".dedup.tsv.tsv","",id)
  id <- recode_fun(id)
  
  # read in files
  gc_df <- read.table(file.path(mydir,gc_file), quote="\"", comment.char="", header = TRUE)
  head(gc_df)
  
  # add smallest number above zero (otherwise 0 divided by 0 will give a problem)
  my_min <- min(gc_df[gc_df$fractionOfReads>0,"fractionOfReads"])
  # add it to fraction of reads
  gc_df$fractionOfReads <- (gc_df$fractionOfReads + my_min)
  
  # separate reference from library stats, prepend prefix
  ref <- gc_df %>% dplyr::filter(Sample=="Reference") %>% rename_all( ~ paste0("REF_", .x))
  lib <- gc_df %>% dplyr::filter(!Sample=="Reference") %>% rename_all( ~ paste0("LIB_", .x))
  
  # rename the library 
  lib$LIB_Sample <- paste0(id)
  
  # reunite them
  myDF <- cbind(ref,lib)
  
  # compute differnece between fraction of reads expected (reference) and fraction of reads obtained (library)
  myDF$diff <- myDF$REF_fractionOfReads/myDF$LIB_fractionOfReads
  
  # clean the df, rename a col: 
  myDF <- myDF %>% dplyr::select(REF_GCcontent,LIB_Sample,diff)
  colnames(myDF)[colnames(myDF) == 'LIB_Sample'] <- 'library'
  head(myDF)
  
  # get correlation between GC content and ratio, using fraction of reads as weight
  rho <- wtd.cor(myDF$REF_GCcontent,myDF$diff,weight=myDF$LIB_fractionOfReads)
  rho <- rho[1,1]
  myDF$rho <- rho
  
  GC_DF <- rbind(GC_DF, myDF)
  
}


# subset
GC_DF <- GC_DF %>% dplyr::filter(library %in% my_subset)

# re-order libs
GC_DF$library <- as.factor(GC_DF$library)
GC_DF <- reorder_lib_fun(GC_DF)

head(GC_DF)

GC_DF_text <- GC_DF %>% dplyr::select(library,rho) %>% distinct() %>%
  dplyr::arrange(library) %>% # to get the factors order
  dplyr::mutate(label=paste0(library), 
                label_rho=paste0("rho=",round(rho,3)),  
                pos=7.5+seq(1:NROW(.)))

# smooth_GC <- GC_DF %>% 
#   dplyr::filter(diff!=1) %>% # don't show in plot ratio=1; rhos is already calculated 
#   ggplot(.,aes(x=REF_GCcontent,y=diff,color=library))+
#   geom_point(alpha=0.3)+
#   theme_bw() +
#   xlim(0.1,0.9) +
#   ylim(-2.5,15)+
#   geom_smooth(size=0.5, alpha=0) +
#   ylab("ratio observed/expected reads") +
#   theme(legend.position="none") +
#   geom_text(
#     data    = GC_DF_text,
#     mapping = aes(x = 0.50, y = pos, label = label, hjust=0, vjust=1), #vjust=1 hjust=1.0
#     size=3
#   )

straight_GC <- GC_DF %>% 
  #dplyr::filter(diff!=1) %>% # don't show in plot ratio=1; rhos is already calculated 
  ggplot(.,aes(x=REF_GCcontent,y=diff,color=library))+
  geom_point(alpha=0.3)+
  theme_bw() +
  xlim(0.1,0.9) +
  ylim(-2.5,15)+
  stat_smooth(method="lm", se=FALSE, size=0.5) +
  ylab("ratio observed/expected reads") +
  theme(legend.position="none") +
  geom_text(
    data    = GC_DF_text,
    mapping = aes(x = 0.45, y = pos, label = label, hjust=0, vjust=1), #vjust=1 hjust=1.0
    size=3
  ) +
  geom_text(
    data    = GC_DF_text,
    mapping = aes(x = 0.75, y = pos, label = label_rho, hjust=0, vjust=1), #vjust=1 hjust=1.0
    size=3
  )

########################################
########################################

# PHRED scores 

PHRED_files = list.files(phred_dir,pattern=".csv")

datalist = list()
for (PHRED_file in PHRED_files) {
  
  # open file contaning all the PHRED scores:
  PHRED_df <- read_csv(file.path(phred_dir,PHRED_file))
  head(PHRED_df)
  
  PHRED_df$library <- gsub("_001","",PHRED_df$library) # to remove (from barcode libs)
  
  # get last two chars = read direction
  PHRED_df$read <- str_sub(PHRED_df$library, -2, -1)
  
  # remove everything from library string after first underscore: 
  PHRED_df$library <- gsub("\\_.*","",PHRED_df$library)
  
  # re-order libs (if this is NOT a HF-barcode libs)
  if (str_sub(unique(PHRED_df$library), 1,10) != "HF-barcode") {
    
    # subset
    PHRED_df <- PHRED_df %>% 
      dplyr::mutate(library=as.factor(recode_fun(library))) %>%
      dplyr::filter(library %in% my_subset) %>%
      drop.levels()
    
    PHRED_df$library <- as.factor(PHRED_df$library)
    PHRED_df <- reorder_lib_fun(PHRED_df)
    
  }
  
  else {
    
    # subset
    PHRED_df <- PHRED_df %>% 
      dplyr::filter(library %in% my_subset) %>%
      drop.levels()
    
  }
  
  datalist[[PHRED_file]] <- PHRED_df # add it to your list
  
}


phred_data = do.call(rbind, datalist)

PHRED_plot <- phred_data %>% 
  ggplot(.,
         aes(x=read_position, y=PHRED_means, colour=library, shape = read,
             group=interaction(library, read))) + 
  geom_point(alpha=0.8, size=0.8) + 
  geom_line(size=0.1)+
  xlab("read position (bp)") +
  ylab("average PHRED score") +
  theme_bw() +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300), lim = c(0, 300)) +
  scale_y_continuous(breaks = c(30,32,34,36,38), lim = c(29, 38))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="right",
        legend.title = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18)) 
if (NROW(unique(phred_data$library))>10) {
  PHRED_plot <- PHRED_plot + 
    theme(legend.position="none") +
    ggtitle("HF barcode libraries")
}
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
  Read = character(),
  stringsAsFactors = FALSE
)

for (rl_file in rl_files) {
  
  # read in file
  rl_df <- read.table(file.path(mydir,rl_file), quote="\"", comment.char="", header = TRUE)
  
  unique(rl_df$Library)
  head(rl_df)
  
  rl_df <- rl_df %>% 
    dplyr::select(Sample,Readlength,Count, Fraction, Read)
  
  rl_df$library <- str_replace_all(rl_df$Sample, "reduced_trimmed2_trimmed_interleaved_", "")
  rl_df$library <- gsub(".dedup","",rl_df$library)
  rl_df$library <- recode_fun(rl_df$library)
  rl_df$Sample <- NULL
  
  rl_to_fill <- rbind(rl_to_fill, rl_df)
}

# subset
rl_to_fill <- rl_to_fill %>% dplyr::filter(library %in% my_subset)

# re-order libs
rl_to_fill$library <- as.factor(rl_to_fill$library)
rl_to_fill <- reorder_lib_fun(rl_to_fill)

# expand rows based on count (to create a density plot, like the insert size one)
#rl_to_fill_expand <- setDT(expandRows(rl_to_fill, "Count"))
#head(rl_to_fill_expand)

# plot, keep both reads
RL_plot <- ggplot(rl_to_fill, aes(x=Readlength, y=Count,colour=Read, fill=Read)) + 
  scale_x_continuous(limits=c(100,310)) + 
  xlab("length of mapped reads (bp)")+
  geom_bar(alpha=0.5, stat="identity", width = 0.001)+
  facet_grid(rows = vars(library), scales = "free") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=8),
        axis.title=element_text(size=14),
        strip.text.y = element_text(size = 8, 
                                    colour = "black", 
                                    angle = 0))


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
  
  
  me_df$Sample <- sub(".*interleaved_", "",me_df$Sample)
  me_df$Sample <- gsub(".dedup", "",me_df$Sample)
  
  me_df <- me_df %>%
    dplyr::mutate(library=as.factor(recode_fun(Sample))) %>%
    select(library, everything()) # bring at first position
  
  MElist[[row]] <- me_df # add it to your list
  
}


ME_data = do.call(rbind, MElist)

# subset
ME_data <- ME_data %>% dplyr::filter(library %in% my_subset)

# re-order libs
ME_data$library <- as.factor(ME_data$library)
ME_data <- reorder_lib_fun(ME_data)

# store names
n <- ME_data$Sample

# transpose all but the first column (name)
ME_data <- as.data.frame(t(ME_data[,-2]))
colnames(ME_data) <- n
ME_data$Sample <- factor(row.names(ME_data))
rownames(ME_data) <- NULL

# Save as Table : 
fwrite(x=ME_data, file=paste0(mydir,"Alfred_ME_data.csv"))

########################################
########################################

# Contamination: 


kraken_files = list.files(mydir,pattern="^kraken_")


# construct an empty dataframe to build on
df_to_fill_kraken <- data.frame(
  library = character(),
  read_count = numeric(),
  percentage_of_tot_unmapped_reads = numeric(),
  species = character(),
  stringsAsFactors = FALSE
)

for (kraken_file in kraken_files) {
  
  kra <- read_delim(file.path(mydir,kraken_file),"\t", 
                    escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE)
  
  id <- kraken_file
  id <- sub(".*interleaved_", "", id)
  id <- gsub(".dedup","",id)
  id <- recode_fun(id)
  kra$library=paste0(as.character(id))
  
  
  # The output of kraken-report is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:
  #   1. Percentage of reads covered by the clade rooted at this taxon
  # 2. Number of reads covered by the clade rooted at this taxon
  # 3. Number of reads assigned directly to this taxon
  # 4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
  # 5. NCBI taxonomy ID
  # 6. indented scientific name
  
  kra_filtered <- kra %>%
    dplyr::filter(X3>0) %>% 
    dplyr::mutate(sum_reads=sum(X3)) %>%
    dplyr::mutate(percentage_of_tot_unmapped_reads=(X3/sum_reads)*100) %>%
    dplyr::select(library,X3,percentage_of_tot_unmapped_reads,X6)
  
  
  colnames(kra_filtered) <- c("library","read_count","percentage_of_tot_unmapped_reads", "species")
  
  df_to_fill_kraken <- rbind(
    df_to_fill_kraken, 
    kra_filtered
  )
  
}

# subset
df_to_fill_kraken <- df_to_fill_kraken %>% dplyr::filter(library %in% my_subset)

# re-order libs
df_to_fill_kraken$library <- as.factor(df_to_fill_kraken$library)
df_to_fill_kraken <- reorder_lib_fun(df_to_fill_kraken)

# Save as Table : 
fwrite(x=df_to_fill_kraken, file=paste0(mydir,"kraken_contamination.csv"))



################################################################################
################################################################################
################################################################################



# plot
pdf(paste0(mydir,species,'_out.pdf'))
# plot coverage 
cov_plot
# print lowest coverage regions
lowest_cov_plot
# plot correlation between libs based on their coverage
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corr numbers (rho)
corrplot(M, method = "color", col = col(200),
         type = "upper", order = "AOE", 
         tl.col = "black", tl.srt = 90, tl.cex = .8, # Text label color, size, and rotation
         # hide correlation coefficient on the principal diagonal
         diag = FALSE, title = "rho based on coverage \n (from bedgraphs - all contigs)",
         # Combine with significance
         p.mat = res$p, sig.level = c(.001, .01, .05), insig = "blank",
         number.cex = .6, addCoef.col = "black") # Add coefficient of correlation
# significance symbols
corrplot(M, method = "color", col = col(200),
         type = "upper", order = "AOE", 
         tl.col = "black", tl.srt = 90, tl.cex = .8, # Text label color, size, and rotation
         # hide correlation coefficient on the principal diagonal
         diag = FALSE, title = "pvalues of rho based on coverage \n (from bedgraphs - all contigs)",
         # Combine with significance
         p.mat = res$p, sig.level = c(.001, .01, .05), insig = "label_sig", 
         pch.cex = .8, pch.col = "white")
#correlation plot
bed_libs %>% correlate() %>%  
  network_plot(min_cor = 0.01)
# plot insert size 
insert_size_plot 
# plot GC content bias
# ggarrange(smooth_GC, straight_GC, 
#           nrow = 2)
straight_GC
# plot PHRED scores
PHRED_plot
# plot read lengths
RL_plot
dev.off()

################################################################################
################################################################################

# ME data save as pdf : 
kbl(ME_data, format = "html") %>%
  kable_classic() %>%
  kable_styling(font_size = 15, "striped") %>%
  kableExtra::as_image(file=paste0(mydir,species,"_ME_data.png")) #%>%
#save_kable(file = paste0(mydir,species,"_ME_data.pdf"))

# stats libs before, after cleaning and after resizing, print pdf: 
x <- cbind(read_lengths[,1:4],n_reads[,2:4], n_bp[,2:4])
kbl(x, format = "html") %>%
  kable_classic() %>%
  kable_styling(font_size = 15, "striped") %>%
  add_header_above(c(" " = 1, 
                     "avg read length" = 3, 
                     "# reads" = 3,
                     "# bp" = 3)) %>%
  kableExtra::as_image(file=paste0(mydir,species,"_stats.png")) #%>%
#save_kable(file = paste0(mydir,species,"_stats.pdf"))


# concatenating created pdfs: 

# setwd(mydir)
# ## Setwd and Collect the names of the figures to be glued together
# ff <- dir(pattern=".pdf")
# ## The name of the pdf doc that will contain all the figures
# outFileName <- paste0(species,"_all.pdf")
# 
# ## Make a system call to pdftk
# system2(command = "pdftk",
#         args = c(shQuote(ff), "cat output", shQuote(outFileName)))

################################################################################
################################################################################
