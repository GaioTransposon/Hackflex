

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
library(Hmisc)


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
#
# mydir <- "~/Desktop/MG1655/goal_paeruginosa/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Pa.SF_1.B1",
#                "Pa.SF_1:50.B1",
#                "Pa.HF.B2",
#                "Pa.HF_55A.B2",
#                "Pa.HF_55A72E.B2") # all from P. aeruginosa
# # 
mydir <- "~/Desktop/MG1655/goal_saureus/"
phred_dir <- "~/Desktop/MG1655/raw_libs/"
my_subset <- c("Sa.SF_1.B1",
               "Sa.SF_1:50.B1",
               "Sa.HF.B2",
               "Sa.HF_55A.B2") # all from S. aureus

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

# mydir <- "~/Desktop/MG1655/goal_size_selection/saureus/"
# phred_dir <- "~/Desktop/MG1655/raw_libs/"
# my_subset <- c("Sa.HF.B2",
#                "Sa.HF_06x.B3")

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

# add prefix to output files, if the goal of the analysis is to compare HF size selection libs with HF: 
if (grepl("goal_size_selection", mydir, fixed = FALSE)==TRUE) {
  goal <- "SizeSel"
} else {
  goal <- "all"
}

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
  scale_y_continuous(breaks = c(26,28,30,32,34,36,38), lim = c(25, 38))+
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


head(phred_data)

phred_data_stats <- phred_data %>% 
  group_by(library) %>%
  dplyr::summarise(
    bp_count = n(),
    PHRED_mean = mean(PHRED_means, na.rm = TRUE),
    PHRED_sd = sd(PHRED_means, na.rm = TRUE),
    PHRED_variance = var(PHRED_means, na.rm = TRUE)
  ) %>%
  arrange(desc(PHRED_mean))

fwrite(x=phred_data_stats, file=paste0(mydir,paste0(species,"_",goal,"_phred.csv")))

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


x <- cbind(read_lengths[,1:4],n_reads[,2:4], n_bp[,2:4])

# # stats libs before, after cleaning and after resizing, print png: 
# kbl(x, format = "html") %>%
#   kable_classic() %>%
#   kable_styling(font_size = 15, "striped") %>%
#   add_header_above(c(" " = 1, 
#                      "avg read length" = 3, 
#                      "# reads" = 3,
#                      "# bp" = 3)) %>%
#   kableExtra::as_image(file=paste0(mydir,species,"_stats.png")) 

to_paste <- c("", rep("avg read length",3), rep("# reads",3), rep("# bp",3))
colnames(x) <- paste(to_paste, colnames(x), sep = " ")

fwrite(x=x, file=paste0(mydir,paste0(species,"_",goal,"_lib_processing_stats.csv")))

########################################
########################################


# BBDUK cleaning get stats 

bbduk_files = list.files(mydir,pattern="cleaning2")

datalist = list()
for (bbduk_file in bbduk_files) {
  
  # open file contaning all the PHRED scores:
  bb <- read_delim(file.path(mydir,bbduk_file), 
                   ":", escape_double = FALSE, col_names = FALSE, 
                   trim_ws = TRUE, skip = 18)
  
  head(bb)
  
  bb <- bb[2:8,]
  # split to get the node number 
  bb <- cSplit(bb, "X2", "\t", drop = TRUE)
  
  bb <- cSplit(bb, "X2_1", "reads", drop = TRUE, stripWhite = FALSE)
  bb <- cSplit(bb, "X2_2", "bases", drop = TRUE, stripWhite = FALSE)
  
  bb$X2_1_2 <- gsub("[ ()]", "", bb$X2_1_2)
  bb$X2_1_2 <- gsub("[%]", "", bb$X2_1_2)
  bb$X2_2_2 <- gsub("[ ()]", "", bb$X2_2_2)
  bb$X2_2_2 <- gsub("[%]", "", bb$X2_2_2)
  
  colnames(bb) <- c("X1","reads","perc_reads","bases","perc_bases")
  
  
  # clean lib names
  id <- sub(".*interleaved_", "", bbduk_file)
  id <- sub(".bedgraph", "", id)
  id <- recode_fun(id)
  bb$library=paste0(as.character(id))
  
  # subset
  bb <- bb %>%
    dplyr::filter(library %in% my_subset) %>%
    drop.levels()
  
  datalist[[bbduk_file]] <- bb # add it to your list
  
}

bbduk_data = do.call(rbind, datalist)

fwrite(x=bbduk_data, file=paste0(mydir,paste0(species,"_",goal,"_bbduk.csv")))


########################################
########################################


# open dups file:


dups_files = list.files(mydir,pattern="dups_stats.txt")

dups_stats <- read_delim(file.path(mydir,dups_files), 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    trim_ws = TRUE)

# construct an empty dataframe to build on 
dups_df <- data.frame(
  PCR_duplicates = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)


libs <- dups_stats %>% 
  dplyr::filter(str_detect(X1, "^red")) %>%
  dplyr::mutate(library=sub(".*interleaved_", "", X1)) %>%
  dplyr::mutate(library=recode_fun(library)) %>%
  dplyr::select(library)
  
PCR_duplicates <- dups_stats %>% 
  dplyr::filter(str_detect(X1, "^DUPLICATE TOTAL")) %>%
  dplyr::mutate(PCR_duplicates=as.numeric(sub("DUPLICATE TOTAL ", "", X1))) %>%
  dplyr::select(PCR_duplicates)



dups_df <- cbind(libs,PCR_duplicates) %>%
  dplyr::filter(library %in% my_subset) %>%
  drop.levels()

fwrite(x=dups_df, file=paste0(mydir,paste0(species,"_",goal,"_dups.csv")))


########################################
########################################


# Insert size: 


IS_files <- grep(list.files(mydir,pattern="^picardIS"), 
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


fwrite(med_IS,file = paste0(mydir,species,"_",goal,"_insert_size.csv"))

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


# # GC content: 
# 
# gc_files = list.files(mydir,pattern="GC_qc")
# 
# GC_DF <- data.frame(
#   library = character(),
#   Lib_GCcontent = numeric(),
#   Lib_fractionOfReads = numeric(),
#   obs.exp = numeric(),
#   rho = numeric(),
#   pval = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (gc_file in gc_files) {
#   
#   gc <- gc_file
#   
#   id <- sub(".*interleaved_", "",gc)
#   id <- gsub(".dedup.tsv.tsv","",id)
#   id <- recode_fun(id)
#   
#   # read in files
#   gc_df <- read.table(file.path(mydir,gc_file), quote="\"", comment.char="", header = TRUE)
#   head(gc_df)
# 
#   ref <- gc_df[1:102,] %>%
#     dplyr::select(Sample,GCcontent,fractionOfReads) 
#   colnames(ref) <- paste("Ref", colnames(ref), sep = "_")
#   ref <- ref %>%
#     dplyr::mutate(bin=seq(1,102))
#   head(ref)
#   
#   lib <- gc_df[103:204,] %>%
#     dplyr::select(Sample,GCcontent,fractionOfReads) 
#   colnames(lib) <- paste("Lib", colnames(lib), sep = "_")
#   lib <- lib %>%
#     dplyr::mutate(bin=seq(1,102))
#   head(lib)
#   
#   gc_df <- inner_join(ref,lib) %>%
#     dplyr::mutate(obs.exp=Lib_fractionOfReads/Ref_fractionOfReads) %>%
#     filter_all(all_vars(!is.na(.))) # remove NaN (These come from 0/0 where no reads are expected in Reference, and neither are observed)
# 
#   # rename the library
#   gc_df <- gc_df %>%
#     dplyr::mutate(Lib_Sample=paste0(id))
#   head(gc_df)
#   
#   # get correlation between GC content and ratio, using fraction of reads as weight
#   rho <- wtd.cor(gc_df$Lib_GCcontent,
#                  gc_df$obs.exp,
#                  weight=gc_df$LIB_fractionOfReads)
#   pval <- rho[1,4]
#   rho <- rho[1,1]
#   
#   myDF <- gc_df %>%
#     dplyr::mutate(rho=rho,
#                   pval=pval,
#                   library=Lib_Sample) %>%
#     dplyr::select(library,Lib_GCcontent,Lib_fractionOfReads,obs.exp,rho,pval)
#   
#   head(myDF)
#   head(GC_DF)
#   GC_DF <- rbind(GC_DF, myDF)
#   
# }
# 
# 
# # subset
# GC_DF <- GC_DF %>% dplyr::filter(library %in% my_subset)
# 
# # re-order libs
# GC_DF$library <- as.factor(GC_DF$library)
# GC_DF <- reorder_lib_fun(GC_DF)
# 
# head(GC_DF)
# 
# mymax <- max(GC_DF$obs.exp)
# GC_DF_text <- GC_DF %>% dplyr::select(library,rho,pval) %>% 
#   distinct() %>%
#   dplyr::arrange(library) %>% # to get the factors order
#   dplyr::mutate(label=paste0(library), 
#                 label_rho=paste0("rho=",
#                                  round(rho,3)),
#                 label_pval=paste0("p=",
#                                  round(pval,5)),  
#                 pos=max(mymax)-seq(1:NROW(.)))
# 
# head(GC_DF)
# 
# 
# straight_GC <- GC_DF %>% 
#   #dplyr::filter(diff!=1) %>% # don't show in plot ratio=1; rhos is already calculated 
#   ggplot(.,aes(x=Lib_GCcontent,y=obs.exp,color=library))+
#   geom_point(alpha=0.3)+
#   theme_bw() +
#   xlim(0,0.9) +
#   #ylim(-2.5,15)+
#   stat_smooth(method="lm", se=FALSE, size=0.5) +
#   ylab("ratio observed/expected reads") +
#   theme(legend.position="none") +
#   geom_text(
#     data    = GC_DF_text,
#     mapping = aes(x = 0.45, y = pos, label = label, hjust=0, vjust=1), #vjust=1 hjust=1.0
#     size=3
#   ) +
#   geom_text(
#     data    = GC_DF_text,
#     mapping = aes(x = 0.65, y = pos, label = label_rho, hjust=0, vjust=1), #vjust=1 hjust=1.0
#     size=3
#   ) +
#   geom_text(
#     data    = GC_DF_text,
#     mapping = aes(x = 0.75, y = pos, label = label_pval, hjust=0, vjust=1), #vjust=1 hjust=1.0
#     size=3
#   )

########################################
########################################

# GC content - Picard


gc_files <- grep(list.files(mydir,pattern="^picardGC_red"), 
                 pattern='.txt', value=TRUE)


GC_DF <- data.frame(
  GC = numeric(),
  WINDOWS = numeric(),
  MEAN_BASE_QUALITY = numeric(),
  NORMALIZED_COVERAGE = numeric(),
  library = character(),
  rho = numeric(),
  pval = numeric(),
  n_obs = numeric(),
  stringsAsFactors = FALSE
)

for (gc_file in gc_files) {
  
  gc <- gc_file
  
  id <- sub(".*interleaved_", "",gc)
  id <- gsub(".dedup.txt","",id)
  id <- recode_fun(id)
  
  # read in files
  gc_df <- read_delim(file.path(mydir,gc_file),
                      "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)
  
  # rename the library and select cols
  gc_df <- gc_df %>%
    dplyr::mutate(library=paste0(id)) %>%
    dplyr::select(GC,WINDOWS,MEAN_BASE_QUALITY,NORMALIZED_COVERAGE,library) %>%
    dplyr::filter(NORMALIZED_COVERAGE>0)
  
  head(gc_df)
  NROW(gc_df)
  
  # get correlation between GC content and ratio, using fraction of reads as weight
  rho <- wtd.cor(gc_df$GC,
                 gc_df$NORMALIZED_COVERAGE,
                 weight=gc_df$WINDOWS)
  pval <- rho[1,4]
  rho <- rho[1,1]
  
  myDF <- gc_df %>%
    dplyr::mutate(rho=rho,
                  pval=pval,
                  n_obs=NROW(gc_df)) 
  
  head(myDF)
  
  GC_DF <- rbind(GC_DF, myDF)
  
}


head(GC_DF)

pic2 <- GC_DF


mymax <- max(pic2$NORMALIZED_COVERAGE)

pic2_text <- pic2 %>% dplyr::select(library,rho,pval,n_obs) %>%
  distinct() %>%
  dplyr::arrange(library) %>% # to get the factors order
  dplyr::mutate(label=paste0(library),
                label_rho=paste0("R=",
                                 round(rho,3)),
                label_pval=paste0("p=",
                                  round(pval,3)))
                #pos=max(mymax)-seq(1:NROW(.)))

fwrite(pic2_text,file = paste0(mydir,species,"_",goal,"_GC_bias.csv"))

#scaleFactor <- (median(pic2$NORMALIZED_COVERAGE) / median(pic2$MEAN_BASE_QUALITY))*1.5
GC_picard <- pic2 %>%
  ggplot(., aes(x = GC, y = NORMALIZED_COVERAGE)) + 
  geom_point(mapping = aes(colour=WINDOWS), alpha=0.5,size=2, shape=1) + #shape=1, stat = "identity"
  scale_color_gradientn(colours = rev(rainbow(5))) +
  stat_smooth(size=0.3, method = "lm", colour="black", aes(weight= WINDOWS)) +
  #geom_line(mapping = aes(x = GC, y = MEAN_BASE_QUALITY*scaleFactor), size = 0.3, color = "green") +
  # scale_y_continuous(name = "Fraction of normalized coverage", limits = c(0,2),
  #                    sec.axis = sec_axis(~./scaleFactor, name = "Mean base quality")) + 
  ylim(-2,2)+
  xlim(0,90)+
  theme(
    #axis.title.y.right = element_text(color = "green"),
    axis.title.y = element_text(color = "black")) +
  facet_grid(rows = vars(library)) +
  geom_text(
    data    = pic2_text,
    mapping = aes(x = 60, y = Inf, label="R=", fontface=3, hjust=0, vjust=1), #vjust=1 hjust=1.0
    size=3
  ) +
  geom_text(
    data    = pic2_text,
    mapping = aes(x = 64, y = Inf, label=round(rho,3), hjust=0, vjust=1), #vjust=1 hjust=1.0
    size=3
  ) + 
  geom_text(
    data    = pic2_text,
    mapping = aes(x = 75, y = Inf, label = label_pval, hjust=0, vjust=1), #vjust=1 hjust=1.0
    size=3
  ) 

########################################
########################################




# Coverage:

bed_files <- grep(list.files(mydir), 
                  pattern='.bedgraph', value=TRUE)


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
# re-order libs
df_to_fill_bed$library <- as.factor(df_to_fill_bed$library)
df_to_fill_bed <- reorder_lib_fun(df_to_fill_bed)


# keep contigs larger than 1Mbp ; needs to be done for P. aeruginosa only as the assembly is made out of 46 contigs
if (species=="Pa") {  
  
  # if P. aeruginosa....
  df_to_fill_bed$contig_n <- str_sub(df_to_fill_bed$contig, -2, -1)
  df_to_fill_bed$contig_n <- gsub("_","",df_to_fill_bed$contig_n)
  df_to_fill_bed$contig_n <- as.numeric(df_to_fill_bed$contig_n)
  keep_contigs <- df_to_fill_bed %>%
    group_by(contig_n) %>%
    dplyr::summarize(max=max(position)) %>%
    dplyr::filter(max>=100000) %>%
    dplyr::select(contig_n) %>%
    distinct()
  df_to_fill_bed <- df_to_fill_bed %>%
    dplyr::filter(contig_n %in% keep_contigs$contig_n)
  
}


text <- df_to_fill_bed %>%
  group_by(library) %>% 
  dplyr::mutate(mean=mean(coverage),
                sd=sd(coverage),
                zeros=NROW(which(coverage==0))) %>%
  dplyr::select(library,mean,sd,zeros) %>%
  distinct() %>%
  dplyr::mutate(label=paste0("mean=",round(mean,2),
                             "\n sd=",round(sd,2),
                             "\n zero=", zeros))


p1 <- df_to_fill_bed %>% 
  group_by(library,coverage) %>%
  tally() %>%
  ggplot(., aes(x=coverage,y=n, color=library))+
  geom_line() +
  xlim(0,100) +
  theme(legend.position="top",
        legend.title = element_blank())+
  labs(x="coverage per site",
       y="Frequency (# sites)") +
  theme_bw()+
  theme(legend.position="none")

p2 <- df_to_fill_bed %>%
  dplyr::filter(coverage < 4) %>%
  ggplot(., aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=0.5)+
  facet_grid(rows = vars(library)) +
  theme(legend.position="top",
        legend.title = element_blank())+
  xlim(-0.5,5.5) +
  labs(x="coverage per site",
       y="Frequency (# sites)") +
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

# display correlation by coverage:

bed_corr <- df_to_fill_bed
head(bed_corr)

bed_corr <- bed_corr %>%
  dplyr::select(contig,position,library,coverage) %>%
  pivot_wider(names_from=library, values_from=coverage, values_fill = 0)

tail(bed_corr)

bed_corr_libs <- bed_corr[,3:ncol(bed_corr)]


# compute correlation between libs based on their coverage (all contigs)
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res <- rcorr(as.matrix(bed_corr_libs), type = "pearson")
corr_df <- flattenCorrMatrix(res$r, res$P) %>%
  arrange(desc(cor))


fwrite(x=corr_df, file=paste0(mydir,paste0(species,"_",goal,"_coverage_correlation.csv")))


########################################

# display top largest contigs' lowest coverage areas 


# keep n largest contigs
if (species=="Pa") {  
  
  keep_largest_contigs <- df_to_fill_bed %>%
    group_by(contig_n) %>%
    dplyr::summarize(max=max(position)) %>%
    dplyr::filter(max>=500000) %>%
    dplyr::select(contig_n) %>%
    distinct()
  df_to_fill_bed <- df_to_fill_bed %>%
    dplyr::filter(contig_n %in% keep_largest_contigs$contig_n)
  
}


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
  geom_point(aes(shape=library, color=library), size=0.01)+
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

# check sites with zero coverage: do the genomic positions overlap? 
x0 <- df_to_fill_bed %>%
  dplyr::filter(coverage < 1) %>% 
  dplyr::mutate(library=as.factor(recode_fun(library))) %>%
  dplyr::filter(library %in% my_subset) %>%
  drop.levels()

multiple_DFs <- split(x0, list(x0$contig, x0$library), drop = TRUE)

# empty df
lowcovDF <- data.frame(contig=character(),
                       start_position=numeric(),
                       end_position=numeric(),
                       library=character())

for (single_DF in multiple_DFs) {
  
  single <- as.data.frame(single_DF)
  lib <- unique(single$library)
  
  x0 <- single %>%
    dplyr::arrange(position)
  
  x0$bins <- c(0, cumsum(diff(x0$position) != 1))
  
  x0 <- x0 %>%
    drop.levels() %>%
    group_by(contig,bins) %>%
    dplyr::summarise(start_position=min(position),
                     end_position=max(position), .groups='drop') %>%
    dplyr::select(contig,start_position,end_position)
  
  x0$library <- paste0(lib)
  
  lowcovDF <- rbind(lowcovDF,x0)
  
}

lowcovDF <- lowcovDF %>%
  dplyr::arrange(contig,start_position) %>%
  dplyr::mutate(length_bp=end_position-start_position)

# for manuscript: 
lowcovDF %>%
  group_by(library) %>%
  dplyr::summarise(sum=sum(length_bp))

fwrite(x=lowcovDF, file=paste0(mydir,paste0(species,"_",goal,"_zero_cov_sites.csv")))


########################################
########################################


# Base quality stats (reads from bam - BQ computed with ALFRED)

BQ_files = list.files(mydir,pattern="^BQ")

datalist = list()
for (BQ_file in BQ_files) {
  
  # open file contaning all the PHRED scores:
  BQ_df <- read_delim(file.path(mydir,BQ_file), 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
  
  head(BQ_df)
  
  
  id <- BQ_file
  
  id <- sub(".*interleaved_", "",id)
  id <- gsub(".dedup.tsv.tsv","",id)
  id <- recode_fun(id)
  
  BQ_df$library <- paste0(id)
  
  BQ_df <- BQ_df %>%
    dplyr::select(Position,BaseQual,Read,library)
  
  # subset
  BQ_df <- BQ_df %>%
    dplyr::filter(library %in% my_subset) %>%
    drop.levels()
  
  datalist[[BQ_file]] <- BQ_df # add it to your list
  
  
}


BQ_data = do.call(rbind, datalist)

bq_data_stats <- BQ_data %>%
  group_by(library, Read) %>%
  dplyr::summarise(median_BQ=median(BaseQual),
                   mean_BQ=mean(BaseQual),
                   sd_BQ=min(BaseQual),
                   bp_count=n())


fwrite(x=bq_data_stats, file=paste0(mydir,paste0(species,"_",goal,"_base_quality.csv")))

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


# # ME data save as png : 
# kbl(ME_data, format = "html") %>%
#   kable_classic() %>%
#   kable_styling(font_size = 15, "striped") %>%
#   kableExtra::as_image(file=paste0(mydir,species,"_ALFRED_data.png")) 


# Save as Table : 
fwrite(x=ME_data, file=paste0(mydir,paste0(species,"_",goal,"_mapping.csv")))

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
  
  if (NROW(kra)>1) {
    
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
  
  
}

# subset
df_to_fill_kraken <- df_to_fill_kraken %>% dplyr::filter(library %in% my_subset)

# re-order libs
df_to_fill_kraken$library <- as.factor(df_to_fill_kraken$library)
df_to_fill_kraken <- reorder_lib_fun(df_to_fill_kraken)

# Save as Table : 
fwrite(x=df_to_fill_kraken, file=paste0(mydir,species,"_",goal,"_kraken.csv"))


########################################
########################################


# Assembly: comparison of HF vs HF_0.6x ( double left clean up )

# if this is "goal_size_selection" read in the assembly stats, parse and save 
if (grepl("size_selection", mydir, fixed = TRUE)==TRUE) {
  
  print ("We look at assembly stats, as these two libs are compared to assess the effect of size selection on assembly")
  
  assembly_stats <- read_delim(paste0(mydir,"assembly_stats.txt"), 
                               "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
  
  id <- assembly_stats$filename
  
  id <- sub(".*interleaved_", "", id)
  id <- gsub("/contigs.fa","",id)
  id
  
  id <- recode_fun(id)
  assembly_stats$library=paste0(as.character(id))
  
  assembly_stats <- assembly_stats %>%
    dplyr::select(library, everything())
  
  # Save as Table : 
  fwrite(x=assembly_stats, file=paste0(mydir,species,"_",goal,"_assembly.csv"))
  
}


################################################################################
################################################################################
################################################################################



# plot
pdf(paste0(mydir,species,"_",goal,'_out.pdf'))
# plot coverage 
cov_plot
# print lowest coverage regions
lowest_cov_plot
# plot insert size 
insert_size_plot 
#straight_GC
GC_picard
# plot PHRED scores
PHRED_plot
# plot read lengths
RL_plot
dev.off()

################################################################################
################################################################################


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
