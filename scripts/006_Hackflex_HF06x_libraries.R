

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

########################################

mydirs <- "~/Desktop/MG1655/goal_size_selection"
my_subset <- c("Ec.HF.B3",
               "Ec.HF_06x.B3",
               "Pa.HF.B2",
               "Pa.HF_06x.B3",
               "Sa.HF.B2",
               "Sa.HF_06x.B3")
out_dir = "/Users/12705859/Desktop/MG1655/goal_size_selection/"

these_dirs <- list.dirs(mydirs, recursive = FALSE)
these_dirs <- grep("goal_size_selection", these_dirs, value = TRUE)


########################################


recode_fun <- function(chars) {
  
  x <- recode(chars, K13 = "Ec.HF.B3")
  x <- recode(x, K14 = "Ec.HF_06x.B3")
  
  x <- recode(x, K8 = "Pa.HF.B2")
  x <- recode(x, K15 = "Pa.HF_06x.B3")
  
  x <- recode(x, K7 = "Sa.HF.B2")
  x <- recode(x, K16 = "Sa.HF_06x.B3")
  
  return(x)
  
}

# reorder lib factor - function 
reorder_lib_fun <- function(df) {
  
  df$library  = factor(df$library, levels=c("Ec.HF.B3",
                                            "Ec.HF_06x.B3",
                                            "Pa.HF.B2",
                                            "Pa.HF_06x.B3",
                                            "Sa.HF.B2",
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

# add prefix to output files, if the goal of the analysis is to compare HF size selection libs with HF: 
if (grepl("goal_size_selection", mydirs, fixed = FALSE)==TRUE) {
  goal <- "SizeSel"
} else {
  goal <- "all"
}

########################################
########################################

# Insert size


filenames <- list.files(these_dirs, pattern="picardIS", full.names=TRUE, recursive = TRUE)
IS_files <- grep(filenames, pattern='.txt', value=TRUE)


# construct an empty dataframe to build on 
df_to_fill_insert_size <- data.frame(
  insert_size = numeric(),
  library = character(),
  stringsAsFactors = FALSE
)

for (IS_file in IS_files) {
  
  # read in file 
  IS_df <- read_delim(file.path(IS_file), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 10)
  
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

head(df_to_fill_insert_size)

df_to_fill_insert_size$spp <- substr(df_to_fill_insert_size$library, 1, 2)
df_to_fill_insert_size$type <- gsub("^.+?\\.(.+?)\\..*$", "\\1", df_to_fill_insert_size$library)

df_to_fill_insert_size <- df_to_fill_insert_size%>%
  dplyr::arrange(library)

pdf(paste0(out_dir,"size_sel_IS.pdf"))
df_to_fill_insert_size %>% 
  dplyr::select(library,insert_size,spp, type) %>%
  ggplot(., aes(insert_size, group=type)) + 
  facet_grid(rows = vars(spp), scales = "free") +
  geom_line(stat='density', aes(linetype = type, color=type), size = 1.2) +
  scale_color_manual(values = c("darkgreen", "green"))+
  geom_line(stat='density', size = 0.2, color="darkgreen") +
  scale_linetype_manual(values=c("dashed", "dotted")) +
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
                                    angle = 0)) +
  labs(x="insert size (bp)")
dev.off()


########################################
########################################


# Assembly stats: 


filenames <- list.files(these_dirs, pattern="assembly_stats", full.names=TRUE, recursive = TRUE)
assembly_files <- grep(filenames, pattern='.txt', value=TRUE)


datalist = list()
for (assembly_file in assembly_files) {
  
  # read in file 
  ass_df <- read_delim(file.path(assembly_file), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
  
  downsampling <- sub(".*shovillsub_", "", ass_df$filename)
  downsampling <- gsub("_R1.*","",downsampling)
  
  id <- ass_df$filename
  id <- sub(".*interleaved_", "", id)
  id <- gsub("/contigs.fa","",id)
  id
  
  id <- recode_fun(id)
  ass_df$library=paste0(as.character(id))
  ass_df$downsampling <- as.numeric(downsampling)
  
  ass_df <- ass_df %>%
    dplyr::select(library, everything())
  
  datalist[[assembly_file]] <- ass_df # add it to your list
  
  
}


assembly_data = do.call(rbind, datalist)
assembly_data$spp <- substr(assembly_data$library, 1, 2)
assembly_data$type <- gsub("^.+?\\.(.+?)\\..*$", "\\1", assembly_data$library)
assembly_data <- assembly_data %>%
  dplyr::arrange(library)

head(assembly_data)

z <- assembly_data %>%
  dplyr::select(downsampling,library,scaf_L50,scaf_N50) %>%
  dplyr::filter(downsampling<60)


contig_bp <- assembly_data %>%
  dplyr::filter(downsampling<50) %>%
  ggplot(., aes(x=downsampling, y= contig_bp)) + 
  scale_color_manual(values = c("blue", 
                                "magenta1", 
                                "green"))+ 
  geom_line(aes(colour = factor(spp),
                linetype = factor(type))) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20), 
        strip.text.y = element_text(size = 15, 
                                    colour = "black", 
                                    angle = 270)) +
  labs(x="downsampling size (coverage depth)",
       y="assembly size (bp)")


L50 <- assembly_data %>%
  dplyr::filter(downsampling<75) %>%
  ggplot(., aes(x=downsampling, y= ctg_N50)) + 
  scale_color_manual(values = c("blue", 
                                "magenta1", 
                                "green"))+ 
  geom_line(aes(colour = factor(spp),
                linetype = factor(type))) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20), 
        strip.text.y = element_text(size = 15, 
                                    colour = "black", 
                                    angle = 270)) +
  labs(x="downsampling size (coverage depth)",
       y="assembly L50")

N50 <- assembly_data %>%
  dplyr::filter(downsampling<200) %>%
  ggplot(., aes(x=downsampling, y= ctg_L50)) + 
  scale_color_manual(values = c("blue", 
                                "magenta1", 
                                "green"))+ 
  geom_line(aes(colour = factor(spp),
                linetype = factor(type))) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=15),
        axis.title=element_text(size=20), 
        strip.text.y = element_text(size = 15, 
                                    colour = "black", 
                                    angle = 270)) +
  labs(x="downsampling size (coverage depth)",
       y="assembly N50")



theme_smaller <- theme(axis.text=element_text(size=9),
                       axis.title=element_text(size=11), 
                       legend.position="none")

all <- ggarrange(contig_bp + theme_smaller + rremove("xlab"), 
                 N50 + theme_smaller + rremove("xlab"), 
                 L50 + theme_smaller + rremove("xlab"), 
                 nrow = 1, ncol=3, 
                 labels=c("B","C","D"))

all2 <- annotate_figure(all,
                        bottom = textGrob("downsampling depth", gp = gpar(cex = 1.1)))

all2


unique(df_to_fill_insert_size$spp)


IS_plot <- df_to_fill_insert_size %>% 
  dplyr::select(library,insert_size,spp, type) %>%
  ggplot(., aes(insert_size)) + 
  geom_line(stat='density', aes(linetype = type, color=spp)) +
  scale_color_manual(values = c("blue", 
                                "magenta1", 
                                "green"))+ 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text=element_text(size=10),
        axis.title=element_text(size=14), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(x="insert size (bp)")

IS_plot <- ggarrange(IS_plot, labels = c("A"))

pdf(paste0(out_dir,"size_sel_final.pdf"))
ggarrange(IS_plot,all2,nrow = 2, common.legend = TRUE)
dev.off()







# some info for manuscript: 

assembly_data %>%
  dplyr::filter(spp=="Pa") %>%
  dplyr::filter(downsampling=="50") %>%
  dplyr::select(library,ctg_L50,ctg_N50)

assembly_data %>%
  dplyr::filter(spp=="Sa") %>%
  dplyr::filter(downsampling=="10"|downsampling=="20"|downsampling=="50") %>%
  dplyr::select(library,ctg_L50,ctg_N50,downsampling)

assembly_data %>%
  dplyr::filter(spp=="Ec") %>%
  dplyr::filter(downsampling=="10"|downsampling=="20") %>%
  dplyr::select(library,ctg_L50,ctg_N50,downsampling)
