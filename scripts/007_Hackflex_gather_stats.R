# find all .csv files in all dirs and sub dirs :
# goal_ecoli, goal_paeruginosa, goal_saureus, goal_barcode


# save all csv files as a single stats.xlsx file


###########################################################################################

library(readr)
library(openxlsx)
library(gtools)
library(dplyr)
library(stringr)


source_dir = "/Users/12705859/Desktop/MG1655"
out_dir = "/Users/12705859/Desktop/MG1655/"

these_dirs <- list.dirs(source_dir, recursive = FALSE)
these_dirs <- grep("goal", these_dirs, value = TRUE)


###########################################################################################

# 002. Table 2: 

my_subset <- c("Ec.SF.B1",
               "Ec.SF_1:50.B1",
               "Ec.SF.B2",
               "Ec.SF_1:50.B2",
               "Ec.HF.B3", # sel from E. coli
               "Pa.SF.B1",
               "Pa.SF_1:50.B1",
               "Pa.HF.B2", # sel from P. aeruginosa
               "Sa.SF.B1",
               "Sa.SF_1:50.B1",
               "Sa.HF.B2") # sel from S. aureus

# reorder lib factor - function 
myord <-  c("Ec.SF.B1",
            "Ec.SF_1:50.B1",
            "Ec.SF.B2",
            "Ec.SF_1:50.B2",
            "Ec.HF.B3",
            "Pa.SF.B1",
            "Pa.SF_1:50.B1",
            "Pa.HF.B2",
            "Sa.SF.B1",
            "Sa.SF_1:50.B1",
            "Sa.HF.B2")


## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="phred.csv", full.names=TRUE, recursive = TRUE)

# do not include the size selection libs (HF vs HF0.6x) : 
filenames <- filenames[!str_detect(filenames,pattern="SizeSel")]

ld <- data.frame(
  library = character(),
  PHRED_mean = numeric(),
  PHRED_sd = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  ldf <- ldf %>% dplyr::filter(library %in% my_subset) 
  
  ld <- rbind(ld,ldf)
  ld$library <- as.factor(ld$library)
  
  phred <- ld %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) %>%
    dplyr::select(library, PHRED_mean, PHRED_sd) %>%
    distinct()
  
  
}

phred




## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="bbduk.csv", full.names=TRUE, recursive = TRUE)

# do not include the size selection libs (HF vs HF0.6x) : 
filenames <- filenames[!str_detect(filenames,pattern="SizeSel")]

ld <- data.frame(
  library = character(),
  tot_gen_reads = numeric(),
  perc_removed_reads = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  ldf <- ldf %>% dplyr::filter(library %in% my_subset)
  
  ldf1 <- ldf %>% dplyr::filter(X1 %in% c("Total Removed")) %>%
    dplyr::select(perc_reads,library) %>%
    dplyr::rename(., perc_removed_reads = perc_reads)
  
  ldf2 <- ldf %>% dplyr::filter(X1 %in% c("Result","Total Removed")) %>%
    group_by(library) %>%
    dplyr::summarise(tot_gen_reads=sum(reads))
  
  ldf3 <- inner_join(ldf1,ldf2) %>%
    dplyr::select(library, everything())
  
  ld <- rbind(ld,ldf3)
  ld$library <- as.factor(ld$library)
  
  bbduk <- ld %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) %>%
    dplyr::select(library,tot_gen_reads,perc_removed_reads) %>%
    distinct()
  
  
}

bbduk



## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="lib_proc_stats.csv", full.names=TRUE, recursive = TRUE)

# do not include the size selection libs (HF vs HF0.6x) : 
filenames <- filenames[!str_detect(filenames,pattern="SizeSel")]

ld <- data.frame(
  library = character(),
  bp_after_resizing = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  ldf <- ldf %>% dplyr::filter(library %in% my_subset)
  
  ldf$library <- as.factor(ldf$library)
  
  ldf1 <- ldf %>%
    dplyr::select(library,X..bp.after_resizing) %>%
    dplyr::rename(., bp_after_resizing = X..bp.after_resizing)
  
  
  ld <- rbind(ld,ldf1)
  ld$library <- as.factor(ld$library)
  
  cleaning_stats <- ld %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) %>%
    dplyr::select(library,bp_after_resizing) %>%
    distinct()
  
  
  
  
}

cleaning_stats

Table2 <- inner_join(inner_join(phred,bbduk), cleaning_stats) %>%
  dplyr::mutate(PHRED_mean=round(PHRED_mean,2),
                PHRED_sd=round(PHRED_sd,2)) %>%
  distinct()

fwrite(x = Table2, file = paste0(out_dir,"Table_2.csv"))



###########################################################################################

# 003. Table 3: 


## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="all_mapping.csv", full.names=TRUE, recursive = TRUE)

ld <- data.frame(
  library = character(),
  X.Unmapped = numeric(),
  MappedFraction = numeric(),
  MismatchRate = numeric(),
  MedianCoverage = numeric(),
  SDCoverage = numeric(),
  stringsAsFactors = FALSE
)

# row to colnames function
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = FALSE))
  
  # filter rows based on last column
  ldf <- ldf %>%
    dplyr::filter(ldf[,ncol(ldf)] %in% c("library","X.Unmapped","MappedFraction","MismatchRate")) 

  ldf <- header.true(ldf)
  
  rownames(ldf) <- ldf$library
  ldf$library <- NULL
  
  ldf <- as.data.frame(x = t(ldf), stringsAsFactors = FALSE)
  ldf$library <- rownames(ldf)
  rownames(ldf) <- NULL
  
  ldf <- ldf %>% dplyr::select(library, everything())
  ld <- rbind(ld, ldf)

}

# class: character to numeric
cols.num <- c("X.Unmapped","MappedFraction","MismatchRate")
ld[cols.num] <- sapply(ld[cols.num],as.numeric)

# rounding and ordering
mapping_table <- ld %>% dplyr::filter(library %in% my_subset) %>%
  mutate(library =  factor(library, levels = myord)) %>%
  arrange(library) %>%
  dplyr::mutate(MappedFraction=round(MappedFraction,3),
                MismatchRate=round(MismatchRate,3))

#####


# add PCR duplicates info: 
dups_files <- list.files(these_dirs, pattern="all_dups.csv", full.names=TRUE, recursive = TRUE)

x0 <- data.frame(
  library = character(),
  PCR_duplicates = numeric(),
  stringsAsFactors = FALSE
)

for (dup in dups_files) {
  
  x <- as.data.frame(lapply(dup, read.csv, header = TRUE))
  
  x <- x %>% dplyr::filter(library %in% my_subset) 
  
  x0 <- rbind(x0,x)
  
}

dups <- x0 
head(dups)


#####


## mean and sd coverage (these have been extracted from bed files, keeping contigs larger than 100,000 bp)
filenames <- list.files(these_dirs, pattern="all_coverage.csv", full.names=TRUE, recursive = TRUE)

ld1 <- data.frame(
  library = character(),
  mean = numeric(),
  sd = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  
  ldf <- ldf %>% 
    dplyr::select(library,mean,sd) %>%
    dplyr::filter(library %in% my_subset) 
  
  ldf$library <- as.factor(ldf$library)
  
  
  ld1 <- rbind(ld1,ldf) %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) %>%
    dplyr::select(library, mean, sd) %>%
    dplyr::mutate(mean=round(mean,2),
                  sd=round(sd,2))
  
}

colnames(ld1) <- c("library","mean_coverage","sd_coverage")

mean_sd_cov <- ld1


Table_3 <- inner_join(inner_join(mapping_table,dups),mean_sd_cov) %>%
  dplyr::mutate(library =  factor(library, levels = myord)) %>%
  dplyr::arrange(library)


fwrite(x = Table_3, file = paste0(out_dir,"Table_3.csv"))


###########################################################################################

# 004. Table 4: 


## Insert size
filenames <- list.files(these_dirs, pattern="all_insert_size.csv", full.names=TRUE, recursive = TRUE)

ld1 <- data.frame(
  library = character(),
  median = numeric(),
  mean = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  

  ldf <- ldf %>% 
    dplyr::filter(library %in% my_subset) 
  
  ldf$library <- as.factor(ldf$library)
  
  
  ld1 <- rbind(ld1,ldf) %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) %>%
    dplyr::select(library, median, mean) 
  
}

colnames(ld1) <- c("library","median_IS","mean_IS")


## GC bias
filenames <- list.files(these_dirs, pattern="all_GC_bias.csv", full.names=TRUE, recursive = TRUE)

ld2 <- data.frame(
  library = character(),
  rho = numeric(),
  pval = numeric(),
  stringsAsFactors = FALSE
)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))

  ldf <- ldf %>% 
    dplyr::filter(library %in% my_subset) %>%
    dplyr::select(library, rho, pval)
  
  ldf$library <- as.factor(ldf$library)
  
  ld2 <- rbind(ld2,ldf) %>%
    mutate(library =  factor(library, levels = myord)) %>%
    arrange(library) 
  
}

colnames(ld2) <- c("library","GC_rho","GC_pval")

ld2 <- ld2 %>%
  dplyr::mutate(GC_rho = round(GC_rho,3),
                GC_pval = stars.pval(GC_pval))

Table_4 <- inner_join(ld1,ld2)

fwrite(x = Table_4, file = paste0(out_dir,"Table_4.csv"))


###########################################################################################


# Supplementary Table 2: 



## 1. create workbook 
wb <- createWorkbook()


#####


# bbduk output

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_bbduk")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "bbduk")
writeData(wb, sheet = "bbduk", df, colNames = TRUE)


#####

# phred scores 

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="phred")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "phred_scores")
writeData(wb, sheet = "phred_scores", df, colNames = TRUE)

#####
#####

# library processing stats

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_lib_proc_stats")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "lib_proc_stats")
writeData(wb, sheet = "lib_proc_stats", df, colNames = TRUE)

#####
#####

# base quality

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_base_quality")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "base_quality")
writeData(wb, sheet = "base_quality", df, colNames = TRUE)

#####
#####

# coverage

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_coverage")]
filenames <- filenames[!str_detect(filenames,pattern="correl")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "coverage")
writeData(wb, sheet = "coverage", df, colNames = TRUE)

#####
#####

# coverage

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_coverage_corr")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "coverage_corr")
writeData(wb, sheet = "coverage_corr", df, colNames = TRUE)

#####
#####

# PCR duplicates

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_dups")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "PCR_duplicates")
writeData(wb, sheet = "PCR_duplicates", df, colNames = TRUE)

#####
#####

# zero coverage sites

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_zero_cov_sites")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "zero_cov_sites")
writeData(wb, sheet = "zero_cov_sites", df, colNames = TRUE)

#####
#####

# contaminants

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_kraken")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "contaminants")
writeData(wb, sheet = "contaminants", df, colNames = TRUE)

#####
#####

# GC_bias

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_GC_bias")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "GC_bias")
writeData(wb, sheet = "GC_bias", df, colNames = TRUE)

#####
#####

# Insert size

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="all_insert_size")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "insert_size")
writeData(wb, sheet = "insert_size", df, colNames = TRUE)

#####
#####

# Assembly HF vs HF0.6x (i.e.: double left clean up)

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="assembly")]

datalist = list()

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  datalist[[each_filename]] <- ldf # add it to your list
  
}

df <- do.call(rbind, datalist)
rownames(df) <- NULL

addWorksheet(wb, sheetName = "assembly")
writeData(wb, sheet = "assembly", df, colNames = TRUE)

#####
#####

# Barcodes: 

filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)
filenames <- filenames[str_detect(filenames,pattern="goal_barcode")]
filenames <- filenames[!str_detect(filenames,pattern="all_possibilities")]

# open each and save as sheets of the workbook
for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  clean_name <- basename(each_filename)
  clean_name <- gsub(".csv", "", clean_name)
  
  addWorksheet(wb, sheetName = clean_name)
  writeData(wb, sheet = clean_name, ldf, colNames = TRUE)
  
}

#####

# save workbook 
saveWorkbook(wb, paste0(out_dir,"stats.xlsx"), overwrite=TRUE)

