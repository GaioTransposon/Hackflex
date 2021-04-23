# find all .csv files in all dirs and sub dirs :
# goal_ecoli, goal_paeruginosa, goal_saureus, goal_barcode


# save all csv files as a single stats.xlsx file


###########################################################################################

library(readr)
library(openxlsx)
library(gtools)



source_dir = "/Users/12705859/Desktop/MG1655"
out_dir = "/Users/12705859/Desktop/MG1655/"

these_dirs <- list.dirs(source_dir, recursive = FALSE)
these_dirs <- grep("goal", these_dirs, value = TRUE)


###########################################################################################

# 001. all the stats: 


## 1. create workbook 
wb <- createWorkbook()


## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern=".csv", full.names=TRUE, recursive = TRUE)

for (each_filename in filenames) {
  
  ldf <- as.data.frame(lapply(each_filename, read.csv, header = TRUE))
  
  clean_name <- basename(each_filename)
  clean_name <- gsub(".csv", "", clean_name)
  
  addWorksheet(wb, clean_name)
  writeData(wb, sheet = clean_name, ldf, colNames = TRUE)
  
}


## 3. save workbook 
saveWorkbook(wb, paste0(out_dir,"stats.xlsx"), overwrite=TRUE)


###########################################################################################

# 002. Table 1: 

my_subset <- c("Ec.SF_1.B1",
               "Ec.SF_1:50.B1",
               "Ec.SF_1.B2",
               "Ec.SF_1:50.B2",
               "Ec.HF.B3", # sel from E. coli
               "Pa.SF_1.B1",
               "Pa.SF_1:50.B1",
               "Pa.HF.B2", # sel from P. aeruginosa
               "Sa.SF_1.B1",
               "Sa.SF_1:50.B1",
               "Sa.HF.B2") # sel from S. aureus

# reorder lib factor - function 
myord <-  c("Ec.SF_1.B1",
            "Ec.SF_1:50.B1",
            "Ec.SF_1.B2",
            "Ec.SF_1:50.B2",
            "Ec.HF.B3",
            "Pa.SF_1.B1",
            "Pa.SF_1:50.B1",
            "Pa.HF.B2",
            "Sa.SF_1.B1",
            "Sa.SF_1:50.B1",
            "Sa.HF.B2")


## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="phred.csv", full.names=TRUE, recursive = TRUE)

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
    dplyr::select(library, PHRED_mean, PHRED_sd) 
  
  
}

phred




## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="bbduk.csv", full.names=TRUE, recursive = TRUE)

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
    dplyr::select(library,tot_gen_reads,perc_removed_reads)
  
  
}

bbduk





## 2. open csv and save as sheets of the workbook
filenames <- list.files(these_dirs, pattern="all_lib_processing_stats.csv", full.names=TRUE, recursive = TRUE)

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
    dplyr::select(library,bp_after_resizing)
  
  
  
  
}

cleaning_stats

Table1 <- inner_join(inner_join(phred,bbduk), cleaning_stats) %>%
  dplyr::mutate(PHRED_mean=round(PHRED_mean,2),
                PHRED_sd=round(PHRED_sd,2))

fwrite(x = Table1, file = paste0(out_dir,"Table_1.csv"))



###########################################################################################

# 003. Table 2: 


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
    dplyr::filter(ldf[,ncol(ldf)] %in% c("library","X.Unmapped","MappedFraction","MismatchRate","MedianCoverage","SDCoverage")) 

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
cols.num <- c("X.Unmapped","MappedFraction","MismatchRate","MedianCoverage","SDCoverage")
ld[cols.num] <- sapply(ld[cols.num],as.numeric)

# rounding and ordering
Table_2 <- ld %>% dplyr::filter(library %in% my_subset) %>%
  mutate(library =  factor(library, levels = myord)) %>%
  arrange(library) %>%
  dplyr::mutate(MappedFraction=round(MappedFraction,3),
                MismatchRate=round(MismatchRate,3),
                SDCoverage=round(SDCoverage,2))

fwrite(x = Table_2, file = paste0(out_dir,"Table_2.csv"))



###########################################################################################

# 004. Table 3: 


## Insert size
filenames <- list.files(these_dirs, pattern="insert_size.csv", full.names=TRUE, recursive = TRUE)

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
filenames <- list.files(these_dirs, pattern="GC_bias.csv", full.names=TRUE, recursive = TRUE)

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

Table_3 <- inner_join(ld1,ld2)

fwrite(x = Table_3, file = paste0(out_dir,"Table_3.csv"))


###########################################################################################



