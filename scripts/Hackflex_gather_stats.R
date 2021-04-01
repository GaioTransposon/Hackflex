# find all .csv files in all dirs and sub dirs :
# goal_ecoli, goal_paeruginosa, goal_saureus, goal_barcode


# save all csv files as a single stats.xlsx file


###########################################################################################

library(readr)
library(openxlsx)



source_dir = "/Users/12705859/Desktop/MG1655"
out_dir = "/Users/12705859/Desktop/MG1655/"

these_dirs <- list.dirs(source_dir, recursive = FALSE)
these_dirs <- grep("goal", these_dirs, value = TRUE)





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


