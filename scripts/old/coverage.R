library(readr)

# set directories
middle_dir <- "~/Desktop/MG1655/"
out_dir <- "~/Desktop/MG1655/"

########################################
########################################

cov_files <- grep(list.files(middle_dir,pattern="^reduced_trimmed2_trimmed_interleaved2_"), 
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
  
  coverage <-read.table(file.path(middle_dir,cov_file))
  
  colnames(coverage) <- c("position","coverage")
  
  id <- str_replace_all(cov_file, "reduced_trimmed2_trimmed_interleaved2_", "")
  coverage$library=paste0(as.character(id))
  
  df_to_fill <- rbind(df_to_fill,coverage)
  
}
  
head(df_to_fill)
NROW(df_to_fill)

fun_cov_plot <- function(x) {
  hist(x$coverage[x$coverage<100], 
       breaks=100,
       xlim=c(0,100), 
       main=as.character(unique(x$library)),
       xlab = "coverage per site", 
       cex.lab = 1.5)
}

fun_low_cov_plot <- function(x) {
  hist(x$coverage[x$coverage<=4], 
       breaks=4,
       xlim=c(0,4), 
       main=NULL,
       xlab = NULL,
       ylab=NULL,
       cex.axis = 0.8)
}

#plot
pdf(paste0(out_dir,'coverage.pdf'))
par(mfrow=c(3,1))
fun_cov_plot(subset(df_to_fill, library == "K1_R1.dedup.tsv"))
fun_cov_plot(subset(df_to_fill, library == "K3_R1.dedup.tsv"))
fun_cov_plot(subset(df_to_fill, library == "K6_R1.dedup.tsv"))
par(fig=c(0.61, 0.99, 0.75, 1), new = T)
fun_low_cov_plot(subset(df_to_fill, library == "K1_R1.dedup.tsv"))
par(fig=c(0.61, 0.99, 0.4, 0.65), new = T)
fun_low_cov_plot(subset(df_to_fill, library == "K3_R1.dedup.tsv"))
par(fig=c(0.61, 0.99, 0.1, 0.35), new = T)
fun_low_cov_plot(subset(df_to_fill, library == "K6_R1.dedup.tsv"))
ggplot(df_to_fill, aes(x=coverage, color=library)) +
  geom_histogram(fill="white", alpha=0.5,
                 position="identity",
                 binwidth=1)+
  facet_grid(rows = vars(library)) +
  xlim(0,100) +
  labs(x="coverage per site",
       y="Frequency")
dev.off()



