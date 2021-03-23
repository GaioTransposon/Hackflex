library(reshape2)
library(textshape)


R1_unpaired_headers <- read_delim("Desktop/MG1655/HF-barcode_undetermined/R1_undet_headers",":", 
                                  escape_double = FALSE, col_names = FALSE,
                                  trim_ws = TRUE)
R2_unpaired_headers <- read_delim("Desktop/MG1655/HF-barcode_undetermined/R2_undet_headers",":", 
                                  escape_double = FALSE, col_names = FALSE,
                                  trim_ws = TRUE)

head(R1_unpaired_headers)

R1_unpaired <- R1_unpaired_headers %>%
  dplyr::select(X10) %>%
  dplyr::mutate(i5 = sub("\\+.*", "", X10)) %>%
  dplyr::mutate(i7 = str_extract(X10, '\\b[^+]+$')) %>%
  dplyr::select(i5,i7)

R2_unpaired <- R2_unpaired_headers %>%
  dplyr::select(X10) %>%
  dplyr::mutate(i5 = sub("\\+.*", "", X10)) %>%
  dplyr::mutate(i7 = str_extract(X10, '\\b[^+]+$')) %>%
  dplyr::select(i5,i7)


NROW(R1_unpaired)
NROW(which(R1_unpaired$i5==R2_unpaired$i5))
NROW(which(R1_unpaired$i7==R2_unpaired$i7))


# keep i5 present in list: 
NROW(R1_unpaired)
unpaired <- R1_unpaired %>% 
  dplyr::filter(i5 %in% HF_barcode_i5_list) %>%
  dplyr::filter(i7 %in% HF_barcode_i7_list) #%>%
  #drop.levels()



u  <- as.data.frame(unpaired)
u$count <- 1

u_sum <- u %>%
  group_by(i5,i7) %>%
  dplyr::summarise(sum=sum(count))

NROW(unique(u_sum[,1:2]))
# it worked. 
u_sum <- as.data.frame(u_sum)
head(u_sum)




NROW(unique(q$i7))

q <- rbind(z3,u_sum) %>%
  group_by(i5,i7) %>%
  dplyr::summarise(sum=sum(sum))

NROW(unique(q[,1:2]))

q

u_summ <- acast(u_sum, i5~i7, value.var="sum", fill = 0)
qq <- acast(q, i5~i7, value.var="sum", fill = 0)

View(qq)

swap <- function(matrixRow,x,y){
  #x is diagonal index
  #y is max of the row
  indexY <- which(matrixRow == y)
  valX <- matrixRow[x]
  matrixRow[x] <- y
  matrixRow[indexY] <- valX
  return(matrixRow)
}

class(qq)

mat <- u_summ
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}




pdf(paste0(mydir,"test.pdf"))
mat <- u_summ
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}
mat %>%
  tidy_matrix('i5','i7') %>%
  ggplot(data = ., aes(x=i5, y=i7, fill=log(value))) +
  geom_tile() + 
  scale_fill_gradient(low = "yellow", high = "red", 
                      na.value = "white") +
  geom_text(aes(label = round(value, 1)), size=1) +
  theme(axis.text.x=element_text(angle=90))+
  theme(
    axis.text.y = element_text(size = 5)   ,
    axis.text.x = element_text(
      size = 5, 
      hjust = 1, 
      vjust = 1, 
      angle = 45
    ),   
    legend.position = 'bottom',
    legend.key.height = grid::unit(.1, 'cm'),
    legend.key.width = grid::unit(.5, 'cm'),
    legend.text = element_text(angle = 90,size=5)
  ) 
# clustered
mat <- u_summ
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}
mat %>%
  cluster_matrix() %>%
  tidy_matrix('i5', 'i7') %>%
  dplyr::mutate(
    i5 = factor(i5, levels = unique(i5)),
    i7 = factor(i7, levels = unique(i7))        
  ) %>%
  group_by(i5) %>%
  ggplot(aes(i5, i7, fill = value)) +
  geom_tile() +
  #geom_text(aes(label = round(value, 1)), size=2) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      na.value = "white") +
  theme(
    axis.text.y = element_text(size = 5)   ,
    axis.text.x = element_text(
      size = 5, 
      hjust = 1, 
      vjust = 1, 
      angle = 45
    ),   
    legend.position = 'bottom',
    legend.key.height = grid::unit(.1, 'cm'),
    legend.key.width = grid::unit(.5, 'cm'),
    legend.text = element_text(angle = 90,size=5)
  ) +
  labs(subtitle = "With Clustering")
# clustered
mat <- qq
for(i in 1:nrow(mat)){
  rowI <- mat[i,]
  y <- max(rowI)
  mat[i,] <- swap(rowI, i, y)
}
mat %>%
  cluster_matrix() %>%
  tidy_matrix('i5', 'i7') %>%
  dplyr::mutate(
    i5 = factor(i5, levels = unique(i5)),
    i7 = factor(i7, levels = unique(i7))        
  ) %>%
  group_by(i5) %>%
  ggplot(aes(i5, i7, fill = value)) +
  geom_tile() +
  #geom_text(aes(label = round(value, 1)), size=2) +
  scale_fill_gradient(low = "yellow", high = "red", 
                      na.value = "white") +
  theme(
    axis.text.y = element_text(size = 5)   ,
    axis.text.x = element_text(
      size = 5, 
      hjust = 1, 
      vjust = 1, 
      angle = 45
    ),   
    legend.position = 'bottom',
    legend.key.height = grid::unit(.1, 'cm'),
    legend.key.width = grid::unit(.5, 'cm'),
    legend.text = element_text(angle = 90,size=5)
  ) +
  labs(subtitle = "With Clustering")
dev.off()



demux_clean <- read_csv("Desktop/MG1655/HF-barcode_undetermined/demux_clean.tsv", col_names = FALSE)
demux_clean <- demux_clean[6:677,]
demux_clean$seq <- rep(seq(1,7), 96)

demux_clean <- as.data.frame(demux_clean)

counts <- demux_clean %>%
  dplyr::filter(seq=="6")

counts$X1 <- substring(counts$X1, 22, str_length(counts$X1))
counts$X1 <- gsub("<.*","",counts$X1)
counts <- as.data.frame(as.numeric(counts$X1))
colnames(counts) <- "counts"


hist(counts$counts)



