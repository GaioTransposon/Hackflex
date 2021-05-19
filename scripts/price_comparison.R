#Price comparison:

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readxl)

Supplementary_Table_1 <- read_excel("~/Hackflex/source_data/Supplementary_Table_1.xlsx", 
                                    sheet = "price_fold_diff")

out_dir <- "~/Desktop/MG1655/"

head(Supplementary_Table_1)


plots_fold <- Supplementary_Table_1 %>%
  dplyr::select(`samples (n)`,SF_vs_Hfv0,SF_vs_Hfv1) %>%
  pivot_longer(cols=c(SF_vs_Hfv0,SF_vs_Hfv1), names_to="comparison",values_to="fold change") %>%
  ggplot(., aes(x=`samples (n)`, y=`fold change`, fill=comparison, color=comparison))+
  geom_point(size=1,aes(shape = comparison)) +
  scale_shape_manual(values=c(3, 17)) +
  geom_line() +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position="top",
        axis.title=element_text(size=18),
        axis.text=element_text(size=14))
  

pdf(paste0(out_dir,"plots_fold.pdf"))
plots_fold
dev.off()