#Price comparison:

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("~/Hackflex/source_data")

out_dir <- "~/Desktop/MG1655/"

price_per_samples_to_R <- read_csv("price_per_sample.csv")
head(price_per_samples_to_R)

plot_fold1 <- price_per_samples_to_R %>%
  dplyr::select(samples,method,`price AUD`) %>%
  pivot_wider(names_from=method,values_from=`price AUD`) %>%
  dplyr::mutate(`Standard Flex vs Hackflex`=Hackflex/`Nextera Flex`,
                `Standard Flex vs Hackflex 1/2 polymerase`=`Hackflex 1/2 polymerase`/`Nextera Flex`) %>%
  dplyr::select(samples,`Standard Flex vs Hackflex`,`Standard Flex vs Hackflex 1/2 polymerase`) %>%
  pivot_longer(cols=c(`Standard Flex vs Hackflex`,
                      `Standard Flex vs Hackflex 1/2 polymerase`), values_to="fold_diff", names_to="comparison") %>%
  ggplot(., aes(x=samples, y=fold_diff, fill=comparison, color=comparison))+
  geom_point(size=1,aes(shape = comparison)) +
  scale_shape_manual(values=c(3, 17)) +
  geom_line() +
  labs(x="samples",
       y="price fold difference") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position="top",
        axis.title=element_text(size=18),
        axis.text=element_text(size=14))+
  xlim(0,5010)+
  ylim(0,1)

plot_fold2 <- price_per_samples_to_R %>%
  dplyr::select(samples,method,`price AUD`) %>%
  pivot_wider(names_from=method,values_from=`price AUD`) %>%
  dplyr::mutate(`Standard Flex vs Hackflex`=`Nextera Flex`/Hackflex,
                `Standard Flex vs Hackflex 1/2 polymerase`=`Nextera Flex`/`Hackflex 1/2 polymerase`) %>%
  dplyr::select(samples,`Standard Flex vs Hackflex`,`Standard Flex vs Hackflex 1/2 polymerase`) %>%
  pivot_longer(cols=c(`Standard Flex vs Hackflex`,
                      `Standard Flex vs Hackflex 1/2 polymerase`), values_to="fold_diff", names_to="comparison") %>%
  ggplot(., aes(x=samples, y=fold_diff, fill=comparison, color=comparison))+
  geom_point(size=1,aes(shape = comparison)) +
  scale_shape_manual(values=c(3, 17)) +
  geom_line() +
  labs(x="samples",
       y="price fold difference") +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position="top",
        axis.title=element_text(size=18),
        axis.text=element_text(size=14))+
  xlim(0,5010)+
  ylim(0,11.5)
  

plots_fold <- ggarrange(plot_fold1, plot_fold2, nrow=2)

pdf(paste0(out_dir,"plots_fold.pdf"))
plots_fold
dev.off()