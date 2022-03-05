########################################################################
# "urssa" Unsupervised Routine Soil Spectral Analysis                  #
# Link: https://github.com/raulpoppiel/urssa                           #
# When using any part of this code, please cite:                       #
# Poppiel, R.R.; Paiva, A.F.S.; Demattê, J.A.M. Bridging the gap       #
# between soil spectroscopy and tradition-al laboratory: insights for  #
# routine implementation. Geoderma, 2022. DOI                          #
########################################################################


## Initial Setting up
options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(ggplot2)
library(ggpubr)

# Adjusting the directory for loading and saving files
setwd("E:/Doutorado FAV-UnB/Parceria/Ariane - ProBASE/Artigo/Topico1/AnaliseDeResultados/Final_code/urssa")
getwd()

# Import outlier data
Observation_outliers <- read.csv("Soil_data_spectraResampled_clustered_outliers.csv", header = TRUE, sep = ";", dec = ".", row.names = "ID", na.strings = "NA",check.names=FALSE)

# Summarize number of outliers by Stratification variable and soil attribute
summarised.Clay <- Observation_outliers %>% select(c(Var_strat, Clay.is.outlier)) %>% filter(Clay.is.outlier %in% "TRUE") %>% 
  group_by(Var_strat) %>% summarise(Clay.is.outlier=n())

summarised.Sand <- Observation_outliers %>% select(c(Var_strat, Sand.is.outlier)) %>% filter(Sand.is.outlier %in% "TRUE") %>% 
  group_by(Var_strat) %>% summarise(Sand.is.outlier=n())

summarised.OM <- Observation_outliers %>% select(c(Var_strat, OM.is.outlier)) %>% filter(OM.is.outlier %in% "TRUE") %>% 
  group_by(Var_strat) %>% summarise(OM.is.outlier=n())

summarised.CEC <- Observation_outliers %>% select(c(Var_strat, CEC.is.outlier)) %>% filter(CEC.is.outlier %in% "TRUE") %>% 
  group_by(Var_strat) %>% summarise(CEC.is.outlier=n())

summarised.V <- Observation_outliers %>% select(c(Var_strat, V.is.outlier)) %>% filter(V.is.outlier %in% "TRUE") %>% 
  group_by(Var_strat) %>% summarise(V.is.outlier=n())

# Merge summarized outliers
summarised.merged <- merge(merge(merge(merge(summarised.Clay,summarised.Sand,by="Var_strat",all=TRUE),
                                       summarised.OM,by="Var_strat",all=TRUE),
                                 summarised.CEC,by="Var_strat",all=TRUE),
                           summarised.V,by="Var_strat",all=TRUE)

summarised.merged <- summarised.merged %>% rename_with(stringr::str_replace, pattern = ".is.outlier", replacement = "", matches(".is.outlier"))

# Replace NAs with zero values (0)  
summarised.merged <- summarised.merged %>% mutate_all(funs(replace_na(.,0)))
colSums(summarised.merged %>% select(Clay:V)) # total number of outliers by soil attribute (adjust the soil attributes selected)
summarised.merged %>% mutate(sum.by.Var_strat = select(., 2:6) %>% rowSums(na.rm = F)) %>% arrange(sum.by.Var_strat) # total number of outliers by Var_strat (laboratory)

# Defining the theme of the plot
themeBasic <- function () { 
  theme_classic2(base_size = 10) %+replace% 
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) + theme_update(plot.title = element_text(hjust = 0.5))
}

# Gathering columns into pairs for plotting 
summarised.merged.gathered <- summarised.merged %>% gather("soil.attr", "n.out", 2:6)

# Order data frame by Rowname and Laboratory
#summarised.merged.gathered$Laboratorio <- str_pad(summarised.merged.gathered$Laboratorio, 2, pad = "0")
summarised.merged.gathered <- summarised.merged.gathered %>% group_by(Var_strat) %>% mutate(orderby=sum(n.out))
summarised.merged.gathered <-summarised.merged.gathered[with(summarised.merged.gathered, order(orderby)),]

# Calculating the percent of total numb of outlier by attribute
SoilAttrib <- "OM"
total <- Observation_outliers %>% select(SoilAttrib) %>% summarise(n=n())
out <- Observation_outliers %>% select(paste0(SoilAttrib,".is.outlier")) %>% 
  rename(var.selected=paste0(SoilAttrib,".is.outlier")) %>%
  filter(var.selected %in% "TRUE") %>% 
  summarise(var.selected=n())
round((out/total)*100,2) # Percent of outlier by attribute

# Calculating the percent of total numb of outlier by Laboratory
total.out.by.lab <- cbind(summarised.merged, total.out = rowSums(summarised.merged))
total.out.by.lab[order(total.out.by.lab$total.out),]
total.out.by.lab[order(total.out.by.lab$V),]

# Grouped barplot
tiff(paste0("Outliers_barplot.tif"), width = 2200, height = 1500, res = 300)
ggplot(summarised.merged.gathered, aes(x=reorder(Var_strat, orderby), y=n.out, fill=soil.attr)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) + scale_y_continuous(expand=c(0,0)) +
  guides(fill=guide_legend("Soil attribute")) +
  ggtitle("") +
  themeBasic() +
  theme_classic(base_size = 12) + 
  xlab("Laboratory (Var_strat)") + 
  ylab("Number of outliers detected")


dev.off()

gc()            # clear the memory
rm(list=ls())   # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console
