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

# Adjusting the directory for loading and saving files
setwd("E:/Doutorado FAV-UnB/Parceria/Ariane - ProBASE/Artigo/Topico1/AnaliseDeResultados/Final_code/urssa")
getwd()



### MODULE 1: Importing and pre-processing data ####

# Importing data, selecting the Laboratory and checking the format
Data_all <- read.csv("Soil_data_spectra.csv", header = TRUE, sep = ";", dec = ".", na.strings = "NA",check.names=FALSE)
Data_all <- Data_all %>% remove_rownames %>% column_to_rownames(var="ID") #Assigning row names from ID column

# Filter rows using a stratification variable (e.g. laboratory, client code, soil type, geology, vegetation cover, etc)
Data_lab <- Data_all %>% filter(Var_strat %in% c(1:2)) # e.g., from lab 1 to 2

# Separating the soil data (observations) and spectra into different (two) tables
Observation_lab <- Data_lab %>% select(c(Var_strat,Sample,Clay,Sand,OM,CEC.ph7,V)) # soil data

wavelength_min <- "350" # adjust the min spectral wavelength
wavelength_max <- "2500" # adjust the max spectral wavelength
Spectra_lab <- Data_lab %>% select(wavelength_min:wavelength_max) # soil spectra (original)
#colnames(Spectra_lab) <-  as.numeric(gsub("X", "", colnames(Spectra_lab))) #Removing "X" from column names

# Resampling the spectral resolution of soil spectra
new_resolution <- 10 # Define the new spectral resolution (in nanometers)
wavelength_resampled <- seq(from = wavelength_min, to = wavelength_max, by = new_resolution) # new band positions

library(prospectr)
spectra_resampled <- as.data.frame(resample(X = Spectra_lab, wav = as.numeric(colnames(Spectra_lab)), 
                                            new.wav = wavelength_resampled, interpol = "spline"))


# Plotting of the first ten original and resampled spectra 
matplot(x = as.numeric(colnames(Spectra_lab)), y = t(Spectra_lab)[,1:10], type = "l",
        lty = 1, lwd = 2, col = 'black', xlab = "Wavelength (nm)", ylab = "Reflectance Factor", main = "Original Spectra")

matplot(x = as.numeric(colnames(spectra_resampled)), y = t(spectra_resampled)[,1:10], type = "l",
        lty = 1, lwd = 2, col = 'red', xlab = "Wavelength (nm)", ylab = "Reflectance Factor", main = "Resampled Spectra")


# Merging the soil data (observations) with the resampled spectra into a single table
Observation_lab_resampled <- merge(Observation_lab,spectra_resampled,by="row.names",all.x=TRUE)
Observation_lab_resampled <- Observation_lab_resampled %>% remove_rownames %>% column_to_rownames(var="Row.names") %>% as.data.frame()


### MODULE 2: Identification of soil attribute outliers ####

## 1) Calculating proximity from URF, 
## 2) Multidimensional Scaling (MDS),
## 3) Clustering samples (K-means)
library(randomUniformForest)

Observation_lab_clustered <- list()

for(i in sort(unique(Observation_lab_resampled$Var_strat))) {
  
  dataset.model <- Observation_lab_resampled %>%
    filter(Var_strat %in% i)
  
  results <- unsupervised.randomUniformForest(object = dataset.model %>% select(wavelength_min:wavelength_max), # select spectra
                                              baseModel = "proximity", # c("proximity", "proximityThenDistance", "importanceThenDistance")
                                              endModel = "MDSkMeans", # c("MDSkMeans", "MDShClust", "MDS", "SpectralkMeans")
                                              endModelMetric = NULL,
                                              samplingMethod = "with bootstrap", # c("uniform univariate sampling","uniform multivariate sampling", "with bootstrap")
                                              MDSmetric = "metricMDS", # c("metricMDS", "nonMetricMDS")
                                              proximityMatrix = NULL,
                                              sparseProximities = FALSE,
                                              outliersFilter = FALSE, 
                                              Xtest = NULL,
                                              predObject = NULL, 
                                              metricDimension = 2, 
                                              coordinates = c(1,2),
                                              bootstrapReplicates = 100,
                                              clusters = NULL,
                                              maxIters = NULL,
                                              importanceObject = NULL,
                                              maxInteractions = 2,
                                              reduceClusters = FALSE, 
                                              maxClusters = 10,
                                              mapAndReduce = FALSE,
                                              OOB = FALSE,
                                              subset = NULL, 
                                              seed = 2014,
                                              uthreads = "auto")
  
  results.final <- data.frame("ID" = rownames(dataset.model),
                              "cluster" = results$unsupervisedModel$cluster,
                              dataset.model)
  
  Observation_lab_clustered[[i]] <- results.final
  
}

# Combine all stratification variables (laboratories) into a single table
Observation_final <- do.call(rbind, Observation_lab_clustered)

# Merging the cluster code with soil data and the RESAMPLED spectra into a single table
Observation_lab_resampled_final <- merge(Observation_final %>% select(cluster),Observation_lab_resampled,by="row.names",all.x=TRUE)
Observation_lab_resampled_final <- Observation_lab_resampled_final %>% remove_rownames %>% column_to_rownames(var="Row.names") %>% as.data.frame()
# Saving the results
write.table(Observation_lab_resampled_final %>% rownames_to_column(var = "ID"),"Soil_data_spectraResampled_clustered.csv", sep = ";", dec = ".",row.names = FALSE)
# Merging the cluster code with soil data and the ORIGINAL spectra into a single table
Observation_lab_original_final <- merge(Observation_final %>% select(cluster),Data_all,by="row.names",all.x=TRUE)
Observation_lab_original_final <- Observation_lab_original_final %>% remove_rownames %>% column_to_rownames(var="Row.names") %>% as.data.frame()
# Saving the results
write.table(Observation_lab_original_final %>% rownames_to_column(var = "ID"),"Soil_data_spectraOriginal_clustered.csv", sep = ";", dec = ".",row.names = FALSE)

gc()            # clear the memory
rm(list=setdiff(ls(), c("wavelength_min","wavelength_max","Observation_lab_resampled","Data_all")))  # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console


# Soil attributes outlier detection using an adjusted boxplot (exponential model) for skewed distributions

library(robustbase)
# Reference: https://www.sciencedirect.com/science/article/pii/S0167947307004434
# https://cran.r-project.org/web/packages/robustbase/robustbase.pdf

Observation_final <- read.csv("Soil_data_spectraResampled_clustered.csv", header = TRUE, sep = ";", dec = ".", row.names = "ID", na.strings = "NA",check.names=FALSE)

Outliers <- Observation_final %>% rownames_to_column(var = "ID") %>%
                              group_by(Var_strat,cluster) %>% 
  mutate(lower.clay = adjboxStats(Clay, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1]) %>% 
  mutate(upper.clay = adjboxStats(Clay, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5]) %>%
  mutate(Clay.is.outlier = ifelse(Clay < lower.clay | Clay > upper.clay, TRUE, FALSE)) %>%
  select(-c(lower.clay,upper.clay))

Outliers <- Outliers %>% group_by(Var_strat,cluster) %>% 
  mutate(lower.Sand = adjboxStats(Sand, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1]) %>% 
  mutate(upper.Sand = adjboxStats(Sand, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5]) %>%
  mutate(Sand.is.outlier = ifelse(Sand < lower.Sand | Sand > upper.Sand, TRUE, FALSE)) %>%
  select(-c(lower.Sand,upper.Sand))

Outliers <- Outliers %>% group_by(Var_strat,cluster) %>% 
  mutate(lower.OM = adjboxStats(OM, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1]) %>% 
  mutate(upper.OM = adjboxStats(OM, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5]) %>%
  mutate(OM.is.outlier = ifelse(OM < lower.OM | OM > upper.OM, TRUE, FALSE)) %>%
  select(-c(lower.OM,upper.OM))

Outliers <- Outliers %>% group_by(Var_strat,cluster) %>% 
  mutate(lower.CEC.ph7 = adjboxStats(CEC.ph7, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1]) %>% 
  mutate(upper.CEC.ph7 = adjboxStats(CEC.ph7, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5]) %>%
  mutate(CEC.is.outlier = ifelse(CEC.ph7 < lower.CEC.ph7 | CEC.ph7 > upper.CEC.ph7, TRUE, FALSE)) %>%
  select(-c(lower.CEC.ph7,upper.CEC.ph7))

Outliers <- Outliers %>% group_by(Var_strat,cluster) %>% 
  mutate(lower.V_perc = adjboxStats(V, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1]) %>% 
  mutate(upper.V_perc = adjboxStats(V, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5]) %>%
  mutate(V.is.outlier = ifelse(V < lower.V_perc | V > upper.V_perc, TRUE, FALSE)) %>%
  select(-c(lower.V_perc,upper.V_perc))

# Summary of the number (n) of outlier identified by stratification variable (laboratory)
Outliers %>% group_by(Var_strat,cluster) %>% tally(Clay.is.outlier=="TRUE") %>% view()
Outliers %>% group_by(Var_strat,cluster) %>% tally(Sand.is.outlier=="TRUE") %>% view()
Outliers %>% group_by(Var_strat,cluster) %>% tally(OM.is.outlier=="TRUE") %>% view()
Outliers %>% group_by(Var_strat,cluster) %>% tally(CEC.is.outlier=="TRUE") %>% view()
Outliers %>% group_by(Var_strat,cluster) %>% tally(V.is.outlier=="TRUE") %>% view()

# Save the soil data that were flagged as outlier
# A column indicating the if the sample is an outlier (FALSE/TRUE) 
# will be created at the end of the table (last columns)
write.table(Outliers,"Soil_data_spectraResampled_clustered_outliers.csv", sep = ";", dec = ".",row.names = FALSE)

# Merging the outliers with the ORIGINAL spectra into a single table
Outliers_originalSpectra <- merge(Outliers %>% column_to_rownames(var="ID") %>%
                                        select(-c(wavelength_min:wavelength_max)),
                                        Data_all %>% select(c(wavelength_min:wavelength_max)),
                                        by="row.names",all.x=TRUE) %>% 
                                        remove_rownames %>% column_to_rownames(var="Row.names") %>% as.data.frame()
# Saving the results
write.table(Outliers_originalSpectra %>% rownames_to_column(var = "ID"),"Soil_data_spectraOriginal_clustered_outliers.csv", sep = ";", dec = ".",row.names = FALSE)


gc()            # clear the memory
rm(list=setdiff(ls(), c("wavelength_min","wavelength_max","Observation_lab_resampled")))  # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console



# Calculating the variable importance from the spectral clustering using URF

#rows<-c(seq(1,nrow(Observation_lab_resampled),3))+2
#length(rows)

matrix <- Observation_lab_resampled %>% select(wavelength_min:wavelength_max)

Spectra_clustered <- unsupervised.randomUniformForest(object = matrix,
                                                      baseModel = "proximity", # c("proximity", "proximityThenDistance", "importanceThenDistance")
                                                      endModel = "MDSkMeans", # c("MDSkMeans", "MDShClust", "MDS", "SpectralkMeans")
                                                      endModelMetric = NULL,
                                                      samplingMethod = "with bootstrap", # c("uniform univariate sampling","uniform multivariate sampling", "with bootstrap")
                                                      MDSmetric = "metricMDS", # c("metricMDS", "nonMetricMDS")
                                                      proximityMatrix = NULL,
                                                      sparseProximities = FALSE,
                                                      outliersFilter = FALSE, 
                                                      Xtest = NULL,
                                                      predObject = NULL, 
                                                      metricDimension = 2, 
                                                      coordinates = c(1,2),
                                                      bootstrapReplicates = 100,
                                                      clusters = NULL,
                                                      maxIters = NULL,
                                                      importanceObject = NULL,
                                                      maxInteractions = 2,
                                                      reduceClusters = FALSE, 
                                                      maxClusters = 10,
                                                      mapAndReduce = FALSE,
                                                      OOB = FALSE,
                                                      subset = NULL, 
                                                      seed = 2014,
                                                      uthreads = "auto")


##################################################################################
Spectra_clustered <- modifyClusters(Spectra_clustered, increaseBy = 1, seed = 2022)

## S3 method for class 'unsupervised'
print(Spectra_clustered)
## S3 method for class 'unsupervised'
plot(Spectra_clustered, importanceObject = NULL, xlim = NULL, ylim = NULL, coordinates = NULL)

# Merging clusters with all data
Observation_lab_clustered <- merge(Spectra_clustered$unsupervisedModel$cluster,Observation_lab_resampled,by="row.names",all.x=TRUE)
rownames(Observation_lab_clustered) <- Observation_lab_clustered$Row.names #Assigning row names from ID column
Observation_lab_clustered$Row.names <- NULL #Removing the ID column
Observation_lab_clustered <- Observation_lab_clustered %>% rename(Cluster = 'x') # Renaming the name of the column of clusters
###################################################################################

## as supervised
spectra.rufUnsup2sup = as.supervised(Spectra_clustered, matrix, mtry = 1, nodesize = 2)
## importance
spectra.ruf.importance = importance(spectra.rufUnsup2sup, Xtest = matrix, maxInteractions = 4)
spectra.ruf.importance.result <- spectra.ruf.importance$globalVariableImportance
write.table(spectra.ruf.importance.result,"spectra.ruf.importance.result.csv", sep = ";", dec = ".", row.names = FALSE)
spectra.ruf.importance.result <- read.csv("spectra.ruf.importance.result.csv", header = TRUE, sep = ";", dec = ".", na.strings = "NA")

spectra.ruf.importance.result$variables2=str_pad(spectra.ruf.importance.result$variables, 4, pad = "0")
spectra.ruf.importance.result=spectra.ruf.importance.result[order(spectra.ruf.importance.result$variables2),c(1,5)]
rownames(spectra.ruf.importance.result)=spectra.ruf.importance.result$variables
spectra.ruf.importance.result2 = spectra.ruf.importance.result %>% select(percent) %>% as.matrix() %>% t()
# Average spectrum 
spectrum <- matrix %>% summarise_all(.funs = c(mean="median"))
spectrum <- spectrum %>% rename_with(stringr::str_replace, 
                                     pattern = "_mean", replacement = "", 
                                     matches("_mean")) %>% as.matrix()

spectra_var.importance <- merge(as.data.frame(spectra.ruf.importance.result2),as.data.frame(spectrum),all=TRUE)
spectra_var.importance[3,]=as.numeric(colnames(spectra_var.importance))
spectra_var.importance=as.data.frame(t(spectra_var.importance))
colnames(spectra_var.importance)= c("espectrum","import","wl")

# Defining the theme of the plot
library(ggpubr)
themeBasic <- function () { 
  theme_classic2(base_size = 10) %+replace% 
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) + theme_update(plot.title = element_text(hjust = 0.5))
}

# Plotting the wavelength importance from URF overlaid by the averaged spectrum
tiff(paste0("spectra_var_import_URF.tif"), width = 2200, height = 2200, res = 300)
ggplot(spectra_var.importance, aes(x = wl)) + 
  geom_col(aes(y = import, fill = import)) +
  scale_fill_gradient2(low = "blue", mid = "violet", high = "red",midpoint = mean(spectra_var.importance$import), name="Importance (%)") +
  geom_line(aes(y = spectrum*150), size = 1.5, color="black", group = 1) +
  xlab("Wavelength (nm)") + ylab("Variable Importane (%) from the Unsupervised RF")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = 70,color = "darkred",linetype = "dashed",size=0.7) +
  theme_classic2(base_size = 16)

dev.off()
########################################################################################





### MODULE 3: Correlation analysis

# Import results
Observation_outliers <- read.csv("Soil_data_spectraResampled_clustered_outliers.csv", header = TRUE, sep = ";", dec = ".", row.names = "ID", na.strings = "NA",check.names=FALSE)

# Selecting the soil attribute to test the normality assumption
#soilAttribute <- 'Sand'
soilAttribute <- 'Clay'
#soilAttribute <- 'OM'
#soilAttribute <- 'CTC.ph7'
#soilAttribute <- 'V_perc'

# Select the laboratories to test the normality assumption (restricted to 5000 samples)
from_Lab <- 1
to_Lab <- from_Lab+1

# Normality test (The data is normal if the p-value is above 0.05)
Shapiro.test <- Observation_outliers %>% as.data.frame() %>% 
  dplyr::filter(!(is.na(soilAttribute))) %>% filter(Var_strat %in% c(from_Lab:to_Lab)) %>% group_by(Var_strat) %>% 
  select(all_of(soilAttribute)) %>%
  summarise_all(.funs = dplyr::funs(statistic = shapiro.test(.)$statistic, 
                                    p.value = shapiro.test(.)$p.value)) %>% ungroup()
Shapiro.test


# Correlation analysis

# Selecting variables for correlation (Stratification variable + Attributes + Spectra)
res.corr.sub <- Observation_outliers %>% select((c(Var_strat,paste(soilAttribute):paste(wavelength_max)))) %>% as.data.frame()

# Calculating correlation coefficients by Laboratory
library(corrr)
library(purrr)

# Select the method to be applied: kendall, spearman or pearson
cor.method <- "pearson"

Observation_lab_clustered <- list()

for(i in sort(unique(res.corr.sub$Var_strat))) {
  
  dataset.model <- res.corr.sub %>%
    filter(Var_strat %in% i) %>% correlate(method = cor.method) %>% mutate_at(vars(Var_strat), ~replace_na(.,i))
  Observation_lab_clustered[[i]] <- dataset.model
}

# Combine correlation matrix of all laboratories into a single matrix
res.corr.byLab.matrix <- do.call(rbind, Observation_lab_clustered)

# Removing values from the diagonal (=1 or =NA) of the correlation matrix
res.corr.byLab.matrix <- res.corr.byLab.matrix %>% filter(term %in% c("Clay","Sand","OM","CEC.ph7","V"))
res.corr.byLab.matrix <- res.corr.byLab.matrix %>% select(-c("Clay","Sand","OM","CEC.ph7","V"))

# Set rownames as "Var_stract-soilAttribute"
rownames(res.corr.byLab.matrix) <- paste0(res.corr.byLab.matrix$Var_strat,"_",res.corr.byLab.matrix$term)
res.corr.byLab.matrix <- res.corr.byLab.matrix %>% as.data.frame()

# Superheat plot -> https://rlbarter.github.io/superheat/adjacent-plots.html#line
library(superheat)
library(RColorBrewer)
library(hrbrthemes)

# Selecting the soil attribute to plot
#soilAttribute <- 'Sand'
soilAttribute <- 'Clay'
#soilAttribute <- 'OM'
#soilAttribute <- 'CTC.ph7'
#soilAttribute <- 'V_perc'

# Selecting the values of corr and ordering by rows for plotting
SoilAttr_corr <- res.corr.byLab.matrix %>% filter(term %in% soilAttribute) %>% 
  select(wavelength_min:wavelength_max) %>%
  apply(.,2,rev) %>% as.matrix()

# Calculating the averaged correlation and ordering by rows for plotting in the right of the Superheat figure
mean_corre <- res.corr.byLab.matrix %>% filter(term %in% soilAttribute) %>%
  select(wavelength_min:wavelength_max) %>% 
  apply(., 1, FUN = median) %>% 
  as.matrix() %>% apply(.,2,rev)

print.noquote(paste0(soilAttribute," -->   ",
                     "Min corr value: ", round(min(mean_corre, na.rm = TRUE),2)," | ",
                     "Max corr value: ", round(max(mean_corre, na.rm = TRUE),2)))

# Average spectrum to plot on the top of the correlation plot
spectrum <- Observation_outliers %>% select(wavelength_min:wavelength_max) %>% summarise_all(.funs = c(mean="median"))
spectrum <- spectrum %>% rename_with(stringr::str_replace,pattern = "_mean", replacement = "",matches("_mean")) %>% as.matrix()

# Plotting correlation values and averaged spectrum
tiff(paste0("Correlation_byLab_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
plot.superheat<- superheat(# Principal plot (correlation values by wavelength)
  X = SoilAttr_corr,
  left.label.col = "white", left.label.size = 0.1, left.label.text.size = 2, left.label.text.alignment = "left",
  row.title = paste0(cor.method," correlation values by Laboratory"), row.title.size = 3,
  column.title = "Wavelength (nm)", column.title.size = 3,
  bottom.label = "variable", bottom.label.size = 0.1,bottom.label.text.angle = 90,bottom.label.text.size = 2, force.bottom.label = T,
  # Tittle of the plot
  title = paste0(cor.method," correlation between spectra and ",soilAttribute, " by Laboratory"), title.size = 3,
  # Palette
  heat.pal = brewer.pal(n = 10, name = "RdBu"), heat.pal.values = c(0,-1, 1), heat.lim = c(-1,1),
  # Legend
  legend.height = 0.10, legend.width = 2, legend.text.size = 10,
  # Top plot (averaged spectrum)
  yt = spectrum, 
  yt.plot.type = "line", yt.line.col = "black", yt.axis = T, yr.axis = T, 
  yt.num.ticks = 3, yt.axis.name = "Averaged Reflectance", yt.axis.size = 10,
  # Right plot (averaged correlation)
  yr = mean_corre,
  yr.plot.type = "scatter", yr.axis.name = "Averaged r", yr.axis.name.size = 12)

dev.off() # Clear plot

gc()            # clear the memory
rm(list=setdiff(ls(), c("wavelength_min","wavelength_max")))  # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console





### MODULE 4: Plotting boxplot and spectra by cluster and Laboratory

# Import results
Observation_outliers <- read.csv("Soil_data_spectraResampled_clustered_outliers.csv", header = TRUE, sep = ";", dec = ".", row.names = "ID", na.strings = "NA",check.names=FALSE)

# Select a stratification variable code for filtering and graphical analysis (boxplot and spectral clusters)
var_strat_code <- 1

# Filter the data using a stratification variable code
spectra_to_plot <- Observation_outliers %>% filter(Var_strat %in% var_strat_code)

# Selecting the soil attribute for analysis
#soilAttribute <- 'Sand'
soilAttribute <- 'Clay'
#soilAttribute <- 'OM'
#soilAttribute <- 'CTC.ph7'
#soilAttribute <- 'V_perc'

# Defining the theme of the plot
themeBasic <- function () { 
  theme_classic2(base_size = 10) %+replace% 
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) + theme_update(plot.title = element_text(hjust = 0.5))
}

# Plot
library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)

plot <- spectra_to_plot %>% rename(Attr =paste0(soilAttribute)) %>% as.data.frame()

# Normality test by cluster (The data is normal if the p-value is above 0.05)
plot  %>% group_by(cluster) %>% select(Attr) %>% 
  summarise_all(.funs = dplyr::funs(statistic = shapiro.test(.)$statistic, 
                                    p.value = shapiro.test(.)$p.value)) %>% ungroup()

# Histogram plot of the selected soil attribute per cluster
tiff(paste0("Frequency_byCluster_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
ggplot(plot, aes(x = Attr)) + 
  geom_histogram(binwidth = 30) + geom_density(aes(y=30 * ..count..), col="red") +
  xlab(bquote(Clay~(g~kg^-1))) +
  #  xlab(paste0(soilAttribute))) +
  ylab(paste0("Frequency from the Lab ",Lab," by cluster")) +
  #  ggtitle(paste0("Histogram from the Lab ",Lab," by cluster")) +
  theme_classic2(base_size = 16) + theme(legend.position="none") + facet_wrap(~cluster)

# Plotting an adjusted boxplot for skewed distributions
library(robustbase)
d <- plot %>% group_by(cluster) %>% 
  summarise(ymin = adjboxStats(Attr, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[1],
            ymax = adjboxStats(Attr, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[5],
            middle = adjboxStats(Attr, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[3],
            lower = adjboxStats(Attr, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[2], 
            upper = adjboxStats(Attr, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)$stats[4],
            n = n(),
            outliers = list(Attr[which(Attr < ymin | Attr > ymax)])) %>% ungroup()


tiff(paste0("Boxplot_byCluster_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
d %>% mutate(clust.ordered = fct_reorder(factor(cluster), middle)) %>%
  ggplot(aes(x= clust.ordered, fill=as.character(clust.ordered))) +
  geom_boxplot(aes(group=clust.ordered,
                   ymin=ymin, ymax=ymax, middle=middle, upper=upper, lower=lower), stat="identity") +
  geom_point(data = . %>% filter(sapply(clust.ordered, length) > 0) %>%
               select(clust.ordered, outliers) %>%
               unnest(cols = c(outliers)), aes(y = unlist(outliers),
                                               alpha=0.6, colour="red", shape = 19,size=3)) +
  scale_shape_identity() + scale_colour_manual(values = c("red")) +
  scale_fill_viridis(discrete = T, alpha=0.6) +
  scale_y_continuous(position = "left") +
  xlab(paste("Spectral Cluster from the Laboratory N° ",Lab)) + 
  ylab(bquote(Clay~(g~kg^-1))) +
  ggtitle(paste("Boxplot for skewed distributions")) +
  theme_classic2(base_size = 17) +
  theme(legend.position="none")


# Creating a hyperSpec Object from a Spectral Matrix and Wavelength Vector
library(hyperSpec)

plot.spectra <- plot %>% rename(Attr.is.outlier =paste0(soilAttribute,".is.outlier")) %>% as.data.frame()

plot.spectra.outlier.false <- plot.spectra %>% filter(Attr.is.outlier=="FALSE")
spc.outlier.false <- new("hyperSpec", spc = select(plot.spectra.outlier.false,c(wavelength_min:wavelength_max)), 
                         wavelength = as.numeric(colnames(select(plot.spectra.outlier.false,c(wavelength_min:wavelength_max)))), 
                         data = select(plot.spectra.outlier.false,c("Sample","cluster","Attr","Attr.is.outlier")))

plot.spectra.outlier.true <- plot.spectra %>% filter(Attr.is.outlier=="TRUE")
spc.outlier.true <- new("hyperSpec", spc = select(plot.spectra.outlier.true,c(wavelength_min:wavelength_max)), 
                        wavelength = as.numeric(colnames(select(plot.spectra.outlier.true,c(wavelength_min:wavelength_max)))), 
                        data = select(plot.spectra.outlier.true,c("Sample","cluster","Attr","Attr.is.outlier")))

# Set axis labels
labels(spc.outlier.false, ".wavelength") <- "Wavelength (nm)"
labels(spc.outlier.false, "spc") <- "Reflectance Factor"

labels(spc.outlier.true, ".wavelength") <- "Wavelength (nm)"
labels(spc.outlier.true, "spc") <- "Reflectance Factor"


# Plotting spectra by cluster 
# samples with soil attributes flagged as outliers have their spectra in red color
tiff(paste0("Spectra_byCluster_",soilAttribute,".tif"), width = 2200, height = 1500, res = 300)
ggplot() +
  geom_line(data = as.long.df(spc.outlier.false), aes(x = .wavelength, y= spc, group= Sample, color="FALSE")) +
  geom_line(data = as.long.df(spc.outlier.true), aes(x = .wavelength, y= spc, group= Sample, color="TRUE"),size=0.5) +
  scale_color_manual(name = "Colors", 
                     values = c("FALSE" = "grey", "TRUE" = "red")) +
  scale_x_continuous(breaks=seq(350, 2500, 750)) +
  xlab(paste0("Wavelength (nm)")) +
  ylab(paste0("Reflectance Spectra from Laboratory N° ",Lab," by cluster")) +
  #    ggtitle(paste0("Spectra by cluster from the Lab ",Lab)) +
  guides(col=guide_legend(paste0(soilAttribute," is outlier"))) +
  theme_classic2(base_size = 12) +
  theme(legend.position="right") + facet_wrap(~cluster, nrow = 2)

gc()            # clear the memory
rm(list=ls())   # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console