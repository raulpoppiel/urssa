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

# Selecting the soil attribute for analysis
soilAttribute <- 'Clay'
#soilAttribute <- 'Sand'
#soilAttribute <- 'OM'
#soilAttribute <- 'CTC.ph7'
#soilAttribute <- 'V_perc'

# Import results and rename the variable for study
Observation_outliers <- read.csv("Soil_data_spectraOriginal_clustered_outliers.csv", header = TRUE, sep = ";", dec = ".", na.strings = "NA")
Observation <- Observation_outliers %>% rename(Attr = paste0(soilAttribute)) %>% as_tibble()

# Removing NAs and outliers (detected in previous script)
Observation <- Observation %>% filter(!is.na(Attr))
Observation <- Observation %>% rename(Attr.is.outlier =paste0(soilAttribute,".is.outlier"))
Observation <- Observation %>% filter(Attr.is.outlier %in% "FALSE")

#######################################################
# Partitioning of data into training and test subsets #

# Group dataset by Laboratory and Cluster and 
# order by the mean spectral reflectance from the most important wavelength (var importance > 70%)
var.importance <- read.csv("spectra.ruf.importance.result.csv", header = TRUE, sep = ";", dec = ".", na.strings = "NA",check.names=FALSE)
importance <- 70 # select the percentage of wavelength importance
wl.most.important <- var.importance %>% filter(percent >= importance) %>% select(variables) %>% arrange(variables)

Observation_ordered <- Observation %>% mutate(mean=rowMeans(select(Observation,c(all_of(paste0("X",wl.most.important$variables)))))) %>%
  group_by(Var_strat,cluster) %>% arrange(mean, .by_group = TRUE) %>% ungroup() #%>% select(-c(mean))

# Obtaining the ID that will be used for modelling by Var_strat (Laboratory) and by Cluster
Observation_ordered <- Observation_ordered %>% mutate(ID_model = paste0(Var_strat,"_",cluster))

# Splitting data into test and train sets by selecting the k_th row/line (sample) of the sequence
k_th <- 3 # Select the amount of data for testing (2=50%, 3=33%, 4=25%, 5=20%, 6=16%...)
k_th_row<-round(c(seq(1,nrow(Observation_ordered),k_th)),0)
train.set <- Observation_ordered %>% slice(-c(k_th_row))
test.set <- Observation_ordered %>% slice(k_th_row)

Perc.to.test<- round(nrow(test.set)/nrow(Observation_ordered),2)*100
print(paste0("Test set size = ",Perc.to.test," % of samples"))
Perc.to.train<- 100-Perc.to.test
print(paste0("Train set size = ",Perc.to.train," % of samples"))

# Use the ID column to verify if the total number of samples is equal to the initial dataset (before splitting)
all.equal(sort(c(train.set$ID,test.set$ID)),sort(c(Observation$ID)))

##################
# Plotting subsets

#Plotting the frequency by Laboratory (Var_strat)
library(ggpubr)
tiff(paste0("Frequency_byLab_trainSet_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
ggplot(train.set, aes(x = Attr)) + 
  geom_histogram(binwidth = 30) + geom_density(aes(y=30 * ..count..), col="red", size=1.5) +
#  xlab(bquote(Clay~(g~kg^-1))) +
  xlab(paste0(soilAttribute," g kg-1")) +
#  ylim(0,30) +
  scale_x_continuous(breaks=seq(0, 1000, 330)) +
  ylab(paste0("Frequency")) +
  ggtitle(paste0("Histogram of the train set (" ,Perc.to.train, "%) by Var_strat (Lab)")) +
  theme_classic2(base_size = 13) + theme(legend.position="none") + facet_wrap(~Var_strat)

tiff(paste0("Frequency_byLab_testSet_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
ggplot(test.set, aes(x = Attr)) + 
  geom_histogram(binwidth = 30) + geom_density(aes(y=30 * ..count..), col="red", size=1.5) +
#  xlab(bquote(Clay~(g~kg^-1))) +
  xlab(paste0(soilAttribute," g kg-1")) +
#  ylim(0,30) +
  scale_x_continuous(breaks=seq(0, 1000, 330)) +
  ylab(paste0("Frequency")) +
  ggtitle(paste0("Histogram of the test set (" ,Perc.to.test, "%) by Var_strat (Lab)")) +
  theme_classic2(base_size = 13) + theme(legend.position="none") + facet_wrap(~Var_strat)

dev.off()

#Plotting the frequency from one Laboratory (Var_strat) and by cluster
# Select a laboratory to analyze their histogram of the train and test sets
var.strat <- 1

# Filter the dataset by the selected laboratory
attribute.plot.train <- train.set %>% filter(Var_strat %in% var.strat)
attribute.plot.test <- test.set %>% filter(Var_strat %in% var.strat)

# Defining the theme of the plot
themeBasic <- function () { 
  theme_classic2(base_size = 10) %+replace% 
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) + theme_update(plot.title = element_text(hjust = 0.5))}

# Plotting histogram by cluster
library(ggplot2)
library(ggpubr)

tiff(paste0("Frequency_byCluster_trainSet_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
ggplot(attribute.plot.train, aes(x = Attr)) + 
  geom_histogram(binwidth = 30) + geom_density(aes(y=30 * ..count..), col="red", size=1.5) +
#  xlab(bquote(Clay~(g~kg^-1))) +
  xlab(paste0(soilAttribute," g kg-1")) +
#  ylim(0,15) +
  scale_x_continuous(breaks=seq(0, 1000, 330)) +
  ylab(paste0("Frequency")) +
  ggtitle(paste0("Histogram of the train set (" ,Perc.to.train, "%) from the Lab ",var.strat," by cluster")) +
  theme_classic2(base_size = 16) + theme(legend.position="none") + facet_wrap(~cluster)

tiff(paste0("Frequency_byCluster_testSet_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".tif"), width = 2200, height = 2200, res = 300)
ggplot(attribute.plot.test, aes(x = Attr)) + 
  geom_histogram(binwidth = 30) + geom_density(aes(y=30 * ..count..), col="red", size=1.5) +
#  xlab(bquote(Clay~(g~kg^-1))) +
  xlab(paste0(soilAttribute," g kg-1")) +
#  ylim(0,15) +
  scale_x_continuous(breaks=seq(0, 1000, 330)) +
  ylab(paste0("Frequency")) +
  ggtitle(paste0("Histogram of the test set (" ,Perc.to.test, "%) from the Lab ",var.strat," by cluster")) +
  theme_classic2(base_size = 16) + theme(legend.position="none") + facet_wrap(~cluster)

dev.off()

####################
# CUBIST ALGORITHM #

# Prediction of soil attributes
# it can be made by: i) Var_strat or ii) ID_model (Var_strat + Cluster)
library(Cubist)
library(chillR)
library(ModelMetrics)

wavelength_min <- "350" # adjust the min spectral wavelength
wavelength_max <- "2500" # adjust the max spectral wavelength
predictors <- colnames(select(train.set,c(paste0("X",wavelength_min):paste0("X",wavelength_max)))) # select only preditors (covariates)

# Select the model using the "#"
# See our paper -> Literature suggests for: i) small (40-70) sample set, use all samples; ii) large (>500) sample sets, stratification/clustering would improve performance.

model <- "Var_strat" # e.g.: Laboratories
#model <- "ID_model" # e.g.: Laboratories + cluster intralaborory

models.val.list.train <- list()
models.val.list.test <- list()

for(i in unique(train.set %>% select(all_of(model)) %>% pull())) {
  
  #Predictors
  dataset.model.predictors <- train.set %>% filter(!!as.symbol(all_of(model)) %in% i) %>% select(all_of(predictors))
  dataset.validation.predictors <- test.set %>% filter(!!as.symbol(all_of(model)) %in% i) %>% select(all_of(predictors))
  
  # Variable to be predicted
  dataset.model.variable <- train.set %>% filter(!!as.symbol(all_of(model)) %in% i) %>% select(c(Attr))
  dataset.validation.variable <- test.set %>% filter(!!as.symbol(all_of(model)) %in% i) %>% select(c(Attr))
  
  # Complete Dataset
  dataset.train <- train.set %>% filter(!!as.symbol(all_of(model)) %in% i)
  dataset.test <- test.set %>% filter(!!as.symbol(all_of(model)) %in% i)
  
  # Training
  model.attr <- cubist(x = dataset.model.predictors,
                       y = dataset.model.variable$Attr,
                       committees = 100)
  
  # Evaluation of the training set
  predicted.values.attr.train <- predict(model.attr, dataset.model.predictors)
  #cor (predicted.values.attr.train, dataset.model.variable$Attr)^2 # R2
  #rmse(predicted.values.attr.train, dataset.model.variable$Attr)   # RMSE
  #RPIQ(predicted.values.attr.train, dataset.model.variable$Attr)   # rpiq 
  
  predicted.val.attr.train <- data.frame("ID" = dataset.train$ID,
                                         "model_used" = dataset.train %>% select(all_of(model)) %>% pull(),
                                         "Var_strat" = dataset.train$Var_strat,
                                         "cluster" = dataset.train$cluster,
                                         "ID_model" = dataset.train$ID_model,
                                         "obs" = dataset.train$Attr,
                                         "pred" = predicted.values.attr.train)
  
  models.val.train <- predicted.val.attr.train %>% dplyr::filter(!(is.na(obs)))
  
  eval.train <- data.frame("RMSE" = sqrt(mean((models.val.train$pred - models.val.train$obs)^2)),
                           "RPIQ" = IQR(models.val.train$obs)/sqrt(mean((models.val.train$pred - models.val.train$obs)^2)),
                           #"R2" = 1-((sum((models.val.train$obs - models.val.train$pred)^2))/(sum((models.val.train$obs - mean(models.val.train$obs))^2))),
                           "R2" = cor(models.val.train$obs,models.val.train$pred)^2,
                           models.val.train)
  
  
  
  # Evaluation of the testing set
  predicted.values.attr.test <- predict(model.attr, dataset.validation.predictors)
  #cor (predicted.values.attr.test, dataset.validation.variable$Attr)^2 # R2
  #rmse(predicted.values.attr.test, dataset.validation.variable$Attr)   # RMSE
  #RPIQ(predicted.values.attr.test, dataset.validation.variable$Attr)   # rpiq 
  
  predicted.val.attr.test <- data.frame("ID" = dataset.test$ID,
                                        "model_used" = dataset.test %>% select(all_of(model)) %>% pull(),
                                        "Var_strat" = dataset.test$Var_strat,
                                        "cluster" = dataset.test$cluster,
                                        "ID_model" = dataset.test$ID_model,
                                        "obs" = dataset.test$Attr,
                                        "pred" = predicted.values.attr.test)
  
  models.val.test <- predicted.val.attr.test %>% dplyr::filter(!(is.na(obs)))
  
  eval.test <- data.frame("RMSE" = sqrt(mean((models.val.test$pred - models.val.test$obs)^2)),
                          "RPIQ" = IQR(models.val.test$obs)/sqrt(mean((models.val.test$pred - models.val.test$obs)^2)),
                          #"R2" = 1-((sum((models.val.test$obs - models.val.test$pred)^2))/(sum((models.val.test$obs - mean(models.val.test$obs))^2))),
                          "R2" = cor(models.val.test$obs,models.val.test$pred)^2,
                          models.val.test)
  
  models.val.list.train[[i]] <- eval.train
  models.val.list.test[[i]] <- eval.test
}

# Combine all laboratories into a single matrix
models.val.results.train <- do.call(rbind, models.val.list.train)
models.val.results.test <- do.call(rbind, models.val.list.test)

# Export validation results
write.table(models.val.results.train,paste0("Cubist_Val_TrainSet_results_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".csv"), sep = ";", dec = ".", row.names = FALSE)
write.table(models.val.results.test,paste0("Cubist_Val_TestSet_results_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,".csv"), sep = ";", dec = ".", row.names = FALSE)

# Calculate and Export mean validation values by Var_strat (Lab)
models.val.results.train.mean <- models.val.results.train %>% group_by_(model) %>% summarise_at(vars(c("RMSE","RPIQ","R2")),funs(mean))
models.val.results.test.mean <- models.val.results.test %>% group_by_(model) %>% summarise_at(vars(c("RMSE","RPIQ","R2")),funs(mean))

write.table(models.val.results.train.mean,paste0("Cubist_Val_TrainSet_results_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,"_mean.csv"), sep = ";", dec = ".", row.names = FALSE)
write.table(models.val.results.test.mean,paste0("Cubist_Val_TestSet_results_",Perc.to.train,"-",Perc.to.test,"_",soilAttribute,"_mean.csv"), sep = ";", dec = ".", row.names = FALSE)

# Assessing results
var.strat=1
p.val <- ggplot(filter(models.val.results.test, Var_strat==var.strat), aes(x = obs, y = pred)) +
  geom_abline(color = "gray", size = 1, linetype = "longdash") + 
  geom_smooth(method=lm , size = 2, linetype = "solid", color="black", fullrange = TRUE, se=FALSE) + 
  geom_point(size = 3, alpha = 0.2) +
  ggtitle(paste0("Predicted vs. observed scatterplot for ",soilAttribute," from the testing set")) +
  xlim(0,1000) + ylim(0,1000) +
  coord_fixed() ; p.val

# Density plot of the distribution of soil attribute values
densityplot(x=filter(models.val.results.train, Var_strat==var.strat)$obs, xlab=paste0("Observed ",soilAttribute," g kg-1")) # Distribution of observed values
densityplot(x=filter(models.val.results.test, Var_strat==var.strat)$pred, xlab=paste0("Predicted ",soilAttribute," g kg-1")) # Distribution of predicted values

# Quantify number of samples used by model
train.set %>% filter(Var_strat==var.strat) %>% group_by_(model) %>% summarise(n()) # number of samples per cluster in the training set
test.set %>% filter(Var_strat==var.strat) %>% group_by_(model) %>% summarise(n())# number of samples per cluster in the testing set

gc()            # clear the memory
rm(list=ls())   # clear the Environment
graphics.off()  # Clear all the plots
shell("cls")    # Clear the Console