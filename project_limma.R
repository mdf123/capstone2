rm(list=ls())
library(tidyverse)
library(caret)
library(matrixStats)
library(limma)
library(Biobase)

#Data: https://www.kaggle.com/datasets/crawford/gene-expression
train_data <- read_csv("data_set_ALL_AML_train.csv", col_names = T)
#View(train_data)
nrow(train_data)

#View(train_data)
test_data <- read_csv("data_set_ALL_AML_independent.csv")
labels <- read_csv("actual.csv")
#View(labels)


#Remove call columns
train_data <- train_data %>% select(-contains("call"))
#View(train_data)

#In train data the columns are not sorted by the patient number
train_data_desc <- train_data[, 1:2]
train_data_data <- train_data[, 3:ncol(train_data)]
train_data_data <- train_data_data[, order(as.numeric(names(train_data_data)))]
#View(train_data_data)

train_data <- cbind(train_data_desc, train_data_data)

#View(train_data)
ncol(train_data)
test_data <- test_data %>% select(-contains("call"))
ncol(test_data)

#Check for missing values
nas <- function(x) ( any(is.na(x)) )
str(train_data)
str(test_data)
print(train_data %>% summarise_all(nas), width = Inf)
print(test_data %>% summarise_all(nas), width = Inf)

#Merging all data
all_data <- train_data %>% left_join(test_data)
ncol(all_data)
nrow(all_data)
View(all_data)
str(all_data)
#View(labels)

#features
f <- all_data[, 1:2]
head(f)
#expression matrix
x <- as.matrix(all_data[, 3:ncol(all_data)])

head(x)
head(labels)

#Colnames for expression set necessary
colnames(x) <- 1:ncol(x)
head(x)

eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(labels),
                      featureData = AnnotatedDataFrame(f))


#Design Matrix
design <- model.matrix(~cancer, data=pData(eset))
head(design)
#Should be the number of observations
nrow(design)
colSums(design)
sum(labels$cancer=="AML")

#fit the model
fit <- lmFit(eset, design)

#Calculate the t-statistics
fit <- eBayes(fit)

results <- decideTests(fit[, "cancerAML"])
summary(results)

#List of differentially expressed genes
topTable(fit, n=70, adjust="fdr")
