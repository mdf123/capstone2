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



#View(train_data)
ncol(train_data)
test_data <- test_data %>% select(-contains("call"))
ncol(test_data)
#View(test_data)

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
#View(all_data)
#str(all_data)
#View(labels)

#Sort patient by patient numbers
all_data_desc <- all_data[, 1:2]
all_data_data <- all_data[, 3:ncol(all_data)]
all_data_data <- all_data_data[, order(as.numeric(names(all_data_data)))]

all_data <- cbind(all_data_desc, all_data_data)
View(all_data)

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
top_DE_genes <- topTable(fit, n=70, adjust="fdr")

str(all_data)
top_data <- all_data %>% filter(`Gene Accession Number` %in% top_DE_genes$Gene.Accession.Number)

View(top_data)

#numerical matrix
x_m <- top_data[, 3:ncol(top_data)]
#View(x_m)


#Transform to make tidy: observations = patients in rows, genes in columns
x_m <- t(x_m)
ncol(x_m) #70 top DE genes
nrow(x_m) ##72 patients

colnames(x_m) <- top_data$`Gene Accession Number`
View(x_m)

#Normalize
x_m_norm <- scale(x_m)
#x_m_norm <- sweep(x_m, 2, colMeans(x_m))
#x_m_norm <- sweep(x_m_norm, 2, colSds(x_m_norm), "/")


#Create Train- and Testdata
indexTrain <- createDataPartition(labels$cancer, p=0.66, list = F)
norm_train_data <- x_m_norm[indexTrain, ]
norm_test_data <- x_m_norm[-indexTrain, ]
train_label <- labels[indexTrain, ]
test_label <- labels[-indexTrain, ]
str(norm_train_data)
nrow(norm_train_data)
nrow(train_label)

#SVN
svmlinear1 <- train(x = norm_train_data, y= train_label$cancer, method = "svmLinear")
svmlinear1

`