rm(list=ls())
library(tidyverse)
library(caret)
library(matrixStats)
library(janitor)

#Data: https://www.kaggle.com/datasets/crawford/gene-expression
train_data <- read_csv("data_set_ALL_AML_train.csv", col_names = T)
View(train_data)
nrow(train_data)

#View(train_data)
test_data <- read_csv("data_set_ALL_AML_independent.csv")
labels <- read_csv("actual.csv")
View(labels)


#Remove call columns
train_data <- train_data %>% select(-contains("call"))
View(train_data)

#In train data the columns are not sorted by the patient number
train_data_desc <- train_data[, 1:2]
train_data_data <- train_data[, 3:ncol(train_data)]
train_data_data <- train_data_data[, order(as.numeric(names(train_data_data)))]
View(train_data_data)

train_data <- cbind(train_data_desc, train_data_data)

View(train_data)
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

#Transform to make tidy
all_data <- t(all_data)


all_data <- all_data %>% row_to_names(row_number = 2)
rownames(all_data) <- paste("patient", 1:nrow(all_data), sep=" ")
View(all_data)

#Numeric matrix

all_data_m <- matrix(as.numeric(unlist(all_data)),nrow=nrow(all_data))
View(all_data_m)
str(all_data_m)
geneMeans <- colMeans(all_data_m)
geneMeans

as.data.frame(geneMeans) %>% ggplot(aes(x=geneMeans)) + geom_histogram()

geneSds <- colSds(all_data_m)
as.data.frame(geneSds) %>% ggplot(aes(x=geneSds)) + geom_histogram() + scale_x_log10()
min(geneSds)

all_data_norm <- sweep(all_data, 2, colMeans(all_data))
all_data_norm <- sweep(all_data_norm, 2, colSds(all_data_m), "/")
colSds(all_data_norm)
colMeans(all_data_norm)
View(all_data_norm)

pca <- prcomp(all_data_norm)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca$x[,1]
pca$x[,2]
rownames(pca$x)
