rm(list=ls())
library(tidyverse)
library(caret)
library(matrixStats)
library(janitor)
set.seed(1, sample.kind="Rounding")

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

#Transform to make tidy
all_data <- t(all_data)


#colnames
all_data <- all_data %>% row_to_names(row_number = 2)


#View(all_data)
str(all_data)
class(all_data)
all_data <- apply(all_data, 2, as.numeric)
str(all_data)

#rownames
rownames(all_data) <- paste("patient", 1:nrow(all_data), sep=" ")
#View(all_data)


#Normalize
all_data_norm <- sweep(all_data, 2, colMeans(all_data))
all_data_norm <- sweep(all_data_norm, 2, colSds(all_data_norm), "/")
colSds(all_data_norm)
colMeans(all_data_norm)
#View(all_data_norm)

#Split data into train and test data

indexTrain <- createDataPartition(labels$cancer, p=0.66, list = F)
norm_train_data <- all_data_norm[indexTrain, ]
norm_test_data <- all_data_norm[-indexTrain, ]
train_label <- labels[indexTrain, ]
test_label <- labels[-indexTrain, ]
str(norm_train_data)
nrow(norm_train_data)
nrow(train_label)





#Support Vector machines
#https://rpubs.com/uky994/593668

svmlinear1 <- train(x = norm_train_data, y= train_label$cancer, method = "svmLinear")
svmlinear1

results <- tibble(model="Linear 1", accuracy=svmlinear1$results[which.max(svmlinear1$results[,2]), 2])
results

#Tune cost factor c
svmlinear2 <- train(x = norm_train_data, y= train_label$cancer, method = "svmLinear", tuneGrid = expand.grid(C=seq(0.1,1,0.2)))
svmlinear2
svmlinear2_result <- c(model="Linear 2", accuracy=svmlinear2$results[which.max(svmlinear2$results[,2]), 2])
results <- rbind(results, svmlinear2_result)
results

#Radial
modelLookup("svmRadial")
svmRadial1 <- train(x = norm_train_data, y= train_label$cancer, method = "svmRadial", tuneLength = 10)
svmRadial1$results
svmRadial1_result <- c(model="Radial 1", accuracy=svmRadial1$results[which.max(svmRadial1$results[,3]), 3])

results <- rbind(results, svmRadial1_result)
results

#Polynomial
modelLookup("svmPoly")
svmPoly1 <- train(x = norm_train_data, y= train_label$cancer, method = "svmPoly", tuneLength = 4)
svmPoly1$results
svmPoly1_result <- c(model="Poly 1", accuracy=svmPoly1$results[which.max(svmPoly1$results[,4]), 4])
results <- rbind(results, svmPoly1_result)
results

svmPoly1$bestTune

#--------- Decision tree ------------------------------

train_rpart <- train(norm_train_data, factor(train_label$cancer),
                     method = "rpart",
                     tuneGrid = data.frame(cp = seq(0.0, 0.1, len = 25)))

train_rpart
plot(train_rpart)


#---------------- PCA ----------------------------------
#Perform pca on train data
pca <- prcomp(norm_train_data)
str(pca)
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

#number of components for 95% variability
cs <- cumsum(pca.var.per)
nc <- length(cs[cs < 95]) + 1
nc

#Select the components
pca_train_data <- pca$x[, 1:nc]
dim(pca_train_data)

train_data_pca <- data.frame(y=train_label$cancer, pca_train_data)
View(train_data_pca)

train_rpart <- train(y ~ .,
                     method = "rpart",
                     tuneGrid = data.frame(cp = seq(0.0, 0.1, len = 25)),
                     data = train_data_pca)
plot(train_rpart)

train_rf_2 <- train(y ~ .,
                    method = "Rborist",
                    tuneGrid = data.frame(predFixed = 2, minNode=seq(1,10, by=1)),
                    data = train_data_pca)
plot(train_rf_2)
