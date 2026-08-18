setwd("path/to/rt_qpcr/")
library("ggplot2")
library("dplyr")
library("plotly")
library("viridis")
library("hrbrthemes")
library("foreign")
library("pheatmap")
library("randomForest")
library("doParallel")
library("caret")
library("e1071")
library("glmnet")
library("doParallel")

#data <- read.table("Machine Learning Input.txt", header = TRUE)
#gr <- c(rep("UCR", 10), rep("UCN", 10))
#write.arff(data, file = "Machine Learning Analysis.arff")

data2 <- read.table("Machine Learning Input - 2.txt", header = TRUE)
gr <- c(rep("UCR", 8), rep("UCN", 10))
write.arff(data, file = "Machine Learning Analysis 2.arff")


#NaiveBayes
library(e1071)
model <- naiveBayes(gr ~ ., data = data)
predict(model, data)
predict(model, data, type = "raw")
pred <- predict(model, data)
table(pred, data$gr)
confusionMatrix(factor(pred),factor(gr), positive = "UCN")

#Lasso Feature Selection
data <- read.table("Machine Learning Input.txt", header = TRUE)
data <- data[,-6]
data <- as.matrix(data)
gr <- c(rep("UCR", 10), rep("UCN", 10))
gr <- as.factor(gr)

cl <- makeCluster(detectCores(), type="PSOCK")
registerDoParallel(cl)

set.seed(100478)
cv.lasso_cont <- cv.glmnet(data, gr, family = "binomial", nfolds = 10 ,
                           alpha = 1,  standardize = TRUE, type.measure = "mse")
Coef_cont <- coef(cv.lasso_cont, s = cv.lasso_cont$lambda.min)
registerDoSEQ()

Coef_cont <- as.matrix(Coef_cont)
Coef1_cont <- data.frame("Gene" = row.names(Coef_cont),"coef" = Coef_cont[,1],
                         "abs" = abs(Coef_cont[,1]))
Coef1_cont <- Coef1_cont[Coef1_cont$coef !=0,]
Coef1_cont <- Coef1_cont[-1,]
write.table(Coef1_cont, "FeaturesLassoValue.txt", sep = "\t")







