setwd("path/to/data/")

library("randomForest")
library("doParallel")
library("caret")
library("e1071")
library("glmnet")
library("doParallel")

#Lasso Feature Selection
DEGs_weka_Expressions <- read.table("DEGs_weka_Expressions.txt")
#DEGs_weka_Expressions <- as.data.frame(scale(DEGs_weka_Expressions))
#data <- data[,-1]
#write.table(data, file = "In Silico Input for Machine Learning.txt")
ran_forest <- t(as.matrix(DEGs_weka_Expressions))
gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))
gr <- as.factor(gr)

cl <- makeCluster(detectCores(), type="PSOCK")
registerDoParallel(cl)

set.seed(100478)
cv.lasso_cont <- cv.glmnet(ran_forest, gr, family = "binomial", nfolds = 10 ,
                   alpha = 1,  standardize = TRUE, type.measure = "mse")
Coef_cont <- coef(cv.lasso_cont, s = cv.lasso_cont$lambda.min)
registerDoSEQ()

Coef_cont <- as.matrix(Coef_cont)
Coef1_cont <- data.frame("Gene" = row.names(Coef_cont),"coef" = Coef_cont[,1],
                         "abs" = abs(Coef_cont[,1]))
Coef1_cont <- Coef1_cont[Coef1_cont$coef !=0,]
Coef1_cont <- Coef1_cont[-1,]
write.table(Coef1_cont, "LassoFeatures.txt", sep = "\t")

###Lasso Features without seting a seed
DEGs_weka_Expressions <- read.table("DEGs_weka_Expressions.txt")
ran_forest <- t(as.matrix(DEGs_weka_Expressions))
gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))
gr <- as.factor(gr)

gene_coef <- c(rep(0, 31))
gene_counter <- c(rep(0, 31))

for (i in seq(1:50)) {
  cl <- makeCluster(detectCores(), type="PSOCK")
  registerDoParallel(cl)
  
  cv.lasso_cont <- cv.glmnet(ran_forest, gr, family = "binomial", nfolds = 10 ,
                             alpha = 1,  standardize = TRUE, type.measure = "mse")
  Coef_cont <- coef(cv.lasso_cont, s = cv.lasso_cont$lambda.min)
  registerDoSEQ()
  
  Coef_cont <- as.matrix(Coef_cont)
  Coef1_cont <- data.frame("Gene" = row.names(Coef_cont),"coef" = Coef_cont[,1],
                           "abs" = abs(Coef_cont[,1]))
  Coef1_cont <- Coef1_cont[Coef1_cont$abs >= 1,]
  Coef1_cont <- Coef1_cont[-1,]
  for (j in seq(1:dim(Coef1_cont)[1])) {
    ndx = which(row.names(DEGs_weka_Expressions)==Coef1_cont[j,1])
    gene_counter[ndx] = gene_counter[ndx] + 1
    gene_coef[ndx] = gene_coef[ndx] + Coef1_cont[j,3]
  }
}

gene_coef = gene_coef/gene_counter
data <- data.frame(row.names(DEGs_weka_Expressions), gene_counter, gene_coef)
write.table(data, "RobustLassoFeatures.txt", sep = "\t")

####Brota Feature Selection
#mli <- read.table("mli.txt", header = TRUE)
#row.names(mli) <- mli$Accession.Code
#mli <- mli[,c(-1,-5)]
#data <- read.table("DEGs_weka_Expressions.txt", header = TRUE)  
#data <- t(data)
#data <- data[order(row.names(data)),]
#mli <- mli[order(row.names(mli)),]
#row.names(mli)[5] == row.names(data)[5]
#row.names(mli)[34] == row.names(data)[34]
#i = 3
#library(Boruta)
#set.seed(5783125)
#y_train = mli[mli$fold != i,]
#y_test = mli[mli$fold == i,]
#data_train <- data[row.names(y_train),]
#data_test <- data[row.names(y_test),] 
#
#y_train$Group <- factor(y_train$Group,levels = c("UCR","UCN"))
#y_train$Group <- relevel(y_train$Group, ref = "UCR")
#d_train <- data.frame(data_train, "Group" = y_train$Group)
#
#boruta_output <- Boruta(factor(Group) ~ ., data=na.omit(d_train), doTrace=0)
#
#boruta_output <- Boruta(x= data_train,y = as.factor(y_train$Group),
#                          mcAdj = TRUE,maxRuns = 1000)
#names(boruta_output)
#boruta_signif <- getSelectedAttributes(boruta_output, withTentative = FALSE)
#print(boruta_signif)
#
#data_train <- data_train[,boruta_signif]
#data_test <- data_test[,boruta_signif]
#
#library(doParallel)
#library(randomForest)
#cl <- makeCluster(detectCores(), type='PSOCK')
#registerDoParallel(cl)
#RF <- randomForest(x = data_train, y = y_train$Group, ntree=10000)
#registerDoSEQ()
#
#pred_RF <- predict(RF,newdata = data_test, type = "class")
#confusionMatrix(factor(pred_RF),factor(y_test$Group), positive = "UCN")
#
#IMP <- as.data.frame(RF$importance)
#IMP$gene <- row.names(IMP)
#IMP <- IMP[order(IMP$MeanDecreaseGini, decreasing = TRUE),]
#
#gene1 <- IMP$gene[1:4]
#
#data_train <- data_train[,gene1]
#data_test <- data_test[,gene1]
#
#cl <- makeCluster(detectCores(), type='PSOCK')
#registerDoParallel(cl)
#RF <- randomForest(x = data_train, y = y_train$Group, ntree=10000)
#registerDoSEQ()
#
#pred_RF <- predict(RF,newdata = data_test, type = "class")
#confusionMatrix(factor(pred_RF),factor(y_test$Group), positive = "UCN")




