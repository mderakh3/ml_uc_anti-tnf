setwd("path/to/roc_curve_analysis/")
library("edgeR")
library("data.table")
library("dplyr")
library("DESeq2")
library("limma")
library("doParallel")
library("randomForest")
library("caret")
library("glmnet")
library("rpart")
library("e1071")# for svm
library("RANN") # for SMOTE
library("class") # for knn
library("MASS") # for lda

#qPCR Data
data <- read.table("qPCR Input for ML.txt", header = TRUE)
gr <- c(rep("UCR", 12), rep("UCN", 10))
colDATA <- read.table("qPCR COLDATA.txt")


#In Silico Data
data <- read.table("In Silico Input for ROC.txt", header = TRUE)
gr <- c(rep("UCN", 23), rep("UCR", 12))
colDATA <- read.table("In Silico COLDATA.txt")


##################ROC Curve and AUC Analysis
ROC_table <- data.frame("response"= "a","predictor"= 0,
                        "fold"= 0,"NO.features"= 0,
                        "method" = "a")
for (i in 1:3) {
  set.seed(5783125)
  y_train = colDATA[colDATA$fold != i,]
  y_test = colDATA[colDATA$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  
  data_train3 <- data_train[,1:3]
  data_train4 <- data_train[,1:4]
  data_train5 <- data_train[,1:5]
  
  data_test3 <- data_test[,1:3]
  data_test4 <- data_test[,1:4]
  data_test5 <- data_test[,1:5]
  
  cl <- makeCluster(detectCores(), type ='PSOCK')
  registerDoParallel(cl)

  #SVM
  
  svm.model3 <- svm(y_train$gr ~ ., data = data_train3, probability = TRUE)
  svm.model4 <- svm(y_train$gr ~ ., data = data_train4, probability = TRUE)
  svm.model5 <- svm(y_train$gr ~ ., data = data_train5, probability = TRUE)
  
  svm.pred3 <- predict(svm.model3, data_test3,probability = TRUE)
  svm.pred4 <- predict(svm.model4, data_test4,probability = TRUE)
  svm.pred5 <- predict(svm.model5, data_test5,probability = TRUE)
  
  svm_pr3  <- as.data.frame(attr(svm.pred3 ,"probabilities"))
  svm_pr4  <- as.data.frame(attr(svm.pred4 ,"probabilities"))
  svm_pr5  <- as.data.frame(attr(svm.pred5 ,"probabilities"))
  
  ROC3 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr3$UCN,
                     "fold"= rep(i,length(svm.pred3)), 
                     "NO.features"= rep(3,length(svm.pred3)),
                     "method" = rep("svm",length(svm.pred3)))
  ROC4 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr4$UCN,
                     "fold"= rep(i,length(svm.pred4)), 
                     "NO.features"= rep(4,length(svm.pred4)),
                     "method" = rep("svm",length(svm.pred4)))
  ROC5 <- data.frame("response"=as.character(y_test$gr),"predictor"= svm_pr5$UCN,
                     "fold"= rep(i,length(svm.pred5)), 
                     "NO.features"= rep(5,length(svm.pred5)),
                     "method" = rep("svm",length(svm.pred5)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
  # LDA
  
  LDA3 <- lda(y_train$gr ~., data = data_train3, probability = TRUE)
  LDA4 <- lda(y_train$gr ~., data = data_train4, probability = TRUE)
  LDA5 <- lda(y_train$gr ~., data = data_train5, probability = TRUE)

  pred3 <- predict(LDA3, as.data.frame(data_test3),probability = TRUE)
  pred4 <- predict(LDA4, as.data.frame(data_test4),probability = TRUE)
  pred5 <- predict(LDA5, as.data.frame(data_test5),probability = TRUE)
  
  LDA_pr3  <- as.data.frame(pred3$posterior)
  LDA_pr4  <- as.data.frame(pred4$posterior)
  LDA_pr5  <- as.data.frame(pred5$posterior)
  
  ROC3 <- data.frame("response"=as.character(y_test$gr),"predictor"= LDA_pr3$UCR,
                     "fold"= rep(i,length(LDA_pr3$UCN)), 
                     "NO.features"= rep(3,length(LDA_pr3$UCN)),
                     "method" = rep("LDA",length(LDA_pr3$UCN)))
  ROC4 <- data.frame("response"=as.character(y_test$gr),"predictor"= LDA_pr4$UCN,
                     "fold"= rep(i,length(LDA_pr4$UCN)), 
                     "NO.features"= rep(4,length(LDA_pr4$UCN)),
                     "method" = rep("LDA",length(LDA_pr4$UCN)))
  ROC5 <- data.frame("response"=as.character(y_test$gr),"predictor"= LDA_pr5$UCN,
                     "fold"= rep(i,length(LDA_pr5$UCN)), 
                     "NO.features"= rep(5,length(LDA_pr5$UCN)),
                     "method" = rep("LDA",length(LDA_pr5$UCN)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
  
  ##logit
  y_train$gr <- as.numeric(y_train$gr)
  y_test$gr <- as.numeric(y_test$gr)
  y_train$gr[y_train$gr == 2]<- 0
  y_test$gr[y_test$gr == 2]<- 0
  
  mylogit3 <- glm(y_train$gr ~ ., data = data_train3, family='binomial')
  mylogit4 <- glm(y_train$gr ~ ., data = data_train4, family='binomial')
  mylogit5 <- glm(y_train$gr ~ ., data = data_train5, family='binomial')
  
  pred3<- predict(mylogit3 ,newdata = data_test3,type = "response", probability = TRUE)
  pred4<- predict(mylogit4 ,newdata = data_test4,type = "response", probability = TRUE)
  pred5<- predict(mylogit5 ,newdata = data_test5,type = "response", probability = TRUE)
  
  ROC3 <- data.frame("response"=y_test$gr,"predictor"= pred3,
                     "fold"= rep(i,length(pred3)), 
                     "NO.features"= rep(3,length(pred3)),
                     "method" = rep("logit",length(pred3)))
  ROC4 <- data.frame("response"=y_test$gr,"predictor"= pred4,
                     "fold"= rep(i,length(pred4)), 
                     "NO.features"= rep(4,length(pred4)),
                     "method" = rep("logit",length(pred4)))
  ROC5 <- data.frame("response"=y_test$gr,"predictor"= pred5,
                     "fold"= rep(i,length(pred5)), 
                     "NO.features"= rep(5,length(pred5)),
                     "method" = rep("logit",length(pred5)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
  #Random forest

  data_train <- as.matrix(data_train)
  data_test <- as.matrix(data_test)
  
  data_train3 <- data_train[,1:3]
  data_train4 <- data_train[,1:4]
  data_train5 <- data_train[,1:5]
  
  data_test3 <- data_test[,1:3]
  data_test4 <- data_test[,1:4]
  data_test5 <- data_test[,1:5]
  
  RF3 <- randomForest(x = data_train3,y = y_train$gr,ntree=100)
  RF4 <- randomForest(x = data_train4,y = y_train$gr,ntree=100)
  RF5 <- randomForest(x = data_train5,y = y_train$gr,ntree=100)
  
  pred3 <- predict(RF3 ,newdata = data_test3 ,type="class")
  pred4 <- predict(RF4 ,newdata = data_test4 ,type="class")
  pred5 <- predict(RF5 ,newdata = data_test5 ,type="class")
  
  ROC3 <- data.frame("response"=y_test$gr,"predictor"= pred3,
                     "fold"= rep(i,length(pred3)), 
                     "NO.features"= rep(3,length(pred3)),
                     "method" = rep("RF",length(pred3)))
  ROC4 <- data.frame("response"=y_test$gr,"predictor"= pred4,
                     "fold"= rep(i,length(pred4)), 
                     "NO.features"= rep(4,length(pred4)),
                     "method" = rep("RF",length(pred4)))
  ROC5 <- data.frame("response"=y_test$gr,"predictor"= pred5,
                     "fold"= rep(i,length(pred5)), 
                     "NO.features"= rep(5,length(pred5)),
                     "method" = rep("RF",length(pred5)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
}

ROC_table <- ROC_table[-1,]
write.table(ROC_table,"ROC Data, dct.txt", sep = "\t", row.names = F)

ME <- unique(ROC_table$method)
NF <- as.numeric(unique(ROC_table$NO.features))
FO <- as.numeric(unique(ROC_table$fold))

ROC_svm <- ROC_table[ROC_table$method == "svm",]
ROC_lda <- ROC_table[ROC_table$method == "LDA",]
ROC_logit <- ROC_table[ROC_table$method == "logit",]
ROC_rf <- ROC_table[ROC_table$method == "RF",]

ROC_svm_3 <- ROC_svm[ROC_svm$NO.features == 3,]
ROC_svm_4 <- ROC_svm[ROC_svm$NO.features == 4,]
ROC_svm_5 <- ROC_svm[ROC_svm$NO.features == 5,]


library(pROC)
roc.svm3 <- roc(response =as.character(ROC_svm_3$response), 
                predictor = as.numeric(ROC_svm_3$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("UCR","UCN")))
roc.svm4 <- roc(response =as.character(ROC_svm_4$response), 
                predictor = as.numeric(ROC_svm_4$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("UCR","UCN")))
roc.svm5 <- roc(response =as.character(ROC_svm_5$response), 
                predictor = as.numeric(ROC_svm_5$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("UCR","UCN")))

LAB <- c(paste("top 3 features: AUC =",round(roc.svm3$auc,2)), 
         paste("top 4 features: AUC =",round(roc.svm4$auc,2)),
         paste("top 5 features: AUC =",round(roc.svm5$auc,2)))

roc.list <- list(roc.svm3,roc.svm4,roc.svm5)
g.list <- ggroc(roc.list,alpha = 0.5, size = 0.85,legacy.axes = TRUE,
                print.auc = TRUE)
g.list + scale_colour_manual(values = c("red","#DF8B00","blue"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")

#LDA
ROC_lda_3 <- ROC_lda[ROC_lda$NO.features == 3,]
ROC_lda_4 <- ROC_lda[ROC_lda$NO.features == 4,]
ROC_lda_5 <- ROC_lda[ROC_lda$NO.features == 5,]


library(pROC)
roc.lda3 <- roc(response =as.character(ROC_lda_3$response), 
                predictor = as.numeric(ROC_lda_3$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("UCR","UCN")))

roc.lda4 <- roc(response =as.character(ROC_lda_4$response), 
                predictor = as.numeric(ROC_lda_4$predictor),
                smooth = FALSE,legacy.axes = TRUE,
                levels = as.factor(c("UCR","UCN")))

roc.lda5 <- roc(response =as.character(ROC_lda_5$response), 
                predictor = as.numeric(ROC_lda_5$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("UCR","UCN")))

LAB <- c(paste("top 3 features: AUC =",round(roc.lda3$auc,2)), 
         paste("top 4 features: AUC =",round(roc.lda4$auc,2)),
         paste("top 5 features: AUC =",round(roc.lda5$auc,2)))

roc.list <- list(roc.lda3,roc.lda4,roc.lda5)
g.list <- ggroc(roc.list,alpha = 0.5, size = 0.85,legacy.axes = TRUE,
                print.auc = TRUE)
g.list + scale_colour_manual(values = c("red","#DF8B00","blue"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")


#logit
ROC_logit_3 <- ROC_logit[ROC_logit$NO.features == 3,]
ROC_logit_4 <- ROC_logit[ROC_logit$NO.features == 4,]
ROC_logit_5 <- ROC_logit[ROC_logit$NO.features == 5,]


library(pROC)
roc.logit3 <- roc(response =as.character(ROC_logit_3$response), 
                  predictor = as.numeric(ROC_logit_3$predictor),
                  smooth = FALSE,legacy.axes = TRUE, 
                  levels = as.factor(c("0","1")))
roc.logit4 <- roc(response =as.character(ROC_logit_4$response), 
                  predictor = as.numeric(ROC_logit_4$predictor),
                  smooth = FALSE,legacy.axes = TRUE,
                  levels = as.factor(c("0","1")))
roc.logit5 <- roc(response =as.character(ROC_logit_5$response), 
                  predictor = as.numeric(ROC_logit_5$predictor),
                  smooth = FALSE,legacy.axes = TRUE, 
                  levels = as.factor(c("0","1")),)

LAB <- c(paste("top 3 features: AUC =",round(roc.logit3$auc,2)), 
         paste("top 4 features: AUC =",round(roc.logit4$auc,2)),
         paste("top 5 features: AUC =",round(roc.logit5$auc,2)))

roc.list <- list(roc.logit3,roc.logit4,roc.logit5)
g.list <- ggroc(roc.list,alpha = 0.5, size = 0.85,legacy.axes = TRUE,
                print.auc = TRUE)
g.list + scale_colour_manual(values = c("red","#DF8B00","blue"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")

#RF
ROC_rf_3 <- ROC_rf[ROC_rf$NO.features == 3,]
ROC_rf_4 <- ROC_rf[ROC_rf$NO.features == 4,]
ROC_rf_5 <- ROC_rf[ROC_rf$NO.features == 5,]


library(pROC)
roc.rf3 <- roc(response =as.character(ROC_rf_3$response), 
                predictor = as.numeric(ROC_rf_3$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("0","1")))

roc.rf4 <- roc(response =as.character(ROC_rf_4$response), 
                predictor = as.numeric(ROC_rf_4$predictor),
                smooth = FALSE,legacy.axes = TRUE,
                levels = as.factor(c("0","1")))

roc.rf5 <- roc(response =as.character(ROC_rf_5$response), 
                predictor = as.numeric(ROC_rf_5$predictor),
                smooth = FALSE,legacy.axes = TRUE, 
                levels = as.factor(c("0","1")))

LAB <- c(paste("top 3 features: AUC =",round(roc.rf3$auc,2)), 
         paste("top 4 features: AUC =",round(roc.rf4$auc,2)),
         paste("top 5 features: AUC =",round(roc.rf5$auc,2)))

roc.list <- list(roc.rf3,roc.rf4,roc.rf5)
g.list <- ggroc(roc.list,alpha = 0.5, size = 0.85,legacy.axes = TRUE,
                print.auc = TRUE)
g.list + scale_colour_manual(values = c("red","#DF8B00","blue"),
                             labels= LAB)+theme_bw()+
  theme(legend.text=element_text(size=19,family = "sans"),
        legend.title = element_text(size=0,family = "sans"),
        axis.title =element_text(size=26,family = "sans"),
        axis.text=element_text(size=12))+
  theme(legend.position = c(0.7, 0.4),
        legend.direction = "vertical")


###ROC Curve by ROCR
#library("ROCR")
#data <- read.table("dct.txt", header = TRUE)
#y_train <- data
#data_train <- as.data.frame(data[,-2])
#row.names(data_train) <- row.names(data)
#colnames(data_train) <- "H19"
#gr <- c(rep("Con", 7), rep("UC", 10))
#LDA <- lda(y_train$gr ~., data = data_train, probability = TRUE)
#pred <- predict(LDA, data_train, type = "probablity")
#pred <- prediction(pred, y_train$gr)
#eval <- performance(pred, "acc")
#plot(eval)
#max <- which.max(slot(eval, "y.values")[[1]])
#acc <- slot(eval, "y.values")[[1]][max]
#cut <- slot(eval, "x.values")[[1]][max]
#print(c(Accuracy = acc, Cutoff = cut))
#pred <- prediction(pred, data$gr)
#roc <- performance(pred, "tpr", "fpr")
#plot(roc, main = "ROC Curve", ylab = "Sensitivity", xlab = "1-Specificity", )




