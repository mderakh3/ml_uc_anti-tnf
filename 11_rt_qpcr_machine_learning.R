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
library("MASS")
library("caret")

###Microarray Input preparation for ML
data <- read.delim("mdata.c.txt", header = T)
data <- tibble::column_to_rownames(data,"X")
features <- c("IL13RA2","MMP1","HCAR3", "CSF3", "INHBA")
data <- data[which(row.names(data) %in% features),]
data <- as.data.frame(t(data))
data <- data[,c(3,5,2,1,4)]
gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))
data <- mutate(data, gr = gr)
write.table(data, "In Silico Input for Machine Learning.txt", sep = "\t")

###Reading Data
#Microarray
data <- read.table("In Silico Input for ML.txt", header = TRUE)

#RT-qPCR
data <- read.table("qPCR Input for ML.txt", header = TRUE)

###Linear Discriminant Analysis
data <- data[,-7]
data1 <- data1[,-7]
y_train = data
y_test = data1
data_train <- data[row.names(y_train),-6]
data_test <- data1[row.names(y_test),-6]
LDA <- lda(y_train$gr ~., data = data_train)
pred <- predict(LDA, data_test,type = "class")
CM <- confusionMatrix(data = factor(pred$class),
                       reference = as.factor(y_test$gr), 
                       positive = "UCN")
fpr_fnr <- data.frame(gr1,pred,row.names(y_test))
write.table(fpr_fnr, file = "fpr and fnr.txt")


###Logistic Regression
y_train = data
y_test = data1
data_train <- data[row.names(y_train),]
data_test <- data1[row.names(y_test),]
data_train <- as.data.frame(data_train[,c(-6)])
data_test <- as.data.frame(data_test[,c(-6)])
mylogit <- glm(y_train$GR ~ ., data = data_train, family='binomial')
pred <- predict(mylogit ,newdata = data_test, type = "response")
CM <- confusionMatrix(data =  as.factor(as.numeric(pred > 0.5)),
                       reference = as.factor(y_test$GR1), 
                       positive = "1")
fpr_fnr <- data.frame(gr1,pred,row.names(data))
write.table(data1, file = "fpr and fnr.txt")


###Support Vector Machine
data <- data[,-7]
y_train = data
y_test = data1
y_train$gr <- as.factor(y_train$gr)
y_test$gr <- as.factor(y_test$gr)
data_train <- data[row.names(y_train),]
data_test <- data1[row.names(y_test),]
data_train <- as.matrix(data_train[,-6])
data_test <- as.matrix(data_test[,-6])
svm.model <- svm(y_train$gr ~ ., data = data_train)
svm.pred <- predict(svm.model, data_test)
CM <- confusionMatrix(data =  factor(svm.pred),
                       reference = as.factor(y_test$gr), 
                       positive = "UCN")

###Random Forest
data <- data[,-7]
data1 <- data1[,-7]
y_train = data
y_test = data1
data_train <- data[row.names(y_train),]
data_test <- data1[row.names(y_test),]
y_train$gr <- as.factor(y_train$gr)
y_test$gr <- as.factor(y_test$gr)
data_train <- as.matrix(data_train[,-6])
data_test <- as.matrix(data_test[,-6])
RF <- randomForest(x = data_train, y = y_train$gr, ntree=100)
pred_test <- predict(RF, newdata = data_test, type="class")
CM <- confusionMatrix(factor(pred_test), factor(y_test$gr), positive = "UCN")

###Cross Fold Validation

CM_table <- data.frame("Sensitivity"=0,"Specificity"=0,"Accuracy"=0,
                       "fold"=0, "NO.features"=0, "method" = "a")

for (i in 1:3) {
  #set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  data_train <- as.matrix(data_train[,c(-6,-7)])
  data_test <- as.matrix(data_test[,c(-6,-7)])
  
  new_train3 <- data_train[,c(1:3)]
  new_train4 <- data_train[,c(1:4)]
  new_train5 <- data_train[,c(1:5)]
  new_test3 <- data_test[,c(1:3)]
  new_test4 <- data_test[,c(1:4)]
  new_test5 <- data_test[,c(1:5)]
  
  rm(data_train,data_test)
  
  ###SVM
  svm.model2 <- svm(y_train$gr ~ ., data = new_train3)
  svm.model3 <- svm(y_train$gr ~ ., data = new_train4)
  svm.model4 <- svm(y_train$gr ~ ., data = new_train5)
  
  svm.pred3 <- predict(svm.model2, new_test3)
  svm.pred4 <- predict(svm.model3, new_test4)
  svm.pred5 <- predict(svm.model4, new_test5)
  
  CM3 <- confusionMatrix(data =  as.factor(svm.pred3),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM4 <- confusionMatrix(data =  as.factor(svm.pred4),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM5 <- confusionMatrix(data =  as.factor(svm.pred5),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"svm"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"svm"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"svm"))

  rm(svm.model2,CM3,svm.pred3,
     svm.model3,CM4,svm.pred4,
     svm.model4,CM5,svm.pred5)
  gc()
  
  ####LDA
  
  LDA3 <- lda(y_train$gr ~., data = as.data.frame(new_train3))
  LDA4 <- lda(y_train$gr ~., data = as.data.frame(new_train4))
  LDA5 <- lda(y_train$gr ~., data = as.data.frame(new_train5))
  
  pred3 <- predict(LDA3, as.data.frame(new_test3),type = "class")
  pred4 <- predict(LDA4, as.data.frame(new_test4),type = "class")
  pred5 <- predict(LDA5, as.data.frame(new_test5),type = "class")

  CM3 <- confusionMatrix(data = as.factor(pred3$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM4 <- confusionMatrix(data = as.factor(pred4$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM5 <- confusionMatrix(data = as.factor(pred5$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"LDA"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"LDA"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"LDA"))
  
  rm(CM3,pred3,LDA3,
     CM4,pred4,LDA4,
     CM5,pred5,LDA5)
  
  gc()
  
  ###Logistic Regression
  
  y_train$gr <- as.numeric(y_train$gr)
  y_test$gr <- as.numeric(y_test$gr)
  
  y_train$gr[y_train$gr == 2]<- 0
  y_test$gr[ y_test$gr == 2]<- 0
  
  new_train3 <- as.data.frame(new_train3)
  new_train4 <- as.data.frame(new_train4)
  new_train5 <- as.data.frame(new_train5)
  
  mylogit3 <- glm(y_train$gr ~ ., data = new_train3, family='binomial')
  mylogit4 <- glm(y_train$gr ~ ., data = new_train4, family='binomial')
  mylogit5 <- glm(y_train$gr ~ ., data = new_train5, family='binomial')
  
  new_test3 <- as.data.frame(new_test3)
  new_test4 <- as.data.frame(new_test4)
  new_test5 <- as.data.frame(new_test5)
  
  pred3<- predict(mylogit3 ,newdata = new_test3, type = "response")
  pred4<- predict(mylogit4 ,newdata = new_test4, type = "response")
  pred5<- predict(mylogit5 ,newdata = new_test5, type = "response")

  CM3 <- confusionMatrix(data =  as.factor(as.numeric(pred3 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  CM4 <- confusionMatrix(data =  as.factor(as.numeric(pred4 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  CM5 <- confusionMatrix(data =  as.factor(as.numeric(pred5 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"logit"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"logit"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"logit"))
  
  
  ###Random Forest
  data <- read.table("In Silico Input for ML.txt",header = TRUE)
  gr <- c(rep("UCN", 23), rep("UCR", 12))
  #set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  data_train <- as.matrix(data_train[,c(-6,-7)])
  data_test <- as.matrix(data_test[,c(-6,-7)])
  
  new_train3 <- data_train[,c(1:3)]
  new_train4 <- data_train[,c(1:4)]
  new_train5 <- data_train[,c(1:5)]
  new_test3 <- data_test[,c(1:3)]
  new_test4 <- data_test[,c(1:4)]
  new_test5 <- data_test[,c(1:5)]
  rm(data_train,data_test)
  
  RF3 <- randomForest(x = new_train3 , 
                      y = y_train$gr,
                      ntree=100)
  RF4 <- randomForest(x = new_train4 , 
                      y = y_train$gr,
                      ntree=100)
  
  RF5 <- randomForest(x = new_train5 , 
                      y = y_train$gr,
                      ntree=100)
  
  pred_test3 <- predict(RF3 ,newdata = new_test3 ,type="class")
  pred_test4 <- predict(RF4 ,newdata = new_test4 ,type="class")
  pred_test5 <- predict(RF5 ,newdata = new_test5 ,type="class")
  
  CM3 <- confusionMatrix(factor(pred_test3),factor(y_test$gr), positive = "UCN")
  CM4 <- confusionMatrix(factor(pred_test4),factor(y_test$gr), positive = "UCN")
  CM5 <- confusionMatrix(factor(pred_test5),factor(y_test$gr), positive = "UCN")

  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"RF"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"RF"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"RF"))
  
}

CM_table <- CM_table[-1,]
ME <- unique(CM_table$method)
NF <- as.numeric(unique(CM_table$NO.features))
CM_all <- data.frame("Sensitivity" = 0,"Specificity" = 0,
                     "Accuracy" = 0, "NO.features" = 0,
                     "method" = "a")

for (j in seq(length(ME))) {
  CM_table_j <- CM_table[CM_table$method == ME[j],]
  for (i in NF) {
    CM_f <- CM_table_j[CM_table_j$NO.features == i, c(1:3,5,6)]
    CM_all <- rbind(CM_all,c(round(mean(as.double(CM_f[,1])), 2),
                             round(mean(as.double(CM_f[,2])), 2),
                             round(mean(as.double(CM_f[,3])), 2),
                             i, ME[j]))
  }}
CM_all <- CM_all[-1,]
write.table(CM_all,"Machine Learning for of the Test In Silico.txt")

###qPCR Classification
data <- read.table("qPCR Input for ML.txt", header = TRUE)

CM_table <- data.frame("Sensitivity"=0,"Specificity"=0,"Accuracy"=0,
                       "fold"=0, "NO.features"=0, "method" = "a")

for (i in 1:2) {
  set.seed(5783125)
  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  data_train <- as.matrix(data_train[,c(-6,-7)])
  data_test <- as.matrix(data_test[,c(-6,-7)])
  
  new_train3 <- data_train[,c(1:3)]
  new_train4 <- data_train[,c(1:4)]
  new_train5 <- data_train[,c(1:5)]
  new_test3 <- data_test[,c(1:3)]
  new_test4 <- data_test[,c(1:4)]
  new_test5 <- data_test[,c(1:5)]
  
  rm(data_train,data_test)
  
  ###SVM
  svm.model2 <- svm(y_train$gr ~ ., data = new_train3)
  svm.model3 <- svm(y_train$gr ~ ., data = new_train4)
  svm.model4 <- svm(y_train$gr ~ ., data = new_train5)
  
  svm.pred3 <- predict(svm.model2, new_test3)
  svm.pred4 <- predict(svm.model3, new_test4)
  svm.pred5 <- predict(svm.model4, new_test5)
  
  CM3 <- confusionMatrix(data =  as.factor(svm.pred3),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM4 <- confusionMatrix(data =  as.factor(svm.pred4),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM5 <- confusionMatrix(data =  as.factor(svm.pred5),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"svm"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"svm"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"svm"))
  
####LDA
  
  LDA3 <- lda(y_train$gr ~., data = as.data.frame(new_train3))
  LDA4 <- lda(y_train$gr ~., data = as.data.frame(new_train4))
  LDA5 <- lda(y_train$gr ~., data = as.data.frame(new_train5))
  
  pred3 <- predict(LDA3, as.data.frame(new_test3),type = "class")
  pred4 <- predict(LDA4, as.data.frame(new_test4),type = "class")
  pred5 <- predict(LDA5, as.data.frame(new_test5),type = "class")
  
  CM3 <- confusionMatrix(data = factor(pred3$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM4 <- confusionMatrix(data = factor(pred4$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  CM5 <- confusionMatrix(data = factor(pred5$class),
                         reference = as.factor(y_test$gr), 
                         positive = "UCN")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"LDA"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"LDA"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"LDA"))
  
  rm(CM3,pred3,LDA3,
     CM4,pred4,LDA4,
     CM5,pred5,LDA5)
  
  gc()
  
  y_train$gr <- as.numeric(y_train$gr)
  y_test$gr <- as.numeric(y_test$gr)
  
  y_train$gr[y_train$gr == 2]<- 0
  y_test$gr[ y_test$gr == 2]<- 0
  
  new_train3 <- as.data.frame(new_train3)
  new_train4 <- as.data.frame(new_train4)
  new_train5 <- as.data.frame(new_train5)
  
  mylogit3 <- glm(y_train$gr ~ ., data = new_train3, family='binomial')
  mylogit4 <- glm(y_train$gr ~ ., data = new_train4, family='binomial')
  mylogit5 <- glm(y_train$gr ~ ., data = new_train5, family='binomial')
  
  new_test3 <- as.data.frame(new_test3)
  new_test4 <- as.data.frame(new_test4)
  new_test5 <- as.data.frame(new_test5)
  
  pred3<- predict(mylogit3 ,newdata = new_test3, type = "response")
  pred4<- predict(mylogit4 ,newdata = new_test4, type = "response")
  pred5<- predict(mylogit5 ,newdata = new_test5, type = "response")
  
  CM3 <- confusionMatrix(data =  as.factor(as.numeric(pred3 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  CM4 <- confusionMatrix(data =  as.factor(as.numeric(pred4 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  CM5 <- confusionMatrix(data =  as.factor(as.numeric(pred5 > 0.5)),
                         reference = as.factor(y_test$gr), 
                         positive = "1")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"logit"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"logit"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"logit"))
  
  
  ###Random Forest
  data <- read.table("qPCR Input for ML.txt",header = TRUE)

  y_train = data[data$fold != i,]
  y_test = data[data$fold == i,]
  data_train <- data[row.names(y_train),]
  data_test <- data[row.names(y_test),]
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  set.seed(100)
  
  y_train$gr <- as.factor(y_train$gr)
  y_test$gr <- as.factor(y_test$gr)
  data_train <- as.matrix(data_train[,c(-6,-7)])
  data_test <- as.matrix(data_test[,c(-6,-7)])
  
  new_train3 <- data_train[,c(1:3)]
  new_train4 <- data_train[,c(1:4)]
  new_train5 <- data_train[,c(1:5)]
  new_test3 <- data_test[,c(1:3)]
  new_test4 <- data_test[,c(1:4)]
  new_test5 <- data_test[,c(1:5)]
  rm(data_train,data_test)
  
  RF3 <- randomForest(x = new_train3 , 
                      y = y_train$gr,
                      ntree=100)
  RF4 <- randomForest(x = new_train4 , 
                      y = y_train$gr,
                      ntree=100)
  
  RF5 <- randomForest(x = new_train5 , 
                      y = y_train$gr,
                      ntree=100)
  
  pred_test3 <- predict(RF3 ,newdata = new_test3 ,type="class")
  pred_test4 <- predict(RF4 ,newdata = new_test4 ,type="class")
  pred_test5 <- predict(RF5 ,newdata = new_test5 ,type="class")
  
  CM3 <- confusionMatrix(factor(pred_test3),factor(y_test$gr), positive = "UCN")
  CM4 <- confusionMatrix(factor(pred_test4),factor(y_test$gr), positive = "UCN")
  CM5 <- confusionMatrix(factor(pred_test5),factor(y_test$gr), positive = "UCN")
  
  CM_table <- rbind(CM_table,c(CM3$byClass[c(1,2)],CM3$overall[1],i,3,"RF"))
  CM_table <- rbind(CM_table,c(CM4$byClass[c(1,2)],CM4$overall[1],i,4,"RF"))
  CM_table <- rbind(CM_table,c(CM5$byClass[c(1,2)],CM5$overall[1],i,5,"RF"))
  
 
}

CM_table <- CM_table[-1,]
ME <- unique(CM_table$method)
NF <- as.numeric(unique(CM_table$NO.features))
CM_all <- data.frame("Sensitivity" = 0,"Specificity" = 0,
                     "Accuracy" = 0, "NO.features" = 0,
                     "method" = "a")

for (j in seq(length(ME))) {
  CM_table_j <- CM_table[CM_table$method == ME[j],]
  for (i in NF) {
    CM_f <- CM_table_j[CM_table_j$NO.features == i, c(1:3,5,6)]
    CM_all <- rbind(CM_all,c(round(mean(as.double(CM_f[,1])), 2),
                             round(mean(as.double(CM_f[,2])), 2),
                             round(mean(as.double(CM_f[,3])), 2),
                             i, ME[j]))
  }}
CM_all <- CM_all[-1,]
write.table(CM_all,"Machine Learning for CinnoRA.txt")


###qPCR ROC Curve
data <- read.table("ML.txt", header = TRUE)
gr <- c(rep("UCR", 12), rep("UCN", 10))

ROC_table <- data.frame("response"= "a","predictor"= 0,
                        "NO.features"= 0,
                        "method" = "a")

  data$gr <- as.factor(data$gr)
  
  data3 <- data[,1:3]
  data4 <- data[,1:4]
  data5 <- data[,1:5]
  
  # LDA
  
  LDA3 <- lda(data$gr ~., data = data3, probability = TRUE)
  LDA4 <- lda(data$gr ~., data = data4, probability = TRUE)
  LDA5 <- lda(data$gr ~., data = data5, probability = TRUE)
  
  pred3 <- predict(LDA3, as.data.frame(data3),probability = TRUE)
  pred4 <- predict(LDA4, as.data.frame(data4),probability = TRUE)
  pred5 <- predict(LDA5, as.data.frame(data5),probability = TRUE)
  
  LDA_pr3  <- as.data.frame(pred3$posterior)
  LDA_pr4  <- as.data.frame(pred4$posterior)
  LDA_pr5  <- as.data.frame(pred5$posterior)
  
  ROC3 <- data.frame("response"=as.character(data$gr),"predictor"= LDA_pr3$UCR,
                     "NO.features"= rep(3,length(LDA_pr3$UCN)),
                     "method" = rep("LDA",length(LDA_pr3$UCN)))
  ROC4 <- data.frame("response"=as.character(data$gr),"predictor"= LDA_pr4$UCN,
                     "NO.features"= rep(4,length(LDA_pr4$UCN)),
                     "method" = rep("LDA",length(LDA_pr4$UCN)))
  ROC5 <- data.frame("response"=as.character(data$gr),"predictor"= LDA_pr5$UCN,
                     "NO.features"= rep(5,length(LDA_pr5$UCN)),
                     "method" = rep("LDA",length(LDA_pr5$UCN)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
  
  ##logit
  data$gr <- as.numeric(data$gr)
  data$gr[data$gr == 2]<- 0
  
  mylogit3 <- glm(data$gr ~ ., data = data3, family='binomial')
  mylogit4 <- glm(data$gr ~ ., data = data4, family='binomial')
  mylogit5 <- glm(data$gr ~ ., data = data5, family='binomial')
  
  pred3<- predict(mylogit3 ,newdata = data3,type = "response", probability = TRUE)
  pred4<- predict(mylogit4 ,newdata = data4,type = "response", probability = TRUE)
  pred5<- predict(mylogit5 ,newdata = data5,type = "response", probability = TRUE)
  
  ROC3 <- data.frame("response"=data$gr,"predictor"= pred3,
                     "NO.features"= rep(3,length(pred3)),
                     "method" = rep("logit",length(pred3)))
  ROC4 <- data.frame("response"=data$gr,"predictor"= pred4,
                     "NO.features"= rep(4,length(pred4)),
                     "method" = rep("logit",length(pred4)))
  ROC5 <- data.frame("response"=data$gr,"predictor"= pred5,
                     "NO.features"= rep(5,length(pred5)),
                     "method" = rep("logit",length(pred5)))
  
  ROC_table <- rbind(ROC_table,ROC3)
  ROC_table <- rbind(ROC_table,ROC4)
  ROC_table <- rbind(ROC_table,ROC5)
  
  
  rm(ROC3,pred3,data3_logit,mylogit3,new_train3,new_test3,pred3,
     ROC4,pred4,data4_logit,mylogit4,new_train4,new_test4,pred4,
     ROC5,pred5,data5_logit,mylogit5,new_train5,new_test5,pred5)
  

ROC_table <- ROC_table[-1,]
write.table(ROC_table,"ROC Data, dct.txt", sep = "\t", row.names = F)
ROC_table <- read.delim("path/to/roc_curve_analysis/ROC Data, dct.txt")

ME <- unique(ROC_table$method)
NF <- as.numeric(unique(ROC_table$NO.features))
FO <- as.numeric(unique(ROC_table$fold))

ROC_lda <- ROC_table[ROC_table$method == "LDA",]
ROC_logit <- ROC_table[ROC_table$method == "logit",]


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














