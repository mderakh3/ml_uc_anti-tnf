##https://www.analyticsvidhya.com/blog/2017/02/introduction-to-ensembling-along-with-implementation-in-r/
##https://topepo.github.io/caret/model-training-and-tuning.html#final


##Best seed for ML for insilico model= 955904


library("caret")
library(caretEnsemble)
library(pROC)

setwd("path/to/data/")
InSilic_data <- read.table("In Silico Input for ML.txt",header = TRUE)
qPCR_data <-read.table("qPCR Input for ML.txt",header = TRUE)
indep_data <- read.table("MLI Zscore -1 UCR sample.txt",header = TRUE)
##############################################


###############################################

algorithmList <- c('lda', 'rf', 'nnet', 'knn')
trainControl <- trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=10,
                             savePredictions='final',
                             classProbs=TRUE)

stackControl <- trainControl(method="repeatedcv",
                             number=10,
                             repeats=10,
                             savePredictions='final',
                             classProbs=TRUE)
seed <- 955904   

Overal_results = data.frame("Method"= rep(c(algorithmList,'glm'),3),
                           "Predictor.No"= rep(0,15),
                           "Accuracy" = rep(0,15),
                           "Kappa" = rep(0,15),
                           "Sensitivity" = rep(0,15),
                           "Specificity" = rep(0,15),
                           "Accuracy_indep" = rep(0,15),
                           "Sensitivity_indep" = rep(0,15),
                           "Specificity_indep" = rep(0,15))
feature_No = 5

qPCR_data <- qPCR_data[,c(1:feature_No,6)]
indep_data <- indep_data[,c(1:feature_No,6)]
InSilic_data <- InSilic_data[,c(1:feature_No,6)]
##############################################
plotname <- "Insilic"
k=0

for (f in c(3:feature_No)) {
    
    data <- InSilic_data[,c(1:f,6)]
    
    testdata <- indep_data[,c(1:f,6)]
    k=k+1;
    set.seed(seed)
    models <- caretList(gr ~.  ,
                        data=data,
                        trControl=trainControl,
                        methodList=algorithmList,
                        metric="Accuracy")
    
    results <- resamples(models)
    su <- summary(results)
    for(j in seq(4))
    {
      Overal_results[k,"Predictor.No"] = f
      Overal_results[k,"Accuracy"] = round(su[["statistics"]][["Accuracy"]][,"Mean"],2)[[j]]
      Overal_results[k,"Kappa"] = round(su[["statistics"]][["Kappa"]][,"Mean"],2)[[j]]
      k=k+1;
    }
    # Box plots to compare models
    #scales <- list(x=list(relation="free"), y=list(relation="free"))
    #bwplot(results, scales=scales)
    
    # Ensemble the predictions of `models` to form a new combined prediction based on glm
    set.seed(seed)
    stack.glm <- caretStack(models,
                            method='glm',
                            metric="Accuracy",
                            trControl=stackControl ) 
    
    #print(stack.glm)
    Overal_results[k,"Predictor.No"] = f
    Overal_results[k,"Accuracy"] = round(stack.glm[["ens_model"]][["results"]][["Accuracy"]],2)[[1]]
    Overal_results[k,"Kappa"] = round(stack.glm[["ens_model"]][["results"]][["Kappa"]],2)[[1]]

    self_CM <- confusionMatrix(data = as.factor(stack.glm[["ens_model"]][["pred"]][["pred"]]),
                reference = as.factor(stack.glm[["ens_model"]][["pred"]][["obs"]]))
    
    Overal_results[k,"Sensitivity"] = round(self_CM[["byClass"]][["Sensitivity"]],2)
    Overal_results[k,"Specificity"] = round(self_CM[["byClass"]][["Specificity"]],2)
    
    stack_predicteds_CM <- predict(stack.glm, newdata=testdata)
    indep_CM<-confusionMatrix(as.factor(indep_data$gr),as.factor(stack_predicteds_CM))
    
    Overal_results[k,"Accuracy_indep"] = round(indep_CM[["overall"]][["Accuracy"]],2)
    Overal_results[k,"Sensitivity_indep"] = round(indep_CM[["byClass"]][["Sensitivity"]],2)
    Overal_results[k,"Specificity_indep"] = round(indep_CM[["byClass"]][["Specificity"]],2)
    
    selfrocCurve <- roc(as.factor(stack.glm[["ens_model"]][["pred"]][["obs"]]),
                        stack.glm[["ens_model"]][["pred"]][["UCN"]],auc = TRUE)

    stack_predicteds <- predict(stack.glm, newdata = indep_data, type = "prob")
    rocCurve <- roc(as.factor(indep_data$gr), stack_predicteds, auc =TRUE)
    if(f==3)
    {
      
      tiff(paste(plotname,".tiff"),width = 2000,height = 2000,res=500)
      plot.roc(selfrocCurve,col=c(4),
               print.auc = TRUE,
               print.auc.y = 0.1,
               print.auc.x = 0.3)
    }
    else
      {
        if (f==4)
          plot.roc(selfrocCurve,col=c(3),
                   print.auc = TRUE,
                   print.auc.y = 0.2,
                   print.auc.x = 0.3,
                   add = TRUE)
        else
          if (f==5)
          {
            plot.roc(selfrocCurve,col=c(2),
                     print.auc = TRUE,
                     print.auc.y = 0.3,
                     print.auc.x = 0.3,
                     add = TRUE)
            dev.off()
        }
      }
}
dev.off()
write.table(Overal_results,paste(plotname,"_Results.txt"),sep="\t")



##########################################
##Accuracy Comparison plot
########################################
{
  library("ggplot2")
  NewAcc <- data.frame("model"=rep(c(algorithmList,'glm'),3),
                       "No.feature"= c(rep(3,length(algorithmList)+1),
                                       rep(4,length(algorithmList)+1),
                                       rep(5,length(algorithmList)+1)),
                       "Accuracy"=c(rep(0,(length(algorithmList)+1)*3)))
 
  
  NewAcc <- Overal_results[,c("Method","Predictor.No","Accuracy")]
 
  tiff(paste(plotname,"333444_Accuracy_Comparision.tiff"),width = 2000,height = 1500,res=500)
  
  ggplot(NewAcc, aes(x = factor(Predictor.No), y = Accuracy, 
                     fill = Method, colour = Method)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(x="Predictor.No - Model", y="Accuracy") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = c(0.9, 0.2)) +
    ylim(0,1) +
    geom_text(aes(label=(NewAcc[,"Method"]), y = 100), position = position_dodge(0.90)))
  
  dev.off()
  
}
