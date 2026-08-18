#################Data Pre-processing, Batch Effect Removal, and Data Integration
setwd("...")

library("limma")
library("plotly")
library("data.table")
library("dplyr")
library("ggplot2")
library("sva")

data1 <- read.delim("GSE12251_data.txt")
data2 <- read.delim("GSE16879_data.txt")
annot <- read.delim("GPL570.txt")
colnames(annot)[1] <- "ID_REF"
data1 <- merge(annot, data1, by = "ID_REF")
data2 <- merge(annot, data2)
data1 <- data1[,-c(1:2)]
data2 <- data2[,-c(1:2)]

for(j in 2:ncol(data1)){
  if(typeof(data1[,j]) != "double"){print(paste("j=",j))}
}

for(j in 2:ncol(data2)){
  if(typeof(data2[,j]) != "double"){print(paste("j=",j))}
}

data1[is.na(data1)]<- 0
data2[is.na(data2)]<- 0
dim(data1)
dim(data2)

data10 <- aggregate(by = list(unique.values = data1$GeneSymbol)
                    , x = data1, FUN = mean)
data20 <- aggregate(by = list(unique.values = data2$GeneSymbol)
                    , x = data2, FUN = mean)

dim(data10)
dim(data20)

data10$GeneSymbol <- NULL
data20$GeneSymbol <- NULL

rownames(data10) <- data10$unique.values
rownames(data20) <- data20$unique.values

data10$unique.values <- NULL
data20$unique.values <- NULL

if(max(data10 > 100)){
  data10 <- log2(data10+1)}

if(max(data20 > 100)){
  data20 <- log2(data20+1)}

max(data10)
min(data10)
max(data20)
min(data20)

data10 <- normalizeQuantiles(data10)
data20 <- normalizeQuantiles(data20)

mdata <- merge(data10, data20, by = 0)

row.names(mdata) <- mdata$Row.names
mdata$Row.names <- NULL
mdata <- mdata[,c(2,6,7,10,1,3,4,5,8,9,11,12:35)]
gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))

colnames(mdata) <- gr
boxplot(mdata, las = 2, col = c(rep("blue", 11), rep("red", 24)),
        main = "Samples' Value Distributions")

Datasets <- factor(c(rep(c("GSE12251","GSE16879"), 
                         times=c(ncol(data10),ncol(data20)))))

pc1 <- prcomp(mdata)
pcr1 <- data.frame(pc1$r[,1:3],Datasets, gr)

ggplot(pcr1,aes(x = PC1, y = PC2, colour=Datasets, shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

mdata.c <- as.data.frame(ComBat(as.matrix(mdata),Datasets))

pc1 <- prcomp(mdata.c)
pcr1 <- data.frame(pc1$r[,1:3],Datasets, gr)
ggplot(pcr1,aes(x = PC1, y = PC2, colour=Datasets,shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))

boxplot(mdata.c, las = 2, col = c(rep("blue", 11), rep("red", 24)),
        main = "Samples' Value Distributions")

write.table(mdata.c, file = "mdata.c.txt", sep = "  ")
save(mdata.c, file = "Metadata")
load("Metadata")

#########################################Gene Expression Analysis
library("limma")
library("GEOquery")
library("ggplot2")
library("reshape2")
library("plyr")
library("gplots")
library("Biobase")

ex <- read.table(file = "mdata.c.txt", header = TRUE)

gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))

dif.ex <- function(ex, gr, con, ann=NULL, number=Inf) {
  if (! is.null(ann)) rownames(ex) <- 1:nrow(ex)
  gr <- factor(gr)
  mat  <- data.frame(G=gr)
  design <- model.matrix(~G+0,mat)
  colnames(design) <- levels(gr)
  fit <- lmFit(ex, design)
  cont.matrix <- makeContrasts(contrasts=con, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01,trend=TRUE,robust=TRUE)
  tT=topTable(fit2, adjust="fdr", number=number)
  if (! is.null(ann)) tT <- cbind(tT, ann[as.numeric(rownames(tT)),,drop=F])
  tT
}

DEGs <- dif.ex(ex = ex, gr = gr, con = "UCN-UCR")
DEGs_down <- DEGs[DEGs$adj.P.Val < 0.05,]
DEGs_down <- DEGs_down[DEGs_down$logFC < -1,]
DEGs_up <- DEGs[DEGs$adj.P.Val < 0.05,]
DEGs_up <- DEGs_up[DEGs_up$logFC > 1,]

write.table(DEGs_down, file = "UCNB vs UCNA.down.txt", sep = " ")
write.table(DEGs_up, file = "UCNB vs UCNA.up.txt", sep = " ")

#########################################Correlation Heatmap
library("pheatmap")
library("matrixStats")
library("data.table")

mdata.c <- read.table("mdata.c.txt")
ddw <- read.table("DEGs_down.txt")
ddu <- read.table("DEGs_up.txt")
ddw <- rownames(ddw)
ddu <- rownames(ddu)
genes <- as.data.frame(append(ddw, ddu))
cohege <- mdata.c[which(row.names(mdata.c) %in% genes$`append(ddw, ddu)`),]
cohege <- as.data.frame(cohege)

gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))

cohege <- cohege[,order(names(cohege))]
gr <- c(rep("UCN", 23), rep("UCR", 12))
colnames(cohege) <- gr

cohege <- as.matrix(cohege)
cohege <- t(scale(t(cohege)))
col <- c(colorRampPalette(c("#0000FF", "#FFFFFF"))(45),
         colorRampPalette(c("#FFFFFF","#FF0000"))(70))
pheatmap(cor(cohege), labels_row = gr, labels_col = gr, color = col
         ,border_color = NA
)

########################################Pathway Bubble Plot 

library("ggplot2")
library("dplyr")
library("plotly")
library("viridis")
library("hrbrthemes")

data <- read.table("pth.txt", header = TRUE)
data <- data[,c(1,2,3)]
names(data)[2] <- paste("adj.P.value")

p = ggplot(data,aes(adj.P.value, Pathway))
p = p + geom_point()  
p = p + geom_point(aes(size = Gene.Ratio))
pbubble = p + geom_point(aes(size = Gene.Ratio, color = adj.P.value))
pr = pbubble + scale_color_gradient(low = "blue", high = "red")
pr = pr + labs(color = expression ("-Log10(adj.P.value)"),size="Gene.Ratio",x="-Log10(adj.P.value)",y="Pathway",title="Pathway Analysis")+theme(plot.title=element_text(hjust=0.5))
pr = pr + theme_bw()

########################################Pearson Correlation Calculation

library("foreign")
library("igraph")
library("dplyr")

x <- read.table("DEGs_down.txt")
y <- read.table("DEGs_up.txt")
mdata.c <- read.table("mdata.c.txt")
ddw <- row.names(x)
duw <- row.names(y)

core_matrix <- as.data.frame(append(ddw, duw))
length(intersect(core_matrix$`append(ddw, duw)` , rownames(mdata.c)))
core_matrix <- mdata.c[which(row.names(mdata.c) %in% core_matrix$`append(ddw, duw)`),]
core_matrix <- t(core_matrix)
dim(core_matrix)

cor_mat <- cor(core_matrix, method = "pearson")
dim(cor_mat)

AA <- graph.adjacency(cor_mat, mode = "undirected", weighted = TRUE, diag = TRUE)
cor_edge_list <- get.data.frame(AA, "edges")

net <- cor_edge_list[cor_edge_list$weight != 1,]

for (i in nrow(net):1) {
  if (net[i,1] == net[i,2]) {net <- net[-i,]}
}

only_sig_pos <- net[abs(net$weight) > 0.9,]
only_sig_neg <- net[abs(net$weight) < -0.9,]

nodes <- unique(c(only_sig$from, only_sig$to))

write.table(only_sig, file = "CorrelationEdge.txt", sep = "\t", row.names = FALSE)

new_g <- graph_from_data_frame(only_sig, F)

only_sig <- read.table("CorrelationEdge.txt", header = TRUE)
ProteinPI <- read.table("ProteinPI.txt", header = TRUE)

for (i in nrow(ProteinPI):1) {
  if (ProteinPI[i,1&2] == only_sig[i,1&2]) {
    new <- only_sig[i, "correlation"]
    ProteinPI[i, ncol(ProteinPI) + 1] <- new
    colnames(ProteinPI)[ncol(ProteinPI)] <- paste0("correlation")
  }
}

###PPI_core_net file Construction
only_sig <- read.table("2.txt", header = TRUE)
ProteinPI <- read.table("1.txt", header = TRUE)
only_sig <- only_sig[,c(-3,-4)]
ProteinPI <- ProteinPI[,c(-3,-4)]
PPI_Core_net <- merge(ProteinPI, only_sig, by = "Interactions")
PPI_Core_net <- distinct(PPI_Core_net)
rownames(PPI_Core_net) <- PPI_Core_net$Interactions
PPI_Core_net <- PPI_Core_net[,-1]
write.table(PPI_Core_net, "PPI_Core_net.txt", sep = "\t")

########################################Network Construction
setwd("path/to/PPI_network/")
library("foreign")
library("igraph")
library("dplyr")

#Nodes Attributes
Proteins <- read.table("Proteins.txt", header = TRUE)
DEGs <- read.table("DEGs.txt", header = TRUE)
nodes <- merge(Proteins, DEGs, by = "Genes")
write.table(nodes, "NodesAttributes1.txt", sep = "\t" ,row.names = FALSE)

#Edge Attributes
PPICoreEdges <- read.table("PPICoreEdges.txt", header = TRUE)
fin_net <- read.table("Venn 100 top genes.txt", header = TRUE)
NET <- PPICoreEdges[PPICoreEdges[,1] %in% fin_net$Genes,] 
NET <- NET[NET[,2] %in% fin_net$Genes,]   

write.table(NET, "fin_net.txt", sep = "\t" ,row.names = FALSE)

NET_genes_logFC <- merge(fin_net, DEGs, by = "Genes")
write.table(NET_genes_logFC, "NET_genes_logFC.txt", sep = "\t" ,row.names = FALSE)

#################################################Lasso Feature Selection

library("randomForest")
library("doParallel")
library("caret")
library("e1071")
library("glmnet")
library("doParallel")

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


####################################################Ensemble Machine Learning
library("caret")
library("caretEnsemble")
library("pROC")

InSilic_data <- read.table("In Silico Input for ML.txt",header = TRUE)
qPCR_data <-read.table("qPCR Input for ML.txt",header = TRUE)
indep_data <- read.table("Independent In Silico Data.txt",header = TRUE)

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

  # Ensemble the predictions of `models` to form a new combined prediction based on glm
  set.seed(seed)
  stack.glm <- caretStack(models,
                          method='glm',
                          metric="Accuracy",
                          trControl=stackControl ) 
  
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


#Accuracy Comparison plot
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
    geom_text(aes(label=(NewAcc[,"Method"]), y = 100), position = position_dodge(0.90))

dev.off()

}
