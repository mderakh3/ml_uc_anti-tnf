setwd("path/to/independent_cohort/")
library("limma")
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
###############Sample Selcetion
#data <- read.delim("GSE73661.txt", header = TRUE)
#samples <- read.delim("Samples.txt", header = TRUE)
#rownames(samples) <- samples$X
#samples$X <- NULL
#rownames(data) <- data$ID_REF
#data$ID_REF <- NULL
#tdata <- as.data.frame(t(data))
#dataf <- merge(samples, tdata, by = "row.names")
#rownames(dataf) <- dataf$Row.names
#dataf$Row.names <- NULL
#dataf$gr <- NULL
#dataf <- as.data.frame(t(dataf))
#write.table(dataf, file = "dataf.txt", sep = "\t")


###########Preparing Annotation File
#library(org.Mm.eg.db) 
#library(annotate)
#library(igraph)
#annot <- read.delim("GPL6244.txt", header = TRUE)
#a <- as.character(annot$ID)
#b <- getSymbols(c(a), data = 'org.Mm.eg.db')
#d <- as.data.frame(b)
#colnames(d) <- "Gene.Symbol"
#d[is.na(d)]<- 0
#d <- subset(d, d$`Gene.Symbol` != 0)

##########Pre Processing
#data <- read.delim("dataf.txt", header = TRUE)
#annot <- read.delim("GPL6244.txt", header = TRUE)
#colnames(annot)[1] <- "ID_REF"
#mdata <- merge(annot, data, by = "ID_REF")
#
#for(j in 2:ncol(mdata)){
#  if(typeof(mdata[,j]) != "double"){print(paste("j=",j))}
#}
#
#mdata <- aggregate(by = list(unique.values = mdata$Gene.symbol)
#                    , x = mdata, FUN = mean)
#mdata$Gene.symbol <- NULL
#rownames(mdata) <- mdata$unique.values
#mdata$unique.values <- NULL
#mdata$ID_REF <- NULL
#
#if(max(mdata > 100)){
#  mdata <- log2(mdata+1)
#}
#
#mdata <- normalizeQuantiles(mdata)
#write.table(mdata, "mdata.txt", sep = "\t")


######PCA UCR Samples Filtering
mdata <- read.delim("PCAGeneQuery.txt", header = TRUE)
genes <- read.table("PCA Genes.txt", header = TRUE)
mdata <- merge(mdata, genes, by = "Genes")
gr <- c(rep("UCN", 3), rep("UCR", 2), rep("UCN", 2), "UCR",
        rep("UCN", 3), rep("UCR", 2), rep("UCN", 2), "UCR", rep("UCN", 2),
        rep("UCR", 2), rep("UCN", 3))
rownames(mdata) <- mdata$Genes
mdata$Genes <- NULL
pc1 <- prcomp(mdata)
pcr1 <- data.frame(pc1$r[,1:3], gr)
ggplot(pcr1,aes(x = PC1, y = PC3, colour=gr,shape = gr))+ geom_point(size=3)+
  theme_bw()+theme(legend.text=element_text(size=13,family = "serif"),
                   legend.title = element_text(size=15,family = "serif"),
                   axis.title =element_text(size=15,family = "serif"))
mdata <- mdata[,-19]
write.table(mdata, "ML Input -1 UCR sample.txt", sep = "\t")

###########Expression Analysis
#library("limma")
#library("GEOquery")
#library("ggplot2")
#library("reshape2")
#library("plyr")
#library("gplots")
#library("Biobase")
#ex <- read.table(file = "mdata.txt", header = TRUE)
#gr <- c(rep("UCN", 3), rep("UCR", 2), rep("UCN", 2), "UCR",
#        rep("UCN", 3), rep("UCR", 2), rep("UCN", 2), "UCR", rep("UCN", 2),
#        rep("UCR", 2), rep("UCN", 3))
#dif.ex <- function(ex, gr, con, ann=NULL, number=Inf) {
#  if (! is.null(ann)) rownames(ex) <- 1:nrow(ex)
#  gr <- factor(gr)
#  mat  <- data.frame(G=gr)
#  design <- model.matrix(~G+0,mat)
#  colnames(design) <- levels(gr)
#  fit <- lmFit(ex, design)
#  cont.matrix <- makeContrasts(contrasts=con, levels=design)
#  fit2 <- contrasts.fit(fit, cont.matrix)
#  fit2 <- eBayes(fit2, 0.01,trend=TRUE,robust=TRUE)
#  tT=topTable(fit2, adjust="fdr", number=number)
#  if (! is.null(ann)) tT <- cbind(tT, ann[as.numeric(rownames(tT)),,drop=F])
#  tT
#}
#DEGs <- dif.ex(ex = ex, gr = gr, con = "UCN-UCR")
#DEGs_down <- DEGs[DEGs$adj.P.Val < 0.05,]
#DEGs_down <- DEGs_down[DEGs_down$logFC < -1,]
#DEGs_up <- DEGs[DEGs$adj.P.Val < 0.05,]
#DEGs_up <- DEGs_up[DEGs_up$logFC > 1,]


############Scaling
a <- read.delim("ML Input -1 UCR sample.txt", header = TRUE)
a <- t(a)
#rownames(a) <- a$X
#a$X <- NULL
#write.table(a, "ML Input.txt", sep = "\t")
#a <- read.delim("ML Input.txt", header = TRUE)
#a <- a[,-6]
a <- scale(a)
a <- t(a)
write.table(a, "MLI Zscore -1 UCR sample.txt", sep = "\t")


###########ML Input Data
data1 <- read.delim("In Silico Input for ML.txt", header = TRUE)
gr1 <- c(rep("UCN", 23), rep("UCR", 12))
rownames(data1) <- data1$X
data1$X <- NULL
data1 <- data1[,-7]
data2 <- read.delim("MLI Zscore -1 UCR sample.txt", header = TRUE)
gr2 <- c(rep("UCN", 15), rep("UCR", 8))
rownames(data2) <- data2$X
data2$X <- NULL


###################Linear Discriminant Analysis
y_train <- data1
y_test <- data2
data_train3 <- data1[,1:3]
data_train4 <- data1[,1:4]
data_train5 <- data1[,1:5]
data_test3 <- data2[,1:3]
data_test4 <- data2[,1:4]
data_test5 <- data2[,1:5]
LDA <- lda(y_train$gr ~., data = data_train5)
pred <- predict(LDA, data_test5,type = "class")
CM <- confusionMatrix(data = factor(pred$class),
                      reference = as.factor(y_test$gr), 
                      positive = "UCN")
fpr_fnr <- data.frame(gr2,pred,row.names(y_test))
write.table(fpr_fnr, file = "fpr and fnr.txt")


##############ROC Curve
ROC_table <- data.frame("response"= "a","predictor"= 0,
                        "NO.features"= 0,
                        "method" = "a")

data1$gr <- as.factor(data1$gr)
data2$gr <- as.factor(data2$gr)

data_train3 <- data1[,1:3]
data_train4 <- data1[,1:4]
data_train5 <- data1[,1:5]

data_test3 <- data2[,1:3]
data_test4 <- data2[,1:4]
data_test5 <- data2[,1:5]

# LDA

LDA3 <- lda(data1$gr ~., data = data_train3, probability = TRUE)
LDA4 <- lda(data1$gr ~., data = data_train4, probability = TRUE)
LDA5 <- lda(data1$gr ~., data = data_train5, probability = TRUE)

pred3 <- predict(LDA3, as.data.frame(data_test3),probability = TRUE)
pred4 <- predict(LDA4, as.data.frame(data_test4),probability = TRUE)
pred5 <- predict(LDA5, as.data.frame(data_test5),probability = TRUE)

LDA_pr3  <- as.data.frame(pred3$posterior)
LDA_pr4  <- as.data.frame(pred4$posterior)
LDA_pr5  <- as.data.frame(pred5$posterior)

ROC3 <- data.frame("response"=as.character(data2$gr),"predictor"= LDA_pr3$UCR,
                   "NO.features"= rep(3,length(LDA_pr3$UCN)),
                   "method" = rep("LDA",length(LDA_pr3$UCN)))
ROC4 <- data.frame("response"=as.character(data2$gr),"predictor"= LDA_pr4$UCN,
                   "NO.features"= rep(4,length(LDA_pr4$UCN)),
                   "method" = rep("LDA",length(LDA_pr4$UCN)))
ROC5 <- data.frame("response"=as.character(data2$gr),"predictor"= LDA_pr5$UCN,
                   "NO.features"= rep(5,length(LDA_pr5$UCN)),
                   "method" = rep("LDA",length(LDA_pr5$UCN)))

ROC_table <- rbind(ROC_table,ROC3)
ROC_table <- rbind(ROC_table,ROC4)
ROC_table <- rbind(ROC_table,ROC5)

rm(ROC3,pred3,data3_logit,mylogit3,new_train3,new_test3,pred3,
   ROC4,pred4,data4_logit,mylogit4,new_train4,new_test4,pred4,
   ROC5,pred5,data5_logit,mylogit5,new_train5,new_test5,pred5)


ROC_table <- ROC_table[-1,]

ME <- unique(ROC_table$method)
NF <- as.numeric(unique(ROC_table$NO.features))
FO <- as.numeric(unique(ROC_table$fold))

ROC_lda <- ROC_table[ROC_table$method == "LDA",]

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



