setRepositories()
library("limma")
library("plotly")
library("data.table")
library("dplyr")
library("ggplot2")
library("sva")

setwd("path/to/data/")

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




