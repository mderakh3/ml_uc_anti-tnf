library("pheatmap")
library("matrixStats")
library("data.table")

x <- read.table("DEGs_down.txt")
y <- read.table("DEGs_up.txt")
mdata.c <- read.table("mdata.c.txt")

ddw <- row.names(x)
duw <- row.names(y)
gene_names <- as.data.frame(append(ddw, duw))
length(intersect(gene_names$`append(ddw, duw)` , rownames(mdata.c)))
DEGs <- mdata.c[which(row.names(mdata.c) %in% gene_names$`append(ddw, duw)`),]
write.table(DEGs, "DEGs.txt", sep = "\t")

#coheda <- read.table("mdata.c.txt")
#coheda$row_var <- rowVars(as.matrix(coheda[,c(1:35)]))
#quantile(coheda$row_var, 0.5)
#coheda <- coheda[coheda$row_var > 0.3061521,]
#coheda$row_var <- NULL
gr <- c(rep("UCR", 4), rep("UCN", 7), rep("UCR", 8), rep("UCN", 16))

#zedScore Calculation
DEGs <- t(scale(t(DEGs)))

col <- c(colorRampPalette(c("#0000FF", "#FFFFFF"))(35),
         colorRampPalette(c("#FFFFFF","#FF0000"))(70))
pheatmap(cor(DEGs), labels_row = gr, labels_col = gr, color = col
         ,border_color = NA)

####Expression Heatmap
#x <- read.table("DEGs_weka_down.txt")
#y <- read.table("DEGs_weka_up.txt")
#ddw <- row.names(DEGs_down)
#duw <- row.names(DEGs_up)
#DEGs_heatmap <- as.data.frame(append(ddw, duw))
#heatmap_data <- data[which(row.names(data) %in% DEGs_heatmap$`append(ddw, duw)`),]
#heatmap_data <- as.matrix(heatmap_data)
#heatmap_data <- t(scale(t(heatmap_data)))

setwd("path/to/data/")
DEGs_weka_Expressions <- read.table("DEGs_weka_Expressions.txt")
DEGs_weka_Expressions <- DEGs_weka_Expressions[,c(1:4, 12:19, 5:11, 20:35)]
gr <- c(rep("UCR", 12), rep("UCN", 23))
DEGs_weka_Expressions <- as.matrix(DEGs_weka_Expressions)
colnames(DEGs_weka_Expressions) <- gr
lasfea <- read.table("Rfeature.txt")
DEGs_weka_Expressions <- DEGs_weka_Expressions[which(row.names(DEGs_weka_Expressions)
                                                   %in% lasfea$V2),]
heatmap.2(DEGs_weka_Expressions, col = bluered(256), scale = "row", trace = "none"
          , Rowv = TRUE, Colv = TRUE)

#Correlation Heatmap for 321 and 33 DEGs
setwd("path/to/data/")
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
         #,cluster_cols = FALSE, cluster_rows = TRUE
          )



