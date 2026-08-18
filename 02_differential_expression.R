setwd("path/to/data/")

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



