setwd("path/to/data/")
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

#cor_g <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = "correlation")

AA <- graph.adjacency(cor_mat, mode = "undirected", weighted = TRUE, diag = TRUE)
cor_edge_list <- get.data.frame(AA, "edges")

#cor_edge_list <- get.data.frame(cor_g, "edges")

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

###Adding adj.P.value and Log(FC) Features
IBD_genes <- read.table("IBD Genes.txt", header = TRUE)
lo_adj <- read.table("Source and Target Nodes Log(FC).txt", header = TRUE)
PPI_Core_net <- read.table("PPI_Core_net.txt", header = TRUE)
names(PPI_Core_net)[1]<-paste("Gene")
PPI_Core_net <- merge(PPI_Core_net, lo_adj, by = "Gene")
names(PPI_Core_net)[1]<-paste("From")
names(PPI_Core_net)[2]<-paste("Gene")
PPI_Core_net <- merge(PPI_Core_net, lo_adj, by = "Gene")
names(PPI_Core_net)[1]<-paste("To")
PPI_Core_net <- PPI_Core_net[,c(2,1,3,4,5,6,7,8)]
names(PPI_Core_net)[5]<-paste("Log(FC)(From)")
names(PPI_Core_net)[7]<-paste("Log(FC)(To)")
names(PPI_Core_net)[6]<- paste("adj.P.value(From)")
names(PPI_Core_net)[8]<- paste("adj.P.value(To)")
write.table(PPI_Core_net, "PPI_Core_net.txt", sep = "\t")







