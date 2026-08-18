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
