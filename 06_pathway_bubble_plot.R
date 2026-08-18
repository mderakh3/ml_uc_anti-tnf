setRepositories()
setwd("path/to/data/")
library(ggplot2)
library(dplyr)
library(plotly)
library(viridis)
library(hrbrthemes)

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
      
                                          