# Load R package
library(Seurat)
library(ggplot2)

######------Figure 2A------######
cors <- c("#8DD3C7","#B3DE69","grey","#F490A9","#FDB499","skyblue1","#84B1ED")
pdf(file=paste0('Figure 2A.pdf'),width = 6,height = 8)
VlnPlot(hip, features = features, 
        stack = TRUE, 
        sort = TRUE, 
        cols = cors,
        split.by =  "Celltype" , 
        flip = TRUE,
        split.plot = TRUE) +
  theme(legend.position = "none")
dev.off()
