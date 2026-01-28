
######------Figure 2M, N and Supplementary Figure 7------######
# Load R package
library(ggplot2)
library(enrichplot)

# drawing
plot <- gseaplot2(kk2,
                  title = "HIF-1 signaling pathway", 
                  "hsa04066", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
ggsave("Figure 2M.pdf", plot, width = 8, height = 8)
ggsave("Figure 2M.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "HIF-1 signaling pathway", 
                  "hsa04066", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
ggsave("Figure 2N.pdf", plot, width = 8, height = 8)
ggsave("Figure 2N.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "IL-17 signaling pathway", 
                  "hsa04657", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7A.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7A.tiff", plot, width = 8, height = 8)


plot <- gseaplot2(kk2,
                  title = "IL-17 signaling pathway", 
                  "hsa04657", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7B.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7B.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "ErbB signaling pathway", 
                  "hsa04012", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7C.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7C.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "ErbB signaling pathway", 
                  "hsa04012", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7D.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7D.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Wnt signaling pathway", 
                  "hsa04310", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7E.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7E.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Wnt signaling pathway", 
                  "hsa04310", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7F.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7F.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "cGMP-PKG signaling pathway", 
                  "hsa04022", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7G.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7G.tiff", plot, width = 8, height = 8)


plot <- gseaplot2(kk2,
                  title = "cGMP-PKG signaling pathway", 
                  "hsa04022", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7H.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7H.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Hippo signaling pathway - multiple species", 
                  "hsa04392", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7I.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7I.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Hippo signaling pathway - multiple species", 
                  "hsa04392", 
                  color="red", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("Supplementary Figure 7J.pdf", plot, width = 8, height = 8)
ggsave("Supplementary Figure 7J.tiff", plot, width = 8, height = 8)

