
# Load R package
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library("stringr")
# Set color
color_1 <- c("green","red")
######------Supplementary Figure 2B------######
p <- ggplot(GO_3, aes(y = GeneRatio, x = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = pvalue)) +
  labs(x = "GO Description", y = "GeneRatio") +
  labs(title = "KEGG signaling pathways (AD vs. Aging)") +
  coord_flip() +
  theme_bw() +
  scale_color_gradient(low = color_1[2], high = color_1[1], labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::number_format(scale = 0.1, decimal.mark = ".", big.mark = ","))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 60))+
  theme(plot.title = element_text(hjust = 0.5,size=12,color = "black", family = "sans",face = "bold"),
        strip.text.y = element_text(size =12,color = "black", family = "sans"), 
        legend.position="right", 
        legend.title = element_text(size=12,color = "black", family = "sans"), 
        legend.text = element_text(size=12,color = "black", family = "sans"), 
        axis.text.x = element_text(size=12,color = "black", family = "sans"), 
        axis.text.y = element_text(size=12,color = "black", family = "sans"),
        axis.title.x = element_text(size=12,color = "black", family = "sans"), 
        axis.title.y = element_text(size=12,color = "black", family = "sans"))
print(p)
ggsave(filename = "Supplementary Figure 2B.pdf", plot = p, width = 6, height = 5, dpi = 600)
ggsave(filename = "Supplementary Figure 2B.png", plot = p, width = 6, height = 5, dpi = 600)


######------Supplementary Figure 2D------######
p <- ggplot(GO_3, aes(y = GeneRatio, x = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = pvalue)) +
  labs(x = "GO Description", y = "GeneRatio") +
  labs(title = "KEGG signaling pathways (Aging vs. Adult)") +
  coord_flip() +
  theme_bw() +
  scale_color_gradient(low = color_1[2], high = color_1[1], labels = scales::scientific_format()) +
  scale_y_continuous(labels = scales::number_format(scale = 0.1, decimal.mark = ".", big.mark = ","))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 60))+
  theme(plot.title = element_text(hjust = 0.5,size=12,color = "black", family = "sans",face = "bold"),
        strip.text.y = element_text(size =12,color = "black", family = "sans"), 
        legend.position="right", 
        legend.title = element_text(size=12,color = "black", family = "sans"), 
        legend.text = element_text(size=12,color = "black", family = "sans"), 
        axis.text.x = element_text(size=12,color = "black", family = "sans"), 
        axis.text.y = element_text(size=12,color = "black", family = "sans"),
        axis.title.x = element_text(size=12,color = "black", family = "sans"), 
        axis.title.y = element_text(size=12,color = "black", family = "sans"))
print(p)
ggsave(filename = "Supplementary Figure 2D.pdf", plot = p, width = 6, height = 5, dpi = 600)
ggsave(filename = "Supplementary Figure 2D.png", plot = p, width = 6, height = 5, dpi = 600)
