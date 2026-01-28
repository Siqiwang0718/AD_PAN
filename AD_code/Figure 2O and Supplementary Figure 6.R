# Load R packages
library(VennDiagram)
library(RColorBrewer)
library(data.table)
library(ggVennDiagram)
library(ggplot2)
######------Supplementary Figure 6------######
# Create a collection and extract the unique values from the PANoptosis_re, Up, and Down columns.
# Check if the PANoptosis_re column contains an empty string or other null values.
data <- Veen_ssGSEA_TP63_TUBA1B_rawdata
TP63_Enrich <- unique(data$TP63_Enrich[!is.na(data$TP63_Enrich) & data$TP63_Enrich != ""])
TUBA1B_Enrich <- unique(data$TUBA1B_Enrich[!is.na(data$TUBA1B_Enrich) & data$TUBA1B_Enrich != ""])
# Create a data list
data.list <- list(TP63_Enrich=TP63_Enrich, TUBA1B_Enrich=TUBA1B_Enrich)
library(ggVennDiagram)
library(RColorBrewer)
# Draw Venn diagram
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   show_percentage = FALSE,  
                   set_color = c("red", "yellowgreen"),
                   label_size = 12,
                   set_size = 8) +  
  scale_fill_gradient(low = "white", high = "white") +  
  theme_void() +  
  theme(
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# Save as PDF file
ggsave("Supplementary Figure 6.pdf", plot = p, width = 8, height = 8)


######------Figure 2O------######
# Create a collection and extract the unique values from the PANoptosis_re, Up, and Down columns.
# Check if the PANoptosis_re column contains an empty string or other null values.
data <- Veen_ssGSEA_PRGs_rawdata
DE.PRGs <- unique(data$DE.PRGs[!is.na(data$DE.PRGs) & data$DE.PRGs != ""])
TP63_HIF.1 <- unique(data$TP63_HIF.1[!is.na(data$TP63_HIF.1) & data$TP63_HIF.1 != ""])
TUBA1B_HIF.1 <- unique(data$TUBA1B_HIF.1[!is.na(data$TUBA1B_HIF.1) & data$TUBA1B_HIF.1 != ""])
# 创建数据列表
data.list <- list(DE.PRGs=DE.PRGs, TP63_HIF.1=TP63_HIF.1, TUBA1B_HIF.1=TUBA1B_HIF.1)
# drawing
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE, 
                   set_color = c("red", "yellowgreen", "blue"),
                   label_size = 12,
                   set_size = 8) +  
  scale_fill_gradient(low = "white", high = "white") + 
  theme_void() +   
  theme(
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
    
  )
print(p)
# Save as PDF file
ggsave("Figure 2O.pdf", plot = p, width = 8, height = 8)
