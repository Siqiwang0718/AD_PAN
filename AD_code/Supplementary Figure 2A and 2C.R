# Load R packages
library(ggVennDiagram)
library(VennDiagram)
library(RColorBrewer)
library(data.table)
library(ggplot2)
######------Supplementary Figure 2A------######
# Create a collection and extract the unique values from the PANoptosis_re, Up, and Down columns.
# Check if the PANoptosis_re column contains an empty string or other null values.
data <- Veen_AD_Aging_rawdata
PAN <- unique(data$PANoptosis_re[!is.na(data$PANoptosis_re) & data$PANoptosis_re != ""])
AD_Aging_Up <- unique(data$up[!is.na(data$up) & data$up != ""])
AD_Aging_Down <- unique(data$down[!is.na(data$down) & data$down != ""])
# Create a data list
data.list <- list(PAN=PAN, AD_Aging_Up=AD_Aging_Up, AD_Aging_Down=AD_Aging_Down)
# Draw Venn diagram
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE, 
                   set_color = c("red", "yellowgreen", "cyan2"),
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
ggsave("Supplementary Figure 2A.pdf", plot = p, width = 8, height = 8)

######------Supplementary Figure 2C------######
# Create a collection and extract the unique values from the PANoptosis_re, Up, and Down columns.
# Check if the PANoptosis_re column contains an empty string or other null values.
data <- Veen_Aging_Adult_rawdata
PAN <- unique(data$PANoptosis_re[!is.na(data$PANoptosis_re) & data$PANoptosis_re != ""])
Aging_Adult_Up <- unique(data$up[!is.na(data$up) & data$up != ""])
Aging_Adult_Down <- unique(data$down[!is.na(data$down) & data$down != ""])


# Create a data list
data.list <- list(PAN=PAN, Aging_Adult_Up=Aging_Adult_Up, Aging_Adult_Down=Aging_Adult_Down)
# Draw Venn diagram
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE,  # 不显示百分比
                   set_color = c("red", "yellowgreen", "cyan2"),
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
ggsave("Supplementary Figure 2C.pdf", plot = p, width = 8, height = 8)

















