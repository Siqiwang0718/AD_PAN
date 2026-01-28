
# Load R packages
library(ggVennDiagram)
library(VennDiagram)
library(RColorBrewer)
library(data.table)
library(ggplot2)
data <- Veen_snRNA_PRGs_rawdata
######------Figure 2G and 2H------######
# Create a collection and extract the unique values from the PANoptosis_re, Up, and Down columns.
# Check if the PANoptosis_re column contains an empty string or other null values.
PRGs <- unique(data$PRGs[!is.na(data$PRGs) & data$PRGs != ""])
ADvsAg.Glu <- unique(data$ADvsAg.Glu[!is.na(data$ADvsAg.Glu) & data$ADvsAg.Glu != ""])
ADvsAg.GABA <- unique(data$ADvsAg.GABA[!is.na(data$ADvsAg.GABA) & data$ADvsAg.GABA != ""])
Ag.vsAd.Glu <- unique(data$Ag.vsAd.Glu[!is.na(data$Ag.vsAd.Glu) & data$Ag.vsAd.Glu != ""])
Ag.vsAd.GABA <- unique(data$Ag.vsAd.GABA[!is.na(data$Ag.vsAd.GABA) & data$Ag.vsAd.GABA != ""])
# Create a data list
data.list1 <- list(PRGs=PRGs, ADvsAg.Glu=ADvsAg.Glu, ADvsAg.GABA=ADvsAg.GABA)
# Draw Venn diagram
p <- ggVennDiagram(data.list1, 
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

# Save as PDF file
ggsave("Figure 2G.pdf", plot = p, width = 8, height = 8)

# Create a data list
data.list2 <- list(PRGs=PRGs, Ag.vsAd.Glu=Ag.vsAd.Glu, Ag.vsAd.GABA=Ag.vsAd.GABA)
# Draw Venn diagram
p <- ggVennDiagram(data.list2, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE, 
                   set_color = c("red", "yellowgreen", "blue"),
                   label_size = 12,
                   set_size = 8) +  
  #scale_fill_distiller(palette = "RdBu") + 
  scale_fill_gradient(low = "white", high = "white") +  
  #scale_color_brewer(palette = "Set1")
  theme_void() +   
  theme(
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
    
  )
# Save as PDF file
ggsave("Figure 2H.pdf", plot = p, width = 8, height = 8)

######------Figure 2I------######
data <- Veen_HIP_Neuron_rawdata
ADvsAg.HIP <- unique(as.character(data$ADvsAg.HIP[!is.na(data$ADvsAg.HIP) & data$ADvsAg.HIP != ""]))
ADvsAg.Neuron <- unique(as.character(data$ADvsAg.Neuron[!is.na(data$ADvsAg.Neuron) & data$ADvsAg.Neuron != ""]))
Ag.vsAd.HIP <- unique(as.character(data$Ag.vsAd.HIP[!is.na(data$Ag.vsAd.HIP) & data$Ag.vsAd.HIP != ""]))
Ag.vsAd.Neuron <- unique(as.character(data$Ag.vsAd.Neuron[!is.na(data$Ag.vsAd.Neuron) & data$Ag.vsAd.Neuron != ""]))

# Create a data list
data.list1 <- list(Ag.vsAd.HIP=Ag.vsAd.HIP,
                   ADvsAg.HIP=ADvsAg.HIP, 
                   ADvsAg.Neuron=ADvsAg.Neuron, 
                   Ag.vsAd.Neurons=Ag.vsAd.Neuron)
# Draw Venn diagram
p <- ggVennDiagram(data.list1, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE, 
                   set_color = c("red", "yellowgreen", "blue","skyblue"),
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
ggsave("Figure 2I.pdf", plot = p, width = 6, height = 6)
ggsave("Figure 2I.png", plot = p, width = 6, height = 6, dpi = 300)


