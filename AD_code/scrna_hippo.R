rm(list = ls())
#install_github("immunogenomics/harmony")
if(!require(BiocManager)) install.packages("BiocManager",update = F,ask = F)
if(!require(Seurat))BiocManager::install("Seurat")
if(!require(Matrix))BiocManager::install("Matrix")
if(!require(reticulate))BiocManager::install("reticulate")
if(!require(sctransform))BiocManager::install("sctransform")
if(!require(viridis))BiocManager::install("viridis")
if(!require(tidyr))BiocManager::install("tidyr")
if(!require(magrittr))BiocManager::install("magrittr")
if(!require(reshape2))BiocManager::install("reshape2")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(cowplot))BiocManager::install("cowplot")
if(!require(dplyr))BiocManager::install("dplyr")
if(!require(purrr))BiocManager::install("purrr")
if(!require(readr))BiocManager::install("readr")
if(!require(readxl))BiocManager::install("readxl")
if(!require(stringr))BiocManager::install("stringr")
library(devtools)
library(harmony)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(Matrix)
library(data.table)
workingDir = "D:/AD/HIP"
setwd(workingDir)
getwd()
options(scipen = 200)

# 读取数据---------------------------------------------
samples <- dir(path="D:/AD/HIP/", pattern="^GSM")
for (i in 1:length(samples)) {
  if (grepl("barcode.txt.gz", samples[i])) {
    # 如果文件名包含 "barcode.txt.gz"，读取条形码文件
    barcode_content <- fread(cmd = paste0("gzip -dc ", samples[i]))
    # 进行后续处理...
  } else if (grepl(".txt.gz", samples[i])) {
    # 如果文件名包含 ".txt.gz"，读取数据文件
    data_content <- fread(cmd = paste0("gzip -dc ", samples[i]))
    # 进行后续处理...
  }
}
length(unique(data_content$V1))
data_content[nrow(data_content),1:4]
hippo=data_content
hippo <- data.frame(hippo)
rownames(hippo) <- as.character(hippo$V1)   
exp=hippo[,2:ncol(hippo)]
dimnames=list(rownames(exp),colnames(exp))
hip=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
hip <- CreateSeuratObject(counts = hip,project = "Hippocampus", min.cells = 3, min.features = 200)#min.cells为基因存在样本最小数，需要根据实际情况选择，min.features = 50基因最小存在细胞数

#使用PercentageFeatureSet函数计算线粒体基因的百分比
nrow(hip)
saveRDS(hip, file = "Seurat_hip.RDS")
hip[["percent.mt"]] <- PercentageFeatureSet(object = hip, pattern = "^MT-")
hip[["percent.hb"]] <- PercentageFeatureSet(object = hip, pattern = "^HBA|^HBB")
hip[["percent.rp"]] <- PercentageFeatureSet(object = hip, pattern = "^RPS|^RPL")
head(hip@meta.data,5)
VlnPlot(object = hip, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
ggsave("VlnPlot.png", path = "D:/AD/HIP", width = 50, height = 30, units = "cm")
#涔嬪悗鍙互鐢荤浉鍏虫€ф暎鐐瑰浘锛?
plot1 <- FeatureScatter(hip,feature1 = "nCount_RNA",
                        feature2 = "percent.mt",
                        pt.size = 1.5)#mt鐧惧垎姣斿緢浣庯紝鎵€浠ュ彲浠ヤ笉鐢ㄧ敾杩欎釜

plot2 <- FeatureScatter(hip,feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        pt.size=1.5)

plot1+plot2 #瀹屾垚璐ㄦ帶锛屼笅闈㈠紑濮嬪彂鐜伴珮鍙樺紓鍩哄洜銆?
ggsave("Scatter.png", path = "D:/AD/HIP", width = 50, height = 30, units = "cm")

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
percent.mt_lower <- 0
percent.mt_upper <- 30
percent.hb_lower <- 0
percent.hb_upper <- 5

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(hip) {
  qc_data <- as_tibble(FetchData(hip, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb", "percent.rp")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = percent.mt))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = percent.hb))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = percent.rp))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), percent.mt, color = nFeature_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), percent.hb, color = nFeature_RNA))) + 
      geom_hline(yintercept = percent.hb_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.hb_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), percent.rp, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), percent.mt, color = nCount_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), percent.hb, color = nCount_RNA))) + 
      geom_hline(yintercept = percent.hb_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.hb_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), percent.rp, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(percent.rp, percent.mt, color = nCount_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(percent.rp, percent.mt, color = nFeature_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}

## Before filtering

hip_unfiltered <- hip

#saveRDS(hip_unfiltered, file = "seurat_objects/all_unfiltered.RDS")

qc_std_plot(hip_unfiltered)
ggsave("SuppFig1A.png", path = "D:/AD/HIP", width = 30, height = 30, units = "cm")


## After filtering

#hip_unfiltered <- readRDS("all_unfiltered.RDS")

hip <- subset(hip_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

qc_std_plot(hip)
ggsave("SuppFig1B.png", path = "D:/AD/HIP", width = 30, height = 30, units = "cm")

View(hip@meta.data)

#### Data normalization ####----------------------------------------------------

hip <- NormalizeData(hip, normalization.method = "LogNormalize",  scale.factor = 10000)
hip <- FindVariableFeatures(hip)
top10 <- head(VariableFeatures(hip), 10)
plot1 <- VariableFeaturePlot(hip,pt.size = 1.5)
plot2 <- LabelPoints(plot = plot1, points = top10,repel = TRUE,pt.size = 1.5)
plot1 + plot2
hip <- ScaleData(hip)
hip <- RunPCA(hip)
print(hip[["pca"]], dims = 1:5,  nfeatures = 5)
VizDimLoadings(object = hip, dims = 1:2, reduction = "pca", nfeatures = 20) #4涓狿C 20涓熀鍥?
DimPlot(hip, reduction = "pca")
DimHeatmap(hip, dims = 1:15, cells = 500, balanced = TRUE, nfeatures = 30, ncol=5)#淇敼dims =10璇曡瘯
hip <- JackStraw(object = hip, num.replicate = 100)
hip <- ScoreJackStraw(object = hip, dims = 1:20)
JackStrawPlot(object = hip, dims = 1:20, reduction = "pca")
ElbowPlot(hip, reduction="pca")#鏍规嵁E1bowP1ot纭畾PC鏁伴噺
JackStrawPlot(object = hip, dims = 1:20, reduction = "pca")+ElbowPlot(hip, reduction="pca")
saveRDS(hip,file = "hip_pca.RDS")

#hip <- hip_pca

#metadata-----------------------------------------------------------------------
# 设置目录路径
directory_path <- "D:/AD/HIP"

# 获取目录下的文件列表
file_list <- list.files(path = directory_path, pattern = "\\.barcode\\.txt\\.gz$", full.names = TRUE)

# 创建一个空的列表用于存储数据框
data_list <- list()

# 逐个读取文件并存储内容
#for (file in file_list) {
  # 读取文件内容
#  file_content <- read.delim(gzfile(file), header = FALSE, stringsAsFactors = FALSE)
  
  # 获取文件名中第一个下划线后面并且第二个下划线前面的内容作为样本名
#  file_name <- basename(file) # 获取文件名
#  sample_name <- gsub("^(?:[^_]*_){1}([^_]*)_(?:[^_]*_).*", "\\1", file_name, perl = TRUE) # 提取样本名
  
  # 获取第二个下划线后且在 ".barcode" 前面的内容作为 group 名
#  group_name <- gsub("^(?:[^_]*_){2}([^_]+)_.*\\.barcode.*", "\\1", file_name, perl = TRUE) # 提取 group 名
  
  # 将文件内容存储为数据框并添加到列表中，同时添加样本名和 group 列
#  file_content$sample <- sample_name
#  file_content$group <- group_name
  
#  data_list[[sample_name]] <- file_content}
#combined_data <- do.call(rbind, data_list)
# 逐个读取文件并存储内容
for (file in file_list) {
  # 读取文件内容
  file_content <- read.delim(gzfile(file), header = FALSE, stringsAsFactors = FALSE)
  
  # 获取文件名作为样本名
  file_name <- basename(file) # 获取文件名
  sample_name <- file_name # 提取样本名
  
  # 获取文件名作为作为 group 名
  group_name <-  file_name # 提取 group 名
  
  # 将文件内容存储为数据框并添加到列表中，同时添加样本名和 group 列
  file_content$sample <- sample_name
  file_content$group <- group_name
  
  data_list[[sample_name]] <- file_content
}
# 合并数据列表中的所有数据框
combined_data <- do.call(rbind, data_list)
#------------------------------------------------------------------------------
cli <- combined_data
names(cli)[1] <- "Barcode"
nrow(cli)
cli <- cli[!grepl("x",cli$Barcode),]
cli <- cli %>% distinct(Barcode, .keep_all = TRUE)
# 获取 meta 数据框的行名
row_names <- rownames(hip@meta.data)
nrow(hip@meta.data)
# 将行名作为新的列添加到 meta 数据框中
hip@meta.data$Barcode <- row_names
nrow(hip@meta.data)
# 调整列的顺序，将 Barcode 移到第一列位置
hip@meta.data <- hip@meta.data[, c("Barcode", names(hip@meta.data)[-which(names(hip@meta.data) == "Barcode")])]
metadata <- FetchData(hip, "Barcode")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = cli, by = "Barcode")
nrow(metadata)
rownames(metadata) <- metadata$cell_id

hip <- AddMetaData(hip, metadata = metadata)
hip@meta.data=hip@meta.data[,-1]
#hip@meta.data$orig.ident <- hip@meta.data$sample




# 将行名作为新的列添加到 meta 数据框中
hip@meta.data$group <- hip@meta.data$sample
nrow(hip@meta.data)
# 提取第一个下划线后面并且第二个下划线前面的内容作为 orig.ident 列的内容
hip@meta.data$sample <- gsub("^(?:[^_]*_){1}([^_]*)_(?:[^_]*_).*", "\\1", hip@meta.data$sample, perl = TRUE)

# 提取第二个下划线后面并且第三个下划线前面的内容作为 Sample 列的内容
hip@meta.data$group <- gsub("^[^_]+_([^_]+)_([^_]+)_.*", "\\2", hip@meta.data$group, perl = TRUE)
hip@meta.data$orig.ident <- hip@meta.data$sample


#-----------------------------------------------------------------------

hip <- FindNeighbors(hip, reduction = "pca",dims=1:20)
hip <- FindClusters(hip, resolution = 0.4)


#tSNE闄嶇淮
hip <- RunTSNE(object = hip, dims = 1:12, reduction = "pca")
print(DimPlot(hip, reduction = "tsne", label = TRUE, pt.size = 0.1)) + labs(title = "Resolution: 0.4")

#方法一：cellmaker法
gene_marker <- c("AQP4","GFAP",    #AST
                 "LHFPL3", #OPC
                 "MBP", "MOBP",      #OLIGO 
                 "CX3CR1","APBB1IP", #MG
                 "SLC17A7",          #Glu_N
                 "GAD2",             #GABA_N
                 "RELN",             #Cajal-Retzius cells
                 "CLIC6",            #Choroid pleus cells
                 "FLT1",             #Endothelial cells
                 "FOLR1"             #Dopa_N  
                 )            

celltype=data.frame(clusterID=0:11,
                    celltype='unknown')

P<-DotPlot(hip, features =c("AQP4","GFAP","MBP","MOBP","SLC17A7","LHFPL3","APBB1IP","FOLR1","GAD2"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
  
p_AST <- DotPlot(hip, features =c("AQP4","GFAP"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(0,6,11),2]="Astro"


p_OLIGO <- DotPlot(hip, features =c("MBP","MOBP"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(1),2]="Oligo"

p_Glu <- DotPlot(hip, features =c("SLC17A7"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(2,4,5),2]="Glu_N"

p_OPC <- DotPlot(hip, features =c("LHFPL3"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(3),2]="OPCs"


p_MG <- DotPlot(hip, features =c("APBB1IP"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(7),2]="Micro"

p_DPN <- DotPlot(hip, features =c("FOLR1"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(9),2]="Dopa_N"

p_GABA <- DotPlot(hip, features =c("GAD2"), group.by = "RNA_snn_res.0.4") + 
  coord_flip() + 
  scale_color_viridis()
celltype[celltype$clusterID %in% c(8,10),2]="GABA_N"
hip_1 <- hip
hip_1@meta.data$celltype="NA"
for (i in 1:nrow(celltype)) {hip_1@meta.data[which(hip_1@meta.data$seurat_clusters == celltype$clusterID[i]),'celltype'] <- celltype$celltype[i]}

Idents(hip_1) <- hip_1@meta.data$celltype
unique(Idents(hip_1))
Glut_N <-  subset(hip_1, idents = "Glu_N")
OPC <- subset(hip_1, idents = "OPCs")
Astro <- subset(hip_1, idents = "Astro")
GABA_N <- subset(hip_1, idents = "GABA_N")
Oligo <- subset(hip_1, idents = "Oligo")
Micro <- subset(hip_1, idents = "Micro")
Dopa_N <- subset(hip_1, idents = "Dopa_N")

Glut_N    <- ScaleData(Glut_N)
OPC       <- ScaleData(OPC)
Astro     <- ScaleData(Astro)
GABA_N    <- ScaleData(GABA_N)
Oligo     <- ScaleData(Oligo)
Micro     <- ScaleData(Micro)
Dopa_N    <- ScaleData(Dopa_N)

hip <- hip_1
# Cell cycle scoring------------------------------------------------------------

### add cell cycle, cc.genes loaded with Seurat

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(hip) {
  hip <- CellCycleScoring(hip, s.genes, g2m.genes)
  hip@meta.data$CC.Diff <- hip@meta.data$S.Score - hip@meta.data$G2M.Score
  return(hip)
}

hip <- score_cc(hip)

FeatureScatter(hip, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)

saveRDS(hip, file = "hip_all.RDS")

#hip <- subset(hip_all, subset = !is.na(group))
hip <- hip_all
#color scheme
unique(hip$celltype)
unique(hip$group)
unique(hip$orig.ident)
library(colorspace)
use_colors <- c(
  AD = "red",
  Control = "green",
  Aging="greenyellow",
  Adult="gray")


####细胞比例####
##################------------未统计分析（按 group 先合并 sample 再算比例）------------##################
library(reshape2)
library(ggplot2)
#define the color
library(ggsci)
cors <- pal_igv()(3) #定义颜色
##################------------统计分析------------##################
hip_subset <- subset(hip, subset = group %in% c("AD", "Aging", "Adult"))
new_order <- c("AD", "Aging", "Adult")
hip_subset$group <- factor(hip_subset$group,levels = new_order, ordered = T)
#准备绘图数据
prop_df1 <- table(hip_subset$celltype,hip_subset$group,hip_subset$sample) %>% reshape2::melt()
colnames(prop_df1) <- c("Cluster","Group","Sample","Number")
prop_df1$Cluster <- factor(prop_df1$Cluster)
prop_df1$Group   <- factor(prop_df1$Group)
library(dplyr)
# 按 sample 计算细胞比例
prop_df1 <- prop_df1 %>% group_by(Group, Sample) %>% filter(sum(Number) > 0) %>% mutate(Proportion = Number / sum(Number)) %>% ungroup()
# 检查：比例是否加和为 1
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_prop = sum(Proportion))
# 检查：是否存在细胞数为 0 的样本
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_cells = sum(Number), .groups = "drop") %>% filter(total_cells == 0)
write.csv(prop_df1, file = "Celltype_Proportion_by_sample2.csv", row.names = FALSE)

#比例图1  尺寸6*3
library(ggplot2)
ggplot(prop_df1, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_manual(values = c("#8DD3C7", "#B3DE69", "grey", "#F490A9", "#FDB499", "skyblue1", "#84B1ED")) +
  theme_classic() +
  labs(x = "Sample", y = "Cell proportion", fill = "Cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file="proportion_barplot_sample.pdf", width=6, height=4)
ggsave(file="proportion_barplott_sample.png", width=6, height=4)
# 每个 group 一个柱（先对 sample 求平均）
prop_group <- prop_df1 %>% group_by(Group, Cluster) %>% summarise(mean_prop = mean(Proportion), .groups = "drop")
ggplot(prop_group,aes(x = mean_prop, y = Group, fill = Cluster)) +
       geom_bar(stat = "identity", width = 0.8, position = "fill") +
       scale_fill_manual(values = c("#8DD3C7", "#B3DE69", "grey", "#F490A9", "#FDB499", "skyblue1", "#84B1ED")) +
       theme_bw() +
       theme(panel.grid = element_blank()) +
       labs(x = "", y = "Cell proportion", fill = "Cell type") +
       theme(
         axis.title.y = element_text(size = 14, colour = "black"),
         axis.title.x = element_text(size = 14, colour = "black"),
         axis.text.y  = element_text(size = 12, colour = "black"),
         axis.text.x  = element_text(size = 12, colour = "black"),
         axis.text.x.bottom = element_text( hjust = 1, vjust = 1, size = 14),
         legend.text  = element_text(size = 12),
         legend.title = element_text(size = 14))
ggsave(file="proportion_barplot_sig.pdf", width=6, height=3)
ggsave(file="proportion_barplot_sig.png", width=6, height=3)


hip_subset <- subset(hip, subset = group %in% c("AD", "Aging", "Adult"))
new_order <- c("Adult", "Aging","AD" )
hip_subset$group <- factor(hip_subset$group,levels = new_order, ordered = T)
#准备绘图数据
prop_df1 <- table(hip_subset$celltype,hip_subset$group,hip_subset$sample) %>% reshape2::melt()
colnames(prop_df1) <- c("Cluster","Group","Sample","Number")
prop_df1$Cluster <- factor(prop_df1$Cluster)
prop_df1$Group   <- factor(prop_df1$Group)
library(dplyr)
# 按 sample 计算细胞比例
prop_df1 <- prop_df1 %>% group_by(Group, Sample) %>% filter(sum(Number) > 0) %>% mutate(Proportion = Number / sum(Number)) %>% ungroup()
# 检查：比例是否加和为 1
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_prop = sum(Proportion))
# 检查：是否存在细胞数为 0 的样本
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_cells = sum(Number), .groups = "drop") %>% filter(total_cells == 0)
write.csv(prop_df1, file = "Celltype_Proportion_by_sample2.csv", row.names = FALSE)
# 先手动计算每个 panel 的 Kruskal p 值（并生成 ns）
library(dplyr)
kw_label <- prop_df1 %>% group_by(Cluster) %>%
  summarise(p = kruskal.test(Proportion ~ Group)$p.value, .groups = "drop") %>%
#  mutate(label = paste0("Kruskal p = ", formatC(p, format = "e", digits = 2)) )
  mutate(label = paste0("Kruskal P = ", signif(p, 3)))
 # 或用 signif(p, 3)
### 绘图---------------
library(ggplot2)
library(ggpubr)
library(dplyr)
kw_label <- prop_df1 %>% group_by(Cluster) %>%
  summarise(p = kruskal.test(Proportion ~ Group)$p.value,
            y_pos = max(Proportion, na.rm = TRUE) * 1.05) %>%
  mutate(label = paste0("Kruskal p = ", signif(p, 3)))
library(ggplot2)
library(ggpubr)
ggplot(prop_df1, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8) +
  facet_wrap(~ Cluster, scales = "free_y") +
  ## Kruskal–Wallis p（右上角，红色，加粗）
  geom_text(data = kw_label,
    aes(x = Inf, y = y_pos, label = label),
    inherit.aes = FALSE, color = "red",
    fontface = "bold", size = 3, hjust = 1.1) +
    ## pairwise（可选，和总检验逻辑不冲突）
    stat_compare_means(method = "wilcox.test",
    label = "p.signif", hide.ns = TRUE) +
    labs(x = "Group", y = "Cell proportion") +
    theme_classic()
ggsave(file="Moyamoya_proportion_barplot_2_sig.pdf", width=9, height=6)
ggsave(file="Moyamoya_proportion_barplot_2_sig.png", width=9, height=6)

 
features = gene_marker
length(features)
#绘图
library(ggsci)
cors <- c("#8DD3C7","#B3DE69","grey","#F490A9","#FDB499","skyblue1","#84B1ED") #定义颜色
pdf(file=paste0('hip_subset_Celltype_vlnplot.pdf'),width = 6,height = 8)
VlnPlot(hip, features = features, 
        stack = TRUE, 
        sort = TRUE, 
        cols = cors,
        split.by =  "celltype" , #每种cluster 一个颜色
        flip = TRUE,
        split.plot = TRUE) +
  theme(legend.position = "none")
dev.off()

##绘制细胞marker气泡图
pdf(file=paste0('hip_Celltype_dotplot.pdf'),width = 7,height = 6)
DotPlot(hip_subset, features = features, assay='RNA') + 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust = 1))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#3C5488FF','#4DBBD5FF','#F39B7FFF','#E64B35FF')) #颜色
dev.off()

#5 差异表达分析----------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(devtools)
#devtools::install_github('immunogenomics/presto')
library(presto)
hip_1 <- hip
#hip_1@meta.data <- hip_1@meta.data[complete.cases(hip_1@meta.data), ]
library(Seurat)
# Filter cells based on desired groups
desired_groups <- c("AD", "Aging", "Adult")
hip_filtered <- hip_1[, hip_1$group %in% desired_groups]
# Reorder group factor levels
new_order <- c("AD", "Aging", "Adult")
hip_filtered$group <- factor(hip_filtered$group, levels = new_order, ordered = TRUE)
# Plot DimPlot
DimPlot(hip_filtered, split.by = 'group', pt.size = 0.8) + 
  theme(text = element_text(size = 14))
hip_1$celltype.group <- paste(hip_1$celltype, hip_1$group, sep = "_")
hip_1$celltype <- Idents(hip_1)
Idents(hip_1) <- "celltype.group"
#saveRDS(hip_1, file = "hip_1.RDS")

hip_1_clean <- na.omit(hip_1)
#OPCs------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'OPCs_AD',ident.2 = 'OPCs_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"OPCs_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'OPCs_Aging',ident.2 = 'OPCs_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"OPCs_Aging_Ault.CSV")

#Astro---------------------------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'Astro_AD',ident.2 = 'Astro_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Astro_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'Astro_Aging',ident.2 = 'Astro_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Astro_Aging_Ault.CSV")
#
#GABA_N---------------------------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'GABA_N_AD',ident.2 = 'GABA_N_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"GABA_N_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'GABA_N_Aging',ident.2 = 'GABA_N_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"GABA_N_Aging_Ault.CSV")

#Oligo---------------------------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'Oligo_AD',ident.2 = 'Oligo_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Oligo_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'Oligo_Aging',ident.2 = 'Oligo_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Oligo_Aging_Ault.CSV")

#Micro----------------------------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'Micro_AD',ident.2 = 'Micro_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Micro_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'Micro_Aging',ident.2 = 'Micro_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Micro_Aging_Ault.CSV")

#Glu_N----------------------------------------------------------------------
mydeg <- FindMarkers(hip_1,ident.1 = 'Glu_N_AD',ident.2 = 'Glu_N_Aging', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Glu_N_AD_Aging.CSV")

mydeg <- FindMarkers(hip_1,ident.1 = 'Glu_N_Aging',ident.2 = 'Glu_N_Adult', verbose = TRUE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
write.csv(mydeg,"Glu_N_Aging_Ault.CSV")

cellfordeg<-levels(hip_1$celltype)
cellfordeg[1]
         
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(hip_1, ident.1 = paste0(cellfordeg[i],"_Aging"), ident.2 = paste0(cellfordeg[i],"_Adult"), verbose = TRUE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".csv"))
}
