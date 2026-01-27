rm(list = ls())
workingDir = "E:/AD_hippoDatabase/GSE48350"
setwd(workingDir)
library(limma)
library(GEOquery)
library(Biobase)
library(Biobase)
library(GEOquery)
library(dplyr)
library(ggplot2)
library(ggrepel)

# 从GEO系列矩阵文件中读取数据
f = "GSE48350"
gset <- getGEO(f, 
               destdir = ".", 
               AnnotGPL = TRUE, 
               getGPL = TRUE)

exp<-exprs(gset[[1]])
cli<-pData(gset[[1]])
GPL<-fData(gset[[1]])
q <- quantile(exp, c(0.01, 0.5, 0.99), na.rm = TRUE)

if (q[3] > 20) {
  message("Raw expression detected. Applying log2.")
  exp <- log2(exp + 1)
} else {
  message("Expression already log-transformed. No further transformation applied.")
}
max(exp)

# 安装并加载dplyr包
library(dplyr)

# 将age (yrs):ch1列转换为数值类型
cli$age_numeric <- as.numeric(cli$`age (yrs):ch1`)

group_list <- sapply(1:nrow(cli), function(i) {
  if (grepl("postcentral gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] != "") {
    "PCG_AD"
  } else if (grepl("postcentral gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] == "") {
    ifelse(cli[i, "age_numeric"] < 60, "PCG_Adult", "PCG_Aging")
  } else if (grepl("post-central gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] != "") {
    "PCG_AD"
  } else if (grepl("post-central gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] == "") {
    ifelse(cli[i, "age_numeric"] < 60, "PCG_Adult", "PCG_Aging")
  } else if (grepl("superior frontal gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] != "") {
    "SFG_AD"
  } else if (grepl("superior frontal gyrus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] == "") {
    ifelse(cli[i, "age_numeric"] < 60, "SFG_Adult", "SFG_Aging")
  } else if (grepl("hippocampus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] != "") {
    "HIP_AD"
  } else if (grepl("hippocampus", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] == "") {
    ifelse(cli[i, "age_numeric"] < 60, "HIP_Adult", "HIP_Aging")
  } else if (grepl("entorhinal cortex", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] != "") {
    "EntorCortex_AD"
  } else if (grepl("entorhinal cortex", cli[i, "brain region:ch1"]) & cli[i, "characteristics_ch1.4"] == "") {
    ifelse(cli[i, "age_numeric"] < 60, "EntorCortex_Adult", "EntorCortex_Aging")
  } else {
    "NA"
  }
})

table(group_list)
#cli$group <- group_list
write.csv(cli,"cli.csv")
###探针注释
library(hgu133plus2.db)
#BiocManager::install("AnnoProbe")
# 使用 toTable 函数获取注释信息
#ids_570<-toTable(hgu133plus2SYMBOL)
ids_570 <- AnnoProbe::idmap('GPL570') # 使用 AnnoProbe 包的 idmap 函数加载 GPL570 平台的注释数据。
load("E:/AD_hippoDatabase/GSE48350/GPL570_bioc.rda")         # 加载存储在路径 "D:/AD/GPL570_bioc.rda" 的下载好的 .rda 文件。
str(GPL570_bioc)                      # 使用 str() 函数显示加载的 GPL570_bioc 对象中数据的结构。
#exprs<-trans_array(exprs,ids_570)  #将 exprs 中的表达数据转换为匹配 ids_570 提供的注释信息，这些信息可能对应于 GPL570 平台。
row<-as.data.frame (rownames (exp))
colnames(row)<-c("probe_id")
expt2<-cbind(row,exp)
expt3<-merge(ids_570,expt2,by="probe_id")
expt3[1:4,1:6]
exp_symbol<-na.omit (expt3) # na.omit() 函数用于去除数据中包含缺失值的行，并将处理后的数据赋值给 exp_symbol 变量。
dim(exp_symbol)
exp_symbol[1:3,1:4]
table(duplicated(exp_symbol[,2]))
library(limma)
exp_unique<-avereps (exp_symbol[,-1],ID=exp_symbol[,2])
dim(exp_unique)
exp_unique[1:3,1:4]
max(exp_unique)

write.table(exp_unique,"E:/AD_hippoDatabase/GSE48350/GSE48350_exp.csv",sep=",",col.names=T, row.names=F,quote=F)

exp_1<-read.csv("E:/AD_hippoDatabase/GSE48350/GSE48350_exp.csv",header=T,row.names=1)
#创建一个分组的矩阵
condition<-factor(group_list,levels = c("PCG_Aging","SFG_Aging","HIP_Aging","EntorCortex_Aging",
                                        "PCG_Adult","SFG_Adult","HIP_Adult","EntorCortex_Adult",
                                        "PCG_AD","SFG_AD","HIP_AD","EntorCortex_AD"),ordered = F)
unique(group_list)
table(condition)
design=model.matrix(~0+condition)
colnames(design)=levels(factor(condition))
ncol(exp_1)
nrow(design)
rownames(design)=colnames(exp_1)
head(design)

#HIP_Aging versus HIP_Adult-------------------------------------------------------
cont.matrix<-makeContrasts('`HIP_Aging`-`HIP_Adult`',levels = design)
#第一步lmfit,其为每一个基因给定一系列的阵列来拟合线性模型
fit<-lmFit(exp_1,design)
#第2步ebayes,其给出一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算
fit2=contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
options(digits = 4)
#topTable(fit2,coef=2,adjust='BH')
tempOutput<-topTable(fit2,coef = 1,n=Inf)
tempOutput<-na.omit(tempOutput)
head(tempOutput)
#tempOutput["MEFV",]
#write.csv(tempOutput,"`HIP_Aging`-`HIP_Adult`.csv")

# 创建新列表示差异表达基因
tempOutput$DEG <- ifelse(tempOutput$adj.P.Val < 0.05 & tempOutput$logFC > 0, "Up",
                         ifelse(tempOutput$adj.P.Val < 0.05 & tempOutput$logFC < 0, "Down", "Stable"))

# 查看tempOutput数据框的列名
names(tempOutput)
# 添加行名为第一列
tempOutput$Gene <- rownames(tempOutput)

# 调整列的顺序，将Gene列放在第一列
tempOutput <- tempOutput[, c("Gene", names(tempOutput)[-which(names(tempOutput) == "Gene")])]
HIP_Aging_Adult <- tempOutput

# 导出为CSV文件
write.csv(HIP_Aging_Adult, file = "HIP_Aging_Adult.csv", row.names = FALSE)

#HIP_AD versus HIP_Aging----------------------------------------------------------
cont.matrix<-makeContrasts('`HIP_AD`-`HIP_Aging`',levels = design)
#第一步lmfit,其为每一个基因给定一系列的阵列来拟合线性模型
fit<-lmFit(exp_1,design)
#第2步ebayes,其给出一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算
fit2=contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
options(digits = 4)
#topTable(fit2,coef=2,adjust='BH')
tempOutput<-topTable(fit2,coef = 1,n=Inf)
tempOutput<-na.omit(tempOutput)
head(tempOutput)
#tempOutput["MEFV",]
#write.csv(tempOutput,"`HIP_AD`-`HIP_Aging`.csv")

# 创建新列表示差异表达基因
tempOutput$DEG <- ifelse(tempOutput$adj.P.Val < 0.05 & tempOutput$logFC > 0, "Up",
                         ifelse(tempOutput$adj.P.Val < 0.05 & tempOutput$logFC < 0, "Down", "Stable"))

# 查看tempOutput数据框的列名
names(tempOutput)
# 添加行名为第一列
tempOutput$Gene <- rownames(tempOutput)

# 调整列的顺序，将Gene列放在第一列
tempOutput <- tempOutput[, c("Gene", names(tempOutput)[-which(names(tempOutput) == "Gene")])]
HIP_AD_Aging <- tempOutput

# 导出为CSV文件
write.csv(HIP_AD_Aging, file = "HIP_AD_Aging.csv", row.names = FALSE)


library(VennDiagram)
library(RColorBrewer)
data <- readRDS("C:/Users/sqyql/Desktop/AD_PDF/Micro-array/Veen_AD_Aging_rawdata.RDS")
# 创建集合，提取 PANoptosis_re、Up 和 Down 列的唯一值
# 检查 PANoptosis_re 列中是否存在空字符串或其他类型的空值
PAN <- unique(data$PANoptosis_re[!is.na(data$PANoptosis_re) & data$PANoptosis_re != ""])
AD_Aging_Up <- unique(data$up[!is.na(data$up) & data$up != ""])
AD_Aging_Down <- unique(data$down[!is.na(data$down) & data$down != ""])

library(VennDiagram)#R包加载
# 创建数据列表
data.list <- list(PAN=PAN, AD_Aging_Up=AD_Aging_Up, AD_Aging_Down=AD_Aging_Down)
# 绘制 Venn 图，设置每个集合的颜色
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE,  # 不显示百分比
                   set_color = c("red", "yellowgreen", "cyan2"),
                   label_size = 12,
                   set_size = 8) +  
  #scale_fill_distiller(palette = "RdBu") + 
  scale_fill_gradient(low = "white", high = "white") +  # 使用渐变效果
  #scale_color_brewer(palette = "Set1")
  theme_void() +  # 使用完全空白的主题
  theme(
    legend.position = "none",  # 隐藏图例
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
    
  )
print(p)
# 保存为PDF
ggsave("Veen_AD_Aging_inter.pdf", plot = p, width = 8, height = 8)

# 确保为每个元素命名
library(data.table)
names(data.list) <- c("PAN", "AD_Aging_Up", "AD_Aging_Down")
inter <- get.venn.partitions(data.list)
for(i in 1:nrow(inter)) 
  inter[i,'values'] <- paste(inter[[i,'..values..']],collapse=',')
# 创建文件并写入 BOM
fileConn <- file("Veen_AD_Aging_inter.csv", open = "w", encoding = "UTF-8")
writeLines("\uFEFF", fileConn)  # 添加 BOM
write.csv(inter[-5], fileConn, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
close(fileConn)


data <- readRDS("C:/Users/sqyql/Desktop/AD_PDF/Micro-array/Veen_Aging_Adult_rawdata.RDS")
# 创建集合，提取 PANoptosis_re、Up 和 Down 列的唯一值
# 检查 PANoptosis_re 列中是否存在空字符串或其他类型的空值
PAN <- unique(data$PANoptosis_re[!is.na(data$PANoptosis_re) & data$PANoptosis_re != ""])
Aging_Adult_Up <- unique(data$up[!is.na(data$up) & data$up != ""])
Aging_Adult_Down <- unique(data$down[!is.na(data$down) & data$down != ""])

library(VennDiagram)#R包加载
# 创建数据列表
data.list <- list(PAN=PAN, Aging_Adult_Up=Aging_Adult_Up, Aging_Adult_Down=Aging_Adult_Down)
# 绘制 Venn 图，设置每个集合的颜色
p <- ggVennDiagram(data.list, 
                   label_alpha = 0, 
                   label = c("count"),
                   label_color = "black",
                   #set_color = "black",
                   show_percentage = FALSE,  # 不显示百分比
                   set_color = c("red", "yellowgreen", "cyan2"),
                   label_size = 12,
                   set_size = 8) +  
  #scale_fill_distiller(palette = "RdBu") + 
  scale_fill_gradient(low = "white", high = "white") +  # 使用渐变效果
  #scale_color_brewer(palette = "Set1")
  theme_void() +  # 使用完全空白的主题
  theme(
    legend.position = "none",  # 隐藏图例
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
    
  )
print(p)
# 保存为PDF
ggsave("Veen_Aging_Adult_inter.pdf", plot = p, width = 8, height = 8)

# 确保为每个元素命名
names(data.list) <- c("PAN", "Aging_Adult_Up", "Aging_Adult_Down")
inter <- get.venn.partitions(data.list)
for(i in 1:nrow(inter)) 
  inter[i,'values'] <- paste(inter[[i,'..values..']],collapse=',')
# 创建文件并写入 BOM
fileConn <- file("Veen_Aging_Adult_inter.csv", open = "w", encoding = "UTF-8")
writeLines("\uFEFF", fileConn)  # 添加 BOM
write.csv(inter[-5], fileConn, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
close(fileConn)
