rm(list = ls())
#rt <- readRDS("C:/Users/sqyql/Desktop/AD_PDF/Micro-array/ssGSEA/AD_Aging_count.RDS")
rt=read.table("C:/Users/sqyql/Desktop/AD_PDF/Micro-array/ssGSEA/AD_Aging_count.txt",sep="\t",header=T,check.names=F)
saveRDS(rt,"AD_Aging_count.RDS")
rt <- readRDS("C:/Users/sqyql/Desktop/AD_PDF/Micro-array/ssGSEA/AD_Aging_count.RDS")
write.csv(rt,"AD_Aging.csv")
#1、单基因表达量分组+差异分析
#加载工具包
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(GSEABase)
library(GSVA)
workingDir = "C:/Users/sqyql/Desktop/AD_PDF/Micro-array/ssGSEA/TP63"
setwd(workingDir)
sgene="TP63"       #输入进行单基因GSEA的基因名称
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
#判断原始数据是否去了log
max(rt)
if(max(rt)>30) rt=log2(rt+1)     #rt最大值大于30则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

#未标准化
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.7)
if(ncol(rt)>40) par(cex = 0.5)   ###设置字体大小
pdf(file = "raw.pdf",width=5,height = 4)
boxplot(rt,las=2,col =cols ) ###绘图
dev.off()

#标准化
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.5)
if(ncol(rt1)>40) par(cex = 0.5)   ###设置字体大小
pdf(file = "nor.pdf",width=5,height = 4.5)
boxplot(rt1,las=2,col =cols ) ###绘图
dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file="norexp_AD_TP63.txt",sep="\t",quote=F,col.names = F)

####差异分析
data=rt1
#按该基因表达情况进行分组
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

#分组之后进行差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
#保存单基因分组之后的差异分析结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("1.","DIFF_all_TP63_AD.xls"),sep="\t",quote=F,col.names=F)

#2、以热图的形式展示差异最大的前30个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]
#分组标签
Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(High="red",Low="blue"  ))

pdf(file=paste0("2.", "DIFF_heatmap_TP63_AD.pdf"),height=7,width=6)
pheatmap(afExp[,rownames(Type1)],                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c("blue","white","red"))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()




#3、以火山图展示差异分析筛选标准
#这里以P=0.05和logFC=0作为差异分析时候的筛选标准
setP=0.05
setlogFC=0
Significant=ifelse((Diff$adj.P.Val<setP & abs(Diff$logFC)>setlogFC), ifelse(Diff$logFC>setlogFC,"Up","Down"), "Not")
#开始绘制
p = ggplot(Diff, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant),size=4)+               #size点的大小
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(setP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=setlogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-setlogFC), colour="gray", linetype="twodash",size=1)
#查看，不添加标记可以直接保存
p
#添加基因点标记，可自行根据差异分析的结果进行标记
point.Pvalue=0.01
point.logFc=1.5
#继续绘制
Diff$symbol=rownames(Diff)
pdf(paste0("3.", "DIFF_vol.pdf"),width=7,height=6)
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & adj.P.Val< point.Pvalue )
p+geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
dev.off()

#4、GSEA富集分析
#GSEA富集分析前数据整理
logFC_t=0
deg$g=ifelse(deg$adj.P.Val>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 1,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)

#将ID转换成symbol
library(tidyverse)
kegg <- kk2@result
# 定义函数：输入一个ENTREZID字符返回对应gene symbol的字符
entrezID2geneSymbol <- function(entrezID){
  # 将字符串切分为向量
  entrezID <-  as.vector(str_split(entrezID,"/",simplify = T))
  # ENTREZID转换为SYMBOL，返回向量
  SYMBOL <- bitr(entrezID,"ENTREZID","SYMBOL",org.Hs.eg.db)[,"SYMBOL"]
  # 再SYMBOL合并为字符串，并返回结果
  SYMBOL_chr <- paste0(SYMBOL,"",collapse ="/")
  return(SYMBOL_chr)
}
# 自定义函数进行转换
kegg$gene_name <- unlist(map(kegg$core_enrichment,entrezID2geneSymbol))
#导出保存
write.csv(kegg, "4_all_GSEA_TP63_AD.csv")


#5、可视化排序后分别取GSEA结果的前10个
num=10
pdf(paste0("4.","GSEA_TP63_AD.pdf"),width = 9,height = 7)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)], col = "black")
dev.off()

library(ggplot2)
plot <- gseaplot2(kk2,
                  title = "HIF-1 signaling pathway", 
                  "hsa04066", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
ggsave("gsea_HIF-1 signaling pathway.pdf", plot, width = 8, height = 8)
ggsave("gsea_HIF-1 signaling pathway.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "cGMP-PKG signaling pathway", 
                  "hsa04022", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("gsea_cGMP-PKG signaling pathway.pdf", plot, width = 8, height = 8)
ggsave("gsea_cGMP-PKG signaling pathway.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "ErbB signaling pathway", 
                  "hsa04012", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("gsea_ErbB signaling pathway.pdf", plot, width = 8, height = 8)
ggsave("gsea_ErbB signaling pathway.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Wnt signaling pathway", 
                  "hsa04310", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("gsea_Wnt signaling pathway.pdf", plot, width = 8, height = 8)
ggsave("gsea_Wnt signaling pathway.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "IL-17 signaling pathway", 
                  "hsa04657", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("gsea_IL-17 signaling pathway.pdf", plot, width = 8, height = 8)
ggsave("gsea_IL-17 signaling pathway.tiff", plot, width = 8, height = 8)

plot <- gseaplot2(kk2,
                  title = "Hippo signaling pathway - multiple species", 
                  "hsa04392", 
                  color="green", 
                  base_size = 20, 
                  subplots = 1:3, 
                  pvalue_table = T)
plot + theme(axis.text.x = element_text(size = 16),
             axis.text.y = element_text(size = 16))
ggsave("gsea_Hippo signaling pathway - multiple species.pdf", plot, width = 8, height = 8)
ggsave("gsea_Hippo signaling pathway - multiple species.tiff", plot, width = 8, height = 8)
