#########Figure 1 A-D###################################################################################
library(WGCNA)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
#参数设置
treatment=c("age", "braak", "mmse", "nft", "pmi") #选择展示的行为学结果
fpkm=read.csv("exp.csv",header=T,row.names=1)
dim(fpkm)
names(fpkm)
# 转置后数据
datExpr0 = as.data.frame(t(fpkm))
############### 读取临床信息 ###############
behavior = read.table("WGCNA_cli.txt",
                      header = TRUE,
                      sep = "\t",
                      check.names = FALSE,
                      row.names = 1)
# 检查数据质量
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK)
{datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]}

behavior = behavior[match(rownames(datExpr0), rownames(behavior)), ]
sex_batch = factor(behavior$sex)
#design = model.matrix(~ age + braak + mmse + nft + pmi,data = behavior)
design = model.matrix(~ age + pmi, data = behavior)
datExpr0_corrected = removeBatchEffect(
  t(datExpr0),
  batch = sex_batch,
  design = design
)

datExpr0 = as.data.frame(t(datExpr0_corrected))
genes <- colnames(datExpr0)
annot <- AnnotationDbi::select(org.Hs.eg.db,keys = genes,
                               columns = c("SYMBOL","GENETYPE"),keytype = "SYMBOL")
head(annot)
protein_genes <- unique(annot$SYMBOL[annot$GENETYPE == "protein-coding"])
datExpr0 <- datExpr0[,colnames(datExpr0) %in% protein_genes]

keep <- colMeans(datExpr0 > 5) > 0.3
datExpr0 <- datExpr0[, keep]

gene_mad <- apply(datExpr0, 2, mad)
datExpr0 <- datExpr0[,gene_mad > quantile(gene_mad, 0.75)]
dim(datExpr0)
summary(as.numeric(as.matrix(datExpr0)))
datExpr <- datExpr0
write.csv(datExpr,file = "WGCNA_expression_matrix.csv",quote = FALSE)

sampleTree = hclust(dist(datExpr), method = "average")
clust = cutreeStatic(sampleTree,cutHeight = 120,minSize = 10)
table(clust)
pdf(file = "1_sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

#Loading clinical trait data
behavior=read.table("WGCNA_cli.txt", header=T, sep="\t", check.names=F, row.names=1)
traitData = behavior[, treatment, drop = FALSE]
traitData = as.data.frame(lapply(traitData, function(x) as.numeric(as.character(x))))
rownames(traitData) = rownames(behavior)
traitData = as.data.frame(scale(traitData))
max(traitData)
dim(traitData)
names(traitData)
allTraits = traitData
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr)
traitSamples =rownames(allTraits)
datTraits = allTraits[match(rownames(datExpr), rownames(allTraits)), ]
all(rownames(datExpr) == rownames(datTraits))
rownames(datTraits) 
sum(is.na(datTraits))
gc()
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
#sizeGrWindow(12,12)
pdf(file="2_Sample dendrogram and trait heatmap.pdf",width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",marAll = c(3,10,3,3))
dev.off()


allowWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:20), seq(from = 22, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower =sft$powerEstimate
softPower
softPower = 4
# Plot the results:
pdf(file="Figure 1A.pdf",width=15,height=9)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower

#adjacency = adjacency(datExpr, power = softPower)
adjacency = adjacency(datExpr, power = softPower, type = "signed")
##### Turn adjacency into topological overlap
#TOM = TOMsimilarity(adjacency);
TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# dynamicMods
#    1    2    3    4    5    6    7 
# 1858  353  301   92   76   64   44 

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# dynamicColors
# black      blue     brown     green       red turquoise    yellow 
#    44       353       301        76        64      1858        92 

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#moduleTraitCor = cor(MEs, datTraits, use = "p")
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Calculate dissimilarity of module eigengenes
# for (trait in names(datTraits)) {
# cor_result <- cor(MEs, datTraits[[trait]], use = "p")
# print(cbind(trait, cor_result))
# }
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
moduleColors = mergedColors

#sizeGrWindow(12, 9)
pdf(file="Figure 1B.pdf", width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Bootstrap-based module assignment stability (FIXED VERSION)
multiExpr <- list()
multiColor <- list()
names(moduleColors) <- colnames(datExpr)
multiExpr[["ref"]] <- list(data = datExpr)
multiColor[["ref"]] <- moduleColors
nBoot <- 10
sampleFrac <- 0.8
set.seed(123)
for (i in 1:nBoot) {
  samp <- sample(
    1:nrow(datExpr),
    floor(sampleFrac * nrow(datExpr))
  )
  dat_boot <- datExpr[samp, , drop = FALSE]
  multiExpr[[paste0("boot", i)]] <- list(data = dat_boot)
  # 使用reference module labels
  multiColor[[paste0("boot", i)]] <- moduleColors
}
mp <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  nPermutations = 50,
  randomSeed = 123,
  verbose = 3
)

table(sapply(multiColor, length))
names(mp$preservation$Z$ref.ref)
pres1 <- mp$preservation$Z$ref.ref$inColumnsAlsoPresentIn.boot1
head(pres1)
bootNames <- grep(
  "boot",
  names(mp$preservation$Z$ref.ref),
  value = TRUE
)

zMat <- sapply(bootNames, function(x){
  
  mp$preservation$Z$ref.ref[[x]][,"Zsummary.pres"]
  
})

summaryDF <- data.frame(
  Module = rownames(mp$preservation$Z$ref.ref[[bootNames[1]]]),
  
  Mean_Zsummary = rowMeans(zMat, na.rm = TRUE),
  
  SD_Zsummary = apply(zMat, 1, sd, na.rm = TRUE)
)

summaryDF <- summaryDF[
  order(summaryDF$Mean_Zsummary, decreasing = TRUE),
]

summaryDF

library(ggplot2)

plotDF <- summaryDF

# 去掉gold和grey
plotDF <- subset(
  plotDF,
  !Module %in% c("gold", "grey")
)

# 保持模块顺序
plotDF$Module <- factor(
  plotDF$Module,
  levels = plotDF$Module
)

p <- ggplot(
  plotDF,
  aes(
    x = Module,
    y = Mean_Zsummary,
    fill = Module
  )
) +
  
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_errorbar(
    aes(
      ymin = Mean_Zsummary - SD_Zsummary,
      ymax = Mean_Zsummary + SD_Zsummary
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    linewidth = 0.7,
    color = "steelblue"
  ) +
  
  geom_hline(
    yintercept = 10,
    linetype = "dashed",
    linewidth = 0.7,
    color = "firebrick"
  ) +
  
  scale_fill_identity() +
  
  coord_flip() +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    
    axis.line = element_line(
      linewidth = 0.8
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  ) +
  
  labs(
    x = NULL,
    y = "Mean Zsummary",
    title = "Bootstrap-based Module Robustness"
  )

# 保存PDF
ggsave("Figure 1C.pdf", plot = p, width = 7, height = 5.2, device = cairo_pdf)




# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(length(unique(moduleColors))))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file="Figure 1D.pdf",width=9,height=8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,font.lab.x = 0.5,font.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

#####
names(datExpr)
probes = names(datExpr)

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = "GS_and_MM.xls",sep="\t",row.names=F)




#########Figure 1E,Figure 1F###################################################################################
library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library("stringr")
library("pathview")
pvalueFilter=0.05         
qvalueFilter=1  
showNum=6

rt=read.table("target_blue.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA) 
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,] 
#rt$entrezID[851] <- '221938'
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="Figure 1E.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bub)
dev.off()

      
showNum=20
keggId="hsa"
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="Figure 1F.pdf",width = 7,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()


#########Supplementary Figure 3###################################################################################
unique(geneInfo$moduleColor)
targetModule <- "blue"
moduleData <- geneInfo[geneInfo$moduleColor == targetModule, ]
moduleData$GS_consensus <- (abs(moduleData$GS.braak) + 
                              abs(moduleData$GS.mmse) + 
                              abs(moduleData$GS.nft)) / 3

hubGenes <- moduleData[moduleData$GS_consensus > 0.25 & abs(moduleData$MMblue) > 0.7, ]
hubGenes <- hubGenes[order(hubGenes$GS_consensus, decreasing = TRUE), ]


#########--------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
moduleData <- geneInfo[geneInfo$moduleColor == targetModule, ]
moduleData <- moduleData[order(-abs(moduleData[[paste0("MM", targetModule)]])),]

# trait----------------mmse
trait <- "GS.mmse"
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])

plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)

plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_mmse_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

# trait----------------mmse
trait <- "GS.mmse"
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])
plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)
plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_mmse_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

# trait----------------braak
trait <- "GS.braak"
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])
plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)
plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_braak_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

# trait----------------nft
trait <- "GS.nft"
# 提取MM和GS
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])
plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)
plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_nft_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

# trait----------------age
trait <- "GS.age"
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])
plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)
plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_age_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

# trait----------------pmi
trait <- "GS.pmi"
MM <- abs(moduleData[[paste0("MM", targetModule)]])
GS <- abs(moduleData[[trait]])
plotData <- data.frame(Gene = rownames(moduleData),MM = MM,GS = GS)
genes_to_label <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
topGenes <- subset(plotData,Gene %in% genes_to_label)
plotData$highlight <- ifelse(plotData$Gene %in% genes_to_label,"Selected","Other")

corVal <- cor(plotData$MM, plotData$GS, use = "p")
corP <- corPvalueStudent(corVal,nrow(plotData))

p <- ggplot(plotData,aes(MM, GS)) +
  geom_point(aes(color = highlight),size = 2.4,alpha = 0.8) +
  scale_color_manual(values = c("Other" = "grey80","Selected" = targetModule)) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 1) +
  geom_text_repel(data = topGenes,aes(label = Gene),size = 4,
                  max.overlaps = Inf,box.padding = 0.4,point.padding = 0.3) +
  labs(title = paste0("Key Genes in ",targetModule," Module"),
       subtitle = paste0("Correlation = ",round(corVal, 3),", P = ",signif(corP, 3)),
       x = paste0("Module Membership (|MM|) in ",targetModule," Module"),
       y = paste0("Gene Significance (|GS|) for ",trait)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",plot.title = element_text(face = "bold",size = 16,hjust = 0.5),
        plot.subtitle = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(face = "bold",size = 14),
        axis.text = element_text(color = "black",size = 12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.1,color = "black"))
ggsave(filename = paste0("Lasso_MM_GS_pmi_",targetModule,".pdf"),
       plot = p,width = 7,height = 6)

