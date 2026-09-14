

if(!"BiocManager" %in% installed.packages()){install.packages('BiocManager')}
if(!"openxlsx" %in% installed.packages()){install.packages('openxlsx')}
if(!"Seurat" %in% installed.packages()){BiocManager::install('Seurat')}
if(!"SeuratObject" %in% installed.packages()){BiocManager::install('SeuratObject')}
if(!"pbapply" %in% installed.packages()){BiocManager::install('pbapply')}
if(!"remotes" %in% installed.packages()){install.packages('remotes')}
if(!"scTenifoldNet" %in% installed.packages()){remotes::install_github('cailab-tamu/scTenifoldNet')}
if(!"pbapply" %in% installed.packages()){remotes::install_github('psolymos/pbapply')}
if(!"scTenifoldKnk" %in% installed.packages()){remotes::install_local('D:/R/scTenifoldKnk.zip',type = 'source', repos = NULL)}
if(!"dplyr" %in% installed.packages()){install.packages('dplyr')}
if(!"ggplot2" %in% installed.packages()){install.packages('ggplot2')}
if(!"ggrepel" %in% installed.packages()){install.packages('ggrepel')}
if(!"parallel" %in% installed.packages()){install.packages('parallel')}
options(warn = -1)
library(Seurat)
library(SeuratObject)
library(openxlsx)
library(scTenifoldKnk)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(scTenifoldNet)

# ==============================================================================
unified_theme <- theme(
  plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 10)),
  axis.title = element_text(face = "bold", size = 11, color = "black"),
  axis.text = element_text(size = 10, color = "black"),
  panel.grid = element_blank(),
  plot.tag = element_text(face = "bold", size = 16, color = "black")
)

hip <- readRDS("GSE268609_seurat_object.rds")
desired_groups <- c("AD", "MCI", "HA")
hip_1 <- hip[, hip$Group %in% desired_groups]
DefaultAssay(hip_1) <- 'RNA'
hip_1$Celltype <- hip_1$Cluster
table(hip_1$Celltype)
Idents(hip_1) <- hip_1$Celltype

core_genes_ens <- c("ENSG00000179915", "ENSG00000157152","ENSG00000064763","ENSG00000152503")
hvg_genes <- Seurat::VariableFeatures(hip_1)[1:5000]
final_features <- unique(c(hvg_genes, core_genes_ens))
hip_1_select_features <- subset(hip_1, features = final_features, idents = "Astrocytes")
table(hip_1$Celltype)
set.seed(666)
head(Cells(hip_1_select_features))
head(rownames(hip_1_select_features@meta.data))
meta_df <- hip_1_select_features@meta.data
meta_df$Group <- trimws(as.character(meta_df$Group))
groups_use <- c("AD", "MCI", "HA")
table(meta_df$Group)
n_min <- min(table(meta_df$Group[meta_df$Group %in% groups_use]))
cells_use <- unlist(lapply(
  groups_use,
  function(g){
    x <- rownames(meta_df)[meta_df$Group == g]
    cat(g, length(x), "\n")
    sample(x, n_min)
  }
))

hip_1_select <- subset(hip_1_select_features, cells = cells_use)
table(hip_1_select$Group)

hip_AD  <- subset(hip_1_select, subset = Group == "AD")
hip_MCI <- subset(hip_1_select, subset = Group == "MCI")
hip_HA  <- subset(hip_1_select, subset = Group == "HA")

sc_AD <- SeuratObject::LayerData(
  hip_AD,
  assay = "RNA",
  layer = "counts"
)

gene_ids <- rownames(sc_AD)


gene_ids_clean <- sub("\\..*$", "", gene_ids)

# ENSG -> SYMBOL
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols[is.na(gene_symbols)] <- gene_ids_clean[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)
rownames(sc_AD) <- gene_symbols
dir_AD <- ""
if(!dir.exists(dir_AD)){
  dir.create(dir_AD, recursive = TRUE)
}
sc_MCI <- SeuratObject::LayerData(hip_MCI, assay = "RNA", layer = "counts")
gene_ids <- rownames(sc_MCI)
gene_ids_clean <- sub("\\..*$", "", gene_ids)
# ENSG -> SYMBOL
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
gene_symbols[is.na(gene_symbols)] <- gene_ids_clean[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)
rownames(sc_MCI) <- gene_symbols
dir_MCI <- ""
if(!dir.exists(dir_MCI)){
  dir.create(dir_MCI, recursive = TRUE)
}

sc_HA <- SeuratObject::LayerData(
  hip_HA,
  assay = "RNA",
  layer = "counts"
)

gene_ids <- rownames(sc_HA)
gene_ids_clean <- sub("\\..*$", "", gene_ids)

# ENSG -> SYMBOL
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols[is.na(gene_symbols)] <- gene_ids_clean[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)
rownames(sc_HA) <- gene_symbols
dir_HA <- ""
if(!dir.exists(dir_HA)){
  dir.create(dir_HA, recursive = TRUE)
}
saveRDS(sc_HA, file = paste0(dir_HA,'/sc_HA.rds'))

c("NRXN1","SYN2","FAR2","TRIM36") %in% rownames(sc_AD)
c("NRXN1","SYN2","FAR2","TRIM36") %in% rownames(sc_MCI)
c("NRXN1","SYN2","FAR2","TRIM36") %in% rownames(sc_HA)


##### AD样本中虚拟过表达NRXN1=================================================
# =========================
# Step 1: Construct GRN
# =========================
gOE <- 'NRXN1' 
set.seed(1234)
# avoid using all CPU cores
n_cores <- max(1, parallel::detectCores() - 2)
cat("Using", n_cores, "cores...\n")
AD_networks <- scTenifoldNet::makeNetworks(
  X       = sc_AD,
  nNet    = 10,
  nc_nCells = min(1000, ncol(sc_AD)),
  nComp   = 3,  # Principal Components Regression
  nCores  = n_cores
)
# tensor decomposition
AD_tensor <- scTenifoldNet::tensorDecomposition(AD_networks, K = 3)
AD_network <- AD_tensor$X

###===================================================================
gOE <- "NRXN1"
AD_OE_network <- AD_network
expr_HA <- median(sc_HA["NRXN1", ])
expr_AD <- median(sc_AD["NRXN1", ])
fc <- expr_HA / expr_AD
cat(fc)
# 1.625
AD_OE_network[gOE, ] <- AD_OE_network[gOE, ] * 1.625
AD_OE_network[, gOE] <- AD_OE_network[, gOE] * 1.625
# =========================
# Step 3: Manifold Alignment
# =========================
MA_results <- scTenifoldNet::manifoldAlignment(
  t(as.matrix(AD_network)),
  t(as.matrix(AD_OE_network)))
# =========================
# Step 4: Calculate distances
# =========================
# obtain gene list from X_ labels
geneList <- rownames(MA_results)[grep("^X_", rownames(MA_results))]
geneList <- gsub("^X_", "", geneList)
# calculate WT vs OE manifold distance
distances <- sapply(geneList, function(gene) {
  wt_name <- paste0("X_", gene)
  oe_name <- paste0("Y_", gene)
  if (!(wt_name %in% rownames(MA_results))) {
    return(NA)
  }
  if (!(oe_name %in% rownames(MA_results))) {
    return(NA)
  }
  wt_vec <- MA_results[wt_name, ]
  oe_vec <- MA_results[oe_name, ]
  as.numeric(dist(rbind(wt_vec, oe_vec)))
})

# remove NA
valid_idx <- !is.na(distances)
distances <- distances[valid_idx]
geneList  <- geneList[valid_idx]

# =========================
# Step 5: Box-Cox transformation
# =========================
positive_distances <- distances[distances > 0]
if (length(positive_distances) > 10) {
  bc <- MASS::boxcox(positive_distances ~ 1, plotit = FALSE)
  lambda_val <- bc$x[which.max(bc$y)]
  cat("Optimal lambda =", lambda_val, "\n")
  if (abs(lambda_val) < 1e-6) {
    transformed_distances <- log(distances)
  } else {
    transformed_distances <-
      (distances^lambda_val - 1) / lambda_val
  }
} else {
  warning("Too few positive distances. Using raw distances.")
  transformed_distances <- distances
}
# =========================
# Step 6: Z-score
# =========================
Z_scores <- as.numeric(scale(transformed_distances))
# =========================
# Step 7: Empirical p-values
# =========================
pValues <- sapply(distances, function(x) {mean(distances >= x)})
pAdjusted <- p.adjust(pValues, method = "fdr")
# =========================
# Step 8: Fold change-like score
# =========================
mean_bg_dist_sq <- mean(distances[geneList != gOE]^2)
FC <- distances^2 / mean_bg_dist_sq

# =========================
# Step 9: Create result table
# =========================
library(dplyr)
library(ggplot2)
library(ggrepel)
dr_results <- data.frame(
  gene      = geneList,
  distance  = as.numeric(distances),
  Z         = Z_scores,
  FC        = FC,
  p.value   = pValues,
  p.adj     = pAdjusted
) %>%
  arrange(desc(Z))

# =========================
# Step 10: Save results
# =========================
xlsx_file <- file.path(".", 
                       paste0(gOE, "_virtual_OE_AD-1.625Fold.xlsx"))
openxlsx::write.xlsx(dr_results, file = xlsx_file, overwrite = TRUE)
# =========================
# Step 11: Visualization
# =========================
dr_results_plot <- dr_results %>%
  mutate(significance = ifelse(p.value < 0.05, "Significant", "Not Significant"),
    logP = -log10(p.value + 1e-300)) %>%
  arrange(p.value) %>%   # 按显著性排序
  mutate(label = if_else(row_number() <= 30, gene, NA_character_))
diff_reg_plot <- ggplot(dr_results_plot, aes(x = Z, y = logP)) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label), max.overlaps = 30, size = 5.5, box.padding = 0.5) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  labs(title = paste0("B. Virtual Restoration of ", gOE),
       subtitle ="Network Perturbation Analysis", x = "Z-score", y = expression(-Log[10]~P)) +
  unified_theme +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(diff_reg_plot)

# =========================
# Step 12: Save plots
# =========================
pdf_file <- file.path("./", "Figure 6B.pdf")
ggsave(filename = pdf_file, plot = diff_reg_plot, width = 12, height = 9)

# =========================
# Step 13: Significant genes
# =========================
head(dr_results)
summary(dr_results$Z)
summary(dr_results$p.value)
table(dr_results$p.value < 0.05)
# sig_genes <- dr_results %>% filter(p.value < 0.05, Z > 2)
sig_genes <- dr_results %>% filter(p.value < 0.05, Z > 1)
cat("\n")
cat("Number of significant genes:", nrow(sig_genes), "\n")
head(sig_genes)
# save significant genes
sig_file <- file.path("./", paste0(gOE, "_significant_genes_AD-1.625Fold.xlsx"))
openxlsx::write.xlsx(sig_genes, file = sig_file, overwrite = TRUE)

# =======================================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(stringr)
top_200 <- dr_results %>%
  arrange(desc(Z)) %>%
  head(200)
gene.df <- bitr(top_200$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", 
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_res <- as.data.frame(ego@result)
head(ego_res[, c("Description","Count","p.adjust")],50)
ekegg <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
ekegg_res <- as.data.frame(ekegg@result)
head(ekegg_res[, c("Description","Count","p.adjust")],50)

### KEGG Pathway Plot===================================================================
selected_pathways <- c(
  "AGE-RAGE signaling pathway in diabetic complications","ECM-receptor interaction",
  "Long-term depression","Platelet activation","PI3K-Akt signaling pathway",
  "Focal adhesion","Glutamatergic synapse","Integrin signaling",
  "Gap junction","Rap1 signaling pathway","Circadian entrainment",
  "Inflammatory mediator regulation of TRP channels","Long-term potentiation",
  "Phospholipase D signaling pathway","Serotonergic synapse","Calcium signaling pathway",
  "Dopaminergic synapse","Phosphatidylinositol signaling system","Apelin signaling pathway",
  "Retrograde endocannabinoid signaling","cGMP-PKG signaling pathway","Wnt signaling pathway",
  "Regulation of actin cytoskeleton","Ether lipid metabolism","Axon guidance",
  "Cytokine-cytokine receptor interaction","MAPK signaling pathway","Neuroactive ligand signaling",
  "Glycerophospholipid metabolism","Arachidonic acid metabolism")

kegg_plot <- ekegg_res %>%
  filter(Description %in% selected_pathways) %>%
  mutate(
    GeneRatio_num = as.numeric(str_split_fixed(GeneRatio, "/", 2)[,1]) /
      as.numeric(str_split_fixed(GeneRatio, "/", 2)[,2]),
    GeneRatio_num = round(GeneRatio_num, 2) 
  ) %>% arrange(desc(-Count))
kegg_plot$Description <- factor(kegg_plot$Description, levels = kegg_plot$Description)
Fig6C <- ggplot(kegg_plot, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_x_continuous(breaks = c(0.03, 0.06, 0.09, 0.12, 0.15)) + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw(base_size = 14) +
  labs(x = "Gene Ratio", y = NULL, color = expression(-Log[10]~P),
       title = "D. KEGG pathway") +
  unified_theme +
  theme(panel.grid = element_blank())
# 保存
ggsave(file="./Fig6C-KEGG Pathway.pdf", Fig6C, width = 8,height = 8)
ggsave(file="./Fig6C-KEGG Pathway.png", Fig6C, width = 8,height = 8,dpi = 600)

### GO Plot========================================================================
library(clusterProfiler)
library(ggplot2)
library(dplyr)
selected_BP <- c("cell-matrix adhesion","integrin-mediated signaling pathway","axonogenesis",
                 "regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction",
                 "regulation of angiogenesis","ether biosynthetic process","icosanoid biosynthetic process",
                 "platelet activation","axon guidance","neuron projection guidance","regulation of endocytosis",
                 "neurotransmitter uptake","icosanoid metabolic process","regulation of cell-matrix adhesion",
                 "regulation of cell-substrate adhesion","regulation of cell junction assembly",
                 "cell aggregation","calcium-mediated signaling","glutamate secretion",
                 "peripheral nervous system axon regeneration","amino acid neurotransmitter reuptake",
                 "retrograde trans-synaptic signaling","regulation of synaptic transmission, GABAergic",
                 "neurotransmitter transport","neuronal action potential propagation",
                 "protein localization to axon","ether lipid biosynthetic process",
                 "glycerol ether biosynthetic process","ether lipid metabolic process",
                 "regulation of sequestering of calcium ion")

go_plot <- ego_res %>%
  filter(Description %in% selected_BP) %>%
  mutate(
    GeneRatio_num = as.numeric(str_split_fixed(GeneRatio, "/", 2)[,1]) /
      as.numeric(str_split_fixed(GeneRatio, "/", 2)[,2]),
    GeneRatio_num = round(GeneRatio_num, 2)
  ) %>%
  arrange(desc(-Count))
go_plot$Description <- factor(go_plot$Description, levels = go_plot$Description)
Fig6D <- ggplot(go_plot, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_x_continuous(breaks = c(0.03, 0.06, 0.09, 0.12, 0.15)) +  
  scale_color_gradient(low = "green", high = "red") +
  theme_bw(base_size = 14) +
  labs(x = "Gene Ratio", y = NULL, color = expression(-Log[10]~P),
       title = "C. Biological Process") +
  unified_theme +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 11))
ggsave(file="./Fig6D-Biological Process.pdf", Fig6D, width = 8,height = 8)


# GO enrichment map
packageVersion("enrichplot")
packageVersion("clusterProfiler")
library(enrichplot)
ego2 <- pairwise_termsim(ego)
Fig6E <- emapplot(ego2, showCategory = 30, layout = "kk") +
  scale_color_gradient(low = "red", high = "green") +
  ggtitle("E. GO Enrichment Network") + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0, face = "bold", size = 16) 
  )
ggsave(file="./Fig6E-Emapplot.pdf", Fig6E, width = 12,height = 10)


fc <- dr_results$Z
names(fc) <- dr_results$gene
Fig6F <- cnetplot(ego, showCategory = 5, foldChange = fc, node_label = "category")+
  scale_color_gradient(low = "green", high = "red") +
  ggtitle("F. Gene–Concept Network") + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0, face = "bold", size = 16))

Fig6F <- Fig6F + 
  ggraph::geom_node_text(aes(label = name), 
                         repel = TRUE,    
                         size = 6)      
ggsave(file="./Fig6F-foldChange-cnetplot.pdf", Fig6F, width = 12,height = 10)



