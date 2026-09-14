# ============================================
# NRXN1-high vs NRXN1-low Astrocytes → Secretome → NicheNet → Neurons
# ============================================
library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(nichenetr)
library(tidyverse)
library(cowplot)
library(pheatmap)
library(ggplot2)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(EnhancedVolcano)
# ============================================
hip_1$Celltype <- hip_1$Cluster
Idents(hip_1) <- "Celltype"
hip_1$Group <- droplevels(hip_1$Group)
DefaultAssay(hip_1) <- "RNA"
astro_sub <- subset(hip_1, subset = Celltype == "Astrocytes" & Group == "AD")
nrxn1_expr <- GetAssayData(astro_sub, layer = "data")["ENSG00000179915", ] %>% as.numeric()
astro_sub$NRXN1_group <- ifelse(nrxn1_expr > median(nrxn1_expr), "High", "Low")
expr <- GetAssayData(astro_sub, layer = "data") %>% as.matrix()
rownames(expr) <- gsub("\\..*", "", rownames(expr))
symbol <- mapIds(org.Hs.eg.db, keys = rownames(expr),
                 keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
expr_symbol <- expr
rownames(expr_symbol) <- symbol
expr_symbol <- expr_symbol[!is.na(rownames(expr_symbol)), ]
expr_symbol <- expr_symbol[!duplicated(rownames(expr_symbol)), ]


# ============================================
DE_NR <- FindMarkers(astro_sub, group.by = "NRXN1_group",
                     ident.1 = "High", ident.2 = "Low",
                     logfc.threshold = 0.25, min.pct = 0.1)
DE_NR$ENSEMBL <- rownames(DE_NR)
DE_NR$SYMBOL <- mapIds(org.Hs.eg.db, keys = DE_NR$ENSEMBL %>% gsub("\\..*", "", .),
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
DE_NR <- DE_NR[!is.na(DE_NR$SYMBOL), ]
head(DE_NR)
de_ensembl <- rownames(DE_NR) %>% gsub("\\..*", "", .)
de_symbol <- mapIds(org.Hs.eg.db, keys = de_ensembl, keytype = "ENSEMBL",
                    column = "SYMBOL", multiVals = "first")
de_symbol <- de_symbol[!is.na(de_symbol)]

# ============================================
# 1. NicheNet lr_network 
data_dir <- "D:/Biofo/GSE268609/metabolism_Astro/NicheNet_Data"
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
lr_network_file <- file.path(data_dir, "lr_network_human_21122021.rds")
ligand_target_matrix_file <- file.path(data_dir, "ligand_target_matrix_human_21122021.rds")
weighted_networks_file <- file.path(data_dir, "weighted_networks_human_21122021.rds")
if(!file.exists(lr_network_file)) {
  cat("[DOWNLOAD] Fetching NicheNet prior knowledge networks from Zenodo...\n")
  download.file("https://zenodo.org/record/3260758/files/lr_network.rds", lr_network_file, mode="wb")
  download.file("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds", ligand_target_matrix_file, mode="wb")
  download.file("https://zenodo.org/record/3260758/files/weighted_networks.rds", weighted_networks_file, mode="wb")}
lr_network <- readRDS(lr_network_file)
ligand_target_matrix <- readRDS(ligand_target_matrix_file)
weighted_networks <- readRDS(weighted_networks_file)
secreted_ligands <- intersect(de_symbol, lr_network$from)
length(secreted_ligands)
secreted_ligands_scores <- rowMeans(expr_symbol[secreted_ligands, , drop=FALSE])
top_ligands <- names(sort(secreted_ligands_scores, decreasing = TRUE)[1:20])
top_ligands

# ============================================
neuron_clusters <- grep("CA_neurons|CA2-4_neurons|GABA_neurons",
                        unique(hip_1$Celltype), value = TRUE)
neurons_list <- lapply(neuron_clusters, function(x){
  subset(hip_1, subset = Celltype == x & Group %in% c("AD","HA"))
})
names(neurons_list) <- neuron_clusters
final_table <- data.frame()
for(subtype in neuron_clusters){
  message("Processing: ", subtype)
  neurons_sub <- neurons_list[[subtype]]
  # -------------------------
  # DE analysis (AD vs HA)
  # -------------------------
  DE_targets <- FindMarkers(neurons_sub, group.by = "Group", ident.1 = "AD",
    ident.2 = "HA", logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
  DE_targets$gene <- rownames(DE_targets)
  # SYMBOL mapping
  de_symbol <- mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", DE_targets$gene),
    keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  DE_targets$symbol <- de_symbol
  DE_targets <- DE_targets[!is.na(DE_targets$symbol), ]
  # -------------------------
  # background genes
  # -------------------------
  bg_symbol <- mapIds( org.Hs.eg.db, keys = gsub("\\..*", "", rownames(neurons_sub)),
    keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  bg_symbol <- bg_symbol[!is.na(bg_symbol)]
  # -------------------------
  # NicheNet ligand activity
  # -------------------------
  ligand_act <- predict_ligand_activities(
    geneset = DE_targets$symbol,
    background_expressed_genes = bg_symbol,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = top_ligands)
  ligand_act <- ligand_act %>% arrange(desc(pearson))
  top_ligands_sel <- head(ligand_act$test_ligand, 20)
  for(lig in top_ligands_sel){
    targets <- intersect(names(ligand_target_matrix[, lig]), DE_targets$symbol)
    if(length(targets) == 0) next
    for(tg in targets){
      tmp <- data.frame(Subtype = subtype, Ligand = lig, Target = tg,
        LigandActivity = ligand_act$pearson[ligand_act$test_ligand == lig],
        Target_logFC = DE_targets$avg_log2FC[match(tg, DE_targets$symbol)],
        Target_padj = DE_targets$p_val_adj[match(tg, DE_targets$symbol)],
        Target_p_val = DE_targets$p_val[match(tg, DE_targets$symbol)],
        Target_pct.1 = DE_targets$pct.1[match(tg, DE_targets$symbol)],
        Target_pct.2 = DE_targets$pct.2[match(tg, DE_targets$symbol)])
      final_table <- rbind(final_table, tmp)}}}

final_table <- final_table %>% arrange(desc(LigandActivity))
write.csv(final_table, "NRXN1_highlow_Astrocytes_Ligand_Receptor_Neuronsub_withActivity.csv", row.names = FALSE)
cat("DONE: Publication-ready table exported!\n")


# ==============================================================================
unified_theme <- theme(
  plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 10)),
  axis.title = element_text(face = "bold", size = 11, color = "black"),
  axis.text = element_text(size = 10, color = "black"),
  panel.grid = element_blank(),
  plot.tag = element_text(face = "bold", size = 16, color = "black")
)
# ============================================
# Fig5A: NRXN1 expression
# ============================================
Fig5A <- VlnPlot(astro_sub, features = "ENSG00000179915", group.by = "NRXN1_group")+
  labs(title = "A. Distribution of NRXN1 expression", x = "Group", y = "Expression Level of NRXN1") +
  unified_theme +
  theme(legend.position="none")
ggsave(file="fig5A_Distribution of NRXN1 expression.pdf", Fig5A, width = 8,height = 6)
ggsave(file="fig5A_Distribution of NRXN1 expression.png", Fig5A, width = 8,height = 6,dpi = 600)

# ============================================
# Fig5B: Astrocyte DE volcano（NRXN1-high vs low）
# ============================================
# top_genes <- DE_NR %>% arrange(p_val_adj) %>% head(50) %>% pull(SYMBOL)
# Fig5B <- EnhancedVolcano(DE_NR, lab = DE_NR$SYMBOL, x = "avg_log2FC", y = "p_val_adj")
Fig5B <- EnhancedVolcano(
  DE_NR, lab = DE_NR$SYMBOL, selectLab = secreted_ligands,
  x = "avg_log2FC", y = "p_val_adj", title = "B. NRXN1-high vs NRXN1-low Astrocytes",
  subtitle = NULL, caption = NULL, pCutoff = 0.05, FCcutoff = 0.25, pointSize = 2, 
  labSize = 5, legendPosition = "top", legendLabSize = 10, legendIconSize = 3) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10, color = "black"))
Fig5B <- Fig5B +
  labs(
    x = expression(bold(Log[2]~Fold~Change)),
    y = expression(bold(-Log[10]~Adjusted~P))
  ) +
  theme(axis.title.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title.y = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "top")
# 保存
ggsave(file="fig5B_NRXN1-high vs NRXN1-low Astrocytes.pdf", Fig5B, width = 8,height = 6)
ggsave(file="fig5B_NRXN1-high vs NRXN1-low Astrocytes.png", Fig5B, width = 8,height = 6,dpi = 600)

# ============================================
# Fig5C: Ligand activity ranking (barplot)（NicheNet核心）
# ============================================
Fig5C <- ligand_act %>% arrange(desc(pearson)) %>% head(20) %>%
  ggplot(aes(x = reorder(test_ligand, pearson), y = pearson, fill = auroc)) +
  geom_col() + coord_flip() +
  scale_fill_gradientn(colours = c("#FDEDEC", "#E74C3C", "#7B241C")) +
  labs(title = "C. Ligand Activity Ranking", fill = "AUROC", y = "Pearson Score", x = NULL) +
  unified_theme
# 保存
ggsave(file="fig5C_Ligand Activity Ranking.pdf", Fig5C, width = 8,height = 6)
ggsave(file="fig5C_Ligand Activity Ranking.png", Fig5C, width = 8,height = 6,dpi = 600)

# ============================================
# Fig5D: Ligand-Receptor network
# ============================================
# 选择 top 20 ligands
top_ligands_sel <- head(ligand_act$test_ligand[order(-ligand_act$pearson)], 20)
# 构建 ligand->receptor edges
edges <- lr_network %>%
  filter(from %in% top_ligands_sel) %>%
  dplyr::select(from, to)
# igraph object
g <- graph_from_data_frame(edges, directed = TRUE)
# network plot
Fig5D <- ggraph(g, layout = "fr") +
  geom_edge_link(arrow = arrow(length = unit(3, "mm")),
                 end_cap = circle(3, "mm"), alpha = 0.5, linewidth = 0.4) +
  geom_node_point(aes(color = ifelse(name %in% top_ligands_sel, "Ligand", "Receptor")), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Ligand"   = "#C0392B", "Receptor" = "#2E86C1")) +
  labs(title = "D. Ligand-Receptor Network", color = NULL) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        legend.position = "top", legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10, color = "black"))
ggsave(file="fig5D_Ligand-Receptor Network.pdf", Fig5D, width = 15,height = 6)
ggsave(file="fig5D_Ligand-Receptor Network.png", Fig5D, width = 15,height = 6,dpi = 600)

# ============================================
# Fig5E: Neuron subtype bubble plot (Ligand activity across subtypes)
# ============================================
bubble_df <- final_table %>% filter(Ligand %in% top_ligands_sel) %>% distinct(Subtype, Ligand, LigandActivity)
Fig5E <- ggplot(bubble_df, aes(x = Ligand, y = Subtype)) +
  geom_point(aes(size = abs(LigandActivity), color = LigandActivity)) +
  scale_color_gradient2(low = "white", mid = "#FDEDEC", high = "#C0392B") +
  scale_size_continuous(range = c(3,10)) +
  theme_bw() +
  unified_theme +
  labs(title = "E. Ligand Activity Across Neuronal Subtypes", color = "Pearson", size = "|Pearson|") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file="fig5E_Ligand Activity Across Neuronal Subtypes.pdf", Fig5E, width = 15,height = 4)
ggsave(file="fig5E_Ligand Activity Across Neuronal Subtypes.png", Fig5E, width = 15,height = 4,dpi = 600)

