# ==============================================================================
# 🌟 INTEGRATED EPIGENETIC-METABOLIC PIPELINE FOR ASTROCYTES (PRODUCTION-READY)
# ==============================================================================
cat("\n[INIT] Loading mandatory environments and visual engines...\n")

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(AUCell)
  library(GSEABase)
  library(scMetabolism)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
  library(ggrepel)
  library(RColorBrewer)
  library(tidyr)
})



my_colors <- c("HA" = "#3A86C8", "MCI" = "#E99B26", "AD" = "#D9383A")
expr_cols <- c("#E5E8E8", "#AED6F1", "#3498DB", "#9B59B6", "#7D3C98")

cat("\n[STEP 2] Loading RDS raw dataset...\n")
hip_1 <- readRDS("D:/Biofo/GSE268609/hip_1.rds")
hip_1$Celltype <- hip_1$Cluster
Idents(hip_1) <- "Celltype"
hip_1$Group <- droplevels(hip_1$Group)
table(hip_1$Group)
Assays(hip_1)
Astrocytes <- subset(hip_1, subset = Celltype == "Astrocytes")
Astrocytes$Celltype <- droplevels(Astrocytes$Celltype)

Astrocytes_sub <- Astrocytes
DefaultDimReduc(Astrocytes_sub) <- "umap"

cat("\n[STEP 4] Normalizing gene naming and Ensembl-to-Symbol mapping...\n")
DefaultAssay(Astrocytes_sub) <- "RNA"
raw_features <- rownames(Astrocytes_sub)
clean_features <- gsub("\\..*", "", raw_features)
gene_symbol <- mapIds(org.Hs.eg.db, keys = clean_features, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

nrxn1_ens <- clean_features[gene_symbol == "NRXN1"] %>% na.omit()
syn2_ens  <- clean_features[gene_symbol == "SYN2"]  %>% na.omit()
far2_ens  <- clean_features[gene_symbol == "FAR2"]  %>% na.omit()
trim36_ens <- clean_features[gene_symbol == "TRIM36"]  %>% na.omit()

if(length(nrxn1_ens) == 0) stop("Fatal: NRXN1 transcripts missing from dataset.")
if(length(syn2_ens) == 0)  stop("Fatal: SYN2 transcripts missing from dataset.")
if(length(far2_ens) == 0)  stop("Fatal: FAR2 transcripts missing from dataset.")
if(length(trim36_ens) == 0)  stop("Fatal: TRIM36 transcripts missing from dataset.")

nrxn1_idx <- which(gene_symbol == "NRXN1")[1]
if(is.na(nrxn1_idx)) stop("Fatal: NRXN1 transcripts missing from dataset.")
nrxn1_ens <- raw_features[nrxn1_idx] 

nrxn1_dt_idx <- grep("NRXN1-DT|NRXN1.DT|NRXN1_DT", gene_symbol, ignore.case = TRUE)[1]
if(is.na(nrxn1_dt_idx)) {
  nrxn1_dt_idx <- grep("NRXN1-DT|NRXN1.DT", raw_features, ignore.case = TRUE)[1]
}
if(is.na(nrxn1_dt_idx)) stop("Fatal: Epigenetic cis-regulator NRXN1-DT missing.")
nrxn1_dt_ens <- raw_features[nrxn1_dt_idx] 


cat("\n[STEP 5] Building SYMBOL matrix and calculating single-cell metabolic AUC...\n")
expr <- GetAssayData(Astrocytes_sub, layer = "data")
rownames(expr) <- gene_symbol
keep_genes <- !is.na(rownames(expr)) & rownames(expr) != ""
expr_s <- expr[keep_genes, ]

unique_symbols <- unique(rownames(expr_s))
match_idx <- match(rownames(expr_s), unique_symbols)
map_mat <- sparseMatrix(i = match_idx, j = 1:nrow(expr_s), x = 1, 
                        dims = c(length(unique_symbols), nrow(expr_s)),
                        dimnames = list(unique_symbols, NULL))
expr_s_final <- map_mat %*% expr_s

rankings <- AUCell_buildRankings(expr_s_final, plotStats = FALSE, verbose = FALSE)
gmt <- getGmt(system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism"))
auc <- AUCell_calcAUC(gmt, rankings, verbose = FALSE)
auc_mat_clean <- getAUC(auc)

auc_mat_clean[is.na(auc_mat_clean)] <- 0
auc_mat_clean[is.infinite(auc_mat_clean)] <- 0
auc_mat_clean <- auc_mat_clean[apply(auc_mat_clean, 1, var) > 0, ]

suppressWarnings({
  Astrocytes_sub[["KEGG_AUC"]] <- CreateAssayObject(counts = auc_mat_clean)
})

common_cells <- intersect(colnames(auc_mat_clean), colnames(Astrocytes_sub))
auc_use <- auc_mat_clean[, common_cells]

# DefaultAssay(Astrocytes_sub) <- "RNA"
nrxn1 <- FetchData(Astrocytes_sub, vars = nrxn1_ens)[, 1] 
nrxn1 <- nrxn1[match(common_cells, colnames(Astrocytes_sub))]
group_use <- Astrocytes_sub$Group[match(common_cells, colnames(Astrocytes_sub))]

pathway_cor <- function(pathway_score, gene_expr){
  keep <- complete.cases(pathway_score, gene_expr)
  if(sum(keep) < 10) return(c(rho=NA, p=NA))
  ct <- suppressWarnings(cor.test(pathway_score[keep], gene_expr[keep], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

calc_group_cor <- function(group_name){
  idx <- group_use == group_name
  res <- t(apply(auc_use[,idx,drop=FALSE], 1, function(x){pathway_cor(pathway_score = x, gene_expr = nrxn1[idx])}))
  res <- as.data.frame(res)
  res$Pathway <- rownames(res)
  res$FDR <- p.adjust(res$p, method = "BH")
  res$Group <- group_name
  rownames(res) <- NULL
  res
}

cor_HA  <- calc_group_cor("HA")
cor_MCI <- calc_group_cor("MCI")
cor_AD  <- calc_group_cor("AD")

all_cor <- bind_rows(cor_HA, cor_MCI, cor_AD)
cor_wide <- all_cor %>%
  dplyr::select(Pathway, Group, rho) %>%
  tidyr::pivot_wider(names_from = Group, values_from = rho, values_fn = mean) %>%
  dplyr::mutate(HA = as.numeric(HA), MCI = as.numeric(MCI), AD = as.numeric(AD))
cor_wide$Delta_AD_HA <- cor_wide$AD - cor_wide$HA
cor_wide <- cor_wide %>% dplyr::filter(is.finite(Delta_AD_HA)) %>% dplyr::arrange(desc(Delta_AD_HA))
write.csv(cor_wide, "NRXN1_pathway_rewiring.csv", row.names = FALSE)
top20 <- cor_wide[order(abs(cor_wide$Delta_AD_HA), decreasing = TRUE), ][1:20, ]
print(head(top20, 20))
"
# A tibble: 20 × 5
   Pathway                                                          HA      MCI      AD Delta_AD_HA
   <chr>                                                         <dbl>    <dbl>   <dbl>       <dbl>
 1 Oxidative phosphorylation                                   0.0879  -0.101   -0.161      -0.249 
 2 Ether lipid metabolism                                      0.0190   0.0521   0.239       0.220 
 3 Phosphonate and phosphinate metabolism                      0.0489   0.0745   0.221       0.173 
 4 Biosynthesis of unsaturated fatty acids                    -0.0568  -0.0404   0.0795      0.136 
 5 Glycerolipid metabolism                                     0.0479   0.116    0.184       0.136 
 6 Histidine metabolism                                       -0.243   -0.141   -0.112       0.131 
 7 Glycosphingolipid biosynthesis - lacto and neolacto series -0.0540  -0.00672  0.0725      0.127 
 8 Glycosphingolipid biosynthesis - globo and isoglobo series -0.00366  0.0150   0.121       0.124 
 9 Purine metabolism                                          -0.0174   0.0321   0.107       0.124 
10 N-Glycan biosynthesis                                      -0.0998  -0.0816  -0.220      -0.120 
11 Glycosaminoglycan degradation                               0.270    0.167    0.151      -0.119 
12 Glycerophospholipid metabolism                              0.145    0.182    0.262       0.117 
13 Phenylalanine metabolism                                   -0.290   -0.198   -0.176       0.114 
14 Pyruvate metabolism                                         0.0901   0.100    0.204       0.114 
15 Steroid biosynthesis                                        0.0539   0.0278  -0.0573     -0.111 
16 Glycosylphosphatidylinositol (GPI)-anchor biosynthesis      0.0599   0.00911 -0.0466     -0.107 
17 Nicotinate and nicotinamide metabolism                      0.00150  0.00652 -0.104      -0.106 
18 Tryptophan metabolism                                      -0.213   -0.117   -0.112       0.101 
19 Cysteine and methionine metabolism                         -0.0352   0.0294   0.0650      0.100 
20 Glycosphingolipid biosynthesis - ganglio series            -0.124   -0.0874  -0.0267      0.0975
"

df <- data.frame(
  NRXN1 = FetchData(Astrocytes_sub, nrxn1_ens[1])[,1],
  SYN2  = FetchData(Astrocytes_sub, syn2_ens[1])[,1],
  FAR2  = FetchData(Astrocytes_sub, far2_ens[1])[,1],
  TRIM36 = FetchData(Astrocytes_sub, trim36_ens[1])[,1],
  Group = Astrocytes_sub$Group
)
df2 <- data.frame(
  NRXN1 = FetchData(Astrocytes_sub, nrxn1_ens[1])[,1],
  FAR2  = FetchData(Astrocytes_sub, far2_ens[1])[,1],
  Group = Astrocytes_sub$Group
)
df3 <- data.frame(
  NRXN1 = FetchData(Astrocytes_sub, nrxn1_ens[1])[,1],
  Group = Astrocytes_sub$Group
)

DefaultAssay(Astrocytes_sub) <- "KEGG_AUC"
df$Ether_Lipid <- FetchData(Astrocytes_sub, "Ether lipid metabolism")[,1]
df2$Ether_Lipid <- FetchData(Astrocytes_sub, "Ether lipid metabolism")[,1]
df3$Ether_Lipid <- FetchData(Astrocytes_sub, "Ether lipid metabolism")[,1]


cat("\n=== Per-Group Pathological Correlation Metrics ===\n")
unique(df$Group)
#df$Group <- droplevels(df$Group)
table(df$Group)
lapply(split(df, df$Group), function(sub_df){
  #cat("\nStage Context:", unique(sub_df$Group), "\n")
  #cat("\nStage Context:", levels(sub_df$Group)[unique(sub_df$Group)], "\n")
  cat("\nStage Context:", as.character(unique(sub_df$Group)), "\n")
  vars <- list(FAR2 = sub_df$FAR2, NRXN1 = sub_df$NRXN1, SYN2 = sub_df$SYN2, 
               TRIM36 = sub_df$TRIM36, Ether_Lipid = sub_df$Ether_Lipid)
  cor_safe <- function(x, y){
    complete <- complete.cases(as.numeric(x), as.numeric(y))
    if(sum(complete) < 2) return(NA)
    cor(as.numeric(x)[complete], as.numeric(y)[complete])
  }
  cat("  -> FAR2 vs SYN2 (Synchronization Index):", cor_safe(vars$FAR2, vars$SYN2), "\n")
  cat("  -> NRXN1 vs Ether_Lipid (Feedback Index):", cor_safe(vars$NRXN1, vars$Ether_Lipid), "\n")
  cat("  -> SYN2 vs Ether_Lipid  (Feedback Index):", cor_safe(vars$SYN2, vars$Ether_Lipid), "\n")
  cat("  -> FAR2 vs Ether_Lipid  (Feedback Index):", cor_safe(vars$FAR2, vars$Ether_Lipid), "\n")
  cat("  -> TRIM36 vs Ether_Lipid  (Feedback Index):", cor_safe(vars$TRIM36, vars$Ether_Lipid), "\n")
  cat("  -> NRXN1 vs SYN2        (Baseline Check) :", cor_safe(vars$NRXN1, vars$SYN2), "\n")
})
'
Stage Context: AD 
-> FAR2 vs SYN2 (Synchronization Index): 0.02396599 
-> NRXN1 vs Ether_Lipid (Feedback Index): 0.2320089 
-> SYN2 vs Ether_Lipid  (Feedback Index): -0.04969269 
-> FAR2 vs Ether_Lipid  (Feedback Index): 0.01937148 
-> TRIM36 vs Ether_Lipid  (Feedback Index): -0.01087594 
-> NRXN1 vs SYN2        (Baseline Check) : -0.06157218 

Stage Context: MCI 
-> FAR2 vs SYN2 (Synchronization Index): 0.009968031 
-> NRXN1 vs Ether_Lipid (Feedback Index): 0.04713819 
-> SYN2 vs Ether_Lipid  (Feedback Index): -0.007318377 
-> FAR2 vs Ether_Lipid  (Feedback Index): 0.02707372 
-> TRIM36 vs Ether_Lipid  (Feedback Index): 0.01890985 
-> NRXN1 vs SYN2        (Baseline Check) : -0.0396581 

Stage Context: HA 
-> FAR2 vs SYN2 (Synchronization Index): 0.03968466 
-> NRXN1 vs Ether_Lipid (Feedback Index): 0.01876344 
-> SYN2 vs Ether_Lipid  (Feedback Index): -0.01334274 
-> FAR2 vs Ether_Lipid  (Feedback Index): -0.004646051 
-> TRIM36 vs Ether_Lipid  (Feedback Index): -0.009323124 
-> NRXN1 vs SYN2        (Baseline Check) : 0.01341089 
$AD
NULL

$MCI
NULL

$HA
NULL
' 

cat("\n=== Multi-variable Linear Regression summary ===\n")
model <- lm(Ether_Lipid ~ NRXN1 + SYN2 + FAR2 + TRIM36 + Group, data = df)
print(summary(model))
# 
# Call:
#   lm(formula = Ether_Lipid ~ NRXN1 + SYN2 + FAR2 + TRIM36 + Group, 
#      data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.037919 -0.012470 -0.000851  0.011276  0.096108 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.983e-02  2.897e-04 102.962  < 2e-16 ***
#   NRXN1      1.287e-03  8.267e-05  15.569  < 2e-16 ***
#   SYN2      -1.052e-03  3.236e-04  -3.252  0.00115 ** 
#   FAR2       1.051e-03  5.100e-04   2.060  0.03941 *  
#   TRIM36    -2.203e-04  8.003e-04  -0.275  0.78315    
# GroupMCI    -3.260e-03  3.422e-04  -9.526  < 2e-16 ***
#   GroupHA   -1.578e-03  2.969e-04  -5.315 1.08e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0172 on 18106 degrees of freedom
# Multiple R-squared:  0.01658,	Adjusted R-squared:  0.01625 
# F-statistic: 50.87 on 6 and 18106 DF,  p-value: < 2.2e-16


### ① NRXN1 → Ether lipid（main）============================================
# Partial Correlation
library(ppcor)
pcor.test(df$NRXN1, df$Ether_Lipid, df[,c("FAR2","SYN2","TRIM36")])
# estimate     p.value statistic     n gp  Method
# 1 0.1039278 1.12852e-44  14.06128 18113  3 pearson


### ② NRXN1 → Ether lipid
df3 %>% group_by(Group) %>% summarise(cor = cor(NRXN1, Ether_Lipid, method="spearman"))
# # A tibble: 3 × 2
# Group    cor
# <fct>  <dbl>
# 1 AD    0.239 
# 2 MCI   0.0521
# 3 HA    0.0190



### ③ NRXN1 ↔ FAR2 =================================================
cor.test(df2$NRXN1, df2$FAR2, method="pearson")
# Pearson's product-moment correlation
# 
# data:  df2$NRXN1 and df2$FAR2
# t = -3.1729, df = 18111, p-value = 0.001512
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.03812058 -0.00901027
# sample estimates:
#         cor 
# -0.02357042

cor.test(df2$NRXN1, df2$FAR2, method="spearman")
# Spearman's rank correlation rho
# 
# data:  df2$NRXN1 and df2$FAR2
# S = 1.0093e+12, p-value = 0.01035
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.01905104 

#  Spearman 
df2 %>% group_by(Group) %>% summarise(cor = cor(NRXN1,FAR2,method="spearman"))
# # A tibble: 3 × 2
# Group     cor
# <fct>   <dbl>
# 1 AD    -0.0260
# 2 MCI   -0.0249
# 3 HA    -0.0235
cor_group <- data.frame(Group=c("HA","MCI","AD"), rho=c(-0.0235,-0.0249,-0.0260))
df2 %>%
  group_by(Group) %>%
  summarise(rho = cor(NRXN1, FAR2, method="spearman"),
    p = cor.test(NRXN1, FAR2, method="spearman")$p.value)
# # A tibble: 3 × 3
# Group     rho      p
# <fct>   <dbl>  <dbl>
# 1 AD    -0.0260 0.0302
# 2 MCI   -0.0249 0.107 
# 3 HA    -0.0235 0.0501

### ④ FAR2 → Ether lipid
pcor.test(df2$FAR2, df2$Ether_Lipid, df2$NRXN1)
#     estimate   p.value statistic    n gp  Method
# 1 0.01421536 0.05573782  1.913203 18113  1 pearson
pcor.test(df$FAR2, df$Ether_Lipid, df[,c("NRXN1","SYN2","TRIM36")])
# estimate   p.value statistic     n gp  Method
# 1 0.01492286 0.0446227  2.008333 18113  3 pearson

model_far2 <- lm(Ether_Lipid ~ FAR2 + NRXN1 + SYN2 + TRIM36 + Group, data=df)
summary(model_far2)
# Call:
#   lm(formula = Ether_Lipid ~ FAR2 + NRXN1 + SYN2 + TRIM36 + Group, 
#      data = df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.037919 -0.012470 -0.000851  0.011276  0.096108 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.983e-02  2.897e-04 102.962  < 2e-16 ***
#   FAR2         1.051e-03  5.100e-04   2.060  0.03941 *  
#   NRXN1        1.287e-03  8.267e-05  15.569  < 2e-16 ***
#   SYN2        -1.052e-03  3.236e-04  -3.252  0.00115 ** 
#   TRIM36      -2.203e-04  8.003e-04  -0.275  0.78315    
# GroupMCI      -3.260e-03  3.422e-04  -9.526  < 2e-16 ***
#   GroupHA     -1.578e-03  2.969e-04  -5.315 1.08e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0172 on 18106 degrees of freedom
# Multiple R-squared:  0.01658,	Adjusted R-squared:  0.01625 
# F-statistic: 50.87 on 6 and 18106 DF,  p-value: < 2.2e-16
confint(model_far2)
# 2.5 %        97.5 %
#   (Intercept)  2.926128e-02  0.0303969984
# FAR2         5.099039e-05  0.0020503653
# NRXN1        1.125093e-03  0.0014491799
# SYN2        -1.686792e-03 -0.0004180425
# TRIM36      -1.788879e-03  0.0013483648
# GroupMCI    -3.930440e-03 -0.0025889492
# GroupHA     -2.159690e-03 -0.0009959224

# 多变量回归进一步支持FAR2对NRXN1的弱负向调控作用
model2 <- lm(FAR2 ~ NRXN1 + Group,data = df2)
print(summary(model2))
# Call:
#   lm(formula = FAR2 ~ NRXN1 + Group, data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.0685 -0.0522 -0.0473 -0.0364  3.4929 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0521522  0.0041499  12.567  < 2e-16 ***
# NRXN1       -0.0043670  0.0012051  -3.624 0.000291 ***
# GroupMCI     0.0007282  0.0049885   0.146 0.883946    
# GroupHA      0.0163478  0.0043282   3.777 0.000159 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.251 on 18109 degrees of freedom
# Multiple R-squared:  0.001511,	Adjusted R-squared:  0.001346 
# F-statistic: 9.137 on 3 and 18109 DF,  p-value: 4.877e-06

model3 <- lm(Ether_Lipid ~ NRXN1*Group + SYN2 + FAR2 + TRIM36, data=df)
summary(model3)
# Call:
#   lm(formula = Ether_Lipid ~ NRXN1 * Group + SYN2 + FAR2 + TRIM36, 
#      data = df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.03978 -0.01247 -0.00081  0.01124  0.09416 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.0273797  0.0003540  77.346  < 2e-16 ***
#   NRXN1           0.0023133  0.0001190  19.433  < 2e-16 ***
#   GroupMCI        0.0011396  0.0007845   1.453  0.14635    
#   GroupHA         0.0040986  0.0005777   7.095 1.35e-12 ***
#   SYN2           -0.0009340  0.0003225  -2.896  0.00379 ** 
#   FAR2            0.0010925  0.0005081   2.150  0.03154 *  
#   TRIM36         -0.0001949  0.0007972  -0.244  0.80690    
# NRXN1:GroupMCI   -0.0016602  0.0002415  -6.874 6.47e-12 ***
#   NRXN1:GroupHA  -0.0020944  0.0001805 -11.604  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01714 on 18104 degrees of freedom
# Multiple R-squared:  0.02442,	Adjusted R-squared:  0.02399 
# F-statistic: 56.64 on 8 and 18104 DF,  p-value: < 2.2e-16
df$Group <- as.character(df$Group)

df$Group <- factor(
  df$Group,
  levels = c("HA", "MCI", "AD")
)

droplevels(df$Group)
df <- droplevels(df)


library(emmeans)
trend_result <- emtrends(model3, pairwise ~ Group, var = "NRXN1")
trend_result$emtrends
# Group NRXN1.trend       SE    df  lower.CL upper.CL
# AD       0.002313 0.000119 18104  0.002080 0.002547
# MCI      0.000653 0.000210 18104  0.000241 0.001065
# HA       0.000219 0.000136 18104 -0.000047 0.000485
# Confidence level used: 0.95 

trend_result$contrasts
# contrast estimate       SE    df t.ratio p.value
# AD - MCI 0.001660 0.000242 18104   6.874 <0.0001
# AD - HA  0.002094 0.000180 18104  11.604 <0.0001
# MCI - HA 0.000434 0.000250 18104   1.735  0.1922
# P value adjustment: tukey method for comparing a family of 3 estimates 


heat_df <- cor_wide %>% dplyr::select(Pathway, HA, MCI, AD)
heat_mat <- heat_df %>% tibble::column_to_rownames("Pathway") %>% as.matrix()
heat_mat <- apply(heat_mat, 2, as.numeric)
rownames(heat_mat) <- heat_df$Pathway


cat("\n[VISUAL] Rendering and saving figures...\n")
all_celltypes <- levels(factor(hip_1$Celltype))
seurat_palette <- scales::hue_pal()(length(all_celltypes))
names(seurat_palette) <- all_celltypes
ca_color_fixed <- seurat_palette["Astrocytes"]


# ==============================================================================
unified_theme <- theme(
  plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 10)),
  axis.title = element_text(face = "bold", size = 11, color = "black"),
  axis.text = element_text(size = 10, color = "black"),
  panel.grid = element_blank(),
  plot.tag = element_text(face = "bold", size = 16, color = "black")
)
# ==============================================================================
#  ROW 1 MODULES: CELLULAR ATLAS & TARGET TRACKING (PANELS A–D)
# ==============================================================================
cat("[MODULE] Constructing Row 1 (Panels A–D)...\n")
all_celltypes <- levels(factor(hip_1$Celltype))
seurat_palette <- scales::hue_pal()(length(all_celltypes))
names(seurat_palette) <- all_celltypes
ca_color_fixed <- seurat_palette["Astrocytes"]
DefaultAssay(Astrocytes_sub) <- "RNA"
fig4A <- DimPlot(Astrocytes_sub, group.by = "Celltype", reduction = "umap", cols = ca_color_fixed) + 
  labs(title = "A. Astrocytes Subset") +
  theme(legend.position = "none") + 
  unified_theme
ggsave(file="fig4A_Astrocytes Subset.pdf", fig4A, width = 8,height = 6)
ggsave(file="fig4A_Astrocytes Subset.png", fig4A, width = 8,height = 6,dpi = 600)

df_nr <- data.frame(NRXN1 = FetchData(Astrocytes_sub, vars = nrxn1_ens)[,1], 
                    Group = Astrocytes_sub$Group)
df_nr$Group <- factor(df_nr$Group, levels = c("HA", "MCI", "AD"))
df_nr %>% group_by(Group) %>% 
  summarise(
    Median_NR = median(NRXN1),
    Mean_NR = mean(NRXN1),
    SD_NR = sd(NRXN1),
    n = n())
# A tibble: 3 × 5
#   Group Median_NR Mean_NR SD_NR     n
#   <fct>     <dbl>   <dbl> <dbl> <int>
# 1 HA         3.60    3.01  1.54  2089
# 2 MCI        3.42    3.02  1.31  1249
# 3 AD         2.98    2.33  1.75  2095


id_rename_vec <- c(nrxn1_ens[1], syn2_ens[1], far2_ens[1], trim36_ens[1])
names(id_rename_vec) <- c("NRXN1", "SYN2", "FAR2","TRIM36")
DefaultAssay(Astrocytes_sub) <- "RNA"
p1d_df <- FetchData(Astrocytes_sub, vars = c(id_rename_vec, "Group"))
p1d_long <- p1d_df %>%
  tidyr::pivot_longer(cols = all_of(as.character(id_rename_vec)), 
                      names_to = "Ensembl_ID", 
                      values_to = "Expression") %>%
  mutate(Gene = names(id_rename_vec)[match(Ensembl_ID, id_rename_vec)]) %>%
  mutate(Gene = factor(Gene, levels = c("NRXN1", "SYN2",  "FAR2","TRIM36"))) %>% 
  mutate(Group = factor(Group, levels = c("HA", "MCI", "AD"))) 
p1d_long <- p1d_long[!is.na(p1d_long$Group), ]
fig4B <- ggplot(p1d_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", alpha = 0.55, color = NA, adjust = 1.2) +
  geom_boxplot(width = 0.14, color = "#2B2D42", outlier.shape = NA, alpha = 0.85, fill = "white", linewidth = 0.4) +
  facet_wrap(~Gene, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = my_colors) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14), 
    legend.position = "none", 
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F4F5F7", color = NA),
    strip.text = element_text(face = "bold.italic", size = 12, color = "#1A1A1A"),
    panel.spacing = unit(0.6, "lines"),
    axis.line = element_line(linewidth = 0.4)) +
  labs(
    title = "B. Gene Expression Levels in Astrocytes",
    y = "Log-Normalized Expression Intensity")+ 
  unified_theme

# 保存
ggsave(file="fig4B_Astrocytes Gene Expression.pdf", fig4B, width = 8,height = 6)
ggsave(file="fig4B_Astrocytes Gene Expression.png", fig4B, width = 8,height = 6,dpi = 600)


fig4C <- FeaturePlot(Astrocytes_sub, features = nrxn1_ens) +
  scale_colour_gradientn(colors = expr_cols) +
  labs(title = "C. NRXN1 Spatial Expression Profile", x = "umap_1", y = "umap_2") +
  theme(legend.position = "right") +
  unified_theme
ggsave(file="fig4C_Astrocytes NRXN1 Spatial Expression Profile.pdf", fig4C, width = 8,height = 6)
ggsave(file="fig4C_Astrocytes NRXN1 Spatial Expression Profile.png", fig4C, width = 8,height = 6,dpi = 600)

df_nr <- data.frame(NRXN1 = FetchData(Astrocytes_sub, vars = nrxn1_ens)[,1], Group = Astrocytes_sub$Group) 
df_nr$Group <- factor(df_nr$Group, levels = c("HA", "MCI", "AD")) 
df_nr %>% group_by(Group) %>% 
  summarise( Median_NR = median(NRXN1), 
             Mean_NR = mean(NRXN1), 
             SD_NR = sd(NRXN1), n = n()) 
# # A tibble: 3 × 5
# Group Median_NR Mean_NR SD_NR     n
# <fct>     <dbl>   <dbl> <dbl> <int>
# 1 HA         3.61    3.04  1.51  6980
# 2 MCI        3.49    3.09  1.26  4183
# 3 AD         3.02    2.37  1.73  6950

cat("[MODULE] Constructing Row 2 (Panels D,E)...\n")
top20_sorted <- top20 %>% dplyr::arrange(desc(Delta_AD_HA))
top_pathways_sorted <- top20_sorted$Pathway
heat_mat_sub <- heat_mat[top_pathways_sorted, c("HA", "MCI", "AD"), drop = FALSE]
p_heat <- pheatmap(
  heat_mat_sub, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  border_color = "white", cellwidth = 32, cellheight = 14, 
  fontsize_row = 9, fontsize_col = 11, angle_col = 0, 
  main = "D. NRXN1-Associated Profiles",
  silent = TRUE)
fig4D <- patchwork::wrap_elements(p_heat$gtable) # 严密封装进入拼图
# 保存
ggsave(file="fig4D_Astrocytes NRXN1-Associated Profiles.pdf", fig4D, width = 8,height = 6)
ggsave(file="fig4D_Astrocytes NRXN1-Associated Profiles.png", fig4D, width = 8,height = 6,dpi = 600)

plot_df_4 <- rbind(head(cor_wide, 10), tail(cor_wide, 10))
plot_df_4$Direction <- ifelse(plot_df_4$Delta_AD_HA > 0, "Accelerated in AD", "Suppressed in AD")
fig4E <- ggplot(plot_df_4, aes(x = reorder(Pathway, Delta_AD_HA), y = Delta_AD_HA, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75, color = "black", size = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Accelerated in AD" = "#C73E1D", "Suppressed in AD" = "#2E86AB")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "E. Metabolic Rewiring Hierarchy", 
       x = "Metabolic Pathway", y = "Delta Spearman Rho (AD - HA)", fill = "Rewiring Status") +
  theme_classic() +
  unified_theme +
  theme(legend.position = "bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 8))
ggsave(file="fig4E_Astrocytes NRXN1-Associated Metabolic Rewiring Hierarchy.pdf", fig4E, width = 8,height = 6)
ggsave(file="fig4E_Astrocytes NRXN1-Associated Metabolic Rewiring Hierarchy.png", fig4E, width = 8,height = 6,dpi = 600)


fig4F <- ggplot(cor_wide, aes(x = HA, y = AD)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", lwd = 0.6) +
  geom_hline(yintercept = 0, color = "grey80", size = 0.4) +
  geom_vline(xintercept = 0, color = "grey80", size = 0.4) +
  geom_point(aes(size = abs(Delta_AD_HA), fill = Delta_AD_HA), shape = 21, alpha = 0.85, color = "black") +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#C73E1D", midpoint = 0, name = "Δ (AD - HA)") +
  scale_size_continuous(range = c(2.5, 6.5), guide = "none") + 
  geom_text_repel(data = top20, aes(label = Pathway), size = 3.5, fontface = "italic", 
                  max.overlaps = 15, box.padding = 0.4, segment.color = "grey30") +
  labs(title = "F. Pathway Rewiring Dynamics", 
       x = "Spearman Rho in HA", y = "Spearman Rho in AD") +
  theme_classic() +
  unified_theme + 
  theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 12))
ggsave(file="fig4F_Astrocytes NRXN1-Associated Metabolic Rewiring Hierarchy.pdf", fig4F, width = 10,height = 6)
ggsave(file="fig4F_Astrocytes NRXN1-Associated Metabolic Rewiring Hierarchy.png", fig4F, width = 10,height = 6,dpi = 600)
 
ether_auc <- GetAssayData(Astrocytes_sub, assay = "KEGG_AUC")["Ether lipid metabolism", ]
df_meta <- data.frame(
  NRXN1 = FetchData(Astrocytes_sub, vars = nrxn1_ens)[,1], 
  Ether = ether_auc, Group = factor(Astrocytes_sub$Group, levels = c("HA", "MCI", "AD")))
fig4G <- ggplot(df_meta, aes(x = nrxn1, y = Ether, color = Group)) +
  geom_point(alpha = 0.15, size = 0.4) +
  geom_smooth(method = "lm", aes(fill = Group), lwd = 1.0, alpha = 0.1) +
  scale_color_manual(values = my_colors, name = "Clinical Stage") +
  scale_fill_manual(values = my_colors, name = "Clinical Stage") +
  labs(title = "G. NRXN1 & Ether Lipid Coupling", 
       x = "NRXN1 Expression Level", y = "Ether Lipid Metabolism (AUC Score)") +
  theme_classic() +    
  unified_theme +      
  theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 8))
ggsave(file="fig4G_Astrocytes NRXN1 & Ether Lipid Coupling.pdf", fig4G, width = 8,height = 6)
ggsave(file="fig4G_Astrocytes NRXN1 & Ether Lipid Coupling.png", fig4G, width = 8,height = 6,dpi = 600)


trend_df <- as.data.frame(trend_result$emtrends) %>%
  transmute(
    Group = factor(Group, levels = c("HA", "MCI", "AD")),
    Estimate = NRXN1.trend,
    SE = SE,
    Lower_CI = lower.CL,
    Upper_CI = upper.CL)
print(trend_df)

fig4H <- ggplot(trend_df, aes(x = Estimate, y = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey45") +
  geom_errorbarh(
    aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.16, linewidth = 0.8, color = "black") +
  geom_point(aes(fill = Group), shape = 21, size = 4.5, stroke = 0.7, color = "black") +
  scale_fill_manual(values = my_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title = "H. Adjusted Stage-Specific NRXN1 Slopes for Ether Lipid Activity",
    x = "Adjusted NRXN1 Slope\n(95% Confidence Interval)",
    y = NULL) +
  theme_classic(base_size = 12) +
  unified_theme +
  theme(axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 20, 10, 10))
fig4H
ggsave(file="fig4H_Astrocytes NRXN1 & Ether Lipid Estimated marginal trends.pdf", fig4H, width = 7,height = 6)
ggsave(file="fig4H_Astrocytes NRXN1 & Ether Lipid Estimated marginal trends.png", fig4H, width = 7,height = 6,dpi = 600)




# Weak association between FAR2 and NRXN1 does not explain ether lipid remodeling
library(tidyverse)
library(patchwork)
library(ggpubr)
pS6A <- ggplot(df2, aes(NRXN1, FAR2)) +
  geom_point(alpha=0.15, size=0.5, color="#4C78A8") +
  geom_smooth(method="lm", color="black", fill="grey80") +
  #stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
           label= "Pearson r = -0.024",
           size=3.5) +
  labs(title="NRXN1–FAR2: Weak Transcriptional Coupling", x="NRXN1 Expression", y="FAR2 Expression") +
  unified_theme

library(dplyr)
set.seed(123)
boot_cor_group <- df2 %>%
  group_by(Group) %>%
  group_modify(~{
    dat <- .x
    B <- 1000
    boot_rho <- numeric(B)
    for(i in 1:B){
      idx <- sample(1:nrow(dat), replace = TRUE)
      boot_rho[i] <- cor(dat$NRXN1[idx], dat$FAR2[idx], method = "spearman")}
    data.frame(rho = cor(dat$NRXN1, dat$FAR2, method="spearman"), 
               lower = quantile(boot_rho, 0.025),
               upper = quantile(boot_rho, 0.975))})
boot_cor_group

# # A tibble: 3 × 4
# # Groups:   Group [3]
# Group     rho   lower    upper
# <fct>   <dbl>   <dbl>    <dbl>
# 1 AD    -0.0260 -0.0505 -0.00347
# 2 MCI   -0.0249 -0.0552  0.00516
# 3 HA    -0.0235 -0.0450 -0.00147

cor_group <- data.frame(
  Group=c("HA","MCI","AD"),
  rho=c(-0.0235,-0.0249,-0.0260),
  lower=c(-0.0450,-0.0552,-0.0505),
  upper=c(-0.00147,0.00516,-0.00347)
)
cor_group$Group <- factor(cor_group$Group, levels = c("HA", "MCI", "AD"))
pS6B <- ggplot(cor_group, aes(Group, rho, fill=Group)) +
  geom_col(width=0.6, color="black") +
  geom_errorbar(aes(ymin=lower,ymax=upper), width=0.15) +
  geom_hline(yintercept=0, linetype=2) +
  scale_fill_manual(values=my_colors) +
  labs(title="Stage-consistent Minimal NRXN1–FAR2 Correlation", y="Spearman rho", x=NULL)+
  unified_theme +
  theme(legend.position="none")

pS6C <- ggplot(df2, aes(FAR2, Ether_Lipid)) +
  geom_point(alpha=0.15, size=0.5, color="#C73E1D") +
  geom_smooth(method="loess", color="black", fill="grey80") +
  geom_smooth(method="lm", linetype=2, color="darkred") +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
           label="Limited association with ether lipid activity",
           size=3.5) +
  labs(title="Modest Association of FAR2 with Ether Lipid Activity",
       x="FAR2 Expression", y="Ether Lipid AUC Score")+
  unified_theme

set.seed(123)
B <- 1000
boot_spearman <- numeric(B)
for(i in 1:B){
  idx <- sample(1:nrow(df2), replace = TRUE)
  boot_spearman[i] <- cor(df2$NRXN1[idx], df2$FAR2[idx], method = "spearman")}
quantile(boot_spearman, c(0.025,0.975))
#           2.5%        97.5% 
#   -0.032659338 -0.005167548
set.seed(123)
B <- 1000
boot_partial <- numeric(B)
for(i in 1:B){
  idx <- sample(1:nrow(df), replace = TRUE)
  tmp <- df[idx,]
  boot_partial[i] <- pcor.test(tmp$FAR2, tmp$Ether_Lipid, tmp[,c("NRXN1","SYN2","TRIM36")])$estimate}
quantile(boot_partial, c(0.025,0.975))
#            2.5%         97.5% 
#   -0.0005703658  0.0309103064 
forest_df <- data.frame(
  Variable=c(
    "Pearson\n(NRXN1-FAR2)",    ### cor.test(df2$NRXN1, df2$FAR2, method="pearson")
    "Spearman\n(overall NRXN1-FAR2)",   ### cor.test(df2$NRXN1, df2$FAR2, method="spearman")
    "Partial correlation\n(FAR2~Ether lipid)",   ### pcor.test(df$FAR2, df$Ether_Lipid, df[,c("NRXN1","SYN2","TRIM36")])
    "Linear Model\n(FAR2-Ether Lipid)"),         ### model_far2 <- lm(Ether_Lipid ~ FAR2 + NRXN1 + SYN2 + TRIM36 + Group, data=df)
  Effect=c(-0.02357, -0.01905, 0.01492, 0.00105),
  Lower=c(-0.03812, -0.03266, -0.00057, 0.00005),
  Upper=c(-0.00901, -0.00517, 0.03091, 0.00205))
pS6D <- ggplot(forest_df, aes(Effect, reorder(Variable, Effect))) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=0.2) +
  geom_point(aes(color=abs(Effect)), size=3) +
  scale_color_gradient(low="grey70", high="#C73E1D") +
  labs(title="Effect Size Estimates Across NRXN1–FAR2 and FAR2–Ether Lipid Analyses", x="Association Estimate", y=NULL) +
  unified_theme +
  theme(legend.position="none")

set.seed(123)
B <- 1000
boot_r <- numeric(B)
# for(i in 1:B){
#   idx <- sample(1:nrow(df2), replace = TRUE)
#   boot_r[i] <- cor(df2$NRXN1[idx], df2$FAR2[idx], method = "pearson")}
for(i in 1:B){
  idx <- sample(
    1:nrow(df2),
    replace=TRUE)
  boot_r[i] <- cor(
    df2$NRXN1[idx],
    df2$FAR2[idx],
    method="pearson")}
quantile(boot_r, c(0.025,0.5,0.975))
boot_df <- data.frame(r = boot_r)
pS6E <- ggplot(boot_df, aes(r)) +
  geom_histogram(bins = 40, fill="#4C78A8", color="white") +
  geom_vline(xintercept=cor(df2$NRXN1,df2$FAR2), linetype=2, color="red") +
  labs(title="Bootstrap Stability of NRXN1–FAR2 Correlation",
       x="Pearson r (bootstrap)", y="Frequency") +
  unified_theme
Supp_Fig6 <- (pS6A | pS6B) / (pS6C | pS6D) / pS6E +
  plot_layout(heights = c(1, 1.15, 0.9)) + 
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, color = "black"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )
ggsave("Supplementary_Figure_S6_FAR2.pdf", Supp_Fig6, width=12, height=12)
ggsave("Supplementary_Figure_S6_FAR2.png", Supp_Fig6, width=12, height=12, dpi=600)


