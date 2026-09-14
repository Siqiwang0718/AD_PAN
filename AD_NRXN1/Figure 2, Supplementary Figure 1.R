
library(glmnet)
library(pROC)
library(ggplot2)
hub_gene_list <- rownames(hub_genes)
x <- as.matrix(datExpr[, colnames(datExpr) %in% hub_gene_list])
x <- scale(x)
cli <- cli[rownames(x), , drop = FALSE]
y <- factor(ifelse(cli$group == "Control", "Control", "AD"), levels = c("Control", "AD"))

set.seed(123)
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = nrow(x))

fit <- glmnet(x, y, family = "binomial", alpha = 1)
pdf("Figure 2A.pdf", width = 6, height = 5)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

pdf("Figure 2B.pdf", width = 6, height = 5)
plot(cv_fit)
dev.off()

coefs <- coef(cv_fit, s = "lambda.min")
final_signature_genes <- rownames(coefs)[which(coefs != 0)][-1]
print(final_signature_genes)


###############################
# Bootstrap internal validation
###############################
set.seed(123)
n_boot <- 1000
bootstrap_selected_genes <- vector("list", n_boot)
bootstrap_auc <- numeric(n_boot)

for(i in 1:n_boot){cat("Bootstrap iteration:", i, "\n")
  # stratified bootstrap
  control_index <- which(y=="Control")
  ad_index <- which(y=="AD")
  
  boot_control <- sample(control_index, length(control_index), replace=TRUE)
  boot_ad <- sample(ad_index, length(ad_index), replace=TRUE)
  boot_index <- c(boot_control, boot_ad)
  
  x_boot <- x[boot_index,]
  y_boot <- y[boot_index]
  cv_boot <- cv.glmnet(x_boot, y_boot, family="binomial", alpha=1, nfolds=5)
  coef_boot <- coef(cv_boot, s="lambda.min")
  
  selected_boot <- rownames(coef_boot)[which(as.numeric(coef_boot)!=0)]
  selected_boot <- selected_boot[selected_boot!="(Intercept)"]
  bootstrap_selected_genes[[i]] <- selected_boot
  prob_boot <- as.numeric(predict(cv_boot, newx=x_boot, s="lambda.min", type="response"))
  roc_boot <- roc(y_boot, prob_boot, quiet=TRUE)
  bootstrap_auc[i] <- as.numeric(auc(roc_boot))
}

gene_frequency <- sort(table(unlist(bootstrap_selected_genes)) / n_boot, decreasing = TRUE)
gene_frequency_df <- data.frame(Gene = names(gene_frequency), Selection_frequency = as.numeric(gene_frequency))
print(gene_frequency_df)

###############################
# 2. Bootstrap AUC
###############################
bootstrap_auc <- bootstrap_auc[!is.na(bootstrap_auc)]
auc_summary <- data.frame(
  Median_AUC = median(bootstrap_auc),
  Mean_AUC = mean(bootstrap_auc),
  Lower95 = quantile(bootstrap_auc,0.025),
  Upper95 = quantile(bootstrap_auc,0.975))
print(auc_summary)

###############################
# Supplementary Figure 1A
# Panel A: Bootstrap selection frequency
# 提取最终4个signature genes
library(dplyr)
signature_genes <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
panelA_df <- gene_frequency_df %>%
  filter(Gene %in% signature_genes) %>%
  mutate(Gene = factor(Gene, levels = c("NRXN1", "SYN2", "TRIM36", "FAR2")),
         Selection_frequency = Selection_frequency * 100)
print(panelA_df)

fig_panelA <- ggplot(panelA_df, aes(x = Gene, y = Selection_frequency)) +
  geom_col(width = 0.65, fill = "#4C72B0") +
  geom_text(aes(label = paste0(round(Selection_frequency,1), "%")),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0,0.12))) +
  labs(title = "A. Bootstrap Selection Stability of LASSO Genes", x = NULL, y = "Selection frequency (%)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0))
fig_panelA
ggsave("Supplementary_Figure_1A.pdf", fig_panelA, width = 5, height = 4)


prob <- as.numeric(predict(cv_fit, newx = x, s = "lambda.min", type = "response"))
roc_full <- roc(y, prob, quiet = TRUE)
auc_val <- round(auc(roc_full), 3)

pdf("Figure 2C.pdf", width = 5, height = 5)
plot(roc_full, col = "#2196F3", lwd = 3,
     main = paste0("Full Model ROC (AUC = ", auc_val, ")"),
     legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")
text(0.4, 0.2, paste0("AUC = ", auc_val), size = 5, font = 2)
dev.off()


tmp_coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
final_df <- data.frame(
  Gene = rownames(tmp_coefs)[which(tmp_coefs != 0)],
  Coefficient = tmp_coefs[which(tmp_coefs != 0)]
)
final_df <- final_df[-1, ]


library(pheatmap)
library(psych) # 用于计算相关性和显著性
my_cols <- c(
  Control = "#4DBBD5",   # 蓝色
  AD      = "#E64B35"    # 红色
)

final_genes <- c("NRXN1", "SYN2", "TRIM36", "FAR2")
gene_expr <- datExpr[, final_genes]

trait_cols <- c("mmse", "braak", "nft", "age", "pmi")
clinical_traits <- cli[rownames(gene_expr), trait_cols]


cor_res <- corr.test(gene_expr, clinical_traits, method = "spearman", adjust = "none")
cor_mat <- cor_res$r  
p_mat <- cor_res$p  

get_sig_text <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}
sig_text <- matrix(sapply(p_mat, get_sig_text), nrow = nrow(p_mat))

pdf("Figure 2H.pdf", width = 5, height = 4)
pheatmap(cor_mat,
         display_numbers = sig_text, # 在方格中显示星号
         number_color = "black",
         fontsize_number = 12,
         color = colorRampPalette(c("#2196F3", "white", "#F44336"))(100), # 蓝白红配色
         cluster_cols = FALSE, # 指标通常按临床逻辑排列，不建议聚类
         cluster_rows = TRUE,  # 基因可以聚类
         main = "Figure 2H",
         angle_col = 45)
dev.off()

library(rms)
final_genes <- c("NRXN1","SYN2","TRIM36","FAR2")
nomo_data <- data.frame(
  AD = ifelse(y=="AD",1,0),
  datExpr[, final_genes]
)
head(nomo_data)
dd <- datadist(nomo_data)
options(datadist="dd")
fit_nom <- lrm(AD ~ NRXN1 + SYN2 + TRIM36 + FAR2, data = nomo_data, x = TRUE, y = TRUE)
fit_nom
summary(fit_nom)
risk_score <- predict(fit_nom,  type = "lp")
nomo_data$RiskScore <- risk_score
nomo_data$Group <- factor(
  ifelse(nomo_data$AD == 1,"AD","Control"),
  levels = c("Control","AD")
)
plot_df <- nomo_data[order(nomo_data$RiskScore),]
plot_df$Sample <- 1:nrow(plot_df)
pdf("Figure 2D.pdf", width = 7,height = 4)
ggplot(plot_df, aes(Sample, RiskScore, color = Group)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(values = my_cols) +
  theme_classic(base_size = 14)
dev.off()

pdf("Figure 2E.pdf",width = 4,height = 5)
ggplot(nomo_data, aes(Group, RiskScore, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15,outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2) +
  scale_fill_manual(values = my_cols) +
  theme_classic(base_size = 14) +
  labs(x = NULL,y = "Risk Score")
dev.off()

library(pheatmap)
ann_colors <- list(Group = my_cols)
expr4 <- t(datExpr[, final_genes])
anno <- data.frame(Group = y)
rownames(anno) <- rownames(datExpr)
head(colnames(expr4))
head(rownames(anno))
all(colnames(expr4)==rownames(anno))
str(expr4)
mode(expr4)

#ann_colors <- list(Group = c(Control = "#4DBBD5", AD = "#E64B35"))
pdf("Figure 2G.pdf",width = 8,height = 4)
pheatmap(expr4, scale = "row",cluster_cols = FALSE, annotation_col = anno,
         annotation_colors = ann_colors, show_colnames = FALSE,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100))
dev.off()
#dev.list()
#while (!is.null(dev.list())) dev.off()


library(ggpubr)
library(reshape2)
plot_df <- datExpr[, final_genes]
plot_df$Group <- y
plot_df <- melt(plot_df, id.vars="Group")
p_exp <- ggplot(plot_df,aes(Group,value,fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15,outlier.shape = NA) +
  geom_jitter(width = 0.08,size = 1.5) +
  stat_compare_means() +
  #facet_wrap(~variable, nrow = 1, scales = "free_y") +
  facet_wrap(~variable,  ncol = 2,scales = "free_y")+
  scale_fill_manual(values = my_cols) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.8),
        axis.line = element_blank(),
        strip.background = element_rect(fill = "grey95",colour = "black"),
        strip.text = element_text(face = "bold"))
# 保存PDF
ggsave(file="Figure 2F.pdf", width = 8,height = 8)




library(limma)
library(GEOquery)
library(Biobase)
library(readxl)
library(GEOquery)
gset <- getGEO("GSE5281",
               destdir = ".", 
               GSEMatrix = TRUE,
               AnnotGPL = TRUE, 
               getGPL = TRUE)
eset <- gset[[1]]
exp<-exprs(eset)
cli<-pData(eset)
GPL<-fData(eset)
cli$age <- ifelse(is.na(cli$`age:ch1`) | cli$`age:ch1` == "", cli$`Age:ch1`,cli$`age:ch1`)
unique(cli$age)
cli$sex <- ifelse(is.na(cli$`sex:ch1`) | cli$`sex:ch1` == "", cli$`Sex:ch1`,cli$`sex:ch1`)
unique(cli$sex)
table(cli$sex)
#female    male  
#    58     103 
cli$group <- ifelse(is.na(cli$`disease state:ch1`) | cli$`disease state:ch1` == "", 
                    cli$`Disease State:ch1`,cli$`disease state:ch1`)
unique(cli$group)

unique(cli$age)
cli$age <- as.numeric(gsub("[^0-9.]", "", cli$age))
cli$age[grepl("days", cli$age)] <- NA
cli$age <- as.numeric(cli$age)
cli$organ <- cli$`Organ Region:ch1`
unique(cli$organ)
# 把 cli_hip 全表“标准化清洗”
cli[] <- lapply(cli, function(x) {
  if (is.character(x)) {
    x <- gsub("\u00A0", "", x)  
    x <- trimws(x)         
  }
  x
})
colnames(cli) <- gsub("\u00A0", "", colnames(cli))
colnames(cli) <- trimws(colnames(cli))
cli[] <- lapply(cli, function(x) {
  if (is.character(x)) {
    x <- gsub("[^[:alnum:] [:space:]\\.,'_-]", "", x)
    x <- trimws(x)
  }
  x
})
unique(cli$age)
unique(cli$sex)
sum(is.na(exp))
exp <- exp[complete.cases(exp), ]
max(exp)
if(max(exp)>30) exp=log2(exp+1)
max(exp)
gpl<-GPL[,c(1,3)]
gpl<-na.omit(gpl)
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),
                              stringsAsFactors = F)[,1]

exp<-as.data.frame(exp)
exp$ID<-rownames(exp)
exp_symbol<-merge(exp,gpl,by="ID", all = TRUE)
exp_2<-na.omit(exp_symbol)
table(duplicated(exp_2$`Gene symbol`))
exp_unique<-avereps(exp_2[,-c(1,ncol(exp_symbol))],
                    ID=exp_2$`Gene symbol`)

all(colnames(exp_unique) == rownames(cli))

unique(cli$characteristics_ch1.4)
cli_hip <- cli[grepl("hippocampus", cli$characteristics_ch1.4, ignore.case = TRUE), ]
unique(cli_hip$characteristics_ch1.4)
cli_hipr <- rownames(cli_hip[grepl("hippocampus", cli_hip$characteristics_ch1.4, ignore.case = TRUE), ])
cli_hipr
exp_hip <- exp_unique[, cli_hipr]
unique(cli_hip$age)
# 保存
write.csv(exp_hip, file = "exp_hip.csv", quote = FALSE)
write.csv(cli_hip, file = "cli_hip.csv", row.names = TRUE)

unique(cli$group)
group_list <- ifelse(cli_hip$group=="normal","Control","AD")
table(group_list)
class(cli_hip$age)

library(limma)
library(pROC)
design_val <- model.matrix(~ age, data = cli_hip)
exp_hip_corrected <- removeBatchEffect(
  exp_hip, 
  batch = factor(cli_hip$sex), 
  design = design_val
)

target_genes <- c("NRXN1", "SYN2", "TRIM36", "FAR2")

# 检查基因是否存在于 exp_hip_corrected 中
present_genes <- intersect(target_genes, rownames(exp_hip_corrected))
cat("当前选中的验证基因:", present_genes, "\n")


val_expr <- as.data.frame(t(exp_hip_corrected[present_genes, ]))
val_expr$Group <- factor(group_list, levels = c("Control", "AD"))

val_model <- glm(Group ~ ., data = val_expr, family = "binomial")
val_probs <- predict(val_model, type = "response")

# 计算 ROC
roc_val <- roc(val_expr$Group, val_probs, quiet = TRUE)
auc_val <- round(auc(roc_val), 3)


pdf("Figure 2I.pdf", width = 5, height = 5)
plot(roc_val, col = "#E91E63", lwd = 3,
     main = paste0("GSE5281 Validation (AUC = ", auc_val, ")"),
     legacy.axes = TRUE,
     print.auc = FALSE) 
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = paste0("AUC = ", auc_val), 
       col = "#E91E63", lwd = 3, bty = "n")
dev.off()
summary(val_model)

library(rms)
dd <- datadist(val_expr)
options(datadist = 'dd')
val_expr$Group_num <- ifelse(val_expr$Group == "AD", 1, 0)
lrm_model <- lrm(Group_num ~ NRXN1 + SYN2 + TRIM36 + FAR2, 
                 data = val_expr, x = TRUE, y = TRUE)


nom <- nomogram(lrm_model, 
                fun = function(x)1/(1+exp(-x)), 
                lp = FALSE, 
                fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                funlabel = "Probability of AD")


pdf("Figure 2J.pdf", width = 10, height = 5)
plot(nom, xfrac = 0.3)
dev.off()























