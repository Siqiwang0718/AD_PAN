
rm(list = ls())
library(dplyr)
library(ggplot2)
library(susieR)
library(coloc)
library(VariantAnnotation)
library(gwasvcf)
library(TwoSampleMR)
library(patchwork)
library(Rfast)

plink <- "E:/software/plink/plink.exe"
bfile <- "D:/1000G_EUR/EUR"
out_prefix <- file.path(workingDir, "NRXN1_region")

vcf_file <- "D:/Biofo/MR/NRXN1/eqtl-a-ENSG00000179915.vcf.gz"
indexTabix(vcf_file, format = "vcf")
vcf <- readVcf(vcf_file, genome = "GRCh37")
dat <- gwasvcf::vcf_to_tibble(vcf) %>% as.data.frame()

dat_cis <- dat %>%
  filter(seqnames == "2", start >= 49145642, start <= 52259674) %>%
  filter(!is.na(ES) & !is.na(SE) & !is.na(LP) & !is.na(AF) & !is.na(SS))

exp <- dat_cis %>%
  mutate(
    SNP = rsid,
    beta = ES,
    se = SE,
    z = ES / SE,
    p = 10^(-LP),
    MAF = pmin(AF, 1 - AF),
    pos = start,
    chr = seqnames,
    SS = as.numeric(SS)
  ) %>%
  filter(MAF > 0.01) %>%
  distinct(SNP, .keep_all = TRUE)

N <- max(exp$SS, na.rm = TRUE)
if(is.na(N) || N <= 1) N <- 31684
cat("eQTL sample size N =", N, "\n")


system(paste(plink, "--bfile", bfile, "--chr 2", "--from-bp 49145642", "--to-bp 52259674",
             "--make-bed --out", out_prefix))
system(paste(plink, "--bfile", out_prefix, "--r square --out", out_prefix))

LD_matrix <- as.matrix(read.table(paste0(out_prefix, ".ld"), header = FALSE))
ld_bim <- read.table(paste0(out_prefix, ".bim"), header = FALSE)
colnames(ld_bim) <- c("CHR","SNP","CM","BP","A1","A2")
rownames(LD_matrix) <- ld_bim$SNP
colnames(LD_matrix) <- ld_bim$SNP

# =========================
common <- intersect(exp$SNP, ld_bim$SNP)
exp2 <- exp %>% filter(SNP %in% common)
exp2 <- exp2[match(common, exp2$SNP), ]
LD_matrix2 <- LD_matrix[common, common]
stopifnot(all(exp2$SNP == rownames(LD_matrix2)))
stopifnot(all(exp2$SNP == colnames(LD_matrix2)))

# =========================
R2 <- LD_matrix2 + diag(1e-3, nrow(LD_matrix2))
susie_fit <- susie_rss(z = exp2$z, R = R2, n = N, L = 3,
                       estimate_residual_variance = TRUE,
                       estimate_prior_variance = TRUE,
                       max_iter = 2000, tol = 1e-3, verbose = TRUE)
library(dplyr)
pip_df <- data.frame(SNP = colnames(susie_fit$alpha), pip = susie_fit$pip) %>%
  dplyr::left_join(exp2 %>% dplyr::select(SNP, pos), by = "SNP") %>%
  dplyr::mutate(cs = "None")

if(!is.null(susie_fit$sets$cs)){
  for(i in seq_along(susie_fit$sets$cs)){
    pip_df$cs[pip_df$SNP %in% susie_fit$sets$cs[[i]]] <- paste0("CS", i)
  }
}

lead_snp <- pip_df$SNP[which.max(pip_df$pip)]
lead_idx <- which(rownames(LD_matrix2) == lead_snp)
pip_df$r2 <- LD_matrix2[, lead_idx][match(pip_df$SNP, rownames(LD_matrix2))]^2

# 绘图
pA <- ggplot(pip_df, aes(x = pos, y = -log10(pip))) +
  geom_point(aes(color = r2, shape = cs), size = 2, alpha = 0.9) +
  geom_point(data = pip_df[pip_df$SNP == lead_snp,], color = "red", size = 3) +
  scale_color_gradient(low = "grey80", high = "firebrick") +
  theme_classic() +
  labs(title = "D. SuSiE fine-mapping", x = "Position", y = "-log10(PIP)")
ggsave("Figure 3D.pdf", pA, width = 8, height = 6, dpi = 600)

# MR + Harmonisation
# =========================
Sys.setenv(OPENGWAS_JWT = "")
outcome_dat <- extract_outcome_data(snps = exp2$SNP,
                                    outcomes = "ebi-a-GCST90027158",
                                    proxies = TRUE)
dat_harmonised <- harmonise_data(
  exposure_dat = exp2 %>%
    transmute(
      SNP = SNP,
      beta.exposure = beta,
      se.exposure = se,
      effect_allele.exposure = ALT,
      other_allele.exposure = REF,
      eaf.exposure = MAF,
      pval.exposure = p,
      samplesize.exposure = SS,
      id.exposure = "NRXN1_eQTL",
      exposure = "NRXN1_expression"
    ),
  outcome_dat = outcome_dat
) %>% filter(mr_keep == TRUE)
saveRDS(dat_harmonised)
# =========================
# 8. Lead SNP QC： R2 & F
# =========================
fallback_N <- 31684
lead_snp_data <- dat_harmonised %>%
  mutate(abs_z = abs(beta.exposure / se.exposure)) %>%
  filter(abs_z == max(abs_z, na.rm = TRUE)) %>%
  distinct(SNP, .keep_all = TRUE)

qc_res <- lead_snp_data %>%
  mutate(
    N_exp = ifelse(!is.na(samplesize.exposure) & samplesize.exposure > 1, samplesize.exposure, fallback_N),
    Lead_R2 = (beta.exposure^2) / (beta.exposure^2 + N_exp * (se.exposure^2)),
    F_statistic = (Lead_R2 * (N_exp - 2)) / (1 - Lead_R2)
  ) %>%
  dplyr::select(exposure, Total_SNPs_In_Region = samplesize.exposure, Lead_SNP = SNP, Lead_R2, F_statistic)

qc_res$Total_SNPs_In_Region <- nrow(dat_harmonised)
write.csv(qc_res, "MR_Instrument_Quality_QC_Fixed.csv", row.names = FALSE)
print(qc_res)
#          exposure Total_SNPs_In_Region   Lead_SNP    Lead_R2 F_statistic
#1 NRXN1_expression                 5939 rs13031157 0.01449366    209.7045

### Coloc
common <- intersect(exp2$SNP, outcome_dat$SNP)
exp_coloc <- exp2 %>% filter(SNP %in% common) %>% arrange(match(SNP, common))
colnames(outcome_dat)
gwas_coloc <- outcome_dat %>%
  filter(SNP %in% common) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  mutate(
    MAF = as.numeric(eaf.outcome),
    MAF = pmin(MAF, 1 - MAF)
  ) %>%
  filter(!is.na(MAF), MAF > 0 & MAF < 1) %>%
  arrange(match(SNP, common))

stopifnot(length(unique(gwas_coloc$SNP)) == nrow(gwas_coloc))
stopifnot(length(unique(exp_coloc$SNP)) == nrow(exp_coloc))
stopifnot(identical(exp_coloc$SNP, gwas_coloc$SNP))
n_cases <- 39106 + 46828 
n_controls <- 401577  
n_total <- n_cases + n_controls
s_prop <- n_cases / n_total
coloc_dat1 <- list(
  beta = exp_coloc$beta,
  varbeta = exp_coloc$se^2,
  snp = exp_coloc$SNP,
  MAF = exp_coloc$MAF,
  N = 31684,
  type = "quant"
)

coloc_dat2 <- list(
  beta = gwas_coloc$beta,
  varbeta = gwas_coloc$se^2,
  snp = gwas_coloc$SNP,
  MAF = gwas_coloc$MAF,
  N = n_total,
  type = "cc",
  s = s_prop
)

coloc_res <- coloc.abf(coloc_dat1, coloc_dat2)
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
#5.15e-41  9.10e-01  4.48e-42  7.93e-02  1.10e-02 
#[1] "PP abf for shared variant: 1.1%"
coloc_res$summary
#       nsnps    PP.H0.abf    PP.H1.abf    PP.H2.abf    PP.H3.abf    PP.H4.abf 
#6.095000e+03 5.146210e-41 9.097064e-01 4.484745e-42 7.926675e-02 1.102688e-02 


### MR分析
mr_res <- mr(dat_harmonised)
mr_res
# id.exposure         id.outcome                                      outcome         exposure                    method nsnp          b           se         pval
# 1  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression                  MR Egger 5939 0.01690481 4.757645e-03 3.835483e-04
# 2  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression           Weighted median 5939 0.02527973 6.092493e-03 3.334581e-05
# 3  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression Inverse variance weighted 5939 0.01467691 3.556737e-03 3.683080e-05
# 4  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression               Simple mode 5939 1.25608406 1.478822e+04 9.999322e-01
# 5  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression             Weighted mode 5939 1.25608406 1.633808e+04 9.999387e-01


# =========================================================
fallback_N <- 31684 

lead_snp_data <- dat_harmonised %>%
  mutate(abs_z = abs(beta.exposure / se.exposure)) %>%
  filter(abs_z == max(abs_z, na.rm = TRUE)) %>%
  distinct(SNP, .keep_all = TRUE)

qc_res <- lead_snp_data %>%
  mutate(
    N_exp = ifelse(!is.na(samplesize.exposure) & samplesize.exposure > 1, samplesize.exposure, fallback_N),
    Lead_R2 = (beta.exposure^2) / (beta.exposure^2 + N_exp * (se.exposure^2)),
    F_statistic = (Lead_R2 * (N_exp - 2)) / (1 - Lead_R2)
  ) %>%
  dplyr::select(exposure, Total_SNPs_In_Region = samplesize.exposure, Lead_SNP = SNP, Lead_R2, F_statistic)


qc_res$Total_SNPs_In_Region <- nrow(dat_harmonised)
print(qc_res)
# exposure Total_SNPs_In_Region   Lead_SNP    Lead_R2 F_statistic
# 1 NRXN1_expression                 5939 rs13031157 0.01449366    209.7045
write.csv(qc_res, file = "MR_Instrument_Quality_QC_Fixed.csv", row.names = FALSE)

mr_heterogeneity(dat_harmonised)
# id.exposure         id.outcome                                      outcome         exposure                    method        Q Q_df Q_pval
# 1  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression                  MR Egger 5015.989 5937      1
# 2  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression Inverse variance weighted 5016.488 5938      1

mr_pleiotropy_test(dat_harmonised)
# id.exposure         id.outcome                                      outcome         exposure egger_intercept           se      pval
# 1  NRXN1_eQTL ebi-a-GCST90027158 Alzheimer's disease || id:ebi-a-GCST90027158 NRXN1_expression   -0.0001361297 0.0001930758 0.4808014

mr_plot <- mr_scatter_plot(mr_res, dat_harmonised)


### Coloc plot
pp4 <- coloc_res$summary["PP.H4.abf"]
pos_map <- exp2[, c("SNP", "pos")]

dat_harmonised$pos <- pos_map$pos[
  match(dat_harmonised$SNP, pos_map$SNP)
]
"pos" %in% colnames(dat_harmonised)
pB <- ggplot(dat_harmonised, aes(x = pos, y = -log10(pval.outcome))) +
  geom_point(color = "steelblue", size = 1.3, alpha = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0("C. Regional GWAS association (PP.H4 = ", round(pp4, 4), ")"),
    x = "Genomic position",
    y = "-log10(P GWAS)"
  )
ggsave("Figure 3C.pdf", pB, width = 8, height = 6)


### MR panel
mr_plot <- mr_scatter_plot(mr_res, dat_harmonised)[[1]]

pC <- mr_plot +
  theme_classic(base_size = 12) +
  labs(title = "B. Mendelian Randomization")
ggsave("Figure 3B.pdf", pC, width = 8, height = 6)


########################################################################################################################################
library(ggplot2)
library(dplyr)

forest_data <- data.frame(
  Gene = c("NRXN1", "NRXN1", "NRXN1", "FAR2", "SYN2", "TRIM36"),
  Method = c("Inverse variance weighted", "Weighted median", "MR Egger", 
             "Inverse variance weighted", "Inverse variance weighted", "Inverse variance weighted"),
  b = c(0.01467691, 0.02527973, 0.01690481, 1.4064289, -0.4592985, 0.4814221),
  se = c(0.003556737, 0.006092493, 0.004757645, 0.8673437, 1.369664, 0.6661026)
)

plot_df <- forest_data %>%
  mutate(
    low_CI = b - 1.96 * se,
    up_CI = b + 1.96 * se,
    Gene = factor(Gene, levels = c("NRXN1", "FAR2", "SYN2", "TRIM36")),
    Method = factor(Method, levels = c("Inverse variance weighted", "Weighted median", "MR Egger")),
    Label = paste0(sprintf("%.3f", b), " (", sprintf("%.3f", low_CI), ", ", sprintf("%.3f", up_CI), ")")
  )

p_forest <- ggplot(plot_df, aes(x = b, y = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.7) +
  geom_errorbarh(aes(xmin = low_CI, xmax = up_CI, color = Gene), height = 0.15, size = 1.0) +
  geom_point(aes(fill = Gene), shape = 21, size = 3.5, color = "black") +
  geom_text(aes(x = max(plot_df$up_CI) * 1.3, label = Label), hjust = 0, size = 3.8, color = "#2D3748") +
  facet_grid(Gene ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("NRXN1" = "#C53030", "FAR2" = "#2B6CB0", "SYN2" = "#4A5568", "TRIM36" = "#718096")) +
  scale_fill_manual(values = c("NRXN1" = "#E53E3E", "FAR2" = "#4299E1", "SYN2" = "#718096", "TRIM36" = "#A0AEC0")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.7))) +
  theme_bw(base_size = 12) +
  labs(
    title = "Mendelian Randomization Analysis of 4 Candidate Hub Genes on AD Risk",
    x = "Causal Effect Size (Beta with 95% CI)",
    y = ""
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "#EDF2F7", color = "grey70"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 15)),
    axis.text.y = element_text(face = "bold", color = "#4A5568", size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )

# =========================================================
print(p_forest)
ggsave("Figure 3A.pdf", p_forest, width = 8.5, height = 5)


res_single <- mr_singlesnp(dat_harmonised)
p_funnel <- mr_funnel_plot(res_single)[[1]]
pleiotropy_p <- 0.4808
heterogeneity_p <- 1.0000
f_stat <- 209.70

# =========================================================
p_funnel_premium <- p_funnel +
  theme_bw(base_size = 12) +
  coord_cartesian(xlim = c(-15, 15)) + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#EDF2F7", size = 0.5),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(color = "#2D3748", size = 10),
    axis.text = element_text(color = "#718096", size = 9),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  ) +
  
  labs(
    title = "Funnel Plot for Genetically Instrumented NRXN1 on AD Risk",
    x = expression(paste("Causal Effect Estimate of Single SNP (", beta[IV], ")")),
    y = expression(paste("Instrument Precision (1 / ", SE[IV], ")")),
    color = "MR Method"
  ) +
  annotate(
    "text", x = 1.5, y = 18, 
    label = paste0(
      "Instrument Strength Check:\n",
      "  • Lead eQTL F-statistic = ", sprintf("%.2f", f_stat), " (F > 10)\n\n",
      "Horizontal Pleiotropy Check:\n",
      "  • Egger intercept P = ", sprintf("%.4f", pleiotropy_p), " (No pleiotropy)\n\n",
      "Heterogeneity Check:\n",
      "  • IVW Cochran's Q P = ", sprintf("%.4f", heterogeneity_p)
    ),
    hjust = 0, size = 3.6, color = "#4A5568", fontface = "plain"
  )
print(p_funnel_premium)
ggsave("Supplementary Figure 4A.pdf", p_funnel_premium, width = 8, height = 6)




# =========================================================
res_leaveoneout <- mr_leaveoneout(dat_harmonised)
overall_ivw_beta <- 0.01467691

if (!exists("res_leaveoneout")) {
  stop("run res_leaveoneout <- mr_leaveoneout(dat_harmonised) ")
}

loo_df <- as.data.frame(res_leaveoneout) %>%
  filter(SNP != "All") %>%
  mutate(
    low_CI = b - 1.96 * se,
    up_CI = b + 1.96 * se,
    y_idx = 1:n()
  )

p_loo_premium <- ggplot(loo_df, aes(x = b, y = y_idx)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "#CBD5E0", size = 0.6) +
  geom_errorbarh(aes(xmin = low_CI, xmax = up_CI), alpha = 0.03, color = "#718096", height = 0) +
  geom_point(alpha = 0.05, color = "#4A5568", size = 0.6) +
  geom_vline(xintercept = 0.01467691, linetype = "dashed", color = "#E53E3E", size = 0.8) +
  coord_cartesian(xlim = c(-0.005, 0.035)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 15)),
    axis.title = element_text(color = "#2D3748", size = 10),
    axis.text.x = element_text(color = "#4A5568", size = 9.5)
  ) +
  
  labs(
    title = "Leave-one-out Sensitivity Analysis for NRXN1 on AD Risk",
    x = expression(paste("Causal Effect Size after Re-estimating (", beta, " with 95% CI)")),
    y = "5,939 Instrumented Variants (Ranked by Genomic Position)"
  ) +
  
  annotate(
    "text", x = 0.018, y = nrow(loo_df) * 0.75,
    label = paste0(
      "Stability Evaluation:\n",
      "  • Total SNPs evaluated = 5,939\n",
      "  • Overall IVW Beta = 0.0147\n\n",
      "Conclusion:\n",
      "  No single eQTL variant disproportionately\n",
      "  drives or invalidates the primary causal\n",
      "  effect, confirming high robust stability."
    ),
    hjust = 0, size = 3.6, color = "#2D3748", fontface = "plain"
  )

print(p_loo_premium)
ggsave("Supplementary Figure 4B.pdf", p_loo_premium, width = 8, height = 6)






