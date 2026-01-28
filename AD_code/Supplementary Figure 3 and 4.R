######------Supplementary Figure 3 and 4------######
# Load R package
library(Seurat)
library(ggplot2)
library(cowplot)
library(viridis)
library(tibble)
library(tidyr)

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
saveRDS(hip_unfiltered, file = "seurat_objects/all_unfiltered.RDS")
qc_std_plot(hip_unfiltered)
ggsave("Supplementary Figure 3.png", path = "D:/AD/HIP", width = 30, height = 30, units = "cm")
## After filtering
hip <- subset(hip_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)
qc_std_plot(hip)
ggsave("Supplementary Figure 4.png", path = "D:/AD/HIP", width = 30, height = 30, units = "cm")

