# Load R package
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

######------Figure 2B------######
hip_subset <- subset(hip, subset = group %in% c("AD", "Aging", "Adult"))
new_order <- c("AD", "Aging", "Adult")
hip_subset$group <- factor(hip_subset$group,levels = new_order, ordered = T)
prop_df1 <- table(hip_subset$celltype,hip_subset$group,hip_subset$sample) %>% reshape2::melt()
colnames(prop_df1) <- c("Cluster","Group","Sample","Number")
prop_df1$Cluster <- factor(prop_df1$Cluster)
prop_df1$Group   <- factor(prop_df1$Group)
# Calculate cell proportion by sample
prop_df1 <- prop_df1 %>% group_by(Group, Sample) %>% filter(sum(Number) > 0) %>% mutate(Proportion = Number / sum(Number)) %>% ungroup()
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_prop = sum(Proportion))
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_cells = sum(Number), .groups = "drop") %>% filter(total_cells == 0)
# One column for each group (average of the sample is calculated first).
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
ggsave(file="Figure 2B.pdf", width=6, height=3)
ggsave(file="Figure 2B.png", width=6, height=3)

######------Supplementary Figure 5------######
hip_subset <- subset(hip, subset = group %in% c("AD", "Aging", "Adult"))
new_order <- c("Adult", "Aging","AD" )
hip_subset$group <- factor(hip_subset$group,levels = new_order, ordered = T)
prop_df1 <- table(hip_subset$celltype,hip_subset$group,hip_subset$sample) %>% reshape2::melt()
colnames(prop_df1) <- c("Cluster","Group","Sample","Number")
prop_df1$Cluster <- factor(prop_df1$Cluster)
prop_df1$Group   <- factor(prop_df1$Group)

# Calculate cell proportion by sample
prop_df1 <- prop_df1 %>% group_by(Group, Sample) %>% filter(sum(Number) > 0) %>% mutate(Proportion = Number / sum(Number)) %>% ungroup()
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_prop = sum(Proportion))
prop_df1 %>% group_by(Group, Sample) %>% summarise(total_cells = sum(Number), .groups = "drop") %>% filter(total_cells == 0)

# Manually calculate the Kruskal p value for each panel (and generate ns).
kw_label <- prop_df1 %>% group_by(Cluster) %>%
  summarise(p = kruskal.test(Proportion ~ Group)$p.value, .groups = "drop") %>%
  mutate(label = paste0("Kruskal P = ", signif(p, 3)))

### drawing---------------
kw_label <- prop_df1 %>% group_by(Cluster) %>%
  summarise(p = kruskal.test(Proportion ~ Group)$p.value,
            y_pos = max(Proportion, na.rm = TRUE) * 1.05) %>%
  mutate(label = paste0("Kruskal p = ", signif(p, 3)))
ggplot(prop_df1, aes(x = Group, y = Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8) +
  facet_wrap(~ Cluster, scales = "free_y") +
  ## Kruskal–Wallis p
  geom_text(data = kw_label,
            aes(x = Inf, y = y_pos, label = label),
            inherit.aes = FALSE, color = "red",
            fontface = "bold", size = 3, hjust = 1.1) +
  ## pairwise
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", hide.ns = TRUE) +
  labs(x = "Group", y = "Cell proportion") +
  theme_classic()
ggsave(file="Supplementary Figure 5.pdf", width=9, height=6)
ggsave(file="Supplementary Figure 5.png", width=9, height=6)
