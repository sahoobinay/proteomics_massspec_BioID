# Adjust assay if needed
DefaultAssay(seurat_obj) <- "RNA"

# Set a threshold (adjust based on expression)
mcherry_pos_cells <- WhichCells(seurat_obj, expression = mCherry > 0.5)

seurat_obj$mcherry_status <- ifelse(colnames(seurat_obj) %in% mcherry_pos_cells, "mCherry_Pos", "mCherry_Neg")

# Combine with sample label if needed
seurat_obj$group <- paste(seurat_obj$sample, seurat_obj$mcherry_status, sep = "_")
table(seurat_obj$group)

Idents(seurat_obj) <- "group"

# Example: compare KI_mCherry_Pos vs WT_mCherry_Neg
markers <- FindMarkers(seurat_obj, ident.1 = "KI_mCherry_Pos", ident.2 = "WT_mCherry_Neg")
head(markers)

markers$gene <- rownames(markers)
markers$significance <- "NS"
markers$significance[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25] <- "Up"
markers$significance[markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.25] <- "Down"

library(ggplot2)
ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  labs(title = "KI mCherry⁺ vs WT mCherry⁻", x = "Log2 Fold Change", y = "-log10(adj p-value)") +
  theme_minimal()

# Top upregulated
top_up <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  head(500)

# Top downregulated
top_down <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.25) %>%
  arrange(avg_log2FC) %>%
  head(100)

top_up$gene
top_down$gene
write.csv(top_up, "top_upregulated_genes.csv", row.names = FALSE)
write.csv(top_down, "top_downregulated_genes.csv", row.names = FALSE)
FeaturePlot(seurat_obj, features = c("Cd48", "mCherry", "Flt3"), cols = c("lightgray", "red"))

VlnPlot(seurat_obj, features = c("Runx1", "Kit", "Sp1"), group.by = "group", pt.size = 0)

markers[c("Bcl2", "Glis2", "Cbfa2t3"), ]
