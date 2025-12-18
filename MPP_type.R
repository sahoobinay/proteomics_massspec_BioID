library(ggplot2)
library(dplyr)
# Marker genes for each MPP subtype
mpp2_genes <- c("Gata1", "Klf1", "Mpl", "Fli1")    # Meg/Ery biased
mpp3_genes <- c("Cebpa", "Elane", "Mpo", "Lyz2")   # Myeloid biased
mpp4_genes <- c("Flt3", "Il7r", "Rag1", "Dntt")    # Lymphoid biased

seurat_obj <- AddModuleScore(seurat_obj, features = list(mpp2_genes), name = "MPP2_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(mpp3_genes), name = "MPP3_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(mpp4_genes), name = "MPP4_Score")

# Extract scores
scores <- seurat_obj@meta.data[, c("MPP2_Score1", "MPP3_Score1", "MPP4_Score1")]

# Assign based on max score
seurat_obj$MPP_type <- apply(scores, 1, function(x) {
  c("MPP2", "MPP3", "MPP4")[which.max(x)]
})

# Stacked bar plot

table_df <- table(seurat_obj$sample, seurat_obj$MPP_type) %>% prop.table(1) %>% as.data.frame()
colnames(table_df) <- c("Sample", "MPP_Type", "Proportion")

ggplot(table_df, aes(x = Sample, y = Proportion, fill = MPP_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "MPP Subtype Composition by Sample")

# Subset to MPP3 cells and compare gene expression
mpp3_obj <- subset(seurat_obj, subset = MPP_type == "MPP3")
Idents(mpp3_obj) <- "sample"
DefaultAssay(mpp3_obj) <- "RNA"
mpp3_diff <- FindMarkers(mpp3_obj, ident.1 = "KI", ident.2 = "WT")

# Repeat for MPP2 / MPP4 if needed

DimPlot(seurat_obj, group.by = "MPP_type", label = TRUE)

