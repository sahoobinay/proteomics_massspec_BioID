# Marker genes for each subtype
LT_HSC_genes <- c("Procr", "Hlf", "Fgd5", "Hoxb5", "Gata2")
ST_HSC_genes <- c("Gata2", "Runx1", "Hlf")
mpp_genes <- c("Cd34", "Flt3", "Cebpa", "Cebpb")
cmp_genes <- c("Gata1", "Gata2")
gmp_genes <- c("Cebpa", "Cebpb", "Prtn3", "Elane")
mep_genes <- c("Gata1", "Klf1", "Epor")
clp_genes <- c("Il7r", "Rag1", "Dntt", "Ly6d")

seurat_obj <- AddModuleScore(seurat_obj, features = list(LT_HSC_genes), name = "LT_HSC_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(ST_HSC_genes), name = "ST_HSC_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(mpp_genes), name = "MPP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(cmp_genes), name = "CMP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(gmp_genes), name = "GMP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(mep_genes), name = "MEP_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = list(clp_genes), name = "CLP_Score")

# Extract scores
# Extract scores
score_cols <- c("LT_HSC_Score1", "ST_HSC_Score1", "MPP_Score1", 
                "CMP_Score1", "GMP_Score1", "MEP_Score1", "CLP_Score1")

# Assign dominant lineage based on max score
seurat_obj$lineage_identity <- apply(seurat_obj@meta.data[, score_cols], 1, function(x) {
  c("LT-HSC", "ST-HSC", "MPP", "CMP", "GMP", "MEP", "CLP")[which.max(x)]
})
# Assign dominant lineage based on max score
seurat_obj$lineage_identity <- apply(seurat_obj@meta.data[, score_cols], 1, function(x) {
  c("LT-HSC", "ST-HSC", "MPP", "CMP", "GMP", "MEP", "CLP")[which.max(x)]
})

DimPlot(seurat_obj, group.by = "lineage_identity", label = TRUE) +
  ggtitle("Hematopoietic Lineages by Gene Signature")

library(dplyr)
library(ggplot2)

table_df <- table(seurat_obj$sample, seurat_obj$lineage_identity) %>%
  prop.table(1) %>% as.data.frame()

colnames(table_df) <- c("Sample", "Lineage", "Proportion")

ggplot(table_df, aes(x = Sample, y = Proportion, fill = Lineage)) +
  geom_bar(stat = "identity") +
  labs(title = "Lineage Composition by Sample", y = "Proportion") +
  theme_minimal()
