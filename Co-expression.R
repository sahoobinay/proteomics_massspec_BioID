c("Cbfa2t3", "Glis2") %in% rownames(seurat_obj)


FeaturePlot(seurat_obj, features = "Cbfa2t3", split.by = "sample")
FeaturePlot(seurat_obj, features = "Glis2", split.by = "sample")
FeaturePlot(seurat_obj, features = "Luciferase", split.by = "sample")
FeaturePlot(seurat_obj, features = "mCherry", split.by = "sample")
# Get RNA expression matrix
DefaultAssay(seurat_obj) <- "RNA"
expr <- GetAssayData(seurat_obj, slot = "data")

# Threshold > 0 for detected expression
coexpressing_cells <- colnames(seurat_obj)[
  expr["Cbfa2t3", ] > 0 & expr["mCherry", ] > 0
 ]

seurat_obj$cbfa2t3_glis2 <- ifelse(colnames(seurat_obj) %in% coexpressing_cells, "Co-expressing", "Other")

DimPlot(seurat_obj, group.by = "cbfa2t3_glis2", cols = c("red", "gray"), split.by = "sample") +
  ggtitle("Glis2,Luciferasa,and mCherry Co-expression on UMAP")

seurat_obj$Cbfa2t3_Glis2_sum <- expr["Cbfa2t3", ] + expr["Glis2", ]
FeaturePlot(seurat_obj, features = "Cbfa2t3_Glis2_sum") +
  ggtitle("Joint Cbfa2t3 + Glis2 Expression")
