# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)

#----------------------------
# 1. Load Cell Ranger Outputs
#----------------------------
# WT
wt_data <- Read10X("data/Wt_Lin-/filtered_feature_bc_matrix")
wt <- CreateSeuratObject(counts = wt_data$`Gene Expression`, project = "WT")
wt[["ADT"]] <- CreateAssayObject(counts = wt_data$`Antibody Capture`)
wt$sample <- "WT"

# KI
ki_data <- Read10X("data/KI_Lin-/filtered_feature_bc_matrix")
ki <- CreateSeuratObject(counts = ki_data$`Gene Expression`, project = "KI")
ki[["ADT"]] <- CreateAssayObject(counts = ki_data$`Antibody Capture`)
ki$sample <- "KI"

#----------------------------
# 2. Merge and QC
#----------------------------
seurat_obj <- merge(wt, y = ki)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

#----------------------------
# 3. Normalize and Scale
#----------------------------
# RNA
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# ADT
seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")
seurat_obj <- ScaleData(seurat_obj, assay = "ADT")

#----------------------------
# 4. Clustering and UMAP
#----------------------------
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

#----------------------------
# 5. Visualization
#----------------------------
DimPlot(seurat_obj, group.by = "sample", reduction = "umap")
DimPlot(seurat_obj, group.by = "seurat_clusters")

# RNA Feature expression
FeaturePlot(seurat_obj, features = c("Glis2",  "Cbfa2t3"), split.by = "sample")

# ADT expression
DefaultAssay(seurat_obj) <- "ADT"

FeaturePlot(seurat_obj, features = c("Ms.CD56", "Ms.CD117"), split.by = "sample")
#----------------------------
# 6. Differential Expression (WT vs KI)
#----------------------------
# RNA
seurat_obj <- JoinLayers(seurat_obj)
Idents(seurat_obj) <- "sample"
DefaultAssay(seurat_obj) <- "RNA"
rna_markers <- FindMarkers(seurat_obj, ident.1 = "KI", ident.2 = "WT")

# ADT
DefaultAssay(seurat_obj) <- "ADT"
adt_markers <- FindMarkers(seurat_obj, ident.1 = "KI", ident.2 = "WT")

# Save
write.csv(rna_markers, "WT_vs_KI_Lin-_RNA_markers.csv")
write.csv(adt_markers, "WT_vs_KI_ADT_markers.csv")


