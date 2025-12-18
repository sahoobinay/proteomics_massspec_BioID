library(ggplot2)


# Extract UMAP embeddings
umap_coords <- as.data.frame(Embeddings(seurat_obj, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Add lineage metadata
umap_df <- cbind(umap_coords, lineage = seurat_obj$lineage_identity)
# Compute lineage centroids
centroids <- umap_df %>%
  group_by(lineage) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# Optional: define differentiation order
lineage_order <- c("LT-HSC", "ST-HSC", "MPP", "CMP", "GMP", "MEP", "CLP")
centroids$lineage <- factor(centroids$lineage, levels = lineage_order)
centroids <- centroids[order(centroids$lineage), ]

# Plot UMAP + differentiation path
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = lineage)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_path(data = centroids, aes(x = UMAP_1, y = UMAP_2), 
            color = "black", arrow = arrow(length = unit(0.2, "inches")), 
            size = 1.2, inherit.aes = FALSE) +
  geom_text(data = centroids, 
            aes(x = UMAP_1, y = UMAP_2, label = lineage), 
            vjust = -1, inherit.aes = FALSE) +
  ggtitle("Hematopoietic Differentiation Trajectory on UMAP") +
  theme_minimal()
