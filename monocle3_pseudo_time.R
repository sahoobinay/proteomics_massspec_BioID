library(monocle3)

root_cells <- WhichCells(seurat_obj, idents = "0")  # change "0" to whatever is LT-HSC
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)
#Save the cds object
rowData(cds)$gene_short_name <- rownames(cds)
plot_genes_in_pseudotime(cds[c("Gata2", "Cebpa", "Klf1"), ],
                         color_cells_by = "pseudotime")

deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
deg_pseudotime <- deg_pseudotime[order(deg_pseudotime$q_value), ]
head(deg_pseudotime)

top_genes <- rownames(deg_pseudotime)[1:10]

plot_genes_in_pseudotime(cds[top_genes, ],
                         color_cells_by = "pseudotime")


"Cbfa2t3" %in% rownames(cds)
plot_genes_in_pseudotime(cds["Luciferase", ],
                         color_cells_by = "pseudotime")

