# 🧲 BioID Proteomics: Proximity Labeling Interactome Analysis using DEP

This repository contains the full R-based workflow and results for analyzing BioID proximity labeling proteomics data using the [DEP](https://bioconductor.org/packages/release/bioc/html/DEP.html) package. The analysis focuses on identifying high-confidence interacting partners of a bait protein through differential enrichment analysis of streptavidin pulldown samples.

---

## 🧬 Project Overview

**BioID (proximity-dependent biotinylation)** allows identification of proximal or interacting proteins in living cells. In this study, we used a **BioID2-tagged bait construct** expressed in [your cell line or system] and compared it to a **BirA* control** to identify specific interactors.

Mass spectrometry was performed after streptavidin pulldown, and data was analyzed using the DEP pipeline to:
- Normalize LFQ intensities
- Filter and impute missing values
- Perform differential enrichment analysis (bait vs control)
- Visualize interactors via volcano plots and heatmaps

---

## 📁 Repository Structure

```
bioid-proteomics/
├── README.md
├── LICENSE
├── data/
│   └── proteinGroups.txt        # MaxQuant output (not uploaded if proprietary)
├── scripts/
│   └── DEP.R                    # DEP pipeline for differential enrichment
├── figures/
│   ├── volcano_plot.pdf         # Volcano plot of enriched interactors
│   ├── heatmap.pdf              # Heatmap of top variable or enriched proteins
│   └── pca_plot.pdf             # PCA of bait vs control samples
```

---

## 🧪 Experimental Design

| Sample Type     | Condition   | Replicates |
|----------------|-------------|------------|
| Control         | BirA* only  | n = X      |
| Experimental    | BioID2-Bait | n = Y      |

Raw `.raw` files were processed using **MaxQuant**, and `proteinGroups.txt` was used as input to DEP.

---

## 🔁 Running the Pipeline

### 1. 📦 Install R Dependencies
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DEP", "SummarizedExperiment", "imputeLCMD", "ComplexHeatmap", "org.Hs.eg.db"))
install.packages(c("ggplot2", "readr", "dplyr", "tidyr", "stringr", "ggrepel"))
```

### 2. ▶️ Run the Analysis
```r
source("scripts/DEP.R")
```

### 3. 🖼 Outputs
The pipeline generates publication-ready:
- **Volcano plot**: Log2 fold change vs –log10 adjusted p-value
- **Heatmap**: Top enriched interactors
- **PCA**: Sample-level separation

---

## 📊 Example Figure Descriptions

- `volcano_plot.pdf`: Highlights significantly enriched interactors (FDR < 0.05, log2FC > 1)
- `heatmap.pdf`: Shows top 30 proteins by variance across samples
- `pca_plot.pdf`: Confirms experimental vs control separation

---

## 📜 License

This project is licensed under the MIT License — see the `LICENSE` file for details.

---

## 🙋 Acknowledgments & Contact

For questions, collaboration, or attribution, contact:

**Dr. Binay K. Sahoo**  
Postdoctoral Fellow, Stanford School of Medicine  
GitHub: [sahoobinay](https://github.com/sahoobinay)

---

## 🔗 Citation

Please cite this repository or associated publication/preprint if used in your work. Citation coming soon.
