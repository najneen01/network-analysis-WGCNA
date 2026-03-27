# network-analysis-WGCNA-GSE247791

This repository contains a complete Weighted Gene Co-expression Network Analysis (WGCNA) pipeline for analyzing transcriptomic changes in response to Abemaciclib treatment. The project identifies highly correlated gene modules and potential biomarkers using data from the GEO dataset **GSE247791**.

## 🧬 Project Overview
The goal of this study was to move beyond simple Differential Gene Expression (DGE) and identify clusters of co-expressed genes (modules) that drive the biological response to CDK4/6 inhibition.

### Key Bioinformatics Workflow:
1.  **Normalization:** VST transformation via `DESeq2`.
2.  **Network Construction:** Scale-free topology selection and blockwise module detection.
3.  **Trait Correlation:** Linking module eigengenes to Abemaciclib treatment.
4.  **Gene Extraction:** Identifying the "Best Module" based on p-value and correlation strength.

---

## 📊 Key Results

### 1. Soft-Threshold Selection
To ensure a scale-free network topology, a soft-thresholding power of **14** was selected, where the $R^2$ reached the 0.8 threshold.

![Soft Threshold Plot](results/figures/Soft_threshold.jpg)

### 2. Module Identification
Gene clustering resulted in several distinct modules. The dendrogram below shows the hierarchical clustering of genes and their assigned module colors.

![Gene Dendrogram](results/figures/Gene_dendrogram_modules.jpg)

### 3. Module–Trait Relationship
The heatmap below illustrates the correlation between identified modules and Abemaciclib treatment. 

* **Top Module:** The **lightgreen** module showed the strongest significant negative correlation ($r = -0.54$, $p = 0.00059$).
* **Biological Significance:** Since Abemaciclib is a cell-cycle inhibitor, this negative correlation suggests the module contains genes involved in proliferation and DNA replication that are suppressed upon treatment.

![Module-Trait Heatmap](results/figures/Module_trait_heatmap.jpg)

---

## 📂 Repository Structure
* `data/`: Contains raw count matrices and metadata.
* `scripts/`: The complete R pipeline (`WGCNA_Analysis.R`).
* `results/`: 
    * `figures/`: JPEG plots of clustering, thresholding, and heatmaps.
    * `tables/`: CSV files including the Module-Trait summary and the specific genes from the Lightgreen module.

## 🛠️ Requirements
* R version 4.x
* Bioconductor packages: `WGCNA`, `DESeq2`, `limma`
* Visualization: `ggplot2`, `pheatmap`
