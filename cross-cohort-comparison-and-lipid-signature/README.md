# Cross-Cohort Comparison and Lipid Signature

Analysis pipeline comparing lipidomic profiles across the Belgian (BECAME) and Canadian (MIRACLE) HFpEF cohorts to identify a minimal lipid signature distinguishing patient clusters.

## Pipeline

Run all steps sequentially with:

```bash
bash runall.sh
```

| Step | Script | Description |
|------|--------|-------------|
| 0 | `0_CleanDataset.R` | Reads raw lipidomic data from both cohorts, merges them, removes batch effects with limma, and exports a normalized expression matrix with metadata. |
| 1 | `1_ClusterClassificationRF.R` | Trains a Random Forest on BECAME patient clusters (B1/B2/B3) and applies it to reassign MIRACLE patients into the same cluster space. |
| 2 | `2_MinimalSignatureLasso.R` | Runs repeated LASSO regressions to identify lipids that most frequently distinguish cluster B1 from the rest of BECAME samples, producing a ranked predictor list. |
| 3 | `3_CorrelationGraph.R` | Builds a lipid–lipid correlation network from the BECAME HFpEF patients, highlighting the frequent LASSO predictors within the graph. |
| 4 | `4_FigureCreation.R` | Generates manuscript figures: cluster heatmaps, PCA plots with RF-predicted clusters, ridge-regression probability scores, and feature boxplots. |

## Configuration

All parameters (paths, thresholds, plot dimensions, colours) are set in `config/config.yaml`.

## Dependencies

R 4.3.2 with the following packages:

| Package | Version | | Package | Version |
|---------|---------|-|---------|---------| 
| car | 3.1-3 | | igraph | 2.0.2 |
| caret | 6.0-94 | | limma | 3.58.1 |
| cowplot | 1.1.3 | | pheatmap | 1.0.12 |
| data.table | 1.15.0 | | plyr | 1.8.9 |
| edgeR | 4.0.16 | | pROC | 1.18.5 |
| ggfortify | 0.4.16 | | randomForest | 4.7-1.2 |
| ggplotify | 0.1.2 | | RColorBrewer | 1.1-3 |
| ggpubr | 0.6.0 | | readxl | 1.4.3 |
| ggrepel | 0.9.5 | | reshape2 | 1.4.4 |
| ggridges | 0.5.6 | | tidyverse | 2.0.0 |
| glmnet | 4.1-8 | | yaml | 2.3.8 |
