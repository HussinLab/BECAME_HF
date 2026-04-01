## Code repository for the paper "Lipidomics Identifies HFpEF Phenogroups and a High-Risk Metabolic Signature" - The BECAME-HF project

### Clustering analysis

The script `clustering.R` performs unsupervised clustering of patients with heart failure with preserved ejection fraction (HFpEF) using lipidomics data and evaluates associated clinical outcomes through survival analysis.

The workflow integrates:
- Principal Component Analysis (PCA)
- Hierarchical Clustering on Principal Components (HCPC)
- Kaplan-Meier survival analysis

Running the script requires to download the Supplemental data on the Mendeley repository (```input_data_clustering.xlsx```). 

#### Dependencies

```r
install.packages(c(
  "ggplot2",
  "readxl",
  "dplyr",
  "FactoMineR",
  "factoextra",
  "survival",
  "survminer"
))
```

For the consensus clustering analysis using ClustOmics, see package [here](https://github.com/galadrielbriere/ClustOmics).


### Cross-Cohort Comparison and Lipid Signature

Analysis pipeline comparing lipidomic profiles across the Belgian (BECAME-HF1) and Canadian (BECAME-HF2) HFpEF cohorts. The second part of the analysis aims at identifying a minimal lipid signature distinguishing the B1 patient cluster from other subjects in the BECAME-HF1 cohort.

#### Pipeline

Running the pipeline requires to download the Supplemental data on the Mendeley repository. 
All steps can be run individually but they all use the processed data from `0_CleanDataset.R`, which would need to be run at least once. `3_CorrelationGraph.R` also needs `2_MinimalSignatureLasso.R` to be run before because it uses the lipids identified for esthetic of the network graph. `4_FigureCreation.R` requires all scripts to have run once, since it gathers the results and creates the final figure. The recommended way is to run them sequentially from the bash script:

```bash
bash runall.sh
```

| Step | Script | Description |
|------|--------|-------------|
| 0 | `0_CleanDataset.R` | Reads raw lipidomic data from both cohorts, merges them, removes batch effects with limma, and exports a normalized expression matrix with metadata. |
| 1 | `1_ClusterClassificationRF.R` | Trains a Random Forest on BECAME-HF1 patient clusters (B1/B2/B3) and applies it to reassign BECAME-HF2 patients into the same cluster space. |
| 2 | `2_MinimalSignatureLasso.R` | Runs repeated LASSO regressions to identify lipids that most frequently distinguish cluster B1 from the rest of BECAME-HF1 samples, producing a ranked predictor list. |
| 3 | `3_CorrelationGraph.R` | Builds a lipid–lipid correlation network from the BECAME-HF1 HFpEF patients, highlighting the frequent LASSO predictors within the graph. |
| 4 | `4_FigureCreation.R` | Generates manuscript figures: cluster heatmaps, PCA plots with RF-predicted clusters, ridge-regression probability scores, and feature boxplots. |

#### Configuration

All parameters (paths, thresholds, plot dimensions, colours) are set in `config/config.yaml`.

#### Dependencies

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
