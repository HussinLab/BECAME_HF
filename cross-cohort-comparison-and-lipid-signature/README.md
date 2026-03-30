# Cross-Cohort Comparison and Lipid Signature

Analysis pipeline comparing lipidomic profiles across the Belgian (BECAME) and Canadian (MIRACLE) HFpEF cohorts to identify a minimal lipid signature distinguishing patient clusters.

## Pipeline

Run all steps sequentially with:

```bash
bash runall.sh
```

| Step | Script | Description |
|------|--------|-------------|
| 0 | `0_CleanDataset.R` | Reads raw lipidomic data from both cohorts, merges them, removes batch effects with edgeR, and exports a normalized expression matrix with metadata. |
| 1 | `1_ClusterClassificationRF.R` | Trains a Random Forest on BECAME patient clusters (B1/B2/B3) and applies it to reassign MIRACLE patients into the same cluster space. |
| 2 | `2_MinimalSignatureLasso.R` | Runs repeated LASSO regressions to identify lipids that most frequently distinguish cluster B1 from the rest, producing a ranked predictor list. |
| 3 | `3_CorrelationGraph.R` | Builds a lipid–lipid correlation network from the BECAME HFpEF patients, highlighting the frequent LASSO predictors within the graph. |
| 4 | `4_FigureCreation.R` | Generates manuscript figures: cluster heatmaps, PCA plots with RF-predicted clusters, ridge-regression probability scores, and feature boxplots. |

## Configuration

All parameters (paths, thresholds, plot dimensions, colours) are set in `config/config.yaml`.
