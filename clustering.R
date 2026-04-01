# ------------------------------------------------------------------------------
# HFpEF Lipidomics Clustering and Survival Analysis
#
# Description:
#   This script performs unsupervised clustering of HFpEF patients based on
#   lipidomics data and evaluates associated clinical outcomes.
#
#   Workflow:
#     1. Import lipidomics and survival datasets from an Excel workbook
#     2. Perform Principal Component Analysis (PCA) on scaled lipid data
#     3. Perform hierarchical clustering on principal components (HCPC)
#     4. Perform Kaplan-Meier survival analyses:
#          - Event-free survival
#          - Overall survival
#
# Input:
#   - Excel file: input_data_clustering.xlsx
#       * Sheet 1: Lipidomics matrix (samples x lipids)
#       * Sheet 2: Clinical/survival data with cluster assignments
#
# Output:
#   - PCA and clustering visualization (fviz_cluster)
#   - Kaplan-Meier survival plots (event-free and overall survival)
# ------------------------------------------------------------------------------

# Package imports ----------------------------------------------------------------
library(ggplot2)
library(readxl)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(survival)
library(survminer)

# Input data ---------------------------------------------------------------------
lipid_data <- read_excel("/path/to/input_data_clustering.xlsx", sheet = 1)
surv_data <- read_excel("/path/to/input_data_clustering.xlsx", sheet = 2)

lipid_data <- as.data.frame(lipid_data)
rownames(lipid_data) <- lipid_data[[1]]
lipid_data <- lipid_data[, -1]

surv_data <- as.data.frame(surv_data)
rownames(surv_data) <- surv_data[[1]]
surv_data <- surv_data[, -1]

# PCA and hierarchical clustering ------------------------------------------------
res.pca <- PCA(scale(lipid_data), ncp = 236, graph = TRUE)
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Plot settings ------------------------------------------------------------------
my_colors <- c("#0072B2", "#E69F00", "#009E73")
# blue   = #0072B2
# orange = #E69F00
# green  = #009E73

# Cluster visualization ----------------------------------------------------------
fviz_cluster(
  res.hcpc,
  labelsize = 0,
  show.clust.cent = FALSE,
  ellipse = FALSE,
  pointsize = 4,
  shape = 16,
  alpha = 0.6,
  palette = my_colors,
  main = "",
  ggtheme = theme_light() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 15),
      panel.border = element_rect(
        color = "darkgray",
        fill = NA,
        size = 0.8
      ),
      legend.position = "none"
    )
)

# Cluster sizes
table(surv_data$clust)
# Expected in original analysis:
# 1  2  3
# 30 51 24

# ------------------------------------------------------------------------------
# Survival analysis 1: Event-free survival
# ------------------------------------------------------------------------------
event_free_fit <- survfit(
  Surv(TimeCombinedMonths, Combined) ~ clust,
  data = surv_data
)

event_free_fit

ggsurvplot(
  event_free_fit,
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  surv.scale = "percent",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  legend.title = "",
  xlab = "Time (Months)",
  ylab = "Event-free survival",
  legend.text = element_text(size = 17),
  palette = c("#0072B2", "#E69F00", "#009E73"),
  risk.table.height = 0.20,
  risk.table.fontsize = 5,
  font.x = c(17),
  font.y = c(17),
  font.tickslab = c(17)
)


# ------------------------------------------------------------------------------
# Survival analysis 2: Overall survival
# ------------------------------------------------------------------------------
overall_survival_fit <- survfit(
  Surv(TimeDeath_Nat_Hist_months, Death) ~ clust,
  data = surv_data
)

overall_survival_fit

ggsurvplot(
  overall_survival_fit,
  conf.int = FALSE,
  pval = TRUE,
  risk.table = TRUE,
  surv.scale = "percent",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  legend.title = "",
  xlab = "Time (Months)",
  ylab = "Overall survival",
  legend.text = element_text(size = 17),
  palette = c("#0072B2", "#E69F00", "#009E73", "#660000"),
  risk.table.height = 0.20,
  risk.table.fontsize = 5,
  font.x = c(17),
  font.y = c(17),
  font.tickslab = c(17)
)
