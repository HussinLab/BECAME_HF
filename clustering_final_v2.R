# ------------------------------------------------------------------------------
# HFpEF lipidomics clustering and survival analysis
#
# Description:
#   This script:
#   1. Loads the clustering input dataset
#   2. Filters the HFpEF subset of interest
#   3. Merges sample metadata with external IDs
#   4. Performs PCA and HCPC clustering on selected lipid variables
#   5. Visualizes clusters
#   6. Runs Kaplan-Meier survival analyses for:
#        - Event-free survival
#        - Overall survival
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
became_short_complete_v2 <- read_excel("./input_clustering.xlsx")
ids <- read.csv("./correspondance_IDs.csv", sep = ",", header = TRUE)


# Filter HFpEF samples -----------------------------------------------------------
hfpEF <- subset(
  became_short_complete_v2,
  filtre_lipido == 1 & groupe_final == 2
)

# Merge with sample correspondence table -----------------------------------------
hfpEF_join <- merge(hfpEF, ids, by = "ID_JTL")

# Set sample names as row names
rownames(hfpEF_join) <- hfpEF_join$Samples

# Select lipid variables ---------------------------------------------------------
lipid_data <- hfpEF_join[, 407:642]

# Remove one feature as in the original script
lipid_data <- lipid_data %>%
  select(-`PE(16:0_22:6)...544`)

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

# Append cluster assignments to original HFpEF data ------------------------------
clustered_data <- cbind(hfpEF, res.hcpc$data.clust)

# Cluster sizes
table(clustered_data$clust)
# Expected in original analysis:
# 1  2  3
# 30 51 24

# ------------------------------------------------------------------------------
# Survival analysis 1: Event-free survival
# ------------------------------------------------------------------------------
event_free_fit <- survfit(
  Surv(TimeCombinedMonths, Combined) ~ clust,
  data = clustered_data
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
  data = clustered_data
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
