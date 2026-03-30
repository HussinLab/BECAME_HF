#!/usr/bin/env Rscript

# Source utility functions
source("utils.R")

# Load configuration
CONFIG <- load_config()

# Load required libraries
required_packages <- c(
    "readxl", "caret", "glmnet", "ggplotify",
    "igraph", "RColorBrewer", "ggfortify", "pROC"
)

for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        load_quietly(package, character.only = TRUE)
    }
}

# Define paths
PROJECT_DIR <- CONFIG$global_vars$project_dir
DATA_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$data_dir)
RESULTS_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$results_dir)
OUTPUT_DIR <- file.path(RESULTS_DIR, "network_analysis")
dir.create(OUTPUT_DIR, showWarnings = FALSE)


read_and_prepare_data <- function(data_excel_path, meta_path) {
    
    # Read metadata
    meta <- read.csv(meta_path)
    
    raw_data_new <- tryCatch(
        data.frame(read_xlsx(data_excel_path, sheet = 2)),
        error = function(e) stop("Error reading Excel file: ", e$message)
    )

    became_meta_new = data.frame(t(raw_data_new[1:4, 8:ncol(raw_data_new)]))
    colnames(became_meta_new) = became_meta_new[1,]
    became_meta_new = became_meta_new[2:nrow(became_meta_new),]
    became_meta_new$samples = rownames(became_meta_new)

    became_num_new = raw_data_new[5:nrow(raw_data_new), ]
    lipid_annotations_new = became_num_new[,!grepl('^X', colnames(became_num_new))]
    rownames(lipid_annotations_new) = paste0('X', rownames(lipid_annotations_new))

    became_num_new = data.frame(t(became_num_new[,grepl('^X', colnames(became_num_new))]))
    became_num_new[] <- lapply(became_num_new, as.numeric)

    return(list(
        became_num = became_num_new,
        became_meta = lipid_annotations_new,
        meta = meta
    ))
}

prepare_network_data <- function(became_num, meta, became_meta) {

    # Subset and prepare data
    became_data <- became_num[rownames(subset(meta, cohort == 'BECAME' & group == 'HFpEF')),]

    became_data_expr = became_data[,grepl('^X', colnames(became_data))]
    mean_lipid_values = colSums(log2(became_data_expr))/nrow(became_data_expr)
    keep_lipids = names(mean_lipid_values[mean_lipid_values>quantile(mean_lipid_values,probs=c(CONFIG$global_vars$quantile_threshold))])
    became_data_expr = became_data_expr[,keep_lipids]

    became_data <- data.frame(scale(log2(became_data_expr)))
    
    # Calculate correlation matrix
    correlation_matrix <- cor(became_data)
    color_coding_new = data.frame(read.csv(CONFIG$CorrelationGraph$color_coding_path))
    
    lipid_ids <- became_meta[rownames(correlation_matrix), c('Compound.Name','Lipid.suclasses','Lipid.ID')]
    lipid_ids$row_id = rownames(lipid_ids)

    color_conversion_df <- merge(lipid_ids, color_coding_new, by = "Lipid.suclasses",sort = FALSE)

    rownames(color_conversion_df) = color_conversion_df$row_id
    color_conversion_df = color_conversion_df[rownames(correlation_matrix),]

    class_colors_df = unique(color_conversion_df[,c('Abbreviations','color')])
    class_colors_vec = class_colors_df$color
    names(class_colors_vec) = class_colors_df$Abbreviations

    node_colors_new = class_colors_vec[color_conversion_df$Abbreviations]

    return(list(
        became_data = became_data,
        correlation_matrix = correlation_matrix,
        col_vector = class_colors_vec,
        node_colors = node_colors_new,
        became_meta = became_meta
    ))
}

create_network_plot <- function(became_data, correlation_matrix, class_colors_vec, node_colors, became_meta, frequent_features_df) {

    rownames(correlation_matrix) = became_meta[rownames(correlation_matrix),'Lipid.ID']
    colnames(correlation_matrix) = became_meta[colnames(correlation_matrix),'Lipid.ID']
    colnames(became_data) = became_meta[colnames(became_data),'Lipid.ID']

    frequent_features <- became_meta[rownames(subset(frequent_features_df, 
                                       se_lambda_predictor_frequency > CONFIG$CorrelationGraph$frequency_threshold)),'Lipid.ID']
    
    very_frequent_features <- became_meta[rownames(subset(frequent_features_df, 
                                       se_lambda_predictor_frequency > .5)),'Lipid.ID']
                                           # Prepare network data
    correlation_matrix[correlation_matrix < CONFIG$CorrelationGraph$correlation_threshold] <- 0
    network <- graph_from_adjacency_matrix(abs(correlation_matrix), weighted = T, mode = "undirected", diag = F)
    
    # Color scale for edges
    c_scale <- colorRamp(rev(c('darkred', 'red', 'darkorange1', 'yellow', 'grey90', 'grey90')))
    edge_colors <- apply(c_scale(E(network)$weight), 1, 
                        function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
    
    # Layout coordinates
    pca_coords <- as.matrix(prcomp(t(became_data[,rownames(correlation_matrix)]))$x[,1:2])
    co <- layout_with_fr(network, coords = pca_coords)
    
    # Create plot
    pdf(file.path(OUTPUT_DIR, "network_plot_with_labels.pdf"), width = 15, height = 10)
    par(mfrow = c(1, 2))
    border_col = c('2'='red','1'='black','0'='grey60')
    border_col_vec = rep(0, nrow(correlation_matrix))
    border_col_vec = border_col_vec + as.numeric(rownames(correlation_matrix)%in%frequent_features)
    border_col_vec = border_col_vec + as.numeric(rownames(correlation_matrix)%in%very_frequent_features)
    border_col_vec = border_col[as.character(border_col_vec)]

    vertex_size_vec = frequent_features_df$se_lambda_predictor_frequency
    names(vertex_size_vec) = frequent_features
    vertex_size_vec = vertex_size_vec[rownames(correlation_matrix)]
    vertex_size_vec[is.na(vertex_size_vec)] = 0.001

    label_size = c('2'=1,'1'=1,'0'=1e-10)
    label_col_vec = rep(0, nrow(correlation_matrix))
    label_col_vec = label_col_vec + as.numeric(rownames(correlation_matrix)%in%frequent_features)
    label_col_vec = label_col_vec + as.numeric(rownames(correlation_matrix)%in%very_frequent_features)
    label_size_vec = label_size[as.character(label_col_vec)]

    # Main network plot
    plot(network,
         layout = co,
         vertex.color = adjustcolor(node_colors, alpha.f = .7),  # Use node_colors instead of col_vector
         vertex.frame.color = border_col_vec,
         vertex.shape = "circle",
         vertex.size = vertex_size_vec*5 + 2,
         vertex.label.color = 'black',
         vertex.label.font = 1,
         vertex.label.cex = label_size_vec,
         vertex.label.degree = 1,
         vertex.label.dist = .7,
         edge.width = scale(1-edge.betweenness(network)),
         edge.color = adjustcolor(edge_colors, alpha.f = .3),
         edge.curved = 0.3
    )
    
    # Legend plot
    par(mar = c(5, 5, 2, 5))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("top", 
           legend = names(class_colors_vec), 
           fill = class_colors_vec, 
           title = "Categories",
           ncol = 3, 
           cex = 0.6, 
           xpd = TRUE)
    
    dev.off()

    pdf(file.path(OUTPUT_DIR, "network_plot_no_labels.pdf"), width = 15, height = 10)
    par(mfrow = c(1, 2))

    label_size = c('2'=1e-10,'1'=1e-10,'0'=1e-10)
    label_col_vec = rep(0, nrow(correlation_matrix))
    label_col_vec = label_col_vec + as.numeric(rownames(correlation_matrix)%in%frequent_features)
    label_col_vec = label_col_vec + as.numeric(rownames(correlation_matrix)%in%very_frequent_features)
    label_size_vec = label_size[as.character(label_col_vec)]

    # Main network plot
    plot(network,
         layout = co,
         vertex.color = adjustcolor(node_colors, alpha.f = .7),  # Use node_colors instead of class_colors_vec
         vertex.frame.color = border_col_vec,
         vertex.shape = "circle",
         vertex.size = vertex_size_vec*5 + 2,
         vertex.label.color = 'black',
         vertex.label.font = 1,
         vertex.label.cex = label_size_vec,
         vertex.label.degree = 1,
         vertex.label.dist = .7,
         edge.width = scale(1-edge.betweenness(network)),
         edge.color = adjustcolor(edge_colors, alpha.f = .3),
         edge.curved = 0.3
    )
    
    # Legend plot
    par(mar = c(5, 5, 2, 5))
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("top", 
           legend = names(class_colors_vec), 
           fill = class_colors_vec, 
           title = "Categories",
           ncol = 3, 
           cex = 0.6, 
           xpd = TRUE)
    
    dev.off()
}

main <- function() {
    # Read data
    data <- read_and_prepare_data(
        file.path(DATA_DIR, CONFIG$global_vars$lipidomic_data_path),
        file.path(RESULTS_DIR, "preprocess/merged_dataset_meta.csv")
    )
    
    # Prepare network data
    network_data <- prepare_network_data(data$became_num, data$meta, data$became_meta)
    
    # Read frequent features
    frequent_features_df <- read.csv(
        file.path(RESULTS_DIR, "lasso_analysis/metrics/se_lambda_predictors.csv"),
        row.names = 1
    )
    
    # Create and save network plot
    create_network_plot(
        network_data$became_data,
        network_data$correlation_matrix, 
        network_data$col_vector,
        network_data$node_colors,
        network_data$became_meta, 
        frequent_features_df
    )
    
    return(network_data)
}

###################
# Run Analysis
###################

if (!interactive()) {
    set.seed(1234)
    results <- main()
    cat("\nNetwork analysis completed. Results saved in:", OUTPUT_DIR, "\n")
}