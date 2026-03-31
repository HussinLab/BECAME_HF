# Source utility functions
source("utils.R")

# Load configuration
CONFIG <- load_config()


required_packages <- c(
    "readxl", "glmnet", "ggplotify", "randomForest",
    "ggridges", "RColorBrewer", "ggplotify", "pheatmap",
    "cowplot", "ggplot2", "reshape2", "plyr", "ggpubr"
)

for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        load_quietly(package, character.only = TRUE)
    }
}

PROJECT_DIR <- CONFIG$global_vars$project_dir
OUTPUT_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$results_dir, "figures")
DATA_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$data_dir)
dir.create(OUTPUT_DIR)

load_and_prepare_lipidomic_data <- function(
    excel_path,
    metadata
) {

    message("Reading Excel file...")

    raw_data_new <- tryCatch(
        data.frame(read_xlsx(excel_path, sheet = 2)),
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
        lipidomic_data = became_num_new,
        sample_metadata = lipid_annotations_new,
        study_metadata = metadata,
        lipid_annotations = lipid_annotations_new
    ))
}

prepare_heatmap_data <- function(normalized_data, metadata) {
    
    normalized_data$cluster <- metadata[rownames(normalized_data), 'cluster']
    
    melted_data <- melt(normalized_data)
    aggregated_means <- aggregate(
        melted_data$value,
        by = list(
            cluster = melted_data$cluster,
            variable = melted_data$variable
        ),
        mean
    )
    
    cluster_matrix <- dcast(aggregated_means, cluster ~ variable)
    rownames(cluster_matrix) <- cluster_matrix[, 1]
    cluster_matrix <- cluster_matrix[, 2:ncol(cluster_matrix)]
    
    return(cluster_matrix)
}

create_heatmaps <- function(
    cluster_matrix, 
    filename = '4A_heatmaps.pdf',
    plot_params = list(width = CONFIG$plot_params$heatmap_4A$width, height = CONFIG$plot_params$heatmap_4A$height, dpi = 300)
) {
    # Create metabolite heatmap
    cluster_heatmap <- as.ggplot(pheatmap(
        cluster_matrix,
        clustering_method = 'ward.D',
        border_col = 'white',
        cutree_row = 8,
        fontsize_col = 1e-10,
        fontsize_row = 20,
        treeheight_col = 0,
        silent=T
    ))
    
    # Create correlation heatmap
    correlation_matrix <- cor(t(cluster_matrix))
    correlation_heatmap <- as.ggplot(pheatmap(
        correlation_matrix,
        border_col = 'white',
        clustering_method = 'ward.D',
        cutree_row = 4,
        cutree_col = 4,
        fontsize_col = 20,
        fontsize_row = 20,
        silent=T
    ))
    
    # Combine and save plots
    combined_plot <- plot_grid(
        cluster_heatmap,
        correlation_heatmap,
        rel_widths = c(2, 1)
    )
    
    ggsave(
        file.path(OUTPUT_DIR, filename),
        combined_plot,
        width = plot_params$width,
        height = plot_params$height,
        dpi = plot_params$dpi
    )
    
    return(list(
        cluster_heatmap = cluster_heatmap,
        correlation_heatmap = correlation_heatmap,
        combined_plot = combined_plot
    ))
}

predict_clusters_random_forest <- function(
    normalized_data, 
    metadata,
    model_path
) {
    # Load model
    rf_model <- tryCatch(
        readRDS(model_path),
        error = function(e) stop("Error loading random forest model: ", e$message)
    )
    
    # Filter and predict
    patient_indices <- rownames(subset(metadata, group %in% c('FEP', 'HFpEF')))
    predicted_clusters <- predict(
        rf_model,
        normalized_data[patient_indices, ],
        type = 'response'
    )
    
    # Update metadata and prepare data
    metadata$cohort=as.character(CONFIG$pretty_labels[metadata$cohort])
    metadata$cluster=as.character(CONFIG$pretty_labels[metadata$cluster])
    metadata$predicted_clusters <- predicted_clusters[rownames(metadata)]
    patient_data <- normalized_data[patient_indices, ]
    patient_metadata <- metadata[patient_indices, ]
    
    # Perform PCA
    pca_results <- prcomp(patient_data)
    pca_data <- data.frame(pca_results$x)
    pca_data$cohort <- patient_metadata$cohort
    pca_data$predicted_clusters <- paste0('B',patient_metadata$predicted_clusters)
    pca_data$original_clusters <- patient_metadata$cluster

    # Create PCA plots
    pca_plot_predicted <- ggplot(pca_data, aes(PC1, PC2, color = predicted_clusters)) +
        geom_point(size = 2) +
        facet_wrap(~cohort) +
        theme_bw() +
        labs(title = "PCA by Predicted Clusters",
             color = "Predicted Clusters")+
        scale_color_manual(values = unlist(CONFIG$plot_colors))+
        theme(strip.background = element_blank(),strip.text = element_text(size = 20))+scale_y_reverse()
    
    pca_plot_original <- ggplot(pca_data, aes(PC1, PC2, color = original_clusters)) +
        geom_point(size = 2) +
        facet_wrap(~cohort) +
        theme_bw() +
        labs(title = "PCA by Original Clusters",
             color = "Original Clusters")+
        scale_color_manual(values = unlist(CONFIG$plot_colors))+
        theme(strip.background = element_blank(),strip.text = element_text(size = 10))+scale_y_reverse()
    
    # Create frequency distribution plot
    cluster_frequencies <- data.frame(table(
        predicted_clusters = paste0('B',patient_metadata$predicted_clusters),
        cohort = patient_metadata$cohort
    ))
    print(cluster_frequencies)
    frequency_plot <- ggplot(cluster_frequencies, 
                           aes(predicted_clusters, Freq, fill = predicted_clusters)) +
        geom_col() +
        facet_wrap(~cohort) +
        theme_bw() +
        labs(title = "Cluster Distribution by Cohort",
             x = "Predicted Clusters",
             y = "Frequency",
             fill = "Clusters")+
        scale_fill_manual(values = unlist(CONFIG$plot_colors))+
        theme(strip.background = element_blank(),strip.text = element_text(size = 10))
    
    return(list(
        predicted_clusters = predicted_clusters,
        pca_results = pca_results,
        plots = list(
            pca_predicted = pca_plot_predicted,
            pca_original = pca_plot_original,
            frequency = frequency_plot
        ),
        cluster_frequencies = table(
            predicted_clusters = patient_metadata$predicted_clusters,
            cohort = patient_metadata$cohort
        )
    ))
}

process_became_data <- function(lipidomic_data, lipidomic_meta, metadata) {
    # Filter BECAME cohort
    became_metadata <- subset(metadata, 
                            cohort == 'BECAME' & cluster %in% c('B1', 'B2', 'B3'))
    
    # Subset and transform data
    filtered_data <- lipidomic_data[rownames(became_metadata), ]
    processed_data <- data.frame(
        scale(log2(filtered_data[, grepl('^X', colnames(filtered_data))])
    ))

    became_metadata_full = subset(metadata,cohort == 'BECAME') 
    became_data_full = lipidomic_data[rownames(became_metadata_full),]
    became_data_full = data.frame(scale(log2(became_data_full[, grepl('^X', colnames(became_data_full))])))
    
    return(list(
        data = processed_data,
        metadata = became_metadata,
        lipidomic_meta = lipidomic_meta,
        full_data = became_data_full,
        full_metadata = became_metadata_full
    ))
}

perform_ridge_analysis <- function(lipidomic_data) {
   
    # Load and prepare data
    message("Loading and preparing data...")
    
    # Process BECAME cohort data
    processed_data <- process_became_data(
        lipidomic_data$lipidomic_data,
        lipidomic_data$lipid_annotations,
        lipidomic_data$study_metadata
    )
    
    # Get frequent features
    predictor_features <- get_frequent_predictors(
        CONFIG$FigureCreation$frequency_threshold
    )
    
    # Fit ridge model
    ridge_results <- fit_ridge_model(
        processed_data$data,
        processed_data$metadata,
        predictor_features
    )
    
    # Create visualizations
    plots <- create_analysis_plots(
        processed_data,
        ridge_results,
        predictor_features
    )
    
    return(list(
        model = ridge_results$model,
        predictions = ridge_results$predictions,
        plots = plots,
        processed_data = processed_data,
        features = predictor_features
    ))
}

get_frequent_predictors <- function(threshold = 0.5) {
    predictor_data <- read.csv(file.path(CONFIG$global_vars$results_dir,'/lasso_analysis/metrics/se_lambda_predictors.csv'), row.names = 1)
    return(rownames(subset(predictor_data, se_lambda_predictor_frequency > threshold)))
}

fit_ridge_model <- function(processed_data, metadata, predictor_features) {
    # Prepare response variable
    processed_data$response_var <- as.factor(
        as.numeric(metadata[rownames(processed_data), 'cluster'] == 'B1')
    )
    
    # Prepare model matrices
    feature_matrix <- as.matrix(processed_data[, predictor_features])
    response_vector <- processed_data$response_var
    
    # Fit ridge model
    ridge_model <- cv.glmnet(
        feature_matrix,
        response_vector,
        alpha = 0,
        family = 'binomial'
    )
    
    # Get predictions
    lambda_min <- ridge_model$lambda.min
    predictions <- predict(
        ridge_model,
        s = lambda_min,
        newx = feature_matrix,
        type = 'response'
    )
    
    return(list(
        model = ridge_model,
        predictions = predictions,
        lambda_min = lambda_min
    ))
}

create_analysis_plots <- function(processed_data, ridge_results, predictor_features) {
    # PCA plot
    pca_plot <- create_pca_plot(
        processed_data$data,
        processed_data$metadata,
        predictor_features,
        ridge_results$predictions
    )
    
    # Boxplot
    boxplot <- create_feature_boxplot(
        processed_data$full_data,
        processed_data$full_metadata,
        processed_data$lipidomic_meta,
        predictor_features
    )
    
    # Correlation heatmap
    correlation_plot <- create_correlation_heatmap(
        processed_data$data,
        processed_data$lipidomic_meta,
        predictor_features
    )
    
    return(list(
        pca = pca_plot,
        boxplot = boxplot,
        correlation = correlation_plot
    ))
}

create_pca_plot <- function(data, metadata, features, predictions) {
    pca_results <- prcomp(data[,features])
    pca_dataframe <- data.frame(pca_results$x)
    pca_dataframe$cluster <- metadata[rownames(pca_dataframe), 'cluster']
    pca_dataframe$cluster <- as.character(CONFIG$pretty_labels[pca_dataframe$cluster]) 
    pca_dataframe$score <- predictions[rownames(pca_dataframe), 's1']
    
    # Calculate convex hulls
    hulls <- ddply(pca_dataframe, "cluster", function(df) df[chull(df$PC1, df$PC2), ])
    
    # Set low scores to NA
    # pca_dataframe$score[pca_dataframe$score < 0.5] <- NA

    ggplot(pca_dataframe, aes(PC1, PC2)) +
        theme_bw() +
        geom_polygon(data = hulls, aes(fill = cluster), alpha = 0.2) +
        geom_point(aes(color = score), size = 2) +
        scale_color_gradientn(
            colors = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
        ) +
        scale_fill_manual(values = unlist(CONFIG$plot_colors))+
        labs(color = 'Probability', fill = 'Cluster')
}


create_feature_boxplot <- function(data, metadata, lipidomic_meta, features) {
    # Prepare data for plotting
    plot_data <- data.frame(
        cbind(data[, features], 
              cluster = metadata[rownames(data), 'cluster'])
    )
    
    # Melt data for ggplot
    melted_data <- melt(plot_data)

    # new_lipid_name_df = data.frame(read_xlsx(CONFIG$new_lipid_names),row.names=2)
    melted_data$lipid_id <- lipidomic_meta[as.character(melted_data$variable), 'Compound.Name']
    melted_data$lipid_name <- lipidomic_meta[as.character(melted_data$variable), 'Lipid.ID']
    melted_data$cluster <- as.character(CONFIG$pretty_labels[melted_data$cluster])
    
    # Create boxplot
    ggplot(melted_data, aes(cluster, value, fill = cluster)) +
        geom_point(aes(color = cluster), 
                  position = position_jitter(width = 0.1)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        facet_wrap(~lipid_name, scales = 'free', ncol = 5) +
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(size = 15)
        ) +
        stat_compare_means(
            comparisons = list(c('B1','B2'), c('B1','B3'))
        ) +
        labs(x = "Cluster", y = "Scaled Expression", fill = "Cluster")+
        scale_fill_manual(values = unlist(CONFIG$plot_colors))+
        scale_color_manual(values = unlist(CONFIG$plot_colors))
}


create_correlation_heatmap <- function(data, lipidomic_meta, features) {
    # Calculate correlation matrix

    correlation_matrix <- cor(data[, features])
    
    # Prepare lipid IDs
    lipidomic_meta$Lipid.ID <- make.unique(lipidomic_meta$Lipid.ID)

    rownames(correlation_matrix) <- lipidomic_meta[rownames(correlation_matrix), 'Lipid.ID']
    colnames(correlation_matrix) <- lipidomic_meta[colnames(correlation_matrix), 'Lipid.ID']
    
    # Prepare display numbers
    display_numbers <- round(correlation_matrix, digits = 2)
    display_numbers[display_numbers == 1 | display_numbers <= 0.5] <- ''
    
    # Create heatmap
    pheatmap(
        correlation_matrix,
        border_col = 'white',
        clustering_method = 'ward.D2',
        display_numbers = display_numbers,
        main = "Correlation Heatmap of Selected Features",
        silent=T
    )
}

main <- function() {
    # Create output directory if it doesn't exist
    if (!dir.exists(OUTPUT_DIR)) {
        dir.create(OUTPUT_DIR, recursive = TRUE)
    }
    
    # 1. Load merged dataset
    message("Loading merged dataset...")
    normalized_data <- read.csv(file.path(PROJECT_DIR, CONFIG$global_vars$results_dir,'/preprocess/merged_dataset.csv'))
    metadata <- read.csv(file.path(PROJECT_DIR, CONFIG$global_vars$results_dir,'/preprocess/merged_dataset_meta.csv'))

    lipidomic_data <- load_and_prepare_lipidomic_data(
        file.path(DATA_DIR, CONFIG$global_vars$lipidomic_data_path),
        metadata
    )
    
    # 2. Create and save heatmaps
    message("Creating heatmaps...")
    cluster_matrix <- prepare_heatmap_data(normalized_data, metadata)
    
    rownames(cluster_matrix) <- as.character(CONFIG$pretty_labels[rownames(cluster_matrix)])
    
    heatmap_results <- create_heatmaps(
        cluster_matrix,
        filename = '4A_heatmaps.pdf'
    )
    
    # 3. Perform random forest prediction and visualization
    message("Performing random forest analysis...")
    rf_results <- predict_clusters_random_forest(normalized_data, metadata, file.path(PROJECT_DIR, CONFIG$global_vars$results_dir, CONFIG$FigureCreation$rf_model))
    
    # Save random forest plots
    ggsave(
        file.path(OUTPUT_DIR, '4B_pca_predicted.pdf'),
        rf_results$plots$pca_predicted,
        height = CONFIG$plot_params$pcaplot_4B$height, width = CONFIG$plot_params$pcaplot_4B$width 
    )
    ggsave(
        file.path(OUTPUT_DIR, '4B_pca_original.pdf'),
        rf_results$plots$pca_original,
        height = CONFIG$plot_params$pcaplot_4B$height, width = CONFIG$plot_params$pcaplot_4B$width 
    )
    ggsave(
        file.path(OUTPUT_DIR, '4B_cluster_frequencies.pdf'),
        rf_results$plots$frequency,
        height = CONFIG$plot_params$freqplot_4B$height, width = CONFIG$plot_params$freqplot_4B$width
    )
    
    # 4. Perform ridge regression analysis for visualization 
    message("Performing ridge regression analysis...")
    ridge_results <- perform_ridge_analysis(lipidomic_data)
    
    # Save ridge regression plots
    ggsave(
        file.path(OUTPUT_DIR, '5C_pca_probability.pdf'),
        ridge_results$plots$pca,
        height=5,width=8
    )
    ggsave(
        file.path(OUTPUT_DIR, '5D_feature_boxplot.pdf'),
        ridge_results$plots$boxplot,
        height=8, width=15
    )
    pdf(file.path(OUTPUT_DIR, '5A_correlation_heatmap.pdf'))
    print(ridge_results$plots$correlation)
    dev.off()
    
    # Return all results
    return(list(
        heatmap_results = heatmap_results,
        rf_results = rf_results,
        ridge_results = ridge_results
    ))
}

if (!interactive()) {
    cat("\nStarting figure creation...\n")
    results <- main()
    cat("\nFigure creation completed. Results saved in:", OUTPUT_DIR, "\n")
}