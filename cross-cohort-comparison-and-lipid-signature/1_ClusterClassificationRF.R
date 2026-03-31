#!/usr/bin/env Rscript

source("utils.R")
CONFIG <- load_config()

required_packages <- c(
    "readxl", "ggfortify", "edgeR", "ggpubr", "ggrepel",
    "randomForest", "ggplotify", "RColorBrewer", "pROC",
    "data.table", "caret", "tidyverse", "pheatmap", "cowplot"
)

for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        load_quietly(package, character.only = TRUE)
    }
}

PROJECT_DIR <- CONFIG$global_vars$project_dir
RESULTS_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$results_dir, 'random_forest')

PLOTS_DIR <- file.path(RESULTS_DIR, "plots")
METRICS_DIR <- file.path(RESULTS_DIR, "metrics")
MODELS_DIR <- file.path(RESULTS_DIR, "models")

for (dir in c(RESULTS_DIR, PLOTS_DIR, METRICS_DIR, MODELS_DIR)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}


load_and_prepare_data <- function(scaled_path, meta_path) {
    scaled_df <- read.csv(scaled_path)
    meta <- read.csv(meta_path)
    
    # Split by cohort and remove controls
    split_cohort <- function(cohort_name) {
        cohort_data <- scaled_df[meta[rownames(scaled_df),]$cohort == cohort_name,]
        cohort_meta <- meta[rownames(cohort_data),]
        patients_only <- cohort_data[rownames(cohort_meta[which(!cohort_meta$group %in% c('non-HF','Non HFpEF')),]),]
        return(list(data = patients_only, meta = cohort_meta[which(!cohort_meta$group %in% c('non-HF','Non HFpEF')),]))
    }
    
    became <- split_cohort('BECAME')
    miracle <- split_cohort('MIRACLE')

    return(list(became = became, miracle = miracle, meta = meta))
}

prepare_rf_data <- function(data, meta) {
    df <- data.frame(data)
    df[] <- lapply(df, function(x) if(is.character(x)) as.numeric(x) else x)
    df$group <- meta[rownames(df),'cluster']
    df$group <- as.factor(gsub('B','', df$group))
    return(df)
}

prettify_matrix <- function(numbers_confusion_matrix){
    frequencies_confusion_matrix <- numbers_confusion_matrix
    for(i in 1:ncol(frequencies_confusion_matrix)){
        frequencies_confusion_matrix[,i] = frequencies_confusion_matrix[,i]/sum(frequencies_confusion_matrix[,i])
    }
    final_confusion_matrix <- matrix(nrow = 3, ncol = 3)
    for (i in 1:3) {
        for (j in 1:3) {
            final_confusion_matrix[i, j] <- paste0(round(frequencies_confusion_matrix[i, j]*100,digits=2), "% (",numbers_confusion_matrix[i, j],")")
        }
    }
    return(final_confusion_matrix) 
}

norm_matrix <- function(confusion_matrix, cohort){
    for(i in 1:ncol(confusion_matrix)){
        confusion_matrix[,i]=confusion_matrix[,i]/sum(confusion_matrix[,i])
    } 
    colnames(confusion_matrix) = paste0(cohort, colnames(confusion_matrix))  
    return(confusion_matrix)
}


train_rf_model <- function(data, formula, ntrees = 2000) {
   
    set.seed(1234)

    train_idx <- unlist(
        lapply(split(seq_len(nrow(data)), data$group), function(idx) {
            sample(idx, size = floor(0.7 * length(idx)))
        })
    )

    train_data <- data[train_idx,]
    model <- randomForest(formula, 
                         data = train_data, 
                         ntree = ntrees,
                         importance = TRUE)
    
    return(list(model = model, train_idx = train_idx))
}

evaluate_model <- function(model, data, true_labels, dataset_name) {
    predictions <- predict(model, data, type = "response")
    conf_matrix <- confusionMatrix(predictions, true_labels)
    
    print(dataset_name)
    print(conf_matrix)
    
    # Calculate class-specific metrics
    class_metrics <- data.frame(
        Class = levels(true_labels),
        Sensitivity = conf_matrix$byClass[,"Sensitivity"],
        Specificity = conf_matrix$byClass[,"Specificity"],
        Precision = conf_matrix$byClass[,"Pos Pred Value"],
        F1_Score = conf_matrix$byClass[,"F1"]
    )
    
    # Get variable importance
    variable_importance <- importance(model)
    top_features <- head(sort(variable_importance[,1], decreasing=TRUE), 10)
    
    return(list(
        dataset = dataset_name,
        predictions = predictions,
        confusion_matrix = conf_matrix,
        class_metrics = class_metrics,
        top_features = top_features,
        true_labels = as.factor(true_labels)
    ))
}

create_visualization <- function(became_results, miracle_results) {
    
    became_contingency_table <- table(became_results$predictions, became_results$true_labels)
    miracle_contingency_table <- table(miracle_results$predictions, miracle_results$true_labels)

    became_confusion_matrix <- as.ggplot(pheatmap(
        norm_matrix(became_contingency_table, "B"),
        border_col='black',
        treeheight_col=0,
        treeheight_row=0,
        border_width = 10,
        number_color='black',
        col=brewer.pal(9,'Reds')[-c(8:11)],
        display_numbers=prettify_matrix(became_contingency_table),
        fontsize=20,
        cluster_row=F,
        cluster_col=F,
        silent=T,
        main="BELGIAN Cohort Results"
    ))
    
    miracle_confusion_matrix <- as.ggplot(pheatmap(
        norm_matrix(miracle_contingency_table, "C"),
        border_col='black',
        treeheight_col=0,
        treeheight_row=0,
        border_width = 10,
        number_color='black',
        col=brewer.pal(9,'PuRd')[-c(8:11)],
        display_numbers=prettify_matrix(miracle_contingency_table),
        fontsize=20,
        cluster_row=F,
        cluster_col=F,
        silent=T,
        main="CANADIAN Cohort Results"
    ))
    
    return(plot_grid(became_confusion_matrix, miracle_confusion_matrix, ncol=2))
}

save_results <- function(results) {

    pdf(file.path(PLOTS_DIR, "classification_results.pdf"), width=12, height=5)
    print(results$plots)
    dev.off()
    
    write.csv(
        results$report$became_class_metrics,
        file.path(METRICS_DIR, "became_metrics.csv")
    )
    
    write.csv(
        results$report$miracle_class_metrics,
        file.path(METRICS_DIR, "miracle_metrics.csv")
    )
    
    write.csv(
        data.frame(
            Feature = names(results$report$top_features),
            Importance = results$report$top_features
        ),
        file.path(METRICS_DIR, "top_features.csv")
    )
    
    # Save model
    saveRDS(results$model, file.path(MODELS_DIR, "random_forest_model.rds"))
    
    # Save full results object
    saveRDS(results, file.path(RESULTS_DIR, "full_results.rds"))
    
}

main <- function() {

    data <- load_and_prepare_data(
        file.path(PROJECT_DIR, CONFIG$global_vars$results_dir,'/preprocess/merged_dataset.csv'),
        file.path(PROJECT_DIR, CONFIG$global_vars$results_dir,'/preprocess/merged_dataset_meta.csv')
    )
    
    became_rf_data <- prepare_rf_data(data$became$data, data$became$meta)
    
    rf_results <- train_rf_model(became_rf_data, group~., ntrees = CONFIG$ClusterClassificationRF$ntrees)
    
    became_results <- evaluate_model(
        rf_results$model,
        became_rf_data[-rf_results$train_idx,],
        became_rf_data$group[-rf_results$train_idx],
        "BECAME"
    )
    
    # MIRACLE clusters were done separately and are not expected to match the clusters of BECAME 
    # The random forest is used later simply to reassign MIRACLE cluster to the BECAME clusters
    miracle_results <- evaluate_model(
        rf_results$model,
        data$miracle$data,
        as.factor(gsub('M','', data$miracle$meta$cluster)),
        "MIRACLE"
    )
    
    plots <- create_visualization(became_results, miracle_results)
    
    return(list(
        model = rf_results$model,
        became_results = became_results,
        miracle_results = miracle_results,
        plots = plots
    ))
}


if (!interactive()) {
    
    results <- main()
    save_results(results)
    
    cat("\nResults have been saved to the following directories:\n")
    cat("Plots:", PLOTS_DIR, "\n")
    cat("Metrics:", METRICS_DIR, "\n")
    cat("Models:", MODELS_DIR, "\n")
}