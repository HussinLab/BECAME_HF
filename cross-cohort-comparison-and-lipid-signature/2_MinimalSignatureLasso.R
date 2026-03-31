#!/usr/bin/env Rscript

source("utils.R")
CONFIG <- load_config()

required_packages <- c(
    "readxl", "caret", "glmnet", "ggplotify", "randomForest",
    "RColorBrewer", "ggfortify", "pROC"
)

for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        load_quietly(package, character.only = TRUE)
    }
}

PROJECT_DIR <- CONFIG$global_vars$project_dir
RESULTS_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$results_dir, 'lasso_analysis')
DATA_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$data_dir)

PLOTS_DIR <- file.path(RESULTS_DIR, "plots")
METRICS_DIR <- file.path(RESULTS_DIR, "metrics")
MODELS_DIR <- file.path(RESULTS_DIR, "models")

for (dir in c(RESULTS_DIR, PLOTS_DIR, METRICS_DIR, MODELS_DIR)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}


load_and_prepare_data <- function(excel_path) {

    became_lipidomic_data <- data.frame(read_xlsx(excel_path, sheet = 2))

    became_lipidomic_meta = data.frame(t(became_lipidomic_data[1:4, 8:ncol(became_lipidomic_data)]))

    colnames(became_lipidomic_meta) = became_lipidomic_meta[1,]
    became_lipidomic_meta = became_lipidomic_meta[2:nrow(became_lipidomic_meta),]
    became_lipidomic_meta$samples = rownames(became_lipidomic_meta)

    became_lipidomic_expr = became_lipidomic_data[5:nrow(became_lipidomic_data), ]

    became_lipid_metadata = became_lipidomic_expr[,1:4]
    rownames(became_lipid_metadata) = paste0('X', rownames(became_lipid_metadata))
    became_lipidomic_expr = data.frame(t(became_lipidomic_expr[,grepl('^X', colnames(became_lipidomic_expr))]))

    became_lipidomic_expr[] <- lapply(became_lipidomic_expr, as.numeric)

    mean_lipid_values = colSums(log2(became_lipidomic_expr))/nrow(became_lipidomic_expr)
    keep_lipids = names(mean_lipid_values[mean_lipid_values>quantile(mean_lipid_values,probs=c(CONFIG$global_vars$quantile_threshold))])
    
    print("Removed lowly expressed lipids:")
    print(became_lipid_metadata[names(mean_lipid_values)[!names(mean_lipid_values)%in%keep_lipids],])

    return(list(
        became_lipidomic_expr = became_lipidomic_expr[,keep_lipids],
        became_lipidomic_meta = became_lipidomic_meta,
        lipid_meta = became_lipid_metadata[keep_lipids,],
        removed_lipids = became_lipid_metadata[!rownames(became_lipid_metadata)%in%keep_lipids,]
    ))
}

prepare_model_data <- function(became_lipidomic_expr, meta_path, cohort_name = 'BECAME') {
    
    meta <- read.csv(meta_path)

    became_data <- became_lipidomic_expr[rownames(subset(meta, cohort == cohort_name)),]

    became_data <- data.frame(scale(log2(became_data[,grepl('^X', colnames(became_data))])))
    became_data$response_var <- as.factor(as.numeric(meta[rownames(became_data), 'cluster'] == 'B1'))
    
    return(became_data)
}

train_and_evaluate_model <- function(data, train_indices) {

    # Prepare training data
    x <- data[train_indices, grepl('^X', colnames(data))]
    y <- data[train_indices,]$response_var
    
    # Train models
    cv_model <- cv.glmnet(as.matrix(x), y, alpha = 1, family = 'binomial')
    
    # Get lambdas
    lambda_min <- cv_model$lambda.min
    lambda_1se <- cv_model$lambda.1se
    
    # Fit models with different lambdas
    min_model <- glmnet(as.matrix(x), y, alpha = 1, lambda = lambda_min, family = 'binomial')
    se_model <- glmnet(as.matrix(x), y, alpha = 1, lambda = lambda_1se, family = 'binomial')
    
    # Get predictors
    min_coef <- coef(min_model)[,1]
    se_coef <- coef(se_model)[,1]
    
    # Remove intercepts
    lambda_min_predictors <- names(min_coef[abs(min_coef) > 0])[-1]  
    lambda_se_predictors <- names(se_coef[abs(se_coef) > 0])[-1] 
    
    # Make predictions
    test_x <- as.matrix(data[-train_indices, grepl('^X', colnames(data))])
    preds_min <- predict(min_model, s = lambda_min, newx = test_x, type = 'response')
    preds_se <- predict(se_model, s = lambda_1se, newx = test_x, type = 'response')
    
    # Evaluate predictions
    results <- list()
    
    if (length(unique(round(preds_min))) > 1 && length(unique(round(preds_se))) > 1) {
        cm_min <- confusionMatrix(data[-train_indices,]$response_var, as.factor(round(preds_min)))
        cm_se <- confusionMatrix(data[-train_indices,]$response_var, as.factor(round(preds_se)))
        
        results <- list(
            min_lambda = list(
                accuracy = cm_min$byClass['Balanced Accuracy'],
                predictors = lambda_min_predictors
            ),
            se_lambda = list(
                accuracy = cm_se$byClass['Balanced Accuracy'],
                predictors = lambda_se_predictors
            )
        )
    }
    
    return(results)
}

summarize_results <- function(results_list,lipid_meta,n_iter) {
    
    valid_results <- results_list[!sapply(results_list, is.null)]
    
    # Calculate average accuracies
    min_accuracies <- sapply(valid_results, function(x) x$min_lambda$accuracy)
    se_accuracies <- sapply(valid_results, function(x) x$se_lambda$accuracy)
    
    # Get common predictors
    min_predictors <- table(unlist(lapply(valid_results, function(x) x$min_lambda$predictors)))
    se_predictors <- table(unlist(lapply(valid_results, function(x) x$se_lambda$predictors)))

    min_predictor_frequency = sort(min_predictors, decreasing = TRUE)/n_iter
    se_predictor_frequency = sort(se_predictors, decreasing = TRUE)/n_iter
    lipid_meta[names(min_predictor_frequency),'min_lambda_predictor_frequency'] = min_predictor_frequency
    lipid_meta[names(se_predictor_frequency),'se_lambda_predictor_frequency'] = se_predictor_frequency
    
    return(list(
        min_lambda = list(
            mean_accuracy = mean(unlist(min_accuracies)),
            sd_accuracy = sd(unlist(min_accuracies)),
            predictor_frequency = lipid_meta[order(-lipid_meta$min_lambda_predictor_frequency),c('Compound.Name','Lipid.ID','min_lambda_predictor_frequency')]
        ),
        se_lambda = list(
            mean_accuracy = mean(unlist(se_accuracies)),
            sd_accuracy = sd(unlist(se_accuracies)),
            predictor_frequency = lipid_meta[order(-lipid_meta$se_lambda_predictor_frequency),c('Compound.Name','Lipid.ID','se_lambda_predictor_frequency')]
        )
    ))
}


save_model_results <- function(results) {
    # Save summary metrics
    write.csv(
        data.frame(
            Metric = c("Mean Accuracy", "SD Accuracy"),
            Min_Lambda = c(results$summary$min_lambda$mean_accuracy,
                         results$summary$min_lambda$sd_accuracy),
            SE_Lambda = c(results$summary$se_lambda$mean_accuracy,
                         results$summary$se_lambda$sd_accuracy)
        ),
        file.path(METRICS_DIR, "accuracy_summary.csv")
    )
    
    # Save predictor frequencies
    write.csv(
        results$summary$min_lambda$predictor_frequency,
        file.path(METRICS_DIR, "min_lambda_predictors.csv")
    )
    
    write.csv(
        results$summary$se_lambda$predictor_frequency,
        file.path(METRICS_DIR, "se_lambda_predictors.csv")
    )
    
    # Save individual iteration results
    saveRDS(results$individual_results, 
            file.path(MODELS_DIR, "iteration_results.rds"))
    
    # Create and save visualization of predictor frequencies
    pdf(file.path(PLOTS_DIR, "predictor_frequencies.pdf"))
    par(mfrow=c(2,1))
    barplot(head(results$summary$min_lambda$predictor_frequency$min_lambda_predictor_frequency, 10),
            main="Top 10 Predictors (Min Lambda)",
            las=2)
    barplot(head(results$summary$se_lambda$predictor_frequency$se_lambda_predictor_frequency, 10),
            main="Top 10 Predictors (SE Lambda)",
            las=2)
    dev.off()
    
    # Save full results object
    saveRDS(results, file.path(RESULTS_DIR, "full_results.rds"))
}

main <- function() {

    set.seed(1234) 
    n_iter <- CONFIG$MinimalSignatureLasso$n_iter
    
    # Load data
    data <- load_and_prepare_data(
        file.path(DATA_DIR, CONFIG$global_vars$lipidomic_data_path)
    )
    
    # Prepare model data
    model_data <- prepare_model_data(
        data$became_lipidomic_expr, 
        file.path(PROJECT_DIR, CONFIG$global_vars$results_dir,'/preprocess/merged_dataset_meta.csv')
    )
    
    # Create cross-validation partitions
    cv_indices <- createDataPartition(model_data$response_var, p = 0.7, 
                                    list = FALSE, times = n_iter)
    
    # Train and evaluate models
    results <- vector("list", ncol(cv_indices))
    for(i in 1:ncol(cv_indices)) {
        results[[i]] <- train_and_evaluate_model(model_data, cv_indices[,i])
    }
    
    # Summarize results
    summary <- summarize_results(results, data$lipid_meta, n_iter)
    
    return(list(
        individual_results = results,
        summary = summary,
        removed_lipids = data$removed_lipids
    ))
}

if (!interactive()) {
    # Run analysis
    results <- main()
    
    # Save all results
    save_model_results(results)
    
    # Print information about saved results
    cat("\nAnalysis completed successfully!\n")
    cat("Results have been saved to:", RESULTS_DIR, "\n")
    cat("\nDirectory structure:\n")
    cat("- Metrics:", METRICS_DIR, "\n")
    cat("- Plots:", PLOTS_DIR, "\n")
    cat("- Models:", MODELS_DIR, "\n")

    # Print summary results
    cat("\nRemoved Lipids\n")
    cat("=======================\n\n")
    print(results$removed_lipids$Lipid.ID)

    cat("\nModel Performance Summary\n")
    cat("=======================\n\n")
    
    cat("\nSE Lambda Results:\n")
    cat(sprintf("Mean Accuracy: %.3f (±%.3f)\n", 
                results$summary$se_lambda$mean_accuracy,
                results$summary$se_lambda$sd_accuracy))
    cat("Top 10 Most Frequent Predictors:\n")
    print(head(results$summary$se_lambda$predictor_frequency, 10))
}