#!/usr/bin/env Rscript

source("utils.R")
CONFIG <- load_config()

required_packages <- c(
    "readxl", "ggfortify", "edgeR", "car", "cowplot", 
    "limma", "tidyverse", "ggrepel", "reshape2", "yaml"
)

for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        load_quietly(package, character.only = TRUE)
    }
}

PROJECT_DIR <- CONFIG$global_vars$project_dir
DATA_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$data_dir)
OUTPUT_DIR <- file.path(PROJECT_DIR, CONFIG$global_vars$results_dir, 'preprocess')
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_data <- function(data_excel_path) {
    tryCatch({
        
        became_lipidomic_data <- data.frame(read_xlsx(data_excel_path, sheet = 2))
        miracle_lipidomic_data <- data.frame(read_xlsx(data_excel_path, sheet = 4))
        
        miracle_lipidomic_meta = data.frame(t(miracle_lipidomic_data[1:2, 7:ncol(miracle_lipidomic_data)]))
        became_lipidomic_meta = data.frame(t(became_lipidomic_data[1:4, 8:ncol(became_lipidomic_data)]))
        
        colnames(miracle_lipidomic_meta) = miracle_lipidomic_meta[1,]
        miracle_lipidomic_meta = miracle_lipidomic_meta[2:nrow(miracle_lipidomic_meta),]
        miracle_lipidomic_meta$samples = rownames(miracle_lipidomic_meta)
        
        colnames(became_lipidomic_meta) = became_lipidomic_meta[1,]
        became_lipidomic_meta = became_lipidomic_meta[2:nrow(became_lipidomic_meta),]
        became_lipidomic_meta$samples = rownames(became_lipidomic_meta)

        miracle_lipidomic_expr = miracle_lipidomic_data[3:nrow(miracle_lipidomic_data),]
        became_lipidomic_expr = became_lipidomic_data[5:nrow(became_lipidomic_data), ]

        rownames(became_lipidomic_expr) = became_lipidomic_expr$Compound.Name

        miracle_lipidomic_meta$cohort = 'MIRACLE'
        became_lipidomic_meta$cohort = 'BECAME'

        became_lipidomic_meta$Cluster[is.na(became_lipidomic_meta$Cluster)] = 'Ctrl'

        miracle_lipidomic_meta$Cluster = paste0('M', miracle_lipidomic_meta$Cluster)
        became_lipidomic_meta$Cluster = paste0('B', became_lipidomic_meta$Cluster)
        
        colnames(became_lipidomic_meta) = c('group','age','sex','cluster','id','cohort')
        colnames(miracle_lipidomic_meta) = c('group','cluster','id','cohort')

        combined_datasets_meta = rbind(
            became_lipidomic_meta[,c('id','group','cluster','cohort')], 
            miracle_lipidomic_meta[,c('id','group','cluster','cohort')]
        )

        lipids_became_order = miracle_lipidomic_expr$Compound.Name..Belgian.cohort.
        miracle_lipids = miracle_lipidomic_expr$Compound.Name..Canadian.cohort.

        miracle_lipidomic_expr = miracle_lipidomic_expr[,grepl('^X',colnames(miracle_lipidomic_expr))]
        print(dim(became_lipidomic_expr))
        # Order BECAME lipids to match MIRACLE lipids
        became_lipidomic_expr = became_lipidomic_expr[lipids_became_order, grepl('^X',colnames(became_lipidomic_expr))]
        became_lipids = rownames(became_lipidomic_expr)
        print(dim(became_lipidomic_expr))
        became_lipidomic_expr = as.data.frame(t(apply(became_lipidomic_expr,2,FUN=function(x){return(as.numeric(x))})))
        miracle_lipidomic_expr = as.data.frame(t(apply(miracle_lipidomic_expr,2,FUN=function(x){return(as.numeric(x))})))
        
        combined_datasets_expr = data.frame(rbind(miracle_lipidomic_expr,became_lipidomic_expr))

        combined_datasets_meta$cluster[combined_datasets_meta$cluster == 'MNA'] <- 'M-CTL'
        combined_datasets_meta$cluster[combined_datasets_meta$cluster == 'BCtrl'] <- 'B-CTL'
        combined_datasets_meta = combined_datasets_meta[rownames(combined_datasets_expr),]

        lipids_correspondance = data.frame(
            cbind(
                miracle_lipids,
                became_lipids,
                colnames(combined_datasets_expr)
            )
        )
        write.csv(lipids_correspondance, file.path(OUTPUT_DIR, "lipids_correspondance.csv"), quote = FALSE, row.names = FALSE)

        return(list(
            became_lipidomic_expr = became_lipidomic_expr,
            became_lipidomic_meta = became_lipidomic_meta,
            miracle_lipidomic_expr = miracle_lipidomic_expr,
            miracle_lipidomic_meta = miracle_lipidomic_meta,
            combined_datasets_expr = combined_datasets_expr,
            combined_datasets_meta = combined_datasets_meta
        ))
    }, error = function(e) {
        stop("Error reading data files: ", e$message)
    })
}

create_separate_pca_plots <- function(miracle_lipidomic_expr, became_lipidomic_expr, combined_datasets_meta) {

    # Function to calculate PCA summary stats
    get_pca_stats <- function(pca) {
        var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
        loadings <- data.frame(pca$rotation)
        top_loadings <- head(sort(abs(loadings$PC1), decreasing = TRUE), 10)
        
        return(list(
            var_explained = var_explained,
            loadings = loadings,
            top_loadings = top_loadings
        ))
    }
    
    # MIRACLE PCA
    miracle_lipidomic_meta <- combined_datasets_meta[rownames(miracle_lipidomic_expr),]
    miracle_lipidomic_expr <- apply(miracle_lipidomic_expr, 2, as.numeric)

    pca_miracle <- prcomp(scale(log2(miracle_lipidomic_expr)))
    miracle_stats <- get_pca_stats(pca_miracle)

    miracle_lipidomic_meta$pretty_clusters <- factor(unlist(CONFIG$pretty_labels[miracle_lipidomic_meta$cluster]))
    
    plot_pca_miracle <- autoplot(pca_miracle, 
                  data = miracle_lipidomic_meta,
                  col = 'pretty_clusters') +
          theme_bw() + 
          scale_color_manual(values = CONFIG$plot_colors) +
          labs(title = 'CANADIAN Cohort',
               x = sprintf("PC1 (%.1f%%)", miracle_stats$var_explained[1]),
               y = sprintf("PC2 (%.1f%%)", miracle_stats$var_explained[2]),
               color = "MIRACLE Clusters")
    
    # BECAME PCA
    became_lipidomic_meta <- combined_datasets_meta[rownames(became_lipidomic_expr),]
    became_lipidomic_expr <- became_lipidomic_expr[rownames(became_lipidomic_meta),]

    became_lipidomic_expr <- apply(became_lipidomic_expr, 2, as.numeric)
    pca_became <- prcomp(scale(log2(became_lipidomic_expr)))
    became_stats <- get_pca_stats(pca_became)

    became_lipidomic_meta$cluster[became_lipidomic_meta$cluster == 'B-CTL'] <- 'B-Non-HF'
    plot_pca_became <- autoplot(pca_became,
                  became_lipidomic_meta,
                  col = 'cluster') +
          scale_color_manual(values = CONFIG$plot_colors) +
          theme_bw() +
          labs(title = 'BELGIAN Cohort',
               x = sprintf("PC1 (%.1f%%)", became_stats$var_explained[1]),
               y = sprintf("PC2 (%.1f%%)", became_stats$var_explained[2]),
               color = "BECAME Clusters")+
        scale_y_reverse()
    
    return(list(
        plots = list(plot_pca_miracle = plot_pca_miracle, plot_pca_became = plot_pca_became),
        stats = list(
            miracle = miracle_stats,
            became = became_stats
        )
    ))
}

process_merged_data <- function(merged, meta) {
    merge_scaled <- log2(merged)
    m <- removeBatchEffect(t(merge_scaled), batch = meta$cohort)
    
    meta$group[meta$group == 'non-HF'] <- 'CONTROL'
    meta$group[meta$group == 'Non HFpEF'] <- 'CONTROL'
    meta$group[meta$group == 'FEP'] <- 'HFpEF'
    
    scaled_df <- data.frame(scale(t(m)))
    rows <- rownames(meta)[rownames(meta) %in% rownames(scaled_df)]
    
    # Function to create PCA plot and get stats
    create_pca_result <- function(data, meta_data, title, color_col, tag="None") {
        pca_result <- prcomp(data)
        var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
        
        if(color_col == 'cohort_label'){
            label_col = 'Cohort'
        }else{
            label_col = 'Disease'
        }

        if(tag == 'Combined'){
            meta_data$group[meta_data$group == 'CONTROL'] <- 'Non-HF'
            plot <- autoplot(pca_result, 
                            data = meta_data, 
                            col = color_col,
                            shape = "cohort_label") + 
                theme_bw() + 
                labs(title = title,
                     x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                     y = sprintf("PC2 (%.1f%%)", var_explained[2]))+
                scale_color_manual(values = CONFIG$plot_colors) + 
                labs(color = label_col,shape = "Cohort")+
                theme(legend.position = "top", legend.direction = "vertical")+scale_shape_manual(values=c(20,3))+scale_y_reverse()
        } else {
            plot <- autoplot(pca_result, 
                            data = meta_data, 
                            col = color_col) + 
                theme_bw() + 
                labs(title = title,
                     x = sprintf("PC1 (%.1f%%)", var_explained[1]),
                     y = sprintf("PC2 (%.1f%%)", var_explained[2]))+
                scale_color_manual(values = CONFIG$plot_colors) + 
                labs(color = label_col)+
                theme(legend.position = "top", legend.direction = "horizontal")
        }
        return(
            list(
                    plot = plot, 
                    pca = pca_result, 
                    var_explained = var_explained
                )
            )
    }
    
    # Create all PCA results
    meta$cohort_label <- unlist(CONFIG$pretty_labels)[meta$cohort]

    pca_concatenated_raw <- create_pca_result(
        merged, 
        meta, 
        'Concatenated', 
        'cohort_label'
    )

    pca_concatenated_log <- create_pca_result(
        merge_scaled, 
        meta, 
        'Log-Transformed', 
        'cohort_label'
    )

    pca_concatenated_batch_effect_removed <- create_pca_result(
        scaled_df, 
        meta, 
        'Removed Batch Effects with EdgeR + Scaling', 
        'cohort_label'
    )

    pca_concatenated_colored_by_groups <- create_pca_result(
        scaled_df, 
        meta, 
        'Colored By Groups', 
        'group'
    )

    pca_concatenated_combined <- create_pca_result(
        scaled_df, 
        meta, 
        'Colored By Groups', 
        'group', 
        tag='Combined'
    )

    plots_list <- list(
        pca_concatenated_raw = pca_concatenated_raw$plot, 
        pca_concatenated_log = pca_concatenated_log$plot, 
        pca_concatenated_batch_effect_removed = pca_concatenated_batch_effect_removed$plot, 
        pca_concatenated_colored_by_groups = pca_concatenated_colored_by_groups$plot,
        pca_concatenated_combined = pca_concatenated_combined$plot
    )
    
    return(list(
        scaled_df = scaled_df,
        rows = rows,
        plots = plots_list,
        pca_results = list(
            pca_concatenated_raw = pca_concatenated_raw, 
            pca_concatenated_log = pca_concatenated_log, 
            pca_concatenated_batch_effect_removed = pca_concatenated_batch_effect_removed, 
            pca_concatenated_colored_by_groups = pca_concatenated_colored_by_groups
        )
    ))
}

save_results <- function(processed_data, combined_meta, pca_dataset_concatenation, pca_individual_datasets, pca_combined_datasets) {
    dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
    
    write.table(processed_data$scaled_df[processed_data$rows, ],
                file.path(OUTPUT_DIR, 'merged_dataset.csv'),
                quote = FALSE, sep = ',')
    
    write.table(combined_meta[processed_data$rows, ],
                file.path(OUTPUT_DIR, 'merged_dataset_meta.csv'),
                quote = FALSE, sep = ',')
    
    # Use plot parameters from config
    pdf(file.path(OUTPUT_DIR, "pca_analysis_results.pdf"), 
        width = CONFIG$plot_params$supplementary_plots$width, 
        height = CONFIG$plot_params$supplementary_plots$height)
    print(pca_dataset_concatenation)
    dev.off()

    pdf(file.path(OUTPUT_DIR, "pca_analysis_results_separate.pdf"), 
        width = CONFIG$plot_params$supplementary_plots$width, 
        height = CONFIG$plot_params$supplementary_plots$height/2)
    print(pca_individual_datasets)
    dev.off()

    pdf(file.path(OUTPUT_DIR, "pca_analysis_results_combined.pdf"), 
        width = CONFIG$plot_params$supplementary_plots$width/2.5, 
        height = CONFIG$plot_params$supplementary_plots$height/2)
    print(pca_combined_datasets)
    dev.off()
    
}

main <- function() {

    # Read and prepare data
    data_files <- read_data(
        file.path(DATA_DIR, CONFIG$global_vars$lipidomic_data_path)
    )

    write.csv(data_files$became_lipidomic_expr, file.path(OUTPUT_DIR, "became_processed.csv"))
    write.csv(data_files$miracle_lipidomic_expr, file.path(OUTPUT_DIR, "miracle_processed.csv"))

    processed_data <- process_merged_data(
        data_files$combined_datasets_expr, 
        data_files$combined_datasets_meta
    )

    pca_plots <- create_separate_pca_plots(
        data_files$miracle_lipidomic_expr, 
        data_files$became_lipidomic_expr, 
        data_files$combined_datasets_meta
    )
    
    pca_dataset_concatenation <- plot_grid(
        processed_data$plots$pca_concatenated_raw,
        processed_data$plots$pca_concatenated_log,
        processed_data$plots$pca_concatenated_batch_effect_removed,
        processed_data$plots$pca_concatenated_colored_by_groups,
        ncol = 2
    )

    pca_individual_datasets <- plot_grid( 
        pca_plots$plots$plot_pca_miracle,
        pca_plots$plots$plot_pca_became,
        ncol = 2
    )
    
    save_results(
        processed_data, 
        data_files$combined_datasets_meta, 
        pca_dataset_concatenation, 
        pca_individual_datasets, 
        processed_data$plots$pca_concatenated_combined
    )
    
    return(list(
        processed_data = processed_data,
        combined_plot = pca_dataset_concatenation
    ))
}

if (!interactive()) {

    results <- main()
    
}