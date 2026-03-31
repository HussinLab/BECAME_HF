load_quietly <- function(package) {
    suppressPackageStartupMessages(library(package, character.only = TRUE))
}

load_config <- function(config_dir = "config", config_file = "config.yaml") {
    config_path <- file.path(config_dir, config_file)
    
    if (!dir.exists(config_dir)) {
        stop("Configuration directory not found: ", config_dir)
    }
    
    if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
    }
    
    tryCatch({
        config <- yaml::read_yaml(config_path)        
        return(config)
    }, error = function(e) {
        stop("Error in configuration: ", e$message)
    })
}
