# config.R - Configuration Management
# Load and validate pipeline configuration

#' Load Configuration from YAML
#' 
#' @param config_file Path to YAML config file
#' @return List containing all configuration parameters
#' @export
load_config <- function(config_file = "config/config.yaml") {
  if (!require("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with: install.packages('yaml')")
  }
  
  if (!file.exists(config_file)) {
    stop(sprintf("Configuration file not found: %s", config_file))
  }
  
  config <- yaml::yaml.load_file(config_file)
  
  # Validate required sections
  required_sections <- c("project", "input", "database", "params", "output")
  missing <- required_sections[!required_sections %in% names(config)]
  if (length(missing) > 0) {
    stop(sprintf("Missing config sections: %s", paste(missing, collapse = ", ")))
  }
  
  # Convert relative paths to absolute
  config <- normalize_paths(config, dirname(config_file))
  
  # Validate file existence
  validate_inputs(config)
  
  return(config)
}

#' Normalize Paths in Configuration
#' 
#' @param config Configuration list
#' @param base_dir Base directory for relative paths
#' @return Config with absolute paths
#' @keywords internal
normalize_paths <- function(config, base_dir) {
  # Helper to make paths absolute
  make_abs <- function(x) {
    if (is.character(x) && !grepl("^/|^~", x)) {
      return(normalizePath(file.path(base_dir, "..", x), mustWork = FALSE))
    }
    return(x)
  }
  
  # Process input files
  config$input <- lapply(config$input, make_abs)
  
  # Process database files
  config$database <- lapply(config$database, make_abs)
  
  # Process output directories
  config$output <- lapply(config$output, make_abs)
  
  return(config)
}

#' Validate Input Files Exist
#' 
#' @param config Configuration list
#' @export
validate_inputs <- function(config) {
  # Check input files
  input_files <- unlist(config$input)
  missing <- input_files[!file.exists(input_files)]
  if (length(missing) > 0) {
    msg <- sprintf("Missing input files:\n%s", paste("  -", missing, collapse = "\n"))
    warning(msg, immediate. = TRUE)
  }
  
  # Check database files
  db_files <- unlist(config$database)
  missing <- db_files[!file.exists(db_files)]
  if (length(missing) > 0) {
    msg <- sprintf("Missing database files:\n%s", paste("  -", missing, collapse = "\n"))
    stop(msg)
  }
  
  invisible(TRUE)
}

#' Get File Path from Config
#' 
#' Safe accessor for config file paths
#' @param config Configuration list
#' @param type Type of file: "input", "database", "output"
#' @param name Name of the file in the config
#' @return File path string
#' @export
get_file_path <- function(config, type, name) {
  if (!type %in% c("input", "database", "output")) {
    stop("Type must be 'input', 'database', or 'output'")
  }
  
  path <- config[[type]][[name]]
  if (is.null(path)) {
    stop(sprintf("Config entry not found: %s.%s", type, name))
  }
  
  return(path)
}

#' Create Output Directories
#' 
#' @param config Configuration list
#' @export
create_output_dirs <- function(config) {
  dirs <- c(
    config$output$results_dir,
    config$output$go_dir
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created directory: %s", dir))
    }
  }
  
  # Create GO subdirectories
  for (onto in config$params$go$ontologies) {
    onto_dir <- file.path(config$output$go_dir, onto)
    if (!dir.exists(onto_dir)) {
      dir.create(onto_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  invisible(TRUE)
}