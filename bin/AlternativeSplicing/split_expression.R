#!/usr/bin/env Rscript

# split_expression.R
# Split expression matrix into multiple parts for parallel processing

library(Seurat)
library(dplyr)
library(data.table)
library(optparse)

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("--seurat_file"), type="character", default=NULL,
                help="Path to Seurat object (.rds file)", metavar="character"),
    make_option(c("--split_size"), type="integer", default=100,
                help="Number of cells per split file", metavar="integer"),
    make_option(c("--output_dir"), type="character", default="split_exp",
                help="Output directory for split files", metavar="character")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  opts <- parse_args(opt_parser)
  
  # Check required arguments
  if (is.null(opts$seurat_file)) {
    print_help(opt_parser)
    stop("--seurat_file is required", call.=FALSE)
  }
  
  return(opts)
}

main <- function() {
  # Parse arguments
  opts <- parse_arguments()
  
  cat("Alternative Splicing - Expression Matrix Splitting\n")
  cat("=================================================\n")
  cat("Seurat file:", opts$seurat_file, "\n")
  cat("Split size:", opts$split_size, "\n")
  cat("Output directory:", opts$output_dir, "\n")
  
  # Create output directory
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load Seurat object
  cat("\nLoading Seurat object...\n")
  if (!file.exists(opts$seurat_file)) {
    stop("Seurat file not found: ", opts$seurat_file)
  }
  
  sce.tr <- readRDS(opts$seurat_file)
  cat("Loaded Seurat object\n")
  cat("Matrix dimensions:", dim(sce.tr), "\n")
  cat("Number of cells:", ncol(sce.tr), "\n")
  cat("Number of transcripts:", nrow(sce.tr), "\n")
  
  # Check if the object has normalized data
  cat("\nChecking available assay layers...\n")
  if ("data" %in% names(sce.tr@assays$RNA@layers)) {
    cat("Using 'data' layer (normalized counts)\n")
    mat_data <- sce.tr@assays$RNA@layers$data
  } else if ("counts" %in% names(sce.tr@assays$RNA@layers)) {
    cat("Using 'counts' layer (raw counts)\n")
    mat_data <- sce.tr@assays$RNA@layers$counts
  } else {
    cat("Warning: No standard layers found, trying @counts\n")
    mat_data <- sce.tr@assays$RNA@counts
  }
  
  # Ensure matrix has row and column names
  if (is.null(rownames(mat_data))) {
    rownames(mat_data) <- rownames(sce.tr)
  }
  if (is.null(colnames(mat_data))) {
    colnames(mat_data) <- colnames(sce.tr)
  }
  
  num_cols <- ncol(mat_data)
  cat("Total cells in matrix:", num_cols, "\n")
  
  # Calculate number of parts
  parts <- as.integer(num_cols / opts$split_size)
  if (num_cols %% opts$split_size > 0) {
    parts <- parts + 1
  }
  
  cols_per_part <- ceiling(num_cols / parts)
  
  cat("\nSplitting parameters:\n")
  cat("  Number of parts:", parts, "\n")
  cat("  Cells per part:", cols_per_part, "\n")
  cat("  Total cells:", num_cols, "\n")
  
  # Split and save files
  cat("\nSplitting matrix and saving files...\n")
  part_files <- character(0)
  
  for (i in seq_len(parts)) {
    start_col <- (i - 1) * cols_per_part + 1
    end_col <- min(start_col + cols_per_part - 1, num_cols)
    
    if (start_col > num_cols) break
    
    # Extract subset
    sub <- mat_data[, start_col:end_col, drop = FALSE]
    
    # Convert to matrix if needed
    if (!is.matrix(sub)) {
      sub <- as.matrix(sub)
    }
    
    # Create filename
    fn <- file.path(opts$output_dir, paste0("part_", i, ".txt"))
    
    # Save to file
    cat("  Creating part", i, "with", ncol(sub), "cells: ", fn, "\n")
    
    # Use fwrite for efficient writing
    # Include both row names (transcript IDs) and column names (cell barcodes)
    fwrite(as.data.frame(sub), file = fn, sep = "\t", 
           row.names = TRUE, col.names = TRUE, quote = FALSE)
    cat("\nRemoving empty cell from header line...\n")
    system(paste("sed -i '1s/^\\t//'", shQuote(fn)))
    part_files <- c(part_files, fn)
  }
  
  # Save split info
  info_file <- file.path(opts$output_dir, "../split_info.txt")
  info_lines <- c(
    paste("Seurat file:", opts$seurat_file),
    paste("Split size:", opts$split_size),
    paste("Total cells:", num_cols),
    paste("Number of parts:", parts),
    paste("Cells per part:", cols_per_part),
    paste("Part files created:", length(part_files)),
    "",
    "File list:",
    paste("-", part_files, collapse="\n")
  )
  
  writeLines(info_lines, info_file)
  
  cat("\nSplit completed successfully!\n")
  cat("Created", length(part_files), "part files\n")
  cat("Split info saved to:", info_file, "\n")
  cat("All files saved in:", opts$output_dir, "\n")
  
  
  # Return success
  invisible(TRUE)
}

# Run main function with error handling
if (!interactive()) {
  tryCatch({
    main()
  }, error = function(e) {
    cat("ERROR in split_expression.R:", e$message, "\n")
    quit(save = "no", status = 1)
  })
}
