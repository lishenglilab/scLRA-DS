#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(optparse)

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("--gene_matrix"), type="character", default=NULL, 
                help="Path to gene matrix directory", metavar="character"),
    make_option(c("--transcript_matrix"), type="character", default=NULL,
                help="Path to transcript matrix directory", metavar="character"),
    make_option(c("--output_dir"), type="character", default=".",
                help="Output directory", metavar="character"),
    make_option(c("--species"), type="character", default="Human",
                help="Species (for mitochondrial gene pattern)", metavar="character"),
    make_option(c("--tissue"), type="character", default="Skin",
                help="Tissue type", metavar="character"),
    make_option(c("--api_key"), type="character", default=NULL,
                help="Deprecated: use the LLM_API_KEY environment variable", metavar="character"),
    make_option(c("--project_dir"), type="character", default=".",
                help="Project directory path (for finding Python scripts)", metavar="character"),
    make_option(c("--min_cells"), type="integer", default=10,
                help="Minimum cells expressing a feature [default: %default]"),
    make_option(c("--min_features"), type="integer", default=200,
                help="Minimum detected features per cell [default: %default]"),
    make_option(c("--llm_enabled"), type="character",
                default=ifelse(nzchar(Sys.getenv("LLM_API_KEY")), "true", "false"),
                help="Whether optional LLM annotation is enabled [default: %default]"),
    make_option(c("--llm_base_url"), type="character", default=Sys.getenv("LLM_BASE_URL"),
                help="OpenAI-compatible API base URL"),
    make_option(c("--llm_model"), type="character", default=Sys.getenv("LLM_MODEL", "deepseek-reasoner"),
                help="Model name exposed by the configured API")
  )
  
  opt_parser <- OptionParser(option_list=option_list)
  opts <- parse_args(opt_parser)
  
  # Check required arguments
  if (is.null(opts$gene_matrix) || is.null(opts$transcript_matrix)) {
    print_help(opt_parser)
    stop("Both --gene_matrix and --transcript_matrix are required", call.=FALSE)
  }
  
  return(opts)
}

# Main processing function
process_seurat <- function(obj, matrix_type, species) {
  # Set mitochondrial pattern based on species
  if (tolower(species) == "human") {
    mt_pattern <- "^MT-"
  } else {
    mt_pattern <- "^Mt-"
  }
  
  # Helper function to filter cells
  filter_cells <- function(obj) {
    nCount <- obj@meta.data$nCount_RNA
    nFeature <- obj@meta.data$nFeature_RNA
    
    obj <- AddMetaData(obj, log10(nFeature), col.name = "log_nFeature")
    obj <- AddMetaData(obj, log10(nCount), col.name = "log_nCount")
    
    obj <- subset(
      obj,
      subset = log_nFeature >= median(log10(nFeature)) - 3 * mad(log10(nFeature)) &
        log_nCount >= median(log10(nCount)) - 3 * mad(log10(nCount))
    )
    
    return(obj)
  }
  
  # Standard analysis pipeline
  standard_analysis <- function(obj) {
    # Split by sample if multiple samples
    if (length(unique(obj$orig.ident)) > 1) {
      obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
      
      obj <- NormalizeData(obj) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()
      
      # Try Harmony integration, fall back to regular if not available
      if (require(harmony)) {
        obj <- IntegrateLayers(
          object = obj,
          method = HarmonyIntegration,
          orig.reduction = "pca",
          new.reduction = "harmony",
          verbose = FALSE
        )
        reduction_use <- "harmony"
      } else {
        cat("Harmony not available, using regular PCA\n")
        reduction_use <- "pca"
      }
    } else {
      obj <- NormalizeData(obj) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()
      reduction_use <- "pca"
    }
    
    obj <- FindNeighbors(obj, reduction = reduction_use, dims = 1:30) %>%
      FindClusters(resolution = 1) %>%
      RunUMAP(reduction = reduction_use, dims = 1:30)
    
    # Join layers if split
    if ("RNA" %in% names(obj@assays) && is.list(obj@assays$RNA@layers)) {
      obj <- JoinLayers(obj)
    }
    
    return(obj)
  }
  
  # Process based on matrix type
  if (matrix_type == 'gene') {
    # Calculate mitochondrial percentage
    obj <- PercentageFeatureSet(obj, pattern = mt_pattern, col.name = "percent.mt")
    
    # Run doublet detection if scrubletR is available
    doublet_rate <- NA_real_
    if (require(scrubletR, quietly = TRUE)) {
      tryCatch({
        obj <- scrublet_R(obj)
        doublet_rate <- sum(obj$predicted_doublets) / length(obj$predicted_doublets)
      }, error = function(e) {
        cat("Doublet detection failed:", e$message, "\n")
        obj$predicted_doublets <- FALSE
      })
    } else {
      obj$predicted_doublets <- FALSE
      cat("scrubletR not available, skipping doublet detection\n")
    }
    
    # Filter cells
    nCount <- obj@meta.data$nCount_RNA
    nFeature <- obj@meta.data$nFeature_RNA
    mt <- obj@meta.data$percent.mt
    
    obj <- AddMetaData(obj, log10(nFeature), col.name = "log_nFeature")
    obj <- AddMetaData(obj, log10(nCount), col.name = "log_nCount")
    
    obj <- subset(
      obj,
      subset = log_nFeature >= median(log10(nFeature)) - 3 * mad(log10(nFeature)) &
        log_nCount >= median(log10(nCount)) - 3 * mad(log10(nCount)) &
        percent.mt <= median(mt) + 3 * mad(mt) &
        !predicted_doublets
    )
    
    # Standard analysis
    obj <- standard_analysis(obj)
    obj$doublet_rate <- doublet_rate
    
  } else if (matrix_type == 'transcript') {
    # Filter cells
    obj <- filter_cells(obj)
    
    # Standard analysis
    obj <- standard_analysis(obj)
    
  } else {
    warning("Unknown matrix_type: ", matrix_type, ". Please use 'gene' or 'transcript'.")
  }
  
  return(obj)
}

# Function to run cell annotation with an OpenAI-compatible LLM endpoint.
annotate_cells_llm <- function(markers_file, output_csv, species, tissue,
                               project_dir, base_url, model) {
  if (!nzchar(Sys.getenv("LLM_API_KEY"))) {
    cat("LLM_API_KEY is not set; skipping LLM cell annotation\n")
    return(NULL)
  }
  
  # Set Python script path
  python_script <- file.path(project_dir, "bin", "llm", "cell_annotation.py")
  
  # Check if Python script exists
  if (!file.exists(python_script)) {
    cat("Python script not found at:", python_script, "\n")
    return(NULL)
  }
  
  args <- c(
    shQuote(python_script),
    "--csv_file_path", shQuote(markers_file),
    "--output_csv", shQuote(output_csv),
    "--species_type", shQuote(species),
    "--tissue_type", shQuote(tissue),
    "--base_url", shQuote(base_url),
    "--model", shQuote(model)
  )

  cat("Running LLM cell annotation...\n")
  status <- system2("python", args = args)

  if (length(status) != 1 || is.na(status) || status != 0) {
    cat("LLM annotation command failed with exit status:", status, "\n")
    return(NULL)
  }
  
  if (file.exists(output_csv)) {
    return(read.csv(output_csv))
  } else {
    cat("LLM annotation failed or no output file created\n")
    return(NULL)
  }
}

# Function to calculate QC metrics
calculate_qc_metrics <- function(obj.ge, obj.tr) {
  # Gene UMI stats
  total_gene_umi_per_cell <- colSums(obj.ge@assays$RNA@layers$counts)
  ge_umi_stats <- c(
    mean = mean(total_gene_umi_per_cell),
    min = min(total_gene_umi_per_cell),
    max = max(total_gene_umi_per_cell),
    sd = sd(total_gene_umi_per_cell)
  )
  
  # Genes detected stats
  genes_detected_per_cell <- colSums(obj.ge@assays$RNA@layers$counts > 0)
  genes_detected_stats <- c(
    mean = mean(genes_detected_per_cell),
    min = min(genes_detected_per_cell),
    max = max(genes_detected_per_cell),
    sd = sd(genes_detected_per_cell)
  )
  
  # Transcript UMI stats
  total_tr_umi_per_cell <- colSums(obj.tr@assays$RNA@layers$counts)
  tr_umi_stats <- c(
    mean = mean(total_tr_umi_per_cell),
    min = min(total_tr_umi_per_cell),
    max = max(total_tr_umi_per_cell),
    sd = sd(total_tr_umi_per_cell)
  )
  
  # Transcripts detected stats
  transcripts_detected_per_cell <- colSums(obj.tr@assays$RNA@layers$counts > 0)
  transcripts_detected_stats <- c(
    mean = mean(transcripts_detected_per_cell),
    min = min(transcripts_detected_per_cell),
    max = max(transcripts_detected_per_cell),
    sd = sd(transcripts_detected_per_cell)
  )
  
  # Mitochondrial percentage stats (if available)
  if ("percent.mt" %in% colnames(obj.ge@meta.data)) {
    mito_percentage_per_cell <- obj.ge$percent.mt
    mito_percentage_stats <- c(
      mean = mean(mito_percentage_per_cell),
      min = min(mito_percentage_per_cell),
      max = max(mito_percentage_per_cell),
      sd = sd(mito_percentage_per_cell)
    )
  } else {
    mito_percentage_stats <- c(mean = NA, min = NA, max = NA, sd = NA)
  }
  
  # Doublet rate (if available)
  if ("doublet_rate" %in% colnames(obj.ge@meta.data)) {
    doublet_rate <- unique(obj.ge$doublet_rate)[1]
  } else {
    doublet_rate <- NA
  }
  
  # Gene-transcript correlation
  common_cells <- intersect(colnames(obj.ge), colnames(obj.tr))
  if (length(common_cells) > 0) {
    
    transcripts_detected_common <- colSums(subset(obj.tr,cells = common_cells)@assays$RNA@layers$counts > 0)
    genes_detected_common <- colSums(subset(obj.ge,cells = common_cells)@assays$RNA@layers$counts > 0)
    correlation <- cor(transcripts_detected_common, genes_detected_common)
  } else {
    correlation <- NA
  }
  
  # Create QC data frame with correct column names
  qc_data <- data.frame(
    total_gene_umi_mean = ge_umi_stats['mean'],
    total_gene_umi_min = ge_umi_stats['min'],
    total_gene_umi_max = ge_umi_stats['max'],
    total_gene_umi_std = ge_umi_stats['sd'],
    
    genes_detected_mean = genes_detected_stats['mean'],
    genes_detected_min = genes_detected_stats['min'],
    genes_detected_max = genes_detected_stats['max'],
    genes_detected_std = genes_detected_stats['sd'],
    
    total_tr_umi_mean = tr_umi_stats['mean'],
    total_tr_umi_min = tr_umi_stats['min'],
    total_tr_umi_max = tr_umi_stats['max'],
    total_tr_umi_std = tr_umi_stats['sd'],
    
    transcripts_detected_mean = transcripts_detected_stats['mean'],
    transcripts_detected_min = transcripts_detected_stats['min'],
    transcripts_detected_max = transcripts_detected_stats['max'],
    transcripts_detected_sd = transcripts_detected_stats['sd'],
    
    mt_rna_mean = ifelse(is.na(mito_percentage_stats['mean']), NA, mito_percentage_stats['mean']),
    mt_rna_min = ifelse(is.na(mito_percentage_stats['min']), NA, mito_percentage_stats['min']),
    mt_rna_max = ifelse(is.na(mito_percentage_stats['max']), NA, mito_percentage_stats['max']),
    mt_rna_std = ifelse(is.na(mito_percentage_stats['sd']), NA, mito_percentage_stats['sd']),
    
    gene_transcript_corr = correlation,
    doublet_rate = doublet_rate
  )
  
  # Convert to numeric (some values might be named vectors)
  qc_data[] <- lapply(qc_data, function(x) if(is.numeric(x)) x else as.numeric(x))
  
  return(qc_data)
}

# Main function
main <- function() {
  # Parse arguments
  opts <- parse_arguments()

  llm_enabled_value <- tolower(trimws(opts$llm_enabled))
  if (!llm_enabled_value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--llm_enabled must be true or false", call. = FALSE)
  }
  llm_enabled <- llm_enabled_value %in% c("true", "1", "yes")

  # Backward compatibility for direct callers while keeping the key out of
  # child-process command lines.
  if (!is.null(opts$api_key) && nzchar(opts$api_key) && !nzchar(Sys.getenv("LLM_API_KEY"))) {
    warning("--api_key is deprecated; set LLM_API_KEY instead")
    Sys.setenv(LLM_API_KEY = opts$api_key)
  }
  
  # Create output directory
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Try to load harmony
  if (!require(harmony)) {
    cat("Harmony not available, will use regular integration\n")
  }
  
  # Try to load scrubletR
  if (!require(scrubletR, quietly = TRUE)) {
    cat("scrubletR not available, doublet detection will be skipped\n")
  }
  
  # Load data
  cat("Loading gene matrix from:", opts$gene_matrix, "\n")
  sce.gene.raw <- CreateSeuratObject(
    counts = Read10X(data.dir = opts$gene_matrix),
    min.cells = opts$min_cells,
    min.features = opts$min_features
  )
  
  cat("Loading transcript matrix from:", opts$transcript_matrix, "\n")
  sce.tr.raw <- CreateSeuratObject(
    counts = Read10X(data.dir = opts$transcript_matrix),
    min.cells = opts$min_cells,
    min.features = opts$min_features
  )
  
  # Keep only common cells from raw data
  common_cells_raw <- intersect(colnames(sce.gene.raw), colnames(sce.tr.raw))
  cat("Number of common cells in raw data:", length(common_cells_raw), "\n")
  
  if (length(common_cells_raw) == 0) {
    stop("No common cells found between gene and transcript matrices")
  }
  
  sce.gene.raw <- subset(sce.gene.raw, cells = common_cells_raw)
  sce.tr.raw <- subset(sce.tr.raw, cells = common_cells_raw)
  
  # Gene-level analysis
  cat("\n=== Gene-level analysis ===\n")
  sce.gene.raw$orig.ident <- sub('_[^_]+$', '', colnames(sce.gene.raw))
  sce.gene <- process_seurat(sce.gene.raw, matrix_type = 'gene', species = opts$species)
  
  # Transcript-level analysis
  cat("\n=== Transcript-level analysis ===\n")
  sce.tr.raw$orig.ident <- sub('_[^_]+$', '', colnames(sce.tr.raw))
  sce.tr <- process_seurat(sce.tr.raw, matrix_type = 'transcript', species = opts$species)
  
  # Get common cells after processing
  common_cells_processed <- intersect(colnames(sce.gene), colnames(sce.tr))
  cat("Number of common cells after processing:", length(common_cells_processed), "\n")
  
  if (length(common_cells_processed) == 0) {
    stop("No common cells remaining after processing gene and transcript data")
  }
  
  # Subset both objects to common cells
  cat("Subsetting to common cells for downstream analysis...\n")
  sce.gene <- subset(sce.gene, cells = common_cells_processed)
  sce.tr <- subset(sce.tr, cells = common_cells_processed)
  
  # Find markers for gene-level clusters
  all_cluster_markers <- FindAllMarkers(sce.gene, only.pos = TRUE, logfc.threshold = 0.25)
  write.table(all_cluster_markers, 
              file.path(opts$output_dir, "gene_all_cluster_markers.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Top markers for visualization
  top_markers <- all_cluster_markers %>%
    filter(!grepl('novel', gene)) %>%
    group_by(cluster) %>%
    top_n(3, wt = avg_log2FC) %>%
    dplyr::select(cluster, gene)
  
  # Create gene-level dot plot
  p <- DotPlot(object = sce.gene, features = unique(top_markers$gene)) +
    scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(file.path(opts$output_dir, "gene_cluster_marker_dotplot.pdf"),
         plot = p, width = 20, height = 15)
  
  # Prepare markers for optional LLM annotation
  anno_ge <- all_cluster_markers %>%
    filter(!grepl('novel', gene)) %>%
    group_by(cluster) %>%
    top_n(30, wt = avg_log2FC) %>%
    dplyr::select(cluster, gene)
  
  anno_file <- file.path(opts$output_dir, "anno_ge_for_llm.csv")
  write.table(anno_ge, anno_file, sep = ",", 
              row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # Run LLM annotation if endpoint and API key are provided.
  annotation_applied <- FALSE
  if (llm_enabled && nzchar(Sys.getenv("LLM_API_KEY")) && nzchar(opts$llm_base_url)) {
    cat("Running LLM cell annotation...\n")
    annotation_file <- file.path(opts$output_dir, "cell_ann_by_llm.csv")
    cell_ann <- annotate_cells_llm(
      markers_file = anno_file,
      output_csv = annotation_file,
      species = opts$species,
      tissue = opts$tissue,
      project_dir = opts$project_dir,
      base_url = opts$llm_base_url,
      model = opts$llm_model
    )
    
    if (!is.null(cell_ann)) {
      # Clean annotation
      cell_ann$Cell.Type <- gsub('\\[.*', '', cell_ann$Cell.Type)
      cell_ann$Cluster <- gsub('Cluster', '', cell_ann$Cluster)
      
      # Add annotation to gene-level Seurat object
      sce.gene$cell_anno_llm <- 'Unknown'
      for (ct in cell_ann$Cluster) {
        if (ct %in% as.character(sce.gene$seurat_clusters)) {
          sce.gene$cell_anno_llm[as.character(sce.gene$seurat_clusters) == ct] <- 
            cell_ann$Cell.Type[cell_ann$Cluster == ct]
        }
      }
      
      # Also add annotation to transcript-level Seurat object (using cell barcodes)
      # First, make sure the order of cells matches
      common_cells <- intersect(colnames(sce.gene), colnames(sce.tr))
      sce.tr$cell_anno_llm <- 'Unknown'
      
      if (length(common_cells) > 0) {
        for (cell in common_cells) {
          if (cell %in% colnames(sce.gene)) {
            sce.tr$cell_anno_llm[colnames(sce.tr) == cell] <- sce.gene$cell_anno_llm[colnames(sce.gene) == cell]
          }
        }
      }
      
      # Create UMAP plots for gene-level
      p1 <- DimPlot(sce.gene, group.by = "seurat_clusters")
      p2 <- DimPlot(sce.gene, group.by = "cell_anno_llm")
      p <- p1 / p2
      
      ggsave(file.path(opts$output_dir, "gene_cluster_umap.pdf"),
             plot = p, width = 20, height = 15)
      
      # Find cell type markers (gene-level)
      Idents(sce.gene) <- sce.gene$cell_anno_llm
      celltype_markers <- FindAllMarkers(sce.gene, only.pos = TRUE, logfc.threshold = 0.25)
      write.table(celltype_markers,
                  file.path(opts$output_dir, "gene_celltype_markers.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Top cell type markers (gene-level)
      top_celltype_markers <- celltype_markers %>%
        filter(!grepl('novel', gene)) %>%
        group_by(cluster) %>%
        top_n(5, wt = avg_log2FC) %>%
        dplyr::select(cluster, gene)
      
      p <- DotPlot(object = sce.gene, features = unique(top_celltype_markers$gene)) +
        scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      ggsave(file.path(opts$output_dir, "gene_celltype_markers_dotplot.pdf"),
             plot = p, width = 20, height = 15)
      annotation_applied <- TRUE
    }
  }

  if (!annotation_applied) {
    # Just create basic UMAP for gene-level
    p <- DimPlot(sce.gene)
    ggsave(file.path(opts$output_dir, "gene_cluster_umap.pdf"),
           plot = p, width = 20, height = 15)
  }
  
  # Find transcript-level markers (by cluster)
  all_tr_markers <- FindAllMarkers(sce.tr, only.pos = TRUE, logfc.threshold = 0.25)
  write.table(all_tr_markers,
              file.path(opts$output_dir, "tr_all_cluster_markers.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Top transcript markers (by cluster)
  top_tr_markers <- all_tr_markers %>%
    filter(!grepl('transcript', gene)) %>%
    group_by(cluster) %>%
    top_n(3, wt = avg_log2FC) %>%
    dplyr::select(cluster, gene)
  
  if (nrow(top_tr_markers) > 0) {
    p <- DotPlot(object = sce.tr, features = unique(top_tr_markers$gene)) +
      scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    ggsave(file.path(opts$output_dir, "tr_cluster_marker_dotplot.pdf"),
           plot = p, width = 20, height = 15)
  }
  
  # Find cell type markers for transcript-level (if cell annotation exists)
  if ("cell_anno_llm" %in% colnames(sce.tr@meta.data)) {
    Idents(sce.tr) <- sce.tr$cell_anno_llm
    
    # Skip if all cells have 'Unknown' annotation
    if (length(unique(sce.tr$cell_anno_llm)) > 1) {
      celltype_tr_markers <- FindAllMarkers(sce.tr, only.pos = TRUE, logfc.threshold = 0.25)
      write.table(celltype_tr_markers,
                  file.path(opts$output_dir, "tr_celltype_markers.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      top_celltype_tr_markers <- celltype_tr_markers %>%
        filter(!grepl('transcript', gene)) %>%
        group_by(cluster) %>%
        top_n(5, wt = avg_log2FC) %>%
        dplyr::select(cluster, gene)
      
      if (nrow(top_celltype_tr_markers) > 0) {
        p <- DotPlot(object = sce.tr, features = unique(top_celltype_tr_markers$gene)) +
          scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(file.path(opts$output_dir, "tr_celltype_markers_dotplot.pdf"),
               plot = p, width = 20, height = 15)
      }
    }
  }
  
  # Transcript UMAP
  p <- DimPlot(sce.tr)
  ggsave(file.path(opts$output_dir, "tr_cluster_umap.pdf"),
         plot = p, width = 20, height = 15)
  
  # If cell annotation exists, also create UMAP colored by cell type for transcript-level
  if ("cell_anno_llm" %in% colnames(sce.tr@meta.data) && 
      length(unique(sce.tr$cell_anno_llm)) > 1) {
    p <- DimPlot(sce.tr, group.by = "cell_anno_llm")
    ggsave(file.path(opts$output_dir, "tr_celltype_umap.pdf"),
           plot = p, width = 20, height = 15)
  }
  
  # Calculate QC metrics
  cat("\n=== Calculating QC metrics ===\n")
  qc_data <- calculate_qc_metrics(sce.gene, sce.tr)
  write.table(qc_data, 
              file.path(opts$output_dir, "qc_metrics.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save Seurat objects
  cat("\n=== Saving results ===\n")
  
  # Use absolute paths to ensure correct filenames
  gene_rds_path <- file.path(opts$output_dir, "sce.gene.rds")
  tr_rds_path <- file.path(opts$output_dir, "sce.tr.rds")
  
  saveRDS(sce.gene, gene_rds_path)
  saveRDS(sce.tr, tr_rds_path)
  
  cat("Saved Seurat objects:\n")
  cat("  Gene-level: ", gene_rds_path, "\n")
  cat("  Transcript-level: ", tr_rds_path, "\n")
  cat("  Cells in gene object: ", ncol(sce.gene), "\n")
  cat("  Cells in transcript object: ", ncol(sce.tr), "\n")
  cat("  Common cells: ", length(intersect(colnames(sce.gene), colnames(sce.tr))), "\n")
  
  cat("Analysis completed successfully!\n")
  cat("Output directory:", opts$output_dir, "\n")
}

# Run main function
if (!interactive()) {
  main()
}
