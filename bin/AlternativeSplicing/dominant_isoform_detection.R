#!/usr/bin/Rscript
# Dominant Isoform Detection and Analysis (Multi-group comparison support)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
  library(optparse)
  library(doParallel)
  library(ggpubr)
})

# Define command line options
option_list <- list(
  make_option(c("-s", "--seurat_rds"), type="character", default=NULL,
              help="Path to Seurat RDS file", metavar="FILE"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL,
              help="Path to transcript models GTF file", metavar="FILE"),
  make_option(c("-o", "--output_dir"), type="character", default="./DominantIsoform",
              help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("-m", "--stage_mapping"), type="character", default=NULL,
              help="Path to stage mapping file (TSV: sample, stage)", metavar="FILE"),
  make_option(c("-c", "--comparisons"), type="character", default=NULL,
              help="Path to comparisons file (TSV: group1, group2, name)", metavar="FILE"),
  make_option(c("--min_cells"), type="integer", default=10,
              help="Minimum number of cells for gene expression [default: %default]",
              metavar="NUM"),
  make_option(c("--n_cores"), type="integer", default=15,
              help="Number of cores for parallel processing [default: %default]",
              metavar="NUM"),
  make_option(c("--logfc_threshold"), type="numeric", default=1.0,
              help="Log2 fold change threshold for differential usage [default: %default]",
              metavar="NUM"),
  make_option(c("--generate_all_plots"), action="store_true", default=FALSE,
              help="Generate individual plots for each comparison [default: %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$seurat_rds) || is.null(opt$gtf_file)) {
  print_help(opt_parser)
  stop("Missing required arguments: --seurat_rds and --gtf_file")
}

# Assign variables
seurat_rds <- opt$seurat_rds
gtf_file <- opt$gtf_file
output_dir <- opt$output_dir
stage_mapping_file <- opt$stage_mapping
comparisons_file <- opt$comparisons
min_cells <- opt$min_cells
n_cores <- opt$n_cores
logfc_threshold <- opt$logfc_threshold
generate_all_plots <- opt$generate_all_plots

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to log messages
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0("[", timestamp, "] ", msg))
  write(paste0("[", timestamp, "] ", msg), 
        file = file.path(output_dir, "analysis.log"), append = TRUE)
}

# Set output file paths
transcript_usage_file <- file.path(output_dir, "transcript_usage.txt")
dominant_transcripts_file <- file.path(output_dir, "dominant_transcripts_per_cell.txt")
dominant_pct_file <- file.path(output_dir, "dominant_transcripts_pct_comparison.RDS")
dominant_stats_file <- file.path(output_dir, "dominant_transcripts_statistics.txt")
combined_plot_file <- file.path(output_dir, "dominant_transcripts_combined_comparison.pdf")

log_message("Starting dominant isoform analysis...")
log_message(paste("Input Seurat RDS:", seurat_rds))
log_message(paste("Input GTF file:", gtf_file))
log_message(paste("Stage mapping file:", ifelse(is.null(stage_mapping_file), "Not provided", stage_mapping_file)))
log_message(paste("Comparisons file:", ifelse(is.null(comparisons_file), "Not provided", comparisons_file)))
log_message(paste("Output directory:", output_dir))
log_message(paste("Min cells threshold:", min_cells))
log_message(paste("Number of cores:", n_cores))
log_message(paste("LogFC threshold:", logfc_threshold))
log_message(paste("Generate all plots:", generate_all_plots))

# ----------------------------------------------------------------------
# Part 1: Calculate transcript usage and dominant isoforms
# ----------------------------------------------------------------------
log_message("Loading Seurat object...")
sce.tr <- readRDS(seurat_rds)

# Extract expression matrix
log_message("Extracting expression matrix...")
exp <- as.data.frame(sce.tr@assays$RNA@layers$counts)
colnames(exp) <- colnames(sce.tr)
rownames(exp) <- rownames(sce.tr)
exp <- rownames_to_column(exp, 'transcript_id')

log_message("Loading GTF annotation...")
gtf <- as.data.frame(import(gtf_file))
gtf <- gtf[gtf$type == 'transcript', ]
id_map <- unique(gtf[, c('gene_id', 'gene_name', 'transcript_id')])
id_map <- id_map %>%
  group_by(gene_id) %>%
  mutate(gene_name = case_when(
    is.na(gene_name) & all(is.na(gene_name)) ~ gene_id,  
    is.na(gene_name) ~ first(gene_name[!is.na(gene_name)]),  
    TRUE ~ gene_name  
  )) %>%
  ungroup()
id_map <- id_map[id_map$transcript_id %in% exp$transcript_id, ]

# Match transcript IDs
rownames(id_map) <- id_map$transcript_id
id_map <- id_map[exp$transcript_id, ]

# Count isoforms per gene
log_message("Counting isoforms per gene...")
genes <- id_map %>% 
  dplyr::select(transcript_id, gene_id, gene_name) %>% 
  group_by(gene_id) %>% 
  summarise(numIso = n_distinct(transcript_id))

colnames(genes)[1] <- 'annot_gene_id'

if (!identical(exp$transcript_id, id_map$transcript_id)) {
  stop("Transcript ID mismatch between expression matrix and annotation")
}

exp$gene_id <- id_map$gene_id

# Calculate transcript usage
log_message("Calculating transcript usage...")
isocount <- data.frame(
  annot_gene_id = exp$gene_id, 
  annot_transcript_id = exp$transcript_id, 
  isocount = rowMeans(exp %>% dplyr::select(-c('gene_id', 'transcript_id')))
)

genes <- genes %>% 
  left_join(isocount %>% 
              group_by(annot_gene_id) %>% 
              summarise(geneCount = sum(isocount)))

# Filter for genes with multiple isoforms
genes_with_multiple_isoforms <- genes %>% 
  filter(numIso > 1) %>% 
  dplyr::select(annot_gene_id) %>% 
  pull()

exp <- exp[exp$gene_id %in% genes_with_multiple_isoforms, ]
genes <- genes[genes$annot_gene_id %in% genes_with_multiple_isoforms, ]

# Calculate gene-level expression
log_message("Calculating gene-level expression...")
exp_ge <- exp %>% 
  group_by(gene_id) %>% 
  summarise(across(where(is.numeric), sum))

exp <- dplyr::select(exp, transcript_id, gene_id, everything())

# Merge transcript and gene expression
log_message("Merging transcript and gene expression...")
merged_matrix <- exp %>%
  left_join(exp_ge, by = "gene_id", suffix = c("_tr", "_ge"))

sample_names <- colnames(exp)[-(1:2)]
sample_names_tr <- paste0(sample_names, "_tr")
sample_names_ge <- paste0(sample_names, "_ge")

# Calculate transcript usage proportion
log_message("Calculating transcript usage proportions...")
usage_matrix <- merged_matrix %>%
  mutate(across(all_of(sample_names_tr), 
                ~ . / get(paste0(sub("_tr$", "", cur_column()), "_ge")),
                .names = "{col}_usage"))

usage_matrix <- usage_matrix %>%
  dplyr::select(transcript_id, gene_id, ends_with("_usage"))
usage_matrix[is.na(usage_matrix)] <- 0
colnames(usage_matrix)[3:ncol(usage_matrix)] <- gsub('_tr_usage', '', 
                                                     colnames(usage_matrix)[3:ncol(usage_matrix)])

# Save transcript usage matrix
log_message(paste("Saving transcript usage matrix to:", transcript_usage_file))
fwrite(usage_matrix, transcript_usage_file, sep = '\t', row.names = FALSE, quote = FALSE)

# ----------------------------------------------------------------------
# Part 2: Identify dominant transcript per cell
# ----------------------------------------------------------------------
log_message("Identifying dominant transcript per cell...")

# Read transcript usage data
dt <- fread(transcript_usage_file)

# Get metadata column names and cell column names
id_cols <- colnames(dt)[1:2]          # transcript_id, gene_id
cell_cols <- colnames(dt)[-(1:2)]     # All Cell IDs

# Function to process single cell column
get_dominant <- function(usage_vec, transcript_vec) {
  if (all(usage_vec == 0)) return(NA_character_)
  return(transcript_vec[which.max(usage_vec)])
}

# Process each cell in parallel
log_message(paste("Processing", length(cell_cols), "cells..."))

results_list <- lapply(cell_cols, function(cell_id) {
  # Extract current cell's usage
  usage_values <- dt[[cell_id]]
  
  # Filter rows with usage > 0
  temp_dt <- dt[usage_values > 0, .(
    transcript_id = get(id_cols[1]), 
    gene_id = get(id_cols[2]), 
    usage = usage_values[usage_values > 0]
  )]
  
  # Find transcript with max usage per gene
  if (nrow(temp_dt) > 0) {
    dom_dt <- temp_dt[, .(DominantTranscript = transcript_id[which.max(usage)]), 
                      by = .(gene_id)]
    dom_dt[, Cellid := cell_id]
    return(dom_dt)
  } else {
    return(NULL)
  }
})

# Combine results
final_results <- rbindlist(results_list[!sapply(results_list, is.null)])

# Reorder columns and save
setcolorder(final_results, c("gene_id", "Cellid", "DominantTranscript"))

log_message(paste("Saving dominant transcripts to:", dominant_transcripts_file))
fwrite(final_results, dominant_transcripts_file, sep = "\t", quote = FALSE)

# ----------------------------------------------------------------------
# Part 3: Process metadata and stage assignment
# ----------------------------------------------------------------------
log_message("Processing metadata and stage assignment...")

# Reload Seurat object for metadata
sce.tr <- readRDS(seurat_rds)

# Get metadata
meta <- sce.tr@meta.data
meta$cell_id <- rownames(meta)


# Extract sample name from cell_id (remove suffix after underscore)
meta$sample <- sub("_[^_]+$", "", meta$cell_id)

# Load stage mapping if provided
if (!is.null(stage_mapping_file) && file.exists(stage_mapping_file)) {
  log_message(paste("Loading stage mapping from:", stage_mapping_file))
  stage_mapping <- fread(stage_mapping_file)
  
  # Check required columns
  if (!all(c("sample", "stage") %in% colnames(stage_mapping))) {
    stop("Stage mapping file must contain 'sample' and 'stage' columns")
  }
  
  # Merge stage information
  meta <- meta %>%
    left_join(stage_mapping, by = "sample")
  
  # Check for unmapped samples
  unmapped_samples <- unique(meta$sample[is.na(meta$stage)])
  if (length(unmapped_samples) > 0) {
    log_message(paste("Warning: Unmapped samples found:", paste(unmapped_samples, collapse = ", ")))
    log_message("Assigning default stage 'Unknown' to unmapped samples")
    meta$stage[is.na(meta$stage)] <- "Unknown"
  }
}

# Get unique stages
all_stages <- unique(meta$stage[meta$stage != "Unknown"])
log_message(paste("All stages found:", paste(all_stages, collapse = ", ")))

# Save metadata with stage information
meta_file <- file.path(output_dir, "cell_metadata_with_stage.txt")
fwrite(meta, meta_file, sep = "\t")
log_message(paste("Saved metadata to:", meta_file))

# ----------------------------------------------------------------------
# Part 4: Load comparisons configuration
# ----------------------------------------------------------------------
if (!is.null(comparisons_file) && file.exists(comparisons_file)) {
  log_message(paste("Loading comparisons from:", comparisons_file))
  comparisons_df <- fread(comparisons_file)
  
  # Check required columns
  if (!all(c("group1", "group2", "name") %in% colnames(comparisons_df))) {
    stop("Comparisons file must contain 'group1', 'group2', and 'name' columns")
  }
  if (nrow(comparisons_df) == 0) {
    stop("Comparisons file does not contain any comparison rows")
  }
  
  # Validate comparisons
  valid_comparisons <- list()
  for (i in 1:nrow(comparisons_df)) {
    group1 <- comparisons_df$group1[i]
    group2 <- comparisons_df$group2[i]
    name <- comparisons_df$name[i]
    
    if (group1 %in% all_stages && group2 %in% all_stages) {
      valid_comparisons[[name]] <- list(group1 = group1, group2 = group2)
      log_message(paste("Valid comparison:", name, "(", group1, "vs", group2, ")"))
    } else {
      stop(paste(
        "Comparison", name, "(", group1, "vs", group2,
        ") contains a group absent from the stage mapping. Available groups:",
        paste(all_stages, collapse = ", ")
      ))
    }
  }
} else {
  log_message("No comparisons file provided or file not found.")
  log_message("Using all pairwise combinations of stages...")
  
  # Generate all pairwise comparisons
  valid_comparisons <- list()
  if (length(all_stages) >= 2) {
    for (i in 1:(length(all_stages)-1)) {
      for (j in (i+1):length(all_stages)) {
        name <- paste0(all_stages[i], "_vs_", all_stages[j])
        valid_comparisons[[name]] <- list(group1 = all_stages[i], group2 = all_stages[j])
        log_message(paste("Auto-generated comparison:", name))
      }
    }
  } else {
    log_message("Need at least 2 valid stages for comparison. Exiting comparison step.")
    valid_comparisons <- list()
  }
}

# Exit if no valid comparisons
if (length(valid_comparisons) == 0) {
  stop("No valid comparisons to analyze")
} else {
  # ----------------------------------------------------------------------
  # Part 5: Calculate transcript usage percentages for all comparisons
  # ----------------------------------------------------------------------
  log_message("Calculating transcript usage percentages for all comparisons...")
  
  # Read dominant transcripts
  df <- fread(dominant_transcripts_file)
  res <- unique(df[, c('gene_id', 'DominantTranscript')])
  
  # Get normalized expression matrix
  exp_norm <- as.data.frame(sce.tr@assays$RNA@layers$data)
  rownames(exp_norm) <- rownames(sce.tr)
  colnames(exp_norm) <- colnames(sce.tr)
  
  # Set up parallel processing for gene-level calculations
  log_message(paste("Setting up parallel processing with", n_cores, "cores..."))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Calculate percentages for each gene and transcript
  log_message("Processing each gene...")
  all_res_data <- foreach(ge = unique(res$gene_id), .combine = rbind, .packages = "dplyr") %dopar% {
    trs <- res$DominantTranscript[res$gene_id == ge]
    exp_sub <- exp_norm[trs, , drop = FALSE]
    
    gene_res <- data.frame()
    
    for (tr in trs) {
      # Calculate expression percentage for each stage
      stage_pcts <- sapply(all_stages, function(stage) {
        stage_cells <- meta$cell_id[meta$stage == stage]
        if (length(stage_cells) == 0) return(0)
        
        # Cells expressing the gene in this stage
        gene_exp_cells <- colnames(exp_sub)[apply(exp_sub[, stage_cells, drop = FALSE] > 0, 2, any)]
        n_gene_stage <- length(gene_exp_cells)
        
        if (n_gene_stage == 0) return(0)
        
        # Cells expressing this transcript in this stage
        n_tr_stage <- sum(exp_sub[tr, stage_cells] > 0)
        return(n_tr_stage / n_gene_stage)
      })
      
      # Add to results
      gene_res <- rbind(gene_res, data.frame(
        gene_id = ge,
        transcript_id = tr,
        matrix(stage_pcts, nrow = 1, dimnames = list(NULL, paste0(all_stages, "_pct")))
      ))
    }
    
    return(gene_res)
  }
  
  # Stop parallel cluster
  stopCluster(cl)
  
  # Add gene names from annotation
  all_res_data <- all_res_data %>%
    left_join(id_map %>% 
                dplyr::select(gene_id, gene_name) %>% 
                distinct(), 
              by = c("gene_id" = "gene_id"))
  
  # Save full results
  log_message(paste("Saving full percentage comparison results to:", dominant_pct_file))
  saveRDS(all_res_data, dominant_pct_file)
  
  # ----------------------------------------------------------------------
  # Part 6: Perform each comparison and generate results
  # ----------------------------------------------------------------------
  log_message("Performing comparisons and generating results...")
  
  all_stats <- list()
  all_plots <- list()
  
  for (comp_name in names(valid_comparisons)) {
    comp <- valid_comparisons[[comp_name]]
    group1 <- comp$group1
    group2 <- comp$group2
    
    log_message(paste("Processing comparison:", comp_name))
    
    # Create comparison-specific output directory
    comp_dir <- file.path(output_dir, comp_name)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Prepare comparison-specific data
    col1 <- paste0(group1, "_pct")
    col2 <- paste0(group2, "_pct")
    
    # Check if both columns exist
    if (!(col1 %in% colnames(all_res_data) && col2 %in% colnames(all_res_data))) {
      log_message(paste("Warning: Missing data for comparison", comp_name, ". Skipping."))
      next
    }
    
    # Extract relevant data
    comp_data <- all_res_data[, c("gene_id", "transcript_id", "gene_name", col1, col2)]
    
    # Calculate log2 fold change with pseudocount
    comp_data$lgFC <- log2((comp_data[[col2]] + 1e-6) / (comp_data[[col1]] + 1e-6))
    
    # Assign types based on threshold
    comp_data$type <- ifelse(
      comp_data$lgFC > logfc_threshold, group2,
      ifelse(comp_data$lgFC < -logfc_threshold, group1, "notSig")
    )
    
    # Save comparison-specific statistics
    comp_stats_file <- file.path(comp_dir, paste0("statistics_", comp_name, ".txt"))
    write.table(comp_data, comp_stats_file, sep = '\t', row.names = FALSE, quote = FALSE)
    
    # Count types
    type_counts <- table(comp_data$type)
    log_message(paste("  Transcript type counts for", comp_name, ":"))
    for (type in names(type_counts)) {
      log_message(paste("    ", type, ":", type_counts[type]))
    }
    
    # Store statistics
    all_stats[[comp_name]] <- list(
      data = comp_data,
      counts = type_counts,
      total_genes = nrow(comp_data)
    )
    
    # Generate comparison-specific plot if requested
    if (generate_all_plots) {
      comp_plot_file <- file.path(comp_dir, paste0("comparison_plot_", comp_name, ".pdf"))
      
      p <- ggplot(comp_data, aes_string(x = col1, y = col2, color = "type")) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c(
          setNames(c("red", "blue", "grey"), c(group1, group2, "notSig"))
        )) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        labs(
          x = paste(group1, "Percentage"),
          y = paste(group2, "Percentage"),
          title = paste("Dominant Transcript Usage:", comp_name),
          subtitle = paste("LogFC threshold =", logfc_threshold),
          color = "Type"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "right"
        )
      
      ggsave(plot = p, filename = comp_plot_file, width = 10, height = 8)
      log_message(paste("  Saved plot to:", comp_plot_file))
    }
    
    # Create a smaller plot for combined figure
    p_small <- ggplot(comp_data, aes_string(x = col1, y = col2, color = "type")) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c(
        setNames(c("red", "blue", "grey"), c(group1, group2, "notSig"))
      )) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      labs(
        x = group1,
        y = group2,
        title = comp_name,
        color = "Type"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "none"
      )
    
    all_plots[[comp_name]] <- p_small
  }
  
  # ----------------------------------------------------------------------
  # Part 7: Create combined summary plot
  # ----------------------------------------------------------------------
  if (length(all_plots) > 0) {
    log_message("Creating combined summary plot...")
    
    # Create summary statistics table
    summary_stats <- data.frame()
    for (comp_name in names(all_stats)) {
      stats <- all_stats[[comp_name]]
      for (type in names(stats$counts)) {
        summary_stats <- rbind(summary_stats, data.frame(
          Comparison = comp_name,
          Type = type,
          Count = as.numeric(stats$counts[type]),
          Percentage = round(as.numeric(stats$counts[type]) / stats$total_genes * 100, 2)
        ))
      }
    }
    
    # Create summary bar plot
    summary_bar <- ggplot(summary_stats, aes(x = Comparison, y = Percentage, fill = Type)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("AD" = "red", "Normal" = "blue", "KO" = "green", "notSig" = "grey")) +
      labs(
        title = "Dominant Transcript Type Distribution",
        x = "Comparison",
        y = "Percentage (%)",
        fill = "Type"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Arrange plots
    if (length(all_plots) == 1) {
      combined_plot <- all_plots[[1]] / summary_bar
    } else if (length(all_plots) == 2) {
      combined_plot <- (all_plots[[1]] | all_plots[[2]]) / summary_bar
    } else if (length(all_plots) == 3) {
      combined_plot <- (all_plots[[1]] | all_plots[[2]] | all_plots[[3]]) / summary_bar
    } else {
      # For more than 3 comparisons, arrange in a grid
      ncol <- ceiling(sqrt(length(all_plots)))
      plot_grid <- wrap_plots(all_plots, ncol = ncol)
      combined_plot <- plot_grid / summary_bar
    }
    
    # Save combined plot
    ggsave(plot = combined_plot, filename = combined_plot_file, 
           width = 12, height = 10)
    log_message(paste("Saved combined plot to:", combined_plot_file))
  }
  
  # ----------------------------------------------------------------------
  # Part 8: Save comprehensive statistics
  # ----------------------------------------------------------------------
  log_message("Saving comprehensive statistics...")
  
  # Create summary table
  summary_table <- data.frame()
  for (comp_name in names(all_stats)) {
    comp <- valid_comparisons[[comp_name]]
    stats <- all_stats[[comp_name]]
    
    for (type in names(stats$counts)) {
      summary_table <- rbind(summary_table, data.frame(
        Comparison = comp_name,
        Group1 = comp$group1,
        Group2 = comp$group2,
        Type = type,
        Count = as.numeric(stats$counts[type]),
        Percentage = round(as.numeric(stats$counts[type]) / stats$total_genes * 100, 2)
      ))
    }
  }
  
  # Save summary statistics
  write.table(summary_table, dominant_stats_file, sep = '\t', row.names = FALSE, quote = FALSE)
  log_message(paste("Saved comprehensive statistics to:", dominant_stats_file))
}

# ----------------------------------------------------------------------
# Save session info
# ----------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "session_info.txt"))

writeLines(capture.output(opt), 
           file.path(output_dir, "command_args.txt"))

log_message("Dominant isoform analysis completed successfully!")
