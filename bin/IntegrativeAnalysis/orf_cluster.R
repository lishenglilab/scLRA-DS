#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(tidyverse)
library(Biostrings)
library(data.table)
library(ggplot2)

# Define command line arguments
option_list <- list(
  make_option(c("--seurat_rds"), type="character", default=NULL,
              help="Path to Seurat RDS file (transcript level)", metavar="FILE"),
  make_option(c("--orf_fasta"), type="character", default=NULL,
              help="Path to ORF FASTA file", metavar="FILE"),
  make_option(c("--gtf_file"), type="character", default=NULL,
              help="Path to GTF file", metavar="FILE"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory", metavar="DIR"),
  make_option(c("--resolution"), type="numeric", default=0.4,
              help="Clustering resolution", metavar="NUM"),
  make_option(c("--dims"), type="character", default="1:30",
              help="Dimensions to use for clustering", metavar="STRING")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$seurat_rds) || is.null(opt$orf_fasta) || is.null(opt$gtf_file)) {
  print_help(opt_parser)
  stop("Required arguments: --seurat_rds, --orf_fasta, and --gtf_file must be provided", call.=FALSE)
}

cat("Starting ORF clustering analysis...\n")
cat("Parameters:\n")
cat(paste("  Seurat RDS:", opt$seurat_rds, "\n"))
cat(paste("  ORF FASTA:", opt$orf_fasta, "\n"))
cat(paste("  GTF file:", opt$gtf_file, "\n"))
cat(paste("  Output dir:", opt$output_dir, "\n"))
cat(paste("  Resolution:", opt$resolution, "\n"))
cat(paste("  Dimensions:", opt$dims, "\n"))

# Create output directory
output_dir <- opt$output_dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Load data
cat("Loading Seurat object...\n")
sce.tr <- readRDS(opt$seurat_rds)

cat("Loading ORF sequences...\n")
fasta_sequences <- readAAStringSet(opt$orf_fasta)

cat("Loading GTF file...\n")
gtf <- as.data.frame(rtracklayer::import(opt$gtf_file))
gtf <- gtf[gtf$type == 'transcript', ]
id_map <- unique(gtf[, c('gene_id', 'gene_name', 'transcript_id')])

# 2. Process missing gene names
cat("Processing gene names...\n")
id_map <- id_map %>%
  group_by(gene_id) %>%
  mutate(gene_name = case_when(
    is.na(gene_name) & all(is.na(gene_name)) ~ gene_id,  
    is.na(gene_name) ~ dplyr::first(gene_name[!is.na(gene_name)]),  
    TRUE ~ gene_name  
  )) %>%
  ungroup()

# 3. Process ORF data
cat("Processing ORF data...\n")
headers <- names(fasta_sequences)
sequences <- as.character(fasta_sequences)

orf_df <- data.frame(Header = headers, Sequence = sequences, stringsAsFactors = FALSE)
orf_df$transcript_id <- gsub('lcl\\||:.*', '', orf_df$Header)
orf_df$orf_len <- nchar(orf_df$Sequence)

# Select the longest ORF for each transcript
orf_df <- orf_df %>% 
  group_by(transcript_id) %>% 
  slice_max(orf_len, with_ties = FALSE) %>%
  ungroup()

# Normalize transcript ID format
orf_df$transcript_id <- gsub('_', '-', orf_df$transcript_id)

# Filter transcripts that exist in Seurat object
orf_df <- orf_df[orf_df$transcript_id %in% rownames(sce.tr), ]
orf_df <- orf_df[!duplicated(orf_df$transcript_id), ]
if (nrow(orf_df) == 0) {
  stop("No ORF-bearing transcripts matched the transcript Seurat object")
}

# Assign names to unique ORF sequences
orf_name <- data.frame(Sequence = unique(orf_df$Sequence))
orf_name$orf_name <- paste('ORF', seq_len(nrow(orf_name)), sep = '_')
orf_df <- left_join(orf_df, orf_name, by = "Sequence")
orf_df <- left_join(orf_df, id_map, by = "transcript_id")

# Save ORF mapping file
output_id_map <- file.path(output_dir, "ORF_id_map.txt")
cat("Saving ORF ID mapping to:", output_id_map, "\n")
fwrite(orf_df, output_id_map, sep = '\t', row.names = FALSE)

# 4. Aggregate expression matrix by ORF
cat("Aggregating expression matrix by ORF...\n")
exp <- as.data.frame(sce.tr@assays$RNA@layers$counts)
rownames(exp) <- rownames(sce.tr)
colnames(exp) <- colnames(sce.tr)
exp$transcript_id <- rownames(exp)

merged_df <- exp %>%
  inner_join(orf_df[, c('transcript_id', 'orf_name')], by = "transcript_id")
merged_df <- select(merged_df, transcript_id, orf_name, everything())

setDT(merged_df)
aggregated_dt <- merged_df[, lapply(.SD, sum), 
                           by = orf_name, 
                           .SDcols = !c("transcript_id")]
setDF(aggregated_dt)

# Save aggregated expression matrix
output_orf_exp <- file.path(output_dir, "orf_aggregated_expression.csv")
cat("Saving aggregated ORF expression to:", output_orf_exp, "\n")
fwrite(aggregated_dt, output_orf_exp, sep = ',', row.names = FALSE, quote = FALSE)

# 5. Create ORF-level Seurat object
cat("Creating ORF-level Seurat object...\n")
orf_exp <- aggregated_dt
rownames(orf_exp) <- orf_exp$orf_name
orf_exp$orf_name <- NULL

# Ensure column names match sce.tr
orf_exp <- orf_exp[, colnames(sce.tr)]

sce.orf <- CreateSeuratObject(orf_exp, meta.data = sce.tr@meta.data)

# Layer normalization
cat("Performing layer normalization and integration...\n")
sce.orf[["RNA"]] <- split(sce.orf[["RNA"]], f = sce.orf$orig.ident)

sce.orf <- NormalizeData(sce.orf) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Harmony integration
cat("Running Harmony integration...\n")
library(harmony)
sce.orf <- IntegrateLayers(
  object = sce.orf,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

# Clustering analysis
cat("Performing clustering...\n")
reduction_use <- "harmony"
parse_dims <- function(value) {
  value <- trimws(value)
  if (grepl("^[0-9]+\\s*:\\s*[0-9]+$", value)) {
    bounds <- as.integer(strsplit(gsub("[[:space:]]", "", value), ":", fixed = TRUE)[[1]])
    dims <- seq.int(bounds[1], bounds[2])
  } else if (grepl("^[0-9]+([[:space:]]*,[[:space:]]*[0-9]+)*$", value)) {
    dims <- as.integer(strsplit(gsub("[[:space:]]", "", value), ",", fixed = TRUE)[[1]])
  } else {
    stop("--dims must be a positive range such as '1:30' or a comma-separated list")
  }

  if (length(dims) == 0 || any(is.na(dims)) || any(dims < 1)) {
    stop("--dims must contain positive integers")
  }
  unique(dims)
}
dims_range <- parse_dims(opt$dims)

sce.orf <- FindNeighbors(sce.orf, reduction = reduction_use, dims = dims_range) %>%
  FindClusters(resolution = opt$resolution) %>%
  RunUMAP(reduction = reduction_use, dims = dims_range)

# 6. Generate DimPlots
cat("Generating DimPlots...\n")

# DimPlot by cluster
dimplot_cluster <- DimPlot(sce.orf, label = TRUE, pt.size = 0.5, repel = TRUE) +
  ggtitle(paste("ORF Clusters (Resolution:", opt$resolution, ")")) +
  theme_minimal()

output_dimplot_cluster <- file.path(output_dir, "ORF_DimPlot_clusters.pdf")
pdf(output_dimplot_cluster, width = 10, height = 8)
print(dimplot_cluster)
dev.off()

# DimPlot by sample (if multiple samples exist)
if (length(unique(sce.orf$orig.ident)) > 1) {
  dimplot_sample <- DimPlot(sce.orf, group.by = "orig.ident", pt.size = 0.5) +
    ggtitle("ORF UMAP by Sample") +
    theme_minimal()
  
  output_dimplot_sample <- file.path(output_dir, "ORF_DimPlot_samples.pdf")
  pdf(output_dimplot_sample, width = 10, height = 8)
  print(dimplot_sample)
  dev.off()
  
  # Also save as PNG for quick viewing
  png(file.path(output_dir, "ORF_DimPlot_samples.png"), width = 1000, height = 800)
  print(dimplot_sample)
  dev.off()
}

# Also save cluster plot as PNG
png(file.path(output_dir, "ORF_DimPlot_clusters.png"), width = 1000, height = 800)
print(dimplot_cluster)
dev.off()

# 7. Save results
output_seurat <- file.path(output_dir, "sce.orf.rds")
cat("Saving ORF Seurat object to:", output_seurat, "\n")
saveRDS(sce.orf, output_seurat)

# 8. Generate summary file
cat("Generating summary...\n")
summary_file <- file.path(output_dir, "orf_cluster_summary.txt")
sink(summary_file)
cat("ORF Clustering Analysis Summary\n")
cat("===============================\n")
cat(paste("Timestamp:", Sys.time(), "\n"))
cat(paste("Input Seurat RDS:", opt$seurat_rds, "\n"))
cat(paste("Input ORF FASTA:", opt$orf_fasta, "\n"))
cat(paste("Input GTF:", opt$gtf_file, "\n"))
cat(paste("Number of unique ORFs:", length(unique(orf_df$orf_name)), "\n"))
cat(paste("Number of transcripts with ORF:", nrow(orf_df), "\n"))
cat(paste("ORF Seurat object dimensions:", paste(dim(sce.orf), collapse = " x "), "\n"))
cat(paste("Clusters found:", length(unique(Idents(sce.orf))), "\n"))
cat(paste("Clustering resolution:", opt$resolution, "\n"))
cat(paste("UMAP dimensions used:", opt$dims, "\n"))
cat("\nOutput files:\n")
cat(paste("  ORF ID mapping:", output_id_map, "\n"))
cat(paste("  Aggregated expression:", output_orf_exp, "\n"))
cat(paste("  ORF Seurat object:", output_seurat, "\n"))
cat(paste("  Cluster DimPlot (PDF):", output_dimplot_cluster, "\n"))
cat(paste("  Cluster DimPlot (PNG):", file.path(output_dir, "ORF_DimPlot_clusters.png"), "\n"))
if (length(unique(sce.orf$orig.ident)) > 1) {
  cat(paste("  Sample DimPlot (PDF):", output_dimplot_sample, "\n"))
  cat(paste("  Sample DimPlot (PNG):", file.path(output_dir, "ORF_DimPlot_samples.png"), "\n"))
}
sink()

cat("ORF clustering analysis completed successfully!\n")
